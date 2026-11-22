//! Solids for item Line 8 mounting recipes. Access bores deliberately
//! cross interrupted stock, so they use a stepped cutter with retained proof.
use crate::{GeometryBackend, ModelError, error};
use multicad_core::{
    document::{Document, EntityId},
    model::{
        self, CountersunkScrewHeadParameters, CrownedScrewHeadParameters,
        ITEM_LINE_8_GROOVE_LIP_DEPTH_MM, MechanicalPartParameters, PROFILE_ACCESS_FROM_END_MM,
        PROFILE_ACCESS_RADIUS_MM, PROFILE_BARREL_DEPTH_MM, PROFILE_PILOT_DEPTH_MM,
        PROFILE_PILOT_RADIUS_MM, PathSegment, ProfileAttachment, ProfileConnectionFrames,
        ProfileConnectionKind, ProfileConnectionNut, ProfileOutputKind, ProfilePartRole,
        ResolvedMechanicalPart, ResolvedThread, ScrewDrive, ScrewHead, ScrewLengthInput,
        ThreadMode, ThreadParameters, length,
    },
    parameters::EvaluatedParameters,
    references::{AnalyticGeometry, Origin, ReferenceFrame, TopologySnapshot},
};
use multicad_kernel::{
    BooleanOperation, BoreStep, CircularPocketTool, CurvedEdge2d, CurvedProfile, EdgeTreatment,
    GeometryPolicy, Line8NutCut, OrthogonalCountersink, PlaneFrame, RevolutionAxis, RigidTransform,
    Solid, TransverseBoreTool, TreatmentEdge,
};
use std::{collections::BTreeMap, f64::consts::TAU};

const SHEET_MM: f64 = 5.0;
const NUT_THREAD_MM: f64 = 4.1;
const BORE_RADIUS_MM: f64 = 4.9;
const BORE_MOUTH_RADIUS_MM: f64 = 8.04;
const BORE_OVERTRAVEL_MM: f64 = 0.02;
const PLATE_PITCH_MM: f64 = 40.;
const PLATE_WIDTH_MM: f64 = 80.;
const PLATE_DEPTH_MM: f64 = 40.;
const HARDWARE_CLEARANCE_VOLUME_MM3: f64 = 0.001;
// The universal fastener's screw passage and head seat meet at the barrel axis.
const BARREL_PASSAGE_RADIUS_MM: f64 = 4.1;
const BARREL_SEAT_RADIUS_MM: f64 = 7.04;
const BARREL_PASSAGE_LENGTH_MM: f64 = 12.0;
const BARREL_PASSAGE_OVERTRAVEL_MM: f64 = 0.02;
// Assembly engineering clearance, independent of catalog and ISO fit classes.
const THREAD_RADIAL_CLEARANCE_MM: f64 = 0.08;
const THREAD_AXIAL_CLEARANCE_MM: f64 = 0.04;

type Point = [f64; 3];
fn add(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn scale(v: Point, factor: f64) -> Point {
    v.map(|x| x * factor)
}
fn dot(a: Point, b: Point) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn magnitude([x, y, z]: Point) -> f64 {
    x.hypot(y).hypot(z)
}
fn shift(frame: ReferenceFrame, distance: f64) -> ReferenceFrame {
    ReferenceFrame {
        origin_mm: add(frame.origin_mm, scale(frame.normal(), distance)),
        ..frame
    }
}
fn plane(frame: ReferenceFrame) -> PlaneFrame {
    PlaneFrame {
        origin: frame.origin_mm,
        u: frame.u,
        v: frame.v,
    }
}
fn placement(frame: ReferenceFrame) -> RigidTransform {
    let [u, v, n] = [frame.u, frame.v, frame.normal()];
    RigidTransform {
        rotation: std::array::from_fn(|i| [u[i], v[i], n[i]]),
        translation_mm: frame.origin_mm,
    }
}
fn polygon(points: &[[f64; 2]]) -> Vec<CurvedEdge2d> {
    points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .map(|(&a, &b)| CurvedEdge2d::Line { a, b })
        .collect()
}
fn rotated(frame: ReferenceFrame, angle: f64) -> ReferenceFrame {
    let (s, c) = angle.sin_cos();
    ReferenceFrame {
        u: add(scale(frame.u, c), scale(frame.v, s)),
        v: add(scale(frame.u, -s), scale(frame.v, c)),
        ..frame
    }
}

fn joining_plate_frame([a, b]: [ReferenceFrame; 2]) -> ReferenceFrame {
    let u = scale(sub(b.origin_mm, a.origin_mm), 1. / PLATE_PITCH_MM);
    let v = cross(a.normal(), u);
    ReferenceFrame {
        origin_mm: sub(
            sub(a.origin_mm, scale(u, PLATE_PITCH_MM / 2.)),
            scale(v, PLATE_DEPTH_MM / 2.),
        ),
        u,
        v,
    }
}

pub(crate) struct Inputs<'a> {
    pub document: &'a Document,
    pub parameters: &'a EvaluatedParameters,
    pub revision: u64,
    pub snapshots: &'a BTreeMap<EntityId, TopologySnapshot>,
    pub paths: &'a BTreeMap<EntityId, PathSegment>,
    pub solids: &'a BTreeMap<EntityId, Solid>,
    pub policy: GeometryPolicy,
}

#[derive(Clone, Copy)]
enum ReceivingHardware {
    Nut,
    /// A profiled screw or cross-bored barrel against members that may be
    /// cut by transverse access bores.
    Round,
    Screw,
    SeatedStock,
}
struct Builder<'a> {
    inputs: Inputs<'a>,
    feature: EntityId,
}
impl Builder<'_> {
    fn checked<T, E: std::fmt::Display>(&self, value: Result<T, E>) -> Result<T, ModelError> {
        value.map_err(|e| error(Some(self.feature), e))
    }
    fn source_checked<T, E: std::fmt::Debug>(&self, value: Result<T, E>) -> Result<T, ModelError> {
        value.map_err(|e| self.fail(format!("Retained source qualification: {e:?}")))
    }
    fn fail(&self, message: impl std::fmt::Display) -> ModelError {
        error(Some(self.feature), message)
    }
    fn solid(&self, producer: EntityId) -> Result<&Solid, ModelError> {
        self.inputs.solids.get(&producer).ok_or_else(|| {
            self.fail("A profile predecessor has no current native geometry. Repair the member reference.")
        })
    }
    fn member(&self, site: &ProfileAttachment) -> Result<model::ResolvedMember, ModelError> {
        let source = self.checked(model::read_member(
            self.inputs.document,
            site.member.producer,
        ))?;
        let path = self
            .inputs
            .paths
            .get(&site.member.producer)
            .ok_or_else(|| {
                self.fail(
                    "The retained profile member path is unavailable. Repair the member source.",
                )
            })?;
        self.checked(model::resolve_member(
            self.inputs.document,
            self.inputs.parameters,
            &source,
            *path,
        ))
    }
    fn cylinder(&self, frame: ReferenceFrame, r: f64, h: f64) -> Result<Solid, ModelError> {
        self.checked(Solid::cylinder(plane(frame), r, h, self.inputs.policy))
    }
    fn cuboid(&self, frame: ReferenceFrame, size: Point) -> Result<Solid, ModelError> {
        self.checked(Solid::box_solid(size, self.inputs.policy))?
            .transformed(placement(frame))
            .map_err(|e| self.fail(e))
    }
    fn subtract(&self, stock: &Solid, tool: &Solid, label: &str) -> Result<Solid, ModelError> {
        self.checked(stock.boolean(tool, BooleanOperation::Subtract, self.inputs.policy))?
            .ok_or_else(|| self.fail(label))
    }
    fn check_clear(&self, a: &Solid, b: &Solid, message: &str) -> Result<(), ModelError> {
        if let Some(overlap) =
            self.checked(a.boolean(b, BooleanOperation::Intersect, self.inputs.policy))?
            && self.checked(overlap.volume_mm3())? > HARDWARE_CLEARANCE_VOLUME_MM3
        {
            return Err(self.fail(message));
        }
        Ok(())
    }

    fn countersink(&self, frame: ReferenceFrame) -> Result<Solid, ModelError> {
        let top = shift(frame, -SHEET_MM - BORE_OVERTRAVEL_MM);
        let radius = BORE_MOUTH_RADIUS_MM;
        let taper_end = radius - BORE_RADIUS_MM;
        let depth = SHEET_MM + 2. * BORE_OVERTRAVEL_MM;
        let outline = [
            [0.0, 0.0],
            [radius, 0.0],
            [BORE_RADIUS_MM, taper_end],
            [BORE_RADIUS_MM, depth],
            [0.0, depth],
        ];
        self.checked(Solid::revolve_curved_profile(
            &CurvedProfile {
                frame: PlaneFrame {
                    origin: top.origin_mm,
                    u: top.u,
                    v: top.normal(),
                },
                outer: polygon(&outline),
                holes: vec![],
            },
            RevolutionAxis {
                origin_mm: top.origin_mm,
                direction: top.normal(),
            },
            TAU,
            self.inputs.policy,
        ))
    }

    fn plate(
        &self,
        kind: ProfileConnectionKind,
        mount: ProfileConnectionFrames,
    ) -> Result<Solid, ModelError> {
        let [a, b] = mount.sites;
        let mut plate = if kind == ProfileConnectionKind::AngleBracket {
            let corner = sub(a.origin_mm, scale(b.normal(), 20.0));
            let leg = |face: ReferenceFrame, other: ReferenceFrame| {
                let u = other.normal();
                let v = cross(face.normal(), u);
                self.cuboid(
                    ReferenceFrame {
                        origin_mm: sub(corner, scale(v, 20.0)),
                        u,
                        v,
                    },
                    [40.0, 40.0, SHEET_MM],
                )
            };
            self.checked(leg(a, b)?.boolean(
                &leg(b, a)?,
                BooleanOperation::Union,
                self.inputs.policy,
            ))?
            .ok_or_else(|| self.fail("Angle bracket legs do not meet at their bend"))?
        } else {
            self.cuboid(
                joining_plate_frame(mount.sites),
                [PLATE_WIDTH_MM, PLATE_DEPTH_MM, SHEET_MM],
            )?
        };
        for frame in mount.fasteners {
            plate = self.subtract(
                &plate,
                &self.countersink(frame)?,
                "Countersink removed the joining part",
            )?;
        }
        Ok(plate)
    }

    fn rust_angle_bracket(&self, mount: ProfileConnectionFrames) -> Result<Solid, ModelError> {
        use multicad_kernel::rust_kernel::{
            orthogonal::{self, Box3, Operation},
            orthogonal_bore_seats::Seat,
        };
        let [a, b] = mount.sites;
        let corner = sub(a.origin_mm, scale(b.normal(), 20.));
        let u = cross(b.normal(), a.normal());
        let frame = ReferenceFrame {
            origin_mm: sub(corner, scale(u, 20.)),
            u,
            v: b.normal(),
        };
        self.checked(frame.validate())?;
        let geometry = self
            .source_checked(orthogonal::boolean(
                &[
                    Box3 {
                        min: [0.; 3],
                        max: [40., 40., SHEET_MM],
                    },
                    Box3 {
                        min: [0.; 3],
                        max: [40., SHEET_MM, 40.],
                    },
                ],
                &[],
                Operation::Union,
                self.inputs.policy.linear_tolerance_mm,
            ))?
            .ok_or_else(|| self.fail("The original angle stock has no physical material"))?;
        let stock = self.checked(Solid::rust_orthogonal_stock_named(
            self.feature,
            geometry,
            self.inputs.policy,
        ))?;
        let travel = SHEET_MM + 2. * BORE_OVERTRAVEL_MM;
        let tools = [
            (
                [20., 20., SHEET_MM + BORE_OVERTRAVEL_MM],
                [0., 0., -travel],
                [1., 0., 0.],
                [0., -1., 0.],
            ),
            (
                [20., SHEET_MM + BORE_OVERTRAVEL_MM, 20.],
                [0., -travel, 0.],
                [1., 0., 0.],
                [0., 0., 1.],
            ),
        ]
        .into_iter()
        .enumerate()
        .map(|(i, (origin, displacement_mm, u, v))| CircularPocketTool {
            frame: PlaneFrame { origin, u, v },
            center_uv: [0.; 2],
            radius_mm: BORE_RADIUS_MM,
            displacement_mm,
            origin: Origin {
                source: self.feature,
                role: format!("angle_bracket.bore.{i}"),
            },
        })
        .collect::<Vec<_>>();
        let drilled = self.checked(stock.rust_drill_circular_machining(
            &tools,
            &[true, true],
            self.inputs.policy,
        ))?;
        let radius = BORE_MOUTH_RADIUS_MM - BORE_OVERTRAVEL_MM;
        let seats = [0, 1].map(|bore| OrthogonalCountersink {
            geometry: Seat {
                bore,
                end: false,
                radius_mm: radius,
                depth_mm: radius - BORE_RADIUS_MM,
            },
            origin: Origin {
                source: self.feature,
                role: format!("angle_bracket.countersink.{bore}"),
            },
        });
        let local =
            self.checked(drilled.rust_countersink_orthogonal(&seats, self.inputs.policy, None))?;
        let bracket = self.checked(local.transformed(placement(frame)))?;
        self.rust_angle_bore_correspondence(&bracket, mount)?;
        Ok(bracket)
    }
    fn rust_angle_bore_correspondence(
        &self,
        bracket: &Solid,
        mount: ProfileConnectionFrames,
    ) -> Result<(), ModelError> {
        use multicad_kernel::rust_kernel::axial_pipe_thread::authored::StockClipPlacement;
        let (source, frame) = bracket
            .rust_orthogonal_countersink_source()
            .ok_or_else(|| self.fail("The angle lost its original finite countersinks"))?;
        let axes = [frame.u, frame.v, frame.normal()];
        let matrix = std::array::from_fn(|i| axes.map(|v| v[i]));
        let placed = self.source_checked(StockClipPlacement::new(matrix, frame.origin))?;
        let vector = self.source_checked(StockClipPlacement::new(matrix, [0.; 3]))?;
        for (i, axis) in mount.fasteners.into_iter().enumerate() {
            let bore = source.bores()[i];
            let position =
                self.source_checked(placed.point_bounds(bore.centers_mm[1].map(|v| [v, v])))?;
            let normal = source.boundary().faces[bore.faces[0]].normal.map(|v| -v);
            let direction = self.source_checked(vector.point_bounds(normal.map(|v| [v, v])))?;
            let target_axes = [axis.u, axis.v, axis.normal()];
            let target = self.source_checked(StockClipPlacement::new(
                std::array::from_fn(|j| target_axes.map(|a| a[j])),
                [0.; 3],
            ))?;
            let difference = |bounds: [[f64; 2]; 3], point: Point| {
                (0..3).fold(0_f64, |sum, j| {
                    (sum + (bounds[j][0] - point[j])
                        .abs()
                        .next_up()
                        .max((bounds[j][1] - point[j]).abs().next_up()))
                    .next_up()
                })
            };
            let residual = (difference(position, axis.origin_mm)
                + 80.
                    * (difference(direction, axis.normal())
                        + placed.metric_error()
                        + target.metric_error())
                    .next_up())
            .next_up();
            if !residual.is_finite() || residual > self.inputs.policy.linear_tolerance_mm / 16. {
                return Err(self.fail("The actual angle bore axes do not match both retained mounting frames within the source construction budget"));
            }
        }
        Ok(())
    }
    fn rust_joining_plate(&self, mount: ProfileConnectionFrames) -> Result<Solid, ModelError> {
        let frame = joining_plate_frame(mount.sites);
        self.checked(frame.validate())?;
        let policy = self.inputs.policy;
        let normal = frame.normal();
        let outer = [
            [0., 0.],
            [PLATE_WIDTH_MM, 0.],
            [PLATE_WIDTH_MM, PLATE_DEPTH_MM],
            [0., PLATE_DEPTH_MM],
        ];
        // A Line 8 groove datum is an outward support plane of the full
        // member section. A corner check therefore proves clearance of the
        // plate enclosure, including for rigidly placed profile occurrences.
        let tolerance = policy.linear_tolerance_mm;
        let penetration =
            tolerance.min(HARDWARE_CLEARANCE_VOLUME_MM3 / (2. * PLATE_WIDTH_MM * PLATE_DEPTH_MM));
        for site in mount.sites {
            for point in outer.map(|uv| plane(frame).point(uv)) {
                for point in [point, add(point, scale(normal, SHEET_MM))] {
                    if dot(sub(point, site.origin_mm), site.normal()) < -penetration {
                        return Err(self.fail(
                            "Joining plate cannot be certified clear of the retained profile surfaces",
                        ));
                    }
                }
            }
        }

        let mut tools = Vec::with_capacity(2);
        let mut countersinks = Vec::with_capacity(2);
        for (i, axis) in mount.fasteners.into_iter().enumerate() {
            let inward = axis.normal();
            if dot(normal, inward) >= 0. || magnitude(cross(normal, inward)) > 1e-12 {
                return Err(self.fail(
                    "Rust joining plate countersinks require axes perpendicular to the plate",
                ));
            }
            let top = shift(axis, -SHEET_MM - BORE_OVERTRAVEL_MM);
            let entry = dot(sub(top.origin_mm, frame.origin_mm), normal) - SHEET_MM;
            let chamfer = BORE_MOUTH_RADIUS_MM - BORE_RADIUS_MM - entry;
            if !(0. ..2. * BORE_OVERTRAVEL_MM).contains(&entry)
                || chamfer <= 0.
                || chamfer >= SHEET_MM
            {
                return Err(self.fail(
                    "The joining plate countersink tool must cut through both plate surfaces",
                ));
            }
            tools.push(CircularPocketTool {
                frame: plane(top),
                center_uv: [0., 0.],
                radius_mm: BORE_RADIUS_MM,
                displacement_mm: scale(inward, SHEET_MM + 2. * BORE_OVERTRAVEL_MM),
                origin: Origin {
                    source: self.feature,
                    role: format!("joining_plate.bore.{i}"),
                },
            });
            countersinks.push((add(top.origin_mm, scale(inward, entry)), chamfer));
        }
        let profile = CurvedProfile {
            frame: plane(frame),
            outer: polygon(&outer),
            holes: vec![],
        };
        let origins: Vec<_> = (0..profile.outer.len())
            .map(|i| Origin {
                source: self.feature,
                role: format!("joining_plate.outline.{i}"),
            })
            .collect();
        let blank =
            crate::rust_backend::stock_profile(self.feature, &profile, &origins, SHEET_MM, policy)?;
        let mut plate =
            self.checked(blank.rust_drill_circular_machining(&tools, &[true; 2], policy))?;
        for (i, (mouth, distance_mm)) in countersinks.into_iter().enumerate() {
            let snapshot =
                self.checked(plate.topology_snapshot(self.feature, self.inputs.revision))?;
            let edge = snapshot
                .elements
                .iter()
                .find(|element| {
                    matches!(
                        element.geometry,
                        Some(AnalyticGeometry::CircularArc { frame, radius_mm, .. })
                            if magnitude(sub(frame.origin_mm, mouth)) <= tolerance
                                && (radius_mm - BORE_RADIUS_MM).abs() <= tolerance
                    )
                })
                .ok_or_else(|| self.fail("Joining plate lost the intended countersink mouth"))?;
            plate = self.checked(plate.treat_edges_named(
                &[TreatmentEdge {
                    key: edge.key,
                    origin: Origin {
                        source: self.feature,
                        role: format!("joining_plate.countersink.{i}"),
                    },
                }],
                EdgeTreatment::Chamfer { distance_mm },
                self.feature,
                policy,
            ))?;
        }
        Ok(plate)
    }

    fn stepped_bores(
        &self,
        stock: &Solid,
        frames: [ReferenceFrame; 2],
    ) -> Result<Solid, ModelError> {
        let mut result = self.checked(stock.transformed(RigidTransform::default()))?;
        let initial_volume = self.checked(result.volume_mm3())?;
        for frame in frames {
            // Start a little outside the actual face while retaining the
            // exact face-relative endpoints of the manufacturer operation.
            let outside = 0.02;
            let start = shift(frame, -outside);
            let counterbore = self.cylinder(
                start,
                PROFILE_ACCESS_RADIUS_MM,
                PROFILE_BARREL_DEPTH_MM + outside,
            )?;
            let pilot = self.cylinder(
                start,
                PROFILE_PILOT_RADIUS_MM,
                PROFILE_PILOT_DEPTH_MM + outside,
            )?;
            let before = self.checked(result.volume_mm3())?;
            result = self.subtract(
                &result,
                &counterbore,
                "Access counterbore removed the member",
            )?;
            if before - self.checked(result.volume_mm3())? < 1.0 {
                return Err(self.fail("The access counterbore missed the profile material. Repair its end or groove orientation."));
            }
            result = self.subtract(&result, &pilot, "Access drill removed the member")?;
        }
        if self.checked(result.solid_count())? != 1
            || self.checked(result.volume_mm3())? >= initial_volume
        {
            return Err(self.fail(
                "Stepped access bores must remove material and preserve one receiving solid",
            ));
        }
        Ok(result)
    }

    /// One coaxial cutter for both entries: each a counterbore over a
    /// deeper pilot, started just outside its face.
    fn rust_stepped_bores(
        &self,
        stock: &Solid,
        frames: [ReferenceFrame; 2],
    ) -> Result<Solid, ModelError> {
        let outside = 0.02;
        let axis = frames[0].normal();
        let origin = frames[0].origin_mm;
        let mut steps = Vec::new();
        for frame in frames {
            let normal = frame.normal();
            let along = dot(normal, axis);
            let offset = sub(frame.origin_mm, origin);
            let station = dot(offset, axis);
            let across = magnitude(sub(offset, scale(axis, station)));
            if (along.abs() - 1.).abs() > 1e-12 || across > self.inputs.policy.linear_tolerance_mm {
                return Err(
                    self.fail("Both access bores must share one axis through the receiving member")
                );
            }
            for (radius_mm, depth) in [
                (PROFILE_ACCESS_RADIUS_MM, PROFILE_BARREL_DEPTH_MM),
                (PROFILE_PILOT_RADIUS_MM, PROFILE_PILOT_DEPTH_MM),
            ] {
                let [start_mm, end_mm] = if along > 0. {
                    [station - outside, station + depth]
                } else {
                    [station - depth, station + outside]
                };
                steps.push(BoreStep {
                    start_mm,
                    end_mm,
                    radius_mm,
                });
            }
        }
        let tool = TransverseBoreTool {
            origin_mm: origin,
            axis,
            steps,
            wall: Origin {
                source: self.feature,
                role: "access_bore.wall".into(),
            },
            shoulder: Origin {
                source: self.feature,
                role: "access_bore.shoulder".into(),
            },
        };
        let before = self.checked(stock.volume_mm3())?;
        let result = self.checked(stock.rust_transverse_bore(&tool, self.inputs.policy))?;
        if self.checked(result.solid_count())? != 1 || self.checked(result.volume_mm3())? >= before
        {
            return Err(self.fail(
                "Stepped access bores must remove material and preserve one receiving solid",
            ));
        }
        Ok(result)
    }

    /// The universal fastener barrel with its screw passage and head seat
    /// crossing at the barrel axis.
    fn rust_barrel(
        &self,
        access: ReferenceFrame,
        screw: ReferenceFrame,
    ) -> Result<Solid, ModelError> {
        let barrel = self.cylinder(
            shift(access, 0.01),
            PROFILE_ACCESS_RADIUS_MM - 0.02,
            PROFILE_BARREL_DEPTH_MM - 0.02,
        )?;
        let station = shift(screw, -PROFILE_ACCESS_FROM_END_MM);
        let tool = TransverseBoreTool {
            origin_mm: station.origin_mm,
            axis: station.normal(),
            steps: vec![
                BoreStep {
                    start_mm: -BARREL_PASSAGE_OVERTRAVEL_MM,
                    end_mm: BARREL_PASSAGE_LENGTH_MM - BARREL_PASSAGE_OVERTRAVEL_MM,
                    radius_mm: BARREL_PASSAGE_RADIUS_MM,
                },
                BoreStep {
                    start_mm: -BARREL_PASSAGE_LENGTH_MM,
                    end_mm: 0.,
                    radius_mm: BARREL_SEAT_RADIUS_MM,
                },
            ],
            wall: Origin {
                source: self.feature,
                role: "universal_fastener.passage".into(),
            },
            shoulder: Origin {
                source: self.feature,
                role: "universal_fastener.seat".into(),
            },
        };
        self.checked(barrel.rust_transverse_bore(&tool, self.inputs.policy))
    }

    fn barrel(&self, access: ReferenceFrame, screw: ReferenceFrame) -> Result<Solid, ModelError> {
        let barrel = self.cylinder(
            shift(access, 0.01),
            PROFILE_ACCESS_RADIUS_MM - 0.02,
            PROFILE_BARREL_DEPTH_MM - 0.02,
        )?;
        let station = shift(screw, -PROFILE_ACCESS_FROM_END_MM);
        let passage = self.cylinder(
            shift(station, -BARREL_PASSAGE_OVERTRAVEL_MM),
            BARREL_PASSAGE_RADIUS_MM,
            BARREL_PASSAGE_LENGTH_MM,
        )?;
        let drive = self.cylinder(
            shift(station, -BARREL_PASSAGE_LENGTH_MM),
            BARREL_SEAT_RADIUS_MM,
            BARREL_PASSAGE_LENGTH_MM,
        )?;
        let body = self.subtract(
            &barrel,
            &passage,
            "Universal fastener lost its screw passage",
        )?;
        self.subtract(&body, &drive, "Universal fastener lost its head seat")
    }

    fn screw(
        &self,
        kind: ProfileConnectionKind,
        frame: ReferenceFrame,
        backend: GeometryBackend,
    ) -> Result<(Solid, ResolvedThread, ReferenceFrame), ModelError> {
        let universal = kind == ProfileConnectionKind::UniversalEnd;
        let height = 4.4;
        let (head, depth, underhead, runout, neck) = if universal {
            (
                ScrewHead::Button {
                    profile: Box::new(CrownedScrewHeadParameters {
                        diameter: length(14.0),
                        height: length(height),
                        crown_radius: length(7.7),
                        top_land_diameter: length(7.7),
                        neck_radius: length(0.7),
                    }),
                },
                2.6,
                30.0,
                1.95,
                -PROFILE_ACCESS_FROM_END_MM,
            )
        } else {
            (
                ScrewHead::Countersunk {
                    profile: Box::new(CountersunkScrewHeadParameters {
                        diameter: length(16.0),
                        height: length(height),
                        included_angle: model::angle(90.0),
                        edge_chamfer: length(0.05),
                        neck_radius: length(0.4),
                    }),
                },
                3.2,
                14.0 - height,
                1.25,
                -SHEET_MM + height,
            )
        };
        let authored = if universal {
            ScrewLengthInput::Underhead {
                value: length(underhead),
            }
        } else {
            ScrewLengthInput::Overall {
                value: length(14.0),
            }
        };
        let part = MechanicalPartParameters::MachineScrew {
            underhead_length: None,
            length_input: Some(authored),
            shank_diameter: length(8.0),
            head,
            drive: ScrewDrive::HexSocket {
                across_flats: length(5.08),
                depth: length(depth),
            },
            thread: Box::new(ThreadParameters {
                length: length(underhead - runout),
                radial_allowance: length(THREAD_RADIAL_CLEARANCE_MM),
                axial_allowance: length(THREAD_AXIAL_CLEARANCE_MM),
                ..Default::default()
            }),
        };
        let resolved = self.checked(part.resolve(self.inputs.parameters))?;
        let ResolvedMechanicalPart::MachineScrew { thread, .. } = &resolved else {
            unreachable!()
        };
        let thread = thread.as_ref().clone();
        let solid = match backend {
            GeometryBackend::Rust => crate::mechanical_parts::build_rust_resolved(
                self.feature,
                resolved,
                self.inputs.policy,
            ),
            GeometryBackend::OpenCascade => {
                crate::mechanical_parts::build_resolved(self.feature, resolved, self.inputs.policy)
            }
        }?;
        let placed = self.checked(solid.transformed(placement(shift(frame, neck))))?;
        Ok((placed, thread, shift(frame, neck + runout)))
    }

    fn nut(
        &self,
        kind: ProfileConnectionKind,
        variant: ProfileConnectionNut,
        frame: ReferenceFrame,
        backend: GeometryBackend,
    ) -> Result<Solid, ModelError> {
        let size = variant.dimensions();
        let half_width = size.width_mm / 2.;
        let seat = shift(frame, ITEM_LINE_8_GROOVE_LIP_DEPTH_MM);
        let section = [
            [-half_width, 0.0],
            [half_width, 0.0],
            [6.5, 3.0],
            [5.4, 5.0],
            [3.3, size.depth_mm],
            [-3.3, size.depth_mm],
            [-5.4, 5.0],
            [-6.5, 3.0],
        ];
        let blank = if backend == GeometryBackend::Rust {
            // Author the original straight local section directly. Projecting
            // a rounded world displacement back into an arbitrary profile
            // frame would create a different, slightly sheared receiving stock.
            let owner = EntityId::new();
            let origins = (0..section.len())
                .map(|i| Origin {
                    source: owner,
                    role: format!("line8.stock.side.{i}"),
                })
                .collect::<Vec<_>>();
            self.checked(Solid::rust_extrude_linear_named(
                &CurvedProfile {
                    frame: PlaneFrame {
                        origin: [-size.bore_from_end_mm, 0., 0.],
                        u: [0., 1., 0.],
                        v: [0., 0., 1.],
                    },
                    outer: polygon(&section),
                    holes: vec![],
                },
                [size.length_mm, 0., 0.],
                0.,
                owner,
                &origins,
                self.inputs.policy,
            ))?
        } else {
            self.checked(Solid::extrude_curved_profile(
                &CurvedProfile {
                    frame: PlaneFrame {
                        origin: sub(seat.origin_mm, scale(seat.u, size.bore_from_end_mm)),
                        u: seat.v,
                        v: seat.normal(),
                    },
                    outer: polygon(&section),
                    holes: vec![],
                },
                size.length_mm,
                self.inputs.policy,
            ))?
        };
        let (screw, external, external_start) = self.screw(kind, frame, backend)?;
        let params = ThreadParameters {
            length: length(NUT_THREAD_MM),
            radial_allowance: length(THREAD_RADIAL_CLEARANCE_MM),
            axial_allowance: length(THREAD_AXIAL_CLEARANCE_MM),
            minimum_wall: length(0.5),
            ..Default::default()
        };
        let mut internal =
            self.checked(params.resolve(self.inputs.parameters, ThreadMode::CutInternal))?;
        let phase = self.checked(model::straight_thread_phase(
            &external,
            external_start,
            &internal,
            seat,
            self.inputs.policy.linear_tolerance_mm,
        ))?;
        let start = phase.axial_station_mm;
        let engagement = external.length_mm.min(start + internal.length_mm) - start.max(0.0);
        if engagement < 2.0 * external.pitch_mm {
            return Err(self.fail(
                "The M8 hardware needs at least two complete modeled turns inside the T-slot nut",
            ));
        }
        let inner = rotated(seat, phase.roll_rad);
        internal.frame = multicad_core::sketch::Frame {
            origin: inner.origin_mm,
            u: inner.u,
            v: inner.v,
            n: inner.normal(),
        };
        if backend == GeometryBackend::Rust {
            let (sin, cos) = phase.roll_rad.sin_cos();
            internal.frame = multicad_core::sketch::Frame {
                origin: [0.; 3],
                u: [cos, sin, 0.],
                v: [-sin, cos, 0.],
                n: [0., 0., 1.],
            };
            let (profile, root, local_frame) =
                crate::threads::rust_profile(self.feature, &internal)?;
            if root.is_some() {
                return Err(
                    self.fail("Rust Line8 nuts require the original flat internal thread profile")
                );
            }
            let local = self.checked(blank.rust_cut_line8_nut(
                Line8NutCut {
                    profile,
                    frame: local_frame,
                    bore_from_start_mm: size.bore_from_end_mm,
                    relief: multicad_kernel::rust_kernel::line8_nut::Relief {
                        start_mm: NUT_THREAD_MM - 0.001,
                        radius_mm: 4.10,
                        tool_length_mm: size.depth_mm - NUT_THREAD_MM + 0.02,
                    },
                    minimum_wall_mm: internal.minimum_wall_mm,
                },
                self.feature,
                self.inputs.policy,
                None,
            ))?;
            let nut = self.checked(local.transformed(placement(seat)))?;
            self.checked(screw.rust_line8_nut_mating(&nut, None))?;
            return Ok(nut);
        }
        let tapped = crate::threads::build_resolved(
            self.feature,
            internal,
            ThreadMode::CutInternal,
            Some(&blank),
            self.inputs.policy,
        )?;
        let relief = self.cylinder(
            shift(seat, NUT_THREAD_MM - 0.001),
            4.10,
            size.depth_mm - NUT_THREAD_MM + 0.02,
        )?;
        let nut = self.subtract(&tapped, &relief, "T-slot nut crown relief removed the part")?;
        self.check_clear(
            &screw,
            &nut,
            "The supported M8 screw and T-slot nut threads interfere",
        )?;
        Ok(nut)
    }
    fn rust_hardware_receivers(
        &self,
        hardware: &Solid,
        kind: ReceivingHardware,
        request: &model::ProfileConnectionRequest,
        access_bores: Option<EntityId>,
    ) -> Result<(), ModelError> {
        for (site, producer) in [
            (
                &request.first,
                access_bores.unwrap_or(request.first.member.producer),
            ),
            (&request.second, request.second.member.producer),
        ] {
            let pose = self.checked(model::resolve_occurrence_pose(
                self.inputs.document,
                &site.occurrence,
                site.member.body,
                self.inputs.parameters,
                self.inputs.revision,
                self.inputs.snapshots,
            ))?;
            let occurrence = RigidTransform {
                rotation: pose.rotation,
                translation_mm: pose.translation_mm,
            };
            match kind {
                ReceivingHardware::Screw => {
                    self.checked(hardware.rust_profiled_screw_receiver_clearance(
                        self.solid(producer)?,
                        occurrence,
                        HARDWARE_CLEARANCE_VOLUME_MM3,
                        None,
                    ))?;
                }
                ReceivingHardware::Nut => {
                    self.checked(hardware.rust_line8_receiver_clearance(
                        self.solid(producer)?,
                        occurrence,
                        HARDWARE_CLEARANCE_VOLUME_MM3,
                        None,
                    ))?;
                }
                ReceivingHardware::Round => {
                    self.checked(hardware.rust_round_hardware_receiver_clearance(
                        self.solid(producer)?,
                        occurrence,
                        HARDWARE_CLEARANCE_VOLUME_MM3,
                        None,
                    ))?;
                }
                ReceivingHardware::SeatedStock => {
                    self.checked(hardware.rust_seated_stock_receiver_clearance(
                        self.solid(producer)?,
                        occurrence,
                        HARDWARE_CLEARANCE_VOLUME_MM3,
                        None,
                    ))?;
                }
            }
        }
        if self.checked(hardware.solid_count())? != 1 {
            return Err(self.fail("Each retained connector component must be one physical solid"));
        }
        Ok(())
    }
}

fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(crate) fn build(
    inputs: Inputs<'_>,
    feature: EntityId,
    backend: GeometryBackend,
) -> Result<Solid, ModelError> {
    let builder = Builder { inputs, feature };
    let ctx = &builder.inputs;
    let output = builder.checked(model::read_profile_connection_output(ctx.document, feature))?;
    let request = builder.checked(model::read_profile_connection_request(
        ctx.document,
        output.connection,
    ))?;
    let first = builder.member(&request.first)?;
    let second = builder.member(&request.second)?;
    let frames = builder.checked(model::resolve_profile_connection_frames(
        ctx.document,
        &request,
        [&first, &second],
        ctx.parameters,
        ctx.revision,
        ctx.snapshots,
    ))?;
    if backend == GeometryBackend::Rust {
        return match output.geometry {
            ProfileOutputKind::AccessBores { target } => {
                let access = frames.access.ok_or_else(|| {
                    builder.fail("The end connection lost its access bore datums")
                })?;
                builder.rust_stepped_bores(builder.solid(target.producer)?, access)
            }
            ProfileOutputKind::Part { index } => match request.kind.roles()[index] {
                ProfilePartRole::Connector(_) => match request.kind {
                    ProfileConnectionKind::JoiningPlate => builder.rust_joining_plate(frames),
                    ProfileConnectionKind::AngleBracket => {
                        let bracket = builder.rust_angle_bracket(frames)?;
                        builder.rust_hardware_receivers(
                            &bracket,
                            ReceivingHardware::SeatedStock,
                            &request,
                            output.access_bores,
                        )?;
                        Ok(bracket)
                    }
                    ProfileConnectionKind::UniversalEnd => {
                        let ProfilePartRole::Connector(i) = request.kind.roles()[index] else {
                            unreachable!("a connector role")
                        };
                        let access = frames.access.ok_or_else(|| {
                            builder.fail("The end connection lost its access bore datums")
                        })?;
                        let barrel = builder.rust_barrel(access[i], frames.fasteners[i])?;
                        builder.rust_hardware_receivers(
                            &barrel,
                            ReceivingHardware::Round,
                            &request,
                            output.access_bores,
                        )?;
                        Ok(barrel)
                    }
                },
                ProfilePartRole::Screw(index) => {
                    let screw = builder
                        .screw(request.kind, frames.fasteners[index], backend)?
                        .0;
                    let kind = if output.access_bores.is_some() {
                        ReceivingHardware::Round
                    } else {
                        ReceivingHardware::Screw
                    };
                    builder.rust_hardware_receivers(&screw, kind, &request, output.access_bores)?;
                    Ok(screw)
                }
                ProfilePartRole::Nut(index) => {
                    let nut =
                        builder.nut(request.kind, request.nut, frames.fasteners[index], backend)?;
                    builder.rust_hardware_receivers(
                        &nut,
                        ReceivingHardware::Nut,
                        &request,
                        output.access_bores,
                    )?;
                    Ok(nut)
                }
            },
        };
    }
    let ProfileOutputKind::Part { index } = output.geometry else {
        let ProfileOutputKind::AccessBores { target } = output.geometry else {
            unreachable!()
        };
        let access = frames
            .access
            .ok_or_else(|| builder.fail("The end connection lost its access bore datums"))?;
        return builder.stepped_bores(builder.solid(target.producer)?, access);
    };
    let solid = match request.kind.roles()[index] {
        ProfilePartRole::Connector(index) => {
            if let Some(access) = frames.access {
                builder.barrel(access[index], frames.fasteners[index])?
            } else {
                builder.plate(request.kind, frames)?
            }
        }
        ProfilePartRole::Screw(index) => {
            builder
                .screw(request.kind, frames.fasteners[index], backend)?
                .0
        }
        ProfilePartRole::Nut(index) => {
            builder.nut(request.kind, request.nut, frames.fasteners[index], backend)?
        }
    };
    for (site, producer) in [
        (
            &request.first,
            output.access_bores.unwrap_or(request.first.member.producer),
        ),
        (&request.second, request.second.member.producer),
    ] {
        let pose = builder.checked(model::resolve_occurrence_pose(
            ctx.document,
            &site.occurrence,
            site.member.body,
            ctx.parameters,
            ctx.revision,
            ctx.snapshots,
        ))?;
        let receiving = builder.checked(builder.solid(producer)?.transformed(RigidTransform {
            rotation: pose.rotation,
            translation_mm: pose.translation_mm,
        }))?;
        builder.check_clear(&solid, &receiving, "Connector hardware interferes with the retained profile section. Check both stations, groove sides and the end orientation.")?;
    }
    if builder.checked(solid.solid_count())? != 1 {
        return Err(builder.fail("Each retained connector component must be one physical solid"));
    }
    Ok(solid)
}
