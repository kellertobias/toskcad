//! Actual profiled heads with an explicitly cosmetic smooth shaft.
use super::*;
use crate::rust_prism_topology::position;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, TopologyElement, TopologyKey, TopologyKind,
};
use multicad_rust_kernel::{cosmetic_profiled_screw::Envelope, profiled_screw::Face as Role};
use std::collections::BTreeSet;
use std::f64::consts::{PI, TAU};
type Point = [f64; 3];
pub(crate) struct CosmeticProfiledBody {
    pub(super) feature: EntityId,
    geometry: Envelope,
    properties: Vec<(Point, f64)>,
}
enum Edge {
    Circle { radius: f64, z: f64 },
    Mouth(usize),
    Floor(usize),
    Post(usize),
}
struct Boundary {
    edge: Edge,
    faces: [usize; 2],
    name: String,
}
impl Solid {
    /// Preserve a qualified physical head/recess, replacing only its thread
    /// interval by a distinct smooth cylindrical envelope. No manufacturing
    /// or STEP admission follows from this conversion.
    pub fn cosmetic_profiled_screw_envelope(&self) -> Result<Self, GeometryError> {
        let unsupported = || {
            error(
                "cosmetic profiled screw",
                GeometryErrorKind::Unsupported,
                "A cosmetic curved-drive envelope requires an authored physical profiled screw",
            )
        };
        let (source, frame) = self
            .rust_profiled_screw_geometry()
            .ok_or_else(unsupported)?;
        if source.recess().hexalobular_geometry().is_none() {
            return Err(unsupported());
        }
        if let multicad_rust_kernel::screw_recess::Recess::ConicalHexalobular(recess) =
            source.recess()
            && recess.finite_binding().is_some()
        {
            return Err(unsupported());
        }
        let SolidBody::Helical(original, _) = &self.body else {
            return Err(unsupported());
        };
        let geometry = Envelope::new(source.clone());
        let heads = helical_body::profiled_screw::torx::head_face_properties(source, self.policy)?;
        let old = source.faces();
        let first = 2 + 4 * source.thread().spec().starts as usize;
        let [a, b] = geometry.shaft_interval_mm();
        let r = geometry.shaft_radius_mm();
        let properties = geometry
            .faces()
            .iter()
            .map(|role| match role {
                Role::ThreadShoulder => (
                    [0., 0., a],
                    PI * (source.profile().end()[0] - r) * (source.profile().end()[0] + r),
                ),
                Role::ThreadTip => ([0., 0., b], PI * r * r),
                Role::ThreadWall(0) => ([0., 0., a.midpoint(b)], TAU * r * (b - a)),
                role => {
                    heads[old
                        .iter()
                        .position(|f| f == role)
                        .expect("retained head role")
                        - first]
                }
            })
            .collect();
        let body = Self::cosmetic_profiled_body(
            CosmeticProfiledBody {
                feature: original.feature,
                geometry,
                properties,
            },
            frame,
            self.policy,
        )?;
        Ok(body)
    }
    fn cosmetic_profiled_body(
        body: CosmeticProfiledBody,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), policy)?;
        let bounds = body.bounds(frame);
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(Self {
            body: SolidBody::CosmeticProfiled(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn cosmetic_profiled_body_ref(&self) -> Option<(&CosmeticProfiledBody, PlaneFrame)> {
        match &self.body {
            SolidBody::CosmeticProfiled(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
}
impl CosmeticProfiledBody {
    fn id(&self, role: Role) -> usize {
        self.geometry
            .faces()
            .iter()
            .position(|r| *r == role)
            .expect("cosmetic physical face")
    }
    pub(super) fn volume_mm3(&self) -> f64 {
        self.geometry.volume_mm3()
    }
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Bounds {
        let axes = [frame.u, frame.v, frame.normal()];
        let support = |i, sign| self.geometry.support(axes.map(|a| sign * a[i]));
        Bounds {
            min: std::array::from_fn(|i| frame.origin[i] - support(i, -1.)),
            max: std::array::from_fn(|i| frame.origin[i] + support(i, 1.)),
        }
    }
    fn plane(&self, role: Role, frame: PlaneFrame) -> Option<ReferenceFrame> {
        let [a, b] = self.geometry.shaft_interval_mm();
        match role {
            Role::ThreadShoulder | Role::ThreadTip => Some(ReferenceFrame {
                origin_mm: position(frame, [0., 0., if role == Role::ThreadTip { b } else { a }]),
                u: frame.u,
                v: frame.v,
            }),
            _ => helical_body::profiled_screw::plane(self.geometry.head_source(), role, frame),
        }
    }
    pub(super) fn faces(&self, frame: PlaneFrame) -> Vec<Face> {
        self.geometry
            .faces()
            .iter()
            .zip(&self.properties)
            .enumerate()
            .map(|(i, (role, (p, area)))| Face {
                index: i as u32,
                plane: self.plane(*role, frame).map(|f| (f.origin_mm, f.normal())),
                center_mm: position(frame, *p),
                area_mm2: *area,
            })
            .collect()
    }
    fn boundaries(&self) -> Vec<Boundary> {
        let source = self.geometry.head_source();
        let [a, b] = self.geometry.shaft_interval_mm();
        let r = self.geometry.shaft_radius_mm();
        let profile: Vec<_> = source.profile().segments().collect();
        let last = Role::ProfileWall(profile.len() as u32 - 1);
        let mut result = vec![Boundary {
            edge: Edge::Circle { radius: r, z: b },
            name: "cosmetic.shaft/tip".into(),
            faces: [self.id(Role::ThreadTip), self.id(Role::ThreadWall(0))],
        }];
        let shoulder = self.geometry.faces().contains(&Role::ThreadShoulder);
        result.push(Boundary {
            edge: Edge::Circle { radius: r, z: a },
            name: "cosmetic.shaft/start".into(),
            faces: [
                self.id(Role::ThreadWall(0)),
                self.id(if shoulder { Role::ThreadShoulder } else { last }),
            ],
        });
        if shoulder {
            result.push(Boundary {
                edge: Edge::Circle {
                    radius: source.profile().end()[0],
                    z: a,
                },
                name: "screw.head/thread.junction/outer".into(),
                faces: [self.id(Role::ThreadShoulder), self.id(last)],
            });
        }
        for (i, segment) in profile.iter().enumerate() {
            result.push(Boundary {
                edge: Edge::Circle {
                    radius: segment.start[0],
                    z: segment.start[1] - source.spec().shank_length_mm,
                },
                name: format!("screw.head/{}/start", segment.role.source_role()),
                faces: [
                    self.id(if i == 0 {
                        Role::HeadTop
                    } else {
                        Role::ProfileWall(i as u32 - 1)
                    }),
                    self.id(Role::ProfileWall(i as u32)),
                ],
            });
        }
        for i in 0..12 {
            let wall = self.id(Role::SocketWall(i as u32));
            result.extend([
                Boundary {
                    edge: Edge::Mouth(i),
                    name: format!("screw.drive/rim/0/{i}"),
                    faces: [self.id(Role::HeadTop), wall],
                },
                Boundary {
                    edge: Edge::Floor(i),
                    name: format!("screw.drive/rim/1/{i}"),
                    faces: [self.id(Role::SocketFloor), wall],
                },
                Boundary {
                    edge: Edge::Post(i),
                    name: format!("screw.drive/post/{i}"),
                    faces: [self.id(Role::SocketWall(((i + 11) % 12) as u32)), wall],
                },
            ]);
        }
        result
    }
    fn edge_geometry(&self, edge: &Edge, frame: PlaneFrame) -> Option<AnalyticGeometry> {
        let source = self.geometry.head_source();
        let top = source.head_top_z();
        let drive = source
            .recess()
            .hexalobular_geometry()
            .expect("qualified curved drive");
        let circle = |center, radius, start, sweep| AnalyticGeometry::CircularArc {
            frame: ReferenceFrame {
                origin_mm: position(frame, center),
                u: frame.u,
                v: frame.v,
            },
            radius_mm: radius,
            start_rad: start,
            sweep_rad: sweep,
        };
        match *edge {
            Edge::Circle { radius, z } => Some(circle([0., 0., z], radius, 0., TAU)),
            Edge::Mouth(i) => {
                let arc = drive.arcs()[i];
                Some(circle(
                    [arc.center_mm[0], arc.center_mm[1], top],
                    arc.radius_mm,
                    arc.start_angle_rad,
                    arc.sweep_rad,
                ))
            }
            Edge::Floor(_) => None,
            Edge::Post(i) => {
                let [x, y] = drive.arcs()[i].start_mm;
                Some(AnalyticGeometry::Line {
                    a_mm: position(frame, [x, y, top]),
                    b_mm: position(
                        frame,
                        [
                            x,
                            y,
                            top + drive.cone().tip_depth_mm
                                - x.portable_hypot(y) * drive.cone().depth_per_radius,
                        ],
                    ),
                })
            }
        }
    }
    fn origin(&self, name: impl Into<String>) -> Origin {
        Origin {
            source: self.feature,
            role: name.into(),
        }
    }
    fn face_origin(&self, role: Role) -> Origin {
        self.origin(match role {
            Role::ThreadShoulder => "cap.start".into(),
            Role::ThreadTip => "cap.end".into(),
            Role::ThreadWall(0) => "thread.cosmetic.envelope".into(),
            Role::HeadTop => "screw.head/top.land".into(),
            Role::ProfileWall(i) => format!(
                "screw.head/{}",
                self.geometry
                    .head_source()
                    .profile()
                    .segments()
                    .nth(i as usize)
                    .expect("retained profile face")
                    .role
                    .source_role()
            ),
            Role::SocketFloor => "screw.drive/end".into(),
            Role::SocketWall(i) => format!("screw.drive/edge/{i}"),
            Role::ThreadWall(_) => unreachable!("one cosmetic shaft wall"),
        })
    }
    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> TopologySnapshot {
        let mut elements: Vec<_> = self
            .geometry
            .faces()
            .iter()
            .enumerate()
            .map(|(i, &role)| {
                let origin = self.face_origin(role);
                let plane = self.plane(role, frame);
                TopologyElement {
                    key: TopologyKey {
                        kind: TopologyKind::Face,
                        ordinal: i as u32,
                    },
                    origins: BTreeSet::from([origin.clone()]),
                    boundary: BTreeSet::new(),
                    frames: plane
                        .into_iter()
                        .map(|frame| AnchoredFrame {
                            origin: origin.clone(),
                            frame,
                        })
                        .collect(),
                    geometry: plane.map(|frame| AnalyticGeometry::Plane { frame }),
                }
            })
            .collect();
        for (index, edge) in self.boundaries().iter().enumerate() {
            let origin = self.origin(edge.name.clone());
            for &face in &edge.faces {
                elements[face].boundary.insert(origin.clone());
            }
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([origin]),
                boundary: BTreeSet::new(),
                frames: vec![],
                geometry: self.edge_geometry(&edge.edge, frame),
            });
        }
        let source = self.geometry.head_source();
        let drive = source
            .recess()
            .hexalobular_geometry()
            .expect("qualified curved drive");
        for (i, arc) in drive.arcs().iter().enumerate() {
            let [x, y] = arc.start_mm;
            for (end, z) in [
                (0, source.head_top_z()),
                (
                    1,
                    source.head_top_z() + drive.cone().tip_depth_mm
                        - x.portable_hypot(y) * drive.cone().depth_per_radius,
                ),
            ] {
                elements.push(TopologyElement {
                    key: TopologyKey {
                        kind: TopologyKind::Vertex,
                        ordinal: (2 * i + end) as u32,
                    },
                    origins: BTreeSet::from([self.origin(format!("screw.drive/corner/{end}/{i}"))]),
                    boundary: BTreeSet::new(),
                    frames: vec![],
                    geometry: Some(AnalyticGeometry::Point {
                        position_mm: position(frame, [x, y, z]),
                    }),
                });
            }
        }
        TopologySnapshot {
            owner,
            revision,
            elements,
        }
    }
    pub(crate) fn planar_loops(&self, face: u32) -> Vec<Vec<(u32, bool)>> {
        let Some(&role) = self.geometry.faces().get(face as usize) else {
            return vec![];
        };
        let Some(plane) = self.plane(role, PlaneFrame::default()) else {
            return vec![];
        };
        let edges = self.boundaries();
        let mut circles = Vec::new();
        let mut mouth = Vec::new();
        for (i, e) in edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.faces.contains(&(face as usize)))
        {
            match e.edge {
                Edge::Circle { radius, .. } => circles.push((radius, i as u32)),
                Edge::Mouth(_) => mouth.push((i as u32, false)),
                _ => {}
            }
        }
        circles.sort_by(|a, b| b.0.total_cmp(&a.0));
        let reverse = plane.normal()[2] < 0.;
        let mut loops: Vec<_> = circles
            .into_iter()
            .enumerate()
            .map(|(i, (_, edge))| vec![(edge, if i == 0 { reverse } else { !reverse })])
            .collect();
        if !mouth.is_empty() {
            loops.push(mouth);
        }
        loops
    }
    pub(super) fn mesh(
        &self,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        crate::rust_display::placed_mesh(self.geometry.mesh(options), frame, false)
    }
    pub(super) fn sample_edge(
        &self,
        frame: PlaneFrame,
        ordinal: u32,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Option<Result<Vec<Point>, multicad_rust_kernel::display::MeshError>> {
        use multicad_rust_kernel::{display::MeshError, hexalobular_drive::Error};
        let edges = self.boundaries();
        let Edge::Floor(arc) = edges.get(ordinal as usize)?.edge else {
            return None;
        };
        Some((|| {
            if !options.deflection.is_finite()
                || options.deflection <= 0.
                || !(0.01..=1.).contains(&options.angle)
                || !(1..=5_000_000).contains(&options.triangle_limit)
            {
                return Err(MeshError::InvalidOptions);
            }
            let source = self.geometry.head_source();
            let drive = source
                .recess()
                .hexalobular_geometry()
                .expect("qualified curved drive");
            let scale = [frame.u, frame.v, frame.normal()]
                .into_iter()
                .flatten()
                .fold(0f64, |a, b| (a + b.abs()).next_up());
            let request = (options.deflection / (4. * scale).next_up()).next_down();
            let local = drive
                .intersection_polyline(arc, request, options.triangle_limit.min(100_000))
                .map_err(|cause| {
                    if cause == Error::Limit {
                        MeshError::BudgetExceeded
                    } else {
                        MeshError::UnresolvedBoundary
                    }
                })?;
            let extent = local.iter().flatten().map(|p| p.abs()).fold(0f64, f64::max)
                + source.head_top_z().abs();
            let world = frame.origin.into_iter().map(f64::abs).fold(0f64, f64::max)
                + (scale * extent).next_up();
            if (64. * f64::EPSILON * world).next_up() > options.deflection / 4. {
                return Err(MeshError::UnresolvedBoundary);
            }
            Ok(local
                .into_iter()
                .map(|mut p| {
                    p[2] += source.head_top_z();
                    position(frame, p)
                })
                .collect())
        })())
    }
}

#[cfg(test)]
mod tests;
