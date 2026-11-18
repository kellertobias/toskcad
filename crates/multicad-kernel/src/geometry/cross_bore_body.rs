//! A stepped shaft with one perpendicular cross-bore on its retained
//! supports: the shaft's own revolution faces, the tool cylinder and planar
//! blind bottoms. The bored wall meets the tool in the exact quartic rim of
//! the kernel provider; display and manufacturing triangles are sampled from
//! that boundary and never define it.
use super::revolution_body::RevolutionBody;
use super::*;
use crate::rust_prism_topology::{place_frame, position};
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, ReferenceFrame, TopologyElement, TopologyKey, TopologyKind,
    TopologySnapshot,
};
use multicad_rust_kernel::{
    cross_bore::{self, CrossBore, Cut, End, Tool},
    display::{MeshError, Tessellation},
};
use std::collections::BTreeSet;

const OPERATION: &str = "cross-bore cut";

/// Authored identity of one tool face, its frames kept literally and in the
/// shaft chart for later placements.
struct Source {
    origins: Vec<Origin>,
    frames: Vec<(AnchoredFrame, ReferenceFrame)>,
}

pub(crate) struct CrossBoreBody {
    pub(super) shaft: Rc<RevolutionBody>,
    pub(super) geometry: CrossBore,
    pub(super) feature: EntityId,
    /// Tool start cap, end cap and side.
    tool: [Source; 3],
    /// Placement that captured the literal source frames.
    authored: PlaneFrame,
}

fn refusal(kind: GeometryErrorKind, message: &str) -> GeometryError {
    error(OPERATION, kind, message)
}

fn kernel_error(cause: cross_bore::Error) -> GeometryError {
    match cause {
        cross_bore::Error::Invalid => refusal(
            GeometryErrorKind::Input,
            "Invalid shaft or cylinder for a cross-bore",
        ),
        cross_bore::Error::Unsupported => refusal(
            GeometryErrorKind::Unsupported,
            "Rust cross-bores support a tool section inside one cylindrical wall of a solid \
             stepped shaft, through or with a blind bottom inside the shaft; tangency, a \
             chordal cut, a shoulder crossing and an enclosed void are not supported",
        ),
        cross_bore::Error::Unresolved => refusal(
            GeometryErrorKind::Unsupported,
            "Cross-bore classification is unresolved within the certified margin",
        ),
        cross_bore::Error::Cancelled => {
            refusal(GeometryErrorKind::Cancelled, "Cross-bore cut was cancelled")
        }
    }
}

fn direction(frame: PlaneFrame, local: [f64; 3]) -> [f64; 3] {
    position(
        PlaneFrame {
            origin: [0.; 3],
            ..frame
        },
        local,
    )
}

fn sources(snapshot: &TopologySnapshot, ordinal: usize, chart: PlaneFrame) -> Source {
    let face = snapshot
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal as usize == ordinal);
    let axes = [chart.u, chart.v, chart.normal()];
    let local = |vector: [f64; 3]| axes.map(|axis| dot(vector, axis));
    let frames = face.map_or(vec![], |face| {
        face.frames
            .iter()
            .map(|anchor| {
                let charted = ReferenceFrame {
                    origin_mm: local(std::array::from_fn(|i| {
                        anchor.frame.origin_mm[i] - chart.origin[i]
                    })),
                    u: local(anchor.frame.u),
                    v: local(anchor.frame.v),
                };
                (anchor.clone(), charted)
            })
            .collect()
    });
    Source {
        origins: face.map_or(vec![], |face| face.origins.iter().cloned().collect()),
        frames,
    }
}

/// Whether every component is 0 or ±1, so a chart change is exact.
fn cardinal(vectors: &[[f64; 3]]) -> bool {
    vectors.iter().flatten().all(|c| *c == 0. || c.abs() == 1.)
}

impl Solid {
    /// The shaft revolution and its placement for a cross-bore target: a
    /// stepped revolution or a plain axial cylinder.
    fn cross_bore_shaft(&self) -> Option<Result<(Rc<RevolutionBody>, PlaneFrame), GeometryError>> {
        if let SolidBody::Revolved(body, frame) = &self.body {
            return Some(Ok((Rc::clone(body), *frame)));
        }
        let converted = self.cylinder_as_revolution()?;
        Some(converted.map(|solid| match solid.body {
            SolidBody::Revolved(body, frame) => (body, frame),
            _ => unreachable!("Rust revolution provider"),
        }))
    }

    /// `None` when the operands belong to another provider: the target is no
    /// solid stepped shaft or the tool is no cylinder perpendicular to it.
    pub(in crate::geometry) fn cross_bore_cut(
        &self,
        tool: &Self,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        let (cylinder, tool_frame) = tool.axial_cylinder_frame()?;
        if !matches!(
            self.body,
            SolidBody::Revolved(..) | SolidBody::AxialCylinder(..)
        ) {
            return None;
        }
        let (shaft, frame) = match self.cross_bore_shaft()? {
            Ok(found) => found,
            Err(error) => return Some(Err(error)),
        };
        cross_bore::stepped_walls(&shaft.geometry)?;
        let axes = [frame.u, frame.v, frame.normal()];
        let local = |vector: [f64; 3]| axes.map(|axis| dot(vector, axis));
        let height = cylinder.height().as_f64();
        let mut travel = local(tool_frame.normal());
        let scale = shaft.geometry.support([1., 0., 0.]).max(
            shaft
                .geometry
                .support([0., 0., 1.])
                .max(-shaft.geometry.support([0., 0., -1.])),
        );
        // A travel within chart noise of perpendicular is perpendicular; any
        // other inclination belongs to another provider.
        let numeric = 32. * f64::EPSILON * scale.max(1.);
        let snapped = travel[2] != 0.;
        if (travel[2] * height).abs() > numeric {
            return None;
        }
        travel[2] = 0.;
        let norm = travel[0].hypot(travel[1]);
        travel = travel.map(|c| c / norm);
        let base: [f64; 3] = std::array::from_fn(|i| tool_frame.origin[i] - frame.origin[i]);
        let binary = |length: Rational| {
            let (numerator, denominator) = length.ratio();
            denominator.count_ones() == 1 && numerator.unsigned_abs() <= 1 << 53
        };
        let radius = cylinder.circle().radius();
        let literal = !snapped
            && binary(radius)
            && binary(cylinder.height())
            && cardinal(&axes)
            && cardinal(&[tool_frame.normal()])
            && (0..3).all(|i| {
                base[i] + frame.origin[i] == tool_frame.origin[i]
                    && tool_frame.origin[i] - base[i] == frame.origin[i]
            });
        let cutter = Tool {
            base_mm: local(base),
            axis: travel,
            radius_mm: radius.as_f64(),
            height_mm: height,
            literal,
        };
        Some(self.cross_bore(tool, shaft, frame, cutter, interrupt))
    }

    fn cross_bore(
        &self,
        tool: &Self,
        shaft: Rc<RevolutionBody>,
        frame: PlaneFrame,
        cutter: Tool,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<Self>, GeometryError> {
        let geometry =
            match cross_bore::subtract(&shaft.geometry, cutter, interrupt).map_err(kernel_error)? {
                Cut::Unchanged => return self.share().map(Some),
                Cut::Body(body) => *body,
            };
        let owner = EntityId::new();
        let cylinder = tool.topology_snapshot(owner, 0)?;
        let body = CrossBoreBody {
            feature: shaft.feature,
            shaft,
            geometry,
            tool: [0, 1, 2].map(|ordinal| sources(&cylinder, ordinal, frame)),
            authored: frame,
        };
        let bounds = body.bounds(frame);
        self.policy.point(bounds.min)?;
        self.policy.point(bounds.max)?;
        Ok(Some(Self {
            body: SolidBody::CrossBore(Rc::new(body), frame),
            policy: self.policy,
            provenance: vec![],
            rust_origins: self.rust_origins.clone(),
            retained_placement: None,
        }))
    }

    pub(crate) fn cross_bore_body(&self) -> Option<(&CrossBoreBody, PlaneFrame)> {
        match &self.body {
            SolidBody::CrossBore(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
}

impl CrossBoreBody {
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Bounds {
        self.shaft.bounds(frame)
    }

    fn bored_face(&self) -> u32 {
        self.shaft
            .geometry
            .wall_face(0, self.geometry.wall().0)
            .expect("validated bored wall")
    }

    /// Ordinals of the tool wall and each blind bottom.
    fn new_faces(&self) -> (u32, [Option<u32>; 2]) {
        let tool = self.geometry.shaft_faces();
        let mut next = tool + 1;
        let caps = [0, 1].map(|side| {
            (self.geometry.ends()[side] == End::Cap).then(|| {
                next += 1;
                next - 1
            })
        });
        (tool, caps)
    }

    pub(super) fn faces(&self, frame: PlaneFrame) -> Vec<Face> {
        let mut faces = self.shaft.faces(frame);
        let bored = self.bored_face();
        let wall = self.geometry.bored_wall();
        for face in faces.iter_mut().filter(|face| face.index == bored) {
            face.area_mm2 = wall.area_mm2;
            face.center_mm = position(frame, wall.centroid_mm);
        }
        let (tool_face, caps) = self.new_faces();
        let tool_wall = self.geometry.tool_wall();
        faces.push(Face {
            index: tool_face,
            plane: None,
            center_mm: position(frame, tool_wall.centroid_mm),
            area_mm2: tool_wall.area_mm2,
        });
        for (side, index) in caps.into_iter().enumerate() {
            let (Some(index), Some((center, normal, radius))) = (index, self.geometry.cap(side))
            else {
                continue;
            };
            faces.push(Face {
                index,
                plane: Some((position(frame, center), direction(frame, normal))),
                center_mm: position(frame, center),
                area_mm2: std::f64::consts::PI * radius * radius,
            });
        }
        faces
    }

    pub(super) fn mesh(
        &self,
        frame: PlaneFrame,
        options: Tessellation,
    ) -> Result<Mesh, GeometryError> {
        crate::rust_display::placed_mesh(self.geometry.mesh(options), frame, false)
    }

    /// Edge ordinals after the shaft's: each through rim, then each bottom
    /// circle, by side.
    fn new_edges(&self, shaft_edges: u32) -> [Option<u32>; 2] {
        let mut next = shaft_edges;
        [0, 1].map(|_| {
            next += 1;
            Some(next - 1)
        })
    }

    fn shaft_edges(&self, frame: PlaneFrame) -> u32 {
        self.shaft
            .topology(frame, EntityId::new(), 0)
            .elements
            .iter()
            .filter(|element| element.key.kind == TopologyKind::Edge)
            .count() as u32
    }

    /// Samples of a bored-wall rim; bottom circles are analytic edges.
    pub(super) fn sample_edge(
        &self,
        frame: PlaneFrame,
        ordinal: u32,
        options: Tessellation,
    ) -> Option<Result<Vec<[f64; 3]>, MeshError>> {
        let edges = self.new_edges(self.shaft_edges(frame));
        let side = edges.iter().position(|edge| *edge == Some(ordinal))?;
        if self.geometry.ends()[side] != End::Through {
            return None;
        }
        Some(self.geometry.sample_rim(side, options).map(|points| {
            points
                .into_iter()
                .map(|point| position(frame, point))
                .collect()
        }))
    }

    /// The bored wall rim or bottom circle as one closed planar loop.
    pub(crate) fn planar_loops(&self, face: u32) -> Vec<Vec<(u32, bool)>> {
        let (_, caps) = self.new_faces();
        if let Some(side) = caps.iter().position(|cap| *cap == Some(face)) {
            let edges = self.new_edges(self.shaft_edges(PlaneFrame::default()));
            return vec![vec![(edges[side].expect("bottom edge"), false)]];
        }
        if face < self.geometry.shaft_faces() {
            return self.shaft.planar_loops(face);
        }
        Vec::new()
    }

    fn placed_sources(&self, source: &Source, frame: PlaneFrame) -> Vec<AnchoredFrame> {
        source
            .frames
            .iter()
            .map(|(literal, charted)| AnchoredFrame {
                origin: literal.origin.clone(),
                frame: if frame == self.authored {
                    literal.frame
                } else {
                    place_frame(frame, *charted)
                },
            })
            .collect()
    }

    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, String> {
        let mut snapshot = self.shaft.topology(frame, owner, revision);
        let shaft_edges = snapshot
            .elements
            .iter()
            .filter(|element| element.key.kind == TopologyKind::Edge)
            .count() as u32;
        let bored = self.bored_face();
        let bored_origins = snapshot
            .elements
            .iter()
            .find(|element| {
                element.key
                    == TopologyKey {
                        kind: TopologyKind::Face,
                        ordinal: bored,
                    }
            })
            .map(|element| element.origins.clone())
            .ok_or("Shaft topology does not publish the bored wall")?;
        let (tool_face, caps) = self.new_faces();
        let tool_origins: BTreeSet<Origin> = self.tool[2].origins.iter().cloned().collect();
        let cap_origins: [BTreeSet<Origin>; 2] =
            [0, 1].map(|side| self.tool[side].origins.iter().cloned().collect());
        let through: Vec<usize> = (0..2)
            .filter(|&side| self.geometry.ends()[side] == End::Through)
            .collect();
        let mut faces = Vec::new();
        let mut tool_boundary = BTreeSet::new();
        if !through.is_empty() {
            tool_boundary.extend(bored_origins.iter().cloned());
        }
        for side in 0..2 {
            if caps[side].is_some() {
                tool_boundary.extend(cap_origins[side].iter().cloned());
            }
        }
        faces.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Face,
                ordinal: tool_face,
            },
            origins: tool_origins.clone(),
            boundary: tool_boundary,
            frames: self.placed_sources(&self.tool[2], frame),
            geometry: None,
        });
        for side in 0..2 {
            let (Some(index), Some((center, normal, _))) = (caps[side], self.geometry.cap(side))
            else {
                continue;
            };
            // Normal = u × v for a normal perpendicular to the shaft axis.
            let plane = ReferenceFrame {
                origin_mm: center,
                u: cross([0., 0., 1.], normal),
                v: [0., 0., 1.],
            };
            faces.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: index,
                },
                origins: cap_origins[side].clone(),
                boundary: tool_origins.clone(),
                frames: self.placed_sources(&self.tool[side], frame),
                geometry: Some(AnalyticGeometry::Plane {
                    frame: place_frame(frame, plane),
                }),
            });
        }
        let edges = self.new_edges(shaft_edges);
        let mut boundaries = Vec::new();
        for side in 0..2 {
            let ordinal = edges[side].expect("allocated edge");
            let (geometry, partner) = match self.geometry.ends()[side] {
                End::Through => (None, &bored_origins),
                End::Cap => {
                    let (center, normal, radius) = self.geometry.cap(side).expect("blind bottom");
                    let circle = ReferenceFrame {
                        origin_mm: center,
                        u: cross([0., 0., 1.], normal),
                        v: [0., 0., 1.],
                    };
                    (
                        Some(AnalyticGeometry::CircularArc {
                            frame: place_frame(frame, circle),
                            radius_mm: radius,
                            start_rad: 0.,
                            sweep_rad: std::f64::consts::TAU,
                        }),
                        &cap_origins[side],
                    )
                }
            };
            let origins: BTreeSet<Origin> = partner.iter().chain(&tool_origins).cloned().collect();
            boundaries.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal,
                },
                origins: origins.clone(),
                boundary: origins,
                frames: vec![],
                geometry,
            });
        }
        if let Some(face) = snapshot.elements.iter_mut().find(|element| {
            element.key
                == TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: bored,
                }
        }) && !through.is_empty()
        {
            face.boundary.extend(tool_origins.iter().cloned());
        }
        let face_count = snapshot
            .elements
            .iter()
            .take_while(|element| element.key.kind == TopologyKind::Face)
            .count();
        let tail = snapshot.elements.split_off(face_count);
        snapshot.elements.extend(faces);
        snapshot.elements.extend(tail);
        snapshot.elements.extend(boundaries);
        snapshot.validate()?;
        Ok(snapshot)
    }

    /// Chord loss of the shared tessellation: the shaft's own bound plus the
    /// sagitta of the bore's rims over its curved faces.
    pub(super) fn volume_error_bound(&self, options: Tessellation) -> f64 {
        let (_, radius, [low, high]) = self.geometry.wall();
        let holes =
            std::f64::consts::TAU * radius * (high - low) - self.geometry.bored_wall().area_mm2;
        let curved = self.geometry.tool_wall().area_mm2 + holes;
        multicad_rust_kernel::revolution_mesh::volume_error_bound(self.geometry.shaft(), options)
            + 2. * options.deflection * curved
    }

    pub(super) fn contains(&self, frame: PlaneFrame, point: [f64; 3], tolerance: f64) -> bool {
        let delta: [f64; 3] = std::array::from_fn(|i| point[i] - frame.origin[i]);
        let local = [frame.u, frame.v, frame.normal()].map(|axis| dot(axis, delta));
        self.geometry.contains(local, tolerance)
    }
}
