//! Retained all-ring spline cap material; display never defines its geometry.
//! Initial source admission requires signed-cardinal unchanged stock/tool charts.
//! Repeated relief composition and general crossing curved Booleans are unqualified.
use super::*;
use crate::rust_prism_topology::{Builder, place_frame, position};
use multicad_core::{
    project_file::ContentHash,
    references::{AnalyticGeometry, AnchoredFrame, TopologyKind},
};
use multicad_rust_kernel::{
    normal_sweep::Edge as AnalyticEdge,
    spline_cap_relief::{Curve, FaceRole, Relief, Surface},
    spline_sweep::Edge,
};
use std::collections::BTreeSet;
mod imported;
mod source;
mod source_precision;
pub(super) use imported::signature as imported_signature;
#[derive(Clone)]
pub(super) struct Labels {
    origins: BTreeSet<Origin>,
    anchors: Vec<AnchoredFrame>,
}
pub(crate) struct SplineCapReliefBody {
    pub(super) geometry: Relief,
    pub(super) feature: EntityId,
    labels: Vec<Labels>,
    edge_map: Vec<Option<u32>>,
    areas: Vec<f64>,
    source_enclosure: source_precision::Receipt,
}
fn vector(frame: PlaneFrame, point: [f64; 3]) -> [f64; 3] {
    let axes = [frame.u, frame.v, frame.normal()];
    std::array::from_fn(|i| (0..3).map(|j| axes[j][i] * point[j]).sum())
}
fn fail(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error("spline cap Boolean", kind, message)
}
impl Solid {
    /// Authoritative exact stock, all motif counter loops and current placement.
    pub fn rust_spline_cap_relief_data(&self) -> Option<(&Relief, PlaneFrame)> {
        self.spline_cap_relief_body().map(|(b, f)| (&b.geometry, f))
    }
    /// Additional retained document-conversion/remap distance (mm) and normal
    /// angle (radians), separate from the pure graph's join/height reserve.
    pub fn rust_spline_cap_relief_source_enclosure(&self) -> Option<(f64, f64)> {
        self.spline_cap_relief_body()
            .map(|(b, _)| (b.source_enclosure.distance_mm, b.source_enclosure.angle_rad))
    }
    pub(crate) fn spline_cap_relief_body(&self) -> Option<(&SplineCapReliefBody, PlaneFrame)> {
        if let SolidBody::SplineCapRelief(body, frame) = &self.body {
            Some((body, *frame))
        } else {
            None
        }
    }
    pub fn boolean_spline_cap_relief(
        operands: &[&Self],
        operation: BooleanOperation,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<Self, GeometryError>> {
        source::build(operands, operation, feature, policy, interrupt)
    }
    pub fn rust_spline_cap_relief_contains(&self, point: [f64; 3]) -> Result<bool, GeometryError> {
        self.policy.point(point)?;
        let (body, frame) = self.spline_cap_relief_body().ok_or_else(|| {
            fail(
                GeometryErrorKind::Unsupported,
                "Not a retained spline cap relief",
            )
        })?;
        let [u, v, n] = [frame.u, frame.v, frame.normal()];
        let determinant = dot(u, cross(v, n));
        let delta = std::array::from_fn(|i| point[i] - frame.origin[i]);
        let local =
            [cross(v, n), cross(n, u), cross(u, v)].map(|row| dot(row, delta) / determinant);
        body.geometry
            .contains(local, self.policy.linear_tolerance_mm)
            .map_err(|e| {
                use multicad_rust_kernel::spline_sweep::ConstructionError as E;
                fail(
                    match e {
                        E::Limit => GeometryErrorKind::Limit,
                        E::Cancelled => GeometryErrorKind::Cancelled,
                        E::Invalid => GeometryErrorKind::Input,
                        E::Unresolved => GeometryErrorKind::Unsupported,
                    },
                    format!("The finite material query is unresolved: {e:?}"),
                )
            })
    }
}
impl SplineCapReliefBody {
    fn new(
        geometry: Relief,
        feature: EntityId,
        labels: Vec<Labels>,
        policy: GeometryPolicy,
        source_enclosure: source_precision::Receipt,
    ) -> Result<Self, GeometryError> {
        let areas = geometry
            .boundary()
            .faces
            .iter()
            .enumerate()
            .map(|(i, face)| {
                if let Some((area, _)) = geometry.cap_properties(i) {
                    return Ok(area);
                }
                let Surface::Extrusion { edge, span_mm, .. } = &face.surface else {
                    unreachable!("cap properties")
                };
                let bound = edge
                    .length_bounds_mm(policy.linear_tolerance_mm / 16.)
                    .map_err(|_| {
                        fail(
                            GeometryErrorKind::Limit,
                            "Spline wall length certification exhausted its budget",
                        )
                    })?;
                Ok(bound[0].midpoint(bound[1]) * (span_mm[1] - span_mm[0]))
            })
            .collect::<Result<Vec<_>, GeometryError>>()?;
        let graph = geometry.boundary();
        let mut incident = vec![vec![]; graph.edges.len()];
        for (f, face) in graph.faces.iter().enumerate() {
            for e in face.loops.iter().flatten() {
                incident[e.edge].push((f, e.forward));
            }
        }
        let mut ordinal = 0;
        let edge_map = incident
            .iter()
            .map(|uses| {
                let [(a, forward), (b, reverse)] = uses[..] else {
                    return Err(fail(
                        GeometryErrorKind::InvalidResult,
                        "Invalid shared boundary incidence",
                    ));
                };
                if forward == reverse {
                    return Err(fail(
                        GeometryErrorKind::InvalidResult,
                        "Shared boundary winding is inconsistent",
                    ));
                }
                if a == b {
                    return Ok(None);
                }
                let id = ordinal;
                ordinal += 1;
                Ok(Some(id))
            })
            .collect::<Result<Vec<_>, GeometryError>>()?;
        Ok(Self {
            geometry,
            feature,
            labels,
            edge_map,
            areas,
            source_enclosure,
        })
    }
    pub(super) fn volume_mm3(&self) -> f64 {
        self.geometry.volume_mm3()
    }
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Bounds {
        let n = frame.normal();
        Bounds {
            min: std::array::from_fn(|i| {
                frame.origin[i] - self.geometry.support([-frame.u[i], -frame.v[i], -n[i]])
            }),
            max: std::array::from_fn(|i| {
                frame.origin[i] + self.geometry.support([frame.u[i], frame.v[i], n[i]])
            }),
        }
    }
    fn plane(&self, index: usize, frame: PlaneFrame) -> Option<ReferenceFrame> {
        match &self.geometry.boundary().faces[index].surface {
            Surface::Plane { z_mm, upward } => Some(ReferenceFrame {
                origin_mm: position(frame, [0., 0., *z_mm]),
                u: frame.u,
                v: frame.v.map(|x| if *upward { x } else { -x }),
            }),
            Surface::Extrusion {
                edge: Edge::Analytic(AnalyticEdge::Line { a, b }),
                span_mm,
                inward,
            } => {
                let length = (b[0] - a[0]).hypot(b[1] - a[1]);
                Some(ReferenceFrame {
                    origin_mm: position(frame, [a[0], a[1], span_mm[0]]),
                    u: vector(frame, [(b[0] - a[0]) / length, (b[1] - a[1]) / length, 0.]),
                    v: frame.normal().map(|x| if *inward { -x } else { x }),
                })
            }
            _ => None,
        }
    }
    pub(super) fn faces(&self, frame: PlaneFrame) -> Vec<Face> {
        self.geometry
            .boundary()
            .faces
            .iter()
            .enumerate()
            .map(|(i, face)| {
                let center = self
                    .geometry
                    .cap_properties(i)
                    .map(|(_, p)| p)
                    .unwrap_or_else(|| {
                        let Surface::Extrusion { edge, span_mm, .. } = &face.surface else {
                            unreachable!("cap properties")
                        };
                        let p = match edge {
                            Edge::Analytic(a) => a.centroid(),
                            Edge::Cubic(_) | Edge::Rational(_) => edge.point(0.5),
                        };
                        [p[0], p[1], span_mm[0].midpoint(span_mm[1])]
                    });
                Face {
                    index: i as u32,
                    plane: self.plane(i, frame).map(|p| (p.origin_mm, p.normal())),
                    center_mm: position(frame, center),
                    area_mm2: self.areas[i],
                }
            })
            .collect()
    }
    pub(super) fn mesh(
        &self,
        frame: PlaneFrame,
        mut options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        options.deflection =
            source_precision::remaining(options.deflection, self.source_enclosure.distance_mm);
        options.angle = source_precision::remaining(options.angle, self.source_enclosure.angle_rad);
        if options.deflection <= 0. || options.angle <= 0. {
            return Err(fail(
                GeometryErrorKind::InvalidResult,
                "The requested mesh is below the retained source enclosure",
            ));
        }
        crate::rust_display::placed_mesh(self.geometry.mesh(options), frame, false)
    }
    pub(super) fn append_step(
        &self,
        writer: &mut multicad_rust_kernel::step::Writer,
        frame: PlaneFrame,
    ) -> Result<(), multicad_rust_kernel::step::StepError> {
        writer.append_spline_cap_relief(&self.geometry, frame.origin, frame.u, frame.v)
    }
    pub(crate) fn planar_loops(&self, index: u32) -> Vec<Vec<(u32, bool)>> {
        let Some(face) = self.geometry.boundary().faces.get(index as usize) else {
            return vec![];
        };
        if self.plane(index as usize, PlaneFrame::default()).is_none() {
            return vec![];
        }
        face.loops
            .iter()
            .map(|ring| {
                ring.iter()
                    .map(|e| {
                        (
                            self.edge_map[e.edge].expect("physical planar rim"),
                            !e.forward,
                        )
                    })
                    .collect()
            })
            .collect()
    }
    pub(super) fn sample_rim(
        &self,
        frame: PlaneFrame,
        ordinal: u32,
        mut options: multicad_rust_kernel::display::Tessellation,
    ) -> Option<Result<Vec<[f64; 3]>, multicad_rust_kernel::display::MeshError>> {
        let index = self.edge_map.iter().position(|id| *id == Some(ordinal))?;
        let Curve::Profile {
            edge: edge @ Edge::Cubic(_),
            z_mm,
        } = &self.geometry.boundary().edges[index].curve
        else {
            return None;
        };
        options.deflection =
            source_precision::remaining(options.deflection, self.source_enclosure.distance_mm);
        options.angle = source_precision::remaining(options.angle, self.source_enclosure.angle_rad);
        if options.deflection <= 0. || options.angle <= 0. {
            return Some(Err(
                multicad_rust_kernel::display::MeshError::UnresolvedBoundary,
            ));
        }
        Some(edge.sample(options).map(|points| {
            points
                .into_iter()
                .map(|p| position(frame, [p[0], p[1], *z_mm]))
                .collect()
        }))
    }
}

impl SplineCapReliefBody {
    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> TopologySnapshot {
        let boundary = self.geometry.boundary();
        let mut builder = Builder::new(
            self.labels
                .iter()
                .map(|labels| labels.origins.iter().cloned().collect())
                .collect(),
        );
        let mut edge_faces = vec![BTreeSet::new(); boundary.edges.len()];
        let mut vertex_faces = vec![BTreeSet::new(); boundary.vertices.len()];
        for (face, data) in boundary.faces.iter().enumerate() {
            for use_ in data.loops.iter().flatten() {
                edge_faces[use_.edge].insert(face);
                for vertex in boundary.edges[use_.edge].vertices {
                    vertex_faces[vertex].insert(face);
                }
            }
        }
        for (index, edge) in boundary.edges.iter().enumerate() {
            if self.edge_map[index].is_none() {
                continue;
            }
            if matches!(
                &edge.curve,
                Curve::Profile {
                    edge: Edge::Cubic(_) | Edge::Rational(_),
                    ..
                }
            ) {
                builder.non_analytic_edge(&edge_faces[index].iter().copied().collect::<Vec<_>>());
                continue;
            }
            let geometry = match &edge.curve {
                Curve::Line { a_mm, b_mm } => AnalyticGeometry::Line {
                    a_mm: position(frame, *a_mm),
                    b_mm: position(frame, *b_mm),
                },
                Curve::Profile {
                    edge: Edge::Analytic(edge),
                    z_mm,
                } => super::normal_sweep_body::rim_in_frame(
                    *edge,
                    ReferenceFrame {
                        origin_mm: position(frame, [0., 0., *z_mm]),
                        u: frame.u,
                        v: frame.v,
                    },
                ),
                Curve::Profile {
                    edge: Edge::Cubic(_) | Edge::Rational(_),
                    ..
                } => unreachable!("handled spline"),
            };
            builder.non_face(
                geometry,
                &edge_faces[index].iter().copied().collect::<Vec<_>>(),
            );
        }
        // A periodic circle's private STEP seam is not a physical CAD vertex.
        let mut physical_vertices = vec![false; boundary.vertices.len()];
        for (index, edge) in boundary.edges.iter().enumerate() {
            if self.edge_map[index].is_some()
                && !matches!(
                    &edge.curve,
                    Curve::Profile {
                        edge: Edge::Analytic(AnalyticEdge::Circle { .. }),
                        ..
                    }
                )
            {
                for vertex in edge.vertices {
                    physical_vertices[vertex] = true;
                }
            }
        }
        for (index, point) in boundary.vertices.iter().enumerate() {
            if physical_vertices[index] {
                builder.non_face(
                    AnalyticGeometry::Point {
                        position_mm: position(frame, *point),
                    },
                    &vertex_faces[index].iter().copied().collect::<Vec<_>>(),
                );
            }
        }
        for (index, face) in boundary.faces.iter().enumerate() {
            let plane = self.plane(index, frame);
            builder.face(
                index,
                plane,
                plane,
                face.loops
                    .iter()
                    .flatten()
                    .flat_map(|edge| edge_faces[edge.edge].iter().copied())
                    .filter(|neighbor| *neighbor != index),
            );
        }
        let mut elements = builder.into_elements();
        for element in &mut elements {
            if element.key.kind == TopologyKind::Face {
                let labels = &self.labels[element.key.ordinal as usize];
                for anchor in &labels.anchors {
                    element
                        .frames
                        .retain(|candidate| candidate.origin != anchor.origin);
                    element.frames.push(AnchoredFrame {
                        origin: anchor.origin.clone(),
                        frame: place_frame(frame, anchor.frame),
                    });
                }
            }
        }
        TopologySnapshot {
            owner,
            revision,
            elements,
        }
    }
}
