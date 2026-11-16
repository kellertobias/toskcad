//! An imported solid that keeps its exchange boundary: support surfaces,
//! trim loops, edge curves and vertices, exactly as the file states them.
//! Volume, bounds, meshes, sections and STEP export read that boundary.
use super::*;
use multicad_core::references::{
    AnalyticGeometry, ReferenceFrame, TopologyElement, TopologyKey, TopologyKind, TopologySnapshot,
};
use multicad_rust_kernel::brep::{self, Curve, Surface};
use multicad_rust_kernel::display::{MeshError, Tessellation};
use multicad_rust_kernel::step_read;
use std::collections::BTreeMap;

const OPERATION: &str = "imported boundary";

/// The most topology elements one retained solid may publish.
const MAX_ELEMENTS: usize = 100_000;

pub(crate) struct ImportedBrepBody {
    pub(super) solid: brep::Solid,
    pub(super) feature: EntityId,
    /// Stable role of every face, edge and vertex, taken from the source
    /// geometry so that a placed copy keeps its references.
    roles: Rc<[Vec<String>; 3]>,
}

pub(super) fn failure(cause: brep::Error) -> GeometryError {
    let (kind, message) = match cause {
        brep::Error::Invalid => (
            GeometryErrorKind::InvalidResult,
            "Imported boundary is not a valid closed solid",
        ),
        brep::Error::Unsupported => (
            GeometryErrorKind::Unsupported,
            "Imported boundary is outside the retained Rust families",
        ),
        brep::Error::Limit => (
            GeometryErrorKind::Limit,
            "Imported boundary exceeds the Rust budget",
        ),
        brep::Error::Cancelled => (GeometryErrorKind::Cancelled, "Imported boundary superseded"),
    };
    error(OPERATION, kind, message)
}

fn coordinates(points: impl IntoIterator<Item = [f64; 3]>) -> Vec<f64> {
    points.into_iter().flatten().collect()
}

/// Roles from geometry alone; equal geometry is told apart by occurrence.
fn roles(solid: &brep::Solid) -> [Vec<String>; 3] {
    let boundary = solid.boundary();
    let mut seen = BTreeMap::<String, u32>::new();
    let mut role = |kind: &str, values: Vec<f64>| {
        let base = step_import::feature_origin(EntityId::new(), kind, values).role;
        let count = seen.entry(base.clone()).or_default();
        *count += 1;
        if *count == 1 {
            base
        } else {
            format!("{base}#{count}")
        }
    };
    let faces = solid
        .faces()
        .iter()
        .map(|face| {
            role(
                "face",
                coordinates([face.centroid_mm, [face.area_mm2, 0., 0.]]),
            )
        })
        .collect();
    let edges = boundary
        .edges
        .iter()
        .map(|edge| {
            let middle = edge.curve.evaluate(0.5 * (edge.range[0] + edge.range[1])).0;
            let ends = edge.vertices.map(|vertex| boundary.vertices[vertex]);
            role("edge", coordinates([ends[0], middle, ends[1]]))
        })
        .collect();
    let vertices = boundary
        .vertices
        .iter()
        .map(|&point| role("vertex", coordinates([point])))
        .collect();
    [faces, edges, vertices]
}

impl Solid {
    /// Import every solid of a flat STEP file by retaining its boundary.
    /// Several solids form one compound in file order.
    pub(super) fn import_rust_brep(
        input: &str,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let solids = step_read::brep::read_solids(input, policy.linear_tolerance_mm, None)
            .map_err(failure)?;
        let single = solids.len() == 1;
        let mut members = Vec::with_capacity(solids.len());
        for (index, solid) in solids.into_iter().enumerate() {
            let child = if single {
                feature
            } else {
                let mut bytes = b"multicad.step-boundary/1".to_vec();
                bytes.extend(feature.uuid_bits().to_le_bytes());
                bytes.extend((index as u64).to_le_bytes());
                let hash = multicad_core::project_file::ContentHash::of(&bytes);
                EntityId::from_uuid_bits(u128::from_le_bytes(hash.0[..16].try_into().unwrap()))
            };
            members.push(Self::from_imported_brep(solid, child, policy)?);
        }
        if single {
            return Ok(members.remove(0));
        }
        let members: Vec<_> = members.iter().collect();
        Self::make_rust_compound(&members, feature, policy, None)
    }

    fn from_imported_brep(
        solid: brep::Solid,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let boundary = solid.boundary();
        let elements = boundary.faces.len() + boundary.edges.len() + boundary.vertices.len();
        if elements > MAX_ELEMENTS {
            return Err(failure(brep::Error::Limit));
        }
        let [low, high] = solid.bounds_mm();
        policy.point(low)?;
        policy.point(high)?;
        let roles = Rc::new(roles(&solid));
        Ok(Self {
            body: SolidBody::ImportedBrep(Rc::new(ImportedBrepBody {
                solid,
                feature,
                roles,
            })),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }

    /// Verify a STEP text written from this retained boundary, or `None`
    /// when the solid does not retain imported boundaries.
    pub fn verify_imported_brep_step(
        &self,
        source: &str,
    ) -> Option<Result<f64, step_read::ReadError>> {
        let members: Vec<&brep::Solid> = match self.imported_brep_body() {
            Some(body) => vec![&body.solid],
            None => self
                .rust_compound()?
                .members()
                .map(|member| member.imported_brep_body().map(|body| &body.solid))
                .collect::<Option<_>>()?,
        };
        Some(step_read::brep::verify_retained(
            source,
            &members,
            self.policy.linear_tolerance_mm,
        ))
    }

    /// The solid's exact STEP export read back as a retained boundary, or
    /// `None` when it has no exact export or the export does not read back.
    pub(crate) fn retained_boundary(&self) -> Option<Solid> {
        if self.imported_brep_body().is_some() {
            return self.share().ok();
        }
        if !self.supports_step_export() {
            return None;
        }
        let mut writer = multicad_rust_kernel::step::Writer::new();
        self.append_step(&mut writer).ok()?;
        let text = writer.finish().ok()?;
        Self::import_rust_brep(&text, EntityId::new(), self.policy).ok()
    }

    pub(crate) fn imported_brep_body(&self) -> Option<&ImportedBrepBody> {
        match &self.body {
            SolidBody::ImportedBrep(body) => Some(body),
            _ => None,
        }
    }
}

impl ImportedBrepBody {
    pub(super) fn transformed(&self, transform: RigidTransform) -> Result<Self, GeometryError> {
        let linear = |v: [f64; 3]| transform.vector(v);
        let map = |p: [f64; 3]| {
            let turned = transform.vector(p);
            std::array::from_fn(|i| turned[i] + transform.translation_mm[i])
        };
        Ok(Self {
            solid: self.solid.moved(&map, &linear).map_err(failure)?,
            feature: self.feature,
            roles: Rc::clone(&self.roles),
        })
    }

    pub(super) fn bounds(&self) -> Bounds {
        let [min, max] = self.solid.bounds_mm();
        Bounds { min, max }
    }

    pub(super) fn faces(&self) -> Vec<Face> {
        let boundary = self.solid.boundary();
        let signs = boundary
            .faces
            .iter()
            .map(|face| if face.same_sense { 1. } else { -1. });
        self.solid
            .faces()
            .iter()
            .zip(&boundary.faces)
            .zip(signs)
            .enumerate()
            .map(|(index, ((properties, face), sign))| Face {
                index: index as u32,
                plane: match &face.surface {
                    Surface::Plane(frame) => Some((frame.origin, frame.z.map(|x| x * sign))),
                    _ => None,
                },
                center_mm: properties.centroid_mm,
                area_mm2: properties.area_mm2,
            })
            .collect()
    }

    pub(super) fn mesh(&self, options: Tessellation) -> Result<Mesh, GeometryError> {
        let options = brep::mesh::Options {
            deflection_mm: options.deflection,
            angle_rad: options.angle,
            triangle_limit: options.triangle_limit,
        };
        let mesh = brep::mesh::tessellate(&self.solid, options, None).map_err(failure)?;
        crate::rust_display::mesh(Ok(mesh.display()))
    }

    /// Plane section loops with exact lines and ellipses where the plane
    /// meets a face in one, chords elsewhere.
    pub(crate) fn section_edges(
        &self,
        origin: [f64; 3],
        normal: [f64; 3],
        deflection_mm: f64,
    ) -> Result<Vec<Vec<multicad_rust_kernel::warped_union::SectionEdge>>, GeometryError> {
        let options = brep::mesh::Options {
            deflection_mm,
            angle_rad: 0.3,
            triangle_limit: 5_000_000,
        };
        brep::section::section_edges(&self.solid, origin, normal, options, None).map_err(failure)
    }

    /// Sample an edge whose curve has no analytic display form.
    pub(super) fn sample_edge(
        &self,
        ordinal: u32,
        options: Tessellation,
    ) -> Option<Result<Vec<[f64; 3]>, MeshError>> {
        let edge = self.solid.boundary().edges.get(ordinal as usize)?;
        Some(
            edge.polyline(options.deflection, options.angle)
                .map_err(|_| MeshError::BudgetExceeded),
        )
    }

    pub(crate) fn topology(&self, owner: EntityId, revision: u64) -> TopologySnapshot {
        let boundary = self.solid.boundary();
        let origin = |kind: usize, index: usize| Origin {
            source: self.feature,
            role: self.roles[kind][index].clone(),
        };
        // The faces around every edge and vertex.
        let mut edge_faces = vec![std::collections::BTreeSet::new(); boundary.edges.len()];
        let mut vertex_faces = vec![std::collections::BTreeSet::new(); boundary.vertices.len()];
        for (index, face) in boundary.faces.iter().enumerate() {
            for trim in &face.loops {
                if let Some(vertex) = trim.vertex {
                    vertex_faces[vertex].insert(origin(0, index));
                }
                for &(edge, _) in &trim.uses {
                    edge_faces[edge].insert(origin(0, index));
                    for vertex in boundary.edges[edge].vertices {
                        vertex_faces[vertex].insert(origin(0, index));
                    }
                }
            }
        }
        let frame = |frame: &brep::Frame, flip: bool| ReferenceFrame {
            origin_mm: frame.origin,
            u: frame.x,
            v: if flip { frame.y.map(|x| -x) } else { frame.y },
        };
        let mut elements = Vec::new();
        for (index, face) in boundary.faces.iter().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: index as u32,
                },
                origins: [origin(0, index)].into(),
                boundary: Default::default(),
                frames: vec![],
                geometry: match &face.surface {
                    Surface::Plane(plane) => Some(AnalyticGeometry::Plane {
                        frame: frame(plane, !face.same_sense),
                    }),
                    _ => None,
                },
            });
        }
        for (index, edge) in boundary.edges.iter().enumerate() {
            let [start, end] = edge.range;
            let geometry = match &edge.curve {
                Curve::Line { .. } => {
                    let [a_mm, b_mm] = edge.vertices.map(|vertex| boundary.vertices[vertex]);
                    Some(AnalyticGeometry::Line { a_mm, b_mm })
                }
                Curve::Circle {
                    frame: placement,
                    radius,
                } => Some(AnalyticGeometry::CircularArc {
                    frame: frame(placement, false),
                    radius_mm: *radius,
                    start_rad: start,
                    sweep_rad: end - start,
                }),
                Curve::Ellipse {
                    frame: placement,
                    semiaxes,
                } => Some(AnalyticGeometry::EllipticArc {
                    frame: frame(placement, false),
                    semiaxes_mm: *semiaxes,
                    start_rad: start,
                    sweep_rad: end - start,
                }),
                Curve::Hyperbola { .. } | Curve::Parabola { .. } | Curve::Nurbs { .. } => None,
            };
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: index as u32,
                },
                origins: edge_faces[index]
                    .iter()
                    .cloned()
                    .chain([origin(1, index)])
                    .collect(),
                boundary: edge_faces[index].clone(),
                frames: vec![],
                geometry,
            });
        }
        for (index, &position_mm) in boundary.vertices.iter().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Vertex,
                    ordinal: index as u32,
                },
                origins: vertex_faces[index]
                    .iter()
                    .cloned()
                    .chain([origin(2, index)])
                    .collect(),
                boundary: vertex_faces[index].clone(),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Point { position_mm }),
            });
        }
        TopologySnapshot {
            owner,
            revision,
            elements,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const FRUSTUM: &str =
        include_str!("../../../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step");
    const VOID: &str =
        include_str!("../../../multicad-app/tests/fixtures/exchange/corpus/occt-box-void.step");

    fn frustum() -> Solid {
        Solid::import_rust_brep(FRUSTUM, EntityId::new(), GeometryPolicy::default()).unwrap()
    }

    #[test]
    fn imported_frustum_reports_its_retained_measurements() {
        let solid = frustum();
        assert!(solid.imported_brep_body().is_some());
        let exact = PI * 10. / 3. * 112.;
        assert!((solid.volume_mm3().unwrap() - exact).abs() < 1e-9 * exact);
        assert_eq!(solid.solid_count().unwrap(), 1);
        let faces = solid.faces().unwrap();
        assert_eq!(faces.len(), 3);
        assert_eq!(faces.iter().filter(|face| face.plane.is_some()).count(), 2);
        let bounds = solid.bounds().unwrap();
        assert!((bounds.max[2] - 10.).abs() < 1e-9 && (bounds.min[0] + 8.).abs() < 1e-9);
        let receipt = solid
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
        assert!(receipt.edge_pairs > 0);
        let topology = solid.rust_topology().unwrap().unwrap();
        topology.validate().unwrap();
        assert_eq!(
            topology
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Face)
                .count(),
            3
        );
        assert!(!solid.display_edges(0.05).unwrap().is_empty());
        assert!(solid.supports_step_export());
    }

    #[test]
    fn placement_keeps_volume_and_references_and_moves_the_boundary() {
        let solid = frustum();
        let transform = RigidTransform {
            rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            translation_mm: [5., 0., 2.],
        };
        let moved = solid.transformed(transform).unwrap();
        assert_eq!(moved.volume_mm3().unwrap(), solid.volume_mm3().unwrap());
        let bounds = moved.bounds().unwrap();
        assert!((bounds.min[2] - 2.).abs() < 1e-9 && (bounds.max[0] - 13.).abs() < 1e-9);
        let roles = |solid: &Solid| {
            let topology = solid.rust_topology().unwrap().unwrap();
            topology
                .elements
                .into_iter()
                .map(|e| e.origins)
                .collect::<Vec<_>>()
        };
        assert_eq!(roles(&moved), roles(&solid));
    }

    #[test]
    fn export_reads_back_and_downstream_operations_refuse() {
        let solid = frustum();
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("frustum.step");
        solid.write_step(&path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let back =
            Solid::import_rust_step(&text, EntityId::new(), GeometryPolicy::default()).unwrap();
        assert!((back.volume_mm3().unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-9);
        let other = Solid::box_solid([4., 4., 4.], GeometryPolicy::default()).unwrap();
        for operation in [
            BooleanOperation::Union,
            BooleanOperation::Subtract,
            BooleanOperation::Intersect,
        ] {
            let refused = solid.boolean(&other, operation, GeometryPolicy::default());
            assert!(refused.is_err(), "{operation:?} must refuse");
        }
        let plane = PlaneFrame {
            origin: [0., 0., 5.],
            ..Default::default()
        };
        assert!(
            solid
                .split_by_plane(plane, GeometryPolicy::default())
                .is_err()
        );
    }

    #[test]
    fn written_boundaries_verify_and_altered_ones_do_not() {
        let solid = frustum();
        let mut writer = multicad_rust_kernel::step::Writer::new();
        writer
            .append_brep(&solid.imported_brep_body().unwrap().solid)
            .unwrap();
        let text = writer.finish().unwrap();
        let volume = solid.verify_imported_brep_step(&text).unwrap().unwrap();
        assert!((volume - solid.volume_mm3().unwrap()).abs() < 1e-9);
        let other =
            Solid::import_rust_brep(VOID, EntityId::new(), GeometryPolicy::default()).unwrap();
        assert!(other.verify_imported_brep_step(&text).unwrap().is_err());
        let box_solid = Solid::box_solid([4., 4., 4.], GeometryPolicy::default()).unwrap();
        assert!(box_solid.verify_imported_brep_step(&text).is_none());
    }

    #[test]
    fn voids_import_and_refusals_keep_their_kind() {
        let hollow =
            Solid::import_rust_brep(VOID, EntityId::new(), GeometryPolicy::default()).unwrap();
        assert!((hollow.volume_mm3().unwrap() - 5760.).abs() < 1e-9);
        let open = FRUSTUM.replacen(
            "CLOSED_SHELL('',(#17,#54,#109))",
            "CLOSED_SHELL('',(#17,#54))",
            1,
        );
        let refused = Solid::import_rust_brep(&open, EntityId::new(), GeometryPolicy::default());
        assert_eq!(refused.err().unwrap().kind, GeometryErrorKind::Unsupported);
        assert_eq!(
            failure(brep::Error::Cancelled).kind,
            GeometryErrorKind::Cancelled
        );
        assert_eq!(failure(brep::Error::Limit).kind, GeometryErrorKind::Limit);
    }
}
