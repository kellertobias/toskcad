//! An imported solid that keeps its exchange boundary: support surfaces,
//! trim loops, edge curves and vertices, exactly as the file states them.
//! Volume, bounds, meshes, sections and STEP export read that boundary.
use super::*;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, ReferenceFrame, TopologyElement, TopologyKey, TopologyKind,
    TopologySnapshot,
};
use multicad_rust_kernel::brep::{self, Curve, Surface};
use multicad_rust_kernel::display::{MeshError, Tessellation};
use multicad_rust_kernel::step_read;
use std::collections::BTreeMap;

mod affine;
mod coped_tube;
mod rounded_box;
mod transverse_bore;

pub use coped_tube::{CopedTube, CopedTubeNames};
pub(crate) use transverse_bore::TransverseRecord;
pub use transverse_bore::{BoreStep, TransverseBoreTool};

const OPERATION: &str = "imported boundary";

/// The most topology elements one retained solid may publish.
const MAX_ELEMENTS: usize = 100_000;

pub(crate) struct ImportedBrepBody {
    pub(super) solid: brep::Solid,
    pub(super) feature: EntityId,
    names: Names,
    /// Parity tessellations of this boundary by deflection, built on use.
    parity: std::cell::RefCell<Vec<(f64, Rc<brep::mesh::Mesh>)>>,
    /// A proved outside cylinder: the material lies outside it.
    excluded: Option<coped_tube::ExcludedCylinder>,
    /// The original prism and stepped tool of a transverse bore.
    pub(super) transverse: Option<Rc<transverse_bore::TransverseRecord>>,
}

/// How a retained boundary names its topology.
enum Names {
    /// Stable role of every face, edge and vertex, taken from the source
    /// geometry so that a placed copy keeps its references.
    Imported(Rc<[Vec<String>; 3]>),
    /// Authored face identities of a boundary this kernel constructed.
    Authored(Rc<super::authored_sweep::AuthoredNaming>),
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

/// The product structure of an imported STEP assembly: every definition
/// once, in its own coordinates, and the placement of its roots.
pub(crate) struct ImportedAssembly {
    pub(super) structure: Rc<step_read::brep::Assembly<brep::Solid>>,
    pub(super) placement: step_read::brep::Placement,
}

impl ImportedAssembly {
    pub(super) fn transformed(&self, transform: RigidTransform) -> Self {
        let motion = step_read::brep::Placement {
            origin: transform.translation_mm,
            axes: std::array::from_fn(|j| std::array::from_fn(|i| transform.rotation[i][j])),
        };
        Self {
            structure: Rc::clone(&self.structure),
            placement: motion.then(&self.placement),
        }
    }
}

/// One placed solid of an imported STEP assembly.
#[derive(Clone, Debug, PartialEq)]
pub struct AssemblyOccurrence {
    /// The root product's name, then the name of every instance down to
    /// this occurrence.
    pub path: Vec<String>,
    /// The shared definition this occurrence places, and its name.
    pub definition: usize,
    pub definition_name: String,
    /// The solid's index within its definition.
    pub solid: usize,
    /// Where the definition's coordinates lie in the world.
    pub origin_mm: [f64; 3],
    pub axes: [[f64; 3]; 3],
    pub volume_mm3: f64,
    pub bounds: Bounds,
}

/// A stable identity for one occurrence: the import, the path of instance
/// indices from its root and the solid's index in its definition.
fn occurrence_feature(feature: EntityId, occurrence: &step_read::brep::Occurrence) -> EntityId {
    let mut bytes = b"multicad.step-occurrence/1".to_vec();
    bytes.extend(feature.uuid_bits().to_le_bytes());
    bytes.extend((occurrence.path.len() as u64).to_le_bytes());
    for &index in occurrence.path.iter().chain([&occurrence.solid]) {
        bytes.extend((index as u64).to_le_bytes());
    }
    let hash = multicad_core::project_file::ContentHash::of(&bytes);
    EntityId::from_uuid_bits(u128::from_le_bytes(hash.0[..16].try_into().unwrap()))
}

impl Solid {
    /// Import every solid of a STEP file by retaining its boundary. Several
    /// solids of a flat file form one compound in file order; an assembly
    /// becomes one compound of its placed occurrences, depth first, that
    /// keeps its product structure.
    pub(super) fn import_rust_brep(
        input: &str,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let verified = super::step_import::pcurve_validation::validate(input, policy)?;
        Self::import_rust_brep_validated(&verified, feature, policy)
    }

    pub(super) fn import_rust_brep_validated(
        verified: &step_read::brep::PcurveValidatedInput<'_>,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if verified.tolerance() != policy.linear_tolerance_mm {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Input,
                "STEP trim verification tolerance does not match the import policy",
            ));
        }
        let mut assembly =
            step_read::brep::read_assembly_validated(verified, None).map_err(failure)?;
        if assembly.nodes.len() == 1 && assembly.nodes[0].children.is_empty() {
            return Self::import_flat_brep(assembly.nodes.remove(0).solids, feature, policy);
        }
        let roles: Vec<Vec<_>> = assembly
            .nodes
            .iter()
            .map(|node| {
                node.solids
                    .iter()
                    .map(|solid| Rc::new(roles(solid)))
                    .collect()
            })
            .collect();
        let mut members = Vec::new();
        for occurrence in assembly.occurrences() {
            members.push(Self::from_imported_brep(
                assembly.place(&occurrence).map_err(failure)?,
                occurrence_feature(feature, &occurrence),
                Rc::clone(&roles[occurrence.node][occurrence.solid]),
                policy,
            )?);
        }
        let members: Vec<_> = members.iter().collect();
        let compound = Self::make_rust_compound(&members, feature, policy, None)?;
        compound.with_imported_assembly(ImportedAssembly {
            structure: Rc::new(assembly),
            placement: step_read::brep::Placement::IDENTITY,
        })
    }

    fn import_flat_brep(
        solids: Vec<brep::Solid>,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
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
            let roles = Rc::new(roles(&solid));
            members.push(Self::from_imported_brep(solid, child, roles, policy)?);
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
        roles: Rc<[Vec<String>; 3]>,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        Self::retained(solid, feature, Names::Imported(roles), policy)
    }

    /// Retain a boundary this kernel assembled from exact exchange
    /// boundaries, naming its elements from their geometry alone.
    pub(super) fn from_assembled_brep(
        solid: brep::Solid,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let roles = Rc::new(roles(&solid));
        Self::from_imported_brep(solid, feature, roles, policy)
    }

    /// Retain a boundary this kernel constructed, with authored face names.
    pub(super) fn from_retained_brep(
        solid: brep::Solid,
        feature: EntityId,
        naming: Rc<super::authored_sweep::AuthoredNaming>,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if naming.faces.len() != solid.boundary().faces.len() {
            return Err(error(
                OPERATION,
                GeometryErrorKind::InvalidResult,
                "Authored names must cover every retained face",
            ));
        }
        Self::retained(solid, feature, Names::Authored(naming), policy)
    }

    fn retained(
        solid: brep::Solid,
        feature: EntityId,
        names: Names,
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
        Ok(Self {
            body: SolidBody::ImportedBrep(Rc::new(ImportedBrepBody {
                solid,
                feature,
                names,
                parity: Default::default(),
                excluded: None,
                transverse: None,
            })),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }

    /// The placed solids of an imported STEP assembly in member order, or
    /// `None` for any other solid.
    pub fn assembly_occurrences(&self) -> Option<Vec<AssemblyOccurrence>> {
        let compound = self.rust_compound()?;
        let assembly = compound.imported_assembly()?;
        let structure = &assembly.structure;
        let occurrences = structure.occurrences();
        Some(
            occurrences
                .iter()
                .zip(compound.members())
                .map(|(occurrence, member)| {
                    let mut path = vec![
                        structure.nodes[structure.roots[occurrence.path[0]]]
                            .name
                            .clone(),
                    ];
                    let mut node = structure.roots[occurrence.path[0]];
                    for &index in &occurrence.path[1..] {
                        let child = &structure.nodes[node].children[index];
                        path.push(child.name.clone());
                        node = child.node;
                    }
                    let placement = assembly.placement.then(&occurrence.placement);
                    let body = &member.imported_brep_body().expect("imported member").solid;
                    let [min, max] = body.bounds_mm();
                    AssemblyOccurrence {
                        path,
                        definition: occurrence.node,
                        definition_name: structure.nodes[occurrence.node].name.clone(),
                        solid: occurrence.solid,
                        origin_mm: placement.origin,
                        axes: placement.axes,
                        volume_mm3: body.volume_mm3(),
                        bounds: Bounds { min, max },
                    }
                })
                .collect(),
        )
    }

    /// The STEP text of an imported assembly with its product structure, or
    /// `None` for any other solid.
    pub(super) fn imported_assembly_step(
        &self,
    ) -> Option<Result<String, multicad_rust_kernel::step::StepError>> {
        let assembly = self.rust_compound()?.imported_assembly()?;
        Some(multicad_rust_kernel::step::write_assembly(
            &assembly.structure,
            assembly.placement,
        ))
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
            names: match &self.names {
                Names::Imported(roles) => Names::Imported(Rc::clone(roles)),
                Names::Authored(naming) => Names::Authored(Rc::new(naming.moved(transform))),
            },
            parity: Default::default(),
            // A rigid motion carries the outside cylinder with the material.
            excluded: self.excluded.map(|cylinder| coped_tube::ExcludedCylinder {
                point_mm: map(cylinder.point_mm),
                axis: linear(cylinder.axis),
                radius_mm: cylinder.radius_mm,
            }),
            transverse: match &self.transverse {
                Some(record) => Some(Rc::new(record.transformed(transform)?)),
                None => None,
            },
        })
    }

    /// Classify a whole world box against this boundary: exact clearance
    /// from every face support, then tessellation parity at a deflection
    /// below that clearance.
    /// `deviation_mm` widens the margin by how far this boundary may lie
    /// from the exact source it stands for.
    pub(super) fn classify_box(
        &self,
        bounds: [[f64; 3]; 2],
        deviation_mm: f64,
    ) -> Result<brep::BoxSide, GeometryError> {
        let budget = |stage: &'static str| {
            move |cause: brep::Error| match cause {
                brep::Error::Limit => error(
                    "Boolean",
                    GeometryErrorKind::Limit,
                    format!("Rust imported-boundary containment exceeded its {stage} budget"),
                ),
                other => failure(other),
            }
        };
        let margin = self.solid.tolerance() + deviation_mm;
        let [low, high] = self.solid.bounds_mm();
        let [lo, hi] = bounds;
        if (0..3).any(|i| hi[i] < low[i] - margin || lo[i] > high[i] + margin) {
            return Ok(brep::BoxSide::Outside);
        }
        let Some(clearance) = self
            .solid
            .box_clearance(bounds, margin, None)
            .map_err(budget("hull"))?
        else {
            return Ok(brep::BoxSide::Unresolved);
        };
        // Parity comes from the coarsest of four halving tessellations whose
        // deflection the clearance clears four times over. A box nearer the
        // boundary first walks its centre away within certified clearance
        // balls, which never cross the boundary.
        let diagonal = (0..3)
            .map(|i| high[i] - low[i])
            .fold(0., f64::portable_hypot);
        let coarse = (diagonal * 1.25e-4).max(10. * margin);
        let levels = [1., 0.5, 0.25, 0.125].map(|f| (coarse * f).max(10. * margin));
        let (point, free) = if clearance >= 4. * coarse {
            (std::array::from_fn(|i| 0.5 * (lo[i] + hi[i])), clearance)
        } else {
            let center = std::array::from_fn(|i| 0.5 * (lo[i] + hi[i]));
            self.solid
                .escape(center, clearance, 4. * coarse, margin, None)
                .map_err(budget("escape"))?
        };
        let Some(&deflection) = levels.iter().find(|&&d| free >= 4. * d) else {
            return Ok(brep::BoxSide::Unresolved);
        };
        let cached = self
            .parity
            .borrow()
            .iter()
            .find(|(d, _)| *d == deflection)
            .map(|(_, mesh)| Rc::clone(mesh));
        let mesh = match cached {
            Some(mesh) => mesh,
            None => {
                let mesh = Rc::new(
                    brep::mesh::tessellate(
                        &self.solid,
                        brep::mesh::Options {
                            deflection_mm: deflection,
                            angle_rad: 0.5,
                            triangle_limit: 4_000_000,
                        },
                        None,
                    )
                    .map_err(budget("tessellation"))?,
                );
                self.parity
                    .borrow_mut()
                    .push((deflection, Rc::clone(&mesh)));
                mesh
            }
        };
        // The box centre's clearance is at least the box's own.
        Ok(self
            .solid
            .box_parity([point, point], free, &mesh, deflection))
    }

    /// Authored face names, when this kernel constructed the boundary.
    pub(super) fn authored_naming(&self) -> Option<&super::authored_sweep::AuthoredNaming> {
        match &self.names {
            Names::Authored(naming) => Some(naming),
            Names::Imported(_) => None,
        }
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

    /// The loops of a planar face as (edge ordinal, reversed) uses, the
    /// outer loop first: the loop whose bounds contain every other loop's.
    pub(crate) fn planar_loops(&self, face: u32) -> Vec<Vec<(u32, bool)>> {
        let boundary = self.solid.boundary();
        let Some(face) = boundary.faces.get(face as usize) else {
            return vec![];
        };
        if !matches!(face.surface, brep::Surface::Plane(_)) {
            return vec![];
        }
        let extent = |uses: &[(usize, bool)]| -> [[f64; 3]; 2] {
            let mut bounds = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
            for &(e, _) in uses {
                let edge = &boundary.edges[e];
                let [lows, highs] = &mut bounds;
                for (axis, (low, high)) in lows.iter_mut().zip(highs.iter_mut()).enumerate() {
                    let [a, b] = edge.curve.extent(edge.range[0], edge.range[1], axis);
                    *low = low.min(a);
                    *high = high.max(b);
                }
            }
            bounds
        };
        let mut loops: Vec<_> = face
            .loops
            .iter()
            .map(|l| {
                let bounds = extent(&l.uses);
                let size: f64 = (0..3).map(|i| bounds[1][i] - bounds[0][i]).sum();
                (
                    size,
                    l.uses
                        .iter()
                        .map(|&(e, forward)| (e as u32, !forward))
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        loops.sort_by(|a, b| b.0.total_cmp(&a.0));
        loops.into_iter().map(|(_, uses)| uses).collect()
    }

    pub(crate) fn topology(&self, owner: EntityId, revision: u64) -> TopologySnapshot {
        let boundary = self.solid.boundary();
        // Every name of a face, edge (kind 1) or vertex (kind 2). Authored
        // edges and vertices take the names of the faces they bound.
        let names = |kind: usize, index: usize| -> Vec<Origin> {
            match &self.names {
                Names::Imported(roles) => vec![Origin {
                    source: self.feature,
                    role: roles[kind][index].clone(),
                }],
                Names::Authored(naming) if kind == 0 => naming.faces[index].origins.clone(),
                Names::Authored(_) => vec![],
            }
        };
        // The faces around every edge and vertex.
        let mut edge_faces = vec![std::collections::BTreeSet::new(); boundary.edges.len()];
        let mut vertex_faces = vec![std::collections::BTreeSet::new(); boundary.vertices.len()];
        for (index, face) in boundary.faces.iter().enumerate() {
            for trim in &face.loops {
                if let Some(vertex) = trim.vertex {
                    vertex_faces[vertex].extend(names(0, index));
                }
                for &(edge, _) in &trim.uses {
                    edge_faces[edge].extend(names(0, index));
                    for vertex in boundary.edges[edge].vertices {
                        vertex_faces[vertex].extend(names(0, index));
                    }
                }
            }
        }
        let authored = self.authored_naming();
        let mut edge_owners = vec![Vec::new(); boundary.edges.len()];
        for (index, face) in boundary.faces.iter().enumerate() {
            for trim in &face.loops {
                for &(edge, _) in &trim.uses {
                    edge_owners[edge].push(index);
                }
            }
        }
        // The names of the faces across each edge of a face.
        let neighbors = |index: usize| {
            let mut around = std::collections::BTreeSet::new();
            for trim in &boundary.faces[index].loops {
                for &(edge, _) in &trim.uses {
                    for &other in &edge_owners[edge] {
                        if other != index {
                            around.extend(names(0, other));
                        }
                    }
                }
            }
            around
        };
        let frame = |frame: &brep::Frame, flip: bool| ReferenceFrame {
            origin_mm: frame.origin,
            u: frame.x,
            v: if flip { frame.y.map(|x| -x) } else { frame.y },
        };
        let mut elements = Vec::new();
        for (index, face) in boundary.faces.iter().enumerate() {
            let geometry = match &face.surface {
                Surface::Plane(plane) => Some(AnalyticGeometry::Plane {
                    frame: frame(plane, !face.same_sense),
                }),
                _ => None,
            };
            let element = match authored {
                None => TopologyElement {
                    key: TopologyKey {
                        kind: TopologyKind::Face,
                        ordinal: index as u32,
                    },
                    origins: names(0, index).into_iter().collect(),
                    boundary: Default::default(),
                    frames: vec![],
                    geometry,
                },
                Some(naming) => {
                    let name = &naming.faces[index];
                    TopologyElement {
                        key: TopologyKey {
                            kind: TopologyKind::Face,
                            ordinal: index as u32,
                        },
                        origins: name.origins.iter().cloned().collect(),
                        boundary: neighbors(index),
                        frames: name
                            .anchor
                            .zip(name.origins.first())
                            .map(|(frame, origin)| AnchoredFrame {
                                origin: origin.clone(),
                                frame,
                            })
                            .into_iter()
                            .collect(),
                        geometry: name
                            .support
                            .map(|frame| AnalyticGeometry::Plane { frame })
                            .or(geometry),
                    }
                }
            };
            elements.push(element);
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
                    .chain(names(1, index))
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
                    .chain(names(2, index))
                    .collect(),
                boundary: vertex_faces[index].clone(),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Point { position_mm }),
            });
        }
        TopologySnapshot {
            owner,
            revision,
            elements: elements.into(),
        }
    }
}

impl Solid {
    /// A finite rectangular or polygon prism wholly inside or wholly outside
    /// a retained imported boundary reuses one unchanged operand.
    pub(in crate::geometry) fn contained_imported_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        use super::helical_body::contained_boolean::{probe_bounds, whole_operand_result};
        use multicad_rust_kernel::helical_thread::BoxClassification;
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        // A retained boundary in world coordinates, or a bevel's verified
        // encoding in its local frame with the probe box charted into it.
        let retained = |solid: &Self| -> Option<(Rc<Self>, PlaneFrame)> {
            let (body, frame) = solid.bevel_body()?;
            Some((body.retained(policy)?, frame))
        };
        let (probe, source_is_left, proxy) = match (&self.body, &other.body) {
            (SolidBody::ImportedBrep(_), _) => (other, true, None),
            (_, SolidBody::ImportedBrep(_)) => (self, false, None),
            (SolidBody::Bevel(..), _) => (other, true, Some(retained(self)?)),
            (_, SolidBody::Bevel(..)) => (self, false, Some(retained(other)?)),
            _ => return None,
        };
        let world = probe_bounds(probe)?;
        let source = if source_is_left { self } else { other };
        let (boundary, bounds, deviation) = match &proxy {
            None => (source, [world.map(|b| b.lo), world.map(|b| b.hi)], 0.),
            Some((retained, frame)) => {
                let local = super::helical_body::contained_boolean::chart_bounds(
                    world,
                    frame.origin,
                    [frame.u, frame.v, frame.normal()],
                );
                // The verified encoding lies within the tolerance of the
                // exact bevel.
                (
                    retained.as_ref(),
                    [local.map(|[lo, _]| lo), local.map(|[_, hi]| hi)],
                    policy.linear_tolerance_mm,
                )
            }
        };
        let SolidBody::ImportedBrep(body) = &boundary.body else {
            return None;
        };
        let relation = match body.classify_box(bounds, deviation) {
            Ok(brep::BoxSide::Inside) => BoxClassification::Inside,
            Ok(brep::BoxSide::Outside) => BoxClassification::Outside,
            Ok(brep::BoxSide::Unresolved) => return None,
            Err(error) if error.kind == GeometryErrorKind::Limit => return Some(Err(error)),
            Err(_) => return None,
        };
        whole_operand_result(self, other, source_is_left, relation, operation)
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
    fn an_imported_assembly_moves_and_exports_with_its_structure() {
        const TWO: &str = include_str!(
            "../../../multicad-app/tests/fixtures/exchange/corpus/occt-two-bodies.step"
        );
        let policy = GeometryPolicy::default();
        let solid = Solid::import_rust_brep(TWO, EntityId::new(), policy).unwrap();
        let transform = RigidTransform {
            rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            translation_mm: [5., 0., 2.],
        };
        let moved = solid.transformed(transform).unwrap();
        let occurrences = moved.assembly_occurrences().unwrap();
        assert_eq!(occurrences.len(), 2);
        assert_eq!(occurrences[0].origin_mm, [5., 0., 2.]);
        assert_eq!(occurrences[0].axes[0], [0., 1., 0.]);
        // The moved structure is written under one product that places it.
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("moved.step");
        moved.write_step(&path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let back = Solid::import_rust_step(&text, EntityId::new(), policy).unwrap();
        let read = back.assembly_occurrences().unwrap();
        assert_eq!(read.len(), 2);
        assert_eq!(read[1].path.len(), 3);
        for (a, b) in read.iter().zip(&occurrences) {
            for (x, y) in [a.bounds.min, a.bounds.max]
                .iter()
                .flatten()
                .zip([b.bounds.min, b.bounds.max].iter().flatten())
            {
                assert!((x - y).abs() < 1e-9);
            }
        }
        // Only an unshared compound takes a structure.
        let structure = || ImportedAssembly {
            structure: Rc::new(step_read::brep::read_assembly(TWO, 1e-7, None).unwrap()),
            placement: step_read::brep::Placement::IDENTITY,
        };
        assert!(frustum().assembly_occurrences().is_none());
        assert!(frustum().imported_assembly_step().is_none());
        assert!(frustum().with_imported_assembly(structure()).is_err());
        let shared = solid.share().unwrap();
        assert!(solid.with_imported_assembly(structure()).is_err());
        drop(shared);
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
