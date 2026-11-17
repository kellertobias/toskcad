//! Actual retained STEP C11 geometry. No authored screw recipe is reconstructed.
use super::*;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, TopologyElement, TopologyKey, TopologyKind,
};
use multicad_rust_kernel::{
    display::{MeshError, Tessellation},
    step_read::npt_stock as read,
};
use std::{
    cell::{Cell, RefCell},
    collections::BTreeSet,
};
mod placement;
mod step;
#[derive(Clone, Copy)]
struct OriginalFace {
    area: [f64; 2],
    center: [[f64; 2]; 3],
}
const BYTES: usize = 64 * 1024 * 1024;
const OVERHEAD: usize = 1024 * 1024;
#[derive(Clone, Copy)]
struct MeshAccounting {
    work: usize,
    actual: usize,
    peak: usize,
}
struct BuiltMesh {
    mesh: Mesh,
    proof: multicad_rust_kernel::manufacturing::MeshProof,
    accounting: MeshAccounting,
}
struct CachedDisplay {
    key: [u64; 3],
    mesh: Mesh,
    work: usize,
    peak: usize,
}
/// Aggregate mesh work includes its immutable proof prerequisites ONCE;
/// actual work distinguishes a fresh query from a cache-copy request. Do not add the
/// preparation field again when aggregating this receipt. Document bytes
/// describe the original neutral parse scope, not current resident storage.
pub struct LiteralNptRuntimeEvidence {
    pub source_query_work: usize,
    pub preparation_work: usize,
    pub neutral_document_budget_bytes: usize,
    pub neutral_document_work: usize,
    pub resident_budget_bytes: usize,
    pub last_mesh_work: Option<usize>,
    pub last_mesh_actual_work: Option<usize>,
    pub last_mesh_peak_bytes: Option<usize>,
    /// Successfully published fresh display/manufacturing results. Refused
    /// requests preserve this count and the last successful receipt.
    pub mesh_builds: usize,
}
pub(crate) struct LiteralNptBody {
    pub(super) feature: EntityId,
    source: read::OwnedMaterial,
    plan: Rc<read::LiteralGeometryPlan>,
    occurrence: Option<read::LiteralOccurrencePlacement>,
    vertices: Vec<[f64; 3]>,
    policy: GeometryPolicy,
    faces: Rc<RefCell<Option<Vec<OriginalFace>>>>,
    display: RefCell<Option<CachedDisplay>>,
    last_mesh: Cell<Option<(usize, usize)>>,
    last_actual: Cell<Option<usize>>,
    mesh_builds: Cell<usize>,
    step_ready: Cell<Option<(usize, bool)>>,
    last_step: Cell<Option<multicad_rust_kernel::step::literal_npt::Receipt>>,
}
pub(in crate::geometry) fn source_signature(
    source: &read::OwnedMaterial,
) -> Result<multicad_core::project_file::ContentHash, GeometryError> {
    // Original numeric leaves and their actual incidence identify geometry.
    // Different control assignments do not become the same part merely
    // because their point tables coincide. Original record IDs are excluded.
    // This pre-existing identity stage has its own 64 MiB arena, separate
    // from manufacturing/query work. Mark physical point ownership once;
    // support-chart anchors do not become physical just by being in the file.
    let (physical, physical_count) = identity_physical_points(source)?;
    let mut bytes = Vec::new();
    let limit = BYTES
        .checked_sub(source.retained_budget_bytes())
        .and_then(|n| n.checked_sub(physical.capacity() * std::mem::size_of::<usize>()))
        .ok_or_else(|| failure(read::Error::Limit))?;
    let mut append = |value: &[u8]| -> Result<(), GeometryError> {
        let end = bytes
            .len()
            .checked_add(value.len())
            .filter(|n| *n <= limit)
            .ok_or_else(|| failure(read::Error::Limit))?;
        if end > bytes.capacity() {
            let additional = (end - bytes.len()).max(4096).min(limit - bytes.len());
            bytes
                .try_reserve_exact(additional)
                .map_err(|_| failure(read::Error::Limit))?;
        }
        if bytes.capacity() > limit {
            return Err(failure(read::Error::Limit));
        }
        bytes.extend_from_slice(value);
        Ok(())
    };
    append(b"multicad.literal-npt.identity/3")?;
    append(&(physical_count as u64).to_le_bytes())?;
    for (index, ordinal) in physical.iter().enumerate() {
        if *ordinal == usize::MAX {
            continue;
        }
        let coordinates = source
            .point(index)
            .expect("retained physical point owner")
            .coordinates();
        append(&(coordinates.len() as u64).to_le_bytes())?;
        for scalar in coordinates {
            append(&(scalar.spelling().len() as u64).to_le_bytes())?;
            append(scalar.spelling().as_bytes())?;
        }
    }
    append(b"vertices")?;
    append(&(source.vertices().len() as u64).to_le_bytes())?;
    for vertex in source.vertices() {
        append(&(physical[vertex.point_index()] as u64).to_le_bytes())?;
    }
    append(b"curves/degree11")?;
    append(&(source.curves().len() as u64).to_le_bytes())?;
    for curve in source.curves() {
        append(&(curve.spans() as u64).to_le_bytes())?;
        append(&(curve.point_indices().len() as u64).to_le_bytes())?;
        for &p in curve.point_indices() {
            append(&(physical[p] as u64).to_le_bytes())?;
        }
    }
    append(b"tensors/degree11/11")?;
    append(&(source.tensors().len() as u64).to_le_bytes())?;
    for tensor in source.tensors() {
        for n in tensor.spans() {
            append(&(n as u64).to_le_bytes())?;
        }
        append(&(tensor.point_indices().len() as u64).to_le_bytes())?;
        for &p in tensor.point_indices() {
            append(&(physical[p] as u64).to_le_bytes())?;
        }
    }
    append(b"edges")?;
    append(&(source.edges().len() as u64).to_le_bytes())?;
    for edge in source.edges() {
        append(&(edge.curve_index() as u64).to_le_bytes())?;
        for v in edge.vertex_indices() {
            append(&(v as u64).to_le_bytes())?;
        }
    }
    append(b"faces")?;
    append(&(source.faces().len() as u64).to_le_bytes())?;
    for face in source.faces() {
        match face.support() {
            read::Support::Tensor(index) => {
                append(&[0, u8::from(face.forward())])?;
                append(&(index as u64).to_le_bytes())?;
            }
            read::Support::Plane { height } => {
                append(&[1, u8::from(face.forward())])?;
                let z = source
                    .point(height)
                    .expect("retained plane owner")
                    .coordinates()[2];
                let spelling = if z.bounds() == [0., 0.] {
                    "0"
                } else {
                    z.spelling()
                };
                append(&(spelling.len() as u64).to_le_bytes())?;
                append(spelling.as_bytes())?;
            }
        }
        append(&(face.coedges().len() as u64).to_le_bytes())?;
        for coedge in face.coedges() {
            append(&(coedge.edge_index() as u64).to_le_bytes())?;
            append(&[u8::from(coedge.forward())])?;
        }
    }
    append(b"original-step-mapping")?;
    append(&(source.original_mapping_identity().len() as u64).to_le_bytes())?;
    append(source.original_mapping_identity())?;
    Ok(multicad_core::project_file::ContentHash::of(&bytes))
}

fn identity_physical_points(
    source: &read::OwnedMaterial,
) -> Result<(Vec<usize>, usize), GeometryError> {
    // Count logical visits before each bounded pass, not unbounded rescans per
    // point. Identity has no caller query ledger or cancellation API.
    let mut remaining = 16_000_000usize;
    let mut charge = |n: usize| -> Result<(), GeometryError> {
        remaining = remaining
            .checked_sub(n)
            .ok_or_else(|| failure(read::Error::Limit))?;
        Ok(())
    };
    charge(source.point_count())?;
    let available = BYTES
        .checked_sub(source.retained_budget_bytes())
        .ok_or_else(|| failure(read::Error::Limit))?;
    let requested = source
        .point_count()
        .checked_mul(std::mem::size_of::<usize>())
        .ok_or_else(|| failure(read::Error::Limit))?;
    if requested > available {
        return Err(failure(read::Error::Limit));
    }
    let mut physical = Vec::new();
    physical
        .try_reserve_exact(source.point_count())
        .map_err(|_| failure(read::Error::Limit))?;
    if physical
        .capacity()
        .checked_mul(std::mem::size_of::<usize>())
        .is_none_or(|n| n > available)
    {
        return Err(failure(read::Error::Limit));
    }
    physical.resize(source.point_count(), usize::MAX);
    charge(source.vertices().len())?;
    for vertex in source.vertices() {
        physical[vertex.point_index()] = 0;
    }
    charge(source.curves().len())?;
    for curve in source.curves() {
        charge(curve.point_indices().len())?;
        for index in curve.point_indices() {
            physical[*index] = 0;
        }
    }
    charge(source.tensors().len())?;
    for tensor in source.tensors() {
        charge(tensor.point_indices().len())?;
        for index in tensor.point_indices() {
            physical[*index] = 0;
        }
    }
    charge(physical.len())?;
    let mut count = 0;
    for ordinal in &mut physical {
        if *ordinal != usize::MAX {
            *ordinal = count;
            count += 1;
        }
    }
    Ok((physical, count))
}

fn failure(cause: read::Error) -> GeometryError {
    error(
        "literal NPT source",
        match cause {
            read::Error::Invalid => GeometryErrorKind::Input,
            read::Error::Unsupported => GeometryErrorKind::Unsupported,
            read::Error::Limit => GeometryErrorKind::Limit,
            read::Error::Cancelled => GeometryErrorKind::Cancelled,
        },
        format!("Retained original geometry query failed: {cause:?}"),
    )
}
fn quality(q: Tessellation) -> Result<(), MeshError> {
    if !q.deflection.is_finite()
        || q.deflection <= 0.
        || !(0.01..=1.).contains(&q.angle)
        || !(1..=5_000_000).contains(&q.triangle_limit)
    {
        Err(MeshError::InvalidOptions)
    } else {
        Ok(())
    }
}
impl Solid {
    pub fn from_literal_npt_source(
        source: read::OwnedMaterial,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        if feature.is_nil() {
            return Err(input("Literal import requires a source identity"));
        }
        // Original STEP placement is already part of these WORLD source
        // APIs. A NEW Native occurrence transform remains a separate gate.
        let b = source.world_enclosure_mm();
        policy.point(b.map(|v| v[0]))?;
        policy.point(b.map(|v| v[1]))?;
        let plan = source
            .geometry_plan(
                read::LiteralGeometryOptions {
                    maximum_bytes: BYTES - OVERHEAD,
                    ..Default::default()
                },
                interrupt,
            )
            .map_err(failure)?;
        let available = (BYTES - OVERHEAD)
            .checked_sub(plan.retained_budget_bytes() - source.retained_budget_bytes())
            .ok_or_else(|| failure(read::Error::Limit))?;
        let mut vertices = Vec::with_capacity(source.vertices().len());
        for index in 0..source.vertices().len() {
            let sample = source
                .sample_vertex(
                    index,
                    read::Options {
                        maximum_bytes: available,
                        ..Default::default()
                    },
                    interrupt,
                )
                .map_err(failure)?;
            let p = sample.position_mm();
            let error = sample
                .position_bounds_mm()
                .into_iter()
                .zip(p)
                .fold(0f64, |sum, (b, x)| {
                    (sum + (x - b[0]).abs().max((b[1] - x).abs()).next_up()).next_up()
                });
            if error > policy.linear_tolerance_mm / 8. {
                return Err(failure(read::Error::Unsupported));
            }
            policy.point(p)?;
            vertices.push(p);
        }
        Ok(Self {
            body: SolidBody::LiteralNpt(Rc::new(LiteralNptBody {
                feature,
                source,
                plan: Rc::new(plan),
                occurrence: None,
                vertices,
                policy,
                faces: Rc::new(RefCell::new(None)),
                display: RefCell::new(None),
                last_mesh: Cell::new(None),
                last_actual: Cell::new(None),
                mesh_builds: Cell::new(0),
                step_ready: Cell::new(None),
                last_step: Cell::new(None),
            })),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn literal_npt_body(&self) -> Option<&LiteralNptBody> {
        match &self.body {
            SolidBody::LiteralNpt(body) => Some(body),
            _ => None,
        }
    }
    /// Additional Native transform provenance; the source getter below keeps
    /// its original literal STEP geometry and original STEP operator unchanged.
    pub fn rust_literal_npt_occurrence(&self) -> Option<&read::LiteralOccurrencePlacement> {
        self.literal_npt_body()
            .and_then(|body| body.occurrence.as_ref())
    }
    pub fn rust_literal_npt_source(&self) -> Option<&read::OwnedMaterial> {
        self.literal_npt_body().map(|b| &b.source)
    }
    pub fn rust_literal_npt_runtime_evidence(&self) -> Option<LiteralNptRuntimeEvidence> {
        self.literal_npt_body()
            .map(|body| LiteralNptRuntimeEvidence {
                source_query_work: body.source.work_used(),
                preparation_work: body.plan.preparation_work_used()
                    + body.occurrence.as_ref().map_or(0, |p| p.work_used()),
                neutral_document_budget_bytes: body.source.original_document_budget_bytes(),
                neutral_document_work: body.source.original_document_work_used(),
                resident_budget_bytes: body.resident_bytes(),
                last_mesh_work: body.last_mesh.get().map(|v| v.0),
                last_mesh_actual_work: body.last_actual.get(),
                last_mesh_peak_bytes: body.last_mesh.get().map(|v| v.1),
                mesh_builds: body.mesh_builds.get(),
            })
    }
    pub fn rust_literal_npt_contains(&self, point: [f64; 3]) -> Result<bool, GeometryError> {
        self.policy.point(point)?;
        let body = self
            .literal_npt_body()
            .ok_or_else(|| failure(read::Error::Unsupported))?;
        let query = match &body.occurrence {
            Some(p) => p.inverse_box_mm(point.map(|v| [v, v])).map_err(failure)?,
            None => point.map(|v| [v, v]),
        };
        let result = body
            .source
            .classify_box(
                query,
                read::QueryOptions {
                    maximum_bytes: body.available_query_bytes()?,
                    ..Default::default()
                },
                None,
            )
            .map_err(failure)?;
        match result.location() {
            read::Location::Inside => Ok(true),
            read::Location::Outside => Ok(false),
            read::Location::Unresolved => Err(failure(read::Error::Unsupported)),
        }
    }
}
impl LiteralNptBody {
    fn available_query_bytes(&self) -> Result<usize, GeometryError> {
        (BYTES - OVERHEAD)
            .checked_sub(
                self.display
                    .borrow()
                    .as_ref()
                    .map_or(0, |c| Self::mesh_bytes(&c.mesh)),
            )
            .and_then(|v| {
                v.checked_sub(
                    self.occurrence
                        .as_ref()
                        .map_or(0, |p| p.retained_budget_bytes()),
                )
            })
            .and_then(|v| {
                v.checked_sub(
                    self.plan.retained_budget_bytes() - self.source.retained_budget_bytes(),
                )
            })
            .ok_or_else(|| failure(read::Error::Limit))
    }
    pub(super) fn volume_mm3(&self) -> f64 {
        let v = self.world_volume();
        v[0].midpoint(v[1])
    }
    pub(super) fn bounds(&self) -> Bounds {
        let b = self.world_bounds();
        Bounds {
            min: b.map(|v| v[0]),
            max: b.map(|v| v[1]),
        }
    }
    fn plane(&self, index: usize) -> Option<ReferenceFrame> {
        let owner = self.source.faces().get(index)?;
        let read::Support::Plane { height } = owner.support() else {
            return None;
        };
        // A plane's arbitrary in-plane chart origin is not a physical cap
        // anchor. Intersect the actual support with the earned thread axis.
        let z = self.source.point(height)?.coordinates()[2].bounds();
        let local = [[0., 0.], [0., 0.], z];
        let mut origin = self.source.map_box_mm(local).ok()?;
        let mut columns = self.source.linear_columns_bounds();
        if let Some(pose) = &self.occurrence {
            origin = pose.map_box_mm(origin).ok()?;
            for column in &mut columns {
                *column = pose.vector_bounds(*column).ok()?;
            }
        }
        let midpoint = |b: [f64; 2]| b[0].midpoint(b[1]);
        let o = origin.map(midpoint);
        let raw = [columns[0].map(midpoint), columns[1].map(midpoint)];
        let scales = raw.map(|a| a.into_iter().map(|x| x * x).sum::<f64>().sqrt());
        if scales.into_iter().any(|s| !s.is_finite() || s <= 0.) {
            return None;
        }
        // Normalization proposes a frame only. Full-cap correspondence below
        // independently accounts its divisions and scaled-product rounding.
        let unit = std::array::from_fn::<_, 2, _>(|j| raw[j].map(|x| x / scales[j]));
        let u = unit[0];
        let sign = if owner.forward() { 1. } else { -1. };
        let v = unit[1].map(|x| sign * x);
        let uncertainty = |ranges: [[f64; 2]; 3], values: [f64; 3]| {
            ranges.into_iter().zip(values).fold(0f64, |sum, (b, x)| {
                (sum + (x - b[0]).abs().max((b[1] - x).abs()).next_up()).next_up()
            })
        };
        let local_box = self.source.local_enclosure_mm();
        // Bound displacement from the actual plane/thread-axis intersection.
        // Directed subtraction includes its axial uncertainty before the
        // column approximation is multiplied.
        let radius = |axis: usize| {
            (local_box[axis][0] - local[axis][1])
                .next_down()
                .abs()
                .max((local_box[axis][1] - local[axis][0]).next_up().abs())
        };
        let mut error = uncertainty(origin, o);
        for (axis, column) in columns.iter().copied().enumerate().take(2) {
            let delta = column
                .into_iter()
                .zip(unit[axis])
                .fold(0f64, |sum, (b, x)| {
                    let product = x * scales[axis];
                    let lo = (b[0] - product.next_up()).next_down();
                    let hi = (b[1] - product.next_down()).next_up();
                    (sum + lo.abs().max(hi.abs())).next_up()
                });
            error = (error + (radius(axis) * delta).next_up()).next_up();
        }
        let dot = |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        if !error.is_finite()
            || error > self.policy.linear_tolerance_mm / 8.
            || (dot(u, u) - 1.).abs() > 1e-10
            || (dot(v, v) - 1.).abs() > 1e-10
            || dot(u, v).abs() > 1e-10
        {
            return None;
        }
        Some(ReferenceFrame { origin_mm: o, u, v })
    }
    pub(super) fn face_plane(&self, index: u32) -> Result<Option<super::FacePlane>, GeometryError> {
        let index = index as usize;
        let face = self
            .source
            .faces()
            .get(index)
            .ok_or_else(|| input("Face index is out of range"))?;
        if matches!(face.support(), read::Support::Tensor(_)) {
            return Ok(None);
        }
        self.plane(index)
            .map(|frame| Some((frame.origin_mm, frame.normal())))
            .ok_or_else(|| failure(read::Error::Unsupported))
    }
    fn mesh_bytes(mesh: &Mesh) -> usize {
        mesh.vertices.capacity() * 8 + (mesh.indices.capacity() + mesh.face_indices.capacity()) * 4
    }
    fn base_resident_bytes(&self) -> usize {
        self.plan.retained_budget_bytes()
            + self
                .occurrence
                .as_ref()
                .map_or(0, |p| p.retained_budget_bytes())
            + std::mem::size_of::<Self>()
            // Self Rc, shared-plan Rc (whose payload is already in the plan
            // receipt), and the face-cache Rc each retain two counters.
            + 6 * std::mem::size_of::<usize>()
            + std::mem::size_of::<RefCell<Option<Vec<OriginalFace>>>>()
            + self.vertices.capacity() * 24
            + self
                .faces
                .borrow()
                .as_ref()
                .map_or(0, |v| v.capacity() * std::mem::size_of::<OriginalFace>())
    }
    fn resident_bytes(&self) -> usize {
        self.base_resident_bytes()
            + self
                .display
                .borrow()
                .as_ref()
                .map_or(0, |c| Self::mesh_bytes(&c.mesh))
    }
    fn copy_display(
        &self,
        mesh: &Mesh,
        proof_work: usize,
        proof_peak: usize,
        actual: usize,
        other_mesh_bytes: usize,
    ) -> Result<(Mesh, MeshAccounting), GeometryError> {
        let requested =
            mesh.vertices.len() * 8 + (mesh.indices.len() + mesh.face_indices.len()) * 4;
        let peak =
            self.base_resident_bytes() + other_mesh_bytes + Self::mesh_bytes(mesh) + requested;
        if peak > BYTES {
            return Err(failure(read::Error::Limit));
        }
        let work = mesh.vertices.len() + mesh.indices.len() + mesh.face_indices.len();
        let total = proof_work
            .checked_add(work)
            .ok_or_else(|| failure(read::Error::Limit))?;
        if total > 1_600_000_000 {
            return Err(failure(read::Error::Limit));
        }
        Ok((
            Mesh {
                vertices: mesh.vertices.clone(),
                indices: mesh.indices.clone(),
                face_indices: mesh.face_indices.clone(),
            },
            MeshAccounting {
                work: total,
                actual: actual + work,
                peak: proof_peak.max(peak),
            },
        ))
    }
    pub(super) fn faces(&self) -> Result<Vec<Face>, GeometryError> {
        if self.faces.borrow().is_none() {
            let mut faces = Vec::with_capacity(self.source.faces().len());
            let available = (BYTES - OVERHEAD)
                .checked_sub(self.cached_bytes())
                .and_then(|v| {
                    v.checked_sub(
                        self.occurrence
                            .as_ref()
                            .map_or(0, |p| p.retained_budget_bytes()),
                    )
                })
                .ok_or_else(|| failure(read::Error::Limit))?;
            for index in 0..self.source.faces().len() {
                let measure = self
                    .plan
                    .face_measure(
                        index,
                        read::FaceMeasureOptions {
                            centroid_width_mm: self.policy.linear_tolerance_mm,
                            maximum_bytes: available,
                            ..Default::default()
                        },
                        None,
                    )
                    .map_err(failure)?;
                faces.push(OriginalFace {
                    area: measure.area_mm2(),
                    center: measure.centroid_mm(),
                });
            }
            *self.faces.borrow_mut() = Some(faces);
        }
        self.faces
            .borrow()
            .as_ref()
            .expect("qualified original face measures")
            .iter()
            .copied()
            .enumerate()
            .map(|(index, original)| {
                let measure = self.face_properties(original)?;
                let center_mm = measure.center.map(|b| b[0].midpoint(b[1]));
                self.policy.point(center_mm)?;
                Ok(Face {
                    index: index as u32,
                    plane: self.plane(index).map(|p| (p.origin_mm, p.normal())),
                    center_mm,
                    area_mm2: measure.area[0].midpoint(measure.area[1]),
                })
            })
            .collect()
    }
    fn origin(&self, kind: &str, entity: u32) -> Origin {
        Origin {
            source: self.feature,
            role: format!("step.literal-npt/{kind}/#{entity}"),
        }
    }
    pub(crate) fn topology(&self, owner: EntityId, revision: u64) -> TopologySnapshot {
        let mut elements = Vec::new();
        for (index, face) in self.source.faces().iter().enumerate() {
            let origin = self.origin("face", face.entity());
            let plane = self.plane(index);
            let mut boundary = BTreeSet::new();
            for coedge in face.coedges() {
                let edge = &self.source.edges()[coedge.edge_index()];
                boundary.insert(self.origin("edge", edge.entity()));
                for vertex in edge.vertex_indices() {
                    boundary.insert(self.origin("vertex", self.source.vertices()[vertex].entity()));
                }
            }
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([origin.clone()]),
                boundary,
                frames: plane
                    .into_iter()
                    .map(|frame| AnchoredFrame {
                        origin: origin.clone(),
                        frame,
                    })
                    .collect(),
                geometry: plane.map(|frame| AnalyticGeometry::Plane { frame }),
            });
        }
        for (index, edge) in self.source.edges().iter().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([self.origin("edge", edge.entity())]),
                boundary: edge
                    .vertex_indices()
                    .map(|v| self.origin("vertex", self.source.vertices()[v].entity()))
                    .into_iter()
                    .collect(),
                frames: vec![],
                geometry: None,
            });
        }
        for (index, vertex) in self.source.vertices().iter().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Vertex,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([self.origin("vertex", vertex.entity())]),
                boundary: BTreeSet::new(),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Point {
                    position_mm: self.vertices[index],
                }),
            });
        }
        TopologySnapshot {
            owner,
            revision,
            elements,
        }
    }
    pub(crate) fn planar_loops(&self, index: u32) -> Vec<Vec<(u32, bool)>> {
        let Some(face) = self.source.faces().get(index as usize) else {
            return vec![];
        };
        if !matches!(face.support(), read::Support::Plane { .. }) {
            return vec![];
        }
        vec![
            face.coedges()
                .iter()
                .map(|c| (c.edge_index() as u32, !c.forward()))
                .collect(),
        ]
    }
    pub(super) fn sample_edge(
        &self,
        index: u32,
        q: Tessellation,
    ) -> Option<Result<Vec<[f64; 3]>, MeshError>> {
        self.source.edges().get(index as usize)?;
        Some((|| {
            quality(q)?;
            let available = self
                .available_query_bytes()
                .map_err(|_| MeshError::BudgetExceeded)?;
            let limit = q.triangle_limit.min(131072);
            let available = available
                .checked_sub(limit * 24)
                .ok_or(MeshError::BudgetExceeded)?;
            let sampled = self
                .source
                .sample_edge(
                    index as usize,
                    self.occurrence
                        .as_ref()
                        .map_or(q.deflection, |p| q.deflection / p.linear_l1_norm() * 0.5),
                    q.triangle_limit.min(131072),
                    read::Options {
                        maximum_bytes: available,
                        ..Default::default()
                    },
                    None,
                )
                .map_err(|cause| match cause {
                    read::Error::Invalid => MeshError::InvalidOptions,
                    read::Error::Limit => MeshError::BudgetExceeded,
                    _ => MeshError::UnresolvedBoundary,
                })?;
            let per_point = if self.occurrence.is_some() { 259 } else { 3 };
            if sampled.points_mm().len() * per_point > sampled.remaining_work() {
                return Err(MeshError::BudgetExceeded);
            }
            self.place_polyline(&sampled, q)
        })())
    }
    fn record_mesh(&self, accounting: MeshAccounting, fresh: bool) {
        self.last_mesh.set(Some((accounting.work, accounting.peak)));
        self.last_actual.set(Some(accounting.actual));
        if fresh {
            self.mesh_builds.set(self.mesh_builds.get() + 1);
        }
    }
    fn cached_bytes(&self) -> usize {
        self.display
            .borrow()
            .as_ref()
            .map_or(0, |cache| Self::mesh_bytes(&cache.mesh))
    }
    fn proved(&self, q: Tessellation) -> Result<BuiltMesh, GeometryError> {
        quality(q).map_err(|_| failure(read::Error::Invalid))?;
        // A fresh query must coexist with the last valid display. Failure
        // leaves that cache and all successful query evidence untouched.
        let cached = self.cached_bytes();
        let available = (BYTES - OVERHEAD)
            .checked_sub(cached)
            .ok_or_else(|| failure(read::Error::Limit))?;
        let options = read::LiteralMeshOptions {
            tessellation: q,
            maximum_bytes: available,
            ..Default::default()
        };
        let (source, proof, accounting) = if let Some(placement) = &self.occurrence {
            let qualified = self
                .plan
                .placed_mesh(placement.clone(), options, None)
                .map_err(failure)?;
            let proof = qualified.embedding().connectivity();
            let accounting = MeshAccounting {
                work: qualified.work_used(),
                actual: qualified.query_work_used(),
                peak: qualified.accounted_peak_bytes() + OVERHEAD + cached,
            };
            (qualified.into_mesh(), proof, accounting)
        } else {
            let qualified = self
                .plan
                .mesh_candidate(options, None)
                .map_err(failure)?
                .certify_embedding(None)
                .map_err(failure)?;
            let proof = qualified.embedding().connectivity();
            let accounting = MeshAccounting {
                work: qualified.candidate().work_used(),
                actual: qualified.candidate().query_work_used(),
                peak: qualified.candidate().accounted_peak_bytes() + OVERHEAD + cached,
            };
            (qualified.into_mesh(), proof, accounting)
        };
        if accounting.peak > BYTES {
            return Err(failure(read::Error::Limit));
        }
        Ok(BuiltMesh {
            mesh: Mesh {
                vertices: source.vertices,
                indices: source.indices,
                face_indices: source.face_indices,
            },
            proof,
            accounting,
        })
    }
    pub(super) fn mesh(&self, q: Tessellation) -> Result<Mesh, GeometryError> {
        let key = [
            q.deflection.to_bits(),
            q.angle.to_bits(),
            q.triangle_limit as u64,
        ];
        if let Some(cached) = self.display.borrow().as_ref().filter(|c| c.key == key) {
            let (copy, accounting) =
                self.copy_display(&cached.mesh, cached.work, cached.peak, 0, 0)?;
            self.record_mesh(accounting, false);
            return Ok(copy);
        }
        let built = self.proved(q)?;
        let (copy, accounting) = self.copy_display(
            &built.mesh,
            built.accounting.work,
            built.accounting.peak,
            built.accounting.actual,
            self.cached_bytes(),
        )?;
        *self.display.borrow_mut() = Some(CachedDisplay {
            key,
            mesh: built.mesh,
            work: built.accounting.work,
            peak: built.accounting.peak,
        });
        self.record_mesh(accounting, true);
        Ok(copy)
    }
    pub(super) fn manufacturing(
        &self,
        q: Tessellation,
    ) -> Result<crate::ManufacturingMesh, GeometryError> {
        let built = self.proved(q)?;
        self.record_mesh(built.accounting, true);
        Ok(crate::ManufacturingMesh {
            mesh: built.mesh,
            maximum_displacement_mm: 0.,
            edge_pairs: built.proof.edge_pairs,
            vertex_links: built.proof.vertex_links,
            removed_native_pole_triangles: 0,
            removed_exact_zero_triangles: 0,
            exact_vertex_pairs: 0,
            internal_constraint_edges: 0,
            removed_meshing_internal_edges: 0,
            removed_meshing_constraint_wire_containers: 0,
            removed_meshing_internal_vertices: 0,
            remeshed_without_internal_constraints: false,
        })
    }
}
#[cfg(test)]
mod tests;
