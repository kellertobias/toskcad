//! Explicit Native provider of independently occupied actual mapped STEP head
//! geometry. No authored screw recipe or OpenCascade geometry is reconstructed.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
use multicad_rust_kernel::{
    display::{MeshError, Tessellation},
    step_read::npt_stock as read,
};
use read::full_head::{
    measure::{HeadGeometryOptions, HeadGeometryPlan, HeadMeasureOptions},
    transport::MappedActualHead,
};
use std::{
    cell::{Cell, RefCell},
    sync::Arc,
};
mod measure;
mod mesh;
mod planes;
mod topology;
const BYTES: usize = 64 * 1024 * 1024;
// Fixed24vertices,17optional scalar descriptors, owner/cache headers, source
// transforms and all fixed Native arithmetic temporaries are bounded here.
// Escaping topology snapshots are a separate query output, not a body cache.
const BODY_SCOPE: usize = 64 * 1024;
#[derive(Clone, Copy)]
struct LocalFace {
    area: [f64; 2],
    center: [[f64; 2]; 3],
}
struct CachedDisplay {
    key: [u64; 3],
    mesh: Mesh,
    record: MeshRecord,
}
#[derive(Clone, Copy)]
struct MeshRecord {
    work: usize,
    fresh: usize,
    peak: usize,
    source_error: f64,
    normal: f64,
    pairs: usize,
}
pub struct ActualClippedHeadRuntimeEvidence {
    pub native_preparation_work: usize,
    pub resident_budget_bytes: usize,
    pub measured_faces: usize,
    pub face_measure_builds: usize,
    pub mesh_builds: usize,
    pub last_mesh_work: Option<usize>,
    pub last_mesh_fresh_work: Option<usize>,
    pub last_mesh_peak_bytes: Option<usize>,
    pub last_mesh_source_error_mm: Option<f64>,
    pub last_mesh_normal_tangent_cone: Option<f64>,
    pub last_embedding_pairs: Option<usize>,
}
type FaceCacheState = (Rc<RefCell<[Option<LocalFace>; 17]>>, Rc<Cell<usize>>);
pub(crate) struct ActualClippedHeadBody {
    pub(super) feature: EntityId,
    source: Arc<MappedActualHead>,
    occurrence: Option<read::LiteralOccurrencePlacement>,
    policy: GeometryPolicy,
    vertices: Vec<[f64; 3]>,
    local_faces: Rc<RefCell<[Option<LocalFace>; 17]>>,
    face_builds: Rc<Cell<usize>>,
    display: RefCell<Option<CachedDisplay>>,
    last_mesh: Cell<Option<MeshRecord>>,
    mesh_builds: Cell<usize>,
    preparation: usize,
}
fn failure(cause: read::Error) -> GeometryError {
    error(
        "actual clipped head",
        match cause {
            read::Error::Invalid => GeometryErrorKind::Input,
            read::Error::Unsupported => GeometryErrorKind::Unsupported,
            read::Error::Limit => GeometryErrorKind::Limit,
            read::Error::Cancelled => GeometryErrorKind::Cancelled,
        },
        format!("Actual retained geometry operation refused: {cause:?}"),
    )
}
fn tick(
    used: &mut usize,
    amount: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        return Err(failure(read::Error::Cancelled));
    }
    *used = used
        .checked_add(amount)
        .filter(|n| *n <= 16_000_000)
        .ok_or_else(|| failure(read::Error::Limit))?;
    Ok(())
}
fn ranges(p: [[f64; 2]; 3]) -> [I; 3] {
    p.map(|v| I::new(v[0], v[1]))
}
fn output(p: [I; 3]) -> [[f64; 2]; 3] {
    p.map(|v| [v.lo, v.hi])
}
fn mesh_bytes(m: &Mesh) -> usize {
    m.vertices.capacity() * 8 + (m.indices.capacity() + m.face_indices.capacity()) * 4
}
impl Solid {
    pub fn from_actual_clipped_head_source(
        source: Arc<MappedActualHead>,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        if feature.is_nil() {
            return Err(input("Actual head import requires a source identity"));
        }
        let body = ActualClippedHeadBody::build(source, feature, policy, None, 0, None, interrupt)?;
        Ok(Self {
            body: SolidBody::ActualClippedHead(Rc::new(body)),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn actual_clipped_head_body(&self) -> Option<&ActualClippedHeadBody> {
        match &self.body {
            SolidBody::ActualClippedHead(v) => Some(v),
            _ => None,
        }
    }
    pub fn rust_actual_clipped_head_source(&self) -> Option<&MappedActualHead> {
        self.actual_clipped_head_body().map(|v| v.source.as_ref())
    }
    pub fn rust_actual_clipped_head_occurrence(&self) -> Option<&read::LiteralOccurrencePlacement> {
        self.actual_clipped_head_body()
            .and_then(|v| v.occurrence.as_ref())
    }
    pub fn rust_actual_clipped_head_runtime_evidence(
        &self,
    ) -> Option<ActualClippedHeadRuntimeEvidence> {
        self.actual_clipped_head_body().map(|v| {
            let r = v.last_mesh.get();
            ActualClippedHeadRuntimeEvidence {
                native_preparation_work: v.preparation,
                resident_budget_bytes: v.resident_bytes(),
                measured_faces: v
                    .local_faces
                    .borrow()
                    .iter()
                    .filter(|v| v.is_some())
                    .count(),
                face_measure_builds: v.face_builds.get(),
                mesh_builds: v.mesh_builds.get(),
                last_mesh_work: r.map(|v| v.work),
                last_mesh_fresh_work: r.map(|v| v.fresh),
                last_mesh_peak_bytes: r.map(|v| v.peak),
                last_mesh_source_error_mm: r.map(|v| v.source_error),
                last_mesh_normal_tangent_cone: r.map(|v| v.normal),
                last_embedding_pairs: r.map(|v| v.pairs),
            }
        })
    }
    pub fn rust_actual_clipped_head_contains(
        &self,
        point: [f64; 3],
    ) -> Result<bool, GeometryError> {
        self.policy.point(point)?;
        let b = self
            .actual_clipped_head_body()
            .ok_or_else(|| failure(read::Error::Unsupported))?;
        let options = read::full_head::query::HeadQueryOptions {
            additional_resident_bytes: b.external_bytes(),
            ..Default::default()
        };
        let input = point.map(|v| [v, v]);
        let q = if let Some(a) = &b.occurrence {
            b.source.classify_placed_box(input, a, options, None)
        } else {
            b.source.classify_box(input, options, None)
        }
        .map_err(failure)?;
        match q.location() {
            read::Location::Inside => Ok(true),
            read::Location::Outside => Ok(false),
            read::Location::Unresolved => Err(failure(read::Error::Unsupported)),
        }
    }
}
impl ActualClippedHeadBody {
    fn build(
        source: Arc<MappedActualHead>,
        feature: EntityId,
        policy: GeometryPolicy,
        occurrence: Option<read::LiteralOccurrencePlacement>,
        inherited: usize,
        cache: Option<FaceCacheState>,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        let mut work = inherited;
        tick(&mut work, 256, interrupt)?;
        let mut result = Self {
            feature,
            source,
            occurrence,
            policy,
            vertices: Vec::new(),
            local_faces: cache
                .as_ref()
                .map_or_else(|| Rc::new(RefCell::new([None; 17])), |v| Rc::clone(&v.0)),
            face_builds: cache
                .as_ref()
                .map_or_else(|| Rc::new(Cell::new(0)), |v| Rc::clone(&v.1)),
            display: RefCell::new(None),
            last_mesh: Cell::new(None),
            mesh_builds: Cell::new(0),
            preparation: 0,
        };
        if result.resident_bytes() > BYTES {
            return Err(failure(read::Error::Limit));
        }
        let b = result.world_box()?;
        policy.point(b.map(|v| v[0]))?;
        policy.point(b.map(|v| v[1]))?;
        result.centroid_shift_guard()?;
        let local = result.source.local_geometry();
        let aux = local.auxiliary_thread();
        result
            .vertices
            .try_reserve_exact(aux.vertices().len() + local.head().vertices.len())
            .map_err(|_| failure(read::Error::Limit))?;
        for vertex in aux.vertices() {
            tick(&mut work, 256, interrupt)?;
            let p = aux
                .point(vertex.point_index())
                .ok_or_else(|| failure(read::Error::Invalid))?
                .coordinates()
                .map(|v| v.bounds());
            let p = local
                .axial_sewing_transfer()
                .map_raw_box_mm(p)
                .map_err(failure)?;
            result
                .vertices
                .push(result.literal_point(result.map_local_box(p)?)?);
        }
        for &p in &local.head().vertices {
            tick(&mut work, 256, interrupt)?;
            result
                .vertices
                .push(result.literal_point(result.map_local_box(p.map(|v| [v, v]))?)?);
        }
        tick(&mut work, 0, interrupt)?;
        result.preparation = work;
        if result.resident_bytes() > BYTES {
            return Err(failure(read::Error::Limit));
        }
        Ok(result)
    }
    fn map_local_box(&self, p: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3], GeometryError> {
        let p = self.source.actual_frame().map_box_mm(p).map_err(failure)?;
        self.occurrence
            .as_ref()
            .map_or(Ok(p), |a| a.map_box_mm(p).map_err(failure))
    }
    fn literal_point(&self, p: [[f64; 2]; 3]) -> Result<[f64; 3], GeometryError> {
        let point = p.map(|v| v[0].midpoint(v[1]));
        let error = ranges(p).into_iter().zip(point).fold(I::ZERO, |s, (v, x)| {
            s + I::point((v - I::point(x)).magnitude())
        });
        if !error.valid() || error.hi > self.policy.linear_tolerance_mm / 8. {
            return Err(failure(read::Error::Unsupported));
        }
        self.policy.point(point)?;
        Ok(point)
    }
    fn external_bytes(&self) -> usize {
        BODY_SCOPE + self.cached_bytes()
    }
    fn resident_bytes(&self) -> usize {
        self.source.receipt().retained_bytes
            + self.external_bytes()
            + self
                .occurrence
                .as_ref()
                .map_or(0, |a| a.retained_budget_bytes())
    }
    fn cached_bytes(&self) -> usize {
        self.display
            .borrow()
            .as_ref()
            .map_or(0, |v| mesh_bytes(&v.mesh))
    }
    fn world_box(&self) -> Result<[[f64; 2]; 3], GeometryError> {
        let p = self.source.world_sewn_enclosure_mm();
        self.occurrence
            .as_ref()
            .map_or(Ok(p), |a| a.map_box_mm(p).map_err(failure))
    }
    pub(super) fn bounds(&self) -> Result<Bounds, GeometryError> {
        let p = self.world_box()?;
        Ok(Bounds {
            min: p.map(|v| v[0]),
            max: p.map(|v| v[1]),
        })
    }
    pub(super) fn volume_mm3(&self) -> f64 {
        let v = self.source.world_sewn_volume_mm3();
        let v = self.occurrence.as_ref().map_or(I::new(v[0], v[1]), |a| {
            let d = a.determinant_bounds();
            I::new(v[0], v[1]) * I::new(d[0], d[1])
        });
        v.lo.midpoint(v.hi)
    }
    fn centroid_shift_guard(&self) -> Result<(), GeometryError> {
        let Some(a) = &self.occurrence else {
            return Ok(());
        };
        let b = self.source.world_sewn_enclosure_mm();
        let span = b
            .into_iter()
            .fold(I::ZERO, |s, v| s + I::point(v[1]) - I::point(v[0]));
        let delta = I::point(a.metric_error());
        let shift = I::point(a.linear_l1_norm()) * span * delta / (I::ONE - delta);
        if !shift.valid() || shift.hi > self.policy.linear_tolerance_mm / 8. {
            return Err(failure(read::Error::Unsupported));
        }
        Ok(())
    }
    pub(in crate::geometry) fn transformed(
        &self,
        t: RigidTransform,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        t.validate(policy)?;
        if self
            .resident_bytes()
            .checked_add(1024 * 1024)
            .is_none_or(|n| n > BYTES)
        {
            return Err(failure(read::Error::Limit));
        }
        let options = read::Options {
            maximum_bytes: 1024 * 1024,
            ..Default::default()
        };
        let occurrence = match &self.occurrence {
            Some(a) => a.transformed(t.rotation, t.translation_mm, options, None),
            None => {
                read::LiteralOccurrencePlacement::new(t.rotation, t.translation_mm, options, None)
            }
        }
        .map_err(failure)?;
        Self::build(
            Arc::clone(&self.source),
            self.feature,
            policy,
            Some(occurrence),
            self.preparation,
            Some((Rc::clone(&self.local_faces), Rc::clone(&self.face_builds))),
            None,
        )
    }
}
