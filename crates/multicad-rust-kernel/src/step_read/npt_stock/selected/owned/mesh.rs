//! Literal C11 facets. Source correspondence and indexed connectivity are
//! earned here; global literal embedding remains a distinct consumer gate.
use super::*;
use crate::display::{DisplayMesh, Tessellation};
mod angular;
mod cells;
mod embedding;
pub use embedding::LiteralEmbeddingProof;
#[cfg(test)]
mod endpoint_tests;
mod hessian;
mod normal;
mod occurrence;
pub use occurrence::PlacedLiteralMesh;
mod plan;
pub use plan::{LiteralGeometryOptions, LiteralGeometryPlan};
mod sheared;
// New literal-mesh aggregate ceiling: measured both-size candidates use
// about1.1B charged operations, including the immutable preparation once.
// Existing reader/material/query ceilings are unchanged.
const WORK: usize = 1_600_000_000;
const BYTES: usize = 64 * 1024 * 1024;
const SCRATCH: usize = 1024 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct LiteralMeshOptions {
    pub tessellation: Tessellation,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for LiteralMeshOptions {
    fn default() -> Self {
        Self {
            tessellation: Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
impl LiteralMeshOptions {
    fn check(self) -> Result<()> {
        self.tessellation.check().map_err(|_| Error::Invalid)?;
        if self.maximum_work == 0
            || self.maximum_work > WORK
            || self.maximum_bytes == 0
            || self.maximum_bytes > BYTES
        {
            return Err(Error::Invalid);
        }
        Ok(())
    }
}
/// Sealed original edge/span parameters of one emitted cap chord. This is
/// query metadata, not an occupied/embedding certificate or a new constructor.
#[derive(Clone, Copy)]
pub(crate) struct LiteralCapEdge {
    pub(crate) edge: usize,
    pub(crate) span: usize,
    pub(crate) parameters: [f64; 2],
    pub(crate) vertices: [u32; 2],
}
pub(crate) struct LiteralCapTrace<'a> {
    pub(crate) face: usize,
    pub(crate) apex: u32,
    pub(crate) rim: &'a [u32],
    pub(crate) edges: &'a [LiteralCapEdge],
}
pub struct LiteralMeshCandidate {
    source: OwnedMaterial,
    mesh: DisplayMesh,
    connectivity: crate::manufacturing::MeshProof,
    position: f64,
    angle: f64,
    work: usize,
    preparation: usize,
    remaining: usize,
    peak: usize,
    root_triangles: Vec<Vec<usize>>,
    caps: [embedding::Cap; 2],
    retained_phase_bytes: usize,
    byte_limit: usize,
}
pub struct LiteralMesh {
    candidate: LiteralMeshCandidate,
    embedding: LiteralEmbeddingProof,
}
impl LiteralMesh {
    pub fn into_mesh(self) -> DisplayMesh {
        self.candidate.mesh
    }
    pub fn candidate(&self) -> &LiteralMeshCandidate {
        &self.candidate
    }
    pub fn mesh(&self) -> &DisplayMesh {
        self.candidate.mesh()
    }
    pub fn source(&self) -> &OwnedMaterial {
        self.candidate.source()
    }
    pub fn embedding(&self) -> LiteralEmbeddingProof {
        self.embedding
    }
}
impl LiteralMeshCandidate {
    pub(crate) fn cap_traces(&self) -> [LiteralCapTrace<'_>; 2] {
        self.caps.each_ref().map(|c| LiteralCapTrace {
            face: c.face,
            apex: c.apex,
            rim: &c.rim,
            edges: &c.edges,
        })
    }
    /// Actual buffers still retained after the construction plan/registry
    /// drop. SAME source Data allocation is shared with the caller's owner.
    pub(crate) fn retained_facet_bytes(&self) -> usize {
        size_of::<Self>()
            + self.mesh.vertices.capacity() * 8
            + (self.mesh.indices.capacity() + self.mesh.face_indices.capacity()) * 4
            + self.root_triangles.capacity() * size_of::<Vec<usize>>()
            + self
                .root_triangles
                .iter()
                .map(|r| r.capacity() * size_of::<usize>())
                .sum::<usize>()
            + self
                .caps
                .iter()
                .map(|c| c.rim.capacity() * 4 + c.edges.capacity() * size_of::<LiteralCapEdge>())
                .sum::<usize>()
    }
    /// Independently earns embedding of these exact rounded vertices. It
    /// continues the original aggregate work/storage ledger, without repair.
    pub fn certify_embedding(
        mut self,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LiteralMesh> {
        let mut budget = Budget {
            used: self.work,
            remaining: self.remaining,
            bytes: self.retained_phase_bytes,
            peak: self.peak,
            limit: self.byte_limit,
            interrupt,
        };
        budget.step(1)?;
        let embedding = embedding::certify(&self.mesh, &self.caps, &mut budget)?;
        self.work = budget.used;
        self.remaining = budget.remaining;
        self.peak = budget.peak;
        Ok(LiteralMesh {
            candidate: self,
            embedding,
        })
    }
    pub fn source(&self) -> &OwnedMaterial {
        &self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn connectivity(&self) -> crate::manufacturing::MeshProof {
        self.connectivity
    }
    pub fn position_error_mm(&self) -> f64 {
        self.position
    }
    pub fn normal_error_radians(&self) -> f64 {
        self.angle
    }
    /// Aggregate receipt: inherited immutable preparation once plus this query.
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn preparation_work_used(&self) -> usize {
        self.preparation
    }
    pub fn query_work_used(&self) -> usize {
        self.work - self.preparation
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
    pub fn original_root_triangle_counts(&self) -> &[Vec<usize>] {
        &self.root_triangles
    }
}
fn interval_point(source: &OwnedMaterial, index: usize) -> [I; 3] {
    source.0.points[index]
        .coordinates
        .each_ref()
        .map(|s| interval(s.bounds))
}
fn basis<const N: usize>(t: f64) -> [I; N] {
    let t = I::point(t);
    let s = I::point(1.) - t;
    let mut basis = [I::point(0.); N];
    basis[0] = I::point(1.);
    for degree in 1..N {
        for i in (0..=degree).rev() {
            basis[i] = if i == 0 {
                s * basis[i]
            } else if i == degree {
                t * basis[i - 1]
            } else {
                s * basis[i] + t * basis[i - 1]
            };
        }
    }
    basis
}
fn side_edges(face: &Face) -> [usize; 4] {
    let sides = face
        .tensor_sides()
        .expect("tensor support has earned pcurve side bindings");
    let mut edges = [usize::MAX; 4];
    for (i, side) in sides.into_iter().enumerate() {
        edges[side] = face.coedges[i].edge;
    }
    edges
}
impl OwnedMaterial {
    /// Facets of the retained ORIGINAL decimal surfaces, in the actual STEP
    /// world frame. This earns source error and indexed connectivity only;
    /// callers must earn global literal embedding before manufacturing use.
    pub fn mesh_candidate(
        &self,
        options: LiteralMeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LiteralMeshCandidate> {
        options.check()?;
        let plan = self.geometry_plan(
            LiteralGeometryOptions {
                maximum_work: options.maximum_work.min(512_000_000),
                maximum_bytes: options.maximum_bytes,
            },
            interrupt,
        )?;
        plan.mesh_candidate(options, interrupt)
    }
}
impl LiteralGeometryPlan {
    pub fn mesh_candidate(
        &self,
        options: LiteralMeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LiteralMeshCandidate> {
        options.check()?;
        let bytes = self.retained_budget_bytes() + SCRATCH;
        if self.preparation_work_used() > options.maximum_work {
            return Err(Error::Limit);
        }
        if bytes > options.maximum_bytes || self.construction_peak_bytes() > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut budget = Budget {
            used: self.preparation_work_used(),
            remaining: options.maximum_work - self.preparation_work_used(),
            bytes,
            peak: bytes.max(self.construction_peak_bytes()),
            limit: options.maximum_bytes,
            interrupt,
        };
        budget.step(1)?;
        let result = angular::build(self, options, &mut budget)?;
        Ok(LiteralMeshCandidate {
            source: self.source().clone(),
            mesh: result.mesh,
            connectivity: result.connectivity,
            position: result.position,
            angle: result.angle,
            work: budget.used,
            preparation: self.preparation_work_used(),
            remaining: budget.remaining,
            peak: budget.peak,
            root_triangles: result.root_triangles,
            caps: result.caps,
            retained_phase_bytes: budget.bytes,
            byte_limit: options.maximum_bytes,
        })
    }
}
#[cfg(test)]
mod tests;
