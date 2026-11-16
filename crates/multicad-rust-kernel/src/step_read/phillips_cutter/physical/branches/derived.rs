//! A separately named data-derived cavity source and measured raw correspondence.
//! No exact raw closure, finite-face filling, embedding, material or mass receipt.
use super::*;
use crate::phillips_drive::{CurveKey, Cutter, FaceKey, VertexKey};
mod closed;
mod fit;
pub use closed::ClosedDerivedCavityCertificate;
#[cfg(test)]
mod tests;
// Only the 49V/76E/29F rounded-inner family is proposed. Construction/topology
// have four fixed corner orbits, bounded section polygons and no adaptive loops.
// Reserve a million source-stage work units before invoking their fixed algebra;
// the separate variable imported spline proof charges 256 units/span as before.
const SOURCE_WORK: usize = 1_000_000;
// geometry::spline schema: <=4097 knots, <=12289 POINT refs/positions,
// <=4096 CubicTrimSpans. Fresh Vec collections under pinned Rust1.99.1
// grow geometrically (never beyond 2*maximum logical length, including the
// minimum4-element allocation). Charge FOUR times each maximum payload to
// include old+new storage during reallocation, plus64KiB headers/small buffers.
// Only one edge is decoded at a time; the next replaces/drops these buffers.
// Shell::point(&self) does not cache poles (unlike vertex(&mut self), unused).
const RAW_SPLINE_SCRATCH: usize = 4
    * (4097 * std::mem::size_of::<f64>()
        + 12289 * (std::mem::size_of::<u32>() + std::mem::size_of::<Point>())
        + 4096 * std::mem::size_of::<CubicTrimSpan>())
    + 64 * 1024;
// Separate fixed orbit topology, views/graph/maps/links and stack scratch.
// <=128 source vertices/edges/faces; imported graph49V/76E/29F and152 uses.
// No raw AST clone or accumulating pole cache exists in this path.
const SCRATCH_BYTES: usize = RAW_SPLINE_SCRATCH + 4 * 1024 * 1024;
const DEFAULT_TRANSFER: f64 = 1e-7 / 16.;
#[derive(Clone, Copy)]
pub struct DerivedCavityOptions {
    /// Measured correspondence limit, independent of the proposed dimensions.
    /// The default preserves the existing source-comparison reserve (1e-7/16).
    pub maximum_transfer_mm: f64,
    /// Includes prerequisite work plus the conservative fixed source-stage charge.
    pub maximum_work: usize,
    /// Includes the original Body once and the live prerequisite/workspace peak.
    pub maximum_bytes: usize,
}
impl Default for DerivedCavityOptions {
    fn default() -> Self {
        Self {
            maximum_transfer_mm: DEFAULT_TRANSFER,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct DerivedVertexCorrespondence {
    pub source_role: VertexKey,
    pub raw_vertex_entity: u32,
    pub raw_point_entity: u32,
}
#[derive(Clone, Copy, Debug)]
pub struct DerivedEdgeCorrespondence {
    pub source_role: CurveKey,
    pub raw_edge_entity: u32,
    pub raw_curve_entity: u32,
    pub raw_edge_reversed: bool,
    pub raw_curve_forward: bool,
}
#[derive(Clone, Copy, Debug)]
pub struct DerivedFaceCorrespondence {
    pub source_role: FaceKey,
    pub raw_material_entity: Option<u32>,
    pub raw_surface_entity: u32,
    pub raw_material_forward: bool,
    pub raw_cavity_forward: bool,
}
/// Immutable NEW coherent cavity, fitted from actual imported supports and
/// checked against every original vertex, edge, oriented loop and support.
/// Its points are not original POINT leaves; the original positive raw gaps
/// remain facts about Body. The measured transfer does not clear raw uncertainty
/// or prove filling, embedding, outer containment, manufacturing or file mass.
/// All accessors are constant-time borrows/copies with no hidden geometry work.
pub struct DerivedCavityRepresentative<'a, 'input> {
    body: &'a Body<'input>,
    source: Cutter,
    frame: [Point; 4],
    vertices: [DerivedVertexCorrespondence; 49],
    edges: [DerivedEdgeCorrespondence; 76],
    faces: [DerivedFaceCorrespondence; 29],
    transfer: PhillipsSourceKnownProof,
    sewing_collar_mm: f64,
    work_used: usize,
    retained_bytes: usize,
}
impl<'a, 'input> DerivedCavityRepresentative<'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.body
    }
    /// Distinct derived source, never a recovered authoring recipe or raw BRep.
    pub fn derived_cavity_source(&self) -> &Cutter {
        &self.source
    }
    /// Raw affine placement of the derived source: origin, U and V. The normal
    /// is their actual cross product; no rigid mass/area invariance is asserted.
    pub fn derived_frame(&self) -> [Point; 3] {
        [self.frame[0], self.frame[1], self.frame[2]]
    }
    /// Actual stored origin,U,V,N columns used in source comparison.
    pub fn derived_affine_columns(&self) -> [Point; 4] {
        self.frame
    }
    pub fn vertices(&self) -> &[DerivedVertexCorrespondence; 49] {
        &self.vertices
    }
    pub fn edges(&self) -> &[DerivedEdgeCorrespondence; 76] {
        &self.edges
    }
    pub fn faces(&self) -> &[DerivedFaceCorrespondence; 29] {
        &self.faces
    }
    /// Source-known comparison to this independently proposed derived source;
    /// explicitly not exact raw incidence or independent decoded mass.
    pub fn transfer_bounds(&self) -> PhillipsSourceKnownProof {
        self.transfer
    }
    /// Actual earned endpoint-correction plus raw-owner enclosure on each
    /// prerequisite sewn edge. This is not the caller's requested tolerance.
    pub fn inherited_sewing_collar_mm(&self) -> f64 {
        self.sewing_collar_mm
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained_bytes
    }
}
impl<'input> Body<'input> {
    /// Consumes no source and publishes atomically. The earned prerequisite
    /// supplies exact schema/POINT authority and cannot refer to a foreign arena.
    /// Its owned payload is live during construction but need not be retained
    /// afterward: the original immutable Body and all new proof data remain.
    pub fn derived_cavity_representative<'a>(
        &'a self,
        prerequisite: &SewnEdgeSources<'_, 'input>,
        options: DerivedCavityOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> BranchResult<DerivedCavityRepresentative<'a, 'input>> {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: std::cell::Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !options.maximum_transfer_mm.is_finite()
                || options.maximum_transfer_mm <= 0.
                || options.maximum_transfer_mm > DEFAULT_TRANSFER
                || options.maximum_work > WORK
                || options.maximum_bytes == 0
                || options.maximum_bytes > BYTES
                || !std::ptr::eq(self, prerequisite.body())
            {
                return Err(ReadError::Invalid);
            }
            work.step(prerequisite.work_used())?;
            let own = std::mem::size_of::<DerivedCavityRepresentative<'_, '_>>();
            // Saturation is an honest refusal. This check precedes even the
            // proposal's small cones Vec and all source/topology allocations.
            let peak = prerequisite
                .retained_capacity_bytes()
                // Prerequisite getters exclude the borrowed original arena.
                // Body includes the shared Boundary/keys/raw input exactly once.
                .saturating_add(self.retained_budget_bytes())
                .saturating_add(own)
                .saturating_add(SCRATCH_BYTES);
            if peak > options.maximum_bytes {
                return Err(ReadError::Limit);
            }
            work.step(SOURCE_WORK)?;
            work.step(prerequisite.edges().len() * 16)?;
            let mut sewing_collar_mm = 0f64;
            for edge in prerequisite.edges() {
                let correction = edge.position_transfer_mm();
                let owners = edge.raw_owner_distances_mm();
                let collar = sewn_edge_collar(correction, owners)?;
                sewing_collar_mm = sewing_collar_mm.max(collar);
            }
            let candidate = fit::candidate(self.cavity(), &mut work)?;
            let topology = candidate
                .source
                .boundary_topology_with_cancel(&|| work.interrupted())
                .map_err(cutter_error)?;
            let boundary = self.cavity();
            let views = boundary
                .faces()
                .iter()
                .map(|f| SourceFaceView {
                    surface: f.surface_entity,
                    forward: f.cavity_forward,
                    uses: f
                        .uses
                        .iter()
                        .map(|&(e, s)| (boundary.edges()[e].entity, s))
                        .collect(),
                })
                .collect::<Vec<_>>();
            let mut shell = Shell::new(
                self.source_entities(),
                options.maximum_transfer_mm * 16.,
                true,
            );
            let graph = graph::read_views(
                &mut shell,
                &views,
                &topology,
                &candidate.pose,
                options.maximum_transfer_mm * 16.,
                &mut work,
            )?;
            let transfer = geometry::verify(
                &shell,
                &graph,
                &topology,
                &candidate.pose,
                options.maximum_transfer_mm * 16.,
                &mut work,
            )?;
            // Graph::read_views established equal lengths and an actual unique
            // bijection, including endpoint reversal and cyclic oriented loops.
            let vertices = std::array::from_fn(|i| DerivedVertexCorrespondence {
                source_role: topology.vertices()[i].source_identity(),
                raw_vertex_entity: graph.vertices[i].0,
                raw_point_entity: graph.vertices[i].1,
            });
            let edges = std::array::from_fn(|i| DerivedEdgeCorrespondence {
                source_role: topology.edges()[i].source_identity(),
                raw_edge_entity: graph.edges[i].entity,
                raw_curve_entity: graph.edges[i].geometry,
                raw_edge_reversed: graph.edges[i].reversed,
                raw_curve_forward: graph.edges[i].forward,
            });
            let faces = std::array::from_fn(|i| {
                let f = &boundary.faces()[graph.faces[i].input_view];
                DerivedFaceCorrespondence {
                    source_role: topology.faces()[i].source_identity(),
                    raw_material_entity: f.material_entity,
                    raw_surface_entity: f.surface_entity,
                    raw_material_forward: f.material_forward,
                    raw_cavity_forward: f.cavity_forward,
                }
            });
            drop(topology);
            // Source and mappings have no heap payload. Borrowed Body budget is
            // cached under its original reader work and includes its raw input.
            let retained_bytes = self.retained_budget_bytes().saturating_add(own);
            work.step(1)?;
            Ok(DerivedCavityRepresentative {
                body: self,
                source: candidate.source,
                frame: [
                    candidate.pose.origin,
                    candidate.pose.axes[0],
                    candidate.pose.axes[1],
                    candidate.pose.axes[2],
                ],
                vertices,
                edges,
                faces,
                transfer,
                sewing_collar_mm,
                work_used: options.maximum_work - work.remaining,
                retained_bytes,
            })
        })();
        result.map_err(|e| {
            if work.cancelled.get() {
                PhillipsReadError::Cancelled
            } else {
                e.into()
            }
        })
    }
}

#[cfg(test)]
mod outer_probe;

// Compose the actual earned position corrections, rejecting nonfinite input or
// arithmetic overflow before publishing a world-space enclosure.
fn sewn_edge_collar(correction: f64, owners: [f64; 2]) -> Result<f64> {
    require(
        std::iter::once(correction)
            .chain(owners)
            .all(|v| v.is_finite() && v >= 0.),
    )?;
    let collar = (crate::tapered_thread::interval::I::point(correction)
        + crate::tapered_thread::interval::I::point(owners[0].max(owners[1])))
    .hi;
    require(collar.is_finite() && collar >= 0.)?;
    Ok(collar)
}
