//! Cap coedges of the SAME source-owned, certified literal mesh. No point
//! matching, inverse fitting, or replacement source law is involved.
use super::*;
use crate::display::Tessellation;

/// Directed original cap boundary and its two literal mesh vertex indices.
/// Edges are stored in source-cell traversal order, not cyclic ring order.
/// `vertices` follows the outward cap coedge (start clockwise, end CCW).
#[derive(Clone, Copy, Debug)]
pub struct StockClipCapTraceEdge {
    pub vertices: [u32; 2],
    pub owners: [StockClipPointOwner; 2],
    pub station: StockClipStation,
    pub patch: StockClipPatch,
    /// Original unwrapped (s,z) endpoint enclosures, not rounded UV proxies.
    pub parameters: [[[f64; 2]; 2]; 2],
    pub source_positions_mm: [[[f64; 2]; 3]; 2],
    /// Whole original finite curve to its actual outgoing literal chord.
    /// Also proves the strict common halfplane needed for the paired-loop
    /// homotopy; this is not a whole composite embedding receipt.
    pub chord_error_mm: f64,
}

#[derive(Debug)]
pub struct StockClipTracedMesh<'a> {
    source: &'a StockClipCellPlan,
    candidate: StockClipMeshCandidate,
    edges: Vec<StockClipCapTraceEdge>,
    work: usize,
    peak: usize,
}
impl StockClipTracedMesh<'_> {
    pub fn source(&self) -> &StockClipCellPlan {
        self.source
    }
    pub fn candidate(&self) -> &StockClipMeshCandidate {
        &self.candidate
    }
    pub fn cap_edges(&self) -> &[StockClipCapTraceEdge] {
        &self.edges
    }
    /// Includes the unchanged candidate's inherited preparation work once.
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Source plan, candidate construction peak, and live trace allocation.
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
impl StockClipCellPlan {
    pub fn mesh_candidate_with_cap_trace(
        &self,
        quality: Tessellation,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipTracedMesh<'_>, Error> {
        quality.check().map_err(|_| Error::Invalid)?;
        if maximum_work == 0
            || maximum_work > 512_000_000
            || maximum_bytes == 0
            || maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let metadata = size_of::<StockClipTracedMesh<'_>>()
            + self.caps.len() * size_of::<StockClipCapTraceEdge>();
        let room = maximum_bytes.checked_sub(metadata).ok_or(Error::Limit)?;
        if room == 0 {
            return Err(Error::Limit);
        }
        let candidate = self.mesh_candidate(quality, maximum_work, room, interrupt)?;
        let mut work = Work {
            used: candidate.work_used(),
            limit: maximum_work,
            bytes: maximum_bytes,
            peak: candidate.retained_budget_bytes() + metadata,
            interrupt,
        };
        work.charge(self.caps.len())?;
        let mut edges = Vec::new();
        edges
            .try_reserve_exact(self.caps.len())
            .map_err(|_| Error::Limit)?;
        work.memory(
            candidate.retained_budget_bytes()
                + size_of::<StockClipTracedMesh<'_>>()
                + edges.capacity() * size_of::<StockClipCapTraceEdge>(),
        )?;
        for edge in &self.caps {
            work.charge(512)?;
            let chord_error_mm = mesh::cap_edge(self, edge, candidate.mesh(), quality)?;
            edges.push(StockClipCapTraceEdge {
                vertices: edge.ids,
                owners: edge.ids.map(|id| self.vertices[id as usize].owner),
                station: edge.station,
                patch: edge.patch,
                parameters: edge.parameters.map(|p| p.map(bounds)),
                source_positions_mm: edge
                    .ids
                    .map(|id| self.vertices[id as usize].position.map(bounds)),
                chord_error_mm,
            });
        }
        Ok(StockClipTracedMesh {
            source: self,
            candidate,
            edges,
            work: work.used,
            peak: work.peak,
        })
    }
}

#[cfg(test)]
mod tests;
