//! ORIGINAL emitted ring ownership, only after all existing mesh certificates.
use super::*;
const BYTES: usize = 32 * 1024 * 1024;
#[derive(Clone, Copy)]
pub struct CapTraceOptions {
    /// Additional retained trace metadata only. Source and mesh storage retain
    /// their existing contracts; this is not a whole-request/RSS byte bound.
    pub maximum_bytes: usize,
}
impl Default for CapTraceOptions {
    fn default() -> Self {
        Self {
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CapRole {
    AxialStart,
    AxialEnd,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapRimFootpoint {
    pub vertex: u32,
    pub side: usize,
    pub side_face: usize,
    pub side_surface: usize,
    pub cap_face: usize,
    pub cap_surface: usize,
    pub cap_role: CapRole,
    /// Ordered side-chart coordinates, not raw tensor control coordinates.
    /// Endpoint v=0/1 owns the actual cap curve (including declared sewing).
    pub uv: [f64; 2],
    /// Actual outgoing cap EDGE_CURVE/coedge ownership. edge_parameter is the
    /// ordered finite cap-use parameter, not a full CIRCLE support angle.
    pub source_vertex: Option<usize>,
    /// At owned knots the declared POINT is authoritative; curve sewing error
    /// remains in the original source proof, never silently welded here.
    pub edge: usize,
    pub edge_forward: bool,
    pub edge_parameter: [f64; 2],
    /// Fresh actual cap-source evaluation enclosure; the emitted literal point
    /// is its midpoint, exactly as in the original mesher. No inverse fitting.
    pub source_position_mm: [[f64; 2]; 3],
}
/// Ring order is the ORIGINAL side cycle. The start cap's emitted fan uses
/// the reverse cycle, the end cap uses this cycle. The borrowed immutable
/// source binds all IDs/footpoints; no manufacturing composite is admitted.
pub struct TracedMesh<'s> {
    source: &'s RadialBoundary,
    mesh: DisplayMesh,
    rims: [Vec<CapRimFootpoint>; 2],
    work: usize,
}
impl<'s> TracedMesh<'s> {
    pub fn source(&self) -> &'s RadialBoundary {
        self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn cap_rims(&self) -> [&[CapRimFootpoint]; 2] {
        [&self.rims[0], &self.rims[1]]
    }
    /// Actual logical provider charges, including all original source proofs,
    /// literal embedding and fresh trace work. Not process instructions/RSS.
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_metadata_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self
                .rims
                .iter()
                .map(|r| r.capacity() * std::mem::size_of::<CapRimFootpoint>())
                .sum::<usize>()
    }
}
fn capacity(bytes: usize, maximum: usize) -> Result<(), Error> {
    if bytes > maximum {
        Err(Error::Limit)
    } else {
        Ok(())
    }
}
impl RadialBoundary {
    pub fn tessellate_with_cap_trace(
        &self,
        options: Tessellation,
        trace: CapTraceOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<TracedMesh<'_>, Error> {
        self.tessellate_with_cap_trace_budget(options, trace, DEFAULT_LIMITS.work, interrupt)
    }
    pub(crate) fn tessellate_with_cap_trace_budget(
        &self,
        options: Tessellation,
        trace: CapTraceOptions,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<TracedMesh<'_>, Error> {
        if maximum_work > DEFAULT_LIMITS.work || trace.maximum_bytes > BYTES {
            return Err(Error::Invalid);
        }
        capacity(std::mem::size_of::<TracedMesh<'_>>(), trace.maximum_bytes)?;
        build_with_output(
            self,
            options,
            interrupt,
            Limits {
                work: maximum_work,
                ..DEFAULT_LIMITS
            },
            |prepared, grid, ring, mesh, work| {
                finish(
                    self,
                    prepared,
                    grid,
                    ring,
                    mesh,
                    trace,
                    (maximum_work, work),
                )
            },
        )
    }
}
fn finish<'s>(
    source: &'s RadialBoundary,
    prepared: &Prepared<'_>,
    grid: &Grid,
    ring: &Ring,
    mesh: DisplayMesh,
    options: CapTraceOptions,
    budget: (usize, &mut Work<'_>),
) -> Result<TracedMesh<'s>, Error> {
    let (initial_work, work) = budget;
    // counts() already bounds n<=2M. Both exact-capacity requests are checked
    // against the caller's lower budget BEFORE allocation, then actual Vec
    // capacities are checked before filling. No additional grid/mesh candidate.
    let n = ring.entries.len();
    let requested =
        std::mem::size_of::<TracedMesh<'_>>() + 2 * n * std::mem::size_of::<CapRimFootpoint>();
    capacity(requested, options.maximum_bytes)?;
    work.charge(2 * n)?;
    let mut result = TracedMesh {
        source,
        mesh,
        rims: std::array::from_fn(|_| Vec::with_capacity(n)),
        work: 0,
    };
    capacity(result.retained_metadata_bytes(), options.maximum_bytes)?;
    for end in 0..2 {
        let row = if end == 0 { 0 } else { grid.v.len() - 1 };
        let cap = &prepared.caps[end];
        for (index, &(side, u_index)) in ring.entries.iter().enumerate() {
            work.charge(64)?;
            let u = grid.u[side][u_index];
            let bounds = cap.point(source, end, side, u, work)?;
            let group = &source.ordered.caps.rims[end][side];
            let (segment, t) = if group.len() == 1 {
                (0, I::point(u))
            } else if u >= 0.5 {
                (1, I::point(2.) * I::point(u) - I::point(1.))
            } else {
                (0, I::point(2.) * I::point(u))
            };
            let (edge, edge_forward) = group[segment];
            let side_face = source.side_faces()[side];
            result.rims[end].push(CapRimFootpoint {
                vertex: (row * n + index) as u32,
                source_vertex: cap.source_vertex(source, end, side, u),
                side,
                side_face,
                side_surface: source.input().faces[side_face].surface,
                cap_face: cap.face,
                cap_surface: source.input().faces[cap.face].surface,
                cap_role: if end == 0 {
                    CapRole::AxialStart
                } else {
                    CapRole::AxialEnd
                },
                uv: [u, end as f64],
                edge,
                edge_forward,
                edge_parameter: [t.lo, t.hi],
                source_position_mm: bounds.map(|p| [p.lo, p.hi]),
            });
        }
    }
    work.step()?;
    result.work = initial_work - work.remaining;
    Ok(result)
}
#[cfg(test)]
mod tests;
