//! Full actual mouth wire from oriented source refinement incidence, never XY fitting.
use super::*;
const BYTES: usize = 32 * 1024 * 1024;
const WORK: usize = 512_000_000;
const SCRATCH: usize = 16384;
#[derive(Clone, Copy)]
pub struct MouthTraceOptions {
    /// New trace bookkeeping/receipt peak only. Existing source and mesh
    /// storage retain their old contracts; this is not whole-request RSS.
    pub maximum_bytes: usize,
    /// New traced-path logical source/refinement query ledger. The original
    /// untraced mesher has no numeric work ledger and keeps its existing limits.
    pub maximum_work: usize,
}
impl Default for MouthTraceOptions {
    fn default() -> Self {
        Self {
            maximum_bytes: BYTES,
            maximum_work: WORK,
        }
    }
}
/// Original source-keyed point on an outgoing finite boundary curve. At a
/// source Anchor, its owned node is authoritative (not a proximity weld).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MouthFootpoint {
    pub vertex: u32,
    pub source_point: MeshPointKey,
    pub curve: MeshCurveKey,
    pub curve_forward: bool,
    pub parameter: f64,
    pub source_position_mm: [[f64; 2]; 3],
}
/// Full shaped source-mouth cycle, in the original Cutter cap winding.
/// This retains the unchanged source-certified boundary mesh, not a globally
/// embedded material/manufacturing claim. Material must reverse Cutter sides.
pub struct TracedMouthMesh<'s> {
    source: &'s Cutter,
    boundary: BoundaryMesh,
    mouth: Vec<MouthFootpoint>,
    work: usize,
    peak: usize,
}
impl<'s> TracedMouthMesh<'s> {
    pub fn source(&self) -> &'s Cutter {
        self.source
    }
    pub fn boundary_mesh(&self) -> &BoundaryMesh {
        &self.boundary
    }
    pub fn mouth(&self) -> &[MouthFootpoint] {
        &self.mouth
    }
    /// Actual proposed cells (including skipped cells), prospective normal
    /// proof quads, fixed graph units, connectivity index units and trace work.
    /// This is a logical query budget, not an instruction count or timing.
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn peak_metadata_bytes(&self) -> usize {
        self.peak
    }
    pub fn retained_metadata_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + self.mouth.capacity() * std::mem::size_of::<MouthFootpoint>()
    }
}
struct Counted<'a> {
    remaining: usize,
    maximum: usize,
    bytes: usize,
    peak: usize,
    cancel: &'a dyn Fn() -> bool,
}
impl Counted<'_> {
    fn memory(&mut self, bytes: usize) -> Result<(), Error> {
        self.charge(0)?;
        capacity(bytes, self.bytes)?;
        self.peak = self.peak.max(bytes);
        Ok(())
    }
}
impl Progress for Counted<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        checkpoint(self.cancel)?;
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
fn capacity(n: usize, maximum: usize) -> Result<(), Error> {
    if n > maximum {
        Err(Error::Limit)
    } else {
        Ok(())
    }
}
fn ordered_wire(
    ends: &[[usize; 2]],
    uses: &[(usize, bool)],
    work: &mut Counted<'_>,
    prefix: usize,
) -> Result<Vec<(usize, bool)>, Error> {
    work.memory(prefix + 2 * uses.len() * std::mem::size_of::<(usize, bool)>())?;
    let mut pending = uses.to_vec();
    let mut ordered = Vec::with_capacity(uses.len());
    work.memory(
        prefix + (pending.capacity() + ordered.capacity()) * std::mem::size_of::<(usize, bool)>(),
    )?;
    let first = *pending.first().ok_or(Error::Unresolved)?;
    let start = ends[first.0][usize::from(!first.1)];
    let mut node = start;
    while !pending.is_empty() {
        work.charge(pending.len())?;
        let mut next = None;
        for (i, &(edge, forward)) in pending.iter().enumerate() {
            if ends[edge][usize::from(!forward)] == node && next.replace(i).is_some() {
                return Err(Error::Unresolved);
            }
        }
        let i = next.ok_or(Error::Unresolved)?;
        let use_ = pending.remove(i);
        node = ends[use_.0][usize::from(use_.1)];
        ordered.push(use_);
        if node == start && !pending.is_empty() {
            return Err(Error::Unresolved);
        }
    }
    if node != start {
        return Err(Error::Unresolved);
    }
    Ok(ordered)
}
impl<'s> PatchAtlas<'s> {
    /// Add source-mouth ownership without changing the original emitted mesh.
    /// The traced path alone counts bounded logical source-query work; the
    /// original path retains its triangle/cancellation contract unchanged.
    pub fn boundary_mesh_with_mouth_trace(
        &self,
        options: Tessellation,
        trace: MouthTraceOptions,
        cancel: &dyn Fn() -> bool,
    ) -> Result<TracedMouthMesh<'s>, Error> {
        if trace.maximum_bytes > BYTES || trace.maximum_work > WORK {
            return Err(Error::Invalid);
        }
        let mut work = Counted {
            remaining: trace.maximum_work,
            maximum: trace.maximum_work,
            bytes: trace.maximum_bytes,
            peak: 0,
            cancel,
        };
        work.memory(std::mem::size_of::<TracedMouthMesh<'_>>() + SCRATCH)?;
        self.refine_with_output(
            options,
            cancel,
            [1; 8],
            &mut work,
            |builder, connectivity, distance, normals, work| {
                finish(builder, connectivity, distance, normals, work)
            },
        )
    }
}
fn finish<'s>(
    builder: Builder<'_, 's>,
    connectivity: MeshProof,
    distance: f64,
    normals: f64,
    work: &mut Counted<'_>,
) -> Result<TracedMouthMesh<'s>, Error> {
    let edges = &builder.complex.edges;
    let nodes = builder.complex.points.len();
    let n = edges.len();
    // Scratch owns signed uses, ordered/pending pairs, endpoint/rank tables and
    // candidate sort records. All sizes derive from already bounded live graph.
    let fixed = SCRATCH
        + std::mem::size_of::<TracedMouthMesh<'_>>()
        + n * (1 + 4 * std::mem::size_of::<(usize, bool)>() + 3 * std::mem::size_of::<usize>())
        + nodes * std::mem::size_of::<usize>();
    work.memory(fixed)?;
    work.charge(n + nodes + 1)?;
    let mut sums = vec![0i8; n];
    for (index, (patch, side)) in builder
        .atlas
        .patches
        .iter()
        .zip(&builder.complex.sides)
        .enumerate()
    {
        if builder.atlas.topology.faces()[patch.face].source_identity() != FaceKey::Mouth {
            continue;
        }
        work.charge(4)?;
        let reverse = builder.atlas.patch_winding(index)? < 0;
        for (index, forward) in [(2, true), (1, true), (3, false), (0, false)] {
            if let Side::Curve(edge) = side[index] {
                sums[edge] += if forward != reverse { 1 } else { -1 };
            }
        }
    }
    // Complex::build already certified two opposite uses globally. Summation
    // cancels only identical source-owned internal edges of the Mouth patches.
    let uses = sums
        .iter()
        .enumerate()
        .filter_map(|(edge, &s)| (s != 0).then_some((edge, s > 0)))
        .collect::<Vec<_>>();
    let ends = edges.iter().map(|e| e.ends).collect::<Vec<_>>();
    let before_wire = SCRATCH
        + std::mem::size_of::<TracedMouthMesh<'_>>()
        + sums.capacity()
        + uses.capacity() * std::mem::size_of::<(usize, bool)>()
        + ends.capacity() * std::mem::size_of::<[usize; 2]>();
    let ordered = ordered_wire(&ends, &uses, work, before_wire)?;
    let mut rank = vec![usize::MAX; n];
    let mut node_rank = vec![usize::MAX; nodes];
    for (i, &(edge, forward)) in ordered.iter().enumerate() {
        rank[edge] = i;
        node_rank[ends[edge][usize::from(!forward)]] = i;
    }
    let owner = |key: PointKey| -> Option<(usize, f64)> {
        match key {
            PointKey::Anchor(node) if node_rank[node] != usize::MAX => {
                let r = node_rank[node];
                Some((r, if ordered[r].1 { 0. } else { 1. }))
            }
            PointKey::Edge(edge, t) if rank[edge] != usize::MAX => {
                Some((rank[edge], f64::from_bits(t)))
            }
            _ => None,
        }
    };
    work.charge(builder.points.len())?;
    let count = builder
        .points
        .keys()
        .filter(|&&key| owner(key).is_some())
        .count();
    let workspace = before_wire
        + ordered.capacity() * std::mem::size_of::<(usize, bool)>()
        + rank.capacity() * std::mem::size_of::<usize>()
        + node_rank.capacity() * std::mem::size_of::<usize>();
    let payload = count
        * (std::mem::size_of::<MouthFootpoint>()
            + std::mem::size_of::<(usize, f64, PointKey, usize)>());
    work.memory(workspace + payload)?;
    work.charge(
        builder.points.len()
            + count * (usize::BITS as usize - count.max(1).leading_zeros() as usize + 1),
    )?;
    let mut candidates = Vec::with_capacity(count);
    for (&key, &vertex) in &builder.points {
        if let Some((r, t)) = owner(key) {
            candidates.push((r, t, key, vertex));
        }
    }
    candidates.sort_unstable_by(|a, b| {
        a.0.cmp(&b.0).then_with(|| {
            if ordered[a.0].1 {
                a.1.total_cmp(&b.1)
            } else {
                b.1.total_cmp(&a.1)
            }
        })
    });
    let mut mouth = Vec::with_capacity(count);
    work.memory(
        workspace
            + candidates.capacity() * std::mem::size_of::<(usize, f64, PointKey, usize)>()
            + mouth.capacity() * std::mem::size_of::<MouthFootpoint>(),
    )?;
    for (r, t, key, vertex) in candidates.iter().copied() {
        work.charge(1)?;
        let (edge, forward) = ordered[r];
        let p = match key {
            PointKey::Anchor(node) => builder.complex.points[node],
            _ => builder.complex.edge_point(builder.atlas, edge, t)?,
        };
        mouth.push(MouthFootpoint {
            vertex: vertex as u32,
            source_point: key,
            curve: edges[edge].key,
            curve_forward: forward,
            parameter: t,
            source_position_mm: p.map(|p| [p.lo, p.hi]),
        });
    }
    let result = TracedMouthMesh {
        source: builder.atlas.topology.events().source(),
        boundary: BoundaryMesh {
            mesh: builder.mesh,
            connectivity,
            source_distance_mm: distance,
            normal_chord: normals,
            refinement: builder.complex.certificate,
        },
        mouth,
        work: 0,
        peak: 0,
    };
    // Account actual new capacities while graph/sort workspaces still live.
    let actual = workspace
        + result.mouth.capacity() * std::mem::size_of::<MouthFootpoint>()
        + candidates.capacity() * std::mem::size_of::<(usize, f64, PointKey, usize)>();
    work.memory(actual)?;
    publish(result, work)
}
fn publish<'s>(
    mut result: TracedMouthMesh<'s>,
    work: &mut Counted<'_>,
) -> Result<TracedMouthMesh<'s>, Error> {
    work.charge(1)?;
    result.work = work.maximum - work.remaining;
    result.peak = work.peak;
    Ok(result)
}
#[cfg(test)]
mod tests;
