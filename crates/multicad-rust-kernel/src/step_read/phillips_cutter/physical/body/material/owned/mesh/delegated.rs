//! Conservative sequential provider workspace charges, excluding allocator
//! fragmentation/RSS and unrelated allocations owned by the caller.
use super::*;
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DelegatedWorkspaceBounds {
    pub retained_thread_bytes: usize,
    pub cutter_bytes: usize,
    pub planar_bytes: usize,
    pub connectivity_bytes: usize,
}
fn bound(base: usize, terms: &[(usize, usize)]) -> R<usize> {
    terms
        .iter()
        .try_fold(base, |sum, &(count, size)| {
            sum.checked_add(count.checked_mul(size)?)
        })
        .ok_or(PhillipsReadError::Limit)
}
impl DelegatedWorkspaceBounds {
    pub(super) fn maximum(self) -> usize {
        self.retained_thread_bytes
            .max(self.cutter_bytes)
            .max(self.planar_bytes)
            .max(self.connectivity_bytes)
    }
}
// Prepared controls: degree9 position and degree17 normal grids/row headers,
// optional fixed NormalJet, knots, rails and caps fit1024 bytes/input pole.
// The fixed512MiB additionally covers both128-per-side strip tables (<=32
// sides),131072-point/65536-cell BTree caches, prospective grid/refinement
// sets,32MiB parallel records,8x1MiB stacks and bounded local polynomial scratch.
// Product-grid output has T/2+2 vertices.144 bytes/T covers simultaneous
// output, cap trace metadata, connectivity and proportional grid buffers.
pub(super) fn retained(source: &RadialBoundary, triangles: usize) -> R<usize> {
    let poles = source
        .input()
        .surfaces
        .iter()
        .try_fold(0usize, |sum, s| {
            sum.checked_add(match s {
                crate::radial_tensor::Surface::Tensor(t) => t.poles().len(),
                _ => 0,
            })
        })
        .ok_or(PhillipsReadError::Limit)?;
    bound(512 * 1024 * 1024, &[(poles, 1024), (triangles, 144)])
}
// Raw quad budget precedes construction: at most2T Builder points. The
// PointKey/usize HashMap uses40-byte entries, including pinned hashbrown's
// <=7/8 load and old+new tables during growth. Growing coordinate/index/face
// vectors, trace/wire descriptors and connectivity remain below512 bytes/T;
// one MiB covers the sealed fixed atlas/schema and fixed query scratch.
pub(super) fn cutter(triangles: usize) -> R<usize> {
    bound(1024 * 1024, &[(triangles, 512)])
}
// Pinned Spade2.15.1 try_bulk_load_cdt inserts only the n original vertices.
// Its planar DCEL has <=3n undirected edges/2n faces. Dedup/index vectors,
// hull buckets, skipped sets, one-constraint conflict/legality paths and the
// final flood/cap buffers are linear live storage.8192 bytes/input point
// conservatively includes their payloads, growth and nested Vec headers.
pub(super) fn planar(points: usize) -> R<usize> {
    bound(1024 * 1024 + 8192 * 8, &[(points, 8192)])
}
// certify():3T edge triples(36T),3T fan parents(12T),T shell parents(4T).
// Vertex links(4V) are allocated after edge triples drop; counting both is
// conservative. Input mesh/owner storage is charged separately by Work.
pub(super) fn connectivity(vertices: usize, triangles: usize) -> R<usize> {
    bound(65536, &[(triangles, 52), (vertices, 4)])
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn checked_provider_charges_scale_and_refuse_unrepresentable_counts() {
        let actual = super::super::super::mesh_probe::owned(1., false);
        assert!(retained(actual.thread(), 1_000_000).unwrap() > 512 * 1024 * 1024);
        assert!(cutter(1_000_000).unwrap() > cutter(100_000).unwrap());
        assert!(planar(1184).unwrap() > planar(1128).unwrap());
        assert!(connectivity(39_018, 78_032).unwrap() > 78_032 * 52);
        assert!(matches!(cutter(usize::MAX), Err(PhillipsReadError::Limit)));
        assert!(matches!(planar(usize::MAX), Err(PhillipsReadError::Limit)));
        assert!(matches!(
            connectivity(usize::MAX, 78_032),
            Err(PhillipsReadError::Limit)
        ));
    }
}
