//! Sequential validator allocation requests, not allocator overhead or RSS.
use super::*;

// Audited against pinned rustc 58cf8f9c3 alloc/collections/btree/node.rs and
// alloc/raw_vec/mod.rs. BTree nodes have eleven key/value slots and twelve
// https://github.com/rust-lang/rust/blob/58cf8f9c35b355199727627d07ae34a4e2846d67/library/alloc/src/collections/btree/node.rs
// https://github.com/rust-lang/rust/blob/58cf8f9c35b355199727627d07ae34a4e2846d67/library/alloc/src/raw_vec/mod.rs
// child pointers. The largest embedding map node (u32 -> Vec<usize>) fits
// 512 bytes on the supported 32/64-bit targets. Each insertion-only tree has
// at most T+1 nodes, including the transient split/root allowance. Charging
// all groups/rim/seen trees together takes <1536T+fixed bytes; growing group
// index buffers, including old+new allocation requests, take <48T bytes.
//
// After those trees drop, embedding retains 16T fan bytes, 48T boxes, 8T IDs
// and 2T BVH nodes. A Node has only two 3D points, two optional indices, one
// index and an optional fan ID; even padding each field separately fits
// 256 bytes. Median splits give depth <=20 for T<=1M. Pair DFS descends one
// or both trees each step, so depth <=41 and pending <=83 (three children
// add at most two pending siblings). Its growth, recursive frames and small
// per-pair vectors fit the fixed MiB allowance. The manifold validator's
// 3T edge triples, 3T fan parents, T shell parents and V links take 52T+4V.
// Summing these sequential phases, rather than relying on their disjoint
// lifetimes, fits the deliberately broad 4096T+16V+1MiB charge below.
// Input meshes/source/cache are charged separately by the consumer. This
// does not claim bounds on allocator bookkeeping, fragmentation or RSS.
pub(super) fn validators(vertices: usize, triangles: usize) -> Result<usize, GeometryError> {
    triangles
        .checked_mul(4096)
        .and_then(|n| vertices.checked_mul(16).and_then(|v| n.checked_add(v)))
        .and_then(|n| n.checked_add(1024 * 1024))
        .ok_or_else(|| {
            failure(
                GeometryErrorKind::Limit,
                "Native validator storage bound overflow",
            )
        })
}
pub(super) fn combined(
    producer: usize,
    consumer: usize,
    validator: usize,
) -> Result<usize, GeometryError> {
    producer
        .checked_add(consumer)
        .and_then(|n| n.checked_add(validator))
        .ok_or_else(|| {
            failure(
                GeometryErrorKind::Limit,
                "Combined native storage bound overflow",
            )
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_validator_counts_have_checked_separate_workspace_bounds() {
        let actual = validators(39_018, 78_032).unwrap();
        assert!(actual > 78_032 * 4096);
        assert!(actual > 65_536 + 78_032 * 52 + 39_018 * 4);
        assert!(validators(500_002, 1_000_000).unwrap() > actual);
        assert_eq!(combined(1024, 2048, actual).unwrap(), actual + 3072);
        assert_eq!(
            combined(usize::MAX, 1, actual).unwrap_err().kind,
            GeometryErrorKind::Limit
        );
        assert_eq!(
            combined(1024, 2048, usize::MAX).unwrap_err().kind,
            GeometryErrorKind::Limit
        );
        assert_eq!(
            validators(usize::MAX, 1).unwrap_err().kind,
            GeometryErrorKind::Limit
        );
        assert_eq!(
            validators(1, usize::MAX).unwrap_err().kind,
            GeometryErrorKind::Limit
        );
    }
}
