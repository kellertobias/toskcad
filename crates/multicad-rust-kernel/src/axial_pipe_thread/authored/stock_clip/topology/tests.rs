use super::*;
use std::collections::{BTreeMap, BTreeSet};
#[test]
fn actual_physical_source_graph_has_closed_opposite_incidence_and_manifold_vertex_links() {
    for quarter in [false, true] {
        for scale in [0.1, 1., 3.] {
            let source = super::super::tests::canonical(quarter, scale);
            let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
            let mut uses: BTreeMap<StockClipCurve, Vec<(StockClipFace, bool)>> = BTreeMap::new();
            let mut links: BTreeMap<StockClipVertex, Vec<[StockClipCurve; 2]>> = BTreeMap::new();
            for face in atlas.physical_faces() {
                let wire = atlas.face_loop(face);
                assert_eq!(wire.face(), face);
                let edges = wire.coedges();
                for (i, e) in edges.iter().enumerate() {
                    uses.entry(e.curve).or_default().push((face, e.forward));
                    let mut ends = e.curve.endpoints();
                    if !e.forward {
                        ends.swap(0, 1);
                    }
                    let next = edges[(i + 1) % edges.len()];
                    let mut next_ends = next.curve.endpoints();
                    if !next.forward {
                        next_ends.swap(0, 1);
                    }
                    assert_eq!(ends[1], next_ends[0]);
                    links
                        .entry(ends[1])
                        .or_default()
                        .push([e.curve, next.curve]);
                }
            }
            assert_eq!(
                uses.keys().copied().collect::<BTreeSet<_>>(),
                atlas.physical_curves().into_iter().collect()
            );
            assert_eq!(uses.len(), 17);
            assert_eq!(links.len(), 12);
            assert_eq!(12 - 17 + 7, 2);
            for (edge, owners) in &uses {
                assert_eq!(owners.len(), 2);
                assert_ne!(owners[0].1, owners[1].1);
                assert_eq!(
                    owners.iter().map(|x| x.0).collect::<BTreeSet<_>>(),
                    atlas.curve_faces(*edge).into_iter().collect()
                );
            }
            for cycles in links.values() {
                let mut adjacency: BTreeMap<StockClipCurve, Vec<StockClipCurve>> = BTreeMap::new();
                for [a, b] in cycles {
                    adjacency.entry(*a).or_default().push(*b);
                    adjacency.entry(*b).or_default().push(*a);
                }
                assert!(adjacency.values().all(|v| v.len() == 2));
                let mut stack = vec![*adjacency.keys().next().unwrap()];
                let mut seen = BTreeSet::new();
                while let Some(e) = stack.pop() {
                    if seen.insert(e) {
                        stack.extend(&adjacency[&e]);
                    }
                }
                assert_eq!(seen.len(), adjacency.len());
            }
            let mut seen = BTreeSet::new();
            let mut stack = vec![atlas.physical_faces()[0]];
            while let Some(face) = stack.pop() {
                if seen.insert(face) {
                    for owners in uses.values() {
                        if owners.iter().any(|x| x.0 == face) {
                            stack.extend(owners.iter().map(|x| x.0));
                        }
                    }
                }
            }
            assert_eq!(seen.len(), 7);
            assert_eq!(*atlas.source(), source);
        }
    }
}
