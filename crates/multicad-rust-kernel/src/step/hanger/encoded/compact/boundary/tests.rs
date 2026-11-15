use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
#[test]
fn actual_compact_boundary_covers_original_owner_arenas_and_one_oriented_shell() {
    for (left, nb, metric_starts, wood_starts) in [
        (false, false, 1, 1),
        (true, true, 1, 1),
        (false, true, 2, 3),
    ] {
        let mut spec = crate::authored_hanger_bolt::tests::fixture(left, nb);
        spec.metric = crate::authored_metric_thread::AuthoredMetricThread::new(
            crate::authored_metric_thread::Spec {
                starts: metric_starts,
                ..spec.metric.spec()
            },
        )
        .unwrap();
        spec.wood = crate::authored_flat_thread::AuthoredThread::new(
            crate::helical_thread::Spec {
                starts: wood_starts,
                ..spec.wood.spec()
            },
            crate::authored_flat_thread::SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-6,
            },
        )
        .unwrap();
        let source = Hanger::new(spec, None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
        let closed = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
        let surfaces = CompactSurfaces::prepare(&field, DomainOptions::default(), None).unwrap();
        let compact =
            CompactBoundary::prepare(&closed, &surfaces, DomainOptions::default(), None).unwrap();
        assert!(std::ptr::eq(compact.source(), &closed));
        assert!(std::ptr::eq(compact.surfaces(), &surfaces));
        assert!(compact.retained_bytes() > 0);
        assert_eq!(compact.faces().len(), surfaces.surfaces().len() + 8);
        assert_eq!(compact.original_edges().len(), closed.edges().len());
        assert_eq!(compact.original_faces().len(), closed.faces().len());
        let mut uses = vec![Vec::new(); compact.edges().len()];
        let mut vertices = BTreeSet::new();
        let mut chi = 0isize;
        for (i, face) in compact.faces().iter().enumerate() {
            assert!(!face.loops().is_empty());
            chi += 2 - face.loops().len() as isize;
            for cycle in face.loops() {
                let mut endpoints = Vec::new();
                for c in cycle {
                    let e = &compact.edges()[c.edge];
                    let pair = if c.sense {
                        e.vertices
                    } else {
                        [e.vertices[1], e.vertices[0]]
                    };
                    vertices.extend(pair);
                    endpoints.push(pair);
                    uses[c.edge].push((i, c.sense))
                }
                for j in 0..endpoints.len() {
                    assert_eq!(endpoints[j][1], endpoints[(j + 1) % endpoints.len()][0])
                }
            }
        }
        assert!(uses.iter().all(|u| u.len() == 2 && u[0].1 != u[1].1));
        assert_eq!(
            vertices.len() as isize - compact.edges().len() as isize + chi,
            2
        );
        let mut reached = BTreeSet::from([0]);
        loop {
            let before = reached.len();
            for pair in &uses {
                if reached.contains(&pair[0].0) || reached.contains(&pair[1].0) {
                    reached.insert(pair[0].0);
                    reached.insert(pair[1].0);
                }
            }
            if reached.len() == before {
                break;
            }
        }
        assert_eq!(reached.len(), compact.faces().len());
        for (i, face) in closed.faces().iter().enumerate() {
            let mapped = compact.original_faces()[i];
            assert_eq!(compact.faces()[mapped].same_sense, face.same_sense);
            let expected = match face.kind {
                ClosedFaceKind::ThreadCell(cell) => ChartFaceKind::Thread(
                    surfaces.cells()[domains.cells()[cell].source_cell].surface(),
                ),
                ClosedFaceKind::PlainWall { quarter } => ChartFaceKind::PlainWall { quarter },
                ClosedFaceKind::Plane(key) => ChartFaceKind::Plane(key),
            };
            assert_eq!(compact.faces()[mapped].kind, expected);
        }
        for (i, edge) in compact.edges().iter().enumerate() {
            match edge.kind {
                ChartEdgeKind::Original(original) => {
                    assert_eq!(edge.vertices, closed.edges()[original].vertices);
                    assert_eq!(
                        compact.original_edges()[original],
                        Some(ChartCoedge {
                            edge: i,
                            sense: true
                        })
                    );
                }
                ChartEdgeKind::Helix { source_edge } => {
                    assert_eq!(compact.original_edges()[source_edge].unwrap().edge, i);
                    let mut next = BTreeMap::new();
                    let mut members = 0;
                    for (old, map) in compact.original_edges().iter().enumerate() {
                        if let Some(c) = map
                            && c.edge == i
                        {
                            let pair = closed.edges()[old].vertices;
                            let pair = if c.sense { pair } else { [pair[1], pair[0]] };
                            assert!(next.insert(pair[0], pair[1]).is_none());
                            members += 1;
                        }
                    }
                    let mut current = edge.vertices[0];
                    let mut count = 0;
                    while current != edge.vertices[1] {
                        current = next.remove(&current).unwrap();
                        count += 1;
                    }
                    assert_eq!(count, members);
                    assert!(next.is_empty());
                }
            }
        }
        let mut old_uses = vec![Vec::new(); closed.edges().len()];
        for (i, face) in closed.faces().iter().enumerate() {
            for c in face.loops().iter().flatten() {
                old_uses[c.edge].push(compact.original_faces()[i])
            }
        }
        for (i, pair) in old_uses.iter().enumerate() {
            assert_eq!(compact.original_edges()[i].is_none(), pair[0] == pair[1]);
        }
        let used = compact.work_used();
        assert!(
            CompactBoundary::prepare(
                &closed,
                &surfaces,
                DomainOptions {
                    maximum_work: used,
                    ..Default::default()
                },
                None
            )
            .is_ok()
        );
        assert!(matches!(
            CompactBoundary::prepare(
                &closed,
                &surfaces,
                DomainOptions {
                    maximum_work: used - 1,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            CompactBoundary::prepare(
                &closed,
                &surfaces,
                DomainOptions {
                    maximum_bytes: closed.retained_bytes() + surfaces.retained_bytes(),
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            CompactBoundary::prepare(&closed, &surfaces, DomainOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        let other = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let other = CompactSurfaces::prepare(&other, DomainOptions::default(), None).unwrap();
        assert!(matches!(
            CompactBoundary::prepare(&closed, &other, DomainOptions::default(), None),
            Err(Error::Invalid)
        ));
        eprintln!(
            "compact boundary left={left} nb={nb} starts={metric_starts}/{wood_starts} faces={} edges={} vertices={} work={used} bytes={}",
            compact.faces().len(),
            compact.edges().len(),
            vertices.len(),
            compact.retained_bytes()
        );
    }
}
