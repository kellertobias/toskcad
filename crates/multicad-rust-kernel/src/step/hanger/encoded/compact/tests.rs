use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
#[test]
fn actual_compact_cartesian_surfaces_reproduce_all_original_homogeneous_coefficients() {
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
        let mut w = Work::new(16_000_000, None).unwrap();
        w.charge(1 + field.work_used()).unwrap();
        let mut m = Memory::new(field.retained_bytes(), 64 * 1024 * 1024).unwrap();
        m.add(SCRATCH).unwrap();
        let result = build(&field, &mut m, &mut w);
        eprintln!(
            "actual compact producer left={left} nb={nb} starts={metric_starts}/{wood_starts} work={} result={:?}",
            w.used,
            result.as_ref().map(|r| (r.surfaces.len(), r.bytes))
        );
        let compact = result.unwrap();
        assert!(std::ptr::eq(compact.field(), &field));
        assert_eq!(compact.cells().len(), field.cells.len());
        assert_eq!(compact.work_used(), w.used);
        assert!(compact.retained_bytes() > 0);
        let per_turn = field.angular.coefficients().len() * 10;
        let mut count = 0;
        for (index, cell) in field.cells.iter().enumerate() {
            let map = &compact.cells()[index];
            let surface = &compact.surfaces()[map.surface()];
            assert!(std::ptr::eq(surface.field(), &field));
            assert_eq!(surface.meridian(), cell.meridian);
            let starts = if field.meridians[cell.meridian].component == Component::Metric {
                metric_starts
            } else {
                wood_starts
            } as i32;
            assert_eq!(
                surface.strand(),
                cell.source_key.total_period.rem_euclid(starts)
            );
            assert_eq!(
                map.turn_shift(),
                i64::from(cell.source_key.total_period.div_euclid(starts))
                    - i64::from(surface.first_turn())
            );
            let dimensions = surface.dimensions();
            assert_eq!(surface.degrees(), [dimensions[0] - 1, 10]);
            assert_eq!(surface.points().len(), dimensions[0] * dimensions[1]);
            assert_eq!(surface.row_weights().len(), dimensions[0]);
            assert_eq!(
                surface.source_bounds(),
                field.meridians[cell.meridian].bounds
            );
            assert!(surface.knots().windows(2).all(|p| p[0] < p[1]));
            let mut oracle = Work::new(16_000_000, None).unwrap();
            for u in 0..dimensions[0] {
                for piece in 0..field.angular.coefficients().len() {
                    for v in 0..=10 {
                        let pole = map.turn_shift() as usize * per_turn + piece * 10 + v;
                        let actual = &surface.points()[u * dimensions[1] + pole];
                        let h = field
                            .homogeneous_control(index, u, piece, v, &mut oracle)
                            .unwrap()
                            .1;
                        let weight = surface.row_weights()[u].number(&mut oracle).unwrap();
                        assert_eq!(weight.0, h[3].0);
                        for axis in 0..3 {
                            assert_eq!(
                                actual[axis]
                                    .number(&mut oracle)
                                    .unwrap()
                                    .mul(&weight, &mut oracle)
                                    .unwrap()
                                    .0,
                                h[axis].0
                            )
                        }
                        count += 1;
                    }
                }
            }
        }
        assert!(count > 30_000);
        let used = compact.work_used();
        assert!(
            CompactSurfaces::prepare(
                &field,
                DomainOptions {
                    maximum_work: used,
                    ..Default::default()
                },
                None
            )
            .is_ok()
        );
        assert!(matches!(
            CompactSurfaces::prepare(
                &field,
                DomainOptions {
                    maximum_work: used - 1,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            CompactSurfaces::prepare(
                &field,
                DomainOptions {
                    maximum_bytes: field.retained_bytes(),
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            CompactSurfaces::prepare(&field, DomainOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
    }
}
#[test]
fn compact_surface_topology_cancels_internal_meridians_and_preserves_oriented_shell() {
    use std::collections::{BTreeMap, BTreeSet};
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
        let closed = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
        let compact = CompactSurfaces::prepare(&field, DomainOptions::default(), None).unwrap();
        let mut face_map = Vec::new();
        let mut nfaces = compact.surfaces().len();
        for face in closed.faces() {
            face_map.push(match face.kind {
                ClosedFaceKind::ThreadCell(i) => {
                    compact.cells()[domains.cells()[i].source_cell].surface()
                }
                _ => {
                    let id = nfaces;
                    nfaces += 1;
                    id
                }
            })
        }
        let mut uses = vec![Vec::new(); closed.edges().len()];
        for (i, face) in closed.faces().iter().enumerate() {
            for c in face.loops().iter().flatten() {
                uses[c.edge].push(i)
            }
        }
        assert!(uses.iter().all(|x| x.len() == 2));
        let mut groups: BTreeMap<[usize; 2], Vec<usize>> = BTreeMap::new();
        let mut map = vec![None; closed.edges().len()];
        let mut edges: Vec<[usize; 2]> = Vec::new();
        let mut removed = 0;
        for (i, edge) in closed.edges().iter().enumerate() {
            let mapped = uses[i].iter().map(|&f| face_map[f]).collect::<Vec<_>>();
            if mapped[0] == mapped[1] {
                assert!(matches!(edge.curve, ClosedCurve::Thread(_)));
                removed += 1;
                continue;
            }
            let a = uses[i][0];
            let b = uses[i][1];
            if matches!(closed.faces()[a].kind, ClosedFaceKind::ThreadCell(_))
                && matches!(closed.faces()[b].kind, ClosedFaceKind::ThreadCell(_))
            {
                groups
                    .entry([mapped[0].min(mapped[1]), mapped[0].max(mapped[1])])
                    .or_default()
                    .push(i)
            } else {
                map[i] = Some((edges.len(), true));
                edges.push(edge.vertices)
            }
        }
        for list in groups.values() {
            let mut degrees: BTreeMap<usize, usize> = BTreeMap::new();
            for &i in list {
                for v in closed.edges()[i].vertices {
                    *degrees.entry(v).or_default() += 1
                }
            }
            let ends = degrees
                .iter()
                .filter(|(_, n)| **n == 1)
                .map(|(&v, _)| v)
                .collect::<Vec<_>>();
            assert_eq!(ends.len(), 2);
            assert!(degrees.values().all(|&n| n == 1 || n == 2));
            let target = edges.len();
            edges.push([ends[0], ends[1]]);
            let mut current = ends[0];
            let mut unseen = list.iter().copied().collect::<BTreeSet<_>>();
            while current != ends[1] {
                let edge = *unseen
                    .iter()
                    .find(|&&i| closed.edges()[i].vertices.contains(&current))
                    .unwrap();
                unseen.remove(&edge);
                let old = closed.edges()[edge].vertices;
                let forward = old[0] == current;
                map[edge] = Some((target, forward));
                current = if forward { old[1] } else { old[0] };
            }
            assert!(unseen.is_empty());
        }
        let mut faces = vec![BTreeMap::new(); nfaces];
        for (i, face) in closed.faces().iter().enumerate() {
            for c in face.loops().iter().flatten() {
                if let Some((edge, forward)) = map[c.edge] {
                    let sense = c.sense == forward;
                    if let Some(previous) = faces[face_map[i]].insert(edge, sense) {
                        assert_eq!(previous, sense)
                    }
                }
            }
        }
        let directed = |e: usize, sense: bool| {
            if sense {
                edges[e]
            } else {
                [edges[e][1], edges[e][0]]
            }
        };
        let mut chi_faces = 0isize;
        let mut edge_uses = vec![Vec::new(); edges.len()];
        let mut vertex_set = BTreeSet::new();
        for (face_id, face) in faces.iter().enumerate() {
            let mut remaining = face.clone();
            let mut loops = 0;
            while let Some((&first, &sense)) = remaining.first_key_value() {
                remaining.remove(&first);
                let [start, mut current] = directed(first, sense);
                edge_uses[first].push((face_id, sense));
                vertex_set.insert(start);
                vertex_set.insert(current);
                while current != start {
                    let (&next, &ns) = remaining
                        .iter()
                        .find(|(e, s)| directed(**e, **s)[0] == current)
                        .unwrap();
                    remaining.remove(&next);
                    let ends = directed(next, ns);
                    edge_uses[next].push((face_id, ns));
                    vertex_set.insert(ends[1]);
                    current = ends[1];
                }
                loops += 1;
            }
            assert!(loops >= 1);
            chi_faces += 2 - loops;
        }
        assert!(edge_uses.iter().all(|u| u.len() == 2 && u[0].1 != u[1].1));
        assert_eq!(
            vertex_set.len() as isize - edges.len() as isize + chi_faces,
            2
        );
        eprintln!(
            "compact shell left={left} nb={nb} faces={nfaces} edges={} vertices={} side_chains={} removed_internal={removed}",
            edges.len(),
            vertex_set.len(),
            groups.len()
        );
    }
}
