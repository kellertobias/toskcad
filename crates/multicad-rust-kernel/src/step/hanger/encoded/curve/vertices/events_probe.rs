//! Actual finite cap-event owner prototype, not an emitted file/material.
use super::events::*;
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
#[test]
fn exact_cap_event_tokens_share_graph_and_internal_composite_junctions() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
        let closed = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
        let graph = closed
            .literal_vertices(CurveOptions::default(), None)
            .unwrap();
        let registry = graph.cap_events(DomainOptions::default(), None).unwrap();
        let events = registry.events();
        let segments = registry.segments();
        let mut w = Work::new(16_000_000, None).unwrap();
        w.charge(registry.work_used()).unwrap();
        let mut m = Memory::new(registry.retained_bytes(), 64 * 1024 * 1024).unwrap();
        if !left && !nb {
            use std::cell::Cell;
            let calls = Cell::new(0_usize);
            let count = || {
                calls.set(calls.get() + 1);
                false
            };
            graph
                .cap_events(DomainOptions::default(), Some(&count))
                .unwrap();
            for checkpoint in 1..=calls.get().min(4096) {
                let current = Cell::new(0_usize);
                let cancel = || {
                    current.set(current.get() + 1);
                    current.get() == checkpoint
                };
                assert!(
                    matches!(
                        graph.cap_events(DomainOptions::default(), Some(&cancel)),
                        Err(Error::Cancelled)
                    ),
                    "checkpoint={checkpoint}"
                );
            }
        }
        assert!(registry.coordinates(usize::MAX).is_none());
        let work_used = registry.work_used();
        assert!(
            graph
                .cap_events(
                    DomainOptions {
                        maximum_work: work_used,
                        ..Default::default()
                    },
                    None
                )
                .is_ok()
        );
        assert!(matches!(
            graph.cap_events(
                DomainOptions {
                    maximum_work: work_used - 1,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            graph.cap_events(
                DomainOptions {
                    maximum_bytes: graph.retained_bytes(),
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            graph.cap_events(DomainOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        // Independent source-token reconstruction and actual-field interval
        // oracle for every segment end, including all periodic seam carries.
        let mut oracle = Work::new(16_000_000, None).unwrap();
        let mut om = Memory::new(0, 64 * 1024 * 1024).unwrap();
        om.add(SCRATCH).unwrap();
        let mut ranges = Vec::with_capacity(segments.len());
        let mut direct_max = 0_f64;
        for segment in segments {
            assert_eq!(
                closed.edges()[segment.original_edge()].curve,
                ClosedCurve::Thread(segment.domain_edge())
            );
            let cell = field.cells[segment.source_cell()];
            let row = &field.meridians[cell.meridian];
            let (p0, starts, length) = if row.component == Component::Metric {
                let s = source.metric().spec();
                (s.pitch_mm, s.starts, s.length_mm)
            } else {
                let s = source.wood().spec();
                (s.pitch_mm, s.starts, s.length_mm)
            };
            let pitch = binary(p0, &mut oracle).unwrap();
            let lead = pitch
                .mul(&N::integer(i64::from(starts)), &mut oracle)
                .unwrap();
            let station = if segment.end() == End::Start {
                N::integer(0)
            } else {
                binary(length, &mut oracle).unwrap()
            };
            let shift = station
                .sub(
                    &pitch
                        .mul(
                            &N::integer(i64::from(cell.source_key.total_period)),
                            &mut oracle,
                        )
                        .unwrap(),
                    &mut oracle,
                )
                .unwrap();
            let a = binary(field.angular.knots()[segment.angular_piece()], &mut oracle).unwrap();
            let b = binary(
                field.angular.knots()[segment.angular_piece() + 1],
                &mut oracle,
            )
            .unwrap();
            let mut qrange = std::array::from_fn::<_, 2, _>(|_| N::integer(0));
            for (qout, &id) in qrange.iter_mut().zip(&segment.event_ids()) {
                let event = &events[id];
                let mut v = event.token.turn.number(&mut oracle).unwrap();
                if v.cmp(&a, &mut oracle).unwrap() == Ordering::Less {
                    v = v.add(&N::integer(1), &mut oracle).unwrap();
                }
                assert!(
                    v.cmp(&a, &mut oracle).unwrap() != Ordering::Less
                        && v.cmp(&b, &mut oracle).unwrap() != Ordering::Greater
                );
                let q = shift
                    .sub(&lead.mul(&v, &mut oracle).unwrap(), &mut oracle)
                    .unwrap();
                let before = om.used;
                let rebuilt = canonical(
                    (row.component, cell.source_key.total_period, starts),
                    q.clone(),
                    v.clone(),
                    &pitch,
                    &mut om,
                    &mut oracle,
                )
                .unwrap();
                assert!(equal(&rebuilt, &event.token, &mut oracle).unwrap());
                let qp = EncodedScalar::retain(&q, &mut om, &mut oracle).unwrap();
                let vp = EncodedScalar::retain(&v, &mut om, &mut oracle).unwrap();
                let (xy, z, _) =
                    footpoint_at(&field, segment.source_cell(), &qp, &vp, &mut oracle).unwrap();
                let direct = bind(coordinates(event, &graph), xy, &z, &mut oracle).unwrap();
                assert!(
                    direct <= event.field_incidence_error_mm(),
                    "direct={direct} receipt={}",
                    event.field_incidence_error_mm()
                );
                direct_max = direct_max.max(direct);
                *qout = q;
                drop((rebuilt, qp, vp));
                om.used = before;
            }
            assert!(qrange[0].cmp(&qrange[1], &mut oracle).unwrap() == Ordering::Less);
            ranges.push((segment.domain_edge(), qrange));
        }
        // Exact ordered phase cover per ORIGINAL graph edge, independently of
        // source-event ID sharing and the two-uses ring inventory below.
        for segment in segments {
            let mut own = ranges.iter().filter(|r| r.0 == segment.domain_edge()).rev();
            let mut last = own.next().unwrap().1[1].clone();
            let first = ranges
                .iter()
                .rfind(|r| r.0 == segment.domain_edge())
                .unwrap()
                .1[0]
                .clone();
            for (_, range) in own {
                assert_eq!(last.0, range[0].0);
                last = range[1].clone();
            }
            let before = om.used;
            let trace = field
                .cap_trace_with_ledgers(segment.source_cell(), segment.end(), &mut om, &mut oracle)
                .unwrap()
                .unwrap();
            assert_eq!(
                first.0,
                trace.endpoints[0].phase.number(&mut oracle).unwrap().0
            );
            assert_eq!(
                last.0,
                trace.endpoints[1].phase.number(&mut oracle).unwrap().0
            );
            drop(trace);
            om.used = before;
        }
        let mut uses = vec![0usize; events.len()];
        for segment in segments {
            assert!(segment.angular_piece() < field.angular.coefficients().len());
            assert!(segment.domain_edge() < domains.edges().len());
            assert_ne!(segment.event_ids()[0], segment.event_ids()[1]);
            for &e in &segment.event_ids() {
                uses[e] += 1;
                assert!(
                    events[e].source_incidence_error_mm() < 1e-7
                        && events[e].field_incidence_error_mm() < 1e-7
                );
                assert_eq!(
                    registry.coordinates(e).unwrap(),
                    coordinates(&events[e], &graph)
                );
            }
        }
        assert!(
            uses.iter().all(|&n| n == 2),
            "one source event shared by exactly two cap segments"
        );
        assert_eq!(segments.len(), if nb { 48 } else { 44 });
        assert_eq!(
            events.len(),
            segments.len(),
            "four complete one-winding rings"
        );
        // Actual rounded-root carriers now use the SAME chosen graph/event
        // points. Their complete true-event/source C0/C1 is independently
        // earned; IDs alone never grant endpoint or tensor incidence.
        let mut root_carriers = 0;
        let mut wrong_order_refused = false;
        let mut worst_source = [0_f64; 3];
        for segment in segments {
            let row = &field.meridians[field.cells[segment.source_cell()].meridian];
            if row.degree() != 2 {
                let targets = segment.event_ids().map(|e| coordinates(&events[e], &graph));
                if !wrong_order_refused {
                    let mut negative_work = Work::new(16_000_000, None).unwrap();
                    let mut negative_memory =
                        Memory::new(graph.retained_bytes(), 64 * 1024 * 1024).unwrap();
                    negative_memory.add(SCRATCH).unwrap();
                    assert!(matches!(
                        super::super::affine_probe::bind_vertices(
                            &field,
                            segment.source_cell(),
                            segment.end(),
                            segment.angular_piece(),
                            [targets[1], targets[0]],
                            &mut negative_memory,
                            &mut negative_work
                        ),
                        Err(Error::Unresolved)
                    ));
                    wrong_order_refused = true;
                }
                let proof = super::super::affine_probe::bind_vertices(
                    &field,
                    segment.source_cell(),
                    segment.end(),
                    segment.angular_piece(),
                    targets,
                    &mut m,
                    &mut w,
                )
                .unwrap();
                for (out, value) in worst_source.iter_mut().zip(proof) {
                    *out = out.max(value);
                }
                continue;
            }
            let curve = field
                .cap_curve_piece(
                    segment.source_cell(),
                    segment.end(),
                    segment.angular_piece(),
                    DomainOptions {
                        maximum_work: w.remaining(),
                        ..DomainOptions::default()
                    },
                    None,
                )
                .unwrap()
                .unwrap();
            let targets = segment.event_ids().map(|e| coordinates(&events[e], &graph));
            let proof = super::super::trim_probe::run_to_vertices(
                &curve,
                targets,
                field.work_used() + w.remaining(),
            )
            .unwrap();
            w.charge(proof.work - field.work_used()).unwrap();
            root_carriers += 1;
            worst_source[0] = worst_source[0].max(proof.source_position);
            worst_source[1] = worst_source[1].max(proof.source_derivative);
            worst_source[2] = worst_source[2].max(proof.original_source_derivative);
        }
        assert!(wrong_order_refused);
        assert_eq!(root_carriers, if nb { 5 } else { 4 });
        eprintln!(
            "shared true-event root carriers={root_carriers} source_C0_D1finite_D1original={worst_source:?} inherited_event_plus_root_work={}",
            w.used
        );
        let graph_count = events
            .iter()
            .filter(|e| matches!(e.point, Point::Graph(_)))
            .count();
        eprintln!(
            "cap event owners left={left} nb={nb} segments={} unique={} graph_reused={graph_count} actual_registry_work={} bytes={} direct_every_peer={direct_max}",
            segments.len(),
            events.len(),
            w.used,
            m.used - SCRATCH
        );
    }
}
