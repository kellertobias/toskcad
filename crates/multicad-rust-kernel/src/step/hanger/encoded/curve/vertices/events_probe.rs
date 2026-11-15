//! Actual finite cap-event owner prototype, not an emitted file/material.
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
struct Token {
    component: Component,
    period: i32,
    phase: EncodedScalar,
    turn: EncodedScalar,
}
enum Point {
    Graph(usize),
    Owned([EncodedScalar; 3]),
}
struct Event {
    token: Token,
    point: Point,
    source_error: f64,
    field_error: f64,
}
struct Segment {
    edge: usize,
    source_cell: usize,
    end: End,
    piece: usize,
    events: [usize; 2],
}
fn payload(x: &EncodedScalar) -> usize {
    x.n.len() + x.d.len()
}
fn equal(a: &Token, b: &Token, w: &mut Work<'_>) -> Result<bool, Error> {
    equal_parts(a, b.component, b.period, &b.phase, &b.turn, w)
}
fn equal_parts(
    a: &Token,
    component: Component,
    period: i32,
    phase: &EncodedScalar,
    turn: &EncodedScalar,
    w: &mut Work<'_>,
) -> Result<bool, Error> {
    w.charge(1 + payload(&a.phase) + payload(&a.turn) + payload(phase) + payload(turn))?;
    Ok(a.component == component && a.period == period && a.phase == *phase && a.turn == *turn)
}

fn canonical(
    domain: (Component, i32, u32),
    mut phase: N,
    mut turn: N,
    pitch: &N,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Token, Error> {
    let (component, mut period, starts) = domain;
    if phase.cmp(pitch, w)? == Ordering::Equal {
        phase = N::integer(0);
        period += 1;
    }
    if turn.cmp(&N::integer(1), w)? == Ordering::Equal {
        turn = N::integer(0);
        period += starts as i32;
    }
    Ok(Token {
        component,
        period,
        phase: EncodedScalar::retain(&phase, m, w)?,
        turn: EncodedScalar::retain(&turn, m, w)?,
    })
}
fn coordinates<'a>(
    event: &'a Event,
    graph: &'a LiteralVertices<'_, '_, '_, '_>,
) -> &'a [EncodedScalar; 3] {
    match &event.point {
        Point::Graph(v) => graph.vertices()[*v].coordinates(),
        Point::Owned(p) => p,
    }
}
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
        let mut w = Work::new(16_000_000, None).unwrap();
        w.charge(graph.work_used()).unwrap();
        let mut m = Memory::new(graph.retained_bytes(), 64 * 1024 * 1024).unwrap();
        m.add(SCRATCH + size_of::<CapTrace>()).unwrap();
        let cap_edges = closed
            .faces()
            .iter()
            .filter(|f| matches!(f.kind, ClosedFaceKind::Plane(_)))
            .flat_map(|f| f.loops().iter().flatten())
            .filter(|c| matches!(closed.edges()[c.edge].curve, ClosedCurve::Thread(_)))
            .count();
        let mut events: Vec<Event> = m
            .allocate(cap_edges * (field.angular.knots().len() + 1), &mut w)
            .unwrap();
        let mut segments: Vec<Segment> = m
            .allocate(cap_edges * field.angular.coefficients().len(), &mut w)
            .unwrap();
        for face in closed.faces() {
            let end = match face.kind {
                ClosedFaceKind::Plane(FaceKey::MetricTip | FaceKey::WoodShoulder) => End::Start,
                ClosedFaceKind::Plane(FaceKey::MetricShoulder | FaceKey::WoodTip) => End::End,
                _ => continue,
            };
            for coedge in face.loops().iter().flatten() {
                let ClosedCurve::Thread(edge) = closed.edges()[coedge.edge].curve else {
                    continue;
                };
                let (cell_index, _) = domains.edges()[edge].faces[0].unwrap();
                let index = domains.cells()[cell_index].source_cell;
                let cell = field.cells[index];
                let row = &field.meridians[cell.meridian];
                let (p0, starts, length) = if row.component == Component::Metric {
                    let s = source.metric().spec();
                    (s.pitch_mm, s.starts, s.length_mm)
                } else {
                    let s = source.wood().spec();
                    (s.pitch_mm, s.starts, s.length_mm)
                };
                let pitch = binary(p0, &mut w).unwrap();
                let lead = pitch.mul(&N::integer(i64::from(starts)), &mut w).unwrap();
                let station = if end == End::Start {
                    N::integer(0)
                } else {
                    binary(length, &mut w).unwrap()
                };
                let shift = station
                    .sub(
                        &pitch
                            .mul(&N::integer(i64::from(cell.source_key.total_period)), &mut w)
                            .unwrap(),
                        &mut w,
                    )
                    .unwrap();
                let before = m.used;
                let trace = field
                    .cap_trace_with_ledgers(index, end, &mut m, &mut w)
                    .unwrap()
                    .unwrap();
                let trace_bytes = m.used - before;
                for piece in 0..field.angular.coefficients().len() {
                    let a = binary(field.angular.knots()[piece], &mut w).unwrap();
                    let b = binary(field.angular.knots()[piece + 1], &mut w).unwrap();
                    let lower = shift.sub(&lead.mul(&b, &mut w).unwrap(), &mut w).unwrap();
                    let upper = shift.sub(&lead.mul(&a, &mut w).unwrap(), &mut w).unwrap();
                    let first = trace.endpoints[0].phase.number(&mut w).unwrap();
                    let last = trace.endpoints[1].phase.number(&mut w).unwrap();
                    let lo = if first.cmp(&lower, &mut w).unwrap() == Ordering::Greater {
                        first
                    } else {
                        lower
                    };
                    let hi = if last.cmp(&upper, &mut w).unwrap() == Ordering::Less {
                        last
                    } else {
                        upper
                    };
                    if lo.cmp(&hi, &mut w).unwrap() != Ordering::Less {
                        continue;
                    };
                    let mut ids = [0; 2];
                    for (id, q) in ids.iter_mut().zip([lo, hi]) {
                        let v = shift.sub(&q, &mut w).unwrap().div(&lead, &mut w).unwrap();
                        let token = canonical(
                            (row.component, cell.source_key.total_period, starts),
                            q.clone(),
                            v.clone(),
                            &pitch,
                            &mut m,
                            &mut w,
                        )
                        .unwrap();
                        let mut found = None;
                        for (i, e) in events.iter().enumerate() {
                            if equal(&e.token, &token, &mut w).unwrap() {
                                found = Some(i);
                                break;
                            }
                        }
                        if let Some(i) = found {
                            events[i].field_error = events[i].field_error.max(
                                (p(events[i].source_error) + p(row.bounds.position_error_mm)).hi,
                            );
                            m.used -= payload(&token.phase) + payload(&token.turn);
                            *id = i;
                            continue;
                        }
                        let mut graph_index = None;
                        for (i, v) in domains.vertices().iter().enumerate() {
                            if equal_parts(&token, v.component, v.period, &v.phase, &v.turn, &mut w)
                                .unwrap()
                            {
                                graph_index = Some(i);
                                break;
                            }
                        }
                        let (point, se, fe) = if let Some(i) = graph_index {
                            let point = &graph.vertices()[i];
                            (
                                Point::Graph(i),
                                point.source_error,
                                (p(point.source_error) + p(row.bounds.position_error_mm)).hi,
                            )
                        } else {
                            let temp = m.used;
                            let qp = EncodedScalar::retain(&q, &mut m, &mut w).unwrap();
                            let vp = EncodedScalar::retain(&v, &mut m, &mut w).unwrap();
                            let (xy, z, er) =
                                footpoint_at(&field, index, &qp, &vp, &mut w).unwrap();
                            m.used = temp;
                            drop((qp, vp));
                            let chosen = [
                                proposal(xy[0], &mut w).unwrap(),
                                proposal(xy[1], &mut w).unwrap(),
                                z.clone(),
                            ];
                            let point = [
                                EncodedScalar::retain(&chosen[0], &mut m, &mut w).unwrap(),
                                EncodedScalar::retain(&chosen[1], &mut m, &mut w).unwrap(),
                                EncodedScalar::retain(&chosen[2], &mut m, &mut w).unwrap(),
                            ];
                            let fe = bind(&point, xy, &z, &mut w).unwrap();
                            (Point::Owned(point), (p(fe) + p(er)).hi, fe)
                        };
                        *id = events.len();
                        events.push(Event {
                            token,
                            point,
                            source_error: se,
                            field_error: fe,
                        });
                    }
                    segments.push(Segment {
                        edge,
                        source_cell: index,
                        end,
                        piece,
                        events: ids,
                    });
                }
                drop(trace);
                m.used -= trace_bytes;
            }
        }
        // Independent source-token reconstruction and actual-field interval
        // oracle for every segment end, including all periodic seam carries.
        let mut oracle = Work::new(16_000_000, None).unwrap();
        let mut om = Memory::new(0, 64 * 1024 * 1024).unwrap();
        om.add(SCRATCH).unwrap();
        let mut ranges = Vec::with_capacity(segments.len());
        let mut direct_max = 0_f64;
        for segment in &segments {
            let cell = field.cells[segment.source_cell];
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
            let station = if segment.end == End::Start {
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
            let a = binary(field.angular.knots()[segment.piece], &mut oracle).unwrap();
            let b = binary(field.angular.knots()[segment.piece + 1], &mut oracle).unwrap();
            let mut qrange = std::array::from_fn::<_, 2, _>(|_| N::integer(0));
            for (qout, &id) in qrange.iter_mut().zip(&segment.events) {
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
                    footpoint_at(&field, segment.source_cell, &qp, &vp, &mut oracle).unwrap();
                let direct = bind(coordinates(event, &graph), xy, &z, &mut oracle).unwrap();
                assert!(
                    direct <= event.field_error,
                    "direct={direct} receipt={}",
                    event.field_error
                );
                direct_max = direct_max.max(direct);
                *qout = q;
                drop((rebuilt, qp, vp));
                om.used = before;
            }
            assert!(qrange[0].cmp(&qrange[1], &mut oracle).unwrap() == Ordering::Less);
            ranges.push((segment.edge, qrange));
        }
        // Exact ordered phase cover per ORIGINAL graph edge, independently of
        // source-event ID sharing and the two-uses ring inventory below.
        for segment in &segments {
            let mut own = ranges.iter().filter(|r| r.0 == segment.edge).rev();
            let mut last = own.next().unwrap().1[1].clone();
            let first = ranges.iter().rfind(|r| r.0 == segment.edge).unwrap().1[0].clone();
            for (_, range) in own {
                assert_eq!(last.0, range[0].0);
                last = range[1].clone();
            }
            let before = om.used;
            let trace = field
                .cap_trace_with_ledgers(segment.source_cell, segment.end, &mut om, &mut oracle)
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
        for segment in &segments {
            assert!(segment.piece < field.angular.coefficients().len());
            assert!(segment.edge < domains.edges().len());
            assert_ne!(segment.events[0], segment.events[1]);
            for &e in &segment.events {
                uses[e] += 1;
                assert!(events[e].source_error < 1e-7 && events[e].field_error < 1e-7);
                let _ = coordinates(&events[e], &graph);
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
        for segment in &segments {
            let row = &field.meridians[field.cells[segment.source_cell].meridian];
            if row.degree() != 2 {
                let targets = segment.events.map(|e| coordinates(&events[e], &graph));
                if !wrong_order_refused {
                    let mut negative_work = Work::new(16_000_000, None).unwrap();
                    let mut negative_memory =
                        Memory::new(graph.retained_bytes(), 64 * 1024 * 1024).unwrap();
                    negative_memory.add(SCRATCH).unwrap();
                    assert!(matches!(
                        super::super::affine_probe::bind_vertices(
                            &field,
                            segment.source_cell,
                            segment.end,
                            segment.piece,
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
                    segment.source_cell,
                    segment.end,
                    segment.piece,
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
                    segment.source_cell,
                    segment.end,
                    segment.piece,
                    DomainOptions {
                        maximum_work: w.remaining(),
                        ..DomainOptions::default()
                    },
                    None,
                )
                .unwrap()
                .unwrap();
            let targets = segment.events.map(|e| coordinates(&events[e], &graph));
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
