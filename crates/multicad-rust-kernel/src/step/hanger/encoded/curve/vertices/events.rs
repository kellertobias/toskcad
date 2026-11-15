//! Source-owned cap event registry with exact canonical tokens and literal
//! graph/junction points. No rounded geometric key or proximity merge is used.
use super::*;
pub(super) struct Token {
    pub(super) component: Component,
    pub(super) period: i32,
    pub(super) phase: EncodedScalar,
    pub(super) turn: EncodedScalar,
}
pub(super) enum Point {
    Graph(usize),
    Owned([EncodedScalar; 3]),
}
pub struct Event {
    pub(super) token: Token,
    pub(super) point: Point,
    pub(super) source_error: f64,
    pub(super) field_error: f64,
}
pub struct Segment {
    pub(super) edge: usize,
    pub(super) closed_edge: usize,
    pub(super) source_cell: usize,
    pub(super) end: End,
    pub(super) piece: usize,
    pub(super) events: [usize; 2],
}
pub(super) fn payload(x: &EncodedScalar) -> usize {
    x.n.len() + x.d.len()
}
pub(super) fn equal(a: &Token, b: &Token, w: &mut Work<'_>) -> Result<bool, Error> {
    equal_parts(a, b.component, b.period, &b.phase, &b.turn, w)
}
pub(super) fn equal_parts(
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

pub(super) fn canonical(
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
pub(super) fn coordinates<'a>(
    event: &'a Event,
    graph: &'a LiteralVertices<'_, '_, '_, '_>,
) -> &'a [EncodedScalar; 3] {
    match &event.point {
        Point::Graph(v) => graph.vertices()[*v].coordinates(),
        Point::Owned(p) => p,
    }
}

pub struct CapEvents<'g, 'c, 'd, 'f, 's> {
    graph: &'g LiteralVertices<'c, 'd, 'f, 's>,
    events: Vec<Event>,
    segments: Vec<Segment>,
    work: usize,
    bytes: usize,
}
impl CapEvents<'_, '_, '_, '_, '_> {
    pub fn events(&self) -> &[Event] {
        &self.events
    }
    pub fn segments(&self) -> &[Segment] {
        &self.segments
    }
    pub fn coordinates(&self, index: usize) -> Option<&[EncodedScalar; 3]> {
        self.events
            .get(index)
            .map(|event| coordinates(event, self.graph))
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
impl Event {
    pub fn source_incidence_error_mm(&self) -> f64 {
        self.source_error
    }
    pub fn field_incidence_error_mm(&self) -> f64 {
        self.field_error
    }
}
impl Segment {
    pub fn original_edge(&self) -> usize {
        self.closed_edge
    }
    pub fn domain_edge(&self) -> usize {
        self.edge
    }
    pub fn source_cell(&self) -> usize {
        self.source_cell
    }
    pub fn end(&self) -> End {
        self.end
    }
    pub fn angular_piece(&self) -> usize {
        self.piece
    }
    pub fn event_ids(&self) -> [usize; 2] {
        self.events
    }
}
impl<'c, 'd, 'f, 's> LiteralVertices<'c, 'd, 'f, 's> {
    pub fn cap_events<'g>(
        &'g self,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<CapEvents<'g, 'c, 'd, 'f, 's>, Error> {
        let graph = self;
        let closed = graph.source();
        let domains = closed.domains();
        let field = domains.field();
        let source = field.source.source_boundary().source();
        let (mut w, mut m) =
            crate::step::hanger::domain::ledgers(options, graph.retained_bytes(), interrupt)?;
        w.charge(graph.work_used())?;
        m.add(SCRATCH + size_of::<CapTrace>() + size_of::<CapEvents<'_, '_, '_, '_, '_>>())?;
        let cap_edges = closed
            .faces()
            .iter()
            .filter(|f| matches!(f.kind, ClosedFaceKind::Plane(_)))
            .flat_map(|f| f.loops().iter().flatten())
            .filter(|c| matches!(closed.edges()[c.edge].curve, ClosedCurve::Thread(_)))
            .count();
        let mut events: Vec<Event> =
            m.allocate(cap_edges * (field.angular.knots().len() + 1), &mut w)?;
        let mut segments: Vec<Segment> =
            m.allocate(cap_edges * field.angular.coefficients().len(), &mut w)?;
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
                let (cell_index, _) = domains.edges()[edge].faces[0].ok_or(Error::Unresolved)?;
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
                let pitch = binary(p0, &mut w)?;
                let lead = pitch.mul(&N::integer(i64::from(starts)), &mut w)?;
                let station = if end == End::Start {
                    N::integer(0)
                } else {
                    binary(length, &mut w)?
                };
                let shift = station.sub(
                    &pitch.mul(&N::integer(i64::from(cell.source_key.total_period)), &mut w)?,
                    &mut w,
                )?;
                let before = m.used;
                let trace = field
                    .cap_trace_with_ledgers(index, end, &mut m, &mut w)?
                    .ok_or(Error::Unresolved)?;
                let trace_bytes = m.used - before;
                for piece in 0..field.angular.coefficients().len() {
                    let a = binary(field.angular.knots()[piece], &mut w)?;
                    let b = binary(field.angular.knots()[piece + 1], &mut w)?;
                    let lower = shift.sub(&lead.mul(&b, &mut w)?, &mut w)?;
                    let upper = shift.sub(&lead.mul(&a, &mut w)?, &mut w)?;
                    let first = trace.endpoints[0].phase.number(&mut w)?;
                    let last = trace.endpoints[1].phase.number(&mut w)?;
                    let lo = if first.cmp(&lower, &mut w)? == Ordering::Greater {
                        first
                    } else {
                        lower
                    };
                    let hi = if last.cmp(&upper, &mut w)? == Ordering::Less {
                        last
                    } else {
                        upper
                    };
                    if lo.cmp(&hi, &mut w)? != Ordering::Less {
                        continue;
                    };
                    let mut ids = [0; 2];
                    for (id, q) in ids.iter_mut().zip([lo, hi]) {
                        let v = shift.sub(&q, &mut w)?.div(&lead, &mut w)?;
                        let token = canonical(
                            (row.component, cell.source_key.total_period, starts),
                            q.clone(),
                            v.clone(),
                            &pitch,
                            &mut m,
                            &mut w,
                        )?;
                        let mut found = None;
                        for (i, e) in events.iter().enumerate() {
                            if equal(&e.token, &token, &mut w)? {
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
                            if equal_parts(
                                &token,
                                v.component,
                                v.period,
                                &v.phase,
                                &v.turn,
                                &mut w,
                            )? {
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
                            let qp = EncodedScalar::retain(&q, &mut m, &mut w)?;
                            let vp = EncodedScalar::retain(&v, &mut m, &mut w)?;
                            let (xy, z, er) = footpoint_at(field, index, &qp, &vp, &mut w)?;
                            m.used = temp;
                            drop((qp, vp));
                            let chosen = [
                                proposal(xy[0], &mut w)?,
                                proposal(xy[1], &mut w)?,
                                z.clone(),
                            ];
                            let point = [
                                EncodedScalar::retain(&chosen[0], &mut m, &mut w)?,
                                EncodedScalar::retain(&chosen[1], &mut m, &mut w)?,
                                EncodedScalar::retain(&chosen[2], &mut m, &mut w)?,
                            ];
                            let fe = bind(&point, xy, &z, &mut w)?;
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
                        closed_edge: coedge.edge,
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

        w.charge(1)?;
        Ok(CapEvents {
            graph,
            events,
            segments,
            work: w.used,
            bytes: m.used,
        })
    }
}
