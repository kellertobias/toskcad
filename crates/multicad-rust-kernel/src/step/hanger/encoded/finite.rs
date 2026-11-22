//! Exact finite phase/turn domains and source-token seam ownership.
//! Plane fills/plain wall and literal STEP encoding are NOT admitted here.
use super::*;
use std::cmp::Ordering;
struct Point {
    q: N,
    v: N,
}
pub struct FiniteCorner {
    pub vertex: usize,
    pub phase: EncodedScalar,
    pub turn: EncodedScalar,
}
pub struct FiniteCell {
    pub source_cell: usize,
    corners: Vec<FiniteCorner>,
}
impl FiniteCell {
    pub fn corners(&self) -> &[FiniteCorner] {
        &self.corners
    }
}
pub struct FiniteVertex {
    pub component: Component,
    pub period: i32,
    pub phase: EncodedScalar,
    pub turn: EncodedScalar,
}
pub struct FiniteEdge {
    pub vertices: [usize; 2],
    pub faces: [Option<(usize, bool)>; 2],
}
pub struct FiniteDomains<'f, 's> {
    field: &'f EncodedField<'s>,
    cells: Vec<FiniteCell>,
    vertices: Vec<FiniteVertex>,
    edges: Vec<FiniteEdge>,
    work: usize,
    bytes: usize,
}
fn value(p: &Point, shift: &N, lead: &N, w: &mut Work<'_>) -> Result<N, Error> {
    p.q.add(shift, w)?.add(&lead.mul(&p.v, w)?, w)
}
fn append(out: &mut Vec<Point>, point: Point, w: &mut Work<'_>) -> Result<(), Error> {
    if let Some(last) = out.last()
        && last.q.cmp(&point.q, w)? == Ordering::Equal
        && last.v.cmp(&point.v, w)? == Ordering::Equal
    {
        return Ok(());
    }
    out.push(point);
    Ok(())
}
fn clip(
    input: Vec<Point>,
    shift: &N,
    lead: &N,
    limit: &N,
    upper: bool,
    w: &mut Work<'_>,
) -> Result<Vec<Point>, Error> {
    let mut out = Vec::with_capacity(6);
    if input.is_empty() {
        return Ok(out);
    }
    for i in 0..input.len() {
        let a = &input[i];
        let b = &input[(i + 1) % input.len()];
        let mut da = value(a, shift, lead, w)?.sub(limit, w)?;
        let mut db = value(b, shift, lead, w)?.sub(limit, w)?;
        if upper {
            da = N::integer(0).sub(&da, w)?;
            db = N::integer(0).sub(&db, w)?;
        }
        let sa = da.cmp(&N::integer(0), w)?;
        let sb = db.cmp(&N::integer(0), w)?;
        if sa != Ordering::Less {
            append(
                &mut out,
                Point {
                    q: a.q.clone(),
                    v: a.v.clone(),
                },
                w,
            )?;
        }
        if (sa == Ordering::Less && sb == Ordering::Greater)
            || (sa == Ordering::Greater && sb == Ordering::Less)
        {
            let t = da.div(&da.sub(&db, w)?, w)?;
            append(
                &mut out,
                Point {
                    q: a.q.add(&b.q.sub(&a.q, w)?.mul(&t, w)?, w)?,
                    v: a.v.add(&b.v.sub(&a.v, w)?.mul(&t, w)?, w)?,
                },
                w,
            )?;
        }
    }
    if out.len() > 1 {
        let first = &out[0];
        let last = &out[out.len() - 1];
        if first.q.cmp(&last.q, w)? == Ordering::Equal
            && first.v.cmp(&last.v, w)? == Ordering::Equal
        {
            out.pop();
        }
    }
    Ok(out)
}
fn twice_area(points: &[Point], w: &mut Work<'_>) -> Result<N, Error> {
    let mut result = N::integer(0);
    for i in 0..points.len() {
        let a = &points[i];
        let b = &points[(i + 1) % points.len()];
        result = result.add(&a.q.mul(&b.v, w)?.sub(&b.q.mul(&a.v, w)?, w)?, w)?;
    }
    Ok(result)
}
fn payload(a: &EncodedScalar, b: &EncodedScalar) -> usize {
    a.n.len() + a.d.len() + b.n.len() + b.d.len()
}
fn vertex(
    vertices: &mut Vec<FiniteVertex>,
    mut point: Point,
    owner: (i32, Component, u32),
    pitch: &N,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<usize, Error> {
    let (mut period, component, starts) = owner;
    if point.q.cmp(pitch, w)? == Ordering::Equal {
        point.q = N::integer(0);
        period += 1;
    }
    if point.v.cmp(&N::integer(1), w)? == Ordering::Equal {
        point.v = N::integer(0);
        period += starts as i32;
    }
    let phase = EncodedScalar::retain(&point.q, m, w)?;
    let turn = EncodedScalar::retain(&point.v, m, w)?;
    for (i, old) in vertices.iter().enumerate() {
        w.charge(1)?;
        if old.component != component || old.period != period {
            continue;
        }
        w.charge(1 + payload(&phase, &turn) + payload(&old.phase, &old.turn))?;
        if old.phase == phase && old.turn == turn {
            m.used -= payload(&phase, &turn);
            return Ok(i);
        }
    }
    vertices.push(FiniteVertex {
        component,
        period,
        phase,
        turn,
    });
    Ok(vertices.len() - 1)
}
fn positive_polygon(polygon: &[Point], w: &mut Work<'_>) -> Result<Option<N>, Error> {
    let area = twice_area(polygon, w)?;
    match area.cmp(&N::integer(0), w)? {
        Ordering::Equal => Ok(None),
        Ordering::Less => Err(Error::Unresolved),
        Ordering::Greater => {
            if polygon.len() > 6 {
                Err(Error::Unresolved)
            } else {
                Ok(Some(area))
            }
        }
    }
}
fn add_edge(
    edges: &mut Vec<FiniteEdge>,
    vertices: [usize; 2],
    face: usize,
    w: &mut Work<'_>,
) -> Result<(), Error> {
    let [a, b] = vertices;
    if a == b {
        return Err(Error::Unresolved);
    }
    let pair = [a.min(b), a.max(b)];
    let sense = a < b;
    let mut found = None;
    for (i, edge) in edges.iter().enumerate() {
        w.charge(1)?;
        if edge.vertices == pair {
            found = Some(i);
            break;
        }
    }
    if let Some(i) = found {
        let edge = &mut edges[i];
        if edge.faces[1].is_some() || edge.faces[0].is_some_and(|(_, s)| s == sense) {
            return Err(Error::Unresolved);
        }
        edge.faces[1] = Some((face, sense));
    } else {
        edges.push(FiniteEdge {
            vertices: pair,
            faces: [Some((face, sense)), None],
        });
    }
    Ok(())
}
fn check_totals(field: &EncodedField<'_>, totals: &[N; 2], w: &mut Work<'_>) -> Result<(), Error> {
    let source = field.source.source_boundary().source();
    for (i, length) in [
        source.metric().full_length_mm(),
        source.wood().full_length_mm(),
    ]
    .into_iter()
    .enumerate()
    {
        if totals[i].cmp(&binary(length, w)?.mul(&N::integer(2), w)?, w)? != Ordering::Equal {
            return Err(Error::Unresolved);
        }
    }
    Ok(())
}
impl<'f, 's> FiniteDomains<'f, 's> {
    pub fn prepare(
        field: &'f EncodedField<'s>,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut w = super::super::domain::work(options, field.bytes, interrupt)?;
        w.charge(field.work_used())?;
        let mut memory = Memory::new(field.bytes + size_of::<Self>(), options.maximum_bytes)?;
        memory.add(SCRATCH)?;
        // <=6 corners after clipping a rectangle by two PARALLEL height
        // halfplanes. <=64 live checked exact values + normalization and
        // small polygon Vec/boxing temporaries fit this existing scratch.
        let mut cells = memory.allocate(field.cells.len(), &mut w)?;
        let mut vertices = memory.allocate(6 * field.cells.len(), &mut w)?;
        let mut edges: Vec<FiniteEdge> = memory.allocate(6 * field.cells.len(), &mut w)?;
        let mut totals = [N::integer(0), N::integer(0)];
        let source = field.source.source_boundary().source();
        for (index, cell) in field.cells.iter().enumerate() {
            let row = &field.meridians[cell.meridian];
            let (p, n, l) = if row.component == Component::Metric {
                let s = source.metric().spec();
                (s.pitch_mm, s.starts, s.length_mm)
            } else {
                let s = source.wood().spec();
                (s.pitch_mm, s.starts, s.length_mm)
            };
            let pitch = binary(p, &mut w)?;
            let lead = pitch.mul(&N::integer(i64::from(n)), &mut w)?;
            let length = binary(l, &mut w)?;
            let shift = pitch.mul(&N::integer(i64::from(cell.source_key.total_period)), &mut w)?;
            let a = row.phase[0]
                .number(&mut w)?
                .div(&row.weights[0].number(&mut w)?, &mut w)?;
            let b = row.phase[row.degree()]
                .number(&mut w)?
                .div(&row.weights[row.degree()].number(&mut w)?, &mut w)?;
            let polygon = vec![
                Point {
                    q: a.clone(),
                    v: N::integer(0),
                },
                Point {
                    q: b.clone(),
                    v: N::integer(0),
                },
                Point {
                    q: b,
                    v: N::integer(1),
                },
                Point {
                    q: a,
                    v: N::integer(1),
                },
            ];
            let polygon = clip(
                clip(polygon, &shift, &lead, &N::integer(0), false, &mut w)?,
                &shift,
                &lead,
                &length,
                true,
                &mut w,
            )?;
            let Some(area) = positive_polygon(&polygon, &mut w)? else {
                continue;
            };
            let component_index = usize::from(row.component == Component::Wood);
            totals[component_index] = totals[component_index].add(&area, &mut w)?;
            let face = cells.len();
            let mut corners = memory.allocate(polygon.len(), &mut w)?;
            for point in polygon {
                let phase = EncodedScalar::retain(&point.q, &mut memory, &mut w)?;
                let turn = EncodedScalar::retain(&point.v, &mut memory, &mut w)?;
                let v = vertex(
                    &mut vertices,
                    point,
                    (cell.source_key.total_period, row.component, n),
                    &pitch,
                    &mut memory,
                    &mut w,
                )?;
                corners.push(FiniteCorner {
                    vertex: v,
                    phase,
                    turn,
                });
            }
            for i in 0..corners.len() {
                let a = corners[i].vertex;
                let b = corners[(i + 1) % corners.len()].vertex;
                add_edge(&mut edges, [a, b], face, &mut w)?;
            }
            cells.push(FiniteCell {
                source_cell: index,
                corners,
            });
        }
        check_totals(field, &totals, &mut w)?;
        w.charge(1)?;
        memory.used -= SCRATCH;
        Ok(Self {
            field,
            cells,
            vertices,
            edges,
            work: w.used,
            bytes: memory.used,
        })
    }
    pub fn field(&self) -> &EncodedField<'s> {
        self.field
    }
    pub fn cells(&self) -> &[FiniteCell] {
        &self.cells
    }
    pub fn vertices(&self) -> &[FiniteVertex] {
        &self.vertices
    }
    pub fn edges(&self) -> &[FiniteEdge] {
        &self.edges
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
#[cfg(test)]
mod tests;
