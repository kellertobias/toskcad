use super::*;
pub(super) struct Cell {
    pub planes: Vec<[Q; 3]>,
    pub owners: Vec<usize>,
}
fn point(p: [f64; 2]) -> R<[Q; 2]> {
    Ok([q(p[0])?, q(p[1])?])
}
fn plane(a: &[Q; 2], b: &[Q; 2]) -> R<[Q; 3]> {
    let dx = bounded(&b[0] - &a[0])?;
    let dy = bounded(&b[1] - &a[1])?;
    Ok([-&dy, dx.clone(), bounded(&dy * &a[0] - &dx * &a[1])?])
}
fn evaluate(p: &[Q; 3], x: &[Q; 2]) -> R<Q> {
    bounded(&p[0] * &x[0] + &p[1] * &x[1] + &p[2])
}
/// Exact axis-aligned enclosure of the full retained circular support AND
/// literal finite endpoints. It never replaces an arc by its chord for bounds.
pub(super) fn edge_bounds(edge: Edge, work: &mut Work<'_>) -> R<[[Q; 2]; 2]> {
    work.charge(16)?;
    match edge {
        Edge::Line { a, b } => Ok(std::array::from_fn(|i| {
            [q(a[i].min(b[i])).unwrap(), q(a[i].max(b[i])).unwrap()]
        })),
        Edge::Arc { center, a, b, .. } => {
            let radius = q(edge.radius_mm().ok_or(Error::Unsupported)?)?;
            // Exact finite endpoint/support incidence is required here. The
            // source constructor's joining tolerance alone is insufficient for
            // the closed-boundary winding proof. Original Line8 corner arcs
            // satisfy this identity exactly in their unchanged local section.
            for endpoint in [a, b] {
                work.charge(12)?;
                let delta = [
                    bounded(q(endpoint[0])? - q(center[0])?)?,
                    bounded(q(endpoint[1])? - q(center[1])?)?,
                ];
                if bounded(&delta[0] * &delta[0] + &delta[1] * &delta[1])? != &radius * &radius {
                    return Err(Error::Unsupported);
                }
            }
            let mut result = [[Q::zero(), Q::zero()], [Q::zero(), Q::zero()]];
            for i in 0..2 {
                let c = q(center[i])?;
                result[i] = [
                    bounded(&c - &radius)?.min(q(a[i])?).min(q(b[i])?),
                    bounded(&c + &radius)?.max(q(a[i])?).max(q(b[i])?),
                ];
            }
            Ok(result)
        }
        Edge::Circle {
            center, radius_mm, ..
        } => {
            let r = q(radius_mm)?;
            Ok([
                [bounded(q(center[0])? - &r)?, bounded(q(center[0])? + &r)?],
                [bounded(q(center[1])? - &r)?, bounded(q(center[1])? + &r)?],
            ])
        }
    }
}
impl Cell {
    pub fn from_boundary(outer: &[Edge], start: usize, work: &mut Work<'_>) -> R<Option<Self>> {
        if !matches!(outer[start], Edge::Line { .. }) {
            return Ok(None);
        }
        // Only maximal clockwise line runs are candidates. Sub-runs would
        // repeat whole-ring qualification without improving the containing
        // convex cell of the actual groove. Other sections may conservatively
        // refuse when their maximal run does not earn a convex exterior cell.
        let previous = outer[(start + outer.len() - 1) % outer.len()];
        if let Edge::Line { a, b } = previous {
            work.charge(24)?;
            if b != outer[start].start() {
                return Ok(None);
            }
            let turn = evaluate(&plane(&point(a)?, &point(b)?)?, &point(outer[start].end())?)?;
            if turn <= Q::zero() {
                return Ok(None);
            }
        }
        let mut vertices = vec![point(outer[start].start())?];
        let mut owners = Vec::new();
        for n in 0..outer.len() {
            work.charge(24)?;
            let index = (start + n) % outer.len();
            let Edge::Line { a, b } = outer[index] else {
                break;
            };
            let a = point(a)?;
            let b = point(b)?;
            if vertices.last() != Some(&a) {
                return Ok(None);
            }
            if vertices.len() >= 2 {
                let p = plane(&vertices[vertices.len() - 2], &a)?;
                if evaluate(&p, &b)? > Q::zero() {
                    break;
                }
            }
            if a == b {
                return Ok(None);
            }
            owners.push(index);
            vertices.push(b);
        }
        if owners.len() < 3 || vertices.first() == vertices.last() {
            return Ok(None);
        }
        Self::qualify(outer, vertices, owners, work)
    }
    pub fn mouth(
        outer: &[Edge],
        first: usize,
        second: usize,
        ranges: &[[Q; 2]; 2],
        work: &mut Work<'_>,
    ) -> R<Option<Self>> {
        work.charge(64)?;
        let (Edge::Line { a, b }, Edge::Line { a: c, b: d }) = (outer[first], outer[second]) else {
            return Ok(None);
        };
        let a = point(a)?;
        let b = point(b)?;
        let c = point(c)?;
        let d = point(d)?;
        for along in 0..2 {
            let across = 1 - along;
            if a[across] != b[across]
                || c[across] != d[across]
                || a[along] != d[along]
                || b[along] != c[along]
                || a[along] == b[along]
                || a[across] == c[across]
            {
                continue;
            }
            let low = a[along]
                .clone()
                .min(b[along].clone())
                .min(ranges[along][0].clone());
            let high = a[along]
                .clone()
                .max(b[along].clone())
                .max(ranges[along][1].clone());
            let mut vertices = vec![a.clone(), b.clone(), c.clone(), d.clone()];
            for (i, v) in vertices.iter_mut().enumerate() {
                let lower = if i == 0 || i == 3 {
                    a[along] < b[along]
                } else {
                    b[along] < a[along]
                };
                v[along] = if lower { low.clone() } else { high.clone() };
            }
            return Self::qualify(outer, vertices, vec![first, second], work);
        }
        Ok(None)
    }
    fn qualify(
        outer: &[Edge],
        vertices: Vec<[Q; 2]>,
        owners: Vec<usize>,
        work: &mut Work<'_>,
    ) -> R<Option<Self>> {
        let mut planes = Vec::with_capacity(vertices.len());
        let mut area = Q::zero();
        for i in 0..vertices.len() {
            work.charge(24)?;
            let a = &vertices[i];
            let b = &vertices[(i + 1) % vertices.len()];
            if a == b {
                return Ok(None);
            }
            area = bounded(area + &a[0] * &b[1] - &a[1] * &b[0])?;
            planes.push(plane(a, b)?);
        }
        if area >= Q::zero() {
            return Ok(None);
        }
        for p in &planes {
            for x in &vertices {
                work.charge(8)?;
                if evaluate(p, x)? > Q::zero() {
                    return Ok(None);
                }
            }
        }
        // Every foreign original owner is separated from the cell interior.
        // Arc disks (including finite endpoints) earn a conservative plane
        // separator. This also permits a winding-preserving chord homotopy.
        for edge in outer.iter().copied() {
            work.charge(owners.len() + 1)?;
            if let Edge::Line { a, b } = edge {
                // Exact strict half-plane clipping: boundary contact is allowed.
                // A line's AABB can straddle a corner even when the line does not.
                let a = point(a)?;
                let b = point(b)?;
                let mut lo = Q::zero();
                let mut hi = Q::from_integer(1.into());
                let mut separated = false;
                for p in &planes {
                    work.charge(32)?;
                    let v = evaluate(p, &a)?;
                    let d = bounded(evaluate(p, &b)? - &v)?;
                    if d.is_zero() {
                        if v >= Q::zero() {
                            separated = true;
                            break;
                        }
                    } else {
                        let t = bounded(-v / &d)?;
                        if d > Q::zero() {
                            hi = hi.min(t);
                        } else {
                            lo = lo.max(t);
                        }
                    }
                }
                if !separated && lo < hi {
                    return Ok(None);
                }
                continue;
            }
            let bounds = edge_bounds(edge, work)?;
            let mut separated = false;
            for p in &planes {
                work.charge(8)?;
                let corner =
                    std::array::from_fn(|i| bounds[i][usize::from(p[i] < Q::zero())].clone());
                if evaluate(p, &corner)? >= Q::zero() {
                    separated = true;
                    break;
                }
            }
            if !separated {
                return Ok(None);
            }
        }
        // The vertex mean lies strictly inside a nondegenerate convex cell.
        // All circular supports above are homotopic to their endpoint chords
        // outside this cell. Exact polygon winding therefore proves the actual
        // outer boundary has zero winding throughout the connected cell.
        let mut center = [Q::zero(), Q::zero()];
        for x in &vertices {
            work.charge(8)?;
            for i in 0..2 {
                center[i] = bounded(&center[i] + &x[i])?;
            }
        }
        work.charge(16)?;
        let count = Q::from_integer(vertices.len().into());
        for c in &mut center {
            *c = bounded(&*c / &count)?;
        }
        for p in &planes {
            work.charge(8)?;
            if evaluate(p, &center)? >= Q::zero() {
                return Ok(None);
            }
        }
        let mut winding = 0isize;
        for (i, edge) in outer.iter().copied().enumerate() {
            work.charge(24)?;
            if matches!(edge, Edge::Circle { .. }) {
                return Ok(None);
            }
            if edge.end() != outer[(i + 1) % outer.len()].start() {
                return Ok(None);
            }
            let a = point(edge.start())?;
            let b = point(edge.end())?;
            let side = evaluate(&plane(&a, &b)?, &center)?;
            if a[1] <= center[1] && center[1] < b[1] && side > Q::zero() {
                winding += 1;
            }
            if b[1] <= center[1] && center[1] < a[1] && side < Q::zero() {
                winding -= 1;
            }
        }
        if winding != 0 {
            return Ok(None);
        }
        Ok(Some(Self { planes, owners }))
    }
}
