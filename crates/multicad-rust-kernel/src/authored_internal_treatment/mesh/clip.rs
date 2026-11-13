//! Exact convex source-domain clipping. Constraint identities belong to the
//! authored parameter graph; nearby numerical coordinates are never merged.
use super::{exact::Q, *};
use std::cmp::Ordering;
#[derive(Clone)]
pub(super) struct Point {
    pub q: Q,
    pub v: Q,
    pub active: u8,
}
#[derive(Clone)]
pub(super) struct Line {
    pub a: Q,
    pub b: Q,
    pub c: Q,
    bounds: [I; 3],
}
impl Line {
    pub fn new(a: Q, b: Q, c: Q, w: &mut Work<'_>) -> Result<Self, Error> {
        let bounds = [a.enclosure(w)?.1, b.enclosure(w)?.1, c.enclosure(w)?.1];
        Ok(Self { a, b, c, bounds })
    }
    pub fn value(&self, p: &Point, w: &mut Work<'_>) -> Result<Ordering, Error> {
        w.step()?;
        let quick = self.bounds[0] * p.q.enclosure(w)?.1
            + self.bounds[1] * p.v.enclosure(w)?.1
            + self.bounds[2];
        if quick.lo > 0. {
            return Ok(Ordering::Greater);
        }
        if quick.hi < 0. {
            return Ok(Ordering::Less);
        }
        Ok(self
            .a
            .mul(&p.q, w)?
            .add(&self.b.mul(&p.v, w)?, w)?
            .add(&self.c, w)?
            .sign())
    }
    fn intersection(&self, other: &Self, w: &mut Work<'_>) -> Result<Point, Error> {
        let det = self.a.mul(&other.b, w)?.sub(&other.a.mul(&self.b, w)?, w)?;
        let q = self
            .b
            .mul(&other.c, w)?
            .sub(&other.b.mul(&self.c, w)?, w)?
            .div(&det, w)?;
        let v = self
            .c
            .mul(&other.a, w)?
            .sub(&other.c.mul(&self.a, w)?, w)?
            .div(&det, w)?;
        Ok(Point { q, v, active: 0 })
    }
}
/// Four rectangle support lines followed by at most two physical constraints.
/// Each produced point retains every exact incident constraint, including a
/// physical station landing on an existing grid corner.
pub(super) fn polygon(
    rect: [Q; 4],
    cuts: &[Line],
    w: &mut Work<'_>,
) -> Result<(Vec<Point>, Vec<Line>), Error> {
    if cuts.len() > 2 {
        return Err(Error::Limit);
    }
    let [ql, qh, vl, vh] = rect;
    if ql.cmp(&qh, w)? != Ordering::Less || vl.cmp(&vh, w)? != Ordering::Less {
        return Err(Error::Unresolved);
    }
    let zero = Q::integer(0, w)?;
    let one = Q::integer(1, w)?;
    let mut lines = vec![
        Line::new(one.clone(), zero.clone(), zero.sub(&ql, w)?, w)?,
        Line::new(one.clone(), zero.clone(), zero.sub(&qh, w)?, w)?,
        Line::new(zero.clone(), one.clone(), zero.sub(&vl, w)?, w)?,
        Line::new(zero.clone(), one, zero.sub(&vh, w)?, w)?,
    ];
    lines.extend_from_slice(cuts);
    let mut points = vec![
        Point {
            q: ql.clone(),
            v: vl.clone(),
            active: 0b0101,
        },
        Point {
            q: qh.clone(),
            v: vl,
            active: 0b0110,
        },
        Point {
            q: qh,
            v: vh.clone(),
            active: 0b1010,
        },
        Point {
            q: ql,
            v: vh,
            active: 0b1001,
        },
    ];
    for index in 4..lines.len() {
        w.step()?;
        // This API returns positive-area sampling cells. Tangencies with at
        // most two vertices have zero source area and cannot acquire area by
        // another halfplane intersection. They are not positive radial strips.
        if points.len() < 3 {
            points.clear();
            break;
        }
        let mut out = Vec::with_capacity(6);
        let mut previous = points.last().unwrap().clone();
        let mut before = lines[index].value(&previous, w)?;
        for mut current in points {
            let after = lines[index].value(&current, w)?;
            // The positive side (including its exact boundary) is retained.
            if (before == Ordering::Less && after == Ordering::Greater)
                || (before == Ordering::Greater && after == Ordering::Less)
            {
                let shared = previous.active & current.active;
                // Each polygon edge inherits an original rectangle support
                // or an earlier clipping line, hence shared is nonzero.
                let support = shared.trailing_zeros() as usize;
                let mut crossing = lines[support].intersection(&lines[index], w)?;
                crossing.active = (1 << support) | (1 << index);
                mark(&mut crossing, &lines, w)?;
                out.push(crossing);
            }
            if after != Ordering::Less {
                if after == Ordering::Equal {
                    current.active |= 1 << index;
                }
                out.push(current.clone());
            }
            previous = current;
            before = after;
        }
        // Four distinct rectangle corners plus at most two halfplanes form
        // a convex polygon with <=6 vertices. Strict crossings exclude zero
        // endpoints and every original vertex is retained only once, so no
        // positive-area polygon duplicates endpoints or needs proximity repair.
        points = out;
    }
    if points.len() < 3 {
        points.clear();
    }
    for p in &mut points {
        mark(p, &lines, w)?;
    }
    Ok((points, lines))
}
fn mark(p: &mut Point, lines: &[Line], w: &mut Work<'_>) -> Result<(), Error> {
    for (i, line) in lines.iter().enumerate() {
        if line.value(p, w)? == Ordering::Equal {
            p.active |= 1 << i;
        }
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    fn work() -> Work<'static> {
        Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        }
    }
    fn q(v: f64, w: &mut Work<'_>) -> Q {
        Q::float(v, w).unwrap()
    }
    #[test]
    fn exact_physical_seam_and_corner_incidence_survive_clipping() {
        let w = &mut work();
        // q+v>=.3, q+2v<=1.1; both constants are actual binary64
        // values, not decimal ideals or tolerance bands.
        let a = Line::new(q(1., w), q(1., w), q(-0.3, w), w).unwrap();
        let b = Line::new(q(-1., w), q(-2., w), q(1.1, w), w).unwrap();
        let (p, lines) = polygon([q(0., w), q(1., w), q(0., w), q(1., w)], &[a, b], w).unwrap();
        assert_eq!(p.len(), 5);
        for vertex in &p {
            assert!(lines[4].value(vertex, w).unwrap() != Ordering::Less);
            assert!(lines[5].value(vertex, w).unwrap() != Ordering::Less);
            assert!(vertex.active.count_ones() >= 2);
            for (i, line) in lines.iter().enumerate() {
                assert_eq!(
                    vertex.active & (1 << i) != 0,
                    line.value(vertex, w).unwrap() == Ordering::Equal
                );
            }
        }
        let a = Line::new(q(1., w), q(1., w), q(-1., w), w).unwrap();
        let (p, _) = polygon([q(0., w), q(1., w), q(0., w), q(1., w)], &[a], w).unwrap();
        assert_eq!(p.len(), 3);
        assert!(p.iter().filter(|p| p.active.count_ones() == 3).count() == 2);
    }
    #[test]
    fn positive_subulp_source_strip_is_not_welded_and_cancel_precedes_work() {
        let w = &mut work();
        let gap = f64::from_bits(1.0_f64.to_bits() + 1);
        let a = Line::new(q(1., w), q(0., w), q(-1., w), w).unwrap();
        let b = Line::new(q(-1., w), q(0., w), q(gap, w), w).unwrap();
        let (p, _) = polygon([q(0., w), q(2., w), q(0., w), q(1., w)], &[a, b], w).unwrap();
        assert_eq!(p.len(), 4);
        let lo = p
            .iter()
            .find(|p| p.q.enclosure(w).unwrap().0 == 1.)
            .unwrap();
        let hi = p
            .iter()
            .find(|p| p.q.enclosure(w).unwrap().0 == gap)
            .unwrap();
        assert_eq!(lo.q.cmp(&hi.q, w).unwrap(), Ordering::Less);
        let mut cancelled = Work {
            remaining: 0,
            used: 0,
            interrupt: Some(&|| true),
        };
        assert!(matches!(
            Q::float(1., &mut cancelled),
            Err(Error::Cancelled)
        ));
    }
}
#[cfg(test)]
mod degenerate_tests {
    use super::*;
    #[test]
    fn exact_tangent_cell_has_no_area_after_further_clipping() {
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let mut f = |v| Q::float(v, w).unwrap();
        let rect = [f(0.), f(1.), f(0.), f(1.)];
        let tangent = Line::new(
            Q::integer(1, w).unwrap(),
            Q::integer(0, w).unwrap(),
            Q::integer(-1, w).unwrap(),
            w,
        )
        .unwrap();
        let cross = Line::new(
            Q::integer(0, w).unwrap(),
            Q::integer(1, w).unwrap(),
            Q::float(-0.5, w).unwrap(),
            w,
        )
        .unwrap();
        assert!(polygon(rect, &[tangent, cross], w).unwrap().0.is_empty());
        let zero = Q::integer(0, w).unwrap();
        let one = Q::integer(1, w).unwrap();
        assert!(matches!(
            polygon(
                [zero.clone(), zero.clone(), zero.clone(), one.clone()],
                &[],
                w
            ),
            Err(Error::Unresolved)
        ));
        let line = Line::new(one.clone(), zero.clone(), zero.clone(), w).unwrap();
        assert!(matches!(
            polygon(
                [zero.clone(), one.clone(), zero, one],
                &[line.clone(), line.clone(), line],
                w
            ),
            Err(Error::Limit)
        ));
    }
}
