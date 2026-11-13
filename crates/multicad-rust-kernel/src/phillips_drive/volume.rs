//! Authored shallow-V section volume with finite-domain eligibility proofs.
use super::*;
use crate::tapered_thread::interval::J;
use std::collections::BinaryHeap;
fn j(v: I) -> J {
    J {
        v,
        ..J::constant(0.)
    }
}
fn n(v: f64) -> J {
    J::constant(v)
}
fn sq(v: J) -> J {
    v * v
}
fn root(v: J) -> Option<J> {
    if v.v.lo <= 0. {
        None
    } else {
        v.sqrt().checked().ok()
    }
}
struct Sections<'a> {
    c: &'a Cutter,
    b0: I,
    b1: I,
    offset: I,
}
impl<'a> Sections<'a> {
    fn new(c: &'a Cutter) -> Result<Self, Error> {
        let k = &c.corners[0];
        let SourceInner::ShallowV { sine, cosine } = k.inner else {
            return Err(Error::Unresolved);
        };
        // Section fronts have slopes b0,b1. With both below one, front(t)-t
        // decreases continuously, including a nonconvex shallow-V junction.
        let b0 = cosine / sine / k.cos_beta;
        let b1 = k.cos_theta / k.sin_theta / k.cos_beta;
        let offset = (k.join_s - k.join_y * k.cos_theta / k.sin_theta) / k.cos_beta;
        let e = &c.envelope;
        let tb = k.sin_beta / k.cos_beta;
        let rb = k.b_half + (k.b_depth - e.tip) * tb;
        // Every removed section point has s>=0. Bound its maximum over the
        // WHOLE cavity box, then prove the authored finite closure and both
        // extrusion caps are inactive before using an infinite-front integral.
        let smax = (k.cos_beta * (I::new(-e.entry.hi, e.entry.hi) - k.b_half)
            + k.sin_beta * (I::new(0., e.tip.hi) - k.b_depth))
            .hi
            .max(0.);
        // If adjacent rotated corners overlapped, monotonicity would put an
        // overlap on the diagonal with t<=entry/sqrt(2). Test that worst t and
        // minimum rb together, without losing the correlation f(t)-t.
        let t = e.entry / I::point(2.).sqrt();
        let f = if t.hi < k.join_y.lo {
            t * cosine / sine
        } else if t.lo > k.join_y.hi {
            k.join_s + (t - k.join_y) * k.cos_theta / k.sin_theta
        } else {
            return Err(Error::Unresolved);
        };
        if b0.lo <= 0.
            || b1.lo <= 0.
            || b0.hi >= 1.
            || b1.hi >= 1.
            || rb.lo <= 0.
            || (rb + offset).lo <= 0.
            || (rb + f / k.cos_beta - t).lo <= 0.
            || smax >= k.far_s.lo
            || (-I::point(c.domain.corner_start_depth_mm) - I::point(smax) * k.sin_beta).lo <= 0.
            || e.tip.hi >= c.domain.corner_end_depth_mm
        {
            return Err(Error::Unresolved);
        }
        Ok(Self { c, b0, b1, offset })
    }
    fn area(&self, d: J, branch: usize) -> Option<J> {
        let e = &self.c.envelope;
        let s = self.c.spec;
        let k = &self.c.corners[0];
        let r = match branch {
            0 => {
                j(e.entry)
                    - root(j(I::point(s.mouth_radius_mm).square()) - sq(d - n(s.mouth_radius_mm)))?
            }
            1 => j(I::point(s.mouth_diameter_mm) / I::point(2.)) - d * j(e.sin_a / e.cos_a),
            _ => (j(e.tip) - d) * j(e.cos_c / e.sin_c),
        };
        if r.v.lo <= 0. {
            return None;
        }
        let rb = j(k.b_half) + (j(k.b_depth) - d) * j(k.sin_beta / k.cos_beta);
        let a = part(r, rb, j(self.b0), I::point(0.), Some(k.join_y))?;
        let b = part(r, rb + j(self.offset), j(self.b1), k.join_y, None)?;
        let value = j(I::new(PI.next_down(), PI.next_up())) * sq(r) - n(8.) * (a + b);
        value.checked().ok()
    }
}
fn primitive(r: J, y: J) -> Option<J> {
    if y.v.lo == 0. && y.v.hi == 0. {
        return Some(n(0.));
    }
    let q = root(sq(r) - sq(y))?;
    Some((y * q + sq(r) * (y / q).atan()) / n(2.))
}
fn part(r: J, a: J, b: J, lo: I, upper: Option<I>) -> Option<J> {
    if r.v.hi <= a.v.lo {
        return Some(n(0.));
    }
    if r.v.lo <= a.v.hi {
        return None;
    }
    let disc = sq(r) * (n(1.) + sq(b)) - sq(a);
    let mut hi = (-a * b + root(disc)?) / (n(1.) + sq(b));
    // Whole-cell order tests only. At any unresolved intersection, clamp,
    // tangency, or inactive/active event the caller uses an area range instead
    // of differentiating across a branch change.
    if let Some(u) = upper {
        if hi.v.lo >= u.hi {
            hi = j(u)
        } else if hi.v.hi > u.lo {
            return None;
        }
    }
    if hi.v.hi <= lo.lo {
        return Some(n(0.));
    }
    if hi.v.lo <= lo.hi {
        return None;
    }
    Some(
        primitive(r, hi)?
            - primitive(r, j(lo))?
            - a * (hi - j(lo))
            - b * (sq(hi) - j(lo.square())) / n(2.),
    )
}
#[derive(Clone, Copy)]
struct Cell {
    lo: f64,
    hi: f64,
    v: I,
    numeric: I,
    quadrature: I,
}
fn cell(s: &Sections<'_>, lo: f64, hi: f64) -> Cell {
    let e = &s.c.envelope;
    let h = I::point(hi) - I::point(lo);
    let mid = (I::point(lo) + I::point(hi)) / I::point(2.);
    let branch = if hi <= e.join.lo {
        Some(0)
    } else if lo >= e.join.hi && hi <= e.flank.lo {
        Some(1)
    } else if lo >= e.flank.hi && hi <= e.tip.lo {
        Some(2)
    } else {
        None
    };
    let max = I::new(PI.next_down(), PI.next_up()) * e.entry.square() * h;
    if let Some(b) = branch
        && let (Some(whole), Some(m)) = (s.area(J::d(I::new(lo, hi)), b), s.area(J::d(mid), b))
    {
        return evaluated_cell(lo, hi, m.v * h, whole.dd, h, max);
    }
    range_cell(lo, hi, max)
}
fn range_cell(lo: f64, hi: f64, maximum: I) -> Cell {
    Cell {
        lo,
        hi,
        v: I::new(0., maximum.hi),
        numeric: I::point(0.),
        quadrature: I::new(0., maximum.hi),
    }
}
fn evaluated_cell(lo: f64, hi: f64, numeric: I, second: I, h: I, maximum: I) -> Cell {
    if let Some((v, remainder)) = midpoint_enclosure(numeric, second, h, maximum) {
        Cell {
            lo,
            hi,
            v,
            numeric: I::new(0., numeric.width()),
            quadrature: I::new(0., remainder.width()),
        }
    } else {
        range_cell(lo, hi, maximum)
    }
}
// The midpoint Peano kernel is nonnegative and integrates to h^3/24.
// Caller proves the C2 branch; arithmetic failure returns to whole-area range.
fn midpoint_enclosure(value: I, second: I, h: I, maximum: I) -> Option<(I, I)> {
    let remainder = second * h * h * h / I::point(24.);
    let v = value + remainder;
    let clipped = I::new(v.lo.max(0.), v.hi.min(maximum.hi));
    if v.valid() && clipped.valid() {
        Some((clipped, remainder))
    } else {
        None
    }
}
fn subdivision_point(lo: f64, hi: f64) -> Result<f64, Error> {
    let mid = lo + (hi - lo) / 2.;
    if mid <= lo || mid >= hi {
        Err(Error::Unresolved)
    } else {
        Ok(mid)
    }
}
pub(super) fn compute(
    c: &Cutter,
    width: f64,
    limit: usize,
    cancel: &dyn Fn() -> bool,
) -> Result<VolumeBounds, Error> {
    if !width.is_finite() || width <= 0. {
        return Err(Error::Invalid);
    }
    if cancel() {
        return Err(Error::Cancelled);
    }
    if limit > 1_000_000 {
        return Err(Error::Limit);
    }
    let s = Sections::new(c)?;
    let e = &c.envelope;
    let mut knots = vec![
        0., e.join.lo, e.join.hi, e.flank.lo, e.flank.hi, e.tip.lo, e.tip.hi,
    ];
    knots.sort_by(f64::total_cmp);
    knots.dedup();
    if limit < knots.len() - 1 {
        return Err(Error::Limit);
    }
    let mut cells: Vec<_> = knots.windows(2).map(|w| cell(&s, w[0], w[1])).collect();
    let mut heap: BinaryHeap<_> = cells
        .iter()
        .enumerate()
        .map(|(i, c)| (c.v.width().to_bits(), i))
        .collect();
    let mut work = cells.len();
    let mut estimate: f64 = cells.iter().map(|c| c.v.width()).sum();
    loop {
        if cancel() {
            return Err(Error::Cancelled);
        }
        if estimate <= width || work % 512 == 0 {
            let mut total = I::point(0.);
            for (index, c) in cells.iter().enumerate() {
                if index % 128 == 0 && cancel() {
                    return Err(Error::Cancelled);
                }
                total = total + c.v;
            }
            total.lo = total.lo.max(0.);
            if total.width() <= width {
                return Ok(VolumeBounds {
                    volume_mm3: bounds(total),
                    section_evaluations: work,
                    midpoint_evaluation_width_mm3: cells
                        .iter()
                        .fold(I::point(0.), |v, c| v + c.numeric)
                        .hi,
                    curvature_and_event_width_mm3: cells
                        .iter()
                        .fold(I::point(0.), |v, c| v + c.quadrature)
                        .hi,
                });
            }
            estimate = cells.iter().map(|c| c.v.width()).sum();
        }
        if work + 2 > limit {
            return Err(Error::Limit);
        }
        let (_, index) = heap.pop().ok_or(Error::Unresolved)?;
        let old = cells[index];
        let mid = subdivision_point(old.lo, old.hi)?;
        let left = cell(&s, old.lo, mid);
        let right = cell(&s, mid, old.hi);
        estimate += left.v.width() + right.v.width() - old.v.width();
        cells[index] = left;
        heap.push((left.v.width().to_bits(), index));
        heap.push((right.v.width().to_bits(), cells.len()));
        cells.push(right);
        work += 2;
    }
}

#[cfg(test)]
mod tests;
