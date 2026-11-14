//! Local literal tensor projection only. Closure, fan links, cap winding and
//! occupied material remain SEPARATE obligations; no STEP admission here.
use super::*;
pub(super) mod polynomial;
use polynomial::{P, cross, derivative, dot_xy};
#[cfg(test)]
mod tests;
const DEPTH: usize = 32;
const SCRATCH: usize = 256 * 1024;
#[derive(Clone, Copy)]
struct Cell {
    p: [P; 3],
    depth: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct ProjectionBounds {
    pub component_radius_mm: [[f64; 2]; 2],
    pub minimum_signed_jacobian: f64,
    pub certified_cells: usize,
    pub work_used: usize,
}
impl SurfacePlan<'_, '_> {
    /// Conditional radial-chart proof from the ACTUAL shared literal controls.
    /// Work/bytes include the entire retained row/domain/source packet once.
    pub fn projection_bounds(
        &self,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ProjectionBounds, Error> {
        let mut w = Work::new(maximum_work, interrupt)?;
        w.charge(self.work_used() + 1)?;
        if maximum_bytes == 0 || maximum_bytes > 64 * 1024 * 1024 {
            return Err(Error::Invalid);
        }
        let mut memory = Memory::new(self.retained_bytes(), maximum_bytes)?;
        memory.add(SCRATCH)?;
        let mut stack = memory.allocate(3 * DEPTH + 1, &mut w)?;
        let mut radii = [[f64::INFINITY, 0_f64]; 2];
        let mut minimum = f64::INFINITY;
        let mut cells = 0;
        for patch in self.patches() {
            let slab = &self.domain().slabs()[patch.slab];
            let [a, b] = patch.rows.map(|i| &self.rows()[i]);
            let spans = a.knots().len() - 1;
            let height = self.domain().boundary().source().full_length_mm();
            height_controls(a.poles(), height, &mut w)?;
            height_controls(b.poles(), height, &mut w)?;
            let hand = self
                .domain()
                .boundary()
                .source_boundary(slab.side)
                .source()
                .spec()
                .left_hand;
            for span in 0..spans {
                surfaces::verify_boundaries(slab, a, b, self.vertices())?;
                let reverse = a.endpoints()[1] == b.endpoints()[1] && span + 1 == spans;
                let factor = (a.endpoints()[0] == b.endpoints()[0] && span == 0) || reverse;
                let controls = std::array::from_fn::<_, 2, _>(|row| {
                    std::array::from_fn::<_, 3, _>(|axis| {
                        let poles = if row == 0 { a.poles() } else { b.poles() };
                        let values = std::array::from_fn::<_, 12, _>(|i| {
                            literal(poles[11 * span + if reverse { 11 - i } else { i }][axis])
                        });
                        P::controls(&values)
                    })
                });
                let sign = if hand != reverse { 1. } else { -1. };
                let p = jacobian(controls, factor, &mut w)?;
                let j = [
                    p[0].scale(sign, &mut w)?,
                    p[1].scale(sign, &mut w)?,
                    p[2].scale(sign, &mut w)?,
                ];
                let (bound, n) = positive(j, &mut stack, &mut w)?;
                minimum = minimum.min(bound.lo);
                cells += n;
                let r = radius_squared(controls, &mut w)?;
                let (bound, n) = positive(r, &mut stack, &mut w)?;
                cells += n;
                let bound = bound.sqrt();
                // positive() returns a finite strictly positive enclosure;
                // sqrt preserves finiteness and positivity on that domain.
                let i = if slab.side == Side::Outer { 0 } else { 1 };
                radii[i] = [radii[i][0].min(bound.lo), radii[i][1].max(bound.hi)];
            }
        }
        w.charge(1)?;
        Ok(ProjectionBounds {
            component_radius_mm: radii,
            minimum_signed_jacobian: minimum,
            certified_cells: cells,
            work_used: w.used,
        })
    }
}
fn literal(x: f64) -> I {
    if x == 0. || (x.fract() == 0. && x.abs() <= 9_007_199_254_740_992.) {
        I::point(x)
    } else {
        I::new(x.next_down(), x.next_up())
    }
}
fn jacobian([a, b]: [[P; 3]; 2], factor: bool, w: &mut Work<'_>) -> Result<[P; 3], Error> {
    let mut d = [a[0]; 3];
    for axis in 0..3 {
        d[axis] = b[axis].sub(a[axis], w)?;
    }
    if factor {
        for p in &mut d {
            p.c[0] = I::point(0.);
        }
    }
    let original = d;
    if factor {
        // BOTH full U-row controls reference the SAME owned anchor. Their
        // difference's constant coefficient is an identity, not overlapping I.
        for axis in &mut d {
            axis.c[0] = I::point(0.);
            *axis = axis.remove_t(w)?;
        }
    }
    let n0 = cross(d, derivative(a, w)?, w)?;
    let mut n1 = cross(d, derivative(d, w)?, w)?;
    if factor {
        for p in &mut n1 {
            *p = p.times_t(w)?;
        }
    }
    let f0 = dot_xy(a, n0, w)?;
    let f1 = dot_xy(a, n1, w)?.add(dot_xy(original, n0, w)?, w)?;
    let f2 = dot_xy(original, n1, w)?;
    Ok([
        f0,
        f0.add(f1.scale(0.5, w)?, w)?,
        f0.add(f1, w)?.add(f2, w)?,
    ])
}
fn radius_squared([a, b]: [[P; 3]; 2], w: &mut Work<'_>) -> Result<[P; 3], Error> {
    Ok([dot_xy(a, a, w)?, dot_xy(a, b, w)?, dot_xy(b, b, w)?])
}
fn positive(p: [P; 3], stack: &mut Vec<Cell>, w: &mut Work<'_>) -> Result<(I, usize), Error> {
    stack.push(Cell { p, depth: 0 });
    let mut total = I::new(f64::INFINITY, f64::NEG_INFINITY);
    let mut cells = 0;
    while let Some(cell) = stack.pop() {
        w.charge(3 * (cell.p[0].degree + 1))?;
        let h = cell.p.map(P::hull);
        if h.iter().any(|x| !x.valid()) {
            return Err(Error::Unresolved);
        }
        let lo = h.iter().map(|x| x.lo).fold(f64::INFINITY, f64::min);
        let hi = h.iter().map(|x| x.hi).fold(f64::NEG_INFINITY, f64::max);
        if lo > 0. {
            total.lo = total.lo.min(lo);
            total.hi = total.hi.max(hi);
            cells += 1;
            continue;
        }
        if hi <= 0. || cell.depth == DEPTH {
            return Err(Error::Unresolved);
        }
        let p = cell.p;
        let left = p[0].add(p[1], w)?.scale(0.5, w)?;
        let right = p[1].add(p[2], w)?.scale(0.5, w)?;
        let mid = left.add(right, w)?.scale(0.5, w)?;
        for u in [[p[0], left, mid], [mid, right, p[2]]] {
            let a = u[0].split(w)?;
            let b = u[1].split(w)?;
            let c = u[2].split(w)?;
            for i in 0..2 {
                stack.push(Cell {
                    p: [a[i], b[i], c[i]],
                    depth: cell.depth + 1,
                });
            }
        }
    }
    Ok((total, cells))
}

// The new writer encodes every coefficient and the cap height with the same
// monotone round-tripping decimal formatter. Equal coefficients are identical
// numeric leaves; the literal convex hull therefore remains in these cap planes.
fn height_controls(poles: &[[f64; 3]], height: f64, w: &mut Work<'_>) -> Result<(), Error> {
    for p in poles {
        w.charge(1)?;
        if !p[2].is_finite() || p[2] < 0. || p[2] > height {
            return Err(Error::Unresolved);
        }
    }
    Ok(())
}
