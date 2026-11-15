//! Certified continuation of the zero set of a smooth function on the unit
//! square. Every step is a chart in which one partial derivative keeps its
//! sign, proved from second-derivative bounds, so the curve is a graph over
//! the other coordinate and is evaluated by solving the function itself.
use super::super::*;
use super::solve::{Search, bisect, isolate, limit};

/// Value and gradient at a point of the unit square.
pub(super) type Eval<'a> = &'a dyn Fn(f64, f64) -> (f64, [f64; 2]);

pub(in crate::drawing) struct Bounds {
    /// Bounds of the two first partial derivatives.
    pub first: [f64; 2],
    /// Bounds of the pure first, mixed and pure second partial derivatives.
    pub second: [f64; 3],
}

/// One chart: the coordinate `u` runs over `u[0]..u[1]` and the other lies
/// in the bracket `w`, on which the function is strictly monotone.
#[derive(Clone, Copy, Debug)]
pub(in crate::drawing) struct Step {
    pub by_first: bool,
    pub u: [f64; 2],
    pub w: [f64; 2],
}

fn singular() -> GeometryError {
    err(
        GeometryErrorKind::Unsupported,
        "Rust twisted gear drawing meets a flank tangentially; the exact curve is singular there",
    )
}

impl Step {
    fn value(&self, eval: Eval<'_>, u: f64, w: f64) -> f64 {
        if self.by_first {
            eval(u, w).0
        } else {
            eval(w, u).0
        }
    }

    /// The dependent coordinate at `u`.
    pub(super) fn solve(&self, eval: Eval<'_>, u: f64) -> f64 {
        let [mut a, mut b] = self.w;
        let [mut fa, mut fb] = [self.value(eval, u, a), self.value(eval, u, b)];
        if fa == 0. {
            return a;
        }
        if fb == 0. {
            return b;
        }
        if (fa < 0.) == (fb < 0.) {
            return if fa.abs() < fb.abs() { a } else { b };
        }
        // Illinois false position, falling back to bisection.
        let mut w = a + (b - a) * fa / (fa - fb);
        for _ in 0..200 {
            if !(w > a.min(b) && w < a.max(b)) {
                w = 0.5 * (a + b);
                if !(w > a.min(b) && w < a.max(b)) {
                    break;
                }
            }
            let f = self.value(eval, u, w);
            if f == 0. {
                return w;
            }
            if (f < 0.) == (fa < 0.) {
                a = w;
                fa = f;
                fb *= 0.5;
            } else {
                b = w;
                fb = f;
                fa *= 0.5;
            }
            if (b - a).abs() <= 2. * f64::EPSILON * (a.abs() + b.abs()) {
                break;
            }
            w = a + (b - a) * fa / (fa - fb);
        }
        if fa.abs() < fb.abs() { a } else { b }
    }

    pub(super) fn point(&self, eval: Eval<'_>, fraction: f64) -> [f64; 2] {
        let u = self.u[0] + (self.u[1] - self.u[0]) * fraction;
        let w = self.solve(eval, u);
        if self.by_first { [u, w] } else { [w, u] }
    }
}

/// A connected curve from one boundary point of the square to another.
pub(super) struct Chain {
    pub steps: Vec<Step>,
    pub ends: [[f64; 2]; 2],
}

const MERGE: f64 = 1e-6;

/// Every curve of the zero set that reaches the boundary of the square.
pub(super) fn trace(
    eval: Eval<'_>,
    bounds: &Bounds,
    budget: &mut usize,
) -> Result<Vec<Chain>, GeometryError> {
    let magnitude = bounds.first[0] + bounds.first[1];
    if !(magnitude.is_finite() && magnitude > 0.) {
        return Ok(Vec::new());
    }
    let mut roots: Vec<([f64; 2], bool)> = Vec::new();
    for axis in 0..2 {
        for side in [0., 1.] {
            let at = |x: f64| if axis == 0 { [side, x] } else { [x, side] };
            let mut found = Vec::new();
            isolate(
                &mut |x| {
                    let [s, z] = at(x);
                    eval(s, z).0
                },
                &Search {
                    slope: bounds.first[1 - axis],
                    leaf: 1e-9,
                    noise: 1e-12 * magnitude,
                    second: bounds.second[if axis == 0 { 2 } else { 0 }],
                    lift: 0.,
                    wrapped: false,
                },
                [0., 1.],
                budget,
                &mut found,
            )?;
            for x in found {
                let p = at(x);
                if !roots
                    .iter()
                    .any(|(q, _)| (q[0] - p[0]).hypot(q[1] - p[1]) <= MERGE)
                {
                    roots.push((p, false));
                }
            }
        }
    }
    let mut chains = Vec::new();
    for index in 0..roots.len() {
        if roots[index].1 {
            continue;
        }
        roots[index].1 = true;
        let start = roots[index].0;
        let (_, g) = eval(start[0], start[1]);
        let mut heading = [-g[1], g[0]];
        let size = heading[0].hypot(heading[1]);
        if size == 0. {
            return Err(singular());
        }
        // Into the square: away from every side the start lies on.
        let mut inward = 0_f64;
        let mut blocked = false;
        for k in 0..2 {
            let normal = if start[k] <= MERGE {
                1.
            } else if start[k] >= 1. - MERGE {
                -1.
            } else {
                continue;
            };
            let component = heading[k] * normal / size;
            if component.abs() <= 1e-9 {
                continue;
            }
            if inward == 0. {
                inward = component.signum();
            } else if inward != component.signum() {
                blocked = true;
            }
        }
        if inward == 0. || blocked {
            // The curve only touches the boundary here.
            continue;
        }
        heading = [heading[0] * inward, heading[1] * inward];
        let chain = walk(eval, bounds, start, heading, budget)?;
        let end = chain.ends[1];
        for (p, used) in &mut roots {
            if (p[0] - end[0]).hypot(p[1] - end[1]) <= MERGE {
                *used = true;
            }
        }
        let length: f64 = chain.steps.iter().map(|s| (s.u[1] - s.u[0]).abs()).sum();
        if length > 1e-9 {
            chains.push(chain);
        }
    }
    Ok(chains)
}

fn walk(
    eval: Eval<'_>,
    bounds: &Bounds,
    start: [f64; 2],
    mut heading: [f64; 2],
    budget: &mut usize,
) -> Result<Chain, GeometryError> {
    let mut steps = Vec::new();
    let mut p = start;
    loop {
        *budget = budget.checked_sub(16).ok_or_else(limit)?;
        let (_, g) = eval(p[0], p[1]);
        let by_first = g[1].abs() >= g[0].abs();
        let (iu, iw) = if by_first { (0, 1) } else { (1, 0) };
        let mut tangent = [-g[1], g[0]];
        if tangent[0] * heading[0] + tangent[1] * heading[1] < 0. {
            tangent = [-tangent[0], -tangent[1]];
        }
        if g[iw] == 0. || !g[iw].is_finite() {
            return Err(singular());
        }
        let direction = tangent[iu].signum();
        let available = if direction > 0. { 1. - p[iu] } else { p[iu] };
        if available <= 1e-13 {
            return Ok(Chain {
                steps,
                ends: [start, p],
            });
        }
        let [gu, gw] = [g[iu].abs(), g[iw].abs()];
        let [muu, muw, mww] = [
            bounds.second[if by_first { 0 } else { 2 }],
            bounds.second[1],
            bounds.second[if by_first { 2 } else { 0 }],
        ];
        let mut span = available.min(0.25);
        let reach = loop {
            let reach = 4. * (gu + muu * span) * span / gw + 1e-7;
            let floor = gw - muw * span - mww * reach;
            if muw * span + mww * reach <= 0.5 * gw
                && (gu + muu * span + muw * reach) * span <= floor * reach
            {
                break reach;
            }
            span *= 0.5;
            if span < 1e-10 {
                return Err(singular());
            }
        };
        let last = span == available;
        let target = if last {
            if direction > 0. { 1. } else { 0. }
        } else {
            p[iu] + direction * span
        };
        let mut step = Step {
            by_first,
            u: [p[iu], target],
            w: [p[iw] - reach, p[iw] + reach],
        };
        let w = step.solve(eval, target);
        if !(-1e-13..=1. + 1e-13).contains(&w) {
            // The curve leaves through a side of the dependent coordinate.
            let bound = if w < 0. { 0. } else { 1. };
            let exit = bisect(
                &mut |u| {
                    let inside = (step.solve(eval, u) - bound) * (0.5 - bound);
                    if inside >= 0. { 1. } else { -1. }
                },
                step.u[0],
                1.,
                target,
            );
            step.u[1] = exit;
            let mut end = [0.; 2];
            end[iu] = exit;
            end[iw] = bound;
            if (exit - step.u[0]).abs() > 0. {
                steps.push(step);
            }
            return Ok(Chain {
                steps,
                ends: [start, end],
            });
        }
        steps.push(step);
        p[iu] = target;
        p[iw] = w.clamp(0., 1.);
        heading = tangent;
        if last {
            return Ok(Chain {
                steps,
                ends: [start, p],
            });
        }
        if steps.len() > 100_000 {
            return Err(limit());
        }
    }
}
