//! Sign-consistent isolation of the crossings of a twice differentiable
//! function from its exact value, derivative and a second-derivative bound.
//!
//! The caller supplies the end values, so neighboring intervals that share
//! an end see one sign there and every sign change between two shared ends
//! is reported an odd number of times. A zero counts as positive; a tangent
//! contact changes no sign and reports nothing.

pub(super) fn positive(value: f64) -> bool {
    value >= 0.
}

/// Evaluations spent on one interval before the remaining pieces fall back
/// to the signs of their ends. Only a function that vanishes identically
/// to rounding, such as a plane containing a whole flat face, gets there.
const BUDGET: usize = 2048;

/// Bisect a sign change between `a` (with sign `from`) and `b`.
pub(super) fn bisect(eval: &dyn Fn(f64) -> (f64, f64), mut a: f64, mut b: f64, from: bool) -> f64 {
    for _ in 0..80 {
        let mid = a + (b - a) / 2.;
        if mid <= a || mid >= b {
            break;
        }
        if positive(eval(mid).0) == from {
            a = mid;
        } else {
            b = mid;
        }
    }
    a + (b - a) / 2.
}

/// Ascending parameters in `[a, b]` at which `eval` changes sign, given
/// `|f''| <= second` and the shared end values.
pub(super) fn crossings(
    eval: &dyn Fn(f64) -> (f64, f64),
    second: f64,
    [a, b]: [f64; 2],
    ends: [f64; 2],
    output: &mut Vec<f64>,
) {
    let mut budget = BUDGET;
    split(eval, second, [a, b], ends, &mut budget, output);
}

fn split(
    eval: &dyn Fn(f64) -> (f64, f64),
    second: f64,
    [a, b]: [f64; 2],
    [fa, fb]: [f64; 2],
    budget: &mut usize,
    output: &mut Vec<f64>,
) {
    let width = b - a;
    let mid = a + width / 2.;
    let changed = positive(fa) != positive(fb);
    if !(mid > a && mid < b) {
        if changed {
            output.push(a);
        }
        return;
    }
    if *budget == 0 {
        if changed {
            output.push(bisect(eval, a, b, positive(fa)));
        }
        return;
    }
    *budget -= 1;
    let (fm, dm) = eval(mid);
    // |f(x) - f(m)| <= |f'(m)| |x - m| + second |x - m|^2 / 2.
    if fm.abs() > dm.abs() * width / 2. + second * width * width / 8. {
        // No interior zero: a differing end sign is a zero at that end.
        match (positive(fa) != positive(fm), positive(fb) != positive(fm)) {
            (true, false) => output.push(a),
            (false, true) => output.push(b),
            _ => {}
        }
        return;
    }
    // f' keeps its sign: at most one zero.
    if dm.abs() >= second * width / 2. {
        if changed {
            output.push(bisect(eval, a, b, positive(fa)));
        }
        return;
    }
    split(eval, second, [a, mid], [fa, fm], budget, output);
    split(eval, second, [mid, b], [fm, fb], budget, output);
}

/// Real roots of `a x^2 + 2 b x + c`, ascending. A tangent or nearly
/// tangent line contributes its point of closest approach instead, which
/// only subdivides an interval of the caller.
pub(super) fn quadratic(a: f64, b: f64, c: f64, output: &mut Vec<f64>) {
    let scale = a.abs().max(b.abs()).max(c.abs());
    if scale == 0. {
        return;
    }
    if a.abs() <= 1e-14 * scale {
        if b.abs() > 1e-14 * scale {
            output.push(-c / (2. * b));
        }
        return;
    }
    let discriminant = b * b - a * c;
    if discriminant <= 1e-12 * (b * b).max((a * c).abs()) {
        output.push(-b / a);
        if discriminant <= 0. {
            return;
        }
    }
    let q = -(b + discriminant.sqrt().copysign(b));
    let mut roots = [q / a, if q == 0. { 0. } else { c / q }];
    if roots[0] > roots[1] {
        roots.swap(0, 1);
    }
    output.extend(roots);
}
