//! Isolation of every root of a Lipschitz function on an interval.
use super::super::*;

pub(super) fn limit() -> GeometryError {
    err(
        GeometryErrorKind::Limit,
        "Rust twisted gear drawing exceeded its root isolation budget",
    )
}

/// Bisect a sign change to adjacent floating-point parameters.
pub(super) fn bisect(f: &mut dyn FnMut(f64) -> f64, mut a: f64, fa: f64, mut b: f64) -> f64 {
    if fa == 0. {
        return a;
    }
    for _ in 0..200 {
        let mid = 0.5 * (a + b);
        if mid <= a.min(b) || mid >= a.max(b) {
            break;
        }
        let fm = f(mid);
        if fm == 0. {
            return mid;
        }
        if (fm < 0.) == (fa < 0.) {
            a = mid;
        } else {
            b = mid;
        }
    }
    0.5 * (a + b)
}

pub(super) struct Search {
    /// Bound of the function's slope on the interval.
    pub slope: f64,
    /// Width below which a sign change is one simple root. Roots that leave
    /// the sign unchanged within a leaf are grazing contacts.
    pub leaf: f64,
    /// Magnitude below which the function only grazes zero: a stretch on
    /// which it stays this small contributes no root.
    pub noise: f64,
}

/// Every sign change of `f` in `[lo, hi]`. An interval is discarded when the
/// slope bound proves it free of roots, or when the function stays within
/// the grazing magnitude at its ends and middle without changing sign.
/// Values may be wrapped angles: a wrap is excluded by the same bound.
pub(super) fn isolate(
    f: &mut dyn FnMut(f64) -> f64,
    search: &Search,
    [lo, hi]: [f64; 2],
    budget: &mut usize,
    output: &mut Vec<f64>,
) -> Result<(), GeometryError> {
    if hi <= lo || !(hi - lo).is_finite() {
        return Ok(());
    }
    let mut stack = vec![(lo, f(lo), hi, f(hi))];
    while let Some((a, fa, b, fb)) = stack.pop() {
        *budget = budget.checked_sub(1).ok_or_else(limit)?;
        let mid = 0.5 * (a + b);
        let fm = f(mid);
        let half = 0.5 * (b - a);
        if fm.abs() > search.slope * half * (1. + 1e-9) + search.noise {
            continue;
        }
        if b - a <= search.leaf {
            if fa == 0. || fb == 0. || (fa < 0.) != (fb < 0.) {
                output.push(bisect(f, a, fa, b));
            }
            continue;
        }
        if [fa, fm, fb].iter().all(|v| v.abs() <= search.noise) {
            // The function stays within the grazing magnitude: at most one
            // crossing is resolved here.
            if (fa < 0.) != (fb < 0.) {
                output.push(bisect(f, a, fa, b));
            }
            continue;
        }
        stack.push((mid, fm, b, fb));
        stack.push((a, fa, mid, fm));
    }
    Ok(())
}
