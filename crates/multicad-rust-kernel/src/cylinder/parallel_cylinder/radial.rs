//! Error-free binary products retain squared gaps below a distance ULP.
use super::*;

fn sum(a: f64, b: f64) -> [f64; 2] {
    let value = a + b;
    let split = value - a;
    [value, (a - (value - split)) + (b - split)]
}

pub(super) fn square_difference(xy: [f64; 2], radius: Rational) -> Option<I> {
    let r = radius.as_f64();
    // The radius here is the actual rational source, not its rounded display.
    if binary_fraction(r) != Some(radius) {
        return None;
    }
    // These bounds prevent overflow and keep every product/residual normal.
    // At 1e-100 the least possible exact binary product bit is above 2^-770;
    // FMA therefore retains the complete multiplication residual. Smaller
    // nonzero operands must use another proof or explicitly remain unresolved.
    let values = [xy[0], xy[1], r];
    if values
        .into_iter()
        .any(|v| !v.is_finite() || v.abs() > 2e7 || (v != 0. && v.abs() < 1e-100))
    {
        return None;
    }
    let products = values.map(|v| [v * v, v.mul_add(v, -(v * v))]);
    let [s, se] = sum(products[0][0], products[1][0]);
    let [gap, ge] = sum(s, -products[2][0]);
    let mut result = I::point(gap);
    for residual in [ge, se, products[0][1], products[1][1], -products[2][1]] {
        if residual != 0. {
            result = result + I::point(residual);
        }
    }
    Some(result)
}
