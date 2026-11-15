//! Fixed-size outward normals of actual finite binary64 triangle vertices.
//! TwoDiff and FMA residuals avoid subtracting independently widened products.
use crate::tapered_thread::interval::I;
/// Conservative fixed work allowance for the bounded3-component calculation.
pub const NORMAL_WORK: usize = 512;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unresolved,
}
fn add(a: f64, b: f64) -> f64 {
    (a + b).next_up()
}
fn mul(a: f64, b: f64) -> f64 {
    (a * b).next_up()
}
fn difference(a: f64, b: f64) -> (f64, f64) {
    let value = a - b;
    let z = value - a;
    let tail = (a - (value - z)) - (b + z);
    (value, tail.abs().next_up())
}
/// Exact input vertices define the plane. No analytic source, normalization
/// target or tolerance can replace or adjust these literals. Coordinates are
/// bounded by1e7; every product/difference stays finite. Residual absolute
/// bounds include a subnormal rounding allowance, and unresolved tiny/zero
/// areas refuse rather than being declared positive.
pub fn normal_bounds(points: [[f64; 3]; 3]) -> Result<[[f64; 2]; 3], Error> {
    if points
        .iter()
        .flatten()
        .any(|v| !v.is_finite() || v.abs() > 1e7)
    {
        return Err(Error::Invalid);
    }
    let edges: [[(f64, f64); 3]; 2] =
        [1, 2].map(|j| std::array::from_fn(|i| difference(points[j][i], points[0][i])));
    let cross = std::array::from_fn::<_, 3, _>(|i| {
        let (j, k) = [(1, 2), (2, 0), (0, 1)][i];
        let (a, ea) = edges[0][j];
        let (b, eb) = edges[1][k];
        let (c, ec) = edges[0][k];
        let (d, ed) = edges[1][j];
        let p = a * b;
        let q = c * d;
        let (value, subtraction) = difference(p, q);
        let tails = add(
            a.mul_add(b, -p).abs().next_up(),
            c.mul_add(d, -q).abs().next_up(),
        );
        let first = add(add(mul(a.abs(), eb), mul(b.abs(), ea)), mul(ea, eb));
        let second = add(add(mul(c.abs(), ed), mul(d.abs(), ec)), mul(ec, ed));
        let error = add(add(tails, subtraction), add(first, second));
        I::new((value - error).next_down(), (value + error).next_up())
    });
    let length = cross
        .into_iter()
        .fold(I::point(0.), |sum, v| sum + v.square())
        .sqrt();
    if !length.lo.is_finite() || length.lo <= 0. || !length.hi.is_finite() {
        return Err(Error::Unresolved);
    }
    let normal = cross.map(|v| v / length);
    // Finite cross components are <1e15. A positive squared-norm lower
    // endpoint is at least the least subnormal, so the positive norm bound
    // exceeds2e-162. Every quotient is therefore finite and ordered.
    Ok(normal.map(|v| [v.lo, v.hi]))
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_tiny_translated_triangle_keeps_its_finite_normal() {
        let e = 2f64.powi(-26);
        for shift in [0., 4., 1_000_000.] {
            let p = [
                [shift, shift, shift],
                [shift + e, shift, shift],
                [shift, shift + e, shift],
            ];
            let n = normal_bounds(p).unwrap();
            assert!(n[2][0] > 0.999999999 && n[2][1] >= 1.);
            assert!(n[0][0] <= 0. && n[0][1] >= 0.);
            let mut reversed = p;
            reversed.swap(1, 2);
            let back = normal_bounds(reversed).unwrap();
            assert!(back[2][1] < -0.999999999);
        }
    }
    #[test]
    fn degenerate_subnormal_and_nonfinite_literals_refuse() {
        assert_eq!(normal_bounds([[0.; 3]; 3]), Err(Error::Unresolved));
        assert_eq!(normal_bounds([[f64::NAN; 3]; 3]), Err(Error::Invalid));
        assert_eq!(normal_bounds([[1e8; 3]; 3]), Err(Error::Invalid));
        assert_eq!(
            normal_bounds([
                [0.; 3],
                [f64::MIN_POSITIVE, 0., 0.],
                [0., f64::MIN_POSITIVE, 0.]
            ]),
            Err(Error::Unresolved)
        );
    }
}
#[cfg(test)]
mod exact_tests {
    use super::*;
    use num_rational::BigRational;
    fn rational(v: f64) -> BigRational {
        BigRational::from_float(v).unwrap()
    }
    #[test]
    fn actual_oblique_slivers_enclose_independent_exact_binary_cross_products() {
        for points in [
            [
                [4.08, 0., 4.099],
                [4.080000000000001, 1e-14, 4.099],
                [3.81, 0.31, 4.098],
            ],
            [
                [1e6, 1e6, 1e6],
                [1e6 + 2f64.powi(-20), 1e6 + 2f64.powi(-21), 1e6],
                [1e6, 1e6 + 2f64.powi(-20), 1e6 + 2f64.powi(-22)],
            ],
        ] {
            let bounds = normal_bounds(points).unwrap();
            let a =
                std::array::from_fn::<_, 3, _>(|i| rational(points[1][i]) - rational(points[0][i]));
            let b =
                std::array::from_fn::<_, 3, _>(|i| rational(points[2][i]) - rational(points[0][i]));
            let c = [
                &a[1] * &b[2] - &a[2] * &b[1],
                &a[2] * &b[0] - &a[0] * &b[2],
                &a[0] * &b[1] - &a[1] * &b[0],
            ];
            let sum = c.iter().map(|v| v * v).sum::<BigRational>();
            let zero = BigRational::from_integer(0.into());
            assert!(sum > zero);
            for i in 0..3 {
                let lo = rational(bounds[i][0]);
                let hi = rational(bounds[i][1]);
                let square = &c[i] * &c[i];
                if lo >= zero {
                    assert!(c[i] >= zero && &lo * &lo * &sum <= square)
                } else if c[i] < zero {
                    assert!(&lo * &lo * &sum >= square)
                }
                if hi <= zero {
                    assert!(c[i] <= zero && &hi * &hi * &sum <= square)
                } else if c[i] > zero {
                    assert!(&hi * &hi * &sum >= square)
                }
            }
        }
    }
}
