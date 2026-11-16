use super::*;
use num_rational::BigRational as Q;
fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn choose(n: usize, k: usize) -> i64 {
    if k > n {
        0
    } else {
        (0..k).fold(1i64, |v, i| v * (n - i) as i64 / (i + 1) as i64)
    }
}
fn enclosed(actual: I, expected: &Q) {
    assert!(q(actual.lo) <= *expected && *expected <= q(actual.hi));
}
// Exact rational power-basis substitution, independently of interval
// de Casteljau. It starts from the ACTUAL binary64 control coefficients.
fn oracle(poles: &[P], axis: usize, a: f64, width: f64) -> Vec<Q> {
    let n = poles.len() - 1;
    let mut power = vec![q(0.); n + 1];
    for (i, p) in poles.iter().enumerate() {
        for (m, coefficient) in power.iter_mut().enumerate().skip(i) {
            *coefficient += q(p[axis])
                * Q::from_integer(
                    (choose(n, i) * choose(n - i, m - i) * if (m - i) % 2 == 0 { 1 } else { -1 })
                        .into(),
                );
        }
    }
    let transformed: Vec<_> = (0..=n)
        .map(|k| {
            (k..=n).fold(q(0.), |v, m| {
                v + &power[m]
                    * Q::from_integer(choose(m, k).into())
                    * q(a).pow((m - k) as i32)
                    * q(width).pow(k as i32)
            })
        })
        .collect();
    (0..=n)
        .map(|j| {
            (0..=j).fold(q(0.), |v, k| {
                v + &transformed[k] * Q::new(choose(j, k).into(), choose(n, k).into())
            })
        })
        .collect()
}
fn work() -> Work<'static> {
    Work {
        remaining: MESH_WORK,
        held: 0,
        peak: 0,
        maximum_bytes: MESH_BYTES,
        maximum_triangles: 1_000_000,
        delegated: Default::default(),
        interrupt: None,
    }
}
fn falling(j: usize, k: usize) -> f64 {
    if j < k {
        0.
    } else {
        (0..k).fold(1., |v, i| v * (j - i) as f64)
    }
}
// An independently specified nontrivial polynomial and its analytic derivative.
// Bernstein controls of u^k are (j)_k/(n)_k, independently of de Casteljau.
fn controls<const N: usize>() -> [[f64; 3]; N] {
    std::array::from_fn(|j| {
        let monomial = |k| falling(j, k) / falling(N - 1, k);
        [2. + monomial(1), 1. + monomial(2), 0.25 - monomial(3)]
    })
}
fn polynomial(t: f64) -> P {
    [2. + t, 1. + t * t, 0.25 - t * t * t]
}
#[test]
fn degree_eleven_actual_curve_has_independent_position_and_derivative_oracles() {
    let poles = controls::<12>();
    let curve = AngularCurve::<12> {
        poles: &poles,
        spans: 1,
        reverse: false,
        range: [0.25, 0.75],
    };
    let mut budget = work();
    for i in 0..=16 {
        let t = 0.25 + i as f64 / 32.;
        let b = curve.gamma([t, t], &mut budget).unwrap();
        for (axis, (b, p)) in b.into_iter().zip(polynomial(t)).enumerate() {
            enclosed(b, &oracle(&poles, axis, t, 0.)[0]);
            // Actual literal coefficients, not the ideal rational polynomial,
            // define this curve. Their rounded-control discrepancy is bounded
            // by convexity; the independent polynomial lies within 2e-14.
            assert!(b.lo <= p + 2e-14 && p - 2e-14 <= b.hi);
            assert!(b.hi - b.lo < 2e-12);
        }
    }
    let literal = poles.map(|p| p.map(I::point));
    let (left, _) = curve_split(literal, I::point(0.75));
    let (_, restricted) = curve_split(left, I::point(1.) / I::point(3.));
    let rational: Vec<_> = (0..3).map(|axis| oracle(&poles, axis, 0.25, 0.5)).collect();
    for j in 0..12 {
        let u = j as f64 / 11.;
        let u2 = falling(j, 2) / falling(11, 2);
        let u3 = falling(j, 3) / falling(11, 3);
        let expected = [
            2.25 + 0.5 * u,
            1.0625 + 0.25 * u + 0.25 * u2,
            0.234375 - 0.09375 * u - 0.1875 * u2 - 0.125 * u3,
        ];
        for k in 0..3 {
            enclosed(restricted[j][k], &rational[k][j]);
            assert!((restricted[j][k] - I::point(expected[k])).magnitude() < 2e-12);
        }
    }
    // Derivatives are in the restricted unit interval, including width=1/2.
    for j in 0..11 {
        let u = j as f64 / 10.;
        let u2 = falling(j, 2) / falling(10, 2);
        let expected = [0.5, 0.25 + 0.5 * u, -0.09375 - 0.375 * u - 0.375 * u2];
        for k in 0..3 {
            let derivative = I::point(11.) * (restricted[j + 1][k] - restricted[j][k]);
            enclosed(
                derivative,
                &((&rational[k][j + 1] - &rational[k][j]) * Q::from_integer(11.into())),
            );
            assert!((derivative - I::point(expected[k])).magnitude() < 5e-11);
        }
    }
    let reverse = AngularCurve {
        reverse: true,
        ..curve
    };
    for t in [0., 0.25, 0.75, 1.] {
        let a = curve.gamma([t, t], &mut budget).unwrap();
        let b = reverse.gamma([1. - t; 2], &mut budget).unwrap();
        for k in 0..3 {
            assert!((a[k] - b[k]).magnitude() < 2e-12);
        }
    }
}
#[test]
fn degree_eleven_curve_work_and_cancellation_are_not_geometric_refinement() {
    let poles = controls::<12>();
    let curve = AngularCurve::<12> {
        poles: &poles,
        spans: 1,
        reverse: false,
        range: [0., 1.],
    };
    let mut exact = work();
    let expected = curve.gamma([0.25, 0.75], &mut exact).unwrap();
    let used = MESH_WORK - exact.remaining;
    assert_eq!(used, 1 + 16 * 12 * 12);
    let mut retry = work();
    retry.remaining = used - 1;
    assert!(matches!(
        curve.gamma([0.25, 0.75], &mut retry),
        Err(PhillipsReadError::Limit)
    ));
    retry.remaining = used;
    let actual = curve.gamma([0.25, 0.75], &mut retry).unwrap();
    assert_eq!(actual.map(|i| [i.lo, i.hi]), expected.map(|i| [i.lo, i.hi]));
    let cancel = || true;
    retry.interrupt = Some(&cancel);
    assert!(matches!(
        curve.gamma([0.25, 0.75], &mut retry),
        Err(PhillipsReadError::Cancelled)
    ));
}

#[test]
fn degree_nine_existing_restriction_keeps_identical_interval_bits_and_work() {
    fn original(mut p: [[I; 3]; 10], t: I) -> ([[I; 3]; 10], [[I; 3]; 10]) {
        let mut left = [[I::point(0.); 3]; 10];
        let mut right = left;
        left[0] = p[0];
        right[9] = p[9];
        for level in 1..10 {
            for i in 0..10 - level {
                p[i] = std::array::from_fn(|k| (I::point(1.) - t) * p[i][k] + t * p[i + 1][k]);
            }
            left[level] = p[0];
            right[9 - level] = p[9 - level];
        }
        (left, right)
    }
    let poles = controls::<10>();
    for t in [0., 0.125, 0.3, 0.75, 1.] {
        let p = poles.map(|p| p.map(I::point));
        let a = original(p, I::point(t));
        let b = curve_split(p, I::point(t));
        let bits = |p: [[I; 3]; 10]| p.map(|p| p.map(|i| [i.lo.to_bits(), i.hi.to_bits()]));
        assert_eq!((bits(a.0), bits(a.1)), (bits(b.0), bits(b.1)));
    }
    let curve = Angular {
        poles: &poles,
        spans: 1,
        reverse: false,
        range: [0., 1.],
    };
    let mut budget = work();
    curve.gamma([0.25, 0.75], &mut budget).unwrap();
    assert_eq!(MESH_WORK - budget.remaining, 1025);
}
