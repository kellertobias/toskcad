use super::*;
fn work() -> Work<'static> {
    Work {
        remaining: 1_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    }
}
fn circle() -> crate::helical_thread::RootCircle {
    crate::helical_thread::RootCircle {
        phase: 0.,
        center_radius: 3.,
        radius: 1.,
    }
}
fn exact_contains(v: I, n: i64, d: i64) {
    use num_rational::BigRational as Q;
    let exact = Q::new(n.into(), d.into());
    assert!(
        Q::from_float(v.lo).unwrap() <= exact && exact <= Q::from_float(v.hi).unwrap(),
        "{v:?} expected{n}/{d}"
    );
}
#[test]
fn independent_actual_circle_complex_derivative_coefficients_and_both_hands() {
    for hand in [1., -1.] {
        let phase = [I::point(-0.25), I::point(0.25)];
        let a = value(circle(), phase, I::point(1.), I::point(hand), I::point(0.5)).unwrap();
        exact_contains(a[0].c[0], 2, 1);
        exact_contains(a[1].c[0], 0, 1);
        exact_contains(a[0].c[1], 0, 1);
        exact_contains(a[1].c[1], -(hand as i64), 1);
        exact_contains(a[0].c[2], -1, 8);
        exact_contains(a[1].c[3], hand as i64, 48);
        let factor = bounds(circle(), phase, I::point(1.), I::point(hand), &mut work()).unwrap();
        assert_eq!(factor.poles.len(), 17);
        assert!(factor.radius >= 2. && factor.radius < 2.1);
        assert!(factor.du >= 1. && factor.du < 1.2);
        assert!(
            factor
                .remainder
                .into_iter()
                .all(|x| x >= 0. && x.is_finite())
        );
    }
}
#[test]
fn real_circle_whole_domain_and_resource_refusals() {
    let phase = [I::point(-0.25), I::point(0.25)];
    for (root, p, rate, hand) in [
        (
            circle(),
            [I::point(-2.), I::point(2.)],
            I::point(1.),
            I::point(1.),
        ),
        (circle(), phase, I::point(0.), I::point(1.)),
        (circle(), phase, I::point(1.), I::point(0.)),
        (circle(), [phase[1], phase[0]], I::point(1.), I::point(1.)),
        (
            crate::helical_thread::RootCircle {
                center_radius: 0.5,
                ..circle()
            },
            phase,
            I::point(1.),
            I::point(1.),
        ),
    ] {
        assert!(matches!(
            bounds(root, p, rate, hand, &mut work()),
            Err(ReadError::Unsupported)
        ));
    }
    let mut limited = Work {
        remaining: 131071,
        ..work()
    };
    assert!(matches!(
        bounds(circle(), phase, I::point(1.), I::point(1.), &mut limited),
        Err(ReadError::Limit)
    ));
}

#[test]
fn all_eighteen_angle_chart_coefficients_enclose_independent_exact_formal_series() {
    use num_bigint::BigInt;
    use num_rational::BigRational as Q;
    let zero = Q::from_integer(0.into());
    let one = Q::from_integer(1.into());
    let mut sin = vec![zero.clone(); 18];
    let mut cos = sin.clone();
    let mut factorial = BigInt::from(1u8);
    for k in 0..18 {
        if k > 0 {
            factorial *= k;
        }
        let den = &factorial * (BigInt::from(1u8) << k);
        if k % 2 == 0 {
            cos[k] = Q::new(
                if (k / 2) % 2 == 0 {
                    1.into()
                } else {
                    (-1).into()
                },
                den,
            );
        } else {
            sin[k] = Q::new(
                if ((k - 1) / 2) % 2 == 0 {
                    1.into()
                } else {
                    (-1).into()
                },
                den,
            );
        }
    }
    let mut radius = cos.iter().map(|x| -x).collect::<Vec<_>>();
    radius[0] += Q::from_integer(3.into());
    for hand in [1i64, -1] {
        let mut exp = vec![[zero.clone(), zero.clone()]; 18];
        exp[0][0] = one.clone();
        for n in 1..18 {
            for k in 1..=n {
                let w = &sin[k] * Q::new((-hand * k as i64).into(), n.into());
                let old = exp[n - k].clone();
                exp[n][0] -= &w * &old[1];
                exp[n][1] += &w * &old[0];
            }
        }
        let value = value(
            circle(),
            [I::point(-0.25), I::point(0.25)],
            I::point(1.),
            I::point(hand as f64),
            I::point(0.5),
        )
        .unwrap();
        for n in 0..18 {
            for (axis, actual) in value.iter().enumerate() {
                let exact = (0..=n).fold(zero.clone(), |s, k| s + &radius[k] * &exp[n - k][axis]);
                let interval = actual.c[n];
                assert!(
                    Q::from_float(interval.lo).unwrap() <= exact
                        && exact <= Q::from_float(interval.hi).unwrap(),
                    "hand{hand} n{n} axis{axis}: {interval:?}"
                );
            }
        }
    }
}

#[test]
fn finite_high_frequency_circle_taylor_norm_overflow_is_an_honest_refusal() {
    // Numerical helper request, not an admitted Thread/Body/material. Individual
    // exact-circle Taylor intervals fit binary64, but squaring their norm does
    // not. The final finite enclosure guard must refuse this whole request.
    let angles = [I::point(-0.25), I::point(0.25)];
    let rate = I::point(1e-11);
    let hand = I::point(1.);
    let whole = value(circle(), angles, rate, hand, I::new(0., 1.)).unwrap();
    assert!(
        whole
            .iter()
            .flat_map(|x| x.c)
            .all(|x| x.lo.is_finite() && x.hi.is_finite())
    );
    assert!(whole.iter().any(|x| x.c[17].magnitude() > 1e154));
    assert!(matches!(
        bounds(circle(), angles, rate, hand, &mut work()),
        Err(ReadError::Unsupported)
    ));
}
