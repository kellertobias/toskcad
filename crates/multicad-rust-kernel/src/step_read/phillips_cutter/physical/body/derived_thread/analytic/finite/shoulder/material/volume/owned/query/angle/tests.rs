use super::*;
use num_bigint::BigInt as Z;
use num_rational::BigRational as Q;
fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn atan(r: Q, n: u32) -> [Q; 2] {
    let mut sum = Q::from_integer(Z::from(0));
    for k in 0..n {
        let term = r.pow((2 * k + 1) as i32) / Z::from(2 * k + 1);
        if k % 2 == 0 {
            sum += term;
        } else {
            sum -= term;
        }
    }
    let extra = r.pow((2 * n + 1) as i32) / Z::from(2 * n + 1);
    if n.is_multiple_of(2) {
        [sum.clone(), sum + extra]
    } else {
        [sum.clone() - extra, sum]
    }
}
fn pi() -> [Q; 2] {
    let a = atan(Q::new(Z::from(1), Z::from(5)), 32);
    let b = atan(Q::new(Z::from(1), Z::from(239)), 32);
    [
        &a[0] * Z::from(16) - &b[1] * Z::from(4),
        &a[1] * Z::from(16) - &b[0] * Z::from(4),
    ]
}
fn ledger() -> Work<'static> {
    Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    }
}
#[test]
fn pythagorean_angle_matches_independent_rational_addition_identity_in_quadrants() {
    let pi = pi();
    let a = atan(Q::new(Z::from(1), Z::from(7)), 40);
    let base = [&pi[0] / Z::from(4) + &a[0], &pi[1] / Z::from(4) + &a[1]];
    for (p, turn) in [
        ([3., 4.], 0),
        ([-4., 3.], 1),
        ([-3., -4.], 2),
        ([4., -3.], -1),
    ] {
        let out = enclose(p.map(I::point), &mut ledger()).unwrap().unwrap().0;
        let shift = if turn >= 0 {
            [
                &pi[0] * Z::from(turn) / Z::from(2),
                &pi[1] * Z::from(turn) / Z::from(2),
            ]
        } else {
            [-&pi[1] / Z::from(2), -&pi[0] / Z::from(2)]
        };
        assert!(
            q(out.lo) <= &base[0] + &shift[0] && q(out.hi) >= &base[1] + &shift[1],
            "point{p:?} {out:?}"
        );
    }
}
#[test]
fn whole_rectangle_extrema_include_axes_and_do_not_cross_an_unearned_branch_cut() {
    let pi = pi();
    for (p, expected) in [
        (
            [I::new(-0.25, 0.25), I::new(3., 4.)],
            [&pi[0] / Z::from(2), &pi[1] / Z::from(2)],
        ),
        ([I::new(-4., -3.), I::new(-0.25, 0.25)], pi.clone()),
        (
            [I::new(-0.25, 0.25), I::new(-4., -3.)],
            [
                &pi[0] * Z::from(3) / Z::from(2),
                &pi[1] * Z::from(3) / Z::from(2),
            ],
        ),
    ] {
        let out = enclose(p, &mut ledger()).unwrap().unwrap().0;
        assert!(q(out.lo) <= expected[0] && q(out.hi) >= expected[1]);
    }
    let out = enclose([I::point(3.), I::point(0.)], &mut ledger())
        .unwrap()
        .unwrap()
        .0;
    assert!(out.lo <= 0. && out.hi >= 0.);
    assert!(
        enclose([I::new(-1., 1.); 2], &mut ledger())
            .unwrap()
            .is_none()
    );
    assert!(
        enclose([I::point(1e-310), I::point(2e-310)], &mut ledger())
            .unwrap()
            .is_none()
    );
    for p in [
        [I::point(f64::NAN), I::point(1.)],
        [I::new(2., 1.), I::point(1.)],
    ] {
        assert!(matches!(
            enclose(p, &mut ledger()),
            Err(PhillipsReadError::Invalid)
        ));
    }
    for limit in [0, 256, 4352, 8448] {
        let mut w = ledger();
        w.remaining = limit;
        assert!(matches!(
            enclose([I::point(3.), I::point(4.)], &mut w),
            Err(PhillipsReadError::Limit)
        ));
    }
}
