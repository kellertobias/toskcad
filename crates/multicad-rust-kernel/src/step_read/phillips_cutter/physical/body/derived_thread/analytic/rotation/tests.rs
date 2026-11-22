use super::*;
use num_rational::BigRational as Q;
use num_traits::One;
fn atan_q(x: Q) -> [Q; 2] {
    let square = &x * &x;
    let mut power = x.clone();
    let mut sum = x;
    for k in 1..=31 {
        power *= &square;
        let term = &power / Q::from_integer(Z::from(2 * k + 1));
        if k % 2 == 0 {
            sum += term;
        } else {
            sum -= term;
        }
    }
    let remainder = power * square / Q::from_integer(Z::from(65));
    [sum.clone(), sum + remainder]
}
fn interval(low: Q, high: Q) -> I {
    let mut a = low.to_f64().unwrap();
    let mut b = high.to_f64().unwrap();
    if Q::from_float(a).unwrap() > low {
        a = a.next_down();
    }
    if Q::from_float(b).unwrap() < high {
        b = b.next_up();
    }
    assert!(Q::from_float(a).unwrap() <= low && Q::from_float(b).unwrap() >= high);
    I::new(a, b)
}
#[test]
fn exact_rational_long_series_independently_checks_reduced_trig_enclosures() {
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let trig = Rotation::new(&mut work).unwrap();
    for (n, d) in [
        (0, 1),
        (1, 3),
        (2, 1),
        (3, 1),
        (7, 1),
        (-7, 1),
        (16, 1),
        (-16, 1),
    ] {
        let x = Q::new(Z::from(n), Z::from(d));
        let square = &x * &x;
        let mut cosine = Q::one();
        let mut sine = x.clone();
        let mut cterm = Q::one();
        let mut sterm = x.clone();
        for k in 1..=50 {
            cterm *= -&square / Q::from_integer(Z::from((2 * k - 1) * (2 * k)));
            sterm *= -&square / Q::from_integer(Z::from((2 * k) * (2 * k + 1)));
            cosine += &cterm;
            sine += &sterm;
        }
        let remainder = Q::new(Z::from(1), Z::from(10).pow(30));
        // |x|<=16 makes the omitted degree102/103 terms below1e-30.
        let actual = trig.sin_cos(interval(x.clone(), x), &mut work).unwrap();
        for (bound, oracle) in actual.into_iter().zip([cosine, sine]) {
            if n == 0 {
                assert!(
                    Q::from_float(bound.lo).unwrap() <= oracle
                        && Q::from_float(bound.hi).unwrap() >= oracle
                );
            } else {
                assert!(Q::from_float(bound.lo).unwrap() <= &oracle - &remainder);
                assert!(Q::from_float(bound.hi).unwrap() >= oracle + &remainder);
            }
        }
    }
    for turns in [-500, 500] {
        let six = Q::from_integer(Z::from(6));
        let t = Q::from_integer(Z::from(2 * turns));
        let a = &Q::new(trig.bounds[0].0.clone(), trig.bounds[0].1.clone())
            * (t.clone() + Q::one() / &six);
        let b = &Q::new(trig.bounds[1].0.clone(), trig.bounds[1].1.clone()) * (t + Q::one() / six);
        let value = trig
            .sin_cos(interval(a.clone().min(b.clone()), a.max(b)), &mut work)
            .unwrap();
        assert!(value[1].lo <= 0.5 && value[1].hi >= 0.5);
    }
}

#[test]
fn fixed_pi_constants_equal_exact_alternating_machin_bounds_and_frontiers_refuse() {
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let trig = Rotation::new(&mut work).unwrap();
    let [a, b] = atan_q(Q::new(Z::one(), Z::from(5)));
    let [c, d] = atan_q(Q::new(Z::one(), Z::from(239)));
    assert_eq!(
        Q::new(trig.bounds[0].0.clone(), trig.bounds[0].1.clone()),
        Q::from_integer(Z::from(16)) * a - Q::from_integer(Z::from(4)) * d
    );
    assert_eq!(
        Q::new(trig.bounds[1].0.clone(), trig.bounds[1].1.clone()),
        Q::from_integer(Z::from(16)) * b - Q::from_integer(Z::from(4)) * c
    );
    for x in [-1e7, -5000., 5000., 1e7, f64::from_bits(1)] {
        assert!(trig.sin_cos(I::point(x), &mut work).is_ok());
    }
    for x in [I::new(-10., 10.), I::point(f64::INFINITY), I::point(1e8)] {
        assert!(trig.sin_cos(x, &mut work).is_err());
    }
    work.remaining = 0;
    assert!(trig.sin_cos(I::point(1.), &mut work).is_err());
    assert!(Rotation::new(&mut work).is_err());
}
