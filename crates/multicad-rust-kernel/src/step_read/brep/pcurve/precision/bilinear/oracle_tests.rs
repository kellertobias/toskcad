//! Independent exact coefficient oracles. These test numerical polynomial
//! contracts, never manufacture a Source, material certificate or pointset PASS.
use super::*;
use num_bigint::BigInt;
use num_rational::BigRational as Q;
fn factorial_integer(n: usize) -> BigInt {
    (1..=n).fold(BigInt::from(1), |a, b| a * BigInt::from(b))
}
fn binomial(n: usize, k: usize) -> BigInt {
    factorial_integer(n) / (factorial_integer(k) * factorial_integer(n - k))
}
fn contains(interval: I, exact: &Q) {
    assert!(interval.valid());
    assert!(Q::from_float(interval.lo).unwrap() <= *exact);
    assert!(*exact <= Q::from_float(interval.hi).unwrap());
}
fn polynomial(degree: usize, phase: usize) -> (Poly, Vec<Q>) {
    let mut value = Poly::zero(degree).unwrap();
    let mut exact = Vec::new();
    for (i, coefficient) in value.c.iter_mut().enumerate().take(degree + 1) {
        let numerator = ((i + phase) % 11) as i64 - 5;
        let q = Q::new(numerator.into(), 16.into());
        *coefficient = I::point(numerator as f64 / 16.);
        exact.push(q);
    }
    (value, exact)
}
#[test]
fn products_through_degree_sixty_four_enclose_every_exact_integer_binomial_coefficient() {
    // 43+21 is the actual highest bilinear-degree11/target-degree11 tangent
    // numerator dot product. Binomial(64,32) exceeds exact binary64 integers.
    assert!(binomial(64, 32) > BigInt::from(1_u64 << 53));
    for (m, n) in [(0, 11), (1, 1), (4, 11), (22, 11), (43, 21)] {
        let (a, qa) = polynomial(m, 0);
        let (b, qb) = polynomial(n, 3);
        let actual = product(a, b, &mut Budget::new(512_000_000, None)).unwrap();
        assert_eq!(actual.degree, m + n);
        for k in 0..=m + n {
            let mut sum = Q::from_integer(0.into());
            for (i, qa) in qa.iter().enumerate() {
                if k < i || k - i > n {
                    continue;
                }
                let j = k - i;
                let factor = Q::new(binomial(m, i) * binomial(n, j), binomial(m + n, k));
                sum += factor * qa * &qb[j];
            }
            contains(actual.c[k], &sum);
        }
    }
}
#[test]
fn half_splits_enclose_independent_closed_form_subdivision_coefficients() {
    for degree in [1, 11, 32, 64] {
        let (original, q) = polynomial(degree, 2);
        let [left, right] = original.split(&mut Budget::new(512_000_000, None)).unwrap();
        for i in 0..=degree {
            // Restricting Bernstein form to t=s/2 gives this lower-triangular
            // exact matrix. The right-half matrix is its reversed counterpart;
            // neither oracle runs the production de Casteljau recurrence.
            let mut a = Q::from_integer(0.into());
            for (j, value) in q.iter().enumerate().take(i + 1) {
                a += Q::new(binomial(i, j), BigInt::from(2).pow(i as u32)) * value;
            }
            let mut b = Q::from_integer(0.into());
            for (j, value) in q.iter().enumerate().skip(i) {
                b += Q::new(
                    binomial(degree - i, j - i),
                    BigInt::from(2).pow((degree - i) as u32),
                ) * value;
            }
            contains(left.c[i], &a);
            contains(right.c[i], &b);
        }
    }
}
#[test]
fn complete_subdivision_bound_has_exact_work_frontier_late_cancel_and_retry() {
    let mut x = Poly::zero(2).unwrap();
    x.c[1] = I::point(1.);
    let zero = Poly::zero(2).unwrap();
    let input = [x, zero, zero];
    // Exact scalar polynomial x(t)=2t(1-t) has maximum1/2, while the initial
    // Bernstein coefficient hull is1. At least one real split is necessary.
    let allowed = 0.5001;
    let run = |b: &mut Budget<'_>| bound(input, 1., allowed, b);
    let mut full = Budget::new(512_000_000, None);
    let first = run(&mut full).unwrap();
    assert!((0.5..=allowed).contains(&first));
    let used = full.used();
    assert!(used > 1000);
    let mut exact = Budget::new(used, None);
    assert_eq!(run(&mut exact).unwrap(), first);
    assert_eq!(exact.remaining, 0);
    assert!(matches!(
        run(&mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    let calls = std::cell::Cell::new(0);
    let stop = || {
        calls.set(calls.get() + 1);
        calls.get() == 5
    };
    let mut interrupted = Budget::new(512_000_000, Some(&stop));
    assert!(matches!(run(&mut interrupted), Err(ReadError::Limit)));
    assert!(interrupted.cancelled.get());
    assert_eq!(calls.get(), 5);
    let mut retry = Budget::new(512_000_000, None);
    assert_eq!(run(&mut retry).unwrap(), first);
    assert_eq!(retry.used(), used);
    // The exact midpoint exceeds this allowance. Real depth/work refusal is
    // acceptable; a positive upper certificate is mathematically impossible.
    assert!(matches!(
        bound(input, 1., 0.49, &mut Budget::new(512_000_000, None)),
        Err(ReadError::Unsupported) | Err(ReadError::Limit)
    ));
}
