//! Standalone numerical contracts; these objects are not Source certificates.
use super::*;
use num_bigint::BigInt;
use num_rational::BigRational;

/// Independent Pascal addition, rather than production's multiplicative u128
/// recurrence or any floating binomial coefficient.
fn pascal(n: usize) -> Vec<Vec<BigInt>> {
    let mut rows: Vec<Vec<BigInt>> = Vec::new();
    for degree in 0..=n {
        let mut row = vec![BigInt::from(1); degree + 1];
        for i in 1..degree {
            row[i] = &rows[degree - 1][i - 1] + &rows[degree - 1][i];
        }
        rows.push(row);
    }
    rows
}
fn integer_control(i: usize, offset: usize) -> i64 {
    ((17 * i + offset) % 31) as i64 - 15
}
fn polynomial(degree: usize, denominator: i64, offset: usize) -> Poly {
    let mut p = Poly::zero(degree).unwrap();
    for i in 0..=degree {
        p.c[i] = I::point(integer_control(i, offset) as f64 / denominator as f64);
    }
    p
}
fn encloses(i: I, expected: &BigRational) {
    assert!(i.valid(), "invalid enclosure: {i:?}");
    let lo = BigRational::from_float(i.lo).unwrap();
    let hi = BigRational::from_float(i.hi).unwrap();
    assert!(
        lo <= *expected && *expected <= hi,
        "{i:?} misses {expected}"
    );
}
fn exact_product(degrees: [usize; 2]) {
    let [m, n] = degrees;
    let a = polynomial(m, 64, 0);
    let b = polynomial(n, 128, 9);
    let mut work = Budget::new(512_000_000, None);
    let p = product(a, b, &mut work).unwrap();
    assert_eq!(p.degree, m + n);
    let choose = pascal(m + n);
    for k in 0..=m + n {
        let mut numerator = BigInt::from(0);
        for i in k.saturating_sub(n)..=k.min(m) {
            let j = k - i;
            numerator += &choose[m][i]
                * &choose[n][j]
                * BigInt::from(integer_control(i, 0))
                * BigInt::from(integer_control(j, 9));
        }
        let expected = BigRational::new(numerator, &choose[m + n][k] * BigInt::from(8192));
        encloses(p.c[k], &expected);
        assert!(
            p.c[k].width() < 1e-10,
            "lost all useful signed-control precision"
        );
    }
}
#[test]
fn all_degree77_and_degree80_product_coefficients_enclose_independent_exact_oracles() {
    exact_product([37, 40]);
    exact_product([39, 41]);
}
#[test]
fn degree80_split_and_derivatives_enclose_every_exact_dyadic_coefficient() {
    let p = polynomial(80, 64, 3);
    let mut row: Vec<BigRational> = (0..=80)
        .map(|i| BigRational::new(integer_control(i, 3).into(), 64.into()))
        .collect();
    let mut left = vec![BigRational::from_integer(0.into()); 81];
    let mut right = left.clone();
    left[0] = row[0].clone();
    right[80] = row[80].clone();
    for level in 1..=80 {
        for i in 0..=80 - level {
            row[i] = (&row[i] + &row[i + 1]) / BigRational::from_integer(2.into());
        }
        left[level] = row[0].clone();
        right[80 - level] = row[80 - level].clone();
    }
    let mut work = Budget::new(512_000_000, None);
    let [a, b] = p.split(&mut work).unwrap();
    for i in 0..=80 {
        encloses(a.c[i], &left[i]);
        encloses(b.c[i], &right[i]);
    }
    let d = p.derivative(&mut work).unwrap();
    assert_eq!(d.degree, 79);
    for i in 0..80 {
        let expected = BigRational::new(
            BigInt::from(80 * (integer_control(i + 1, 3) - integer_control(i, 3))),
            64.into(),
        );
        encloses(d.c[i], &expected);
    }
}
#[test]
fn degree5_to80_elevation_encloses_exact_pascal_weights_without_new_geometry() {
    let p = polynomial(5, 64, 7);
    let elevated = p.elevate(80, &mut Budget::new(512_000_000, None)).unwrap();
    let choose = pascal(80);
    for i in 0usize..=80 {
        let mut numerator = BigInt::from(0);
        for j in i.saturating_sub(75)..=i.min(5) {
            numerator += &choose[5][j] * &choose[75][i - j] * BigInt::from(integer_control(j, 7));
        }
        encloses(
            elevated.c[i],
            &BigRational::new(numerator, &choose[80][i] * BigInt::from(64)),
        );
    }
}
#[test]
fn actual_numerical_work_frontiers_observed_cancellation_and_unchanged_retry() {
    let a = polynomial(39, 64, 0);
    let b = polynomial(41, 128, 9);
    let mut full = Budget::new(512_000_000, None);
    let expected = product(a, b, &mut full).unwrap();
    let used = full.used();
    assert!(matches!(
        product(a, b, &mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    let mut exact = Budget::new(used, None);
    let retry = product(a, b, &mut exact).unwrap();
    assert_eq!(exact.remaining, 0);
    for i in 0..=80 {
        assert_eq!(retry.c[i].lo.to_bits(), expected.c[i].lo.to_bits());
        assert_eq!(retry.c[i].hi.to_bits(), expected.c[i].hi.to_bits());
    }
    let polls = std::cell::Cell::new(0);
    let observe = || {
        polls.set(polls.get() + 1);
        false
    };
    product(a, b, &mut Budget::new(512_000_000, Some(&observe))).unwrap();
    let final_poll = polls.get();
    assert!(final_poll >= 4);
    polls.set(0);
    let stop = || {
        polls.set(polls.get() + 1);
        polls.get() == final_poll
    };
    let mut cancelled = Budget::new(512_000_000, Some(&stop));
    assert!(matches!(
        product(a, b, &mut cancelled),
        Err(ReadError::Limit)
    ));
    assert!(cancelled.cancelled.get());
    let retry = product(a, b, &mut Budget::new(512_000_000, None)).unwrap();
    assert_eq!(retry.c[40].lo.to_bits(), expected.c[40].lo.to_bits());
    assert!(matches!(
        product(
            polynomial(40, 64, 0),
            b,
            &mut Budget::new(512_000_000, None)
        ),
        Err(ReadError::Unsupported)
    ));
    assert!(matches!(
        Poly::constant(
            80,
            I::point(f64::INFINITY),
            &mut Budget::new(512_000_000, None)
        ),
        Err(ReadError::Unsupported)
    ));
}
/// Exact `k`th Bernstein control of `polynomial(m, 64, a) * polynomial(n, 128, b)`.
fn exact_coefficient(m: usize, n: usize, [a, b]: [usize; 2], k: usize) -> BigRational {
    let choose = pascal(m + n);
    let mut numerator = BigInt::from(0);
    for i in k.saturating_sub(n)..=k.min(m) {
        numerator += &choose[m][i]
            * &choose[n][k - i]
            * BigInt::from(integer_control(i, a))
            * BigInt::from(integer_control(k - i, b));
    }
    BigRational::new(numerator, &choose[m + n][k] * BigInt::from(8192))
}
#[test]
fn shared_weight_products_and_elevations_enclose_every_factor_exactly() {
    let (m, n) = (39, 41);
    let factors = [0, 5, 11].map(|offset| polynomial(m, 64, offset));
    let b = polynomial(n, 128, 9);
    let mut work = Budget::new(512_000_000, None);
    let products = product_many(factors, b, &mut work).unwrap();
    let used = work.used();
    for (p, offset) in products.iter().zip([0, 5, 11]) {
        assert_eq!(p.degree, m + n);
        for k in 0..=m + n {
            encloses(p.c[k], &exact_coefficient(m, n, [offset, 9], k));
            assert!(p.c[k].width() < 1e-10, "lost signed-control precision");
        }
    }
    // The batch is cheaper than three separate products of the same pairs.
    let mut separate = Budget::new(512_000_000, None);
    for f in factors {
        product(f, b, &mut separate).unwrap();
    }
    assert!(used < separate.used());
    // Its work frontier is exact and a refused attempt keeps no partial result.
    assert!(matches!(
        product_many(factors, b, &mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    let mut exact = Budget::new(used, None);
    product_many(factors, b, &mut exact).unwrap();
    assert_eq!(exact.remaining, 0);

    let low = [polynomial(5, 64, 7), polynomial(5, 64, 2)];
    let elevated = elevate_many(low, 80, &mut Budget::new(512_000_000, None)).unwrap();
    let choose = pascal(80);
    for (p, offset) in elevated.iter().zip([7, 2]) {
        for i in 0usize..=80 {
            let mut numerator = BigInt::from(0);
            for j in i.saturating_sub(75)..=i.min(5) {
                numerator +=
                    &choose[5][j] * &choose[75][i - j] * BigInt::from(integer_control(j, offset));
            }
            encloses(
                p.c[i],
                &BigRational::new(numerator, &choose[80][i] * BigInt::from(64)),
            );
        }
    }
    let same = elevate_many(low, 5, &mut Budget::new(512_000_000, None)).unwrap();
    assert!(same.iter().zip(&low).all(|(a, b)| {
        a.degree == 5
            && a.c
                .iter()
                .zip(&b.c)
                .all(|(x, y)| x.lo == y.lo && x.hi == y.hi)
    }));
    // Mixed degrees, lowering and degrees past the arena are refused.
    let mixed = [polynomial(5, 64, 0), polynomial(6, 64, 0)];
    let mut work = Budget::new(512_000_000, None);
    assert!(matches!(
        product_many(mixed, b, &mut work),
        Err(ReadError::Unsupported)
    ));
    assert!(matches!(
        elevate_many(mixed, 80, &mut work),
        Err(ReadError::Unsupported)
    ));
    assert!(matches!(
        elevate_many(low, 4, &mut work),
        Err(ReadError::Unsupported)
    ));
    assert!(matches!(
        product_many([polynomial(40, 64, 0)], polynomial(41, 64, 0), &mut work),
        Err(ReadError::Unsupported)
    ));
}
