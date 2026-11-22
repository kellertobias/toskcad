use super::*;
fn q(n: i64, d: i64) -> BigRational {
    BigRational::new(n.into(), d.into())
}
fn factorial(n: usize) -> BigInt {
    (1..=n).fold(BigInt::from(1), |v, i| v * BigInt::from(i))
}
fn binomial(n: usize, k: usize) -> BigRational {
    BigRational::new(factorial(n), factorial(k) * factorial(n - k))
}
fn contains(interval: I, value: &BigRational) {
    assert!(interval.valid());
    assert!(BigRational::from_float(interval.lo).unwrap() <= *value);
    assert!(BigRational::from_float(interval.hi).unwrap() >= *value);
}
#[test]
fn every_degree_sixteen_product_weight_matches_exact_factorials() {
    let mut b = Budget::new(512_000_000, None);
    for i in 0..=8 {
        for j in 0..=8 {
            let mut a = [I::point(0.); 9];
            a[i] = I::point(1.);
            let mut c = [I::point(0.); 9];
            c[j] = I::point(1.);
            let a = Polynomial::new(&a, &mut b).unwrap();
            let c = Polynomial::new(&c, &mut b).unwrap();
            let product = a.multiply(&c, &mut b).unwrap();
            assert_eq!(product.degree(), 16);
            for (k, value) in product.controls().iter().enumerate() {
                let expected = if k == i + j {
                    binomial(8, i) * binomial(8, j) / binomial(16, k)
                } else {
                    q(0, 1)
                };
                contains(*value, &expected);
            }
        }
    }
}
#[test]
fn exact_split_elevation_derivative_and_resource_frontiers() {
    fn operation(b: &mut Budget<'_>) -> Result<Polynomial> {
        let p = Polynomial::new(&[I::point(1.), I::point(2.), I::point(4.), I::point(8.)], b)?;
        let elevated = p.elevate(12, b)?;
        let [left, right] = elevated.split(I::point(0.25), b)?;
        let d = p.derivative(b)?;
        let sum = left.add(&right, b)?.scale(I::point(0.5), b)?;
        check(d.hull().lo > 0., ReadError::Unsupported)?;
        Ok(sum)
    }
    let mut b = Budget::new(512_000_000, None);
    let p = Polynomial::new(&[I::point(0.), I::point(1.)], &mut b).unwrap();
    let elevated = p.elevate(16, &mut b).unwrap();
    for (i, c) in elevated.controls().iter().enumerate() {
        contains(*c, &q(i as i64, 16));
    }
    let [left, right] = p.split(I::point(0.25), &mut b).unwrap();
    contains(left.controls()[1], &q(1, 4));
    contains(right.controls()[0], &q(1, 4));
    let d = p.derivative(&mut b).unwrap();
    contains(d.controls()[0], &q(1, 1));
    let mut b = Budget::new(512_000_000, None);
    operation(&mut b).unwrap();
    let used = b.used();
    assert!(matches!(
        operation(&mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    operation(&mut Budget::new(used, None)).unwrap();
    let calls = std::cell::Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    operation(&mut Budget::new(512_000_000, Some(&observe))).unwrap();
    let last = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == last
    };
    assert!(matches!(
        operation(&mut Budget::new(512_000_000, Some(&cancel))),
        Err(ReadError::Limit)
    ));
    operation(&mut Budget::new(512_000_000, None)).unwrap();
}
