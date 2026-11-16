use super::*;
#[test]
fn actual_closed_unit_reduction_retains_cardinals_original_parameter_and_shared_bounds() {
    for left in [false, true] {
        for position in [1e-7, 1e-10] {
            let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
            let mut w = Work::new(16_000_000, None).unwrap();
            let r = fit_reduced_unit(
                left,
                CurveOptions {
                    position_error_mm: position,
                    ..CurveOptions::default()
                },
                &mut memory,
                &mut w,
            )
            .unwrap();
            assert_eq!(r.left_hand(), left);
            assert_eq!(r.knots()[0], 0.);
            assert_eq!(*r.knots().last().unwrap(), 1.);
            assert!(r.position_error() <= position && r.derivative_error() <= 1e-6);
            for pair in r.coefficients().windows(2) {
                assert_eq!(pair[0][10], pair[1][0]);
            }
            assert_eq!(
                r.coefficients().first().unwrap()[0],
                r.coefficients().last().unwrap()[10]
            );
            assert_eq!(
                r.coefficients().first().unwrap()[0],
                [BigRational::one(), BigRational::zero()]
            );
            eprintln!(
                "reduced unit {left}/{position:e}: {}pieces E0/E1 {:?}, {}work {}bytes",
                r.coefficients().len(),
                [r.position_error(), r.derivative_error()],
                w.used,
                memory.used
            );
            r.release_memory(&mut memory);
            drop(r);
            assert_eq!(memory.used, 0);
        }
    }
}
#[test]
fn actual_unit_resource_precision_cancel_and_retry_do_not_publish_partial_receipts() {
    use std::cell::Cell;
    let count = Cell::new(0);
    let poll = || {
        count.set(count.get() + 1);
        false
    };
    let mut w = Work::new(16_000_000, Some(&poll)).unwrap();
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let baseline = fit_reduced_unit(false, CurveOptions::default(), &mut memory, &mut w).unwrap();
    let work = w.used;
    let callbacks = count.get();
    baseline.release_memory(&mut memory);
    drop(baseline);
    for at in [0, callbacks / 2, callbacks - 1] {
        let calls = Cell::new(0);
        let stop = || {
            let n = calls.get();
            calls.set(n + 1);
            n >= at
        };
        let mut w = Work::new(16_000_000, Some(&stop)).unwrap();
        let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
        assert!(matches!(
            fit_reduced_unit(false, CurveOptions::default(), &mut memory, &mut w),
            Err(Error::Cancelled)
        ));
    }
    let mut w = Work::new(work - 1, None).unwrap();
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    assert!(matches!(
        fit_reduced_unit(false, CurveOptions::default(), &mut memory, &mut w),
        Err(Error::Limit)
    ));
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut memory = Memory::new(0, 1).unwrap();
    assert!(matches!(
        fit_reduced_unit(false, CurveOptions::default(), &mut memory, &mut w),
        Err(Error::Limit)
    ));
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    assert!(matches!(
        fit_reduced_unit(
            false,
            CurveOptions {
                position_error_mm: f64::from_bits(1),
                ..CurveOptions::default()
            },
            &mut memory,
            &mut w
        ),
        Err(Error::Unresolved)
    ));
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    assert!(matches!(
        fit_reduced_unit(
            false,
            CurveOptions {
                position_error_mm: 0.,
                ..CurveOptions::default()
            },
            &mut memory,
            &mut w
        ),
        Err(Error::Invalid)
    ));
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    fit_reduced_unit(false, CurveOptions::default(), &mut memory, &mut w).unwrap();
}

#[test]
fn public_unit_query_uses_its_own_real_limits_and_same_immutable_coefficients() {
    let q = ReducedUnitAngular::prepare(false, CurveOptions::default(), None).unwrap();
    assert!(!q.left_hand());
    assert!(q.position_error() <= 1e-7);
    assert!(matches!(
        ReducedUnitAngular::prepare(
            false,
            CurveOptions {
                work: 1,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        ReducedUnitAngular::prepare(
            false,
            CurveOptions {
                work: 0,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
}
fn atan_reciprocal(k: i64) -> [BigRational; 2] {
    let x = BigRational::new(1.into(), k.into());
    let square = &x * &x;
    let mut power = x;
    let mut sum = BigRational::zero();
    for n in 0..80 {
        let term = &power / BigInt::from(2 * n + 1);
        if n % 2 == 0 {
            sum += term;
        } else {
            sum -= term;
        }
        power *= &square;
    }
    // Eighty alternating terms end below atan(x); next term bounds the tail.
    let upper = &sum + power / BigInt::from(161);
    [sum, upper]
}
#[test]
fn independent_exact_rational_eighth_turn_value_and_original_derivative_oracle() {
    let a = atan_reciprocal(5);
    let b = atan_reciprocal(239);
    // Machin's identity follows by tan(4atan(1/5)-atan(1/239))=1 and its
    // principal-angle range. This oracle uses no production trig/evaluation.
    let pi = [
        &a[0] * BigInt::from(16) - &b[1] * BigInt::from(4),
        &a[1] * BigInt::from(16) - &b[0] * BigInt::from(4),
    ];
    let denominator = BigInt::from(10).pow(60);
    let digits = BigInt::parse_bytes(
        b"707106781186547524400844362104849039284835937688474036588339",
        10,
    )
    .unwrap();
    let h = [
        BigRational::new(digits.clone(), denominator.clone()),
        BigRational::new(digits + 1, denominator),
    ];
    assert!(&h[0] * &h[0] * BigInt::from(2) < BigRational::one());
    assert!(&h[1] * &h[1] * BigInt::from(2) > BigRational::one());
    for left in [false, true] {
        let unit = ReducedUnitAngular::prepare(left, CurveOptions::default(), None).unwrap();
        let index = unit.knots().iter().position(|t| *t == 0.125).unwrap();
        let p = &unit.coefficients()[index];
        let span = BigRational::from_float(unit.knots()[index + 1] - unit.knots()[index]).unwrap();
        let mut position_error = BigRational::zero();
        let mut derivative_error = BigRational::zero();
        for (axis, (position, next)) in p[0].iter().zip(&p[1]).enumerate() {
            let negative = axis == 1 && left;
            let value = if negative {
                [-h[1].clone(), -h[0].clone()]
            } else {
                h.clone()
            };
            position_error += (position - &value[0])
                .abs()
                .max((position - &value[1]).abs());
            let magnitude = [
                &pi[0] * &h[0] * BigInt::from(2),
                &pi[1] * &h[1] * BigInt::from(2),
            ];
            let negative = axis == 0 || left;
            let derivative = if negative {
                [-magnitude[1].clone(), -magnitude[0].clone()]
            } else {
                magnitude
            };
            let actual = (next - position) * BigInt::from(10) / &span;
            derivative_error += (&actual - &derivative[0])
                .abs()
                .max((&actual - &derivative[1]).abs());
        }
        assert!(position_error <= BigRational::from_float(unit.position_error()).unwrap());
        assert!(derivative_error <= BigRational::from_float(unit.derivative_error()).unwrap());
    }
}
