use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use crate::step::hanger::encoded::cap::tests::{bezier, number};
use num_traits::Zero;
type Q = BigRational;
fn jet(c: Vec<[Q; 4]>, u: &Q) -> ([Q; 3], [Q; 3]) {
    let degree = Q::from_integer((c.len() - 1).into());
    let d = c
        .windows(2)
        .map(|p| std::array::from_fn(|k| (&p[1][k] - &p[0][k]) * &degree))
        .collect();
    let h = bezier(c, u);
    let dh: [Q; 4] = bezier(d, u);
    (
        std::array::from_fn(|k| &h[k] / &h[3]),
        std::array::from_fn(|k| (&dh[k] * &h[3] - &h[k] * &dh[3]) / (&h[3] * &h[3])),
    )
}
#[test]
fn actual_literal_rational_weights_and_cartesian_poles_bound_position_and_derivative() {
    for left in [false, true] {
        let s = Hanger::new(
            crate::authored_hanger_bolt::tests::fixture(left, true),
            None,
        )
        .unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        for band in [ThreadBandRole::Root, ThreadBandRole::Upper] {
            let index = f
                .cells()
                .iter()
                .enumerate()
                .find_map(|(i, c)| {
                    (matches!(c.source_key.face,FaceKey::MetricThread{band:b,..} if b==band)
                        && c.source_key.total_period == -1)
                        .then_some(i)
                })
                .unwrap();
            let curve = (0..f.angular().coefficients().len())
                .find_map(|p| {
                    f.cap_curve_piece(index, End::Start, p, DomainOptions::default(), None)
                        .unwrap()
                })
                .unwrap();
            let literal = curve.literal(None, CurveOptions::default(), None).unwrap();
            assert!(std::ptr::eq(literal.source(), &curve));
            assert_eq!(literal.poles().len(), curve.controls().len());
            assert!(literal.retained_bytes() > curve.retained_bytes);
            let source = curve
                .controls()
                .iter()
                .map(|p| p.each_ref().map(number))
                .collect::<Vec<_>>();
            let literal_h = literal
                .poles()
                .iter()
                .zip(literal.weights())
                .map(|(p, w)| {
                    let w = number(w);
                    let p = p.each_ref().map(number);
                    [&p[0] * &w, &p[1] * &w, &p[2] * &w, w]
                })
                .collect::<Vec<_>>();
            let ratio = &literal_h[0][3] / &source[0][3];
            for (a, b) in source.iter().zip(&literal_h) {
                assert_eq!(&b[3] / &a[3], ratio);
                assert!(b[3] > Q::zero() && b[3] < Q::one());
            }
            for pole in literal.poles() {
                for x in pole {
                    let den = number(x).denom().magnitude().clone();
                    assert_eq!(&den & (&den - 1_u8), num_bigint::BigUint::zero());
                }
            }
            for k in 0..=16 {
                let u = Q::new(k.into(), 16.into());
                let (a, da) = jet(source.clone(), &u);
                let (b, db) = jet(literal_h.clone(), &u);
                let error = (0..3).map(|i| (&a[i] - &b[i]).abs()).sum::<Q>();
                let derivative = (0..3).map(|i| (&da[i] - &db[i]).abs()).sum::<Q>();
                assert!(error <= Q::from_float(literal.position_error_mm()).unwrap());
                assert!(derivative <= Q::from_float(literal.derivative_error_mm()).unwrap());
            }
            eprintln!(
                "literal cap band={band:?} left={left} work={} position={} derivative={} bytes={}",
                literal.work_used(),
                literal.position_error_mm(),
                literal.derivative_error_mm(),
                literal.retained_bytes()
            );
        }
    }
}
#[test]
fn exact_numeric_quantizer_weights_and_shared_endpoint_overrides_are_measured() {
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut m = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let half = Q::new(BigInt::one(), BigInt::one() << (FRACTION_BITS + 1));
    for x in [
        Q::new(1.into(), 3.into()),
        Q::new((-1).into(), 3.into()),
        half.clone(),
        -half.clone(),
    ] {
        let chosen = quantize(&N(x.clone()), &mut w).unwrap();
        assert!((&chosen.0 - x).abs() <= half);
    }
    assert!(matches!(
        finite_weights(
            &[[N::integer(0), N::integer(0), N::integer(0), N::integer(0)]],
            &mut m,
            &mut w
        ),
        Err(Error::Unresolved)
    ));
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let curve = (0..f.cells().len())
        .find_map(|i| {
            f.cap_curve_piece(i, End::Start, 3, DomainOptions::default(), None)
                .unwrap()
        })
        .unwrap();
    let endpoint = |i: usize, w: &mut Work<'_>, m: &mut Memory, shift: N| {
        let row = &curve.controls()[i];
        let weight = row[3].number(w).unwrap();
        std::array::from_fn(|axis| {
            let x = row[axis].number(w).unwrap().div(&weight, w).unwrap();
            let x = if axis == 0 {
                x.add(&shift, w).unwrap()
            } else {
                x
            };
            EncodedScalar::retain(&x, m, w).unwrap()
        })
    };
    let tiny = N(Q::new(BigInt::one(), BigInt::one() << 80));
    let endpoints = [
        endpoint(0, &mut w, &mut m, tiny),
        endpoint(curve.degree(), &mut w, &mut m, N::integer(0)),
    ];
    let literal = curve
        .literal(Some(&endpoints), CurveOptions::default(), None)
        .unwrap();
    for (end, point) in endpoints.iter().enumerate() {
        let pole = &literal.poles()[if end == 0 { 0 } else { curve.degree() }];
        for (actual, expected) in pole.iter().zip(point) {
            assert_eq!(
                number(actual),
                quantize(&N(number(expected)), &mut w).unwrap().0
            );
        }
    }
    assert!(literal.position_error_mm() > 0.);
    let wrong = [
        endpoint(0, &mut w, &mut m, N::integer(1)),
        endpoint(curve.degree(), &mut w, &mut m, N::integer(0)),
    ];
    assert!(matches!(
        curve.literal(Some(&wrong), CurveOptions::default(), None),
        Err(Error::Unresolved)
    ));
}
#[test]
fn literal_errors_cancellation_precision_and_actual_capacity_limits_are_atomic() {
    use std::cell::Cell;
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let curve = (0..f.cells().len())
        .find_map(|i| {
            f.cap_curve_piece(i, End::Start, 3, DomainOptions::default(), None)
                .unwrap()
        })
        .unwrap();
    for (options, expected) in [
        (
            CurveOptions {
                position_error_mm: 0.,
                ..CurveOptions::default()
            },
            Error::Invalid,
        ),
        (
            CurveOptions {
                position_error_mm: f64::from_bits(1),
                ..CurveOptions::default()
            },
            Error::Unresolved,
        ),
        (
            CurveOptions {
                work: 1,
                ..CurveOptions::default()
            },
            Error::Limit,
        ),
        (
            CurveOptions {
                bytes: 1,
                ..CurveOptions::default()
            },
            Error::Limit,
        ),
    ] {
        assert!(matches!(curve.literal(None,options,None),Err(e) if e==expected));
    }
    let observed = Cell::new(0usize);
    let count = || {
        observed.set(observed.get() + 1);
        false
    };
    let literal = curve
        .literal(None, CurveOptions::default(), Some(&count))
        .unwrap();
    for stop in (1..=32).chain([observed.get() / 2, observed.get()]) {
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        assert!(matches!(
            curve.literal(None, CurveOptions::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    assert!(matches!(
        curve.literal(
            None,
            CurveOptions {
                work: literal.work_used() - 1,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(s.spec(), b.source().spec());
}
#[test]
fn genuine_subnormal_plain_station_is_retained_and_off_plane_endpoint_is_refused() {
    let mut spec = crate::authored_hanger_bolt::tests::fixture(false, false);
    spec.plain_length_mm = f64::from_bits(1);
    let s = Hanger::new(spec, None).unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let index = f
        .cells()
        .iter()
        .position(|c| {
            c.source_key.total_period == 23
                && matches!(
                    c.source_key.face,
                    FaceKey::MetricThread {
                        band: ThreadBandRole::Root,
                        ..
                    }
                )
        })
        .unwrap();
    let curve = f
        .cap_curve_piece(index, End::End, 3, DomainOptions::default(), None)
        .unwrap()
        .unwrap();
    let literal = curve.literal(None, CurveOptions::default(), None).unwrap();
    let plane = Q::from_float(-f64::from_bits(1)).unwrap();
    for p in literal.poles() {
        assert_eq!(number(&p[2]), plane);
    }
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let mut endpoints = [
        literal.poles()[0].clone(),
        literal.poles()[curve.degree()].clone(),
    ];
    endpoints[0][2] = EncodedScalar::retain(&N::integer(0), &mut memory, &mut w).unwrap();
    assert!(matches!(
        curve.literal(Some(&endpoints), CurveOptions::default(), None),
        Err(Error::Unresolved)
    ));
    assert_eq!(s.spec().plain_length_mm, f64::from_bits(1));
}

#[test]
fn numeric_encoder_checkpoints_cover_lcm_and_quantizer_without_partial_receipts() {
    use std::cell::Cell;
    let calls = Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    // Different positive rational weights require genuine LCM/GCD work;
    // these are numeric primitive inputs, never a fabricated source body.
    let h = (0..=11)
        .map(|i| {
            let weight = N(Q::new(1.into(), (3 + 2 * i).into()));
            [weight.clone(), weight.clone(), N::integer(0), weight]
        })
        .collect::<Vec<_>>();
    let mut w = Work::new(16_000_000, Some(&observe)).unwrap();
    let mut m = Memory::new(0, 64 * 1024 * 1024).unwrap();
    encode_cap(&h, None, CurveOptions::default(), &mut m, &mut w).unwrap();
    for stop in 1..=calls.get() {
        let count = Cell::new(0);
        let cancel = || {
            count.set(count.get() + 1);
            count.get() >= stop
        };
        let mut w = Work::new(16_000_000, Some(&cancel)).unwrap();
        let mut m = Memory::new(0, 64 * 1024 * 1024).unwrap();
        assert!(
            matches!(
                encode_cap(&h, None, CurveOptions::default(), &mut m, &mut w),
                Err(Error::Cancelled)
            ),
            "stop={stop}"
        );
    }
}
