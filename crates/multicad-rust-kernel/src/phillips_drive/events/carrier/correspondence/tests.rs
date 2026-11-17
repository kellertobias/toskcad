use super::*;
use num_rational::BigRational as Q;
use num_traits::{One, ToPrimitive, Zero};
use std::cell::Cell;
fn enclose(q: &Q) -> I {
    let v = q.to_f64().unwrap();
    let exact = Q::from_float(v).unwrap();
    match exact.cmp(q) {
        std::cmp::Ordering::Equal => I::point(v),
        std::cmp::Ordering::Less => I::new(v, v.next_up()),
        std::cmp::Ordering::Greater => I::new(v.next_down(), v),
    }
}
fn choose(n: usize, k: usize) -> usize {
    if k > n {
        0
    } else {
        (0..k).fold(1, |p, i| p * (n - i) / (i + 1))
    }
}
fn power(i: I, n: usize) -> I {
    (0..n).fold(I::point(1.), |p, _| p * i)
}
fn scalar(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn bezier(p: &[[Q; 3]; 4], t: &Q, derivative: bool) -> [Q; 3] {
    let one = Q::one() - t;
    std::array::from_fn(|axis| {
        if derivative {
            (0..3).fold(Q::zero(), |s, i| {
                s + Q::from_integer(3.into())
                    * (&p[i + 1][axis] - &p[i][axis])
                    * Q::from_integer(choose(2, i).into())
                    * t.pow(i as i32)
                    * one.pow((2 - i) as i32)
            })
        } else {
            (0..4).fold(Q::zero(), |s, i| {
                s + &p[i][axis]
                    * Q::from_integer(choose(3, i).into())
                    * t.pow(i as i32)
                    * one.pow((3 - i) as i32)
            })
        }
    })
}
#[test]
fn trim_c1_taylor3_bound_encloses_independent_rational_polynomials_and_chain_factors() {
    let columns = [[0.8, 0.6, 0.], [-0.36, 0.48, 0.8], [0.48, -0.64, 0.6]];
    let origin = [11., -7., 23.];
    let scale = [
        Q::from_integer(2.into()),
        Q::from_integer((-3).into()),
        Q::new(1.into(), 2.into()),
    ];
    let offset = [
        Q::from_integer(1_000_000.into()),
        Q::from_integer((-1_000_000).into()),
        Q::new(7.into(), 10.into()),
    ];
    for degree in 0..=4usize {
        for (a, b) in [
            (Q::zero(), Q::one()),
            (Q::new(1.into(), 10.into()), Q::new(3.into(), 10.into())),
            (Q::new(1.into(), 8.into()), Q::new(3.into(), 4.into())),
        ] {
            let h = &b - &a;
            let ai = enclose(&a);
            let bi = enclose(&b);
            let hi = bi - ai;
            let m = (ai + bi) / I::point(2.);
            let midpoint = std::array::from_fn(|axis| Series::<5> {
                c: std::array::from_fn(|k| {
                    let term = if k <= degree {
                        enclose(&scale[axis])
                            * I::point(choose(degree, k) as f64)
                            * power(m, degree - k)
                    } else {
                        I::point(0.)
                    };
                    if k == 0 {
                        term + enclose(&offset[axis])
                    } else {
                        term
                    }
                }),
            });
            let p = polynomial(midpoint, hi);
            let fourth = std::array::from_fn(|axis| {
                if degree == 4 {
                    enclose(&scale[axis])
                } else {
                    I::point(0.)
                }
            });
            let value = |t: &Q, derivative: bool| {
                std::array::from_fn::<_, 3, _>(|axis| {
                    if derivative {
                        if degree == 0 {
                            Q::zero()
                        } else {
                            &scale[axis]
                                * Q::from_integer(degree.into())
                                * t.pow((degree - 1) as i32)
                                * &h
                        }
                    } else {
                        &offset[axis] + &scale[axis] * t.pow(degree as i32)
                    }
                })
            };
            let first = value(&a, false);
            let last = value(&b, false);
            let df = value(&a, true);
            let dl = value(&b, true);
            let local: [[_; 3]; 4] = [
                first.clone(),
                std::array::from_fn(|i| &first[i] + &df[i] / Q::from_integer(3.into())),
                std::array::from_fn(|i| &last[i] - &dl[i] / Q::from_integer(3.into())),
                last,
            ];
            let mut raw: [[Q; 3]; 4] = std::array::from_fn(|j| {
                std::array::from_fn(|i| {
                    scalar(origin[i])
                        + (0..3).fold(Q::zero(), |s, k| s + scalar(columns[k][i]) * &local[j][k])
                })
            });
            raw[2][1] += Q::new(1.into(), 1_000_000_000.into());
            let controls = raw.clone().map(|row| row.map(|q| enclose(&q)));
            let bounds =
                difference_bounds(p, fourth, hi * hi * hi * hi, controls, origin, columns).unwrap();
            for i in 0..=32 {
                let t = Q::new(i.into(), 32.into());
                let at = &a + &h * &t;
                for derivative in [false, true] {
                    let source = value(&at, derivative);
                    let expected: [Q; 3] = std::array::from_fn(|axis| {
                        (if derivative {
                            Q::zero()
                        } else {
                            scalar(origin[axis])
                        }) + (0..3)
                            .fold(Q::zero(), |s, k| s + scalar(columns[k][axis]) * &source[k])
                    });
                    let actual = bezier(&raw, &t, derivative);
                    let square =
                        (0..3).fold(Q::zero(), |s, k| s + (&expected[k] - &actual[k]).pow(2));
                    assert!(
                        square <= scalar(bounds[usize::from(derivative)]).pow(2),
                        "degree{degree} s{t} derivative{derivative} bounds{bounds:?}"
                    );
                }
            }
        }
    }
    let huge = [[[I::point(f64::MAX); 3]; 4]; 1][0];
    assert!(matches!(
        difference_bounds(
            huge,
            [I::point(0.); 3],
            I::point(1.),
            [[I::point(0.); 3]; 4],
            [0.; 3],
            [[2.; 3]; 3]
        ),
        Err(Error::Unresolved)
    ));
}
#[test]
fn trim_c1_source_binding_exact_limits_cancellation_and_retry_refuse_invalid_spans() {
    let (spec, domain) = crate::phillips_drive::test_fixture(1);
    let cutter = Cutter::new(spec, domain).unwrap();
    let events = cutter.boundary_events().unwrap();
    let raw = events.trim_carrier(0, 1e-8, 4096).unwrap();
    let span = &raw[0];
    let parameters = span.parameters.map(|x| [x, x]);
    let count = Cell::new(0usize);
    let poll = || {
        count.set(count.get() + 1);
        false
    };
    let cache = events
        .prepare_trim_span_correspondence(0, parameters, MAX_WORK, MAX_BYTES, Some(&poll))
        .unwrap();
    assert!(std::ptr::eq(cache.events(), &events));
    assert!(std::ptr::eq(cache.events().source(), &cutter));
    assert_eq!(cache.trim_index(), 0);
    assert_eq!(cache.source_parameter_bounds(), parameters);
    assert_eq!(cache.construction_work_used(), TRIM_SPAN_SOURCE_WORK);
    for (work, bytes) in [
        (TRIM_SPAN_SOURCE_WORK - 1, MAX_BYTES),
        (MAX_WORK, cache.construction_peak_bytes() - 1),
    ] {
        assert!(matches!(
            events.prepare_trim_span_correspondence(0, parameters, work, bytes, None),
            Err(Error::Limit)
        ));
    }
    events
        .prepare_trim_span_correspondence(
            0,
            parameters,
            TRIM_SPAN_SOURCE_WORK,
            cache.construction_peak_bytes(),
            None,
        )
        .unwrap();
    for (index, p) in [
        (4, parameters),
        (0, [[-0.1, -0.1], [0.5, 0.5]]),
        (0, [[0.5, 0.5], [0.1, 0.1]]),
        (0, [[f64::NAN, f64::NAN], [1., 1.]]),
    ] {
        assert!(matches!(
            events.prepare_trim_span_correspondence(index, p, MAX_WORK, MAX_BYTES, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        events.prepare_trim_span_correspondence(
            0,
            [[0., 0.5], [0.4, 1.]],
            MAX_WORK,
            MAX_BYTES,
            None
        ),
        Err(Error::Unresolved)
    ));
    let controls = span
        .poles_mm
        .map(|p| p.map(|x| [x.next_down(), x.next_up()]));
    let axes = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let proof = cache
        .compare_controls(&controls, [0.; 3], axes, MAX_WORK, MAX_BYTES, None)
        .unwrap();
    assert!(std::ptr::eq(proof.source(), &cache));
    assert!(std::ptr::eq(proof.original_control_bounds(), &controls));
    assert_eq!(proof.new_work_used(), TRIM_SPAN_COMPARE_WORK);
    assert_eq!(proof.source_origin_mm(), [0.; 3]);
    assert_eq!(proof.source_columns(), axes);
    assert!(matches!(
        events.prepare_trim_span_correspondence(0, parameters, MAX_WORK + 1, MAX_BYTES, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        cache.compare_controls(&controls, [0.; 3], axes, MAX_WORK, MAX_BYTES + 1, None),
        Err(Error::Invalid)
    ));
    for (work, bytes) in [
        (TRIM_SPAN_COMPARE_WORK - 1, MAX_BYTES),
        (MAX_WORK, proof.accounted_peak_bytes() - 1),
    ] {
        assert!(matches!(
            cache.compare_controls(&controls, [0.; 3], axes, work, bytes, None),
            Err(Error::Limit)
        ));
    }
    let exact = cache
        .compare_controls(
            &controls,
            [0.; 3],
            axes,
            TRIM_SPAN_COMPARE_WORK,
            proof.accounted_peak_bytes(),
            None,
        )
        .unwrap();
    assert_eq!(
        exact.maximum_derivative_mm().to_bits(),
        proof.maximum_derivative_mm().to_bits()
    );
    let mut bad = controls;
    bad[0][0] = [1., 0.];
    assert!(matches!(
        cache.compare_controls(&bad, [0.; 3], axes, MAX_WORK, MAX_BYTES, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        cache.compare_controls(&controls, [f64::NAN; 3], axes, MAX_WORK, MAX_BYTES, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        cache.compare_controls(&controls, [1e13; 3], axes, MAX_WORK, MAX_BYTES, None),
        Err(Error::Limit)
    ));
    for stop in [1, count.get()] {
        let calls = Cell::new(0usize);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        assert!(matches!(
            events.prepare_trim_span_correspondence(
                0,
                parameters,
                MAX_WORK,
                MAX_BYTES,
                Some(&cancel)
            ),
            Err(Error::Cancelled)
        ));
    }
    let calls = Cell::new(0usize);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() >= 2
    };
    assert!(matches!(
        cache.compare_controls(&controls, [0.; 3], axes, MAX_WORK, MAX_BYTES, Some(&cancel)),
        Err(Error::Cancelled)
    ));
    let retry = cache
        .compare_controls(&controls, [0.; 3], axes, MAX_WORK, MAX_BYTES, None)
        .unwrap();
    assert_eq!(
        retry.maximum_position_mm().to_bits(),
        proof.maximum_position_mm().to_bits()
    );
    assert_eq!(
        retry.maximum_derivative_mm().to_bits(),
        proof.maximum_derivative_mm().to_bits()
    );
}
