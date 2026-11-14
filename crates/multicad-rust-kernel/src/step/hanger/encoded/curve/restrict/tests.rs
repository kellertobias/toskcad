use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use crate::step::hanger::encoded::cap::tests::{bezier, number};
use num_traits::{One, Signed, Zero};
type Q = BigRational;
fn jet(c: Vec<[Q; 4]>, u: &Q) -> ([Q; 3], [Q; 3]) {
    let n = Q::from_integer((c.len() - 1).into());
    let d: Vec<[Q; 4]> = c
        .windows(2)
        .map(|p| std::array::from_fn(|i| (&p[1][i] - &p[0][i]) * &n))
        .collect();
    let h = bezier(c, u);
    let d = bezier(d, u);
    (
        std::array::from_fn(|i| &h[i] / &h[3]),
        std::array::from_fn(|i| (&d[i] * &h[3] - &h[i] * &d[3]) / (&h[3] * &h[3])),
    )
}
#[test]
fn actual_finite_restrictions_preserve_exact_point_and_derivative_correspondence() {
    for (left, nb) in [(false, false), (true, true)] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let mut count = 0;
        let mut uncertain = 0;
        let mut work = 0;
        let mut incremental = 0;
        for index in 0..f.cells().len() {
            for end in [End::Start, End::End] {
                for piece in 0..f.angular.coefficients().len() {
                    let Some(curve) = f
                        .cap_curve_piece(index, end, piece, DomainOptions::default(), None)
                        .unwrap()
                    else {
                        continue;
                    };
                    let finite = curve
                        .restrict_finite(CurveOptions::default(), None)
                        .unwrap();
                    count += 1;
                    work += finite.work_used();
                    incremental += finite.work_used() - f.work_used() - curve.work_used;
                    assert!(std::ptr::eq(finite.source(), &curve));
                    assert!(finite.retained_bytes() > curve.retained_bytes);
                    let [a, b] = finite.meridian_parameters().each_ref().map(number);
                    assert!(a >= Q::zero() && b <= Q::one() && a < b);
                    let bounds = finite.bounds();
                    assert!(bounds.position_error_mm <= 1e-7 && bounds.derivative_error_mm <= 1e-6);
                    if bounds.meridian_parameter_error != [0., 0.] {
                        uncertain += 1;
                    }
                    let original = curve
                        .controls()
                        .iter()
                        .map(|p| p.each_ref().map(number))
                        .collect::<Vec<_>>();
                    let restricted = finite
                        .controls()
                        .iter()
                        .map(|p| p.each_ref().map(number))
                        .collect::<Vec<_>>();
                    for k in 0..=4 {
                        let t = Q::new(k.into(), 4.into());
                        let u = &a + (&b - &a) * &t;
                        let (x, dx) = jet(original.clone(), &u);
                        let (y, dy) = jet(restricted.clone(), &t);
                        assert_eq!(x, y);
                        assert_eq!(dy, dx.map(|d| d * (&b - &a)));
                    }
                }
            }
        }
        eprintln!(
            "finite cap restriction left={left} nb={nb} pieces={count} uncertain={uncertain} summed_standalone_work={work} restriction_incremental={incremental}"
        );
        assert!(count >= 44);
        if nb {
            assert!(uncertain > 0);
        } else {
            assert_eq!(uncertain, 0);
        }
    }
}
#[test]
fn exact_restriction_scratch_checks_reject_invalid_numeric_domains() {
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut m = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let line = [
        [N::integer(0), N::integer(0), N::integer(0), N::integer(1)],
        [N::integer(1), N::integer(0), N::integer(0), N::integer(1)],
    ];
    for (a, b) in [(0, 0), (1, 0), (-1, 1), (0, 2)] {
        assert!(matches!(
            restrict(&line, &N::integer(a), &N::integer(b), &mut m, &mut w),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        restrict::<4>(&[], &N::integer(0), &N::integer(1), &mut m, &mut w),
        Err(Error::Invalid)
    ));
    let huge = vec![line[0].clone(); 24];
    assert!(matches!(
        restrict(&huge, &N::integer(0), &N::integer(1), &mut m, &mut w),
        Err(Error::Invalid)
    ));
    let half = N(Q::new(1.into(), 2.into()));
    let quarter = N(Q::new(1.into(), 4.into()));
    let data = restrict(&line, &quarter, &half, &mut m, &mut w).unwrap();
    assert_eq!(data[0][0].0, quarter.0);
    assert_eq!(data[1][0].0, half.0);
    let bounds = derivatives(&line, &N::integer(1), &mut w).unwrap();
    assert!(bounds[0] >= 1. && bounds[0] < 1. + 1e-12 && bounds[1] <= 1e-15);
    let zero: [[N; 4]; 2] = std::array::from_fn(|_| std::array::from_fn(|_| N::integer(0)));
    assert!(matches!(
        derivatives(&zero, &N::integer(1), &mut w),
        Err(Error::Unresolved)
    ));
}
#[test]
fn actual_algebraic_endpoint_quality_cancel_and_capacity_failures_stay_explicit() {
    use std::cell::Cell;
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(true, true),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let index = f
        .cells()
        .iter()
        .position(|c| {
            c.source_key.total_period == 25
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
        .cap_curve_piece(index, End::End, 0, DomainOptions::default(), None)
        .unwrap()
        .unwrap();
    let observed = Cell::new(0usize);
    let count = || {
        observed.set(observed.get() + 1);
        false
    };
    let good = curve
        .restrict_finite(CurveOptions::default(), Some(&count))
        .unwrap();
    assert!(good.bounds().meridian_parameter_error[1] > 0.);
    for stop in (1..=256).chain([observed.get() / 2, observed.get()]) {
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        assert!(matches!(
            curve.restrict_finite(CurveOptions::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
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
                position_error_mm: 1e-30,
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
        assert!(matches!(curve.restrict_finite(options,None),Err(e) if e==expected));
    }
    assert!(matches!(
        curve.restrict_finite(
            CurveOptions {
                work: good.work_used() - 1,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(s.spec(), b.source().spec());
}
fn second(c: &[[Q; 4]], u: &Q) -> [Q; 3] {
    let n = Q::from_integer((c.len() - 1).into());
    let d: Vec<[Q; 4]> = c
        .windows(2)
        .map(|p| std::array::from_fn(|i| (&p[1][i] - &p[0][i]) * &n))
        .collect();
    let dd: Vec<[Q; 4]> = d
        .windows(2)
        .map(|p| std::array::from_fn(|i| (&p[1][i] - &p[0][i]) * (&n - Q::one())))
        .collect();
    let h = bezier(c.to_vec(), u);
    let d = bezier(d, u);
    let dd = bezier(dd, u);
    std::array::from_fn(|i| {
        &dd[i] / &h[3]
            - Q::from_integer(2.into()) * &d[i] * &d[3] / (&h[3] * &h[3])
            - &h[i] * &dd[3] / (&h[3] * &h[3])
            + Q::from_integer(2.into()) * &h[i] * &d[3] * &d[3] / (&h[3] * &h[3] * &h[3])
    })
}
#[test]
fn independently_bracketed_root_events_and_actual_curve_jets_check_chain_rule_bound() {
    for left in [false, true] {
        let s = Hanger::new(
            crate::authored_hanger_bolt::tests::fixture(left, true),
            None,
        )
        .unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let index = f
            .cells()
            .iter()
            .position(|c| {
                c.source_key.total_period == 25
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
            .cap_curve_piece(index, End::End, 0, DomainOptions::default(), None)
            .unwrap()
            .unwrap();
        let finite = curve
            .restrict_finite(CurveOptions::default(), None)
            .unwrap();
        let row = &f.meridians[f.cells[index].meridian];
        let qw = (0..=row.degree())
            .map(|i| [number(&row.phase[i]), number(&row.weights[i])])
            .collect::<Vec<_>>();
        let q_at = |u: &Q| {
            let h = bezier(qw.clone(), u);
            &h[0] / &h[1]
        };
        for (end, phase) in curve.phase.iter().enumerate() {
            let q = number(phase);
            let mut low = Q::zero();
            let mut high = Q::one();
            if q_at(&low) == q {
                high = low.clone();
            } else if q_at(&high) == q {
                low = high.clone();
            } else {
                for _ in 0..128 {
                    let mid = (&low + &high) / Q::from_integer(2.into());
                    match q_at(&mid).cmp(&q) {
                        Ordering::Less => low = mid,
                        Ordering::Greater => high = mid,
                        Ordering::Equal => {
                            low = mid.clone();
                            high = mid;
                            break;
                        }
                    }
                }
            }
            assert!(q_at(&low) <= q && q_at(&high) >= q);
            assert!(low >= Q::from_float(curve.parameter_bounds[end][0]).unwrap());
            assert!(high <= Q::from_float(curve.parameter_bounds[end][1]).unwrap());
            let chosen = number(&finite.meridian_parameters()[end]);
            let error = (&chosen - &low).abs().max((&chosen - &high).abs());
            assert!(error <= Q::from_float(finite.bounds().meridian_parameter_error[end]).unwrap());
        }
        let mut w = Work::new(16_000_000, None).unwrap();
        let mut m = Memory::new(0, 64 * 1024 * 1024).unwrap();
        let low = binary(curve.parameter_bounds[0][0], &mut w).unwrap();
        let high = binary(curve.parameter_bounds[1][1], &mut w).unwrap();
        let new = factor_derivatives(&curve, &low, &high, &mut m, &mut w).unwrap();
        let input = curve
            .controls()
            .iter()
            .map(|p| p.each_ref().map(|x| N(number(x))))
            .collect::<Vec<_>>();
        let old = derivatives(&input, &N::integer(1), &mut w).unwrap();
        let data = curve
            .controls()
            .iter()
            .map(|p| p.each_ref().map(number))
            .collect::<Vec<_>>();
        for k in 0..=8 {
            let u = &low.0 + (&high.0 - &low.0) * Q::new(k.into(), 8.into());
            let (_, d) = jet(data.clone(), &u);
            let dd = second(&data, &u);
            let d = d.into_iter().map(|x| x.abs()).sum::<Q>();
            let dd = dd.into_iter().map(|x| x.abs()).sum::<Q>();
            assert!(d <= Q::from_float(new[0]).unwrap() && dd <= Q::from_float(new[1]).unwrap());
            assert!(d <= Q::from_float(old[0]).unwrap() && dd <= Q::from_float(old[1]).unwrap());
        }
        let bad = [
            [N::integer(0), N::integer(0), N::integer(-1)],
            [N::integer(0), N::integer(0), N::integer(-1)],
        ];
        assert!(matches!(
            scalar_jet(&bad, 0, &N::integer(1), &mut w),
            Err(Error::Unresolved)
        ));
        let outside = amplification(I::new(-0.25, 1.5), 10, &mut w).unwrap();
        assert!(outside.lo >= 1. && outside.hi >= 1024.);
    }
}
#[test]
fn numeric_chain_overflow_is_refused_without_a_source_geometry_waiver() {
    let mut w = Work::new(16_000_000, None).unwrap();
    let r = ScalarJet {
        value: p(f64::MAX),
        first: 1.,
        second: 1.,
    };
    let q = ScalarJet {
        value: p(1.),
        first: 1.,
        second: 1.,
    };
    assert!(matches!(
        chain_bound(r, q, p(1.), [p(1.), p(1.), p(1.)], &mut w),
        Err(Error::Unresolved)
    ));
}
#[test]
fn actual_wood_start_chain_and_every_prefix_cancel_keep_original_support() {
    use std::cell::Cell;
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, true),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let curve = f
        .cells()
        .iter()
        .enumerate()
        .filter(|(_, c)| {
            matches!(
                c.source_key.face,
                FaceKey::WoodThread {
                    band: ThreadBandRole::Upper,
                    ..
                }
            )
        })
        .find_map(|(i, _)| {
            (0..f.angular.coefficients().len()).find_map(|p| {
                f.cap_curve_piece(i, End::Start, p, DomainOptions::default(), None)
                    .unwrap()
            })
        })
        .unwrap();
    let low = number(
        &curve
            .restrict_finite(CurveOptions::default(), None)
            .unwrap()
            .meridian_parameters()[0],
    );
    let high = Q::from_float(curve.parameter_bounds[1][1]).unwrap();
    let observed = Cell::new(0usize);
    let count = || {
        observed.set(observed.get() + 1);
        false
    };
    let mut w = Work::new(16_000_000, Some(&count)).unwrap();
    let mut m = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let bounds =
        factor_derivatives(&curve, &N(low.clone()), &N(high.clone()), &mut m, &mut w).unwrap();
    let data = curve
        .controls()
        .iter()
        .map(|p| p.each_ref().map(number))
        .collect::<Vec<_>>();
    for k in 0..=8 {
        let u = &low + (&high - &low) * Q::new(k.into(), 8.into());
        let (_, d) = jet(data.clone(), &u);
        let dd = second(&data, &u);
        assert!(d.into_iter().map(|x| x.abs()).sum::<Q>() <= Q::from_float(bounds[0]).unwrap());
        assert!(dd.into_iter().map(|x| x.abs()).sum::<Q>() <= Q::from_float(bounds[1]).unwrap());
    }
    for stop in (1..=observed.get().min(1024)).chain([observed.get() / 2, observed.get()]) {
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        let mut w = Work::new(16_000_000, Some(&cancel)).unwrap();
        let mut m = Memory::new(0, 64 * 1024 * 1024).unwrap();
        assert!(
            matches!(
                factor_derivatives(&curve, &N(low.clone()), &N(high.clone()), &mut m, &mut w),
                Err(Error::Cancelled)
            ),
            "stop={stop}"
        );
    }
}
