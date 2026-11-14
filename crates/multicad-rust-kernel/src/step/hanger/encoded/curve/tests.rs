use super::super::cap::tests::{bezier, number};
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use num_traits::{One, Zero};
type Q = BigRational;
fn value_derivative<const D: usize>(c: Vec<[Q; D]>, u: &Q) -> ([Q; D], [Q; D]) {
    let n = c.len() - 1;
    let dc = c
        .windows(2)
        .map(|p| std::array::from_fn(|i| (&p[1][i] - &p[0][i]) * Q::from_integer(n.into())))
        .collect();
    (bezier(c, u), bezier(dc, u))
}
fn quotient(h: [Q; 4], d: [Q; 4]) -> ([Q; 3], [Q; 3]) {
    (
        std::array::from_fn(|i| &h[i] / &h[3]),
        std::array::from_fn(|i| (&d[i] * &h[3] - &h[i] * &d[3]) / (&h[3] * &h[3])),
    )
}
#[test]
fn actual_root_and_flat_cap_curves_equal_same_surface_and_first_derivatives() {
    for (left, nb) in [(false, false), (true, true)] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let mut count = 0;
        let mut work = 0;
        for index in 0..f.cells().len() {
            for end in [End::Start, End::End] {
                let Some(trace) = f
                    .cap_trace(index, end, DomainOptions::default(), None)
                    .unwrap()
                else {
                    continue;
                };
                let row = &f.meridians[f.cells[index].meridian];
                let spec = if row.component == Component::Metric {
                    s.metric().spec().pitch_mm
                } else {
                    s.wood().spec().pitch_mm
                };
                let starts = if row.component == Component::Metric {
                    s.metric().spec().starts
                } else {
                    s.wood().spec().starts
                };
                let lead = Q::from_float(spec).unwrap() * Q::from_integer(starts.into());
                let station = if end == End::Start {
                    Q::zero()
                } else {
                    Q::from_float(if row.component == Component::Metric {
                        s.metric().full_length_mm()
                    } else {
                        s.wood().full_length_mm()
                    })
                    .unwrap()
                };
                let shifted = station
                    - Q::from_float(spec).unwrap()
                        * Q::from_integer(f.cells[index].source_key.total_period.into());
                let mut intervals = Vec::new();
                for piece in 0..f.angular.coefficients().len() {
                    let Some(curve) = f
                        .cap_curve_piece(index, end, piece, DomainOptions::default(), None)
                        .unwrap()
                    else {
                        continue;
                    };
                    count += 1;
                    work += curve.work_used;
                    assert_eq!(curve.degree(), 11 * row.degree());
                    assert!(std::ptr::eq(curve.field(), &f));
                    assert_eq!(curve.source_cell(), index);
                    assert_eq!(curve.angular_piece(), piece);
                    assert_eq!(curve.end(), end);
                    assert!(curve.retained_bytes() > 0);
                    assert!(curve.work_used() > 0);
                    assert_eq!(curve.parameter_bounds(), curve.parameter_bounds);
                    assert_eq!(curve.phase(), &curve.phase);
                    intervals.push(curve.phase.each_ref().map(number));
                    let ua = Q::from_float(curve.parameter_bounds[0][1]).unwrap();
                    let ub = Q::from_float(curve.parameter_bounds[1][0]).unwrap();
                    for numerator in [1, 2, 3] {
                        let u = (&ua * Q::from_integer((4 - numerator).into())
                            + &ub * Q::from_integer(numerator.into()))
                            / Q::from_integer(4.into());
                        let (h, dh) = value_derivative(
                            curve
                                .controls()
                                .iter()
                                .map(|v| v.each_ref().map(number))
                                .collect(),
                            &u,
                        );
                        let (point, tangent) = quotient(h, dh);
                        let (rqw, drqw) = value_derivative(
                            (0..=row.degree())
                                .map(|j| {
                                    [
                                        number(&row.radius[j]),
                                        number(&row.phase[j]),
                                        number(&row.weights[j]),
                                    ]
                                })
                                .collect(),
                            &u,
                        );
                        let r = &rqw[0] / &rqw[2];
                        let q = &rqw[1] / &rqw[2];
                        let dr = (&drqw[0] * &rqw[2] - &rqw[0] * &drqw[2]) / (&rqw[2] * &rqw[2]);
                        let dq = (&drqw[1] * &rqw[2] - &rqw[1] * &drqw[2]) / (&rqw[2] * &rqw[2]);
                        let v = (&shifted - q) / &lead;
                        let a = Q::from_float(f.angular.knots()[piece]).unwrap();
                        let width = Q::from_float(f.angular.knots()[piece + 1]).unwrap() - &a;
                        let t = (v - a) / &width;
                        assert!(t >= Q::zero() && t <= Q::one());
                        let (actual, reverse) = if left {
                            (f.angular.coefficients().len() - 1 - piece, true)
                        } else {
                            (piece, false)
                        };
                        let cc = (0..=10)
                            .map(|j| {
                                f.angular.coefficients()[actual][if reverse { 10 - j } else { j }]
                                    .clone()
                            })
                            .collect();
                        let (xy, dxy) = value_derivative(cc, &t);
                        let dv = -dq / &lead;
                        assert_eq!(point, [&r * &xy[0], &r * &xy[1], number(&trace.plane_z)]);
                        assert_eq!(
                            tangent,
                            [
                                &dr * &xy[0] + &r * &dxy[0] * &dv / &width,
                                &dr * &xy[1] + &r * &dxy[1] * &dv / &width,
                                Q::zero()
                            ]
                        );
                    }
                }
                intervals.sort_by(|a, b| a[0].cmp(&b[0]));
                assert_eq!(intervals[0][0], number(&trace.endpoints[0].phase));
                assert_eq!(
                    intervals.last().unwrap()[1],
                    number(&trace.endpoints[1].phase)
                );
                for pair in intervals.windows(2) {
                    assert_eq!(pair[0][1], pair[1][0]);
                }
            }
        }
        eprintln!(
            "actual cap curves left={left} nb={nb} pieces={count} summed_individual_query_work={work}"
        );
        assert!(count > 20);
    }
}
#[test]
fn bounded_polynomial_operations_and_errors_do_not_publish_partial_curves() {
    let mut w = Work::new(16_000_000, None).unwrap();
    assert!(matches!(
        Polynomial::coefficients(&[], &mut w),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        Polynomial::coefficients(&vec![N::integer(1); 24], &mut w),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        Polynomial::zero(22).product(&Polynomial::zero(1), &mut w),
        Err(Error::Limit)
    ));
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    assert!(matches!(
        f.cap_curve_piece(0, End::Start, usize::MAX, DomainOptions::default(), None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        f.cap_curve_piece(usize::MAX, End::Start, 0, DomainOptions::default(), None),
        Err(Error::Invalid)
    ));
    let options = DomainOptions {
        maximum_bytes: 1,
        ..DomainOptions::default()
    };
    assert!(matches!(
        f.cap_curve_piece(0, End::Start, 0, options, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        f.cap_curve_piece(0, End::Start, 0, DomainOptions::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
}
#[test]
fn actual_cap_prefix_cancellation_resource_and_positive_sliver_refusals_are_atomic() {
    use std::cell::Cell;
    let mut spec = crate::authored_hanger_bolt::tests::fixture(false, false);
    spec.metric = AuthoredMetricThread::new(crate::authored_metric_thread::Spec {
        major_radius_mm: 1.46,
        pitch_mm: 0.23,
        length_mm: 0.43125,
        crest_width_mm: 0.023,
        starts: 1,
        left_hand: false,
    })
    .unwrap();
    let s = Hanger::new(spec, None).unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let index = f
        .cells()
        .iter()
        .position(|c| {
            c.source_key.total_period == 1
                && matches!(
                    c.source_key.face,
                    FaceKey::MetricThread {
                        band: ThreadBandRole::Root,
                        ..
                    }
                )
        })
        .unwrap();
    assert!(
        f.cap_trace(index, End::End, DomainOptions::default(), None)
            .unwrap()
            .is_some()
    );
    let outcome = f.cap_curve_piece(index, End::End, 4, DomainOptions::default(), None);
    assert!(
        matches!(outcome, Err(Error::Unresolved)),
        "positive source/root angular sliver is not welded"
    );
    let vacant = f
        .cells()
        .iter()
        .enumerate()
        .find_map(|(i, _)| {
            f.cap_trace(i, End::Start, DomainOptions::default(), None)
                .unwrap()
                .is_none()
                .then_some(i)
        })
        .unwrap();
    assert!(
        f.cap_curve_piece(vacant, End::Start, 0, DomainOptions::default(), None)
            .unwrap()
            .is_none()
    );
    assert!(matches!(
        f.cap_curve_piece(
            index,
            End::End,
            2,
            DomainOptions {
                maximum_work: 1,
                ..DomainOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let observed = Cell::new(0usize);
    let count = || {
        observed.set(observed.get() + 1);
        false
    };
    let successful = f
        .cap_curve_piece(index, End::End, 2, DomainOptions::default(), Some(&count))
        .unwrap()
        .unwrap();
    for stop in (1..=512.min(observed.get())).chain([observed.get() / 2, observed.get()]) {
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        assert!(
            matches!(
                f.cap_curve_piece(index, End::End, 2, DomainOptions::default(), Some(&cancel)),
                Err(Error::Cancelled)
            ),
            "stop={stop}"
        );
    }
    assert_eq!(s.spec(), b.source().spec());
    assert!(successful.work_used > 0);
}
