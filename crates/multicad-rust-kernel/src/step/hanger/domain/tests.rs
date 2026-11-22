use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
fn overlap(a: I, b: I) {
    assert!(
        a.valid() && b.valid() && a.lo <= b.hi && b.lo <= a.hi,
        "{a:?} {b:?}"
    );
}
#[test]
fn actual_graph_period_hand_cap_and_rail_parameters_bind_original_curve_c0_c1_c2() {
    for left in [false, true] {
        for nb in [false, true] {
            for starts in [1, 3] {
                let mut spec = crate::authored_hanger_bolt::tests::fixture(left, nb);
                spec.metric = AuthoredMetricThread::new(crate::authored_metric_thread::Spec {
                    starts,
                    ..spec.metric.spec()
                })
                .unwrap();
                let s = Hanger::new(spec, None).unwrap();
                let b = s.boundary(BoundaryOptions::default(), None).unwrap();
                for start in 0..starts {
                    let face = FaceKey::MetricThread {
                        start,
                        band: ThreadBandRole::Root,
                    };
                    let d =
                        RootFaceDomain::prepare(&b, face, DomainOptions::default(), None).unwrap();
                    assert!(std::ptr::eq(d.source_boundary(), &b));
                    assert_eq!(d.face_key(), face);
                    assert_eq!(d.same_sense(), left);
                    assert!(
                        d.phase_period() < starts
                            && d.work_used() <= 16_000_000
                            && d.retained_bytes() >= b.retained_bytes()
                    );
                    let turns = d.physical_turns();
                    let span = p(turns[1]) - p(turns[0]);
                    let lead = p(s.metric().spec().pitch_mm) * p(f64::from(starts));
                    let rate = I::new(
                        std::f64::consts::TAU.next_down(),
                        std::f64::consts::TAU.next_up(),
                    ) * p(if left { -1. } else { 1. });
                    for curve in b.curves().iter().filter(|c| c.faces.contains(&face)) {
                        for t in [0., 0.25, 0.5, 0.75, 1.] {
                            let uv = d
                                .trim_bounds(curve.key, [t, t], DomainOptions::default(), None)
                                .unwrap();
                            let original = b
                                .curve_bounds(curve.key, uv.original_parameter, None)
                                .unwrap();
                            let root = s
                                .metric()
                                .root_meridian()
                                .bounds(uv.point_uv[0], None)
                                .unwrap();
                            let [q, r] = root.point_mm.map(hull);
                            let [qu, ru] = root.first_derivative_mm.map(hull);
                            let [quu, ruu] = root.second_derivative_mm.map(hull);
                            let [ut, vt] = uv.first_derivative.map(hull);
                            let [utt, vtt] = uv.second_derivative.map(hull);
                            let v = p(turns[0]) + span * hull(uv.point_uv[1]);
                            let theta = rate * v;
                            let tht = rate * span * vt;
                            let thtt = rate * span * vtt;
                            let (c, sin) = (theta.cos(), theta.sin());
                            let rt = ru * ut;
                            let rtt = ruu * ut.square() + ru * utt;
                            let qt = qu * ut;
                            let qtt = quu * ut.square() + qu * utt;
                            let point = [r * c, r * sin, q + lead * v + hull(d.axial_offset_mm())];
                            let first = [
                                rt * c - r * sin * tht,
                                rt * sin + r * c * tht,
                                qt + lead * span * vt,
                            ];
                            let aa = rtt - r * tht.square();
                            let bb = p(2.) * rt * tht + r * thtt;
                            let second = [
                                aa * c - bb * sin,
                                aa * sin + bb * c,
                                qtt + lead * span * vtt,
                            ];
                            let parameter_first = hull(uv.original_parameter_first);
                            let parameter_second = hull(uv.original_parameter_second);
                            for i in 0..3 {
                                overlap(point[i], hull(original.point_mm[i]));
                                overlap(
                                    first[i],
                                    hull(original.first_derivative_mm[i]) * parameter_first,
                                );
                                overlap(
                                    second[i],
                                    hull(original.second_derivative_mm[i])
                                        * parameter_first.square()
                                        + hull(original.first_derivative_mm[i]) * parameter_second,
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
#[test]
fn real_wrong_roles_other_periods_precision_resources_and_cancellation_refuse_without_mutation() {
    use std::cell::Cell;
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(true, true),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let opts = DomainOptions::default();
    let face = FaceKey::MetricThread {
        start: 0,
        band: ThreadBandRole::Root,
    };
    let d = RootFaceDomain::prepare(&b, face, opts, None).unwrap();
    let key = CurveKey::Cap {
        component: Component::Metric,
        end: End::Start,
        period: d.phase_period(),
        band: ThreadBandRole::Root,
    };
    for face in [
        FaceKey::PlainWall,
        FaceKey::MetricThread {
            start: 99,
            band: ThreadBandRole::Root,
        },
    ] {
        assert!(matches!(
            RootFaceDomain::prepare(&b, face, opts, None),
            Err(Error::Invalid)
        ));
    }
    for key in [
        CurveKey::PlainCircle { end: End::Start },
        CurveKey::Rail {
            component: Component::Metric,
            period: 0,
            knot: Knot::CrestStart,
        },
        CurveKey::Cap {
            component: Component::Wood,
            end: End::Start,
            period: 0,
            band: ThreadBandRole::Root,
        },
    ] {
        assert_eq!(
            d.trim_bounds(key, [0., 1.], opts, None),
            Err(Error::Invalid)
        );
    }
    for parameter in [[f64::NAN, 1.], [-1., 0.], [1., 0.], [0., 1_f64.next_up()]] {
        assert_eq!(
            d.trim_bounds(key, parameter, opts, None),
            Err(Error::Invalid)
        );
    }
    assert!(matches!(
        RootFaceDomain::prepare(
            &b,
            face,
            DomainOptions {
                maximum_work: 1,
                ..opts
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        RootFaceDomain::prepare(
            &b,
            face,
            DomainOptions {
                maximum_work: 0,
                ..opts
            },
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        RootFaceDomain::prepare(
            &b,
            face,
            DomainOptions {
                maximum_bytes: 0,
                ..opts
            },
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        RootFaceDomain::prepare(
            &b,
            face,
            DomainOptions {
                maximum_bytes: 1,
                ..opts
            },
            None
        ),
        Err(Error::Limit)
    ));
    let exact = d.retained_bytes() + SCRATCH;
    assert!(
        d.trim_bounds(
            key,
            [0., 1.],
            DomainOptions {
                maximum_bytes: exact,
                ..opts
            },
            None
        )
        .is_ok()
    );
    assert_eq!(
        d.trim_bounds(
            key,
            [0., 1.],
            DomainOptions {
                maximum_bytes: exact - 1,
                ..opts
            },
            None
        ),
        Err(Error::Limit)
    );
    let count = Cell::new(0);
    d.trim_bounds(
        key,
        [0., 1.],
        opts,
        Some(&|| {
            count.set(count.get() + 1);
            false
        }),
    )
    .unwrap();
    for stop in [0, 1, count.get() / 2, count.get() - 1] {
        let n = Cell::new(0);
        assert_eq!(
            d.trim_bounds(
                key,
                [0., 1.],
                opts,
                Some(&|| {
                    let at = n.get();
                    n.set(at + 1);
                    at == stop
                })
            ),
            Err(Error::Cancelled)
        );
    }
    let count = Cell::new(0);
    RootFaceDomain::prepare(
        &b,
        face,
        opts,
        Some(&|| {
            count.set(count.get() + 1);
            false
        }),
    )
    .unwrap();
    for stop in [0, 1, count.get() / 2, count.get() - 1] {
        let n = Cell::new(0);
        assert!(matches!(
            RootFaceDomain::prepare(
                &b,
                face,
                opts,
                Some(&|| {
                    let at = n.get();
                    n.set(at + 1);
                    at == stop
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(*b.source(), s);
}
