use super::*;
use std::cell::Cell;
fn cubic(p: [[f64; 3]; 4], t: f64) -> [f64; 3] {
    let mix = |a: [f64; 3], b: [f64; 3]| std::array::from_fn(|i| a[i] * (1. - t) + b[i] * t);
    let a = mix(p[0], p[1]);
    let b = mix(p[1], p[2]);
    let c = mix(p[2], p[3]);
    mix(mix(a, b), mix(b, c))
}
fn cubic_interval(p: [[f64; 3]; 4], t: f64) -> [I; 3] {
    let t = I::point(t);
    let mix = |a: [I; 3], b: [I; 3]| std::array::from_fn(|i| a[i] * (I::point(1.) - t) + b[i] * t);
    let [a, b, c, d] = p.map(|p| p.map(I::point));
    let [a, b, c] = [mix(a, b), mix(b, c), mix(c, d)];
    mix(mix(a, b), mix(b, c))
}
#[test]
fn metric_trim_carriers_bind_whole_authored_curves_and_shared_events() {
    for fixture in 0..=5 {
        let (s, d) = super::super::super::tests::fixture(fixture);
        let cutter = Cutter::new(s, d).unwrap();
        let graph = cutter.boundary_events().unwrap();
        let now = std::time::Instant::now();
        let mut counts = Vec::new();
        for index in 0..4 {
            let spans = graph.trim_carrier(index, 1e-8, 4096).unwrap();
            counts.push(spans.len());
            assert!(!spans.is_empty());
            let ends = graph.trims[index].ends;
            assert_eq!(
                spans[0].poles_mm[0],
                graph.event(ends[0]).point.map(I::midpoint)
            );
            assert_eq!(
                spans.last().unwrap().poles_mm[3],
                graph.event(ends[1]).point.map(I::midpoint)
            );
            for pair in spans.windows(2) {
                assert_eq!(pair[0].parameters[1], pair[1].parameters[0]);
                assert_eq!(pair[0].poles_mm[3], pair[1].poles_mm[0]);
            }
            for span in &spans {
                let bound = graph.certify_trim_span(index, span).unwrap();
                assert!(bound.source_curve_error_bound_mm() <= 1e-8);
                for j in 0..=8 {
                    let u = j as f64 / 8.;
                    let t = span.parameters[0] + u * (span.parameters[1] - span.parameters[0]);
                    let source = graph.trim_point_mm_bounds(index, t).unwrap().map(b);
                    let point = cubic(span.poles_mm, u);
                    let delta = norm(std::array::from_fn(|i| I::point(point[i]) - source[i]));
                    assert!(
                        delta.hi <= bound.source_curve_error_bound_mm() + 1e-11,
                        "{fixture}:{index} {delta:?} {bound:?}"
                    );
                }
            }
            let mut changed = spans[0];
            changed.poles_mm[1][0] += 0.01;
            let error = graph.certify_trim_span(index, &changed).unwrap();
            assert!(error.control_error_mm > 0.0099);
        }
        eprintln!(
            "Phillips fixture{fixture} carrier spans={counts:?} elapsed={:?}",
            now.elapsed()
        );
    }
}
#[test]
fn literal_span_validation_work_cancel_and_precision_exhaustion_are_observable() {
    let (s, d) = super::super::super::tests::fixture(1);
    let cutter = Cutter::new(s, d).unwrap();
    let graph = cutter.boundary_events().unwrap();
    let before = cutter.clone();
    let valid = graph.trim_carrier(0, 1e-7, 4096).unwrap()[0];
    assert!(matches!(
        graph.certify_trim_span(4, &valid),
        Err(Error::Invalid)
    ));
    for parameters in [[0., 0.], [0.5, 0.25], [-1., 0.5], [0., 2.], [f64::NAN, 1.]] {
        assert!(matches!(
            graph.certify_trim_span(
                0,
                &CubicTrimSpan {
                    parameters,
                    ..valid
                }
            ),
            Err(Error::Invalid)
        ))
    }
    let mut bad = valid;
    bad.poles_mm[2][1] = f64::NAN;
    assert!(matches!(
        graph.certify_trim_span(0, &bad),
        Err(Error::Invalid)
    ));
    bad.poles_mm[2][1] = 1e10;
    assert!(matches!(
        graph.certify_trim_span(0, &bad),
        Err(Error::Limit)
    ));
    for (index, error, limit, expected) in [
        (4, 1e-7, 10, Error::Invalid),
        (0, 0., 10, Error::Invalid),
        (0, f64::NAN, 10, Error::Invalid),
        (0, 1e-7, 0, Error::Limit),
        (0, 1e-7, 1_000_001, Error::Limit),
        (0, 1e-8, 1, Error::Limit),
    ] {
        assert_eq!(graph.trim_carrier(index, error, limit), Err(expected))
    }
    assert_eq!(
        graph.trim_carrier_with_cancel(0, 1e-8, 4096, &|| true),
        Err(Error::Cancelled)
    );
    let calls = Cell::new(0);
    assert_eq!(
        graph.trim_carrier_with_cancel(0, 1e-8, 4096, &|| {
            calls.set(calls.get() + 1);
            calls.get() > 5
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(graph.trim_carrier(0, 1e-30, 2000), Err(Error::Unresolved));
    assert_eq!(cutter, before);
    assert!(graph.trim_carrier(0, 1e-7, 4096).is_ok());
}

#[test]
fn exact_quartic_remainder_and_overflow_have_independent_numeric_controls() {
    // f(t)=t^4 has Hermite endpoint controls[0,0,-1/3,1]. Its error
    // f-H=t²(1-t)² reaches1/16, exactly the fourth-derivative certificate.
    let bound =
        hermite_error(0., [I::point(1.), I::point(0.), I::point(0.)], I::point(1.)).unwrap();
    let literal = cubic(
        [[0., 0., 0.], [0., 0., 0.], [-1. / 3., 0., 0.], [1., 0., 0.]],
        0.5,
    )[0];
    let actual = 0.5_f64.powi(4) - literal;
    assert!((actual - 1. / 16.).abs() < 1e-16);
    assert!(bound.hermite_remainder_mm >= actual);
    assert!(bound.hermite_remainder_mm < actual + 1e-14);
    // A finite fourth coefficient can overflow its Euclidean norm calculation;
    // it must refuse a certificate, even with finite literal controls.
    assert_eq!(
        hermite_error(
            0.,
            [I::point(1e200), I::point(0.), I::point(0.)],
            I::point(1.)
        ),
        Err(Error::Unresolved)
    );
}
#[test]
fn literal_world_controls_bind_actual_rotated_reflected_far_and_affine_source_curves() {
    let identity = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let rolled = [[0.36, 0.48, 0.8], [-0.8, 0.6, 0.], [-0.48, -0.64, 0.6]];
    let frames = [
        ([0.; 3], identity),
        ([7., -4., 2.], rolled),
        ([0.; 3], [[1., 0., 0.], [0., -1., 0.], [0., 0., 1.]]),
        ([1e6, -2e6, 3e6], rolled),
        ([0.; 3], [[1.2, 0.1, 0.], [0.05, 0.9, 0.2], [0., 0.01, 0.8]]),
    ];
    for fixture in [0, 1] {
        let (s, d) = super::super::super::tests::fixture(fixture);
        let source = Cutter::new(s, d).unwrap();
        let events = source.boundary_events().unwrap();
        for index in 0..4 {
            let spans = events.trim_carrier(index, 1e-7, 4096).unwrap();
            for (origin, axes) in frames {
                let place = |p: [I; 3]| {
                    std::array::from_fn::<_, 3, _>(|i| {
                        I::point(origin[i])
                            + (0..3).fold(I::point(0.), |s, j| s + I::point(axes[j][i]) * p[j])
                    })
                };
                let nominal = |p: [f64; 3]| {
                    std::array::from_fn::<_, 3, _>(|i| {
                        origin[i] + axes[0][i] * p[0] + axes[1][i] * p[1] + axes[2][i] * p[2]
                    })
                };
                for span in &spans {
                    let placed = CubicTrimSpan {
                        parameters: span.parameters,
                        poles_mm: span.poles_mm.map(nominal),
                    };
                    let bound = events
                        .certify_trim_span_in_frame(index, &placed, origin, axes)
                        .unwrap();
                    assert!(bound.source_curve_error_bound_mm() < 1e-5);
                    for t in [0., 0.25, 0.5, 0.75, 1.] {
                        let parameter =
                            span.parameters[0] + t * (span.parameters[1] - span.parameters[0]);
                        let exact = place(
                            events
                                .trim_point_mm_bounds(index, parameter)
                                .unwrap()
                                .map(b),
                        );
                        let evaluated = cubic(placed.poles_mm, t);
                        let residual =
                            norm(std::array::from_fn(|i| I::point(evaluated[i]) - exact[i]));
                        // The independent point/de Casteljau evaluation has its
                        // own arithmetic reserve; it is not source uncertainty.
                        let source_radius = norm(std::array::from_fn(|i| {
                            exact[i] - I::point(exact[i].midpoint())
                        }));
                        let polynomial = cubic_interval(placed.poles_mm, t);
                        let evaluation_error = norm(std::array::from_fn(|i| {
                            polynomial[i] - I::point(evaluated[i])
                        }));
                        let threshold = I::point(bound.source_curve_error_bound_mm())
                            + I::point(2.) * source_radius
                            + evaluation_error;
                        assert!(residual.hi <= threshold.hi);
                    }
                    let mut changed = placed;
                    changed.poles_mm[1][0] += 0.001;
                    assert!(
                        events
                            .certify_trim_span_in_frame(index, &changed, origin, axes)
                            .unwrap()
                            .control_error_mm
                            > 0.00099
                    );
                }
            }
        }
    }
}
#[test]
fn posed_trim_certification_refuses_nonfinite_parameters_and_explicit_domain_limits() {
    let (s, d) = super::super::super::tests::fixture(1);
    let source = Cutter::new(s, d).unwrap();
    let events = source.boundary_events().unwrap();
    let span = events.trim_carrier(0, 1e-6, 4096).unwrap()[0];
    let axes = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    for parameters in [[f64::NAN, 1.], [0.5, 0.5], [-0.1, 0.5], [0., 1.1]] {
        assert_eq!(
            events.certify_trim_span_in_frame(
                0,
                &CubicTrimSpan { parameters, ..span },
                [0.; 3],
                axes
            ),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        events.certify_trim_span_in_frame(4, &span, [0.; 3], axes),
        Err(Error::Invalid)
    );
    assert_eq!(
        events.certify_trim_span_in_frame(0, &span, [f64::NAN; 3], axes),
        Err(Error::Invalid)
    );
    assert_eq!(
        events.certify_trim_span_in_frame(0, &span, [0.; 3], [[f64::NAN; 3]; 3]),
        Err(Error::Invalid)
    );
    let mut changed = span;
    changed.poles_mm[1][0] = f64::INFINITY;
    assert_eq!(
        events.certify_trim_span_in_frame(0, &changed, [0.; 3], axes),
        Err(Error::Invalid)
    );
    changed.poles_mm[1][0] = 1e13;
    assert_eq!(
        events.certify_trim_span_in_frame(0, &changed, [0.; 3], axes),
        Err(Error::Limit)
    );
    assert_eq!(
        events.certify_trim_span_in_frame(0, &span, [1e13; 3], axes),
        Err(Error::Limit)
    );
    assert_eq!(
        events.certify_trim_span_in_frame(0, &span, [0.; 3], [[1e10; 3]; 3]),
        Err(Error::Limit)
    );
}
