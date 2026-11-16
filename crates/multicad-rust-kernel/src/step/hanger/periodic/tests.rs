use super::*;
#[test]
fn actual_both_handed_literal_templates_have_one_exact_periodic_xy_owner() {
    for left in [false, true] {
        let q = PeriodicAngular::prepare(left, CurveOptions::default(), None).unwrap();
        assert_eq!(q.left_hand(), left);
        assert_eq!(q.poles()[0], [1., 0., 0.]);
        assert_eq!(*q.poles().last().unwrap(), [1., 0., 1.]);
        assert_eq!(q.knots()[0], 0.);
        assert_eq!(*q.knots().last().unwrap(), 1.);
        assert_eq!(q.poles().len(), 11 * (q.knots().len() - 1) + 1);
        assert!(q.position_error() <= 1e-7 && q.derivative_error() <= 1e-6);
        assert!(q.signed_angular_jacobian_lower() > 6.);
        assert!(q.work_used() <= 16_000_000 && q.retained_bytes() <= 64 * 1024 * 1024);
        // Independent evaluation of actual polynomial coefficients and its
        // derivative, including both ONE-SIDED derivatives of every knot.
        for (span, ab) in q.knots().windows(2).enumerate() {
            for t in [0., 0.25, 0.5, 0.75, 1.] {
                let data = q.poles()[11 * span..11 * span + 12].to_vec();
                let derivative = data
                    .windows(2)
                    .map(|v| std::array::from_fn(|k| 11. * (v[1][k] - v[0][k]) / (ab[1] - ab[0])))
                    .collect();
                let actual = super::super::tests::bezier(data, t);
                let first = super::super::tests::bezier(derivative, t);
                let sign = if left { -1. } else { 1. };
                let jac = sign * (actual[0] * first[1] - actual[1] * first[0]);
                assert!(jac + 1e-10 >= q.signed_angular_jacobian_lower());
                let physical = ab[0] + (ab[1] - ab[0]) * t;
                let theta = sign * std::f64::consts::TAU * physical;
                assert!(
                    (actual[0] - theta.cos()).abs()
                        + (actual[1] - theta.sin()).abs()
                        + (actual[2] - physical).abs()
                        <= q.position_error() + 1e-14
                );
            }
        }
    }
}
#[test]
fn degree_and_projection_predicate_refuses_actual_double_wrap_and_tangent_fold_enclosures() {
    // B(t)=exp(4*pi*i*t): |B-C|<=2, |B'-C'|<=6*pi. This
    // is a genuine closed double wrap, excluded by the homotopy guard.
    assert_eq!(
        radial_proof(2., 3. * std::f64::consts::TAU),
        Err(Error::Unresolved)
    );
    // B=C+0.3*sin(4*tau*t)*i*C. Uniform errors <=0.3 and
    // 1.5*tau are genuine bounds. At t=1/8 its angular derivative
    // is tau*(1-1.2)<0, so C0 closeness alone does not prove a radial graph.
    let t = 0.125;
    let tau = std::f64::consts::TAU;
    let angle_derivative =
        tau + 1.2 * tau * (4. * tau * t).cos() / (1. + 0.09 * (4. * tau * t).sin().powi(2));
    assert!(angle_derivative < 0.);
    assert_eq!(radial_proof(0.3, 1.5 * tau), Err(Error::Unresolved));
    assert!(radial_proof(0., 0.).unwrap() > 6.);
}
#[test]
fn template_resource_cancel_and_precision_failures_never_publish_partial_data() {
    use std::cell::Cell;
    let options = CurveOptions::default();
    assert!(matches!(
        PeriodicAngular::prepare(false, CurveOptions { work: 1, ..options }, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        PeriodicAngular::prepare(
            false,
            CurveOptions {
                bytes: 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        PeriodicAngular::prepare(
            false,
            CurveOptions {
                maximum_spans_per_curve: 0,
                ..options
            },
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        PeriodicAngular::prepare(
            false,
            CurveOptions {
                maximum_spans_per_curve: 1,
                position_error_mm: 1e-30,
                ..options
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    let n = Cell::new(0);
    let good = PeriodicAngular::prepare(
        false,
        options,
        Some(&|| {
            n.set(n.get() + 1);
            false
        }),
    )
    .unwrap();
    for stop in [0, 1, n.get() / 2, n.get() - 1] {
        let at = Cell::new(0);
        assert!(matches!(
            PeriodicAngular::prepare(
                false,
                options,
                Some(&|| {
                    let k = at.get();
                    at.set(k + 1);
                    k == stop
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(
        PeriodicAngular::prepare(false, options, None)
            .unwrap()
            .poles(),
        good.poles()
    );
}
