use super::*;

fn plan() -> StockClipPlan {
    let mut spec = super::super::super::tests::source(1. / 32.).spec();
    spec.crest_radius_limit_mm = Some(0.405 * 25.4 / 2.);
    AuthoredThread::new(spec, 1e-7, 1e-5)
        .unwrap()
        .stock_clip_plan(Default::default(), None)
        .unwrap()
}
fn midpoint(b: [f64; 2]) -> f64 {
    b[0].midpoint(b[1])
}
fn multiply(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

#[test]
fn all_twenty_original_edges_keep_symbolic_event_and_high_order_source_law() {
    let plan = plan();
    let original = *plan.source();
    let contact = plan.contact_height_bounds_mm();
    for curve in plan.curves() {
        let law = plan.curve_law(curve, Default::default(), None).unwrap();
        assert!(std::ptr::eq(law.plan(), &plan));
        assert_eq!(law.curve(), curve);
        assert!(law.work_used() > 0 && law.retained_budget_bytes() > 0);
        let endpoints = law.affine_fields_mm().map(|p| p.map(midpoint));
        let omega = std::f64::consts::TAU / original.spec().pitch_mm;
        let dr = endpoints[1][0] - endpoints[0][0];
        let dz = endpoints[1][2] - endpoints[0][2];
        let dq = endpoints[1][1] - endpoints[0][1];
        let angular_rate = omega * (dz - dq);
        for t in [0., 0.17, 0.5, 0.87, 1.] {
            let jet = law.taylor([t, t], Default::default(), None).unwrap();
            assert!(jet.work_used() > law.work_used());
            assert!(jet.retained_budget_bytes() > law.retained_budget_bytes());
            let c = jet.coefficients_mm();
            let radius = endpoints[0][0] + dr * t;
            let theta = omega * (endpoints[0][2] - endpoints[0][1] + (dz - dq) * t);
            let unit = [theta.cos(), theta.sin()];
            // Independent closed complex-exponential derivative formula,
            // including the affine-radius Leibniz term, through order12.
            let mut previous = [1., 0.];
            let mut power = [1., 0.];
            let mut factorial = 1.;
            for (n, coefficient) in c.into_iter().enumerate() {
                if n > 0 {
                    previous = power;
                    power = multiply(power, [0., angular_rate]);
                    factorial *= n as f64;
                }
                let mut expected = multiply(unit, power).map(|x| radius * x / factorial);
                if n > 0 {
                    let leibniz = multiply(unit, previous);
                    for j in 0..2 {
                        expected[j] += dr * n as f64 * leibniz[j] / factorial;
                    }
                }
                for j in 0..2 {
                    let error = (midpoint(coefficient[j]) - expected[j]).abs();
                    assert!(
                        error < 2e-11 * expected[j].abs().max(1.),
                        "{curve:?} {t} {n} {j}: {error}"
                    );
                }
                assert!(coefficient.into_iter().flatten().all(f64::is_finite));
            }
            let direct = plan
                .curve_bounds(curve, [t, t], Default::default(), None)
                .unwrap()
                .position_mm();
            for j in 0..3 {
                assert!(c[0][j][0] <= direct[j][1] && direct[j][0] <= c[0][j][1]);
            }
            assert!(c[1][2][0] <= dz && dz <= c[1][2][1]);
            assert!(c[2..].iter().all(|c| c[2][0] <= 0. && c[2][1] >= 0.));
        }
        let whole = law
            .taylor([0., 1.], Default::default(), None)
            .unwrap()
            .coefficients_mm();
        for t in [0., 0.25, 0.5, 0.75, 1.] {
            let point = law
                .taylor([t, t], Default::default(), None)
                .unwrap()
                .coefficients_mm();
            for n in 0..13 {
                for j in 0..3 {
                    assert!(whole[n][j][0] <= point[n][j][0] && point[n][j][1] <= whole[n][j][1]);
                }
            }
        }
    }
    assert_eq!(*plan.source(), original);
    assert_eq!(plan.contact_height_bounds_mm(), contact);
}

#[test]
fn source_edge_queries_charge_whole_request_and_refuse_invalid_domains_before_publication() {
    let plan = plan();
    let curve = plan.curves()[0];
    let options = StockClipOptions::default();
    let law = plan.curve_law(curve, options, None).unwrap();
    let jet = law.taylor([0., 1.], options, None).unwrap();
    for domain in [[-0.001, 0.], [0., 1.001], [0.8, 0.2], [f64::NAN, 0.]] {
        assert_eq!(
            law.taylor(domain, options, None).unwrap_err(),
            Error::Invalid
        );
    }
    for maximum_work in [0, options.maximum_work + 1] {
        let o = StockClipOptions {
            maximum_work,
            ..options
        };
        assert!(matches!(
            plan.curve_law(curve, o, None),
            Err(Error::Invalid)
        ));
    }
    for o in [
        StockClipOptions {
            maximum_work: 1,
            ..options
        },
        StockClipOptions {
            maximum_bytes: 1,
            ..options
        },
        StockClipOptions {
            maximum_work: law.work_used() - 1,
            ..options
        },
    ] {
        assert!(matches!(plan.curve_law(curve, o, None), Err(Error::Limit)));
    }
    let exact = plan
        .curve_law(
            curve,
            StockClipOptions {
                maximum_work: law.work_used(),
                ..options
            },
            None,
        )
        .unwrap();
    assert_eq!(exact.work_used(), law.work_used());
    assert_eq!(
        law.taylor(
            [0., 1.],
            StockClipOptions {
                maximum_work: jet.work_used() - 1,
                ..options
            },
            None
        )
        .unwrap_err(),
        Error::Limit
    );
    assert_eq!(
        law.taylor(
            [0., 1.],
            StockClipOptions {
                maximum_work: jet.work_used(),
                ..options
            },
            None
        )
        .unwrap()
        .work_used(),
        jet.work_used()
    );
    assert!(matches!(
        plan.curve_law(curve, options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert_eq!(
        law.taylor([0., 1.], options, Some(&|| true)).unwrap_err(),
        Error::Cancelled
    );
    // Final publication poll, after actual recurrence work, remains observable.
    let polls = std::cell::Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    law.taylor([0., 1.], options, Some(&count)).unwrap();
    let seen = std::cell::Cell::new(0usize);
    let cancel = || {
        let n = seen.get();
        seen.set(n + 1);
        n + 1 == polls.get()
    };
    assert_eq!(
        law.taylor([0., 1.], options, Some(&cancel)).unwrap_err(),
        Error::Cancelled
    );
    assert_eq!(
        law.taylor([0., 1.], options, None)
            .unwrap()
            .coefficients_mm(),
        jet.coefficients_mm()
    );
}
