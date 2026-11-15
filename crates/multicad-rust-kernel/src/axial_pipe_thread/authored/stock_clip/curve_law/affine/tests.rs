use super::*;
fn plan() -> StockClipPlan {
    let mut spec = super::super::super::super::tests::source(1. / 32.).spec();
    spec.crest_radius_limit_mm = Some(0.405 * 25.4 / 2.);
    AuthoredThread::new(spec, 1e-7, 1e-5)
        .unwrap()
        .stock_clip_plan(Default::default(), None)
        .unwrap()
}
#[test]
fn all_original_affine_edges_match_general_series_and_source_jets_without_repaying_preparation() {
    let source = plan();
    let original = *source.source();
    for curve in source.curves() {
        let law = source.curve_law(curve, Default::default(), None).unwrap();
        for parameter in [
            [0., 0.],
            [0.13, 0.13],
            [0.5, 0.5],
            [0.91, 0.91],
            [1., 1.],
            [0., 1.],
        ] {
            let fast = law
                .affine_taylor(parameter, Default::default(), None)
                .unwrap();
            let general = law.taylor(parameter, Default::default(), None).unwrap();
            assert_eq!(fast.work_used(), 1025);
            assert_eq!(fast.preparation_work_used(), law.work_used());
            assert!(fast.retained_budget_bytes() > law.retained_budget_bytes());
            let a = fast.coefficients_mm();
            let b = general.coefficients_mm();
            for n in 0..13 {
                for axis in 0..3 {
                    assert!(a[n][axis][0] <= b[n][axis][1] && b[n][axis][0] <= a[n][axis][1]);
                }
            }
            if parameter[0] == parameter[1] {
                for (n, c) in a.iter().enumerate() {
                    let norm = c.iter().map(|b| b[0].midpoint(b[1]).abs()).sum::<f64>();
                    assert!(norm <= fast.coefficient_l1_bounds_mm()[n] * (1. + 1e-10));
                }
            }
        }
    }
    assert_eq!(*source.source(), original);
}
#[test]
fn affine_queries_preserve_domain_errors_real_work_byte_frontiers_and_cancellation() {
    let source = plan();
    let law = source
        .curve_law(source.curves()[0], Default::default(), None)
        .unwrap();
    let options = StockClipOptions::default();
    let result = law.affine_taylor([0.2, 0.4], options, None).unwrap();
    for parameter in [[f64::NAN, 0.], [-0.1, 0.], [0.8, 0.7]] {
        assert!(matches!(
            law.affine_taylor(parameter, options, None),
            Err(Error::Invalid)
        ));
    }
    for bad in [
        StockClipOptions {
            maximum_work: result.work_used() - 1,
            ..options
        },
        StockClipOptions {
            maximum_bytes: result.retained_budget_bytes() - 1,
            ..options
        },
    ] {
        assert!(matches!(
            law.affine_taylor([0.2, 0.4], bad, None),
            Err(Error::Limit)
        ));
    }
    assert_eq!(
        law.affine_taylor(
            [0.2, 0.4],
            StockClipOptions {
                maximum_work: result.work_used(),
                maximum_bytes: result.retained_budget_bytes()
            },
            None
        )
        .unwrap()
        .work_used(),
        result.work_used()
    );
    assert!(matches!(
        law.affine_taylor([0., 1.], options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let seen = std::cell::Cell::new(0);
    let cancel = || {
        seen.set(seen.get() + 1);
        seen.get() == 2
    };
    assert!(matches!(
        law.affine_taylor([0., 1.], options, Some(&cancel)),
        Err(Error::Cancelled)
    ));
}
