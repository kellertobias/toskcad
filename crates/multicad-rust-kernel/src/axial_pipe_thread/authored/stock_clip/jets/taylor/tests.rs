use super::*;
fn plan() -> StockClipPlan {
    let mut spec = super::super::super::super::tests::source(1. / 32.).spec();
    spec.crest_radius_limit_mm = Some(0.405 * 25.4 / 2.);
    AuthoredThread::new(spec, 1e-7, 1e-5)
        .unwrap()
        .stock_clip_plan(Default::default(), None)
        .unwrap()
}
fn mid(b: [f64; 2]) -> f64 {
    b[0].midpoint(b[1])
}
fn cmul(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn cpow(x: f64, n: usize) -> [f64; 2] {
    (0..n).fold([1., 0.], |v, i| {
        cmul(v, [0., x]).map(|c| c / (i + 1) as f64)
    })
}
fn exp_coefficient(b: f64, c: f64, d: f64, i: usize, j: usize) -> [f64; 2] {
    // Independent combinatorial expansion of exp(iBs)exp(iCt)exp(iDst),
    // rather than the production first-derivative recurrence.
    (0..=i.min(j)).fold([0.; 2], |sum, k| {
        let term = cmul(cmul(cpow(b, i - k), cpow(c, j - k)), cpow(d, k));
        std::array::from_fn(|axis| sum[axis] + term[axis])
    })
}
#[test]
fn original_mixed_orders_include_bilinear_angular_coupling_and_match_independent_expansion() {
    let plan = plan();
    let saved = *plan.source();
    let contact = plan.contact_height_bounds_mm();
    let mut coupled = false;
    for patch in plan.patches() {
        let law = plan.patch_law(patch, Default::default(), None).unwrap();
        let corners = law.bilinear_fields_mm().map(|p| p.map(mid));
        for (u, v) in [(0., 0.), (0.13, 0.29), (0.5, 0.5), (0.91, 0.83), (1., 1.)] {
            let fields: [[f64; 4]; 3] = std::array::from_fn(|axis| {
                let c = corners.map(|p| p[axis]);
                let d = c[3] - c[2] - c[1] + c[0];
                [
                    c[0] + (c[1] - c[0]) * u + (c[2] - c[0]) * v + d * u * v,
                    c[1] - c[0] + d * v,
                    c[2] - c[0] + d * u,
                    d,
                ]
            });
            let omega = std::f64::consts::TAU / saved.spec().pitch_mm;
            let theta: [f64; 4] = std::array::from_fn(|k| omega * (fields[2][k] - fields[1][k]));
            coupled |= theta[3].abs() > 1e-8;
            let phase = [theta[0].cos(), theta[0].sin()];
            let actual = law
                .taylor([u, u], [v, v], Default::default(), None)
                .unwrap();
            let c = actual.coefficients_mm();
            assert!(actual.work_used() > 0);
            assert_eq!(actual.preparation_work_used(), law.work_used());
            assert!(actual.retained_budget_bytes() > law.retained_budget_bytes());
            for (i, row) in c.iter().enumerate() {
                for (j, coefficient) in row.iter().enumerate() {
                    let observed_l1 = coefficient.iter().map(|b| mid(*b).abs()).sum::<f64>();
                    assert!(observed_l1 <= actual.coefficient_l1_bounds_mm()[i][j] * (1. + 1e-8));
                    let mut xy = [0.; 2];
                    for (di, dj, radius) in [
                        (0, 0, fields[0][0]),
                        (1, 0, fields[0][1]),
                        (0, 1, fields[0][2]),
                        (1, 1, fields[0][3]),
                    ] {
                        if i >= di && j >= dj {
                            let e = exp_coefficient(theta[1], theta[2], theta[3], i - di, j - dj);
                            let term = cmul(phase, e);
                            for axis in 0..2 {
                                xy[axis] += radius * term[axis];
                            }
                        }
                    }
                    for axis in 0..2 {
                        let observed = mid(coefficient[axis]);
                        assert!(
                            (observed - xy[axis]).abs() <= 1e-8 * xy[axis].abs().max(1.),
                            "patch={patch:?} i={i} j={j} axis={axis} got={observed} want={}",
                            xy[axis]
                        );
                    }
                }
            }
            let jets = law.jets([u, u], [v, v], Default::default(), None).unwrap();
            for (index, expected, factor) in [
                ([1, 0], jets.du(), 1.),
                ([0, 1], jets.dv(), 1.),
                ([2, 0], jets.duu(), 2.),
                ([1, 1], jets.duv(), 1.),
                ([0, 2], jets.dvv(), 2.),
            ] {
                for axis in 0..3 {
                    let a = I::new(
                        c[index[0]][index[1]][axis][0],
                        c[index[0]][index[1]][axis][1],
                    ) * I::point(factor);
                    assert!(a.lo <= expected[axis][1] && expected[axis][0] <= a.hi);
                }
            }
            let whole = law
                .taylor([0., 1.], [0., 1.], Default::default(), None)
                .unwrap();
            for (i, row) in c.iter().enumerate() {
                for (j, coefficient) in row.iter().enumerate() {
                    let observed_l1 = coefficient.iter().map(|b| mid(*b).abs()).sum::<f64>();
                    assert!(observed_l1 <= whole.coefficient_l1_bounds_mm()[i][j] * (1. + 1e-8));
                    for (axis, component) in coefficient.iter().enumerate() {
                        let b = whole.coefficients_mm()[i][j][axis];
                        assert!(b[0] <= mid(*component) && mid(*component) <= b[1]);
                    }
                }
            }
        }
    }
    assert!(coupled);
    assert_eq!(*plan.source(), saved);
    assert_eq!(plan.contact_height_bounds_mm(), contact);
}
#[test]
fn original_patch_taylor_queries_honor_domains_actual_budgets_and_cancellation() {
    let plan = plan();
    let law = plan
        .patch_law(plan.patches()[0], Default::default(), None)
        .unwrap();
    let options = StockClipOptions::default();
    let result = law.taylor([0.2, 0.4], [0.6, 0.8], options, None).unwrap();
    for (u, v) in [
        ([f64::NAN, 0.], [0., 1.]),
        ([0., 1.], [0.8, 0.7]),
        ([-0.1, 0.], [0., 1.]),
    ] {
        assert!(matches!(
            law.taylor(u, v, options, None),
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
            law.taylor([0.2, 0.4], [0.6, 0.8], bad, None),
            Err(Error::Limit)
        ));
    }
    assert_eq!(
        law.taylor(
            [0.2, 0.4],
            [0.6, 0.8],
            StockClipOptions {
                maximum_work: result.work_used(),
                maximum_bytes: result.retained_budget_bytes(),
            },
            None
        )
        .unwrap()
        .work_used(),
        result.work_used()
    );
    assert!(matches!(
        law.taylor([0., 1.], [0., 1.], options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let calls = std::cell::Cell::new(0);
    let cancel = || {
        let n = calls.get() + 1;
        calls.set(n);
        n == 2
    };
    assert!(matches!(
        law.taylor([0., 1.], [0., 1.], options, Some(&cancel)),
        Err(Error::Cancelled)
    ));
}
