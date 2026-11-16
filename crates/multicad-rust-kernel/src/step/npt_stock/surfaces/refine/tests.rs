use super::*;
#[test]
fn actual_source_curve_refinement_shares_dyadic_owners_and_retains_c0_c1() {
    let source = super::super::super::tests::source(false);
    let options = CurveOptions::default();
    let mut plan = CurvePlan::new(&source, options, None).unwrap();
    let mut c = plan.curves.remove(0);
    let count = c.spans.len() * 2;
    let mut memory = Memory::new(plan.bytes, options.maximum_bytes).unwrap();
    let mut work = Work::new(options.maximum_work, None).unwrap();
    curve(&mut c, count, &mut memory, &mut work).unwrap();
    assert_eq!(c.spans.len(), count);
    for (i, span) in c.spans.iter().enumerate() {
        assert_eq!(
            span.owners,
            [
                curves::owner(c.source, i, count),
                curves::owner(c.source, i + 1, count)
            ]
        );
        let law = source
            .curve_law(c.source, Default::default(), None)
            .unwrap();
        for local in [0., 0.13, 0.5, 0.91, 1.] {
            let t = span.domain[0] + local * (span.domain[1] - span.domain[0]);
            let expected = law
                .taylor([t, t], Default::default(), None)
                .unwrap()
                .coefficients_mm();
            let actual = super::super::super::tests::evaluate(span.poles, local);
            let derivative: [[f64; 3]; 11] = std::array::from_fn(|i| {
                std::array::from_fn(|axis| {
                    11. * (span.poles[i + 1][axis] - span.poles[i][axis])
                        / (span.domain[1] - span.domain[0])
                })
            });
            let actual_d = super::super::super::tests::evaluate(derivative, local);
            for axis in 0..3 {
                assert!(
                    (actual[axis] - expected[0][axis][0].midpoint(expected[0][axis][1])).abs()
                        <= span.position_error + 1e-10
                );
                assert!(
                    (actual_d[axis] - expected[1][axis][0].midpoint(expected[1][axis][1])).abs()
                        <= span.derivative_error + 1e-10
                );
            }
        }
    }
    for pair in c.spans.windows(2) {
        assert_eq!(pair[0].owners[1], pair[1].owners[0]);
        assert_eq!(pair[0].poles[11], pair[1].poles[0]);
    }
    let mut work = Work::new(1, None).unwrap();
    assert!(matches!(
        curve(&mut c, count * 2, &mut memory, &mut work),
        Err(Error::Limit)
    ));
}
