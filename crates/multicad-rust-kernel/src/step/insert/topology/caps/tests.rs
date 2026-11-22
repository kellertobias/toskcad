use super::*;
use multicad_precision::PortableMath as _;
#[test]
fn actual_literal_cap_rims_are_coplanar_and_strict_one_turn_jordan_curves() {
    for kind in 0..4 {
        let body = if kind == 0 {
            super::super::super::tests::metric_graph()
        } else {
            super::super::super::tests::graph(kind == 2, if kind == 3 { 2 } else { 1 })
        };
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let proof = t
            .cap_regularity_bounds(16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        assert!(proof.minimum_signed_angular_derivative > 0.);
        assert!(proof.certified_segments > 0);
        assert!(proof.work_used <= 16_000_000);
        eprintln!(
            "CAP_JORDAN {kind}: angular{} segments{} work{}",
            proof.minimum_signed_angular_derivative, proof.certified_segments, proof.work_used
        );
    }
}
#[test]
fn winding_one_does_not_certify_a_backtracking_literal_rim() {
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut points = [0., 1., 0.5, 2., 4.]
        .map(|t: f64| [t.portable_cos(), t.portable_sin(), 0.])
        .to_vec();
    points.push(points[0]);
    let mut winding = 0;
    for p in points.windows(2) {
        winding += super::super::winding::crossing(p[0], p[1], &mut w).unwrap();
    }
    assert_eq!(winding, 1);
    let back = &points[1..=2];
    let x = P::controls(&back.iter().map(|p| literal(p[0])).collect::<Vec<_>>());
    let y = P::controls(&back.iter().map(|p| literal(p[1])).collect::<Vec<_>>());
    let p = x
        .product(y.derivative(&mut w).unwrap(), &mut w)
        .unwrap()
        .sub(
            y.product(x.derivative(&mut w).unwrap(), &mut w).unwrap(),
            &mut w,
        )
        .unwrap();
    let mut stack = Vec::with_capacity(DEPTH + 1);
    assert!(matches!(
        positive(p, &mut stack, &mut w),
        Err(Error::Unresolved)
    ));
}
#[test]
fn actual_cap_regularization_resources_and_primitive_domains_fail_closed() {
    let body = super::super::super::tests::metric_graph();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    for bytes in [0, 64 * 1024 * 1024 + 1] {
        assert!(matches!(
            t.cap_regularity_bounds(16_000_000, bytes, None),
            Err(Error::Invalid)
        ));
    }
    for bytes in [
        t.retained_bytes(),
        t.retained_bytes() + SCRATCH,
        t.retained_bytes() + SCRATCH + 1024,
    ] {
        assert!(matches!(
            t.cap_regularity_bounds(16_000_000, bytes, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        t.cap_regularity_bounds(1, 64 * 1024 * 1024, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        t.cap_regularity_bounds(16_000_000, 64 * 1024 * 1024, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let mut w = Work::new(16_000_000, None).unwrap();
    assert!(matches!(
        coplanar(&[[0., 0., 1.]], 0., &mut w),
        Err(Error::Unresolved)
    ));
    let mut stack = Vec::with_capacity(DEPTH + 1);
    let p = P::controls(&[I::point(1.), I::point(-0.1), I::point(1.)]);
    assert!(positive(p, &mut stack, &mut w).unwrap() > 0.);
    for values in [[-1., -1., -1.], [1., -1., 1.], [f64::NAN, 1., 1.]] {
        stack.clear();
        let p = P::controls(&values.map(I::point));
        assert!(matches!(
            positive(p, &mut stack, &mut w),
            Err(Error::Unresolved)
        ));
    }
}
