use super::*;
#[test]
fn actual_canonical_and_nonbinary_tensors_have_positive_literal_projection() {
    for kind in 0..3 {
        let body = if kind == 0 {
            super::super::tests::metric_graph()
        } else {
            super::super::tests::graph(kind == 2, 1)
        };
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let p = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let proof = p
            .projection_bounds(16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        assert!(proof.minimum_signed_jacobian > 0.);
        assert!(proof.component_radius_mm[0][0] > proof.component_radius_mm[1][1]);
        assert!(proof.certified_cells >= p.patches().len() * 2);
        assert!(proof.work_used <= 16_000_000);
        eprintln!(
            "literal projection {kind}: {:?}, jac{} cells{} work{}",
            proof.component_radius_mm,
            proof.minimum_signed_jacobian,
            proof.certified_cells,
            proof.work_used
        );
    }
}
#[test]
fn actual_projection_certification_respects_global_resource_and_cancel_limits() {
    let body = super::super::tests::metric_graph();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    for bytes in [0, 64 * 1024 * 1024 + 1] {
        assert!(matches!(
            s.projection_bounds(16_000_000, bytes, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        s.projection_bounds(1, 64 * 1024 * 1024, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        s.projection_bounds(16_000_000, s.retained_bytes(), None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        s.projection_bounds(16_000_000, 64 * 1024 * 1024, Some(&|| true)),
        Err(Error::Cancelled)
    ));
}

#[test]
fn literal_control_hulls_cannot_escape_actual_cap_planes() {
    let mut w = Work::new(16_000_000, None).unwrap();
    height_controls(&[[0., 0., 0.], [0., 0., 4.]], 4., &mut w).unwrap();
    for z in [-f64::MIN_POSITIVE, 4_f64.next_up(), f64::NAN] {
        assert!(matches!(
            height_controls(&[[0., 0., z]], 4., &mut w),
            Err(Error::Unresolved)
        ));
    }
}
#[test]
fn actual_bernstein_positive_and_tangent_controls_subdivide_or_refuse() {
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut stack = Vec::with_capacity(3 * DEPTH + 1);
    let p = P::controls(&[I::point(1.), I::point(-0.1), I::point(1.)]);
    let (bound, cells) = positive([p; 3], &mut stack, &mut w).unwrap();
    assert!(bound.lo > 0. && cells > 1);
    for values in [[-1., -1., -1.], [1., -1., 1.], [f64::NAN, 1., 1.]] {
        stack.clear();
        let p = P::controls(&values.map(I::point));
        assert!(matches!(
            positive([p; 3], &mut stack, &mut w),
            Err(Error::Unresolved)
        ));
    }
    let table = polynomial::test_choices();
    for (n, row) in table.iter().enumerate() {
        let mut exact = 1u64;
        for (k, &value) in row.iter().enumerate().take(n + 1) {
            assert_eq!(value, exact);
            exact = exact * (n - k) as u64 / (k + 1) as u64;
        }
    }
}

#[test]
fn finite_ruled_projection_honors_every_actual_work_frontier() {
    let p = |values: &[f64]| P::controls(&values.iter().map(|v| I::point(*v)).collect::<Vec<_>>());
    // Actual noncollapsed affine patch S=[1+u,t,u], with negative
    // cylindrical projection determinant -(1+u) on the entire square.
    let a = [p(&[1., 1.]), p(&[0., 1.]), p(&[0., 0.])];
    let b = [p(&[2., 2.]), p(&[0., 1.]), p(&[1., 1.])];
    let mut complete = Work::new(16_000_000, None).unwrap();
    let expected = jacobian([a, b], false, &mut complete).unwrap();
    assert!(expected.iter().all(|p| p.hull().hi < 0.));
    for limit in 1..complete.used {
        let mut short = Work::new(limit, None).unwrap();
        assert!(matches!(
            jacobian([a, b], false, &mut short),
            Err(Error::Limit)
        ));
        assert!(short.used <= limit);
    }
    let mut retry = Work::new(complete.used, None).unwrap();
    let actual = jacobian([a, b], false, &mut retry).unwrap();
    for (a, b) in actual.iter().zip(expected) {
        assert_eq!([a.hull().lo, a.hull().hi], [b.hull().lo, b.hull().hi]);
    }
}
