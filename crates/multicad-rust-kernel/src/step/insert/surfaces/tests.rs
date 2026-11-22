use super::*;
fn domain(body: &Boundary) -> DomainPlan<'_> {
    DomainPlan::prepare(body, DomainOptions::default(), None).unwrap()
}
#[test]
fn actual_shared_rows_have_identical_tensor_boundaries_and_owned_event_vertices() {
    for (left, starts) in [(false, 1), (true, 1), (false, 2), (true, 4)] {
        let body = super::super::tests::graph(left, starts);
        let d = domain(&body);
        let surface = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        assert!(std::ptr::eq(surface.domain(), &d));
        assert_eq!(surface.patches().len(), d.slabs().len());
        assert!(surface.rows().len() < surface.patches().len() * 2);
        for row in surface.rows() {
            let ends = row.endpoints();
            assert_eq!(row.poles()[0], surface.vertices()[ends[0]].point_mm());
            assert_eq!(
                *row.poles().last().unwrap(),
                surface.vertices()[ends[1]].point_mm()
            );
            assert!(row.position_error_mm() <= 1e-7);
            assert!(row.derivative_error_mm() <= 1e-6);
            for (i, station) in [row.key().interval, row.key().interval + 1]
                .into_iter()
                .enumerate()
            {
                if let Some(q) = d.stations(row.key().curve.owner())[station as usize].cardinal() {
                    let p = surface.vertices()[ends[i]].point_mm();
                    assert_eq!(p[if q % 2 == 0 { 1 } else { 0 }], 0.);
                }
            }
            if let Some(cap) = body
                .curve(row.key().curve)
                .unwrap()
                .faces()
                .into_iter()
                .find(|&f| f < 2)
            {
                let z = if cap == 0 {
                    0.
                } else {
                    body.source().full_length_mm()
                };
                assert!(row.poles().iter().all(|p| p[2] == z));
            }
        }
        for pair in surface.rows().windows(2) {
            if pair[0].key().curve == pair[1].key().curve
                && pair[0].key().interval + 1 == pair[1].key().interval
            {
                assert_eq!(pair[0].endpoints()[1], pair[1].endpoints()[0]);
            }
        }
        for patch in surface.patches() {
            let slab = &d.slabs()[patch.slab];
            let [a, b] = patch.rows.map(|i| &surface.rows()[i]);
            assert_eq!(a.knots(), b.knots());
            assert_eq!(
                patch.same_sense,
                body.source_boundary(slab.side).source().spec().left_hand
                    != (slab.side == Side::Inner)
            );
            for i in 0..2 {
                if slab.collapsed[i] {
                    assert_eq!(a.endpoints()[i], b.endpoints()[i]);
                }
            }
        }
        assert!(
            surface
                .vertices()
                .iter()
                .any(|p| matches!(p.owner(), VertexOwner::Station { .. }))
        );
        assert!(
            surface
                .vertices()
                .iter()
                .any(|p| matches!(p.owner(), VertexOwner::Original(_)))
        );
        assert!(surface.work_used() <= 16_000_000);
        assert!(surface.retained_bytes() <= 64 * 1024 * 1024);
        eprintln!(
            "shared surfaces {left}/{starts}: {} rows {} patches {} vertices {} work",
            surface.rows().len(),
            surface.patches().len(),
            surface.vertices().len(),
            surface.work_used()
        );
    }
}
#[test]
fn actual_row_refinement_preserves_common_grids_and_failure_budgets() {
    let body = super::super::tests::graph(false, 1);
    let d = domain(&body);
    let tight = CurveOptions {
        derivative_error_mm: 1e-10,
        maximum_spans_per_curve: 32,
        ..CurveOptions::default()
    };
    let plan = SurfacePlan::prepare(&d, tight, None).unwrap();
    assert!(plan.rows().iter().any(|r| r.knots().len() > 2));
    for patch in plan.patches() {
        assert_eq!(
            plan.rows()[patch.rows[0]].knots(),
            plan.rows()[patch.rows[1]].knots()
        );
    }
    assert!(matches!(
        SurfacePlan::prepare(
            &d,
            CurveOptions {
                maximum_spans_per_curve: 1,
                ..tight
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        SurfacePlan::prepare(
            &d,
            CurveOptions {
                work: d.work_used(),
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        SurfacePlan::prepare(
            &d,
            CurveOptions {
                bytes: d.retained_bytes(),
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        SurfacePlan::prepare(
            &d,
            CurveOptions {
                position_error_mm: 0.,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        SurfacePlan::prepare(&d, CurveOptions::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let polls = std::cell::Cell::new(0);
    let stop = || {
        polls.set(polls.get() + 1);
        polls.get() > 1000
    };
    assert!(matches!(
        SurfacePlan::prepare(&d, CurveOptions::default(), Some(&stop)),
        Err(Error::Cancelled)
    ));
}
#[test]
fn cancellation_during_actual_source_restriction_never_publishes_rows() {
    let body = super::super::tests::graph(false, 1);
    let d = domain(&body);
    let polls = std::cell::Cell::new(0);
    let stop = || {
        polls.set(polls.get() + 1);
        polls.get() > 25_000
    };
    assert!(matches!(
        SurfacePlan::prepare(&d, CurveOptions::default(), Some(&stop)),
        Err(Error::Cancelled)
    ));
    assert_eq!(polls.get(), 25_001);
}

#[test]
fn actual_derived_candidate_rejects_lost_collapse_identity_and_boundary_corruption() {
    let body = super::super::tests::graph(false, 1);
    let d = domain(&body);
    let mut plan = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let patch = *plan
        .patches()
        .iter()
        .find(|p| d.slabs()[p.slab].collapsed.contains(&true))
        .unwrap();
    let slab = &d.slabs()[patch.slab];
    let [a, b] = patch.rows;
    let old = plan.rows[b].ends;
    let end = usize::from(!slab.collapsed[0]);
    let duplicate = plan.vertices.len();
    let old_vertex = &plan.vertices[old[end]];
    plan.vertices.push(SurfaceVertex {
        owner: old_vertex.owner,
        point: old_vertex.point,
    });
    plan.rows[b].ends[end] = duplicate;
    assert!(matches!(
        verify_boundaries(slab, &plan.rows[a], &plan.rows[b], &plan.vertices),
        Err(Error::Unresolved)
    ));
    plan.rows[b].ends = old;
    let original = plan.rows[b].fit.poles[0];
    plan.rows[b].fit.poles[0][0] += 0.001;
    assert!(matches!(
        verify_boundaries(slab, &plan.rows[a], &plan.rows[b], &plan.vertices),
        Err(Error::Unresolved)
    ));
    plan.rows[b].fit.poles[0] = original;
    plan.rows[b].fit.knots[0] = 0.001;
    assert!(matches!(
        verify_boundaries(slab, &plan.rows[a], &plan.rows[b], &plan.vertices),
        Err(Error::Unresolved)
    ));
}
#[test]
fn source_rational_coefficients_earn_float_enclosures_instead_of_trusting_proposals() {
    use num_bigint::BigInt;
    use num_rational::BigRational as Q;
    let mut w = Work::new(16_000_000, None).unwrap();
    assert_eq!(N::integer(1).enclosure(&mut w).unwrap().lo, 1.);
    for numerator in [1, -1] {
        let n = N(Q::new(numerator.into(), 3.into()));
        let i = n.enclosure(&mut w).unwrap();
        assert!(Q::from_float(i.lo).unwrap() <= n.0 && Q::from_float(i.hi).unwrap() >= n.0);
        assert!(matches!(
            n.enclose_proposal(0., &mut w),
            Err(Error::Unresolved)
        ));
    }
    let huge = N(Q::from_integer(BigInt::from(1) << 2048));
    assert!(matches!(huge.enclosure(&mut w), Err(Error::Unresolved)));
    assert!(matches!(
        N::integer(1).enclose_proposal(f64::INFINITY, &mut w),
        Err(Error::Unresolved)
    ));
    for sign in [1, -1] {
        let n = N(
            (Q::from_float(f64::MAX).unwrap() + Q::from_integer(1.into()))
                * Q::from_integer(sign.into()),
        );
        assert!(matches!(n.enclosure(&mut w), Err(Error::Unresolved)));
    }
}

#[test]
fn canonical_metric_insert_builds_real_local_coupled_rows() {
    let body = super::super::tests::metric_graph();
    let d = domain(&body);
    let p = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    assert_eq!(p.patches().len(), d.slabs().len());
    assert!(
        p.rows()
            .iter()
            .all(|r| r.position_error_mm() <= 1e-7 && r.derivative_error_mm() <= 1e-6)
    );
    eprintln!(
        "CANONICAL_LOCAL rows{} patches{} vertices{} work{} bytes{}",
        p.rows().len(),
        p.patches().len(),
        p.vertices().len(),
        p.work_used(),
        p.retained_bytes()
    );
}
