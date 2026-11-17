use super::*;
fn fresh_default_budget() -> Budget<'static> {
    // The same valid initial state as the public calculation: no successful
    // debit, observed work or earned receipt is fabricated by this fixture.
    Budget {
        initial: MAX_WORK,
        remaining: MAX_WORK,
        peak: 0,
        evaluated: 0,
        interrupt: None,
        work_by_stage: [0; 10],
        attempted_stage: None,
        attempted_work: 0,
        attempted_bytes: 0,
        axis_samples: [0; 2],
        axis_unresolved: [0; 2],
        axis_remainder_max: [[0.; 4]; 2],
        last_aggregate_widths: None,
    }
}
#[test]
fn actual_distinct_plane_pairs_refuse_disjoint_private_numeric_arguments_and_retry() {
    let plane = |scale| {
        let source = RetainedTensor11::new(
            [1, 1],
            (0..12)
                .flat_map(|i| (0..12).map(move |j| [scale * i as f64, j as f64, 3.]))
                .collect(),
        )
        .unwrap();
        let geometry = source.cell([0, 0]).unwrap();
        let normal = Arc::new(NormalPolynomial11::new(&geometry).unwrap());
        let squared = SquaredNormal11::new(normal).unwrap();
        (geometry, squared)
    };
    let (a, normal_a) = plane(1.);
    let (b, normal_b) = plane(2.);
    assert!(!std::ptr::eq(normal_a.source(), normal_b.source()));
    let options = SurfaceOptions11::default();
    let own_a = enclosure(&a, &normal_a, [0.; 3], options, &mut fresh_default_budget())
        .unwrap()
        .value
        .map(pair);
    let own_b = enclosure(&b, &normal_b, [0.; 3], options, &mut fresh_default_budget())
        .unwrap()
        .value
        .map(pair);
    assert!(own_a[0][0] <= 121. && 121. <= own_a[0][1]);
    assert!(own_b[0][0] <= 242. && 242. <= own_b[0][1]);
    // TWO different genuine numeric sources, not contradictory enclosures
    // claimed for one source. The public Cell route always pairs its own
    // position and normal. This specific disjoint pair must refuse; this is
    // not a universal wrong-owner detector or a public-reach claim.
    assert_eq!(
        enclosure(&a, &normal_b, [0.; 3], options, &mut fresh_default_budget())
            .map(|x| x.value.map(pair)),
        Err(Error::Unresolved)
    );
    assert_eq!(
        enclosure(&a, &normal_a, [0.; 3], options, &mut fresh_default_budget())
            .unwrap()
            .value
            .map(pair),
        own_a
    );
    assert_eq!(
        enclosure(&b, &normal_b, [0.; 3], options, &mut fresh_default_budget())
            .unwrap()
            .value
            .map(pair),
        own_b
    );
}
#[test]
fn actual_d11_point_and_squared_normal_helper_enclose_independent_child_geometry() {
    let source = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| (0..12).map(move |j| [i as f64, j as f64, (i * j) as f64]))
            .collect(),
    )
    .unwrap();
    let root = source.cell([0, 0]).unwrap();
    let [_, child] = root.subdivide(0).unwrap();
    for (cell, original_u, jacobian_squared) in [(&root, 0.25, 1.), (&child, 0.625, 0.25)] {
        let (p, q) = point_at(cell, [I::point(0.25), I::point(0.5)]).unwrap();
        let expected = (14641. + 1771561. * (original_u * original_u + 0.25)) * jacobian_squared;
        assert!(q.lo <= expected && expected <= q.hi);
        for (bound, value) in p
            .into_iter()
            .zip([11. * original_u, 5.5, 121. * original_u * 0.5])
        {
            assert!(bound.lo <= value && value <= bound.hi);
        }
    }
    let (bounds, q) = point_at(&root, [I::new(0.25, 0.75), I::new(0.5, 1.)]).unwrap();
    for u in [0.25, 0.75] {
        for v in [0.5, 1.] {
            for (bound, value) in bounds.into_iter().zip([11. * u, 11. * v, 121. * u * v]) {
                assert!(bound.lo <= value && value <= bound.hi);
            }
            let expected = 14641. + 1771561. * (u * u + v * v);
            assert!(q.lo <= expected && expected <= q.hi);
        }
    }
}
#[test]
fn actual_distinct_plane_area_enclosures_refuse_contradictory_helper_intersection() {
    // An independent interval-helper contract. The public same-source path
    // never pairs these two different actual geometric surfaces.
    let plane = |scale| {
        RetainedTensor11::new(
            [1, 1],
            (0..12)
                .flat_map(|i| (0..12).map(move |j| [scale * i as f64, j as f64, 3.]))
                .collect(),
        )
        .unwrap()
    };
    let a = plane(1.);
    let b = plane(2.);
    let first = a.surface_properties(Default::default(), None).unwrap();
    let second = b.surface_properties(Default::default(), None).unwrap();
    assert!(first.area_mm2[1] < second.area_mm2[0]);
    let x = I::new(first.area_mm2[0], first.area_mm2[1]);
    let y = I::new(second.area_mm2[0], second.area_mm2[1]);
    assert_eq!(intersect(x, y).map(pair), Err(Error::Unresolved));
    assert_eq!(intersect(x, x).map(pair).unwrap(), first.area_mm2);
    assert_eq!(intersect(y, y).map(pair).unwrap(), second.area_mm2);
    assert_eq!(
        a.surface_properties(Default::default(), None).unwrap(),
        first
    );
    assert_eq!(
        b.surface_properties(Default::default(), None).unwrap(),
        second
    );
}
