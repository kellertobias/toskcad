use super::*;
#[test]
fn actual_source_box_queries_distinguish_material_void_and_boundary() {
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let point = |p: [f64; 3]| p.map(|x| [x, x]);
    for p in [[0., 0., 1.], [1., -1., 3.], [2., 2., 4.]] {
        assert_eq!(
            atlas
                .classify_box(point(p), Default::default(), None)
                .unwrap()
                .classification(),
            StockClipClassification::Inside
        );
    }
    for p in [[0., 0., -1.], [0., 0., 10.], [6., 0., 1.], [-6., 0., 1.]] {
        assert_eq!(
            atlas
                .classify_box(point(p), Default::default(), None)
                .unwrap()
                .classification(),
            StockClipClassification::Outside
        );
    }
    assert_eq!(
        atlas
            .classify_box(point([0., 0., 0.]), Default::default(), None)
            .unwrap()
            .classification(),
        StockClipClassification::Undetermined
    );
    for patch in atlas.patches() {
        for (u, v) in [(0.2, 0.2), (0.5, 0.5), (0.8, 0.8)] {
            let p = atlas
                .point_bounds(patch, [u, u], [v, v], Default::default(), None)
                .unwrap()
                .position_mm()
                .map(|v| I::new(v[0], v[1]).midpoint());
            let inside = [p[0] * 0.999, p[1] * 0.999, p[2]];
            let outside = [p[0] * 1.001, p[1] * 1.001, p[2]];
            assert_eq!(
                atlas
                    .classify_box(point(inside), Default::default(), None)
                    .unwrap()
                    .classification(),
                StockClipClassification::Inside
            );
            assert_eq!(
                atlas
                    .classify_box(point(outside), Default::default(), None)
                    .unwrap()
                    .classification(),
                StockClipClassification::Outside
            );
        }
    }
}
#[test]
fn source_queries_keep_budget_cancel_numeric_and_source_contracts() {
    let source = super::super::tests::canonical(true, 0.1);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    assert_eq!(atlas.bounds_mm()[2], [0., source.spec().length_mm]);
    let b = [[0., 0.], [0., 0.], [0.1, 0.2]];
    let earned = atlas.classify_box(b, Default::default(), None).unwrap();
    assert!(earned.work_used() > 0);
    assert!(earned.retained_budget_bytes() > 0);
    for limits in [
        StockClipOptions {
            maximum_work: 0,
            ..Default::default()
        },
        StockClipOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            atlas.classify_box(b, limits, None),
            Err(Error::Invalid)
        ));
    }
    for limits in [
        StockClipOptions {
            maximum_work: 1,
            ..Default::default()
        },
        StockClipOptions {
            maximum_bytes: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            atlas.classify_box(b, limits, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        atlas.classify_box(b, Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        atlas.classify_box([[1e8, 1e8], [0., 0.], [0., 0.]], Default::default(), None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        atlas.classify_box([[1., 0.], [0., 0.], [0., 0.]], Default::default(), None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        atlas.classify_box(
            [[f64::NAN, 0.], [0., 0.], [0., 0.]],
            Default::default(),
            None
        ),
        Err(Error::Invalid)
    ));
    let count = std::cell::Cell::new(0);
    let cancel = || {
        let n = count.get();
        count.set(n + 1);
        n > 5
    };
    let r = atlas.stock_radius_limit_mm() * 0.99;
    assert!(matches!(
        atlas.classify_box(
            [[r, r], [0., 0.], [0.3, 0.3]],
            Default::default(),
            Some(&cancel)
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(*atlas.source(), source);
    assert_eq!(
        atlas
            .classify_box(b, Default::default(), None)
            .unwrap()
            .classification(),
        earned.classification()
    );
}
#[test]
fn angle_quadrants_axis_crossings_and_seam_keep_true_directions() {
    for (x, y) in [
        (1., 2.),
        (1., -2.),
        (-1., 2.),
        (-1., -2.),
        (0., 2.),
        (0., -2.),
        (0., 0.),
    ] {
        let r = angle(I::point(x), I::point(y));
        let a = y.portable_atan2(x);
        assert!(r.lo <= a && a <= r.hi);
    }
    let r = angle(I::new(-1., 1.), I::new(-1., 1.));
    assert!(r.lo < -std::f64::consts::PI && r.hi > std::f64::consts::PI);
}
