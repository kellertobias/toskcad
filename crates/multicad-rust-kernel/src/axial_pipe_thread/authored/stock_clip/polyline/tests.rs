use super::*;
#[test]
fn actual_source_edges_have_shared_endpoint_owned_bounded_chords() {
    for quarter in [false, true] {
        let source = super::super::tests::canonical(quarter, 1.);
        let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
        for curve in atlas.curves() {
            let p = atlas
                .curve_polyline(curve, 0.01, 8192, Default::default(), None)
                .unwrap();
            assert_eq!(p.curve(), curve);
            assert!(p.source_chord_error_mm() <= 0.01);
            assert!(p.work_used() <= 16_000_000);
            assert!(p.retained_budget_bytes() <= 64 * 1024 * 1024);
            for (i, key) in curve.endpoints().into_iter().enumerate() {
                let expected =
                    chosen(atlas.vertex_bounds(key, Default::default(), None).unwrap()).p;
                let actual = if i == 0 {
                    p.points_mm()[0]
                } else {
                    *p.points_mm().last().unwrap()
                };
                assert_eq!(actual, expected);
            }
        }
    }
}
#[test]
fn actual_curve_chords_resources_cancellation_and_tight_refusal_are_transactional() {
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let curve = atlas.curves()[0];
    let good = atlas
        .curve_polyline(curve, 0.01, 8192, Default::default(), None)
        .unwrap();
    for (error, points) in [
        (0., 8192),
        (-1., 8192),
        (f64::NAN, 8192),
        (0.01, 1),
        (0.01, 131073),
    ] {
        assert!(matches!(
            atlas.curve_polyline(curve, error, points, Default::default(), None),
            Err(Error::Invalid)
        ));
    }
    for limits in [
        StockClipOptions {
            maximum_work: good.work_used() - 1,
            ..Default::default()
        },
        StockClipOptions {
            maximum_bytes: good.retained_budget_bytes() - 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            atlas.curve_polyline(curve, 0.01, 8192, limits, None),
            Err(Error::Limit)
        ));
    }
    let limits = StockClipOptions {
        maximum_work: good.work_used(),
        maximum_bytes: good.retained_budget_bytes(),
    };
    assert_eq!(
        atlas
            .curve_polyline(curve, 0.01, 8192, limits, None)
            .unwrap()
            .points_mm(),
        good.points_mm()
    );
    assert!(matches!(
        atlas.curve_polyline(curve, 0.01, 2, Default::default(), None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        atlas.curve_polyline(curve, 0.01, 8192, Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let count = std::cell::Cell::new(0);
    let cancel = || {
        let n = count.get();
        count.set(n + 1);
        n > 100
    };
    assert!(matches!(
        atlas.curve_polyline(curve, 0.01, 8192, Default::default(), Some(&cancel)),
        Err(Error::Cancelled)
    ));
    let patch = atlas.patches()[0];
    let law = atlas.patch_law(patch, Default::default(), None).unwrap();
    let prefix = atlas.work_used() + law.work_used();
    assert!(matches!(
        atlas.curve_polyline(
            curve,
            0.01,
            8192,
            StockClipOptions {
                maximum_work: prefix,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(
        atlas
            .curve_polyline(curve, 1e-30, 8192, Default::default(), None)
            .is_err()
    );
    assert_eq!(*atlas.source(), source);
    assert_eq!(
        atlas
            .curve_polyline(curve, 0.01, 8192, Default::default(), None)
            .unwrap()
            .points_mm(),
        good.points_mm()
    );
}
#[test]
fn real_initial_output_storage_is_charged_before_any_polyline_allocation() {
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let fixed = size_of::<StockClipPlan>() + size_of::<StockClipPolyline>() + 2 * SCRATCH;
    assert!(matches!(
        atlas.curve_polyline(
            atlas.curves()[0],
            0.01,
            8192,
            StockClipOptions {
                maximum_bytes: fixed,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(
        atlas
            .curve_polyline(atlas.curves()[0], 0.01, 8192, Default::default(), None)
            .is_ok()
    );
}
#[test]
fn independent_point_to_segment_samples_check_original_curve_distance() {
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    for curve in atlas.physical_curves() {
        let line = atlas
            .curve_polyline(curve, 0.01, 8192, Default::default(), None)
            .unwrap();
        for i in 0..=128 {
            let t = i as f64 / 128.;
            let source = atlas
                .curve_bounds(curve, [t, t], Default::default(), None)
                .unwrap();
            let p = source.position_mm().map(|v| I::new(v[0], v[1]).midpoint());
            let distance = line
                .points_mm()
                .windows(2)
                .map(|ab| {
                    let d: [f64; 3] = std::array::from_fn(|i| ab[1][i] - ab[0][i]);
                    let v: [f64; 3] = std::array::from_fn(|i| p[i] - ab[0][i]);
                    let norm = d.iter().map(|x| x * x).sum::<f64>();
                    let t = ((0..3).map(|i| v[i] * d[i]).sum::<f64>() / norm).clamp(0., 1.);
                    let delta: [f64; 3] = std::array::from_fn(|i| v[i] - t * d[i]);
                    delta[0].hypot(delta[1]).hypot(delta[2])
                })
                .fold(f64::INFINITY, f64::min);
            assert!(
                distance <= line.source_chord_error_mm() + 1e-10,
                "{curve:?}/{t} {distance}"
            );
        }
    }
}
