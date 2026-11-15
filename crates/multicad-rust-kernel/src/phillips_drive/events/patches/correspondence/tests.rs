use super::*;
#[test]
fn real_active_cutter_fields_compare_without_domain_equality() {
    let original = Cutter::new(
        Spec {
            mouth_diameter_mm: 1.9,
            mouth_radius_mm: 0.3,
            b_width_mm: 0.61,
            g_diameter_mm: 0.81,
            tilt_angle_deg: 7.,
            outer_flank_angle_deg: 26.5,
            floor_angle_deg: 28.,
            outer_corner_angle_deg: 92.,
            inner: Inner::Rounded { radius_mm: 0.25 },
        },
        Domain {
            corner_start_depth_mm: -4.,
            corner_end_depth_mm: 4.,
            corner_reach_mm: 20.,
        },
    )
    .unwrap();
    let mut spec = original.spec();
    spec.b_width_mm = f64::from_bits(spec.b_width_mm.to_bits() + 1);
    let d = original.domain();
    let changed = Cutter::new(
        spec,
        Domain {
            corner_start_depth_mm: d.corner_start_depth_mm * 2.,
            corner_end_depth_mm: d.corner_end_depth_mm * 2.,
            corner_reach_mm: d.corner_reach_mm * 2.,
        },
    )
    .unwrap();
    let a = original.surface_patches().unwrap();
    let b = changed.surface_patches().unwrap();
    let now = std::time::Instant::now();
    let r = a
        .source_correspondence(&b, Default::default(), None)
        .unwrap();
    eprintln!("cavity whole patch C0C1 {r:?} elapsed {:?}", now.elapsed());
    assert!(r.position_mm < 1e-7 && r.du_mm < 1e-7 && r.dv_mm < 1e-7);
}

#[test]
fn real_correspondence_numeric_cell_and_cancellation_frontiers_are_observable() {
    use std::cell::Cell;
    let (spec, domain) = crate::phillips_drive::tests::fixture(0);
    let cutter = Cutter::new(spec, domain).unwrap();
    let atlas = cutter.surface_patches().unwrap();
    let good = atlas
        .source_correspondence(&atlas, Default::default(), None)
        .unwrap();
    for options in [
        PatchCorrespondenceOptions {
            maximum_position_mm: f64::NAN,
            ..Default::default()
        },
        PatchCorrespondenceOptions {
            maximum_position_mm: 0.,
            ..Default::default()
        },
        PatchCorrespondenceOptions {
            maximum_derivative_mm: f64::INFINITY,
            ..Default::default()
        },
        PatchCorrespondenceOptions {
            maximum_derivative_mm: 1.,
            ..Default::default()
        },
        PatchCorrespondenceOptions {
            maximum_cells: CELLS + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            atlas.source_correspondence(&atlas, options, None),
            Err(Error::Invalid)
        ));
    }
    for cells in [0, good.cells_used - 1] {
        assert!(matches!(
            atlas.source_correspondence(
                &atlas,
                PatchCorrespondenceOptions {
                    maximum_cells: cells,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    assert!(
        atlas
            .source_correspondence(
                &atlas,
                PatchCorrespondenceOptions {
                    maximum_cells: good.cells_used,
                    ..Default::default()
                },
                None
            )
            .is_ok()
    );
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    atlas
        .source_correspondence(&atlas, Default::default(), Some(&count))
        .unwrap();
    let total = calls.get();
    for stop in [1, 2, total / 2, total] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert!(matches!(
            atlas.source_correspondence(&atlas, Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    let changed = Cutter::new(
        Spec {
            b_width_mm: spec.b_width_mm + 0.01,
            ..spec
        },
        domain,
    )
    .unwrap();
    let other = changed.surface_patches().unwrap();
    assert!(
        atlas
            .source_correspondence(
                &other,
                PatchCorrespondenceOptions {
                    maximum_cells: 100,
                    ..Default::default()
                },
                None
            )
            .is_err()
    );
    for i in 0..atlas.patches().len() {
        for u in [0., 0.25, 0.5, 0.75, 1.] {
            for v in [0., 0.5, 1.] {
                let a = atlas.patch_bounds(i, [[u, u], [v, v]]).unwrap();
                let b = atlas.patch_bounds(i, [[u, u], [v, v]]).unwrap();
                for (left, right, bound) in [
                    (a.point_mm, b.point_mm, good.position_mm),
                    (a.du_mm, b.du_mm, good.du_mm),
                    (a.dv_mm, b.dv_mm, good.dv_mm),
                ] {
                    for k in 0..3 {
                        assert!(
                            (I::new(left[k][0], left[k][1]) - I::new(right[k][0], right[k][1]))
                                .magnitude()
                                <= bound + 1e-10
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn actual_canonical_cutter_families_retain_distinct_patch_inventories() {
    let mut inventories = Vec::new();
    for index in 0..5 {
        let (spec, domain) = crate::phillips_drive::tests::fixture(index);
        let cutter = Cutter::new(spec, domain).unwrap();
        let atlas = cutter.surface_patches().unwrap();
        atlas
            .source_correspondence(&atlas, Default::default(), None)
            .unwrap();
        inventories.push(cutter);
    }
    for a in &inventories {
        for b in &inventories {
            let a = a.surface_patches().unwrap();
            let b = b.surface_patches().unwrap();
            if a.patches().len() != b.patches().len()
                || a.patches()
                    .iter()
                    .zip(b.patches())
                    .any(|(a, b)| a.source_identity() != b.source_identity())
            {
                assert!(matches!(
                    a.source_correspondence(&b, Default::default(), None),
                    Err(Error::Unresolved)
                ));
            }
        }
    }
}
