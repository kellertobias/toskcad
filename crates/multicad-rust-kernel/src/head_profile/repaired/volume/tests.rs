use super::*;
use crate::head_profile::{Head, Profile};
#[test]
fn actual_repaired_cylinder_and_curved_heads_have_new_source_volume() {
    for head in [
        Head::Cylindrical {
            radius_mm: 2.,
            height_mm: 1.6,
        },
        Head::Crowned {
            radius_mm: 2.,
            height_mm: 1.6,
            crown_radius_mm: 3.2,
            top_land_radius_mm: 1.25,
            neck_radius_mm: 0.1,
        },
        Head::Countersunk {
            radius_mm: 3.,
            height_mm: 1.65,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.1,
            neck_radius_mm: 0.1,
        },
    ] {
        let shank = if matches!(head, Head::Countersunk { .. }) {
            1.5
        } else {
            1.
        };
        let profile = Profile::new_for_thread_union(head, shank, 0.).unwrap();
        let source = profile
            .closed_repaired_meridian(Default::default(), None)
            .unwrap();
        let volume = source.volume_bounds(Default::default(), None).unwrap();
        assert!(std::ptr::eq(volume.source(), &source));
        let [lo, hi] = volume.volume_mm3();
        assert!(lo > 0. && hi >= lo && hi - lo <= 1e-6);
        if let Head::Cylindrical {
            radius_mm: r,
            height_mm: h,
        } = head
        {
            let expected = pi() * I::point(r).square() * I::point(h);
            assert!(lo <= expected.lo && expected.hi <= hi);
        }
        eprintln!(
            "HEADV {head:?} {:?} cells{} work{}",
            volume.volume_mm3(),
            volume.cells(),
            volume.work_used()
        );
    }
}

fn quarter() -> Profile {
    Profile::from_retained_segments(&[Segment {
        role: Role::Crown,
        start: [1., -1.],
        end: [2., 0.],
        curve: Curve::Arc {
            center: [1., 0.],
            radius_mm: 1.,
            sweep: std::f64::consts::FRAC_PI_2,
        },
    }])
    .unwrap()
}
#[test]
fn normalized_quarter_circle_contains_independent_closed_form() {
    let profile = quarter();
    let source = profile
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let result = source.volume_bounds(Default::default(), None).unwrap();
    // r(z)=1+sqrt(1-z²), -1<=z<=0: integral π*r² = 5π/3+π²/2.
    let expected = pi() * I::point(5.) / I::point(3.) + pi().square() / I::point(2.);
    let [lo, hi] = result.volume_mm3();
    assert!(
        lo <= expected.lo && hi >= expected.hi,
        "{lo} {hi} vs {expected:?}"
    );
}
#[test]
fn real_resource_frontiers_cancel_and_numeric_floor_are_atomic() {
    let profile = Profile::new_for_thread_union(
        Head::Cylindrical {
            radius_mm: 2.,
            height_mm: 1.6,
        },
        1.,
        0.,
    )
    .unwrap();
    let source = profile
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let receipt = source.volume_bounds(Default::default(), None).unwrap();
    let work = receipt.work_used();
    let cells = receipt.cells();
    let peak = std::mem::size_of::<RepairedHeadVolume<'_, '_>>()
        + source.retained_budget_bytes()
        + 4 * (cells + 1) * std::mem::size_of::<Cell>()
        + SCRATCH;
    let options = RepairedHeadVolumeOptions {
        maximum_work: work,
        maximum_cells: cells,
        maximum_bytes: peak,
        ..Default::default()
    };
    assert!(source.volume_bounds(options, None).is_ok());
    for options in [
        RepairedHeadVolumeOptions {
            maximum_work: work - 1,
            ..options
        },
        RepairedHeadVolumeOptions {
            maximum_work: source.work - 1,
            ..options
        },
        RepairedHeadVolumeOptions {
            maximum_cells: cells - 1,
            ..options
        },
        RepairedHeadVolumeOptions {
            maximum_bytes: peak - 1,
            ..options
        },
    ] {
        assert!(matches!(
            source.volume_bounds(options, None),
            Err(RepairedHeadError::Limit)
        ));
    }
    for options in [
        RepairedHeadVolumeOptions {
            maximum_width_mm3: f64::NAN,
            ..options
        },
        RepairedHeadVolumeOptions {
            maximum_width_mm3: 0.,
            ..options
        },
        RepairedHeadVolumeOptions {
            maximum_cells: 0,
            ..options
        },
        RepairedHeadVolumeOptions {
            maximum_cells: CELLS + 1,
            ..options
        },
        RepairedHeadVolumeOptions {
            maximum_work: WORK + 1,
            ..options
        },
        RepairedHeadVolumeOptions {
            maximum_bytes: BYTES + 1,
            ..options
        },
    ] {
        assert!(matches!(
            source.volume_bounds(options, None),
            Err(RepairedHeadError::Invalid)
        ));
    }
    let count = std::cell::Cell::new(0);
    source
        .volume_bounds(
            options,
            Some(&|| {
                count.set(count.get() + 1);
                false
            }),
        )
        .unwrap();
    for at in 0..count.get() {
        let n = std::cell::Cell::new(0);
        assert!(matches!(
            source.volume_bounds(
                options,
                Some(&|| {
                    let old = n.get();
                    n.set(old + 1);
                    old == at
                })
            ),
            Err(RepairedHeadError::Cancelled)
        ));
    }
    let large = Profile::new_for_thread_union(
        Head::Cylindrical {
            radius_mm: 9e6,
            height_mm: 1e6,
        },
        8e6,
        0.,
    )
    .unwrap();
    let large = large
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(matches!(
        large.volume_bounds(Default::default(), None),
        Err(RepairedHeadError::Unresolved)
    ));
    assert!(
        large
            .volume_bounds(
                RepairedHeadVolumeOptions {
                    maximum_width_mm3: 1e8,
                    ..Default::default()
                },
                None
            )
            .is_ok()
    );
}
#[test]
fn curved_source_respects_cell_work_and_precision_limits() {
    let profile = quarter();
    let source = profile
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    for options in [
        RepairedHeadVolumeOptions {
            maximum_cells: 1,
            ..Default::default()
        },
        RepairedHeadVolumeOptions {
            maximum_work: source.work + 4095,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            source.volume_bounds(options, None),
            Err(RepairedHeadError::Limit)
        ));
    }
    let tiny = source.volume_bounds(
        RepairedHeadVolumeOptions {
            maximum_width_mm3: 1e-300,
            ..Default::default()
        },
        None,
    );
    assert!(
        matches!(tiny, Err(RepairedHeadError::Unresolved)),
        "{:?}",
        tiny.err()
    );
    let count = std::cell::Cell::new(0);
    source
        .volume_bounds(
            Default::default(),
            Some(&|| {
                count.set(count.get() + 1);
                false
            }),
        )
        .unwrap();
    for at in [
        1,
        2,
        3,
        4,
        5,
        count.get() / 2,
        count.get() - 2,
        count.get() - 1,
    ] {
        let n = std::cell::Cell::new(0);
        assert!(matches!(
            source.volume_bounds(
                Default::default(),
                Some(&|| {
                    let old = n.get();
                    n.set(old + 1);
                    old == at
                })
            ),
            Err(RepairedHeadError::Cancelled)
        ));
    }
}
#[test]
fn actual_subnormal_material_cannot_publish_an_uncertified_positive_volume() {
    let r = 1e-108;
    let profile = Profile::from_retained_segments(&[Segment {
        role: Role::Skirt,
        start: [r, -r],
        end: [r, 0.],
        curve: Curve::Line,
    }])
    .unwrap();
    let source = profile
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(matches!(
        source.volume_bounds(Default::default(), None),
        Err(RepairedHeadError::Unresolved)
    ));
}
