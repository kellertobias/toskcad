use super::*;
use std::cell::Cell;
#[test]
fn actual_hardware_contracts_and_clip_phase_cancellation_remain_bounded() {
    let thread = source(false, true, Some(5.1));
    let kind = Kind::HexHead {
        across_flats_mm: 16.,
        height_mm: 3.,
    };
    for options in [
        Options {
            limits: StockClipOptions {
                maximum_work: 0,
                ..Default::default()
            },
            ..Default::default()
        },
        Options {
            limits: StockClipOptions {
                maximum_work: 16_000_001,
                ..Default::default()
            },
            ..Default::default()
        },
        Options {
            limits: StockClipOptions {
                maximum_bytes: 0,
                ..Default::default()
            },
            ..Default::default()
        },
        Options {
            limits: StockClipOptions {
                maximum_bytes: 64 * 1024 * 1024 + 1,
                ..Default::default()
            },
            ..Default::default()
        },
        Options {
            cap_area_width_mm2: f64::NAN,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            Material::new(thread, kind, options, None),
            Err(Error::Invalid)
        ));
    }
    for kind in [
        Kind::HexHead {
            across_flats_mm: 1e7 + 1.,
            height_mm: 3.,
        },
        Kind::HexHead {
            across_flats_mm: 16.,
            height_mm: 1e7,
        },
    ] {
        assert!(matches!(
            Material::new(thread, kind, Default::default(), None),
            Err(Error::Limit)
        ));
    }
    for stop in [2, 6] {
        let calls = Cell::new(0);
        let interrupt = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        assert!(matches!(
            Material::new(thread, kind, Default::default(), Some(&interrupt)),
            Err(Error::Cancelled)
        ));
    }
    let contact = AuthoredThread::new(
        AuthoredSpec {
            crest_radius_limit_mm: Some(5.),
            ..thread.spec()
        },
        1e-7,
        1e-6,
    )
    .unwrap();
    assert!(matches!(
        Material::new(contact, kind, Default::default(), None),
        Err(Error::Unresolved)
    ));
    let body = Material::new(thread, kind, Default::default(), None).unwrap();
    let box_mm = [[0., 0.], [0., 0.], [3., 3.]];
    for limits in [
        StockClipOptions {
            maximum_work: 0,
            ..Default::default()
        },
        StockClipOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            body.classify_box(box_mm, limits, None),
            Err(Error::Invalid)
        ));
    }
}

#[test]
fn actual_near_terminal_clip_cap_retains_both_valid_area_branches() {
    let base = source(false, true, None).spec();
    let limit = (base.major_radius_mm + base.radial_taper * base.length_mm).next_down();
    let thread = AuthoredThread::new(
        AuthoredSpec {
            crest_radius_limit_mm: Some(limit),
            ..base
        },
        1e-7,
        1e-6,
    )
    .unwrap();
    let body = Material::new(
        thread,
        Kind::HexHead {
            across_flats_mm: 16.,
            height_mm: 3.,
        },
        Default::default(),
        None,
    )
    .unwrap();
    let plan = body.clipped_thread().unwrap();
    assert!(plan.contact_height_bounds_mm()[1] < base.length_mm);
    assert!(base.length_mm - plan.contact_height_bounds_mm()[0] < 1e-10);
    overlap(
        plan.face_measure(StockClipFace::EndCap, Default::default(), None)
            .unwrap()
            .area_mm2(),
        body.original_cap_areas_mm2()[1],
    );
}

fn width(bounds: [f64; 2]) -> f64 {
    I::new(bounds[0], bounds[1]).width()
}
#[test]
fn actual_hex_minus_cap_precision_is_checked_after_each_positive_input() {
    let thread = source(false, false, None);
    // A broad real hex makes subtraction rounding significant at the first
    // cap; ordinary hardware separates the second cap's larger input width.
    for (across, index) in [(1024., 0usize), (16., 1)] {
        let kind = Kind::HexNut {
            across_flats_mm: across,
        };
        let body = Material::new(
            thread,
            kind,
            Options {
                volume_width_mm3: 1e-4,
                ..Default::default()
            },
            None,
        )
        .unwrap();
        let input = body
            .original_cap_areas_mm2()
            .into_iter()
            .map(width)
            .fold(width(body.hex_area_mm2()), f64::max);
        let planes = body.planes().as_slice();
        let lower = if index == 0 {
            input
        } else {
            input.max(width(planes[0].area_mm2))
        };
        let upper = width(planes[index].area_mm2);
        assert!(lower < upper, "across={across}: {lower:e} vs {upper:e}");
        let budget = lower.midpoint(upper);
        assert!(
            body.original_cap_areas_mm2()
                .into_iter()
                .all(|v| width(v) < budget)
        );
        assert!(width(body.hex_area_mm2()) < budget);
        assert!(matches!(
            Material::new(
                thread,
                kind,
                Options {
                    volume_width_mm3: 1e-4,
                    cap_area_width_mm2: budget,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Unresolved)
        ));
    }
}
