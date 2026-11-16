use super::*;
use std::cell::Cell;

#[test]
fn actual_retained_file_queries_reject_invalid_quality_and_budget_contracts() {
    let source = source(false, false);
    let input = text(false);
    for options in [
        Options {
            maximum_work: 0,
            ..Default::default()
        },
        Options {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        Options {
            maximum_bytes: 0,
            ..Default::default()
        },
        Options {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
        Options {
            position_error_mm: f64::NAN,
            ..Default::default()
        },
        Options {
            original_derivative_error_mm: 0.,
            ..Default::default()
        },
        Options {
            normal_chord_error: -1.,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            verify_exact_replay(&source, options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        ReadbackOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        ReadbackOptions {
            maximum_bytes: 128 * 1024 * 1024 + 1,
            ..Default::default()
        },
        ReadbackOptions {
            normal_chord_error: f64::INFINITY,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            verify_readback(&source, None, &input, 0, options, None),
            Err(Error::Invalid)
        ));
    }
    let fresh = verify_readback(&source, None, &input, 0, Default::default(), None).unwrap();
    assert_eq!(
        fresh.receipt().matched_entities(),
        source.source_entity_count()
    );
    assert!(matches!(
        verify_readback(
            &source,
            None,
            &input,
            0,
            ReadbackOptions {
                maximum_work: fresh.receipt().query_work_used() - 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
}

#[test]
fn actual_fresh_readback_cancellation_reaches_document_material_and_binding_phases() {
    let source = source(false, false);
    let input = text(false);
    let calls = Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let before =
        verify_readback(&source, None, &input, 0, Default::default(), Some(&count)).unwrap();
    let total = calls.get();
    assert!(total > 1024);
    for stop in [2, 32, 256, total / 8, total / 2, total - 1] {
        let calls = Cell::new(0usize);
        let interrupt = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        assert!(
            matches!(
                verify_readback(
                    &source,
                    None,
                    &input,
                    0,
                    Default::default(),
                    Some(&interrupt)
                ),
                Err(Error::Cancelled)
            ),
            "stop={stop}/{total}"
        );
    }
    let after = verify_readback(&source, None, &input, 0, Default::default(), None).unwrap();
    assert_eq!(before.actual_volume_mm3(), after.actual_volume_mm3());
    assert_eq!(before.actual_bounds_mm(), after.actual_bounds_mm());
}

#[test]
fn actual_file_units_and_finite_operator_with_overflowing_world_volume_are_refused() {
    let source = source(false, false);
    let input = text(false);
    let wrong_units = input.replacen("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)", 1);
    assert_ne!(wrong_units, input);
    assert!(matches!(
        CompoundDocument::parse(&wrong_units, Default::default(), None),
        Err(crate::step_read::owned_insert::Error::Unsupported)
    ));
    assert!(matches!(
        verify_readback(&source, None, &wrong_units, 0, Default::default(), None),
        Err(Error::Unsupported)
    ));
    let mapped = mapped(&input);
    let overflowing = mapped.replacen("#900014,1.)", "#900014,1.E102)", 1);
    assert_ne!(overflowing, mapped);
    let document = CompoundDocument::parse(&overflowing, Default::default(), None).unwrap();
    let selected = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let pose = selected.placement(Default::default(), None).unwrap();
    let determinant = pose.determinant_bounds();
    assert!(determinant[0] > 0. && determinant[1].is_finite());
    let occupied = selected
        .parsed()
        .projection(Default::default(), None)
        .unwrap()
        .occupied(Default::default(), None)
        .unwrap();
    let local = occupied.volume_mm3();
    assert!(local[0] > 0. && local[1].is_finite());
    assert!(!(I::new(local[0], local[1]) * I::new(determinant[0], determinant[1])).valid());
    assert!(matches!(
        verify_readback(&source, None, &overflowing, 0, Default::default(), None),
        Err(Error::Unsupported)
    ));
}
