use super::*;
#[test]
fn independent_d11_reader_detaches_actual_fitted_source_after_inputs_drop() {
    let text = super::super::derived_thread::actual_m3(true, false);
    let read = read_degree11_phillips_material(&text, Default::default(), None).unwrap();
    assert_eq!(read.reserved_delegated_work(), TOTAL);
    assert!(read.material_proof_work_used() <= 512_000_000);
    assert!(read.material_proof_work_used() > 64_000_000);
    assert_eq!(read.material().original_input(), text);
    let expected = read.material().world_volume_mm3();
    let owned = read.into_material();
    drop(text);
    assert_eq!(owned.world_volume_mm3(), expected);
    assert!(
        owned
            .thread_faces()
            .iter()
            .all(|f| f.source().spans().iter().all(|&n| n > 0))
    );
    assert!(owned.local_volume_mm3()[0] > 0.);
}

#[test]
fn independent_reader_preserves_real_flat_control_and_fitted_both_hand_frames() {
    for fitted in [false, true] {
        for left in [false, true] {
            let text = super::super::derived_thread::actual_m3(fitted, left);
            let read = read_degree11_phillips_material(&text, Default::default(), None).unwrap();
            assert_eq!(read.material().original_input(), text);
            assert_eq!(read.material().thread_faces().len(), 4);
            assert!(
                read.material()
                    .containment_clearances_mm()
                    .iter()
                    .all(|&x| x > 0.)
            );
            assert!(read.material().join_radial_margin_mm() > 0.);
            assert!(read.material().world_volume_mm3()[0] > 0.);
            // These are actual original D11 cells, never a renamed D9 source.
            for face in read.material().thread_faces() {
                let [u, v] = face.source().spans();
                assert_eq!(face.source().dimensions(), [11 * u + 1, 11 * v + 1]);
                assert!(face.transfer().distance_mm < 1e-7);
                assert!(face.transfer().derivative_mm.iter().all(|&x| x < 1e-7));
            }
            eprintln!(
                "D11_READ fitted={fitted} left={left} work={} bytes={} volume={:?}",
                read.material_proof_work_used(),
                read.material().retained_budget_bytes(),
                read.material().world_volume_mm3()
            );
        }
    }
}
#[test]
fn actual_reader_refuses_bad_limits_and_corrupt_file_before_publication() {
    let text = super::super::derived_thread::actual_m3(true, false);
    let options = Degree11ReadOptions::default();
    let before = text.clone();
    for options in [
        Degree11ReadOptions {
            maximum_reserved_work: TOTAL + 1,
            ..options
        },
        Degree11ReadOptions {
            maximum_input_bytes: INPUT_BYTES + 1,
            ..options
        },
        Degree11ReadOptions {
            stage_work_limits: [0; STAGES],
            ..options
        },
    ] {
        let error = read_degree11_phillips_material(&text, options, None)
            .err()
            .unwrap();
        assert_eq!(error, fail(PhillipsReadError::Invalid));
    }
    for options in [
        Degree11ReadOptions {
            maximum_reserved_work: TOTAL - 1,
            ..options
        },
        Degree11ReadOptions {
            maximum_input_bytes: text.len() - 1,
            ..options
        },
    ] {
        let error = read_degree11_phillips_material(&text, options, None)
            .err()
            .unwrap();
        assert_eq!(error, fail(PhillipsReadError::Limit));
    }
    let error = read_degree11_phillips_material(&text[..text.len() / 2], options, None)
        .err()
        .unwrap();
    assert_eq!(error.stage, Some(Degree11ReadStage::Body));
    assert!(matches!(
        error.cause,
        PhillipsReadError::Invalid | PhillipsReadError::Unsupported
    ));
    let error = read_degree11_phillips_material(&text, options, Some(&|| true))
        .err()
        .unwrap();
    assert_eq!(error, fail(PhillipsReadError::Cancelled));
    assert_eq!(text, before);
    assert!(read_degree11_phillips_material(&text, options, None).is_ok());
}
#[test]
fn real_reader_child_work_frontiers_and_staged_cancellation_do_not_substitute_a_recipe() {
    use std::cell::Cell;
    let text = super::super::derived_thread::actual_m3(true, false);
    let options = Degree11ReadOptions::default();
    for stage in [
        Degree11ReadStage::Outer,
        Degree11ReadStage::Carrier,
        Degree11ReadStage::Volume,
        Degree11ReadStage::Owned,
    ] {
        let mut limits = HARD;
        limits[stage as usize] = 1;
        let error = read_degree11_phillips_material(
            &text,
            Degree11ReadOptions {
                stage_work_limits: limits,
                ..options
            },
            None,
        )
        .err()
        .unwrap();
        assert_eq!(
            error,
            Degree11ReadError {
                stage: Some(stage),
                cause: PhillipsReadError::Limit
            }
        );
    }
    for stage in [
        Degree11ReadStage::Body,
        Degree11ReadStage::Head,
        Degree11ReadStage::Owned,
    ] {
        let cancelled = Cell::new(false);
        let progress = |s| {
            if s == stage {
                cancelled.set(true)
            }
        };
        let error = read_degree11_phillips_material_with_progress(
            &text,
            options,
            Some(&|| cancelled.get()),
            Some(&progress),
        )
        .err()
        .unwrap();
        assert_eq!(
            error,
            Degree11ReadError {
                stage: Some(stage),
                cause: PhillipsReadError::Cancelled
            }
        );
    }
    let polls = Cell::new(0usize);
    let observe = || {
        polls.set(polls.get() + 1);
        false
    };
    let mut seen = Vec::new();
    // Progress callbacks are read-only observations; a RefCell stores the
    // observed order without changing a parser/proof/source object.
    let stages = std::cell::RefCell::new(Vec::new());
    let read = read_degree11_phillips_material_with_progress(
        &text,
        options,
        Some(&observe),
        Some(&|s| stages.borrow_mut().push(s)),
    )
    .unwrap();
    seen.extend(stages.into_inner());
    assert_eq!(seen.len(), STAGES);
    assert_eq!(seen[0], Degree11ReadStage::Body);
    assert_eq!(*seen.last().unwrap(), Degree11ReadStage::Owned);
    let stop = polls.get() - 1;
    let polls = Cell::new(0usize);
    let cancel = || {
        let n = polls.get();
        polls.set(n + 1);
        n == stop
    };
    let error = read_degree11_phillips_material(&text, options, Some(&cancel))
        .err()
        .unwrap();
    assert_eq!(
        error,
        Degree11ReadError {
            stage: Some(Degree11ReadStage::Owned),
            cause: PhillipsReadError::Cancelled
        }
    );
    assert_eq!(read.material().original_input(), text);
    assert!(read_degree11_phillips_material(&text, options, None).is_ok());
}
