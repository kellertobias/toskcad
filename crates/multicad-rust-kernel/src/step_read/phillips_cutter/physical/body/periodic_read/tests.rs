use super::*;
use std::cell::{Cell, RefCell};
fn text() -> String {
    let original = crate::step_read::phillips_cutter::physical::tests::fixture(0, false, true, 1);
    let owner = OwnedPeriodicMaterial::from_authored(&original, Default::default(), None).unwrap();
    let proposal = owner.authored_mesh_approximation(16_000_000, None).unwrap();
    crate::step::write_phillips_screw(
        proposal.proposal(),
        [0.; 3],
        [1., 0., 0.],
        [0., 1., 0.],
        1e-7,
    )
    .unwrap()
}
#[test]
fn actual_file_owns_independent_periodic_material_under_an_explicit_staged_contract() {
    let text = text();
    let now = std::time::Instant::now();
    let read = read_periodic_phillips_material(&text, Default::default(), None).unwrap();
    assert_eq!(read.material().original_input(), Some(text.as_str()));
    assert_eq!(read.reserved_delegated_work(), 224_000_000);
    assert!(read.material_proof_work_used() < 64_000_000 && read.material_proof_work_used() > 0);
    assert!(read.material().world_volume_mm3()[0] > 0.);
    let input = read.into_material();
    assert_eq!(input.original_input(), Some(text.as_str()));
    let progress = RefCell::new(Vec::new());
    let record = |stage| progress.borrow_mut().push(stage);
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let repeated = read_periodic_phillips_material_with_progress(
        &text,
        Default::default(),
        Some(&count),
        Some(&record),
    )
    .unwrap();
    let complete = calls.get();
    calls.set(0);
    let late_cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == complete
    };
    assert_eq!(
        read_periodic_phillips_material(&text, Default::default(), Some(&late_cancel))
            .err()
            .unwrap(),
        PeriodicReadError {
            stage: Some(PeriodicReadStage::Owned),
            cause: PhillipsReadError::Cancelled
        }
    );
    assert_eq!(progress.borrow().len(), STAGES);
    for (i, s) in progress.borrow().iter().enumerate() {
        assert_eq!(*s as usize, i);
    }
    assert_eq!(
        repeated.material().world_volume_mm3(),
        input.world_volume_mm3()
    );
    eprintln!(
        "actual owning staged periodic read {:?} reserved{} proof{} retained{}",
        now.elapsed(),
        repeated.reserved_delegated_work(),
        repeated.material_proof_work_used(),
        repeated.material().retained_budget_bytes()
    );
}
#[test]
fn actual_reader_stage_budgets_and_progress_cancellation_never_publish_partial_material() {
    let text = text();
    let good = PeriodicReadOptions::default();
    for bad in [
        PeriodicReadOptions {
            maximum_reserved_work: TOTAL + 1,
            ..good
        },
        PeriodicReadOptions {
            maximum_input_bytes: INPUT_BYTES + 1,
            ..good
        },
    ] {
        assert_eq!(
            read_periodic_phillips_material(&text, bad, None)
                .err()
                .unwrap(),
            fail(PhillipsReadError::Invalid)
        );
    }
    for bad in [
        PeriodicReadOptions {
            maximum_reserved_work: 223_999_999,
            ..good
        },
        PeriodicReadOptions {
            maximum_input_bytes: text.len() - 1,
            ..good
        },
    ] {
        assert_eq!(
            read_periodic_phillips_material(&text, bad, None)
                .err()
                .unwrap(),
            fail(PhillipsReadError::Limit)
        );
    }
    assert_eq!(
        read_periodic_phillips_material(&text, good, Some(&|| true))
            .err()
            .unwrap(),
        fail(PhillipsReadError::Cancelled)
    );
    let polls = Cell::new(0);
    let early = || {
        polls.set(polls.get() + 1);
        polls.get() == 2
    };
    assert_eq!(
        read_periodic_phillips_material(&text, good, Some(&early))
            .err()
            .unwrap(),
        PeriodicReadError {
            stage: Some(PeriodicReadStage::Body),
            cause: PhillipsReadError::Cancelled
        }
    );
    let stages = [
        PeriodicReadStage::Body,
        PeriodicReadStage::Outer,
        PeriodicReadStage::Endpoints,
        PeriodicReadStage::Analytic,
        PeriodicReadStage::Lifts,
        PeriodicReadStage::Cubics,
        PeriodicReadStage::Sewn,
        PeriodicReadStage::Representative,
        PeriodicReadStage::ClosedCavity,
        PeriodicReadStage::CommonMouth,
        PeriodicReadStage::Head,
        PeriodicReadStage::Periodic,
        PeriodicReadStage::Finite,
        PeriodicReadStage::Shoulder,
        PeriodicReadStage::Material,
        PeriodicReadStage::Volume,
        PeriodicReadStage::Owned,
    ];
    for stage in stages {
        let mut bad = good;
        bad.stage_work_limits[stage as usize] = 0;
        assert_eq!(
            read_periodic_phillips_material(&text, bad, None)
                .err()
                .unwrap(),
            fail(PhillipsReadError::Invalid)
        );
        bad.stage_work_limits[stage as usize] = HARD[stage as usize] + 1;
        assert_eq!(
            read_periodic_phillips_material(&text, bad, None)
                .err()
                .unwrap(),
            fail(PhillipsReadError::Invalid)
        );
        bad.stage_work_limits[stage as usize] = 1;
        assert_eq!(
            read_periodic_phillips_material(&text, bad, None)
                .err()
                .unwrap(),
            PeriodicReadError {
                stage: Some(stage),
                cause: PhillipsReadError::Limit
            },
            "stage {stage:?}"
        );
        let stop = Cell::new(false);
        let progress = |s| {
            if s == stage {
                stop.set(true)
            }
        };
        assert_eq!(
            read_periodic_phillips_material_with_progress(
                &text,
                good,
                Some(&|| stop.get()),
                Some(&progress)
            )
            .err()
            .unwrap(),
            PeriodicReadError {
                stage: Some(stage),
                cause: PhillipsReadError::Cancelled
            }
        );
    }
}
