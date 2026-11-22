use super::super::super::degree11::tests::actual_m3_parts;
use super::*;
use crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head;
fn text(source: &PhillipsScrew, [o, u, v]: [Point; 3]) -> String {
    crate::step::write_phillips_screw(source, o, u, v, 1e-7).unwrap()
}
#[test]
fn real_supplied_m3_author_and_frame_have_original_cell_bounds_without_metadata_alias() {
    for fitted in [false, true] {
        for left in [false, true] {
            let (source, frame) = actual_m3_parts(fitted, left);
            let input = text(&source, frame);
            with_head(&input, &mut |head| {
                let before = head.body().inferred_thread_spec();
                let receipt = head
                    .supplied_thread_correspondence(&source, frame, Default::default(), None)
                    .unwrap();
                assert!(std::ptr::eq(receipt.source(), &source));
                assert!(std::ptr::eq(receipt.body(), head.body()));
                assert_eq!(
                    receipt.source().thread().spec().root_radius_mm.is_some(),
                    fitted
                );
                assert_eq!(receipt.frame(), frame);
                let determinant = receipt.frame_determinant_bounds();
                assert!(determinant[0] > 0. && determinant[0] <= 1. && determinant[1] >= 1.);
                assert_eq!(head.body().inferred_thread_spec(), before);
                let mut c0 = 0f64;
                let mut du = 0f64;
                let mut dv = 0f64;
                let mut cells = 0usize;
                for t in receipt.transfers() {
                    c0 = c0.max(t.distance_mm);
                    du = du.max(t.derivative_mm[0]);
                    dv = dv.max(t.derivative_mm[1]);
                    cells += t.original_cells[0] * t.original_cells[1];
                    assert!(
                        t.distance_mm <= 1e-7 && t.derivative_mm.into_iter().all(|d| d <= 1e-7)
                    );
                }
                assert_eq!(receipt.body().literal_input(), input);
                assert_eq!(
                    receipt.total_work_used(),
                    receipt.inherited_head_work_used() + receipt.new_work_used()
                );
                assert!(
                    receipt.new_work_used() <= WORK
                        && receipt.total_work_used() <= TOTAL
                        && receipt.retained_budget_bytes() <= BYTES
                );
                assert!(
                    receipt.retained_budget_bytes()
                        > head.retained_budget_bytes() + resident_bytes(&source)
                );
                if left {
                    assert_ne!(
                        before,
                        source.thread().spec(),
                        "posed inference must not silently stand in for exact author"
                    );
                }
                eprintln!(
                    "SUPPLIED_M3 fitted={fitted} left={left} work={} inherited={} retained={} cells={cells} c0={c0} du={du} dv={dv}",
                    receipt.new_work_used(),
                    receipt.inherited_head_work_used(),
                    receipt.retained_budget_bytes()
                );
            });
        }
    }
}
#[test]
fn actual_different_and_huge_supplied_laws_and_frames_refuse_without_mutating_body() {
    let (source, frame) = actual_m3_parts(true, false);
    let input = text(&source, frame);
    let (different, _) = actual_m3_parts(false, false);
    let mut spec = source.thread().spec();
    spec.length_mm = 60.;
    let huge = PhillipsScrew::new(
        crate::helical_thread::Thread::new(spec).unwrap(),
        *source.profile(),
        source.spec(),
        source.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    with_head(&input, &mut |head| {
        let before = head.body().inferred_thread_spec();
        assert!(matches!(
            head.supplied_thread_correspondence(&different, frame, Default::default(), None),
            Err(PhillipsReadError::Unsupported)
        ));
        assert!(matches!(
            head.supplied_thread_correspondence(&huge, frame, Default::default(), None),
            Err(PhillipsReadError::Limit) | Err(PhillipsReadError::Unsupported)
        ));
        let mut translated = frame;
        translated[0][0] += 1.;
        assert!(matches!(
            head.supplied_thread_correspondence(&source, translated, Default::default(), None),
            Err(PhillipsReadError::Unsupported)
        ));
        for invalid in [
            [[0.; 3]; 3],
            [[f64::NAN, 0., 0.], frame[1], frame[2]],
            [[f64::MAX, 0., 0.], frame[1], frame[2]],
        ] {
            assert!(matches!(
                head.supplied_thread_correspondence(&source, invalid, Default::default(), None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        assert_eq!(head.body().inferred_thread_spec(), before);
        head.supplied_thread_correspondence(&source, frame, Default::default(), None)
            .unwrap();
    });
}
#[test]
fn actual_supplied_law_work_bytes_cancel_and_exact_retry_contract() {
    let (source, frame) = actual_m3_parts(true, false);
    let input = text(&source, frame);
    with_head(&input, &mut |head| {
        let o = AnalyticThreadOptions::default();
        let good = head
            .supplied_thread_correspondence(&source, frame, o, None)
            .unwrap();
        for bad in [
            AnalyticThreadOptions {
                maximum_work: good.new_work_used() - 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_total_work: good.total_work_used() - 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_bytes: head.retained_budget_bytes()
                    + std::mem::size_of::<AnalyticThreadRepair<'_, '_>>()
                    + resident_bytes(&source)
                    + std::mem::size_of::<SuppliedThreadCorrespondence<'_, '_, '_>>()
                    + SCRATCH_ANALYTIC
                    - 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_distance_mm: 1e-20,
                ..o
            },
        ] {
            assert!(
                head.supplied_thread_correspondence(&source, frame, bad, None)
                    .is_err()
            );
        }
        let exact = AnalyticThreadOptions {
            maximum_work: good.new_work_used(),
            maximum_total_work: good.total_work_used(),
            ..o
        };
        let calls = Cell::new(0usize);
        head.supplied_thread_correspondence(
            &source,
            frame,
            exact,
            Some(&|| {
                calls.set(calls.get() + 1);
                false
            }),
        )
        .unwrap();
        for stop in [1, 2, calls.get() / 2, calls.get()] {
            let count = Cell::new(0usize);
            assert!(matches!(
                head.supplied_thread_correspondence(
                    &source,
                    frame,
                    exact,
                    Some(&|| {
                        count.set(count.get() + 1);
                        count.get() == stop
                    })
                ),
                Err(PhillipsReadError::Cancelled)
            ));
        }
        let retry = head
            .supplied_thread_correspondence(&source, frame, exact, None)
            .unwrap();
        assert_eq!(retry.new_work_used(), good.new_work_used());
        assert!(std::ptr::eq(retry.source(), good.source()));
        for (a, b) in retry.transfers().iter().zip(good.transfers()) {
            assert_eq!(
                (
                    a.face_entity,
                    a.distance_mm,
                    a.derivative_mm,
                    a.original_cells
                ),
                (
                    b.face_entity,
                    b.distance_mm,
                    b.derivative_mm,
                    b.original_cells
                )
            );
        }
    });
}

#[test]
fn supplied_frame_requires_a_directed_positive_actual_determinant() {
    assert!(valid_frame([[0.; 3], [1., 0., 0.], [0., 1., 0.]]));
    let uncertain = [[0.; 3], [1., 1., 0.], [1., 1. + f64::EPSILON, 0.]];
    assert!(!valid_frame(uncertain));
}

#[test]
fn actual_multistart_body_refuses_a_different_single_start_author_before_generation() {
    let (source, frame) = actual_m3_parts(true, false);
    let multistart =
        crate::step_read::phillips_cutter::physical::tests::fixture(0, false, false, 2);
    let input = text(&multistart, frame);
    with_head(&input, &mut |head| {
        let before = head.body().inferred_thread_spec();
        assert_eq!(before.starts, 2);
        assert!(matches!(
            head.supplied_thread_correspondence(&source, frame, Default::default(), None),
            Err(PhillipsReadError::Unsupported)
        ));
        assert_eq!(head.body().inferred_thread_spec(), before);
    });
}

#[test]
fn original_entry_rim_rotation_cannot_bypass_a_reduced_whole_source_work_budget() {
    let (source, frame) = actual_m3_parts(true, false);
    let input = text(&source, frame);
    with_head(&input, &mut |head| {
        let before = head.body().inferred_thread_spec();
        let receipt = head
            .supplied_thread_correspondence(&source, frame, Default::default(), None)
            .unwrap();
        let options = AnalyticThreadOptions {
            maximum_work: receipt.new_work_used() - 1000,
            ..Default::default()
        };
        assert!(matches!(
            head.supplied_thread_correspondence(&source, frame, options, None),
            Err(PhillipsReadError::Limit)
        ));
        assert!(
            head.supplied_thread_correspondence(&source, frame, Default::default(), None)
                .is_ok()
        );
        assert_eq!(head.body().inferred_thread_spec(), before);
    });
}
