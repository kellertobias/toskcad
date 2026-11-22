use super::*;
#[test]
fn actual_m2_analytic_thread_has_whole_cell_world_transfer_and_budget_identity() {
    degree11::tests::with_head(&mut |head| {
        let start = std::time::Instant::now();
        let receipt = head
            .analytic_thread_repair(Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(receipt.body(), head.body()));
        assert_eq!(receipt.thread().spec(), head.body().inferred_thread_spec());
        assert_eq!(receipt.frame(), head.body().inferred_thread_frame());
        assert_eq!(receipt.transfers().len(), 4);
        assert_eq!(receipt.inherited_head_work_used(), head.work_used());
        assert_eq!(
            receipt.total_work_used(),
            receipt.new_work_used() + head.work_used()
        );
        let mut cells = 0;
        for t in receipt.transfers() {
            assert!(t.face_entity > 0 && t.surface_entity > 0);
            assert!(t.original_rail_entities.into_iter().all(|id| id > 0));
            assert!(t.distance_mm < 1e-7 && t.derivative_mm.into_iter().all(|d| d < 1e-7));
            cells += t.original_cells[0] * t.original_cells[1];
        }
        assert!(cells > 0 && receipt.retained_budget_bytes() < BYTES);
        eprintln!(
            "ANALYTIC THREAD {:?} new{} inherited{} total{} retained{} bounds{:?}",
            start.elapsed(),
            receipt.new_work_used(),
            head.work_used(),
            receipt.total_work_used(),
            receipt.retained_budget_bytes(),
            receipt.transfers()
        );
    });
}
#[test]
fn actual_body_frontiers_are_atomic_and_total_budget_includes_inherited_head() {
    degree11::tests::with_head(&mut |head| {
        let o = AnalyticThreadOptions::default();
        let good = head.analytic_thread_repair(o, None).unwrap();
        let used = good.new_work_used();
        let total = good.total_work_used();
        for bad in [
            AnalyticThreadOptions {
                maximum_distance_mm: 0.,
                ..o
            },
            AnalyticThreadOptions {
                maximum_cell_derivative_mm: f64::NAN,
                ..o
            },
            AnalyticThreadOptions {
                maximum_distance_mm: 1e-6,
                ..o
            },
            AnalyticThreadOptions {
                maximum_work: WORK + 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_total_work: TOTAL + 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_bytes: BYTES + 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_work: 0,
                ..o
            },
            AnalyticThreadOptions {
                maximum_work: 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_work: used - 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_total_work: head.work_used() - 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_total_work: total - 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_bytes: head.retained_budget_bytes() + SCRATCH_ANALYTIC,
                ..o
            },
            AnalyticThreadOptions {
                maximum_bytes: head.retained_budget_bytes() + SCRATCH_ANALYTIC + 1_000_000,
                ..o
            },
            AnalyticThreadOptions {
                maximum_distance_mm: 1e-12,
                ..o
            },
        ] {
            assert!(head.analytic_thread_repair(bad, None).is_err());
        }
        let exact = head
            .analytic_thread_repair(
                AnalyticThreadOptions {
                    maximum_work: used,
                    maximum_total_work: total,
                    ..o
                },
                None,
            )
            .unwrap();
        assert_eq!(exact.new_work_used(), used);
        let count = Cell::new(0);
        let observe = || {
            count.set(count.get() + 1);
            false
        };
        head.analytic_thread_repair(o, Some(&observe)).unwrap();
        let calls = count.get();
        for frontier in [1, 2, calls / 2, calls] {
            let calls = Cell::new(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() == frontier
            };
            assert!(matches!(
                head.analytic_thread_repair(o, Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ));
        }
        assert_eq!(
            head.analytic_thread_repair(o, None)
                .unwrap()
                .thread()
                .spec(),
            good.thread().spec()
        );
    });
}

#[test]
fn actual_left_handed_posed_packet_qualifies_and_distinct_unsupported_laws_refuse() {
    use crate::{helical_thread::Thread, phillips_drive::PhillipsScrew};
    for (starts, root) in [(1, None), (2, None), (1, Some(0.04))] {
        let original = super::super::super::super::tests::fixture(0, true, true, starts);
        let mut spec = original.thread().spec();
        spec.root_radius_mm = root;
        if let Some(r) = root {
            let a = spec.lower_flank_deg.to_radians();
            let b = spec.upper_flank_deg.to_radians();
            spec.depth_mm = (spec.pitch_mm - spec.crest_width_mm)
                / (a.portable_tan() + b.portable_tan())
                - r * (((b - a) * 0.5).portable_cos() / ((a + b) * 0.5).portable_sin() - 1.);
        }
        let source = PhillipsScrew::new(
            Thread::new(spec).unwrap(),
            *original.profile(),
            original.spec(),
            original.cavity().clone(),
            100_000,
            &|| false,
        )
        .unwrap();
        let text = crate::step::write_phillips_screw(
            &source,
            [10., 20., 30.],
            [0.8, 0.6, 0.],
            [-0.36, 0.48, 0.8],
            1e-7,
        )
        .unwrap();
        super::super::super::common_mouth::tests::derived_thread_probe::with_head(
            &text,
            &mut |head| {
                let result = head.analytic_thread_repair(Default::default(), None);
                if starts == 1 {
                    let receipt =
                        result.unwrap_or_else(|e| panic!("starts{starts} root{root:?}: {e:?}"));
                    assert!(receipt.thread().spec().left_hand);
                    assert_eq!(receipt.frame(), head.body().inferred_thread_frame());
                    assert_eq!(receipt.thread().spec(), head.body().inferred_thread_spec());
                    assert_eq!(
                        receipt.thread().spec().root_radius_mm.is_some(),
                        root.is_some()
                    );
                    for t in receipt.transfers() {
                        assert!(
                            t.distance_mm <= 1e-7 && t.derivative_mm.into_iter().all(|d| d <= 1e-7)
                        );
                        assert!(t.original_cells.into_iter().all(|n| n > 0));
                    }
                    if root.is_some() {
                        assert!(matches!(
                            head.periodic_thread_repair(Default::default(), None),
                            Err(PhillipsReadError::Unsupported)
                        ));
                        assert!(matches!(
                            head.common_chart_periodic_thread_repair(Default::default(), None),
                            Err(PhillipsReadError::Unsupported)
                        ));
                    } else {
                        let periodic = head
                            .periodic_thread_repair(Default::default(), None)
                            .unwrap();
                        assert!(periodic.source().left_hand());
                        let common = head
                            .common_chart_periodic_thread_repair(Default::default(), None)
                            .unwrap();
                        assert_eq!(common.chart(), head.chart());
                        assert!(common.source().left_hand());
                        let finite = head
                            .periodic_thread_finite_curves(&common, Default::default(), None)
                            .unwrap();
                        assert_eq!(finite.curves().len(), 16);
                    }
                } else {
                    assert!(matches!(result, Err(PhillipsReadError::Unsupported)));
                    assert!(matches!(
                        head.periodic_thread_repair(Default::default(), None),
                        Err(PhillipsReadError::Unsupported)
                    ));
                }
            },
        );
    }
}

#[test]
fn actual_m3_flat_and_real_fitted_root_whole_angle_cells_both_hands_and_pose() {
    for fitted in [false, true] {
        for left in [false, true] {
            let text = super::super::degree11::tests::actual_m3(fitted, left);
            super::super::super::common_mouth::tests::derived_thread_probe::with_head(
                &text,
                &mut |head| {
                    let receipt = head
                        .analytic_thread_repair(Default::default(), None)
                        .unwrap();
                    assert_eq!(receipt.thread().spec().root_radius_mm.is_some(), fitted);
                    assert_eq!(receipt.thread().spec(), head.body().inferred_thread_spec());
                    assert!(std::ptr::eq(receipt.body(), head.body()));
                    let mut c0 = 0f64;
                    let mut d = [0f64; 2];
                    let mut cells = 0usize;
                    for t in receipt.transfers() {
                        c0 = c0.max(t.distance_mm);
                        for (a, b) in d.iter_mut().zip(t.derivative_mm) {
                            *a = a.max(b);
                        }
                        cells += t.original_cells[0] * t.original_cells[1];
                        assert!(
                            t.distance_mm <= 1e-7 && t.derivative_mm.into_iter().all(|x| x <= 1e-7)
                        );
                    }
                    eprintln!(
                        "ROUNDED_ANALYTIC fitted={fitted} left={left} work={} inherited={} retained={} originalcells={cells} c0={c0} d={d:?}",
                        receipt.new_work_used(),
                        receipt.inherited_head_work_used(),
                        receipt.retained_budget_bytes()
                    );
                    assert!(receipt.new_work_used() <= WORK && receipt.total_work_used() <= TOTAL);
                    assert_eq!(receipt.body().literal_input(), text);
                    if fitted {
                        assert!(matches!(
                            head.periodic_thread_repair(Default::default(), None),
                            Err(PhillipsReadError::Unsupported)
                        ));
                        assert!(matches!(
                            head.analytic_thread_repair(
                                AnalyticThreadOptions {
                                    maximum_distance_mm: 1e-20,
                                    ..Default::default()
                                },
                                None
                            ),
                            Err(PhillipsReadError::Unsupported)
                        ));
                    }
                },
            );
        }
    }
}

#[test]
fn real_fitted_root_correspondence_preserves_limits_cancel_and_exact_retry() {
    let text = super::super::degree11::tests::actual_m3(true, false);
    super::super::super::common_mouth::tests::derived_thread_probe::with_head(&text, &mut |head| {
        let o = AnalyticThreadOptions::default();
        let good = head.analytic_thread_repair(o, None).unwrap();
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
                    + SCRATCH_ANALYTIC
                    - 1,
                ..o
            },
            AnalyticThreadOptions {
                maximum_cell_derivative_mm: 1e-20,
                ..o
            },
        ] {
            assert!(head.analytic_thread_repair(bad, None).is_err());
        }
        let exact = AnalyticThreadOptions {
            maximum_work: good.new_work_used(),
            maximum_total_work: good.total_work_used(),
            ..o
        };
        let calls = Cell::new(0usize);
        let observed = head
            .analytic_thread_repair(
                exact,
                Some(&|| {
                    calls.set(calls.get() + 1);
                    false
                }),
            )
            .unwrap();
        assert_eq!(observed.new_work_used(), good.new_work_used());
        for stop in [1, 2, calls.get() / 2, calls.get()] {
            let count = Cell::new(0usize);
            assert!(matches!(
                head.analytic_thread_repair(
                    exact,
                    Some(&|| {
                        count.set(count.get() + 1);
                        count.get() == stop
                    })
                ),
                Err(PhillipsReadError::Cancelled)
            ));
        }
        let retry = head.analytic_thread_repair(exact, None).unwrap();
        assert_eq!(retry.thread().spec(), good.thread().spec());
        assert_eq!(retry.new_work_used(), good.new_work_used());
        for (a, b) in retry.transfers().iter().zip(good.transfers()) {
            assert_eq!(
                (
                    a.face_entity,
                    a.surface_entity,
                    a.distance_mm,
                    a.derivative_mm,
                    a.original_cells
                ),
                (
                    b.face_entity,
                    b.surface_entity,
                    b.distance_mm,
                    b.derivative_mm,
                    b.original_cells
                )
            );
        }
    });
}
