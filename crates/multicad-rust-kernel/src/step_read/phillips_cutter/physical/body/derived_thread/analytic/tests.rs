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
            spec.depth_mm = (spec.pitch_mm - spec.crest_width_mm) / (a.tan() + b.tan())
                - r * (((b - a) * 0.5).cos() / ((a + b) * 0.5).sin() - 1.);
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
                if starts == 1 && root.is_none() {
                    let receipt = result.unwrap();
                    assert!(receipt.thread().spec().left_hand);
                    assert_eq!(receipt.frame(), head.body().inferred_thread_frame());
                    let periodic = head
                        .periodic_thread_repair(Default::default(), None)
                        .unwrap();
                    assert!(periodic.source().left_hand());
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
