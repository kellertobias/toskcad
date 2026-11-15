use super::*;
#[test]
fn actual_m2_full_fields_bind_the_same_head_cavity_chart_without_reprojection() {
    degree11::tests::with_head(&mut |head| {
        let start = std::time::Instant::now();
        let receipt = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(receipt.body(), head.body()));
        assert_eq!(receipt.chart(), head.chart());
        assert_eq!(
            receipt.chart().columns,
            head.common_mouth()
                .closed_cavity()
                .representative()
                .derived_affine_columns()
        );
        assert_eq!(
            receipt.source().pitch_mm(),
            head.body().inferred_thread_spec().pitch_mm
        );
        assert_eq!(receipt.transfers().len(), 4);
        assert_eq!(receipt.inherited_head_work_used(), head.work_used());
        assert_eq!(
            receipt.total_work_used(),
            receipt.new_work_used() + head.work_used()
        );
        assert!(receipt.retained_budget_bytes() < BYTES);
        for t in receipt.transfers() {
            assert!(t.distance_mm < 1e-7 && t.derivative_mm.into_iter().all(|x| x < 1e-7));
        }
        eprintln!(
            "COMMON PERIODIC {:?} new{} total{} chart{:?} bounds{:?}",
            start.elapsed(),
            receipt.new_work_used(),
            receipt.total_work_used(),
            receipt.chart(),
            receipt.transfers()
        );
    });
}
#[test]
fn actual_common_chart_producer_honors_public_quota_cancel_and_recovery() {
    degree11::tests::with_head(&mut |head| {
        let o = AnalyticThreadOptions::default();
        let good = head.common_chart_periodic_thread_repair(o, None).unwrap();
        let work = good.new_work_used();
        for limit in [0, 1, work - 1] {
            assert!(
                head.common_chart_periodic_thread_repair(
                    AnalyticThreadOptions {
                        maximum_work: limit,
                        ..o
                    },
                    None
                )
                .is_err()
            );
        }
        assert!(
            head.common_chart_periodic_thread_repair(
                AnalyticThreadOptions {
                    maximum_total_work: good.total_work_used() - 1,
                    ..o
                },
                None
            )
            .is_err()
        );
        assert!(
            head.common_chart_periodic_thread_repair(
                AnalyticThreadOptions {
                    maximum_distance_mm: 1e-12,
                    ..o
                },
                None
            )
            .is_err()
        );
        assert!(matches!(
            head.common_chart_periodic_thread_repair(o, Some(&|| true)),
            Err(PhillipsReadError::Cancelled)
        ));
        let again = head.common_chart_periodic_thread_repair(o, None).unwrap();
        assert_eq!(again.chart(), good.chart());
        assert_eq!(again.new_work_used(), work);
    });
}
