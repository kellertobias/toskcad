use super::*;
#[test]
fn actual_m2_fresh_periodic_source_transfer_retains_shared_node_identity() {
    degree11::tests::with_head(&mut |head| {
        let start = std::time::Instant::now();
        let receipt = head
            .periodic_thread_repair(Default::default(), None)
            .unwrap();
        let source = receipt.source();
        let phase = source.phase_node_expressions_mm();
        let radius = source.radius_nodes_mm();
        assert_eq!(phase[4], [phase[0][0], source.pitch_mm()]);
        assert_eq!(radius[0], radius[4]);
        assert!(radius.into_iter().all(|r| r > 0.));
        assert_eq!(
            source.length_mm(),
            head.body().inferred_thread_spec().length_mm
        );
        assert!(!source.left_hand());
        assert!(std::ptr::eq(receipt.body(), head.body()));
        assert_eq!(receipt.frame(), head.body().inferred_thread_frame());
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
            "PERIODIC M2 {:?} new{} total{} bounds{:?}",
            start.elapsed(),
            receipt.new_work_used(),
            receipt.total_work_used(),
            receipt.transfers()
        );
        for band in 0..4 {
            let a = source.bounds(band, [1.; 2], [0.; 2], WORK, None).unwrap();
            let b = source
                .bounds((band + 1) % 4, [0.; 2], [0.; 2], WORK, None)
                .unwrap();
            // Source equality is the shared-node expression above, not interval overlap.
            // These independently evaluated enclosures must additionally be compatible.
            for k in 0..3 {
                assert!(
                    a.position_mm[k][0] <= b.position_mm[k][1]
                        && b.position_mm[k][0] <= a.position_mm[k][1]
                );
            }
            let whole = source
                .bounds(band, [0., 1.], [0., source.length_mm()], WORK, None)
                .unwrap();
            assert!(whole.position_mm.into_iter().flatten().all(f64::is_finite));
        }
    });
}
#[test]
fn actual_periodic_receipt_and_source_queries_honor_frontiers() {
    degree11::tests::with_head(&mut |head| {
        let o = AnalyticThreadOptions::default();
        let good = head.periodic_thread_repair(o, None).unwrap();
        let used = good.new_work_used();
        for limit in [0, 1, used - 1] {
            assert!(
                head.periodic_thread_repair(
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
            head.periodic_thread_repair(
                AnalyticThreadOptions {
                    maximum_distance_mm: 1e-12,
                    ..o
                },
                None
            )
            .is_err()
        );
        let source = good.source();
        for (band, u, z, budget) in [
            (4, [0.; 2], [0.; 2], WORK),
            (0, [-1., 0.], [0.; 2], WORK),
            (0, [0.; 2], [-1., 0.], WORK),
            (0, [0.; 2], [source.length_mm() + 1.; 2], WORK),
            (0, [0.; 2], [0.; 2], WORK + 1),
            (0, [f64::NAN; 2], [0.; 2], WORK),
            (0, [0.; 2], [0.; 2], 0),
            (0, [0.; 2], [0.; 2], 256),
            (0, [0.; 2], [0.; 2], 4352),
            (0, [0.; 2], [0.; 2], 8448),
        ] {
            assert!(source.bounds(band, u, z, budget, None).is_err());
        }
        assert!(matches!(
            source.bounds(0, [0.; 2], [0.; 2], WORK, Some(&|| true)),
            Err(PhillipsReadError::Cancelled)
        ));
        let count = Cell::new(0);
        let observe = || {
            count.set(count.get() + 1);
            false
        };
        head.periodic_thread_repair(o, Some(&observe)).unwrap();
        let total = count.get();
        for frontier in [1, 2, total / 2, total] {
            let n = Cell::new(0);
            let cancel = || {
                n.set(n.get() + 1);
                n.get() == frontier
            };
            assert!(matches!(
                head.periodic_thread_repair(o, Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ));
        }
    });
}
#[test]
fn independent_axis_point_and_one_sided_radial_derivative_match_shared_nodes() {
    degree11::tests::with_head(&mut |head| {
        let receipt = head
            .periodic_thread_repair(Default::default(), None)
            .unwrap();
        let source = receipt.source();
        let phase = source.phase_node_expressions_mm();
        let radius = source.radius_nodes_mm();
        for band in 0..4 {
            let z = phase[band][0];
            let b = source.bounds(band, [0.; 2], [z; 2], WORK, None).unwrap();
            let contains = |r: [f64; 2], v: f64| r[0] <= v && r[1] >= v;
            assert!(contains(b.position_mm[0], radius[band]));
            assert!(contains(b.position_mm[1], 0.));
            assert!(contains(b.position_mm[2], z));
            assert!(contains(
                b.derivative_mm[0][0],
                radius[band + 1] - radius[band]
            ));
            assert!(contains(b.derivative_mm[0][2], 0.));
            assert!(contains(b.derivative_mm[1][0], 0.));
            assert!(contains(b.derivative_mm[1][2], 1.));
        }
    });
}
#[test]
fn legal_nearly_vanishing_radial_band_can_honestly_refuse_interval_ordering() {
    use crate::helical_thread::{Spec, Thread};
    let mut depth = (1. - 0.1) / (2. * 30_f64.to_radians().tan());
    let mut refused = false;
    for _ in 0..16 {
        depth = depth.next_down();
        if let Ok(thread) = Thread::new(Spec {
            major_radius_mm: 2.,
            pitch_mm: 1.,
            length_mm: 1.,
            depth_mm: depth,
            crest_width_mm: 0.1,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        }) {
            let mut work = Work {
                remaining: WORK,
                interrupt: None,
                cancelled: Cell::new(false),
            };
            refused |= PeriodicThreadSource::new(thread, &mut work).is_err();
        }
    }
    assert!(refused);
}
