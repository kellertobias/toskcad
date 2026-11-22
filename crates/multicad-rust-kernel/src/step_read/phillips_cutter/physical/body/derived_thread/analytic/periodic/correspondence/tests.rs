use super::*;
fn source(spec: crate::helical_thread::Spec) -> PeriodicThreadSource {
    let thread = crate::helical_thread::Thread::new(spec).unwrap();
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    PeriodicThreadSource::new(thread, &mut work).unwrap()
}
#[test]
fn real_metric_periodic_coefficients_and_two_term_caps_compare_over_every_band() {
    let original = crate::step_read::phillips_cutter::physical::tests::fixture(0, false, true, 1);
    for left in [false, true] {
        let spec = crate::helical_thread::Spec {
            left_hand: left,
            ..original.thread().spec()
        };
        let a = source(spec);
        let b = source(crate::helical_thread::Spec {
            pitch_mm: f64::from_bits(spec.pitch_mm.to_bits() + 1),
            ..spec
        });
        let caps = [[0.1, -0.01], [spec.length_mm, 0.]];
        let r = a.finite_correspondence(&b, caps, caps, WORK, None).unwrap();
        assert!(r.position_mm < 1e-7 && r.du_mm < 1e-7 && r.dv_mm < 1e-7);
        eprintln!("periodic finite field source correspondence {r:?}");
        assert!(
            a.finite_correspondence(&a, caps, caps, WORK, None)
                .unwrap()
                .position_mm
                < 1e-7
        );
        let changed_caps = [[0.1, -0.03], [spec.length_mm - 0.1, 0.]];
        let moved = a
            .finite_correspondence(&b, caps, changed_caps, WORK, None)
            .unwrap();
        assert!(moved.position_mm > 0.01 && moved.dv_mm > 0.01);
        for limit in [0, 1, r.work_used - 1] {
            assert!(
                a.finite_correspondence(&b, caps, caps, limit, None)
                    .is_err()
            );
        }
        assert!(
            a.finite_correspondence(&b, caps, caps, r.work_used, None)
                .is_ok()
        );
        for bad in [
            [[f64::NAN, 0.], [1., 0.]],
            [[0., 0.], [f64::INFINITY, 0.]],
            [[0., 0.], [1e8, 0.]],
            [[0., 0.], [0., 0.]],
            [[-1., 0.], [1., 0.]],
            [[0., 0.], [spec.length_mm, 1.]],
        ] {
            assert!(matches!(
                a.finite_correspondence(&b, bad, caps, WORK, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        assert!(matches!(
            a.finite_correspondence(&b, caps, caps, WORK + 1, None),
            Err(PhillipsReadError::Invalid)
        ));
        let opposite = source(crate::helical_thread::Spec {
            left_hand: !left,
            ..spec
        });
        assert!(
            a.finite_correspondence(&opposite, caps, caps, WORK, None)
                .is_err()
        );
        let calls = Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        a.finite_correspondence(&b, caps, caps, WORK, Some(&count))
            .unwrap();
        let total = calls.get();
        for stop in 1..=total {
            calls.set(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() == stop
            };
            assert!(matches!(
                a.finite_correspondence(&b, caps, caps, WORK, Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ));
        }
    }
}
