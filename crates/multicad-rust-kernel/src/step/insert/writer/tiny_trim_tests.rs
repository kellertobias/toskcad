use super::*;
#[test]
fn actual_tiny_caps_earn_exact_finite_direction_without_admitting_axial_overshoot() {
    for left in [false, true] {
        let length = 4_f64.next_up();
        let body = crate::step_read::owned_insert::tests::graph(left, length);
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let cap = t
            .cap_regularity_bounds(16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        assert!(cap.minimum_signed_angular_derivative > 0.);
        assert_eq!(cap.winding, [[-1, 1], [1, -1]]);
        if !left {
            let prerequisite = t
                .winding_bounds(16_000_000, 64 * 1024 * 1024, None)
                .unwrap();
            // A genuine budget shortly beyond the earned winding prefix must
            // stop in finite-derivative validation rather than publish a cap.
            assert!(matches!(
                t.cap_regularity_bounds(prerequisite.work_used + 379, 64 * 1024 * 1024, None),
                Err(Error::Limit)
            ));
            // Enough for winding/interval work, insufficient for the explicitly
            // additional rare-span rational scratch: no uncharged fallback.
            assert!(matches!(
                t.cap_regularity_bounds(
                    16_000_000,
                    t.retained_bytes() + crate::step::insert::exact_cap::SCRATCH,
                    None
                ),
                Err(Error::Limit)
            ));
            t.cap_regularity_bounds(16_000_000, 64 * 1024 * 1024, None)
                .unwrap();
        }
        // A cap-only proof cannot waive the independently observed invalid
        // axial control slab. Full nextafter export remains an explicit gap.
        assert!(
            s.rows()
                .iter()
                .any(|row| row.poles().iter().any(|p| p[2] > length))
        );
        assert!(matches!(
            write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None),
            Err(Error::Unresolved)
        ));
    }
}
