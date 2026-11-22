use super::*;

fn thread(major: f64, rounded: bool, left: bool, starts: u32) -> Thread {
    let [a, b] = [30_f64.to_radians(), 35_f64.to_radians()];
    Thread::new(Spec {
        major_radius_mm: major,
        pitch_mm: 0.7,
        length_mm: 4.2,
        depth_mm: if rounded {
            (0.7 - 0.1) / (a.portable_tan() + b.portable_tan())
                - 0.04 * (((b - a) * 0.5).portable_cos() / ((b + a) * 0.5).portable_sin() - 1.)
        } else {
            0.35
        },
        crest_width_mm: 0.1,
        lower_flank_deg: 30.,
        upper_flank_deg: 35.,
        starts,
        left_hand: left,
        root_radius_mm: rounded.then_some(0.04),
    })
    .unwrap()
}
#[test]
fn whole_period_bounds_retain_actual_flat_and_rounded_owners_hands_and_starts() {
    for left in [false, true] {
        for starts in [1, 3] {
            for outer_round in [false, true] {
                for inner_round in [false, true] {
                    let outer = thread(2., outer_round, left, starts);
                    let inner = thread(2.2, inner_round, left, starts);
                    let proof = outer
                        .prove_profile_containment(&inner, Default::default(), None)
                        .unwrap();
                    assert!(std::ptr::eq(proof.external(), &outer));
                    assert!(std::ptr::eq(proof.cavity(), &inner));
                    assert!(proof.minimum_gap_mm() > 0.);
                    assert!(proof.cells_visited() > 1);
                    for i in 0..2048 {
                        let phase = i as f64 * outer.spec().pitch_mm / 2048.;
                        let gap = inner.radius_at_phase(phase) - outer.radius_at_phase(phase);
                        assert!(gap > 0. && gap + 1e-12 >= proof.minimum_gap_mm());
                    }
                    let exact = ProfileContainmentOptions {
                        maximum_cells: proof.cells_visited(),
                    };
                    assert_eq!(
                        outer
                            .prove_profile_containment(&inner, exact, None)
                            .unwrap()
                            .minimum_gap_mm(),
                        proof.minimum_gap_mm()
                    );
                    assert_eq!(
                        outer
                            .prove_profile_containment(
                                &inner,
                                ProfileContainmentOptions {
                                    maximum_cells: exact.maximum_cells - 1
                                },
                                None
                            )
                            .unwrap_err(),
                        ProfileContainmentError::Limit
                    );
                }
            }
        }
    }
}
#[test]
fn identical_profiles_earn_contact_without_equating_their_finite_lengths() {
    let outer = thread(2., true, true, 3);
    let inner = Thread::new(Spec {
        length_mm: 2.1,
        ..outer.spec()
    })
    .unwrap();
    let proof = outer
        .prove_profile_containment(&inner, ProfileContainmentOptions { maximum_cells: 1 }, None)
        .unwrap();
    assert_eq!(proof.minimum_gap_mm(), 0.);
    assert_eq!(proof.cells_visited(), 0);
    assert_eq!(proof.external().spec().length_mm, 4.2);
    assert_eq!(proof.cavity().spec().length_mm, 2.1);
}
#[test]
fn wrong_phase_overlap_and_bounded_contact_uncertainty_refuse_without_source_changes() {
    use ProfileContainmentError::*;
    let outer = thread(2., false, false, 1);
    let before = outer;
    for spec in [
        Spec {
            pitch_mm: 0.8,
            ..outer.spec()
        },
        Spec {
            starts: 2,
            ..outer.spec()
        },
        Spec {
            left_hand: true,
            ..outer.spec()
        },
    ] {
        let inner = Thread::new(spec).unwrap();
        assert_eq!(
            outer
                .prove_profile_containment(&inner, Default::default(), None)
                .unwrap_err(),
            Invalid
        );
    }
    for maximum_cells in [0, 1_000_001] {
        assert_eq!(
            outer
                .prove_profile_containment(
                    &outer,
                    ProfileContainmentOptions { maximum_cells },
                    None
                )
                .unwrap_err(),
            Invalid
        );
    }
    let too_small = thread(1.8, false, false, 1);
    assert_eq!(
        outer
            .prove_profile_containment(&too_small, Default::default(), None)
            .unwrap_err(),
        Unresolved
    );
    // A different valid profile shares the crest. It has no strictly positive
    // clearance there; the finite adaptive proof must not manufacture a gap.
    let touching = Thread::new(Spec {
        depth_mm: 0.25,
        ..outer.spec()
    })
    .unwrap();
    assert_eq!(
        outer
            .prove_profile_containment(&touching, Default::default(), None)
            .unwrap_err(),
        Limit
    );
    let base = Thread::new(Spec {
        pitch_mm: 1.,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        depth_mm: 0.375,
        ..outer.spec()
    })
    .unwrap();
    let same_floor = Thread::new(Spec {
        major_radius_mm: 2.125,
        depth_mm: 0.5,
        ..base.spec()
    })
    .unwrap();
    assert_eq!(
        base.prove_profile_containment(&same_floor, Default::default(), None)
            .unwrap_err(),
        Unresolved
    );
    assert_eq!(outer, before);
}
#[test]
fn cancellation_is_polled_through_real_subdivision_and_before_publication() {
    use std::cell::Cell;
    let outer = thread(2., true, false, 1);
    let inner = thread(2.2, false, false, 1);
    for cavity in [&outer, &inner] {
        let calls = Cell::new(0);
        outer
            .prove_profile_containment(
                cavity,
                Default::default(),
                Some(&|| {
                    calls.set(calls.get() + 1);
                    false
                }),
            )
            .unwrap();
        let count = calls.get();
        for stop in 0..count {
            calls.set(0);
            let cancelled = || {
                let n = calls.get();
                calls.set(n + 1);
                n == stop
            };
            assert_eq!(
                outer
                    .prove_profile_containment(cavity, Default::default(), Some(&cancelled))
                    .unwrap_err(),
                ProfileContainmentError::Cancelled
            );
        }
        outer
            .prove_profile_containment(cavity, Default::default(), None)
            .unwrap();
    }
}
#[test]
fn shifted_cavities_are_proved_at_their_actual_phase() {
    use ProfileContainmentError::*;
    for starts in [1, 3] {
        let outer = thread(2., true, false, starts);
        // A deep cavity whose root still clears the outer crest at any phase.
        let loose = thread(2.6, true, false, starts);
        let pitch = outer.spec().pitch_mm;
        for shift in [0., 0.13, pitch / 2., -0.4, 5. * pitch + 0.02] {
            let proof = outer
                .prove_shifted_profile_containment(
                    &loose,
                    [shift - 1e-12, shift + 1e-12],
                    Default::default(),
                    None,
                )
                .unwrap();
            assert!(proof.minimum_gap_mm() > 0.);
            for i in 0..512 {
                let phase = i as f64 * pitch / 512.;
                let gap = loose.radius_at_phase((phase - shift).rem_euclid(pitch))
                    - outer.radius_at_phase(phase);
                assert!(gap + 1e-12 >= proof.minimum_gap_mm(), "{shift} {phase}");
            }
        }
        // A close cavity clears only near its own phase.
        let close = thread(2.2, true, false, starts);
        let zero = outer
            .prove_shifted_profile_containment(&close, [0., 0.], Default::default(), None)
            .unwrap();
        let unshifted = outer
            .prove_profile_containment(&close, Default::default(), None)
            .unwrap();
        assert_eq!(zero.minimum_gap_mm(), unshifted.minimum_gap_mm());
        assert_eq!(
            outer
                .prove_shifted_profile_containment(
                    &close,
                    [pitch / 2., pitch / 2.],
                    Default::default(),
                    None
                )
                .unwrap_err(),
            Unresolved
        );
    }
    let outer = thread(2., true, false, 1);
    let cavity = thread(2.2, true, false, 1);
    let pitch = outer.spec().pitch_mm;
    for bad in [[1., 0.], [f64::NAN, 0.], [0., pitch]] {
        assert_eq!(
            outer
                .prove_shifted_profile_containment(&cavity, bad, Default::default(), None)
                .unwrap_err(),
            Invalid
        );
    }
    let left = thread(2.2, true, true, 1);
    assert_eq!(
        outer
            .prove_shifted_profile_containment(&left, [0., 0.], Default::default(), None)
            .unwrap_err(),
        Invalid
    );
}
