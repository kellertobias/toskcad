use super::*;
use crate::{
    Circle, Rational,
    authored_flat_thread::SourceOptions,
    cylinder::AxialCylinder,
    helical_thread::{Spec, Thread},
    thread_end_relief,
};
fn spec(left: bool, major: f64) -> Spec {
    Spec {
        major_radius_mm: major,
        pitch_mm: 0.5,
        length_mm: 4.,
        depth_mm: 0.25,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    }
}
fn male(s: Spec) -> Chamfer {
    Chamfer::new(
        Thread::new(s).unwrap(),
        0.2,
        Some(thread_end_relief::Intent {
            width_mm: 0.5,
            radial_depth_mm: 0.05,
        }),
    )
    .unwrap()
}
fn female(left: bool, end: Rational) -> Material {
    let authored = AuthoredThread::new(
        spec(left, 1.5),
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let cavity = Cavity::new_with_exact_receiving_extent(
        authored,
        Intent {
            entry_depth_mm: Some(0.4),
            terminal: Some(Terminal {
                width_mm: 0.5,
                radial_depth_mm: 0.05,
            }),
        },
        end,
    )
    .unwrap();
    let zero = Rational::new(0, 1).unwrap();
    Material::new(
        AxialCylinder::new(
            7,
            Circle::new([zero; 2], Rational::new(7, 2).unwrap()).unwrap(),
            zero,
            end,
        )
        .unwrap(),
        cavity,
        1.,
    )
    .unwrap()
}
#[test]
fn whole_actual_profiles_and_both_treatments_earn_positive_gap_in_finite_overlap() {
    for left in [false, true] {
        for end in [
            Rational::new(19, 5).unwrap(),
            Rational::new(4, 1).unwrap(),
            Rational::new(41, 10).unwrap(),
        ] {
            let a = male(spec(left, 1.375));
            let b = female(left, end);
            let before = (a, b);
            let proof = b
                .prove_coaxial_chamfer_separation(&a, Default::default(), None)
                .unwrap();
            assert!(std::ptr::eq(proof.external(), &a));
            assert!(std::ptr::eq(proof.receiving(), &b));
            assert!(proof.reference_profile_gap_mm() > proof.radial_transfer_error_mm());
            assert!(
                proof.radial_transfer_error_mm() > 0. && proof.radial_transfer_error_mm() < 1e-7
            );
            assert!(proof.minimum_radial_gap_mm() > 0. && proof.minimum_radial_gap_mm() <= 0.125);
            assert!(proof.terminal_radial_gap_mm().unwrap() > 0.17);
            let overlap = proof.overlap_end_mm();
            let expected = 4_f64.min(end.as_f64());
            assert!(overlap.lower <= expected && overlap.upper >= expected);
            assert!(proof.cells_visited() > 0);
            // Independently evaluate the finite source sections as supporting
            // witnesses; the whole-profile/cylinder proof above is authority.
            for z in [0., 0.1, 1., 3.49, 3.6, 3.79] {
                for theta in [0., 0.43, 2.1, 6.0] {
                    let radius = a.radius(theta, z);
                    assert!(b.cavity().contains_void(
                        [
                            radius * theta.portable_cos(),
                            radius * theta.portable_sin(),
                            z
                        ],
                        0.
                    ));
                }
            }
            assert_eq!((a, b), before);
        }
    }
}
#[test]
fn changed_actual_laws_overlap_contact_limits_and_cancellation_do_not_earn_separation() {
    let b = female(false, Rational::new(4, 1).unwrap());
    let original = b;
    for changed in [
        Spec {
            pitch_mm: 0.6,
            ..spec(false, 1.375)
        },
        Spec {
            starts: 2,
            ..spec(false, 1.375)
        },
        spec(true, 1.375),
        Spec {
            length_mm: 5.,
            ..spec(false, 1.375)
        },
    ] {
        assert_eq!(
            b.prove_coaxial_chamfer_separation(&male(changed), Default::default(), None)
                .err(),
            Some(Error::Invalid)
        );
    }
    for major in [1.5, 1.625] {
        assert!(
            b.prove_coaxial_chamfer_separation(&male(spec(false, major)), Default::default(), None)
                .is_err()
        );
    }
    let a = male(spec(false, 1.375));
    assert_eq!(
        b.prove_coaxial_chamfer_separation(
            &a,
            ProfileContainmentOptions { maximum_cells: 0 },
            None
        )
        .err(),
        Some(Error::Invalid)
    );
    assert_eq!(
        b.prove_coaxial_chamfer_separation(
            &a,
            ProfileContainmentOptions { maximum_cells: 1 },
            None
        )
        .err(),
        Some(Error::Limit)
    );
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let done = b
        .prove_coaxial_chamfer_separation(&a, Default::default(), Some(&count))
        .unwrap();
    let total = calls.get();
    for stop in [0, total / 2, total - 1] {
        calls.set(0);
        let cancel = || {
            let n = calls.get();
            calls.set(n + 1);
            n >= stop
        };
        assert_eq!(
            b.prove_coaxial_chamfer_separation(&a, Default::default(), Some(&cancel))
                .err(),
            Some(Error::Cancelled)
        );
        assert_eq!(b, original);
    }
    let retry = b
        .prove_coaxial_chamfer_separation(&a, Default::default(), None)
        .unwrap();
    assert_eq!(retry.minimum_radial_gap_mm(), done.minimum_radial_gap_mm());
    assert_eq!(retry.cells_visited(), done.cells_visited());
}

#[test]
fn actual_whole_period_support_retains_literal_overshoot_and_independent_flank_endpoints() {
    let overshoot = Thread::new(Spec {
        major_radius_mm: 1e7,
        pitch_mm: 0.0032990690049907053,
        length_mm: 0.03,
        depth_mm: 0.0007390051778989907,
        crest_width_mm: 0.0007390799746356486,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: Some(0.0014780103557979804),
    })
    .unwrap();
    let actual = overshoot.radius_at_phase(overshoot.profile_breaks()[2]);
    assert!(actual > overshoot.spec().major_radius_mm);
    let support = whole_male_radius(&overshoot).unwrap();
    assert!(support.hi >= actual && support.hi > overshoot.spec().major_radius_mm);
    // This is a genuine public Thread support test, not a forged Chamfer:
    // the existing material constructor correctly refuses that raw source.
    assert_eq!(
        Chamfer::new(overshoot, 0.001, None),
        Err(crate::thread_entry_chamfer::Error::Unresolved)
    );
    let q = |v: f64| num_rational::BigRational::from_float(v).unwrap();
    for thread in [overshoot, Thread::new(spec(false, 1.375)).unwrap()] {
        let bound = whole_male_radius(&thread).unwrap();
        let knots = thread.profile_breaks();
        for band in [1, 3] {
            let anchor = if band == 1 { knots[1] } else { knots[4] };
            let slope = thread.slope_at_phase(anchor, band);
            for end in [knots[band], knots[band + 1]] {
                let exact = q(thread.spec().major_radius_mm) + q(slope) * (q(end) - q(anchor));
                assert!(q(bound.lo) <= exact && exact <= q(bound.hi));
            }
        }
    }
}

#[test]
fn entry_only_actual_receiver_earns_whole_pair_gap_without_a_terminal_cylinder() {
    for left in [false, true] {
        let combined = female(left, Rational::new(4, 1).unwrap());
        let cavity = Cavity::new(
            *combined.cavity().source_thread(),
            Intent {
                entry_depth_mm: Some(0.4),
                terminal: None,
            },
        )
        .unwrap();
        let receiving = Material::new(*combined.stock(), cavity, 1.).unwrap();
        let external = male(spec(left, 1.375));
        let proof = receiving
            .prove_coaxial_chamfer_separation(&external, Default::default(), None)
            .unwrap();
        assert!(proof.minimum_radial_gap_mm() > 0.);
        assert_eq!(proof.terminal_radial_gap_mm(), None);
        assert!(proof.overlap_end_mm().lower <= 4. && proof.overlap_end_mm().upper >= 4.);
        // Supporting actual source sections near the former terminal region;
        // the whole-profile certificate above, not these samples, is authority.
        for theta in [0., 0.43, 2.1, 6.] {
            let r = external.radius(theta, 3.8);
            assert!(receiving.cavity().contains_void(
                [r * theta.portable_cos(), r * theta.portable_sin(), 3.8],
                0.
            ));
        }
    }
}
