use super::*;
#[test]
fn proposal_ray_rotation_has_a_no_libm_rational_oracle_on_the_whole_angle_domain() {
    use num_bigint::BigInt as Z;
    use num_rational::BigRational as Q;
    for angle in [-7_f64, -1. / 3., 0., 1. / 3., 7.] {
        let x = Q::from_float(angle).unwrap();
        let mut terms = [Q::from_integer(Z::from(1)), Q::from_integer(Z::from(0))];
        let mut sums = terms.clone();
        for k in 1..=127 {
            terms = [-&terms[1] * &x / Z::from(k), &terms[0] * &x / Z::from(k)];
            for axis in 0..2 {
                sums[axis] += &terms[axis];
            }
        }
        // Independent exact degree127 sum and Lagrange remainder.
        let mut error = Q::from_integer(Z::from(1));
        for k in 1..=128 {
            error = error * Q::from_float(angle.abs()).unwrap() / Z::from(k);
        }
        let got = rotation(I::point(angle)).unwrap();
        for axis in 0..2 {
            assert!(Q::from_float(got[axis].lo).unwrap() <= &sums[axis] - &error);
            assert!(Q::from_float(got[axis].hi).unwrap() >= &sums[axis] + &error);
        }
    }
    for angle in [I::point(f64::NAN), I::point(7.1)] {
        assert!(matches!(
            rotation(angle),
            Err(RepairedHeadError::Unresolved)
        ));
    }
}
#[test]
fn actual_crowned_arcs_have_a_whole_correspondence_with_small_c0_c1_and_normal_error() {
    let p = Profile::new_for_thread_union(
        crate::head_profile::Head::Crowned {
            radius_mm: 2.,
            height_mm: 1.6,
            crown_radius_mm: 3.2,
            top_land_radius_mm: 1.25,
            neck_radius_mm: 0.1,
        },
        1.,
        0.,
    )
    .unwrap();
    let source = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let proof = source.proposal_transfer(WORK, None).unwrap();
    assert!(proof.position_mm > 0. && proof.position_mm < 1e-11);
    assert!(proof.derivative_mm > 0. && proof.derivative_mm < 1e-10);
    assert!(proof.normal_angle_radians > 0. && proof.normal_angle_radians < 1e-8);
    assert_eq!(proof.work_used, 1 + 4096 * source.segment_count());
    assert!(matches!(
        source.proposal_transfer(WORK + 1, None),
        Err(RepairedHeadError::Invalid)
    ));
    assert!(matches!(
        source.proposal_transfer(proof.work_used - 1, None),
        Err(RepairedHeadError::Limit)
    ));
    let calls = std::cell::Cell::new(0);
    for target in 0..=source.segment_count() + 1 {
        calls.set(0);
        let cancel = || {
            let n = calls.get();
            calls.set(n + 1);
            n == target
        };
        assert!(matches!(
            source.proposal_transfer(WORK, Some(&cancel)),
            Err(RepairedHeadError::Cancelled)
        ));
    }
    // Independent finite differences are only a regression oracle. The receipt
    // above comes from complete rays, determinants and cubic coefficients.
    for i in 0..source.segment_count() {
        let c = source.segment_certificate(i).unwrap();
        assert!(proof.position_mm >= c.support_transfer_mm);
    }
}
