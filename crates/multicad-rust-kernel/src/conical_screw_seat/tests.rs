use super::*;
use crate::head_profile::Head;
fn head(radius: f64, angle: f64, neck: f64) -> Profile {
    let shank = radius * 0.5;
    Profile::new(
        Head::Countersunk {
            radius_mm: radius,
            height_mm: (radius - shank) / (angle.to_radians() * 0.5).tan() + 0.4,
            included_angle_deg: angle,
            edge_chamfer_mm: 0.1,
            neck_radius_mm: neck,
        },
        shank,
        1.,
    )
    .unwrap()
}
fn seat(radius: f64, angle: f64) -> Seat {
    let bore = radius * 0.5 + 0.2;
    Seat {
        bore_radius_mm: bore,
        mouth_radius_mm: radius,
        depth_mm: (radius - bore) / (angle.to_radians() * 0.5).tan(),
    }
}
fn options() -> Options {
    Options {
        tolerance_mm: 1e-7,
        minimum_land_mm: 0.1,
    }
}
#[test]
fn actual_metric_profiles_derive_safe_neck_and_positive_radial_bearing_without_recipe_alignment() {
    for radius in [3., 6.] {
        for angle in [60., 82., 90., 100., 120.] {
            let profile = head(radius, angle, 0.1);
            let original = profile;
            let seat = seat(radius, angle);
            let receipt = certify(&profile, seat, options()).unwrap();
            assert!(receipt.safe_neck_z_mm < receipt.ideal_neck_z_mm[0]);
            assert!(receipt.bearing_axial_mm[1] > receipt.bearing_axial_mm[0]);
            assert!(receipt.bearing_radial_mm[1] - receipt.bearing_radial_mm[0] >= 0.1);
            assert!(receipt.maximum_radial_gap_mm > 0. && receipt.maximum_radial_gap_mm <= 1e-7);
            for segment in profile.segments() {
                for i in 0..101 {
                    let [r, z] = segment.point(f64::from(i) / 100.);
                    let z = z + receipt.safe_neck_z_mm;
                    if z >= 0. {
                        let cavity = seat.bore_radius_mm.max(
                            seat.mouth_radius_mm
                                - (seat.mouth_radius_mm - seat.bore_radius_mm) / seat.depth_mm * z,
                        );
                        assert!(r <= cavity, "literal interference r={r} allowed={cavity}");
                    }
                }
            }
            assert_eq!(profile, original);
        }
    }
}
#[test]
fn mismatch_thin_bearing_invalid_precision_and_neck_interference_refuse_without_mutation() {
    let profile = head(3., 90., 0.1);
    let receiving = seat(3., 90.);
    let original = profile;
    assert_eq!(
        certify(
            &profile,
            Seat {
                depth_mm: receiving.depth_mm * 1.01,
                ..receiving
            },
            options()
        ),
        Err(Error::IncompatibleCone)
    );
    assert_eq!(
        certify(
            &profile,
            receiving,
            Options {
                minimum_land_mm: 2.,
                ..options()
            }
        ),
        Err(Error::InsufficientLand)
    );
    for tolerance in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            certify(
                &profile,
                receiving,
                Options {
                    tolerance_mm: tolerance,
                    ..options()
                }
            ),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        certify(
            &profile,
            receiving,
            Options {
                tolerance_mm: 1e-30,
                ..options()
            }
        ),
        Err(Error::Unresolved)
    );
    for bad in [
        Seat {
            bore_radius_mm: 0.,
            ..receiving
        },
        Seat {
            depth_mm: 0.,
            ..receiving
        },
        Seat {
            mouth_radius_mm: receiving.bore_radius_mm,
            ..receiving
        },
    ] {
        assert_eq!(certify(&profile, bad, options()), Err(Error::Invalid));
    }
    let tight = Seat {
        bore_radius_mm: 1.501,
        mouth_radius_mm: 3.,
        depth_mm: 1.499,
    };
    assert!(certify(&profile, tight, options()).is_err());
    assert_eq!(profile, original);
    assert!(certify(&profile, receiving, options()).is_ok());
    let flat = Profile::new(
        Head::Crowned {
            radius_mm: 3.,
            height_mm: 2.,
            crown_radius_mm: 4.,
            top_land_radius_mm: 0.5,
            neck_radius_mm: 0.1,
        },
        1.5,
        1.,
    )
    .unwrap();
    assert_eq!(certify(&flat, receiving, options()), Err(Error::NoCone));
}

#[test]
fn finite_source_and_clearance_proof_controls_never_promote_unresolved_arithmetic() {
    assert_eq!(
        checked(I::new(f64::NEG_INFINITY, f64::INFINITY)).unwrap_err(),
        Error::Unresolved
    );
    assert!(trig(I::point(2.1)).is_err());
    assert!(trig(I::point(f64::NAN)).is_err());
    let profile = head(3., 90., 0.1);
    let cone = profile.segments().find(|s| s.role == Role::Cone).unwrap();
    assert!(arc_bounds(cone).is_err());
    for bad in [
        Segment {
            curve: Curve::Arc {
                center: [0., 0.],
                sweep: 0.1,
            },
            ..cone
        },
        Segment {
            end: cone.start,
            ..cone
        },
        Segment {
            end: [cone.start[0] + 1., cone.end[1]],
            ..cone
        },
    ] {
        assert!(supporting_cone(&[bad]).is_err());
    }
    assert!(supporting_cone(&[cone, cone]).is_err());
    let receiving = seat(3., 90.);
    let slope = I::point(1.);
    assert!(
        clear_segment(
            Segment {
                start: [1., 1.],
                end: [1., 0.],
                curve: Curve::Line,
                role: Role::Shank
            },
            0.,
            0.,
            receiving,
            slope
        )
        .is_err()
    );
    let neck = profile.segments().find(|s| s.role == Role::Neck).unwrap();
    assert!(clear_segment(neck, -100., profile.rounding_error_mm(), receiving, slope).is_ok());
    let interferes = Segment {
        start: [4., 0.],
        end: [4., 1.],
        curve: Curve::Line,
        role: Role::Rim,
    };
    assert_eq!(
        clear_segment(interferes, 0., 0., receiving, slope),
        Err(Error::Interference)
    );
    let touching = Segment {
        start: [1.7, 0.],
        end: [1.7, 1.],
        curve: Curve::Line,
        role: Role::Shank,
    };
    assert_eq!(
        clear_segment(touching, 100., 0., receiving, slope),
        Err(Error::Unresolved)
    );
    let no_bearing = Seat {
        bore_radius_mm: 1.3,
        mouth_radius_mm: 1.4,
        depth_mm: 0.1,
    };
    assert_eq!(
        certify(&profile, no_bearing, options()),
        Err(Error::InsufficientLand)
    );
    let tiny_budget = Options {
        tolerance_mm: profile.rounding_error_mm() * 2.,
        ..options()
    };
    assert_eq!(
        certify(&profile, receiving, tiny_budget),
        Err(Error::Unresolved)
    );
    assert!(certify(&head(3., 90., 0.), receiving, options()).is_ok());
}

#[test]
fn finite_but_unresolvable_receiving_dimensions_never_become_a_supported_seat() {
    let profile = head(3., 90., 0.1);
    let impossible = Seat {
        bore_radius_mm: 1.,
        mouth_radius_mm: 1e308,
        depth_mm: 1e-308,
    };
    assert_eq!(
        certify(&profile, impossible, options()),
        Err(Error::Unresolved)
    );
    assert!(certify(&profile, seat(3., 90.), options()).is_ok());
}

#[test]
fn retained_actual_rim_and_chamfer_material_are_checked_beyond_the_supporting_cone() {
    let line = |start, end, role| Segment {
        start,
        end,
        curve: Curve::Line,
        role,
    };
    let rim = [
        line([5., -3.], [5., -1.], Role::Rim),
        line([5., -1.], [3., -1.], Role::Bearing),
        line([3., -1.], [1.5, 0.5], Role::Cone),
        line([1.5, 0.5], [1.5, 1.], Role::Shank),
    ];
    let profile = Profile::from_retained_segments(&rim).unwrap();
    let narrow = Seat {
        bore_radius_mm: 1.7,
        mouth_radius_mm: 3.,
        depth_mm: 1.3,
    };
    assert!(
        certify(&profile, narrow, options()).is_ok(),
        "upper actual rim is outside receiving stock"
    );
    let wide = Seat {
        bore_radius_mm: 1.7,
        mouth_radius_mm: 5.,
        depth_mm: 3.3,
    };
    assert_eq!(certify(&profile, wide, options()), Err(Error::Interference));
    let chamfer = [
        line([4., -3.], [5., -2.], Role::Chamfer),
        line([5., -2.], [3., -2.], Role::Bearing),
        line([3., -2.], [1.5, -0.5], Role::Cone),
        line([1.5, -0.5], [1.5, 1.], Role::Shank),
    ];
    let profile = Profile::from_retained_segments(&chamfer).unwrap();
    assert_eq!(certify(&profile, wide, options()), Err(Error::Interference));
    let bearing = line([1.6, 0.], [1.5, 0.], Role::Bearing);
    assert!(clear_segment(bearing, 1., 0., narrow, I::point(1.)).is_ok());
    let original = head(3., 90., 0.3);
    let retained =
        Profile::from_retained_segments(&original.segments().collect::<Vec<_>>()).unwrap();
    assert!(
        certify(&retained, narrow, options()).is_ok(),
        "finite neck support must avoid full-circle false refusal"
    );
}

#[test]
fn rounded_final_poses_require_whole_head_clearance_and_a_bounded_real_gap() {
    let profile = head(3., 90., 0.1);
    let receiving = seat(3., 90.);
    let proposed = certify(&profile, receiving, options()).unwrap();
    let safe = proposed.safe_neck_z_mm;
    assert!(proposed.maximum_radial_gap_mm >= options().tolerance_mm / 4.);
    let interval = [safe - 1e-10, safe + 1e-10];
    let actual = certify_at_neck(&profile, receiving, options(), interval, 1e-10).unwrap();
    assert_eq!(actual.certified_neck_z_mm, interval);
    assert_eq!(actual.safe_neck_z_mm, safe);
    assert_eq!(actual.placement_deviation_mm, 1e-10);
    assert!(actual.maximum_radial_gap_mm <= options().tolerance_mm);
    assert!(certify_at_neck(&profile, receiving, options(), proposed.ideal_neck_z_mm, 0.).is_err());
    assert_eq!(
        certify_at_neck(&profile, receiving, options(), [safe + 1e-6; 2], 0.),
        Err(Error::Interference)
    );
    assert_eq!(
        certify_at_neck(&profile, receiving, options(), [safe - 1e-6; 2], 0.),
        Err(Error::Unresolved)
    );
    for interval in [[safe + 1., safe], [f64::NAN, safe], [safe, f64::INFINITY]] {
        assert_eq!(
            certify_at_neck(&profile, receiving, options(), interval, 0.),
            Err(Error::Invalid)
        );
    }
    for extra in [-1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            certify_at_neck(&profile, receiving, options(), [safe; 2], extra),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        certify_at_neck(&profile, receiving, options(), [safe; 2], 1.),
        Err(Error::Unresolved)
    );
    assert!(certify_at_neck(&profile, receiving, options(), [safe; 2], 0.).is_ok());
}
