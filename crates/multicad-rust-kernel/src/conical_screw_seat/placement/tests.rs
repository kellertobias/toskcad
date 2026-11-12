use super::*;
use crate::conical_screw_seat::tests::{head, options, seat};
fn placed(
    profile: &Profile,
    receiving: Seat,
    origin: [f64; 3],
    axes: [[f64; 3]; 3],
    offset: f64,
) -> (Placement, Axis) {
    let safe = certify(profile, receiving, options())
        .unwrap()
        .safe_neck_z_mm;
    let placement = Placement {
        origin_mm: std::array::from_fn(|i| origin[i] + axes[2][i] * (safe - offset)),
        axes,
        axial_offset_mm: offset,
    };
    (
        placement,
        Axis {
            origin_mm: origin,
            direction: axes[2],
        },
    )
}
#[test]
fn actual_metric_ordinary_rotated_and_far_frames_certify_the_literal_placed_head() {
    let profile = head(3., 90., 0.1);
    let receiving = seat(3., 90.);
    for (origin, axes) in [
        ([0.; 3], [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
        ([17., -3., 21.], [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]]),
        (
            [1e6, -2e6, 3e6],
            [[0.8, 0.6, 0.], [-0.36, 0.48, 0.8], [0.48, -0.64, 0.6]],
        ),
    ] {
        for offset in [0., -1., -7.] {
            let (placement, axis) = placed(&profile, receiving, origin, axes, offset);
            let original = placement;
            let receipt =
                certify_placement(&profile, receiving, options(), placement, axis).unwrap();
            assert!(receipt.maximum_point_deviation_mm < options().tolerance_mm);
            assert!(receipt.bearing.maximum_radial_gap_mm <= options().tolerance_mm);
            assert_eq!(receipt.neck_z_mm, receipt.bearing.certified_neck_z_mm);
            assert_eq!(placement, original);
        }
    }
}
#[test]
fn displacement_tilt_scale_reflection_and_tight_ambiguity_cannot_be_repaired_by_admission() {
    let profile = head(3., 90., 0.1);
    let receiving = seat(3., 90.);
    let (mut placement, axis) = placed(
        &profile,
        receiving,
        [0.; 3],
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        -1.,
    );
    let original = placement;
    placement.origin_mm[0] += 0.01;
    assert!(certify_placement(&profile, receiving, options(), placement, axis).is_err());
    for axes in [
        [[1.001, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[1., 0., 0.], [0., 1., 0.], [0., 0., -1.]],
        [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
    ] {
        placement = original;
        placement.axes = axes;
        assert!(certify_placement(&profile, receiving, options(), placement, axis).is_err());
    }
    let (placement, axis) = placed(
        &profile,
        receiving,
        [1e6, -2e6, 3e6],
        [[0.8, 0.6, 0.], [-0.36, 0.48, 0.8], [0.48, -0.64, 0.6]],
        -1.,
    );
    assert!(
        certify_placement(
            &profile,
            receiving,
            Options {
                tolerance_mm: 1e-10,
                ..options()
            },
            placement,
            axis
        )
        .is_err()
    );
    assert!(certify_placement(&profile, receiving, options(), placement, axis).is_ok());
    assert!(
        certify_placement(
            &profile,
            receiving,
            options(),
            placement,
            Axis {
                direction: [0.; 3],
                ..axis
            }
        )
        .is_err()
    );
    assert!(
        certify_placement(
            &profile,
            receiving,
            options(),
            placement,
            Axis {
                origin_mm: [f64::NAN; 3],
                ..axis
            }
        )
        .is_err()
    );
    let mut invalid = placement;
    invalid.axes[0][0] = f64::INFINITY;
    assert!(certify_placement(&profile, receiving, options(), invalid, axis).is_err());
    assert!(unit([I::point(0.); 3]).is_err());
}

#[test]
fn matching_neck_origin_does_not_hide_whole_head_tilt_or_small_affine_scaling() {
    let profile = head(3., 90., 0.1);
    let receiving = seat(3., 90.);
    let (base, axis) = placed(
        &profile,
        receiving,
        [0.; 3],
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        -1.,
    );
    let safe = certify(&profile, receiving, options())
        .unwrap()
        .safe_neck_z_mm;
    let angle = 1e-6_f64;
    let (sin, cos) = angle.sin_cos();
    let axes = [[cos, 0., -sin], [0., 1., 0.], [sin, 0., cos]];
    let tilted = Placement {
        axes,
        origin_mm: [sin, 0., safe + cos],
        ..base
    };
    assert!((tilted.origin_mm[0] + axes[2][0] * tilted.axial_offset_mm).abs() < 1e-15);
    assert!(frame(axes.map(|a| a.map(I::point))).is_ok());
    assert_eq!(
        certify_placement(&profile, receiving, options(), tilted, axis),
        Err(Error::Placement)
    );
    let mut scaled = base;
    scaled.axes[0][0] = 1. + 1e-8;
    assert!(frame(scaled.axes.map(|a| a.map(I::point))).is_ok());
    assert!(certify_placement(&profile, receiving, options(), scaled, axis).is_err());
    for radius in [3., 6.] {
        for angle in [60., 120.] {
            let profile = head(radius, angle, 0.1);
            let receiving = seat(radius, angle);
            let (placement, axis) = placed(
                &profile,
                receiving,
                [1e6, -2e6, 3e6],
                [[0.8, 0.6, 0.], [-0.36, 0.48, 0.8], [0.48, -0.64, 0.6]],
                -1.,
            );
            assert!(certify_placement(&profile, receiving, options(), placement, axis).is_ok());
        }
    }
}

#[test]
fn axis_normalization_keeps_real_magnitude_independence_and_refuses_numeric_exhaustion() {
    let profile = head(3., 90., 0.1);
    let receiving = seat(3., 90.);
    let (placement, axis) = placed(
        &profile,
        receiving,
        [0.; 3],
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        -1.,
    );
    assert!(
        certify_placement(
            &profile,
            receiving,
            options(),
            placement,
            Axis {
                direction: [0., 0., 1e300],
                ..axis
            }
        )
        .is_ok()
    );
    assert_eq!(
        certify_placement(
            &profile,
            receiving,
            options(),
            placement,
            Axis {
                direction: [0., 0., f64::from_bits(1)],
                ..axis
            }
        ),
        Err(Error::Unresolved)
    );
    assert!(certify_placement(&profile, receiving, options(), placement, axis).is_ok());
}

#[test]
fn nearly_proper_frame_cannot_define_a_coaxial_basis_from_a_parallel_radial_column() {
    let profile = head(3., 90., 0.1);
    let receiving = seat(3., 90.);
    let placement = Placement {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [1e-8, 0., 1.]],
        axial_offset_mm: 0.,
    };
    let axis = Axis {
        origin_mm: [0.; 3],
        direction: [1., 0., 0.],
    };
    assert!(frame(placement.axes.map(|a| a.map(I::point))).is_ok());
    assert_eq!(
        certify_placement(&profile, receiving, options(), placement, axis),
        Err(Error::Placement)
    );
}
