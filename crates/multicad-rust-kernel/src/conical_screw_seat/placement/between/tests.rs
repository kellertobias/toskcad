use super::*;
use crate::conical_screw_seat::tests::{head, options, seat};
fn world(frame: Frame, local: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| {
        frame.origin_mm[i]
            + frame.axes[0][i] * local[0]
            + frame.axes[1][i] * local[1]
            + frame.axes[2][i] * local[2]
    })
}
fn pair(frame: Frame, profile: &Profile, seat: Seat) -> Placement {
    let safe = certify(profile, seat, options()).unwrap().safe_neck_z_mm;
    let local = [
        frame.offset_mm[0],
        frame.offset_mm[1],
        frame.offset_mm[2] + safe + 1.,
    ];
    Placement {
        origin_mm: world(frame, local),
        axes: frame.axes,
        axial_offset_mm: -1.,
    }
}
#[test]
fn actual_cardinal_proper_far_and_offset_receiver_charts_keep_physical_bounds() {
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
        let receiver = Frame {
            origin_mm: origin,
            axes,
            offset_mm: [20., -7., 3.],
            source_error_mm: 1e-12,
        };
        let actual_head = pair(receiver, &profile, receiving);
        let copy = receiver;
        let proof =
            certify_between_frames(&profile, receiving, options(), actual_head, receiver).unwrap();
        assert!(proof.maximum_physical_gap_mm <= options().tolerance_mm);
        assert!(proof.minimum_physical_land_mm >= options().minimum_land_mm);
        assert!(
            proof.receiver_operator_norm_upper >= 1. && proof.receiver_inverse_norm_upper >= 1.
        );
        assert_eq!(receiver, copy);
    }
}
#[test]
fn source_uncertainty_shear_and_actual_relative_interference_are_never_normalized_away() {
    let profile = head(3., 90., 0.1);
    let receiving = seat(3., 90.);
    let base = Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        offset_mm: [0.; 3],
        source_error_mm: 0.,
    };
    let actual_head = pair(base, &profile, receiving);
    let mut displaced = actual_head;
    displaced.origin_mm[0] += 0.01;
    assert!(certify_between_frames(&profile, receiving, options(), displaced, base).is_err());
    let mut shear = base;
    shear.axes[1][0] = 1e-9;
    assert!(certify_between_frames(&profile, receiving, options(), actual_head, shear).is_ok());
    shear.axes[1][0] = 1e-7;
    assert!(certify_between_frames(&profile, receiving, options(), actual_head, shear).is_err());
    let transported = pair(shear, &profile, receiving);
    assert!(certify_between_frames(&profile, receiving, options(), transported, shear).is_ok());
    let safe = certify(&profile, receiving, options())
        .unwrap()
        .safe_neck_z_mm;
    assert!(
        certify_between_frames(
            &profile,
            receiving,
            options(),
            actual_head,
            Frame {
                source_error_mm: 1e-10,
                ..base
            }
        )
        .is_ok()
    );
    assert!(
        certify_between_frames(
            &profile,
            receiving,
            options(),
            actual_head,
            Frame {
                source_error_mm: 2e-8,
                ..base
            }
        )
        .is_err()
    );
    assert!(
        certify_between_frames(
            &profile,
            receiving,
            options(),
            actual_head,
            Frame {
                source_error_mm: 1.,
                ..base
            }
        )
        .is_err()
    );
    for error in [-1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            certify_between_frames(
                &profile,
                receiving,
                options(),
                actual_head,
                Frame {
                    source_error_mm: error,
                    ..base
                }
            ),
            Err(Error::Invalid)
        );
    }
    for axes in [
        [[2., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[1., 0., 0.], [0., 1., 0.], [0., 0., -1.]],
        [[1., 0., 0.], [1., 0., 0.], [0., 0., 1.]],
    ] {
        assert!(
            certify_between_frames(
                &profile,
                receiving,
                options(),
                actual_head,
                Frame { axes, ..base }
            )
            .is_err()
        );
    }
    let intruding = Placement {
        origin_mm: [0., 0., safe + 1. + 1e-6],
        ..actual_head
    };
    assert!(certify_between_frames(&profile, receiving, options(), intruding, base).is_err());
    assert!(certify_between_frames(&profile, receiving, options(), actual_head, base).is_ok());
}

#[test]
fn exact_frame_transfer_retains_large_local_offsets_and_simultaneous_proper_transport() {
    let profile = head(6., 90., 0.1);
    let receiving = seat(6., 90.);
    let receiver = Frame {
        origin_mm: [0.; 3],
        axes: [[0.8, 0.6, 0.], [-0.36, 0.48, 0.8], [0.48, -0.64, 0.6]],
        offset_mm: [1e6, -2e6, 3e6],
        source_error_mm: 2e-9,
    };
    let actual_head = pair(receiver, &profile, receiving);
    assert!(certify_between_frames(&profile, receiving, options(), actual_head, receiver).is_ok());
    let q = [[0.6, -0.8, 0.], [0.8, 0.6, 0.], [0., 0., 1.]];
    let rotate =
        |p: [f64; 3]| std::array::from_fn(|i| q[0][i] * p[0] + q[1][i] * p[1] + q[2][i] * p[2]);
    let transport = |p: [f64; 3]| {
        let p = rotate(p);
        [p[0] + 7., p[1] - 9., p[2] + 11.]
    };
    let receiver = Frame {
        origin_mm: transport(receiver.origin_mm),
        axes: receiver.axes.map(rotate),
        ..receiver
    };
    let actual_head = Placement {
        origin_mm: transport(actual_head.origin_mm),
        axes: actual_head.axes.map(rotate),
        ..actual_head
    };
    let receipt =
        certify_between_frames(&profile, receiving, options(), actual_head, receiver).unwrap();
    assert!(receipt.maximum_physical_gap_mm <= options().tolerance_mm);
    assert!(receipt.minimum_physical_land_mm >= options().minimum_land_mm);
    let mut tilt = actual_head;
    tilt.axes[2][0] += 1e-5;
    assert!(certify_between_frames(&profile, receiving, options(), tilt, receiver).is_err());
    assert!(
        certify_between_frames(
            &profile,
            receiving,
            Options {
                tolerance_mm: 1e-10,
                ..options()
            },
            actual_head,
            receiver
        )
        .is_err()
    );
}

#[test]
fn all_six_inward_cardinal_charts_keep_large_exact_offsets_and_source_error() {
    let profile = head(6., 120., 0.1);
    let receiving = seat(6., 120.);
    for axes in [
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[-1., 0., 0.], [0., 1., 0.], [0., 0., -1.]],
        [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]],
        [[0., -1., 0.], [0., 0., 1.], [-1., 0., 0.]],
        [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
        [[0., 0., -1.], [1., 0., 0.], [0., -1., 0.]],
    ] {
        let receiver = Frame {
            origin_mm: [0.; 3],
            axes,
            offset_mm: [1e6, -2e6, 3e6],
            source_error_mm: 2e-9,
        };
        let actual_head = pair(receiver, &profile, receiving);
        let proof =
            certify_between_frames(&profile, receiving, options(), actual_head, receiver).unwrap();
        assert!(proof.maximum_physical_gap_mm <= options().tolerance_mm);
        assert!(proof.minimum_physical_land_mm >= options().minimum_land_mm);
    }
}
#[test]
fn finite_overflow_and_underflow_transfers_refuse_or_enclose_without_poisoning_retry() {
    let profile = head(3., 90., 0.1);
    let receiving = seat(3., 90.);
    let base = Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        offset_mm: [0.; 3],
        source_error_mm: 0.,
    };
    let mut actual_head = pair(base, &profile, receiving);
    actual_head.origin_mm[0] = f64::MAX;
    assert_eq!(
        certify_between_frames(
            &profile,
            receiving,
            options(),
            actual_head,
            Frame {
                origin_mm: [-f64::MAX, 0., 0.],
                ..base
            }
        ),
        Err(Error::Unresolved)
    );
    let overflow = Frame {
        axes: [[1.0000001, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        offset_mm: [f64::MAX, 0., 0.],
        ..base
    };
    assert_eq!(
        certify_between_frames(&profile, receiving, options(), actual_head, overflow),
        Err(Error::Unresolved)
    );
    let tiny = Frame {
        axes: [[0.8, 0.6, 0.], [-0.6, 0.8, 0.], [0., 0., 1.]],
        offset_mm: [f64::from_bits(1), 0., 0.],
        ..base
    };
    assert!(
        certify_between_frames(
            &profile,
            receiving,
            options(),
            pair(tiny, &profile, receiving),
            tiny
        )
        .is_ok()
    );
    let actual_head = pair(base, &profile, receiving);
    assert!(certify_between_frames(&profile, receiving, options(), actual_head, base).is_ok());
}
