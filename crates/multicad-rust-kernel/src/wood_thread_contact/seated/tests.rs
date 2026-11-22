use super::*;
use crate::{
    axial_profile_stack::{ProfileStack, RectangularBlindCounterbore},
    headed_screw::ScrewSpec,
    helical_thread::{Spec, Thread},
};
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn dimensions() -> RectangularBlindCounterbore {
    RectangularBlindCounterbore {
        min_mm: [-10.; 2],
        max_mm: [10.; 2],
        height_mm: 20.,
        bore_radius_mm: 2.,
        bore_depth_mm: 16.,
        seat_radius_mm: 6.,
        seat_depth_mm: 1.,
        minimum_wall_mm: 1.,
        minimum_floor_mm: 1.,
    }
}
fn req() -> Requirements {
    Requirements {
        minimum_wall_mm: 1.,
        minimum_floor_mm: 1.,
        minimum_engagement_mm: 12.,
        tip_clearance_mm: 1.,
        linear_tolerance_mm: 1e-7,
        minimum_feature_mm: 1e-4,
    }
}
fn screw(left: bool, starts: u32) -> Screw {
    Screw::new(
        Thread::new(Spec {
            major_radius_mm: 3.,
            pitch_mm: 2.,
            length_mm: 12.,
            depth_mm: 1.,
            crest_width_mm: 0.4,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts,
            left_hand: left,
            root_radius_mm: None,
        })
        .unwrap(),
        ScrewSpec {
            head: ScrewHead::ExternalHex {
                across_flats_mm: 10.,
                height_mm: 4.,
            },
            shank_length_mm: 0.,
            shank_radius_mm: 3.,
        },
    )
    .unwrap()
}
fn frame() -> Frame {
    Frame {
        origin_mm: [0., 0., 1.],
        ..identity()
    }
}
fn run(
    s: RectangularBlindCounterbore,
    screw: &Screw,
    f: Frame,
    r: Requirements,
    c: &dyn Fn() -> bool,
) -> R<SeatedContact> {
    let stack = ProfileStack::blind_counterbore(s, 1e-7, 1e-4).unwrap();
    certify_seated(
        stack.rectangular_blind_counterbore_source().unwrap(),
        identity(),
        screw,
        f,
        r,
        c,
    )
}
#[test]
fn original_custom_thread_and_actual_shoulder_prove_complete_engagement() {
    for left in [false, true] {
        for starts in [1, 2, 8] {
            let screw = screw(left, starts);
            let proof = run(dimensions(), &screw, frame(), req(), &|| false).unwrap();
            assert!(proof.contact.engagement_mm[0] <= 12. && proof.contact.engagement_mm[1] >= 12.);
            assert!(
                proof.contact.tip_clearance_mm[0] <= 3. && proof.contact.tip_clearance_mm[1] >= 3.
            );
            assert_eq!(proof.seat_depth_mm, 1.);
            assert_eq!(proof.bore_depth_mm, 16.);
            assert_eq!(
                screw
                    .thread()
                    .classify_box(proof.contact.thread_witness_mm)
                    .unwrap(),
                BoxClassification::Inside
            );
            let b = proof.contact.thread_witness_mm;
            let radial = I::new(b[0][0], b[0][1]).square() + I::new(b[1][0], b[1][1]).square();
            assert!(radial.lo > 4. && b[2][0] + 1. > 1. && b[2][1] + 1. < 16.);
            assert!(proof.contact.witness_volume_lower_mm3 > 0.);
            assert!(
                proof.contact.source_chart_bearing_radial_land_mm[1]
                    > proof.contact.source_chart_bearing_radial_land_mm[0]
            );
        }
    }
    let screw = screw(false, 1);
    let stack = ProfileStack::blind_counterbore(dimensions(), 1e-7, 1e-4).unwrap();
    let stock = Frame {
        origin_mm: [0.; 3],
        axes: [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]],
    };
    let placed = Frame {
        origin_mm: [1., 0., 0.],
        ..stock
    };
    assert!(
        certify_seated(
            stack.rectangular_blind_counterbore_source().unwrap(),
            stock,
            &screw,
            placed,
            req(),
            &|| false
        )
        .is_ok()
    );
}
#[test]
fn actual_head_plain_core_tip_and_policy_failures_do_not_borrow_the_seat() {
    let s = screw(false, 1);
    for n in 0..9 {
        let mut d = dimensions();
        let mut f = frame();
        let mut r = req();
        match n {
            0 => f.origin_mm[2] = 0.,
            1 => d.seat_radius_mm = 5.,
            2 => d.bore_radius_mm = 1.9,
            3 => d.bore_radius_mm = 3.,
            4 => r.minimum_engagement_mm = 12_f64.next_up(),
            5 => r.tip_clearance_mm = 3_f64.next_up(),
            6 => r.minimum_wall_mm = 5.,
            7 => r.minimum_floor_mm = 5.,
            _ => f.axes = [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
        }
        assert!(run(d, &s, f, r, &|| false).is_err(), "case {n}");
    }
    let mut spec = s.spec();
    spec.shank_length_mm = 1.;
    let plain = Screw::new(*s.thread(), spec).unwrap();
    let mut f = frame();
    f.origin_mm[2] = 2.;
    assert!(run(dimensions(), &plain, f, req(), &|| false).is_err());
    let mut invalid = req();
    invalid.minimum_engagement_mm = f64::NAN;
    assert!(matches!(
        run(dimensions(), &s, frame(), invalid, &|| false),
        Err(Error::Invalid)
    ));
    let mut r = req();
    r.minimum_feature_mm = 3.;
    assert!(run(dimensions(), &s, frame(), r, &|| false).is_err());
}
#[test]
fn every_actual_cancellation_checkpoint_preserves_atomic_publication() {
    use std::cell::Cell;
    let s = screw(false, 1);
    let calls = Cell::new(0);
    assert!(
        run(dimensions(), &s, frame(), req(), &|| {
            calls.set(calls.get() + 1);
            false
        })
        .is_ok()
    );
    for stop in 0..calls.get() {
        let at = Cell::new(0);
        let result = run(dimensions(), &s, frame(), req(), &|| {
            let n = at.get();
            at.set(n + 1);
            n == stop
        });
        assert!(matches!(result, Err(Error::Cancelled)));
    }
}

#[test]
fn real_source_precision_limits_and_frame_orientation_remain_explicit() {
    let original = screw(false, 1);
    let rounded = Thread::new(Spec {
        root_radius_mm: Some(0.1),
        depth_mm: (2. - 0.4) / (2. * 30_f64.to_radians().portable_tan())
            - 0.1 * (1. / 30_f64.to_radians().portable_sin() - 1.),
        ..original.thread().spec()
    })
    .unwrap();
    let rounded = Screw::new(rounded, original.spec()).unwrap();
    assert!(matches!(
        run(dimensions(), &rounded, frame(), req(), &|| false),
        Err(Error::Unsupported)
    ));
    let mut tiny = dimensions();
    tiny.bore_radius_mm = 3_f64.next_down();
    assert!(matches!(
        run(tiny, &original, frame(), req(), &|| false),
        Err(Error::Unsupported)
    ));
    // A real eccentric placement can defeat this bounded witness proposal.
    // It is refused rather than counting a thread box still inside the pilot.
    let mut eccentric = dimensions();
    eccentric.bore_radius_mm = 2.8;
    eccentric.seat_radius_mm = 8.;
    let mut moved = frame();
    moved.origin_mm[0] = -0.5;
    assert!(matches!(
        run(eccentric, &original, moved, req(), &|| false),
        Err(Error::Unsupported)
    ));
    let socket = Screw::new(
        *original.thread(),
        ScrewSpec {
            head: ScrewHead::SocketCap {
                radius_mm: 5.,
                height_mm: 4.,
                socket_across_flats_mm: 3.,
                socket_depth_mm: 2.,
            },
            ..original.spec()
        },
    )
    .unwrap();
    assert!(run(dimensions(), &socket, frame(), req(), &|| false).is_ok());
    let stack = ProfileStack::blind_counterbore(dimensions(), 1e-7, 1e-4).unwrap();
    for stock_frame in [
        Frame {
            axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., -1.]],
            ..identity()
        },
        Frame {
            axes: [[1., 0., 0.], [0.01, 1., 0.], [0., 0., 1.]],
            ..identity()
        },
        Frame {
            origin_mm: [f64::NAN, 0., 0.],
            ..identity()
        },
    ] {
        assert!(
            certify_seated(
                stack.rectangular_blind_counterbore_source().unwrap(),
                stock_frame,
                &original,
                frame(),
                req(),
                &|| false
            )
            .is_err()
        );
    }
    let stock_frame = Frame {
        axes: [[0.8, 0.6, 0.], [-0.6, 0.8, 0.], [0., 0., 1.]],
        ..identity()
    };
    let screw_frame = Frame {
        origin_mm: [0., 0., 1.],
        ..stock_frame
    };
    assert!(
        certify_seated(
            stack.rectangular_blind_counterbore_source().unwrap(),
            stock_frame,
            &original,
            screw_frame,
            req(),
            &|| false
        )
        .is_ok()
    );
    let mut touching = dimensions();
    touching.bore_depth_mm = 13.;
    let mut zero_gap = req();
    zero_gap.tip_clearance_mm = 0.;
    zero_gap.minimum_floor_mm = 7.;
    let proof = run(touching, &original, frame(), zero_gap, &|| false).unwrap();
    assert_eq!(proof.contact.tip_clearance_mm[0], 0.);
    assert!(proof.contact.tip_clearance_mm[1] >= 0.);
    let mut turned = frame();
    turned.axes = [[0.8, 0.6, 0.], [-0.6, 0.8, 0.], [0., 0., 1.]];
    assert!(run(dimensions(), &original, turned, req(), &|| false).is_err());
}
