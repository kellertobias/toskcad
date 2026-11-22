use super::*;
use crate::HelicalThread;
use multicad_core::references::Origin;
use multicad_rust_kernel::headed_screw::{ScrewHead, ScrewSpec};
pub(in crate::geometry::wood_contact) fn fixture(length: f64) -> (Solid, Solid) {
    let policy = GeometryPolicy::default();
    let stock = Solid::box_solid([20.; 3], policy).unwrap();
    let stock = stock
        .rust_connection_blind_bore(
            &CircularPocketTool {
                frame: PlaneFrame {
                    origin: [10., 10., 20.001],
                    u: [1., 0., 0.],
                    v: [0., -1., 0.],
                },
                center_uv: [0.; 2],
                radius_mm: 4.,
                displacement_mm: [0., 0., -6.001],
                origin: Origin {
                    source: EntityId::new(),
                    role: "actual.stock.pilot".into(),
                },
            },
            1.,
            1.,
            policy,
        )
        .unwrap()
        .0;
    let screw = Solid::rust_machine_screw_at_thread_frame(
        HelicalThread {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: length,
            depth_mm: 0.25,
            crest_width_mm: 0.05,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: false,
        },
        None,
        ScrewSpec {
            head: ScrewHead::SocketCap {
                radius_mm: 2.75,
                height_mm: 3.,
                socket_across_flats_mm: 2.,
                socket_depth_mm: 1.,
            },
            shank_length_mm: 0.,
            shank_radius_mm: 1.5,
        },
        EntityId::new(),
        PlaneFrame {
            origin: [10., 10., 20.],
            u: [1., 0., 0.],
            v: [0., -1., 0.],
        },
        policy,
    )
    .unwrap();
    (stock, screw)
}
#[test]
fn real_five_point_five_mm_tip_fails_real_one_mm_blind_bottom_requirement() {
    for (length, gap) in [(3., 3.), (4.5, 1.5), (5., 1.), (5.5, 0.5)] {
        let (stock, screw) = fixture(length);
        let proof = stock
            .rust_wood_insert_screw_clearance(&screw, 0.4, None)
            .unwrap();
        assert!(proof.tip_clearance_mm[0] <= gap && proof.tip_clearance_mm[1] >= gap);
        if length == 5.5 {
            let e = stock
                .rust_wood_insert_screw_clearance(&screw, 1., None)
                .unwrap_err();
            assert_eq!(e.kind, GeometryErrorKind::Unsupported);
            assert!(e.to_string().contains("blind-bottom clearance"));
        } else {
            assert!(
                stock
                    .rust_wood_insert_screw_clearance(&screw, 1., None)
                    .is_ok()
            );
        }
    }
}
#[test]
fn cancelled_unknown_and_mismatched_policy_queries_publish_no_stock_exclusion() {
    let (stock, screw) = fixture(4.5);
    assert_eq!(
        stock
            .rust_wood_insert_screw_clearance(&screw, 1., Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert_eq!(
        stock
            .rust_wood_insert_screw_clearance(&screw, -1., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    let plain = Solid::box_solid([20.; 3], GeometryPolicy::default()).unwrap();
    assert_eq!(
        plain
            .rust_wood_insert_screw_clearance(&screw, 1., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        screw
            .rust_wood_insert_screw_clearance(&screw, 1., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        stock
            .rust_wood_insert_screw_clearance(&plain, 1., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let mut policy = GeometryPolicy::default();
    policy.linear_tolerance_mm *= 2.;
    let unequal = Solid::box_solid([20.; 3], policy).unwrap();
    assert_eq!(
        stock
            .rust_wood_insert_screw_clearance(&unequal, 1., None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    assert!(
        stock
            .rust_wood_insert_screw_clearance(&screw, 1., None)
            .is_ok()
    );
}

#[test]
fn actual_periodic_phillips_continuation_preserves_blind_floor_and_radial_clearance() {
    for (length, gap) in [(3., 3.), (4.5, 1.5), (5., 1.), (5.5, 0.5)] {
        let (stock, regular) = fixture(length);
        let (original, frame) = regular.rust_machine_screw_geometry().unwrap();
        let screw = super::super::mating::tests::periodic_screw(*original.thread(), 1., frame);
        let receipt = stock
            .rust_wood_insert_screw_clearance(&screw, 0.4, None)
            .unwrap();
        assert!(receipt.tip_clearance_mm[0] <= gap && receipt.tip_clearance_mm[1] >= gap);
        assert!(receipt.minimum_radial_clearance_mm > 2.4);
        if length == 5.5 {
            assert_eq!(
                stock
                    .rust_wood_insert_screw_clearance(&screw, 1., None)
                    .unwrap_err()
                    .kind,
                GeometryErrorKind::Unsupported
            )
        } else {
            assert!(
                stock
                    .rust_wood_insert_screw_clearance(&screw, 1., None)
                    .is_ok()
            )
        }
        assert_eq!(
            stock
                .rust_wood_insert_screw_clearance(&screw, 0.4, Some(&|| true))
                .unwrap_err()
                .kind,
            GeometryErrorKind::Cancelled
        );
    }
}

#[test]
fn periodic_positive_neck_and_exact_two_term_prefix_keep_real_stock_clearance() {
    for plain in [0., 0.01] {
        let (stock, regular) = fixture(6.);
        let (original, mut frame) = regular.rust_machine_screw_geometry().unwrap();
        frame.origin[2] += 0.5;
        let screw = super::super::mating::tests::periodic_screw(*original.thread(), plain, frame);
        let proof = stock
            .rust_wood_insert_screw_clearance(&screw, 0.4, None)
            .unwrap();
        assert!(proof.tip_clearance_mm[0] <= 0.5 && proof.tip_clearance_mm[1] >= 0.5);
        let expected = 0.4 + plain;
        assert!(
            proof.prefix_entry_clearance_mm[0] > expected - 1e-10
                && proof.prefix_entry_clearance_mm[1] < expected + 1e-10
        );
        assert_eq!(
            stock
                .rust_wood_insert_screw_clearance(&screw, 1., None)
                .unwrap_err()
                .kind,
            GeometryErrorKind::Unsupported
        );
    }
}
