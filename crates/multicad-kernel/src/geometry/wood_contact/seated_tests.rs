use super::*;
use multicad_rust_kernel::cylindrical_blind_pocket::Seat;
fn recipient() -> Solid {
    let policy = GeometryPolicy::default();
    let stock =
        crate::geometry::construction::box_solid([20.; 3], EntityId::new(), policy).unwrap();
    stock
        .rust_cut_rectangular_blind_counterbore(
            CylindricalBlindPocketSpec {
                radius_mm: 2.,
                depth_mm: 16.,
                minimum_wall_mm: 1.,
                minimum_floor_mm: 1.,
                seat: Some(Seat {
                    radius_mm: 6.,
                    depth_mm: 1.,
                }),
            },
            PlaneFrame {
                origin: [10., 10., 20.],
                u: [1., 0., 0.],
                v: [0., -1., 0.],
            },
            EntityId::new(),
            policy,
            None,
        )
        .unwrap()
}
#[test]
fn actual_seated_sources_and_occurrence_poses_keep_full_thread_contact() {
    let cut = recipient();
    let (_, screw, req) = super::tests::fixture();
    let screw = screw
        .transformed(RigidTransform {
            translation_mm: [0., 0., -1.],
            ..Default::default()
        })
        .unwrap();
    let proof = cut
        .rust_seated_direct_wood_contact(&screw, req, None)
        .unwrap();
    assert_eq!(proof.seat_radius_mm, 6.);
    assert_eq!(proof.seat_depth_mm, 1.);
    assert!(proof.contact.engagement_mm[0] <= 12. && proof.contact.engagement_mm[1] >= 12.);
    let pose = RigidTransform {
        rotation: [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]],
        translation_mm: [9., -3., 7.],
    };
    let placed = cut.transformed(pose).unwrap();
    let screw = screw.transformed(pose).unwrap();
    assert!(
        placed
            .rust_seated_direct_wood_contact(&screw, req, None)
            .is_ok()
    );
    assert_eq!(
        cut.rust_seated_direct_wood_contact(&screw, req, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
}
#[test]
fn source_policy_identity_and_cancel_are_required_before_contact_publication() {
    let cut = recipient();
    let (pilot, screw, req) = super::tests::fixture();
    assert_eq!(
        pilot
            .rust_seated_direct_wood_contact(&screw, req, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        cut.rust_seated_direct_wood_contact(&cut, req, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        cut.rust_seated_direct_wood_contact(&screw, req, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let screw = screw
        .transformed(RigidTransform {
            translation_mm: [0., 0., -1.],
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        cut.rust_seated_direct_wood_contact(&screw, req, Some(&|| true))
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let mut wrong = req;
    wrong.linear_tolerance_mm = 1e-6;
    assert_eq!(
        cut.rust_seated_direct_wood_contact(&screw, wrong, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
}
