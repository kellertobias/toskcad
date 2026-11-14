use super::*;
fn req(policy: GeometryPolicy) -> PassageRequirements {
    PassageRequirements {
        minimum_radial_clearance_mm: policy.minimum_feature_mm,
        linear_tolerance_mm: policy.linear_tolerance_mm,
        minimum_feature_mm: policy.minimum_feature_mm,
    }
}
fn clamp(policy: GeometryPolicy) -> Solid {
    let stock = Solid::box_solid([20., 20., 4.1], policy)
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: [0., 0., 20.],
            ..Default::default()
        })
        .unwrap();
    stock
        .rust_connection_through_bore(
            &CircularPocketTool {
                frame: PlaneFrame {
                    origin: [10., 10., 24.101],
                    u: [1., 0., 0.],
                    v: [0., -1., 0.],
                },
                center_uv: [0.; 2],
                radius_mm: 1.7,
                displacement_mm: [0., 0., -4.102],
                origin: Origin {
                    source: EntityId::new(),
                    role: "actual.clamp".into(),
                },
            },
            1.,
            policy,
        )
        .unwrap()
        .0
}
fn ring(policy: GeometryPolicy) -> Solid {
    Solid::rust_annulus(
        1.7,
        3.5,
        0.2,
        EntityId::new(),
        PlaneFrame {
            origin: [10., 10., 24.3],
            u: [1., 0., 0.],
            v: [0., -1., 0.],
        },
        policy,
    )
    .unwrap()
}
#[test]
fn actual_clamp_and_washer_sources_clear_hanger_and_reject_real_interface_penetration() {
    let policy = GeometryPolicy::default();
    for posed in [false, true] {
        let (_, mut hanger, _, _) = super::super::tests::fixture(false, false, policy);
        let mut clamp = clamp(policy);
        let mut ring = ring(policy);
        if posed {
            let pose = RigidTransform {
                rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
                translation_mm: [11., -7., 5.],
            };
            hanger = hanger.transformed(pose).unwrap();
            clamp = clamp.transformed(pose).unwrap();
            ring = ring.transformed(pose).unwrap();
        }
        let proof = clamp
            .rust_hanger_through_clearance(&hanger, 1., req(policy), None)
            .unwrap();
        assert!(proof.radial_clearance_lower_mm > 0.19);
        assert!(
            ring.rust_hanger_annulus_clearance(&hanger, req(policy), None)
                .unwrap()
                .radial_clearance_lower_mm
                > 0.19
        );
    }
    let (_, hanger, _, _) = super::super::tests::fixture(false, false, policy);
    let moved = hanger
        .transformed(RigidTransform {
            translation_mm: [0., 0., 1e-10],
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        clamp(policy)
            .rust_hanger_through_clearance(&moved, 1., req(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
}
#[test]
fn material_ownership_policies_and_cancellation_are_not_replaced_by_nominal_holes() {
    let policy = GeometryPolicy::default();
    let (pilot, hanger, _, _) = super::super::tests::fixture(false, false, policy);
    let clamp = clamp(policy);
    let ring = ring(policy);
    let unknown = Solid::box_solid([2.; 3], policy).unwrap();
    let uncut = Solid::box_solid([20.; 3], policy)
        .unwrap()
        .orthogonal_stock()
        .unwrap();
    assert_eq!(
        uncut
            .rust_hanger_through_clearance(&hanger, 1., req(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        unknown
            .rust_hanger_through_clearance(&hanger, 1., req(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        pilot
            .rust_hanger_through_clearance(&hanger, 1., req(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        clamp
            .rust_hanger_through_clearance(&unknown, 1., req(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        unknown
            .rust_hanger_annulus_clearance(&hanger, req(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        ring.rust_hanger_annulus_clearance(&unknown, req(policy), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        clamp
            .rust_hanger_through_clearance(&hanger, 1., req(policy), Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert_eq!(
        ring.rust_hanger_annulus_clearance(&hanger, req(policy), Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let mut bad = req(policy);
    bad.linear_tolerance_mm = 1e-6;
    assert_eq!(
        clamp
            .rust_hanger_through_clearance(&hanger, 1., bad, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        ring.rust_hanger_annulus_clearance(&hanger, bad, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
}
