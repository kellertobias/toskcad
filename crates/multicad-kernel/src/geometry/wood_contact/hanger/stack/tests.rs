use super::*;
fn req(policy: GeometryPolicy) -> StackRequirements {
    StackRequirements {
        wood: Requirements {
            minimum_wall_mm: 1.,
            minimum_floor_mm: 1.,
            minimum_engagement_mm: 4.,
            tip_clearance_mm: 1.,
            linear_tolerance_mm: policy.linear_tolerance_mm,
            minimum_feature_mm: policy.minimum_feature_mm,
        },
        minimum_clamp_wall_mm: 1.,
        passage: source::PassageRequirements {
            minimum_radial_clearance_mm: policy.minimum_feature_mm,
            linear_tolerance_mm: policy.linear_tolerance_mm,
            minimum_feature_mm: policy.minimum_feature_mm,
        },
        mating: HangerMatingRequirements {
            minimum_engagement_mm: 2.,
            ..Default::default()
        },
        minimum_bearing_land_mm: policy.minimum_feature_mm,
        maximum_bearing_gap_mm: policy.linear_tolerance_mm,
    }
}
fn clamp(policy: GeometryPolicy) -> Solid {
    Solid::box_solid([20., 20., 4.1], policy)
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: [0., 0., 20.],
            ..Default::default()
        })
        .unwrap()
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
fn actual_current_native_stack_proves_all_related_materials() {
    let p = GeometryPolicy::default();
    for left in [false, true] {
        let (wood, hanger, nut, _) = super::super::tests::fixture(left, false, p);
        let clamp = clamp(p);
        let ring = ring(p);
        // Independently rounded placements really penetrate at the nut/ring
        // cap. This remains a refusal, even though the gap is below tolerance.
        assert_eq!(
            hanger
                .rust_hanger_stack(&wood, &clamp, &nut, &[&ring], req(p), None)
                .unwrap_err()
                .kind,
            GeometryErrorKind::Unsupported
        );
        // This test physically places the nut clear; the compiler must not
        // silently make this correction to arbitrary current occurrences.
        let nut = nut
            .transformed(RigidTransform {
                translation_mm: [0., 0., 1e-11],
                ..Default::default()
            })
            .unwrap();
        let result = hanger
            .rust_hanger_stack(&wood, &clamp, &nut, &[&ring], req(p), None)
            .unwrap();
        assert_eq!(result.bearings.len(), 3);
        assert!(
            result
                .bearings
                .iter()
                .all(|b| b.radial_land_lower_mm >= p.minimum_feature_mm)
        );
        let into_clamp = ring
            .transformed(RigidTransform {
                translation_mm: [0., 0., -1e-10],
                ..Default::default()
            })
            .unwrap();
        assert_eq!(
            hanger
                .rust_hanger_stack(&wood, &clamp, &nut, &[&into_clamp], req(p), None)
                .unwrap_err()
                .kind,
            GeometryErrorKind::Unsupported
        );
    }
}

#[test]
fn wrong_current_sources_policies_capacity_and_cancellation_never_publish_a_stack() {
    let p = GeometryPolicy::default();
    let (wood, hanger, nut, _) = super::super::tests::fixture(false, false, p);
    let clamp = clamp(p);
    let ring = ring(p);
    let unknown = Solid::box_solid([2.; 3], p).unwrap();
    let uncut = Solid::box_solid([20.; 3], p)
        .unwrap()
        .orthogonal_stock()
        .unwrap();
    for case in 0..7 {
        let result = match case {
            0 => unknown.rust_hanger_stack(&wood, &clamp, &nut, &[&ring], req(p), None),
            1 => hanger.rust_hanger_stack(&unknown, &clamp, &nut, &[&ring], req(p), None),
            2 => hanger.rust_hanger_stack(&uncut, &clamp, &nut, &[&ring], req(p), None),
            3 => hanger.rust_hanger_stack(&wood, &unknown, &nut, &[&ring], req(p), None),
            4 => hanger.rust_hanger_stack(&wood, &wood, &nut, &[&ring], req(p), None),
            5 => hanger.rust_hanger_stack(&wood, &clamp, &unknown, &[&ring], req(p), None),
            _ => hanger.rust_hanger_stack(&wood, &clamp, &nut, &[&unknown], req(p), None),
        };
        assert_eq!(
            result.unwrap_err().kind,
            GeometryErrorKind::Unsupported,
            "case{case}"
        );
    }
    for case in 0..5 {
        let mut requirements = req(p);
        match case {
            0 => requirements.wood.linear_tolerance_mm = 1e-6,
            1 => requirements.wood.minimum_feature_mm = 1e-4,
            2 => requirements.passage.linear_tolerance_mm = 1e-6,
            3 => requirements.passage.minimum_feature_mm = 1e-4,
            _ => requirements.mating.linear_tolerance_mm = 1e-6,
        }
        assert_eq!(
            hanger
                .rust_hanger_stack(&wood, &clamp, &nut, &[&ring], requirements, None)
                .unwrap_err()
                .kind,
            GeometryErrorKind::Input
        );
    }
    assert_eq!(
        hanger
            .rust_hanger_stack(&wood, &clamp, &nut, &vec![&ring; 63], req(p), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Limit
    );
    assert_eq!(
        hanger
            .rust_hanger_stack(&wood, &clamp, &nut, &[&ring], req(p), Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let other = GeometryPolicy {
        linear_tolerance_mm: 1e-6,
        minimum_feature_mm: 1e-4,
        ..p
    };
    let different = ring.transformed(RigidTransform::default()).unwrap();
    let other_solid = Solid::box_solid([2.; 3], other).unwrap();
    assert_eq!(
        hanger
            .rust_hanger_stack(&wood, &clamp, &nut, &[&other_solid], req(p), None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(different.rust_annulus_data(), ring.rust_annulus_data());
}
