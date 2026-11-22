use super::*;
#[test]
fn actual_counterbored_clamp_retains_original_stock_height_and_source_frame() {
    let policy = GeometryPolicy::default();
    let stock = Solid::box_solid([20., 20., 4.1], policy).unwrap();
    let tool = crate::CircularPocketTool {
        frame: PlaneFrame {
            origin: [10., 10., 4.10001],
            u: [1., 0., 0.],
            v: [0., -1., 0.],
        },
        center_uv: [0.; 2],
        radius_mm: 1.7,
        displacement_mm: [0., 0., -4.10002],
        origin: Origin {
            source: EntityId::new(),
            role: "actual.clamp".into(),
        },
    };
    let (clamp, _) = stock
        .rust_connection_counterbore(
            &tool,
            3.5,
            1.,
            1.,
            Origin {
                source: EntityId::new(),
                role: "actual.seat".into(),
            },
            policy,
        )
        .unwrap();
    let placed = clamp.polygon_prism().unwrap();
    let cuts = &placed.bore_treatments.as_ref().unwrap().geometry;
    eprintln!(
        "ACTUAL CLAMP frame={:?}, reverse={}, height={:?}, layers={:?}",
        placed.frame,
        placed.reverse_sweep,
        placed.prism.height().ratio(),
        cuts.stacks()
    );
    assert_eq!(placed.prism.height().ratio(), (41, 10));
    assert_eq!(placed.frame.origin[2], 4.1);
    assert_eq!(placed.frame.normal(), [0., 0., -1.]);
}

#[test]
fn actual_clocked_native_queries_keep_source_identity_policies_and_cancellation() {
    let (original, insert) =
        super::super::mating::tests::fixture(false, 8., 0., PlaneFrame::default());
    let (sin, cos) = (-72_f64.to_radians()).portable_sin_cos();
    let screw = original
        .transformed(RigidTransform {
            rotation: [[cos, -sin, 0.], [sin, cos, 0.], [0., 0., 1.]],
            translation_mm: [0., 0., -3.1],
        })
        .unwrap();
    let proof = screw
        .rust_clocked_wood_insert_mating(&insert, 4., None)
        .unwrap();
    assert!(proof.mating.minimum_radial_clearance_lower_mm > 0.007);
    assert_eq!(
        screw
            .rust_clocked_wood_insert_mating(&insert, 4., Some(&|| true))
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let (stock, _) = super::super::clearance::tests::fixture(4.);
    let screw = screw
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
            translation_mm: [10., 10., 20.],
        })
        .unwrap();
    assert!(
        stock
            .rust_clocked_wood_insert_screw_clearance(&screw, 1., None)
            .unwrap()
            .tip_clearance_mm[0]
            > 1.09
    );
    assert_eq!(
        stock
            .rust_clocked_wood_insert_screw_clearance(&screw, 1., Some(&|| true))
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let box_ = Solid::box_solid([20.; 3], GeometryPolicy::default()).unwrap();
    assert_eq!(
        box_.rust_clocked_wood_insert_mating(&insert, 1., None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        original
            .rust_clocked_wood_insert_mating(&box_, 1., None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        box_.rust_clocked_wood_insert_screw_clearance(&screw, 1., None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        stock
            .rust_clocked_wood_insert_screw_clearance(&box_, 1., None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let other = Solid::box_solid(
        [20.; 3],
        GeometryPolicy {
            linear_tolerance_mm: 1e-6,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        original
            .rust_clocked_wood_insert_mating(&other, 1., None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        other
            .rust_clocked_wood_insert_screw_clearance(&screw, 1., None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
}
