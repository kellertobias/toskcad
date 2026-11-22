//! Intersections and cuts proved empty when an axial body passes through a
//! hole with clearance, and refused when it does not fit.
use multicad_core::document::EntityId;
use multicad_kernel::{BooleanOperation, GeometryPolicy, PlaneFrame, Solid};

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

/// A 20 × 20 × 5 plate with a through bore of `radius` at (10, 10).
fn bored_plate(radius: f64) -> Solid {
    let plate = Solid::box_solid([20., 20., 5.], policy()).unwrap();
    let bore = Solid::cylinder(
        PlaneFrame {
            origin: [10., 10., -1.],
            ..Default::default()
        },
        radius,
        7.,
        policy(),
    )
    .unwrap();
    plate
        .boolean(&bore, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap()
}

/// A hex nut of `across_flats` standing on the plate's floor plane.
fn nut(across_flats: f64) -> Solid {
    Solid::rust_cosmetic_hex_nut(
        across_flats,
        0.5,
        5.,
        0.1,
        EntityId::new(),
        PlaneFrame {
            origin: [10., 10., 0.],
            ..Default::default()
        },
        policy(),
    )
    .unwrap()
}

#[test]
fn a_hex_body_inside_a_through_bore_shares_no_material_with_the_plate() {
    let plate = bored_plate(3.);
    // Corner radius 4 / √3 ≈ 2.31 mm stays inside the 3 mm bore.
    let inside = nut(4.);
    for (a, b) in [(&plate, &inside), (&inside, &plate)] {
        assert!(
            a.boolean(b, BooleanOperation::Intersect, policy())
                .unwrap()
                .is_none()
        );
    }
    let kept = plate
        .boolean(&inside, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap();
    assert_eq!(kept.volume_mm3().unwrap(), plate.volume_mm3().unwrap());
    assert_eq!(kept.bounds().unwrap(), plate.bounds().unwrap());
}

#[test]
fn a_hex_body_wider_than_the_bore_is_never_reported_clear() {
    let plate = bored_plate(3.);
    // Corner radius 6 / √3 ≈ 3.46 mm reaches into the plate.
    let wide = nut(6.);
    for (a, b) in [(&plate, &wide), (&wide, &plate)] {
        assert!(!matches!(
            a.boolean(b, BooleanOperation::Intersect, policy()),
            Ok(None)
        ));
    }
}

/// A 20 × 20 × 5 box with a blind pocket of radius 3 and depth 2 from the top.
fn pocketed_block() -> Solid {
    let block = Solid::box_solid([20., 20., 5.], policy()).unwrap();
    let pocket = Solid::cylinder(
        PlaneFrame {
            origin: [10., 10., 3.],
            ..Default::default()
        },
        3.,
        3.,
        policy(),
    )
    .unwrap();
    block
        .boolean(&pocket, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap()
}

fn pin(radius: f64, bottom: f64) -> Solid {
    Solid::cylinder(
        PlaneFrame {
            origin: [10., 10., bottom],
            ..Default::default()
        },
        radius,
        10.,
        policy(),
    )
    .unwrap()
}

#[test]
fn a_pin_in_a_blind_pocket_and_above_the_block_misses_its_material() {
    let block = pocketed_block();
    // Down to the floor at z = 3 and up beyond the top face.
    let resting = pin(2.5, 3.);
    assert!(
        block
            .boolean(&resting, BooleanOperation::Intersect, policy())
            .unwrap()
            .is_none()
    );
    assert!(
        resting
            .boolean(&block, BooleanOperation::Intersect, policy())
            .unwrap()
            .is_none()
    );
    // Below the floor, or wider than the pocket, the pin reaches material.
    for pin in [pin(2.5, 2.), pin(3.5, 3.)] {
        assert!(!matches!(
            block.boolean(&pin, BooleanOperation::Intersect, policy()),
            Ok(None)
        ));
    }
}
