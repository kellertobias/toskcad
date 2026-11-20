use super::*;
use std::f64::consts::PI;

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-7 * expected.abs().max(1.),
        "{actual} != {expected}"
    );
}

fn rod(x: f64, z: f64, radius: f64, height: f64) -> Solid {
    Solid::cylinder(
        PlaneFrame {
            origin: [x, 28., z],
            ..Default::default()
        },
        radius,
        height,
        policy(),
    )
    .unwrap()
}

fn cut(a: &Solid, b: &Solid) -> Solid {
    a.boolean(b, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap()
}

/// The drill-guide fence: a 120 x 14 x 22 bar below z = 0.
fn fence() -> Solid {
    Solid::box_solid([120., 14., 22.], policy())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: [0., 21., -22.],
            ..Default::default()
        })
        .unwrap()
}

/// The horizontal annular shoulder faces of a solid, as (z, area).
fn shoulders(solid: &Solid, outer: f64, inner: f64) -> Vec<f64> {
    let area = PI * (outer * outer - inner * inner);
    solid
        .faces()
        .unwrap()
        .into_iter()
        .filter(|face| face.plane.is_some() && (face.area_mm2 - area).abs() < 1e-6)
        .map(|face| face.center_mm[2])
        .collect()
}

#[test]
fn insert_pockets_around_clearance_bores_are_exact_counterbores() {
    let mut body = fence();
    for x in [12., 108.] {
        body = cut(&body, &rod(x, -23., 3.3, 24.));
        body = cut(&body, &rod(x, -22.1, 4.1, 6.));
    }
    let bores = 2. * PI * 3.3 * 3.3 * 22.;
    let seats = 2. * PI * (4.1 * 4.1 - 3.3 * 3.3) * 5.9;
    near(body.volume_mm3().unwrap(), 120. * 14. * 22. - bores - seats);
    let placed = body.polygon_prism().unwrap();
    let treatments = placed.bore_treatments.as_ref().unwrap();
    assert_eq!(treatments.geometry.counterbores().len(), 2);
    // Each shoulder sits 5.9 mm above the bottom cap.
    let levels = shoulders(&body, 4.1, 3.3);
    assert_eq!(levels.len(), 2);
    for z in levels {
        near(z, -22. + 5.9);
    }
}

#[test]
fn a_seat_from_the_top_cap_is_measured_from_that_cap() {
    let body = cut(&fence(), &rod(60., -23., 3.3, 24.));
    let seated = cut(&body, &rod(60., -4., 5., 5.));
    near(
        seated.volume_mm3().unwrap(),
        120. * 14. * 22. - PI * 3.3 * 3.3 * 22. - PI * (25. - 3.3 * 3.3) * 4.,
    );
    let levels = shoulders(&seated, 5., 3.3);
    assert_eq!(levels.len(), 1);
    near(levels[0], -4.);
}

#[test]
fn off_axis_narrow_through_or_repeated_seats_are_not_counterbores() {
    let body = cut(&fence(), &rod(60., -23., 3.3, 24.));
    let declined = |tool: &Solid| {
        body.coaxial_counterbore_cut(tool, BooleanOperation::Subtract, policy())
            .is_none()
    };
    // Off the bore axis, no wider than the bore, through both caps, or
    // starting inside the stock.
    assert!(declined(&rod(60.5, -22.1, 4.1, 6.)));
    assert!(declined(&rod(60., -22.1, 3.3, 6.)));
    assert!(declined(&rod(60., -23., 4.1, 24.)));
    assert!(declined(&rod(60., -20., 4.1, 6.)));
    // A second seat on the same end of the same bore.
    let seated = cut(&body, &rod(60., -22.1, 4.1, 6.));
    assert!(
        seated
            .coaxial_counterbore_cut(
                &rod(60., -22.1, 5., 3.),
                BooleanOperation::Subtract,
                policy()
            )
            .is_none()
    );
    // A plain box with no bores has nothing to seat or reapply.
    assert!(
        fence()
            .coaxial_counterbore_cut(
                &rod(60., -23., 3.3, 24.),
                BooleanOperation::Subtract,
                policy()
            )
            .is_none()
    );
    assert!(
        body.coaxial_counterbore_cut(&rod(60., -22.1, 4.1, 6.), BooleanOperation::Union, policy())
            .is_none()
    );
}
