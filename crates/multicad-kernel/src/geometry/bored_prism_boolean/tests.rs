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

fn block(size: [f64; 3], at: [f64; 3]) -> Solid {
    Solid::box_solid(size, policy())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: at,
            ..Default::default()
        })
        .unwrap()
}

/// A 120 x 90 x 18 plate with three 3 mm bores along x = 30, 60, 90.
fn bored_plate() -> Solid {
    [30., 60., 90.]
        .into_iter()
        .fold(block([120., 90., 18.], [0.; 3]), |plate, x| {
            let bore = Solid::cylinder(
                PlaneFrame {
                    origin: [x, 70., -1.],
                    ..Default::default()
                },
                3.,
                20.,
                policy(),
            )
            .unwrap();
            plate
                .boolean(&bore, BooleanOperation::Subtract, policy())
                .unwrap()
                .unwrap()
        })
}

#[test]
fn a_through_slot_clear_of_the_bores_cuts_the_bored_prism_exactly() {
    let plate = bored_plate();
    let (prism, _) = plate.rust_prism_frame().unwrap();
    assert_eq!(prism.inner_loops().len(), 3);
    let slot = block([7., 40., 20.], [8.5, 20., -1.]);
    let cut = plate
        .boolean(&slot, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap();
    near(
        cut.volume_mm3().unwrap(),
        120. * 90. * 18. - 3. * 9. * PI * 18. - 7. * 40. * 18.,
    );
    // Every bore is still a full-height cylindrical wall.
    let walls = cut
        .faces()
        .unwrap()
        .into_iter()
        .filter(|face| face.plane.is_none())
        .count();
    assert_eq!(walls, 3);
    // A blind pocket clear of the bores commutes the same way.
    let pocket = block([10., 10., 5.], [100., 10., 15.]);
    let cut = plate
        .boolean(&pocket, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap();
    near(
        cut.volume_mm3().unwrap(),
        120. * 90. * 18. - 3. * 9. * PI * 18. - 10. * 10. * 3.,
    );
}

#[test]
fn tools_reaching_a_bore_or_other_operations_are_not_commuted() {
    let plate = bored_plate();
    // The slot's bounds overlap the middle bore's bounds.
    let across = block([7., 40., 20.], [58., 50., -1.]);
    // The slot's bounds touch the first bore's bounds within tolerance.
    let touching = block([7., 40., 20.], [33. + 1e-9, 50., -1.]);
    for tool in [&across, &touching] {
        assert!(
            plate
                .bored_prism_subtract(tool, BooleanOperation::Subtract, policy())
                .is_none()
        );
    }
    let clear = block([7., 40., 20.], [8.5, 20., -1.]);
    for operation in [BooleanOperation::Union, BooleanOperation::Intersect] {
        assert!(
            plate
                .bored_prism_subtract(&clear, operation, policy())
                .is_none()
        );
    }
    // A plain box has no bores to commute.
    assert!(
        block([10., 10., 10.], [0.; 3])
            .bored_prism_subtract(&clear, BooleanOperation::Subtract, policy())
            .is_none()
    );
}
