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

fn pin() -> Solid {
    Solid::cylinder(
        PlaneFrame {
            origin: [12., 28., -21.],
            ..Default::default()
        },
        3.,
        50.,
        policy(),
    )
    .unwrap()
}

#[test]
fn an_inscribing_cap_block_flush_with_the_pin_end_joins_into_one_stack() {
    // The 6 mm square cap touches the 3 mm pin on all four sides.
    let cap = block([6., 6., 3.], [9., 25., 26.]);
    let joined = pin()
        .boolean(&cap, BooleanOperation::Union, policy())
        .unwrap()
        .unwrap();
    near(joined.volume_mm3().unwrap(), 9. * PI * 47. + 36. * 3.);
    assert!(joined.axial_stack_body().is_some());
    let bounds = joined.bounds().unwrap();
    near(bounds.min[2], -21.);
    near(bounds.max[2], 29.);
    near(bounds.min[0], 9.);
    // The union is symmetric and repeatable in its identities.
    let swapped = cap
        .boolean(&pin(), BooleanOperation::Union, policy())
        .unwrap()
        .unwrap();
    near(swapped.volume_mm3().unwrap(), joined.volume_mm3().unwrap());
}

#[test]
fn blocks_that_do_not_inscribe_the_pin_or_other_operations_are_not_stacked() {
    // The block is narrower than the pin's circle.
    let narrow = block([4., 4., 3.], [10., 26., 26.]);
    assert!(
        pin()
            .stacked_cylinder_union(&narrow, BooleanOperation::Union, policy())
            .is_none()
    );
    // A cap with clearance on every side would need an annular cap face,
    // which the profile stack does not certify.
    let clear = block([8., 8., 3.], [8., 24., 26.]);
    assert!(
        pin()
            .stacked_cylinder_union(&clear, BooleanOperation::Union, policy())
            .is_none()
    );
    let cap = block([6., 6., 3.], [9., 25., 26.]);
    for operation in [BooleanOperation::Subtract, BooleanOperation::Intersect] {
        assert!(
            pin()
                .stacked_cylinder_union(&cap, operation, policy())
                .is_none()
        );
    }
    // Two blocks have no cylinder to stack on.
    assert!(
        cap.stacked_cylinder_union(&narrow, BooleanOperation::Union, policy())
            .is_none()
    );
}
