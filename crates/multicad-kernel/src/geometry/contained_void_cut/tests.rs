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

fn stock() -> Solid {
    Solid::box_solid([20., 20., 4.], policy()).unwrap()
}

/// The BOM washer: bore 4.3, outside 9, 0.8 thick, centred at (10, 10, 1).
fn washer(z: f64) -> Solid {
    Solid::rust_annulus(
        2.15,
        4.5,
        0.8,
        EntityId::new(),
        PlaneFrame {
            origin: [10., 10., z],
            ..Default::default()
        },
        policy(),
    )
    .unwrap()
}

#[test]
fn a_buried_washer_becomes_one_reversed_void_shell() {
    let cut = stock()
        .boolean(&washer(1.), BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap();
    let ring = PI * (4.5 * 4.5 - 2.15 * 2.15) * 0.8;
    near(cut.volume_mm3().unwrap(), 1600. - ring);
    let body = cut.imported_brep_body().unwrap();
    let shells = &body.solid.boundary().shells;
    assert_eq!(shells.len(), 2);
    assert!(shells[1].reversed);
    let bounds = cut.bounds().unwrap();
    near(bounds.max[0], 20.);
    near(bounds.max[2], 4.);
    // The void's faces join the stock's: 6 box faces plus the washer's 4.
    assert_eq!(cut.faces().unwrap().len(), 10);
}

#[test]
fn tools_reaching_the_stock_boundary_are_not_voids() {
    let declined = |stock: &Solid, tool: &Solid| {
        stock
            .contained_void_cut(tool, BooleanOperation::Subtract, policy())
            .is_none()
    };
    // Touching the floor, breaking through the top, or a cut that is not a
    // subtraction.
    assert!(declined(&stock(), &washer(0.)));
    assert!(declined(&stock(), &washer(3.5)));
    assert!(
        stock()
            .contained_void_cut(&washer(1.), BooleanOperation::Intersect, policy())
            .is_none()
    );
    // Holed stock is not convex.
    let holed = stock()
        .boolean(
            &Solid::cylinder(
                PlaneFrame {
                    origin: [3., 3., -1.],
                    ..Default::default()
                },
                1.,
                6.,
                policy(),
            )
            .unwrap(),
            BooleanOperation::Subtract,
            policy(),
        )
        .unwrap()
        .unwrap();
    assert!(declined(&holed, &washer(1.)));
}
