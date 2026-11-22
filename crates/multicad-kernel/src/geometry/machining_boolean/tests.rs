use super::*;
use crate::{CurvedEdge2d, CurvedProfile};
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

fn rectangle(w: f64, h: f64) -> Vec<CurvedEdge2d> {
    let corners = [[0., 0.], [w, 0.], [w, h], [0., h]];
    (0..4)
        .map(|i| CurvedEdge2d::Line {
            a: corners[i],
            b: corners[(i + 1) % 4],
        })
        .collect()
}

fn circle(frame: PlaneFrame, center: [f64; 2], radius_mm: f64, travel: [f64; 3]) -> Solid {
    Solid::extrude_curved_profile_vector(
        &CurvedProfile {
            frame,
            outer: vec![CurvedEdge2d::Circle { center, radius_mm }],
            holes: vec![],
        },
        travel,
        policy(),
    )
    .unwrap()
}

fn subtract(a: &Solid, b: &Solid) -> Result<Option<Solid>, GeometryError> {
    a.boolean(b, BooleanOperation::Subtract, policy())
}

/// A plain sweep along a tilted frame, with a circular hole.
fn tilted_stock() -> (PlaneFrame, Solid) {
    let frame = PlaneFrame {
        origin: [10., 20., 30.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    };
    let stock = Solid::extrude_curved_profile_vector(
        &CurvedProfile {
            frame,
            outer: rectangle(30., 20.),
            holes: vec![vec![CurvedEdge2d::Circle {
                center: [5., 5.],
                radius_mm: 1.,
            }]],
        },
        [-12., 0., 0.],
        policy(),
    )
    .unwrap();
    assert!(stock.normal_sweep_body().is_some());
    (frame, stock)
}

#[test]
fn blind_cylinder_cuts_a_straight_sweep_as_an_exact_pocket() {
    let (frame, stock) = tilted_stock();
    let tool = circle(
        PlaneFrame {
            origin: [11., 20., 30.],
            ..frame
        },
        [12., 10.],
        2.,
        [-5., 0., 0.],
    );
    let cut = subtract(&stock, &tool).unwrap().unwrap();
    // The pocket enters 1 mm outside the cap and leaves a 4 mm deep floor.
    near(cut.volume_mm3().unwrap(), 7200. - 12. * PI - 16. * PI);
    assert!(cut.polygon_prism().is_some_and(|p| p.pockets.len() == 1));
    // A second pocket from the opposite cap cuts the first Boolean result.
    let opposite = circle(
        PlaneFrame {
            origin: [-3., 20., 30.],
            ..frame
        },
        [12., 10.],
        2.,
        [4., 0., 0.],
    );
    let cut = subtract(&cut, &opposite).unwrap().unwrap();
    near(
        cut.volume_mm3().unwrap(),
        7200. - 12. * PI - 16. * PI - 12. * PI,
    );
    assert!(cut.polygon_prism().is_some_and(|p| p.pockets.len() == 2));
    // Two caps, four walls and the hole, plus a wall and floor per pocket.
    assert_eq!(cut.faces().unwrap().len(), 11);
    let reading = stock.normal_sweep_polygon_stock().unwrap();
    near(reading.volume_mm3().unwrap(), stock.volume_mm3().unwrap());
}

#[test]
fn through_cylinder_on_placed_polygon_stock_becomes_a_polygon_hole() {
    let (frame, sweep) = tilted_stock();
    let stock = sweep.normal_sweep_polygon_stock().unwrap();
    assert!(stock.polygon_prism().is_some_and(|p| p.frame == frame));
    let tool = circle(
        PlaneFrame {
            origin: [11., 20., 30.],
            ..frame
        },
        [20., 10.],
        3.,
        [-14., 0., 0.],
    );
    let cut = subtract(&stock, &tool).unwrap().unwrap();
    near(cut.volume_mm3().unwrap(), 7200. - 12. * PI - 108. * PI);
    let placed = cut.polygon_prism().unwrap();
    assert_eq!(placed.prism.circular_holes().len(), 2);
    assert!(placed.pockets.is_empty());
}

#[test]
fn stadium_sweep_cuts_box_stock_as_a_rounded_slot() {
    let stock = Solid::box_solid([20., 16., 12.], policy()).unwrap();
    let frame = PlaneFrame {
        origin: [5., 6., 13.],
        ..Default::default()
    };
    let [a, b, c, d] = [[0., -2.], [4., -2.], [4., 2.], [0., 2.]];
    let slot = Solid::extrude_curved_profile_vector(
        &CurvedProfile {
            frame,
            outer: vec![
                CurvedEdge2d::Line { a, b },
                CurvedEdge2d::Arc {
                    center: [4., 0.],
                    a: b,
                    b: c,
                    clockwise: false,
                },
                CurvedEdge2d::Line { a: c, b: d },
                CurvedEdge2d::Arc {
                    center: [0., 0.],
                    a: d,
                    b: a,
                    clockwise: false,
                },
            ],
            holes: vec![],
        },
        [0., 0., -5.],
        policy(),
    )
    .unwrap();
    let tool = slot.normal_sweep_slot_tool().unwrap();
    assert_eq!(tool.centers_uv, [[4., 0.], [0., 0.]]);
    assert_eq!(tool.radius_mm, 2.);
    assert_eq!(tool.displacement_mm, [0., 0., -5.]);
    let cut = subtract(&stock, &slot).unwrap().unwrap();
    let area = 16. + 4. * PI;
    near(cut.volume_mm3().unwrap(), 20. * 16. * 12. - 4. * area);
    // A blind round hole beside the slot is then an exact pocket too.
    let hole = circle(frame, [9., 4.], 1., [0., 0., -4.]);
    let cut = subtract(&cut, &hole).unwrap().unwrap();
    near(
        cut.volume_mm3().unwrap(),
        20. * 16. * 12. - 4. * area - 3. * PI,
    );
}

#[test]
fn curved_rings_and_non_stadium_arcs_have_no_machining_reading() {
    let frame = PlaneFrame::default();
    // A quarter-round outline is neither straight stock nor a stadium.
    let fillet = Solid::extrude_curved_profile_vector(
        &CurvedProfile {
            frame,
            outer: vec![
                CurvedEdge2d::Line {
                    a: [0., 0.],
                    b: [2., 0.],
                },
                CurvedEdge2d::Arc {
                    center: [0., 0.],
                    a: [2., 0.],
                    b: [0., 2.],
                    clockwise: false,
                },
                CurvedEdge2d::Line {
                    a: [0., 2.],
                    b: [0., 0.],
                },
            ],
            holes: vec![],
        },
        [0., 0., 3.],
        policy(),
    )
    .unwrap();
    assert!(fillet.normal_sweep_polygon_stock().is_none());
    assert!(fillet.normal_sweep_slot_tool().is_none());
    // Unequal end radii are not a slot.
    let pill = Solid::extrude_curved_profile_vector(
        &CurvedProfile {
            frame,
            outer: vec![
                CurvedEdge2d::Line {
                    a: [0., -2.],
                    b: [4., -1.],
                },
                CurvedEdge2d::Arc {
                    center: [4., 0.],
                    a: [4., -1.],
                    b: [4., 1.],
                    clockwise: false,
                },
                CurvedEdge2d::Line {
                    a: [4., 1.],
                    b: [0., 2.],
                },
                CurvedEdge2d::Arc {
                    center: [0., 0.],
                    a: [0., 2.],
                    b: [0., -2.],
                    clockwise: false,
                },
            ],
            holes: vec![],
        },
        [0., 0., 3.],
        policy(),
    )
    .unwrap();
    assert!(pill.normal_sweep_slot_tool().is_none());
}

#[test]
fn refused_machining_keeps_the_typed_boolean_refusal() {
    let (frame, stock) = tilted_stock();
    // A tool that stops inside the stock without reaching a cap is neither a
    // blind pocket nor a through hole.
    let buried = circle(
        PlaneFrame {
            origin: [5., 20., 30.],
            ..frame
        },
        [20., 10.],
        2.,
        [-3., 0., 0.],
    );
    assert!(
        stock
            .subtract_machining_tool(&buried, BooleanOperation::Subtract, policy())
            .is_none()
    );
    let Err(error) = subtract(&stock, &buried) else {
        panic!("a buried tool has no exact machining result");
    };
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    // Union and intersection are never read as machining.
    let tool = circle(
        PlaneFrame {
            origin: [11., 20., 30.],
            ..frame
        },
        [12., 10.],
        2.,
        [-5., 0., 0.],
    );
    for operation in [BooleanOperation::Union, BooleanOperation::Intersect] {
        assert!(
            stock
                .subtract_machining_tool(&tool, operation, policy())
                .is_none()
        );
    }
}
