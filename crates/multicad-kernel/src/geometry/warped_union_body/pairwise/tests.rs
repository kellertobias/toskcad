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

fn line(a: [f64; 2], b: [f64; 2]) -> CurvedEdge2d {
    CurvedEdge2d::Line { a, b }
}

/// A bare 2 mm quarter-round fillet along +Y, 10 mm long, without a base
/// or back wall: the material below the arc about (0, 2) in XZ.
fn bare_fillet() -> Solid {
    let profile = CurvedProfile {
        frame: PlaneFrame {
            v: [0., 0., 1.],
            ..Default::default()
        },
        outer: vec![
            line([0., 0.], [2., 0.]),
            line([2., 0.], [2., 2.]),
            CurvedEdge2d::Arc {
                a: [2., 2.],
                b: [0., 0.],
                center: [0., 2.],
                clockwise: true,
            },
        ],
        holes: vec![],
    };
    Solid::extrude_curved_profile_vector(&profile, [0., 10., 0.], policy()).unwrap()
}

/// A vertical prism through the whole band over an XY polygon.
fn vertical_prism(points: &[[f64; 2]]) -> Solid {
    let outer = points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .map(|(&a, &b)| line(a, b))
        .collect();
    Solid::extrude_curved_profile_vector(
        &CurvedProfile {
            frame: PlaneFrame {
                origin: [0., 0., -1.],
                ..Default::default()
            },
            outer,
            holes: vec![],
        },
        [0., 0., 4.],
        policy(),
    )
    .unwrap()
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

fn boss(y: f64, radius: f64) -> Solid {
    Solid::cylinder(
        PlaneFrame {
            origin: [9.95, y, 20.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        },
        radius,
        3.55,
        policy(),
    )
    .unwrap()
}

fn apply(a: &Solid, b: &Solid, operation: BooleanOperation) -> Option<Solid> {
    a.boolean(b, operation, policy()).unwrap()
}

/// The fillet clipped by the plane y = 3 + x / 2: integrate
/// (2 - sqrt(4 - x^2)) (7 - x / 2) over 0 <= x <= 2.
const CLIPPED: f64 = 28. - 7. * PI - 2. / 3.;

#[test]
fn an_oblique_vertical_clip_keeps_the_exact_quarter_round() {
    let fillet = bare_fillet();
    near(fillet.volume_mm3().unwrap(), 10. * (4. - PI));
    let shear = vertical_prism(&[[-1., -1.], [3., -1.], [3., 4.5], [-1., 2.5]]);
    let clipped = apply(&fillet, &shear, BooleanOperation::Subtract).unwrap();
    assert!(clipped.warped_union_body().is_some());
    near(clipped.volume_mm3().unwrap(), CLIPPED);
    let bounds = clipped.bounds().unwrap();
    near(bounds.min[1], 3.);
    near(bounds.max[2], 2.);

    // A block beside the trimmed fillet, partly inside the removed region,
    // joins the recovered cut material rather than the uncut stock.
    let joined = apply(
        &clipped,
        &block([1., 4., 2.], [2., 4., 0.]),
        BooleanOperation::Union,
    )
    .unwrap();
    near(joined.volume_mm3().unwrap(), CLIPPED + 8.);
    assert_eq!(joined.solid_count().unwrap(), 1);
    // The cutter is not part of the material: the block's corner below the
    // clip plane stays.
    near(joined.bounds().unwrap().min[1], 3.);
}

#[test]
fn a_cutter_enclosing_the_fillet_removes_everything() {
    let all = vertical_prism(&[[-1., -1.], [3., -1.], [3., 11.], [-1., 11.]]);
    assert!(apply(&bare_fillet(), &all, BooleanOperation::Subtract).is_none());
}

#[test]
fn bosses_join_planar_stock_one_at_a_time_in_any_placement() {
    let stock = block([10., 80., 27.], [0., 0., 11.]);
    let one = apply(&stock, &boss(20., 2.), BooleanOperation::Union).unwrap();
    let outside = |r: f64| PI * r * r * 3.5;
    near(one.volume_mm3().unwrap(), 21600. + outside(2.));
    // The cylinder may be either operand.
    let two = apply(&boss(60., 3.), &one, BooleanOperation::Union).unwrap();
    near(
        two.volume_mm3().unwrap(),
        21600. + outside(2.) + outside(3.),
    );
    let bounds = two.bounds().unwrap();
    near(bounds.max[0], 13.5);

    // A placed retained union takes the placed boss in its own chart.
    let pose = RigidTransform {
        rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
        translation_mm: [13., -7., 42.],
    };
    let placed = one.transformed(pose).unwrap();
    let placed_two = apply(
        &placed,
        &boss(60., 3.).transformed(pose).unwrap(),
        BooleanOperation::Union,
    )
    .unwrap();
    near(placed_two.volume_mm3().unwrap(), two.volume_mm3().unwrap());
}

#[test]
fn unproved_pairs_are_left_to_the_caller() {
    let policy = policy();
    let stock = block([10., 80., 27.], [0., 0., 11.]);
    let declined =
        |a: &Solid, b: &Solid, operation| a.pairwise_warped_boolean(b, operation, policy).is_none();
    // A buried cylinder has no exit face; a crossing one exits twice.
    let buried = Solid::cylinder(
        PlaneFrame {
            origin: [2., 20., 20.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        },
        2.,
        3.,
        policy,
    )
    .unwrap();
    assert!(declined(&stock, &buried, BooleanOperation::Union));
    let crossing = Solid::cylinder(
        PlaneFrame {
            origin: [-1., 20., 20.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        },
        2.,
        12.,
        policy,
    )
    .unwrap();
    assert!(declined(&stock, &crossing, BooleanOperation::Union));
    // Plain planar pairs, cylinder pairs and intersections belong elsewhere.
    let other = block([2., 2., 2.], [20., 0., 0.]);
    assert!(declined(&stock, &other, BooleanOperation::Union));
    assert!(declined(
        &boss(20., 2.),
        &boss(60., 2.),
        BooleanOperation::Union
    ));
    assert!(declined(
        &stock,
        &boss(20., 2.),
        BooleanOperation::Intersect
    ));
    // A retained union with a boss takes no further planar material, and
    // plain stock subtracts through other providers.
    let one = apply(&stock, &boss(20., 2.), BooleanOperation::Union).unwrap();
    assert!(declined(&one, &other, BooleanOperation::Union));
    assert!(declined(&stock, &other, BooleanOperation::Subtract));
    // A fillet clip whose wall leans through the rounded band does not keep
    // its plane in the chart.
    let leaning = block([4., 4., 4.], [-1., -1., -1.])
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [0., 0.8, -0.6], [0., 0.6, 0.8]],
            ..Default::default()
        })
        .unwrap();
    assert!(declined(
        &bare_fillet(),
        &leaning,
        BooleanOperation::Subtract
    ));
}

#[test]
fn open_fillet_recognition_keeps_the_strict_reader_strict() {
    let tolerance = policy().linear_tolerance_mm;
    let fillet = bare_fillet();
    assert!(QuarterFillet::recognize(&fillet, tolerance).is_none());
    let open = QuarterFillet::recognize_open(&fillet, tolerance).unwrap();
    near(open.radius, 2.);
    assert_eq!(open.pieces(2., 2., tolerance).unwrap().len(), 1);
    // A plain box is no fillet.
    assert!(QuarterFillet::recognize_open(&block([1., 1., 1.], [0.; 3]), tolerance).is_none());
}
