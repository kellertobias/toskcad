use super::*;
use std::f64::consts::PI;

fn r(n: i64) -> Rational {
    Rational::new(n, 1).unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.max(1.),
        "{actual} != {expected}"
    );
}

/// An 80 x 60 plate with one circular hole, of `height` swept from `z`.
fn plate(z: f64, height: i64, reverse: bool, hole: i64) -> Solid {
    let feature = EntityId::new();
    let ring = [[0, 0], [80, 0], [80, 60], [0, 60]]
        .map(|p| p.map(r))
        .to_vec();
    let prism = PolygonPrism::with_circular_holes(
        vec![ring],
        vec![Circle::new([r(40), r(30)], r(hole)).unwrap()],
        r(height),
    )
    .unwrap();
    let origins: Vec<_> = (0..5)
        .map(|i| Origin {
            source: feature,
            role: format!("edge/{i}"),
        })
        .collect();
    Solid::from_polygon_prism_sweep(
        prism,
        PlaneFrame {
            origin: [0., 0., z],
            ..Default::default()
        },
        reverse,
        feature,
        &origins,
        &[false; 5],
        GeometryPolicy::default(),
    )
    .unwrap()
}

const AREA: f64 = 4800. - 25. * PI;

fn run(a: &Solid, b: &Solid, operation: BooleanOperation) -> Option<Solid> {
    a.boolean(b, operation, GeometryPolicy::default()).unwrap()
}

fn z_range(solid: &Solid) -> [f64; 2] {
    let bounds = solid.bounds().unwrap();
    [bounds.min[2], bounds.max[2]]
}

#[test]
fn touching_and_overlapping_same_profile_prisms_join_into_one_prism() {
    let base = plate(0., 10, false, 5);
    let top = plate(10., 5, false, 5);
    let joined = run(&base, &top, BooleanOperation::Union).unwrap();
    near(joined.volume_mm3().unwrap(), AREA * 15.);
    assert_eq!(z_range(&joined), [0., 15.]);
    assert!(joined.polygon_prism().is_some());
    // A reversed sweep hanging below the base overlaps it by 2 mm.
    let below = plate(2., 6, true, 5);
    let joined = run(&base, &below, BooleanOperation::Union).unwrap();
    near(joined.volume_mm3().unwrap(), AREA * 14.);
    assert_eq!(z_range(&joined), [-4., 10.]);
}

#[test]
fn end_cuts_and_common_intervals_keep_the_exact_profile() {
    let base = plate(0., 10, false, 5);
    // A face-profile cut 4 mm down from the top.
    let cut = run(&base, &plate(10., 4, true, 5), BooleanOperation::Subtract).unwrap();
    near(cut.volume_mm3().unwrap(), AREA * 6.);
    assert_eq!(z_range(&cut), [0., 6.]);
    let bottom = run(&base, &plate(-3., 5, false, 5), BooleanOperation::Subtract).unwrap();
    assert_eq!(z_range(&bottom), [2., 10.]);
    let common = run(&base, &plate(7., 9, false, 5), BooleanOperation::Intersect).unwrap();
    near(common.volume_mm3().unwrap(), AREA * 3.);
    assert_eq!(z_range(&common), [7., 10.]);
    assert!(run(&base, &plate(-1., 12, false, 5), BooleanOperation::Subtract).is_none());
}

#[test]
fn splitting_cuts_and_different_profiles_are_not_stacked_prisms() {
    let policy = GeometryPolicy::default();
    let base = plate(0., 10, false, 5);
    let middle = plate(3., 2, false, 5);
    assert!(
        base.stacked_prism_boolean(&middle, BooleanOperation::Subtract, policy)
            .is_none()
    );
    let other_hole = plate(10., 5, false, 6);
    assert!(
        base.stacked_prism_boolean(&other_hole, BooleanOperation::Union, policy)
            .is_none()
    );
    let apart = plate(11., 5, false, 5);
    assert!(
        base.stacked_prism_boolean(&apart, BooleanOperation::Union, policy)
            .is_none()
    );
    let shifted = plate(10., 5, false, 5)
        .transformed(RigidTransform {
            translation_mm: [1., 0., 0.],
            ..Default::default()
        })
        .unwrap();
    assert!(
        base.stacked_prism_boolean(&shifted, BooleanOperation::Union, policy)
            .is_none()
    );
}
