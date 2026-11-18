use super::points::*;
use multicad_core::{
    document::Document,
    model::{self, NumericTransform},
    parameters::EvaluatedParameters,
};

fn evaluated() -> EvaluatedParameters {
    multicad_core::parameters::evaluate(&Document::new("Point to point").unwrap()).unwrap()
}

#[test]
fn alignment_translates_a_rotated_staged_pose_so_the_source_lands_on_the_target() {
    let p = evaluated();
    let staged = NumericTransform {
        translation: [model::length(3.), model::length(-2.), model::length(1.)],
        pivot: [model::length(5.), model::length(7.), model::length(0.)],
        axis: [0., 0., 1.],
        angle: model::angle(30.),
    };
    let source = [20., 10., 5.];
    let target = [61.5, 17.25, -3.];
    let result = aligned(&staged, &p, source, target).unwrap();
    let placed = result.resolve(&p).unwrap().point(source);
    for i in 0..3 {
        assert!((placed[i] - target[i]).abs() < 1e-9, "{placed:?}");
    }
    // Only translation changes: orientation, axis and pivot are kept.
    assert_eq!(result.angle, staged.angle);
    assert_eq!(result.axis, staged.axis);
    assert_eq!(result.pivot, staged.pivot);
    let before = staged.resolve(&p).unwrap();
    let after = result.resolve(&p).unwrap();
    assert_eq!(before.rotation, after.rotation);
}

#[test]
fn coincident_points_leave_the_placement_unchanged_and_nonfinite_points_fail() {
    let p = evaluated();
    let staged = NumericTransform {
        translation: [model::length(4.), model::length(0.), model::length(0.)],
        ..Default::default()
    };
    let same = aligned(&staged, &p, [1., 2., 3.], [5., 2., 3.]).unwrap();
    assert_eq!(same, staged);
    assert!(aligned(&staged, &p, [f64::NAN, 0., 0.], [0.; 3]).is_err());
    assert!(aligned(&staged, &p, [0.; 3], [0., f64::INFINITY, 0.]).is_err());
}
