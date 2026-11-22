use super::*;
use crate::{BSpline2d, GeometryPolicy};

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 1e-10, "{actual} != {expected}");
}

fn near_point(actual: Point, expected: Point) {
    near(actual[0], expected[0]);
    near(actual[1], expected[1]);
}

#[test]
fn full_circle_has_one_seam_and_continuous_projections_fail_explicitly() {
    let policy = GeometryPolicy::default();
    let circle = CurvedEdge2d::Circle {
        center: [0.; 2],
        radius_mm: 2.,
    };
    for point in [[3., 0.], [3., -1e-20]] {
        let projected = circle.project_parameters(point, policy).unwrap();
        assert_eq!(projected.len(), 1);
        assert_eq!(projected[0].parameter, 0.);
        assert_eq!(projected[0].point, [2., 0.]);
    }
    let adjacent = circle.project_parameters([3., -1e-6], policy).unwrap();
    assert!((TAU - 1e-6..TAU).contains(&adjacent[0].parameter));
    for point in [[0.; 2], [policy.linear_tolerance_mm * 0.5, 0.]] {
        assert_eq!(
            circle.project_parameters(point, policy).unwrap_err().kind,
            GeometryErrorKind::Input
        );
    }
    let arc = CurvedEdge2d::Arc {
        center: [0.; 2],
        a: [2., 0.],
        b: [0., 2.],
        clockwise: false,
    };
    assert_eq!(
        arc.project_parameters([0.; 2], policy).unwrap_err().kind,
        GeometryErrorKind::Input
    );
    let near_center = arc.project_parameters([1e-9, 1e-9], policy).unwrap();
    assert_eq!(near_center.len(), 3);
}

#[test]
fn splines_keep_their_authored_parameter_domain() {
    let policy = GeometryPolicy::default();
    let spline = CurvedEdge2d::BSpline(BSpline2d {
        degree: 1,
        poles: vec![[0., 0.], [3., 4.]],
        knots: vec![2., 5.],
        multiplicities: vec![2, 2],
        weights: None,
    });
    let sample = spline.evaluate(3.5, policy).unwrap();
    near_point(sample.point, [1.5, 2.]);
    let hits = spline.project_parameters(sample.point, policy).unwrap();
    assert_eq!(hits.len(), 1);
    near(hits[0].parameter, 3.5);
}

#[test]
fn rational_quarter_circle_evaluation_matches_its_analytic_tangent() {
    let policy = GeometryPolicy::default();
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let spline = BSpline2d {
        degree: 2,
        poles: vec![[1., 0.], [1., 1.], [0., 1.]],
        knots: vec![0., 1.],
        multiplicities: vec![3, 3],
        weights: Some(vec![1., w, 1.]),
    };
    let middle = spline.evaluate(0.5, policy).unwrap();
    near_point(middle.point, [w, w]);
    near_point(middle.derivative, [-2. / (1. + w), 2. / (1. + w)]);
    near_point(
        spline.evaluate(0., policy).unwrap().derivative,
        [0., 2. * w],
    );
    near_point(
        spline.evaluate(1., policy).unwrap().derivative,
        [-2. * w, 0.],
    );
}
