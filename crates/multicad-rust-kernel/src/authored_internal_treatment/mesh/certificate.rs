//! Whole actual source patch versus each literal emitted triangle.
use super::*;
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |sum, (a, b)| sum + a * b)
}
pub(super) fn chord(field: &source::Field, errors: [f64; 3]) -> f64 {
    let curvature = field.point.iter().fold(I::point(0.), |sum, p| {
        sum + I::point(p.second[0].magnitude())
            + I::point(2.) * I::point(p.second[1].magnitude())
            + I::point(p.second[2].magnitude())
    });
    // Taylor expansion at the barycentric source point cancels first
    // derivatives. Each coordinate variance is <=1/4; Cauchy bounds the
    // absolute mixed covariance by1/4 too. Thus the remainder factor is1/8.
    // Shared-vertex source-to-literal discrepancy is a convex combination.
    (I::point(0.125) * curvature + I::point(errors.into_iter().fold(0., f64::max))).hi
}
pub(super) fn oriented_triangle(
    points: [[f64; 3]; 3],
    normal: source::Normal,
    angle: f64,
) -> Option<(bool, f64)> {
    let p = points.map(|p| p.map(I::point));
    let edge = |i: usize| std::array::from_fn(|k| p[i][k] - p[0][k]);
    let facet = cross(edge(1), edge(2));
    // Yaw is a finite sampling proposal only. Outward exact-real rotation of
    // BOTH vectors preserves angle/orientation, while evaluating theta-yaw
    // directly avoids cancellation between independent world sine/cosine hulls.
    let yaw = facet[1]
        .lo
        .midpoint(facet[1].hi)
        .portable_atan2(facet[0].lo.midpoint(facet[0].hi));
    // Actual-source vertices are bounded by 1e7 in each coordinate. Their
    // finite unnormalized cross products stay below 1e15; atan2(0,0) is finite.
    let a = I::point(yaw);
    let c = a.cos();
    let s = a.sin();
    let mut facet = [
        facet[0] * c + facet[1] * s,
        -facet[0] * s + facet[1] * c,
        facet[2],
    ];
    let normal = normal.components(yaw);
    let direction = dot(facet, normal);
    let reverse = direction.hi < 0.;
    if reverse {
        facet = facet.map(|v| -v);
    }
    let projection = dot(facet, normal);
    if !projection.valid() || projection.lo <= 0. {
        return None;
    }
    let error = cross(facet, normal);
    let magnitude = dot(error, error).sqrt();
    let bound = (I::new(0., magnitude.hi) / I::point(projection.lo)).atan();
    (bound.valid() && bound.hi <= angle).then_some((reverse, bound.hi))
}
