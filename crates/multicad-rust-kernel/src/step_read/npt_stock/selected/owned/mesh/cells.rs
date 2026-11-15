//! Whole-cell correspondence to actual literal triangles, in the STEP frame.
use super::*;
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn length(a: [I; 3]) -> I {
    a.into_iter()
        .fold(I::point(0.), |sum, x| sum + x.square())
        .sqrt()
}
pub(super) fn l1(a: [I; 3]) -> f64 {
    a.into_iter()
        .fold(I::point(0.), |sum, x| sum + I::point(x.magnitude()))
        .hi
}
pub(super) fn vector(source: &OwnedMaterial, v: [[f64; 2]; 3]) -> [I; 3] {
    let v = v.map(interval);
    std::array::from_fn(|axis| {
        (0..3).fold(I::point(0.), |sum, j| sum + source.0.linear[j][axis] * v[j])
    })
}
fn facet_normal(p: [[f64; 3]; 3]) -> [I; 3] {
    let [a, b, c] = p.map(|p| p.map(I::point));
    cross(
        std::array::from_fn(|k| b[k] - a[k]),
        std::array::from_fn(|k| c[k] - a[k]),
    )
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |sum, (a, b)| sum + a * b)
}
fn cone(corners: [[I; 3]; 4], error: [I; 3], facet: [I; 3]) -> f64 {
    let perpendicular = length(cross(error, facet)).hi;
    let axial = dot(error, facet).magnitude();
    let mut bound = 0f64;
    for point in corners {
        let denominator = dot(point, facet) - I::point(axial);
        if !denominator.valid() || denominator.lo <= 0. {
            return f64::INFINITY;
        }
        let numerator = I::point(length(cross(point, facet)).hi) + I::point(perpendicular);
        let value = numerator / I::point(denominator.lo);
        if !value.valid() {
            return f64::INFINITY;
        }
        bound = bound.max(value.hi);
    }
    // N=sum B_i*N_i+E, B_i>=0,sum B_i=1. Each corner obeys
    // ||N_i x F||+max||E x F|| <= r*(N_i.F-max|E.F|).
    // Convexity gives the same cone for EVERY original source normal.
    // Positive dot earns theta=atan(r)<=r without a libm angle assumption.
    bound
}
pub(super) fn cap_angle(source: &OwnedMaterial, upward: bool, triangle: [[f64; 3]; 3]) -> f64 {
    let u = source.0.linear[0];
    let v = source.0.linear[1];
    let normal = if upward { cross(u, v) } else { cross(v, u) };
    cone([normal; 4], [I::point(0.); 3], facet_normal(triangle))
}
