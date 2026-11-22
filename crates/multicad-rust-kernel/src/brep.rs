//! Retained boundary representation of an imported solid.
//!
//! The model keeps the support surfaces, trim loops, edge curves and
//! vertices exactly as the exchange file states them. Mass properties are
//! integrated over the trimmed parameter domains with the exact surface
//! evaluators; tessellation, sections and bounds read the same boundary.

use multicad_precision::PortableMath as _;
pub mod affine;
mod containment;
pub mod curve;
mod domain;
mod extent;
mod mass;
pub mod mesh;
pub mod nurbs;
mod polygon;
pub mod section;
pub mod solid;
pub mod surface;
pub mod vertex;

pub use affine::Affine;
pub use containment::BoxSide;
pub use curve::Curve;
pub use nurbs::{NurbsCurve, NurbsSurface};
pub use solid::{Boundary, Edge, Face, FaceProperties, Loop, Shell, Solid};
pub use surface::{Frame, Surface};

pub type Point = [f64; 3];

/// Refusal classes of the retained boundary provider.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    /// The boundary is malformed, open, inconsistently oriented or lies
    /// outside the linear tolerance.
    Invalid,
    /// The boundary is well formed but outside the supported families.
    Unsupported,
    /// A resource budget was exhausted.
    Limit,
    /// The caller's interrupt asked the computation to stop.
    Cancelled,
}

pub type Result<T> = std::result::Result<T, Error>;

/// A caller-supplied poll that returns true once the work is superseded.
pub type Interrupt<'a> = Option<&'a dyn Fn() -> bool>;

pub(crate) fn ensure(ok: bool, error: Error) -> Result<()> {
    if ok { Ok(()) } else { Err(error) }
}

pub(crate) fn checkpoint(interrupt: Interrupt<'_>) -> Result<()> {
    ensure(!interrupt.is_some_and(|stop| stop()), Error::Cancelled)
}

pub(crate) fn add(a: Point, b: Point) -> Point {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub(crate) fn sub(a: Point, b: Point) -> Point {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(crate) fn scale(a: Point, s: f64) -> Point {
    [a[0] * s, a[1] * s, a[2] * s]
}

pub(crate) fn dot(a: Point, b: Point) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(crate) fn norm(a: Point) -> f64 {
    dot(a, a).sqrt()
}

/// Angle between two directions, in `[0, pi]`; zero when either vanishes.
pub(crate) fn turn(a: Point, b: Point) -> f64 {
    let (sine, cosine) = (norm(cross(a, b)), dot(a, b));
    if sine == 0. && cosine == 0. {
        0.
    } else {
        sine.portable_atan2(cosine)
    }
}

pub(crate) fn distance(a: Point, b: Point) -> f64 {
    norm(sub(a, b))
}

/// Distance from `p` to the closest point of a triangle.
pub(crate) fn point_triangle(p: Point, [a, b, c]: [Point; 3]) -> f64 {
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    let (va, vb, vc) = (d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2);
    let closest = if d1 <= 0. && d2 <= 0. {
        a
    } else if d3 >= 0. && d4 <= d3 {
        b
    } else if d6 >= 0. && d5 <= d6 {
        c
    } else if vc <= 0. && d1 >= 0. && d3 <= 0. {
        add(a, scale(ab, d1 / (d1 - d3)))
    } else if vb <= 0. && d2 >= 0. && d6 <= 0. {
        add(a, scale(ac, d2 / (d2 - d6)))
    } else if va <= 0. && d4 >= d3 && d5 >= d6 {
        add(b, scale(sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6))))
    } else {
        let total = va + vb + vc;
        add(a, add(scale(ab, vb / total), scale(ac, vc / total)))
    };
    distance(p, closest)
}

/// Whether a ray from `origin` along a fixed skew direction passes through
/// the triangle.
pub(crate) fn crosses_ray(origin: Point, triangle: &[Point; 3]) -> bool {
    const RAY: Point = [
        0.574_352_934_563_6,
        0.577_919_283_817_1,
        0.579_868_302_102_6,
    ];
    let [a, b, c] = *triangle;
    let (ab, ac) = (sub(b, a), sub(c, a));
    let p = cross(RAY, ac);
    let determinant = dot(ab, p);
    if determinant == 0. {
        return false;
    }
    let t = sub(origin, a);
    let u = dot(t, p) / determinant;
    let q = cross(t, ab);
    let v = dot(RAY, q) / determinant;
    u >= 0. && v >= 0. && u + v <= 1. && dot(ac, q) / determinant > 0.
}

/// Eight-point Gauss-Legendre rule on `[-1, 1]`.
pub(crate) const GAUSS: [(f64, f64); 8] = [
    (-0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    (-0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (-0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (-0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];

/// Shift `value` by whole periods to the representative nearest `target`.
pub(crate) fn nearest(value: f64, target: f64, period: Option<f64>) -> f64 {
    match period {
        Some(period) => value + ((target - value) / period).round() * period,
        None => value,
    }
}

#[cfg(test)]
mod tests;
