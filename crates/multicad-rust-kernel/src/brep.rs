//! Retained boundary representation of an imported solid.
//!
//! The model keeps the support surfaces, trim loops, edge curves and
//! vertices exactly as the exchange file states them. Mass properties are
//! integrated over the trimmed parameter domains with the exact surface
//! evaluators; tessellation, sections and bounds read the same boundary.

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
        sine.atan2(cosine)
    }
}

pub(crate) fn distance(a: Point, b: Point) -> f64 {
    norm(sub(a, b))
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
