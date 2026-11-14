//! Exact view-space paths shared by the scene members.
use super::*;
use multicad_rust_kernel::{curved_projection::Ellipse, normal_sweep, spur_gear::Edge};
use std::f64::consts::{PI, TAU};

pub(in crate::drawing) fn add(a: P3, b: P3) -> P3 {
    std::array::from_fn(|i| a[i] + b[i])
}

pub(in crate::drawing) fn scaled(a: P3, s: f64) -> P3 {
    a.map(|v| v * s)
}

pub(in crate::drawing) fn dot3(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(in crate::drawing) fn dot2(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

/// A straight view-space segment.
pub(in crate::drawing) struct Segment(pub P3, pub P3);

impl Path for Segment {
    fn point(&self, t: f64) -> P3 {
        std::array::from_fn(|i| self.0[i] + (self.1[i] - self.0[i]) * t)
    }

    fn flatten(&self, _: f64) -> Vec<f64> {
        vec![0., 1.]
    }

    fn paper(&self, [a, b]: [f64; 2]) -> Option<Curve> {
        let [p, q] = [self.point(a), self.point(b)];
        Some(Curve::Line([[p[0], p[1]], [q[0], q[1]]]))
    }
}

/// Affine placement of a planar chart in view space.
#[derive(Clone, Copy, Debug)]
pub(in crate::drawing) struct Lift {
    pub origin: P3,
    pub x: P3,
    pub y: P3,
}

impl Lift {
    pub(in crate::drawing) fn point(&self, [x, y]: P2) -> P3 {
        std::array::from_fn(|i| self.origin[i] + x * self.x[i] + y * self.y[i])
    }

    fn stretch(&self) -> f64 {
        (dot3(self.x, self.x) + dot3(self.y, self.y)).sqrt().max(1.)
    }
}

/// Center, radius, start angle and signed sweep of a circular profile edge.
pub(in crate::drawing) fn round(edge: normal_sweep::Edge) -> Option<(P2, f64, f64, f64)> {
    match edge {
        normal_sweep::Edge::Line { .. } => None,
        normal_sweep::Edge::Arc { center, a, .. } => {
            let [x, y] = [a[0] - center[0], a[1] - center[1]];
            Some((center, x.hypot(y), y.atan2(x), edge.signed_sweep_rad()?))
        }
        normal_sweep::Edge::Circle {
            center,
            radius_mm,
            clockwise,
        } => Some((center, radius_mm, 0., if clockwise { -TAU } else { TAU })),
    }
}

/// Exact evaluation of a profile edge that stays accurate at any fraction.
pub(in crate::drawing) fn edge_point(edge: Edge, fraction: f64) -> P2 {
    // The authored ends are shared exactly with the neighboring edges.
    if let Edge::Planar(planar) = edge
        && fraction != 0.
        && fraction != 1.
        && let Some((center, radius, start, sweep)) = round(planar)
    {
        let (sin, cos) = (start + fraction * sweep).sin_cos();
        return [center[0] + radius * cos, center[1] + radius * sin];
    }
    edge.point(fraction)
}

/// A range of one profile edge placed by an affine lift.
pub(in crate::drawing) struct ProfilePath {
    pub edge: Edge,
    pub range: [f64; 2],
    pub lift: Lift,
    /// The paper image is a straight segment traversed monotonically.
    pub straight: bool,
}

impl ProfilePath {
    fn fraction(&self, t: f64) -> f64 {
        self.range[0] + (self.range[1] - self.range[0]) * t
    }
}

/// The paper ellipse `center + cos(a) * cosine + sin(a) * sine` over the
/// angles from `start` sweeping by `sweep`, in the published convention: an
/// increasing trim on a counterclockwise orthonormal basis.
pub(in crate::drawing) fn conic(
    center: P2,
    cosine: P2,
    sine: P2,
    start: f64,
    sweep: f64,
) -> Option<Curve> {
    let turn = 0.5 * (2. * dot2(cosine, sine)).atan2(dot2(cosine, cosine) - dot2(sine, sine));
    let (sin, cos) = turn.sin_cos();
    let major = [
        cos * cosine[0] + sin * sine[0],
        cos * cosine[1] + sin * sine[1],
    ];
    let minor = [
        cos * sine[0] - sin * cosine[0],
        cos * sine[1] - sin * cosine[1],
    ];
    let [a, b] = [major, minor].map(|axis| axis[0].hypot(axis[1]));
    if !(a.is_finite() && b > 1e-9 * a.max(1.)) {
        return None;
    }
    let u = [major[0] / a, major[1] / a];
    // Mirror the angle when the minor axis is clockwise from the major axis.
    let mirrored = u[0] * minor[1] - u[1] * minor[0] < 0.;
    let mut from = start - turn;
    let mut sweep = sweep;
    if mirrored {
        from = -from;
        sweep = -sweep;
    }
    if sweep < 0. {
        from += sweep;
        sweep = -sweep;
    }
    let (u, from) = if (a - b).abs() <= 1e-12 * a {
        // A circle rotates its basis so that the trim starts at zero.
        let (sin, cos) = from.sin_cos();
        ([cos * u[0] - sin * u[1], sin * u[0] + cos * u[1]], 0.)
    } else {
        (u, from.rem_euclid(TAU))
    };
    Some(Curve::Ellipse(Ellipse {
        center,
        u,
        v: [-u[1], u[0]],
        major: a,
        minor: b,
        trim: [from, from + sweep.min(TAU)],
    }))
}

impl Path for ProfilePath {
    fn point(&self, t: f64) -> P3 {
        self.lift.point(edge_point(self.edge, self.fraction(t)))
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        let span = (self.range[1] - self.range[0]).abs();
        let tolerance = tolerance_mm / self.lift.stretch();
        let count = match self.edge {
            Edge::Planar(edge) => match round(edge) {
                None => 1.,
                Some((_, radius, _, sweep)) => {
                    let step = if tolerance < radius {
                        2. * (1. - tolerance / radius).acos()
                    } else {
                        PI
                    };
                    (sweep.abs() * span / step.min(PI / 3.)).ceil()
                }
            },
            Edge::Involute(curve) => {
                let [a, b] = curve.interval();
                let second = curve.base_radius_mm() * a.max(b).hypot(1.);
                ((b - a).abs() * span / (8. * tolerance / second).sqrt().min(0.3)).ceil()
            }
        };
        let count = (count as usize).clamp(1, 1 << 20);
        (0..=count).map(|i| i as f64 / count as f64).collect()
    }

    fn paper(&self, [from, to]: [f64; 2]) -> Option<Curve> {
        let ends = [from, to].map(|t| {
            let p = self.point(t);
            [p[0], p[1]]
        });
        let Edge::Planar(edge) = self.edge else {
            return self.straight.then_some(Curve::Line(ends));
        };
        let Some((center, radius, start, sweep)) = round(edge) else {
            return Some(Curve::Line(ends));
        };
        if self.straight {
            return Some(Curve::Line(ends));
        }
        let c = self.lift.point(center);
        let [fa, fb] = [self.fraction(from), self.fraction(to)];
        conic(
            [c[0], c[1]],
            [radius * self.lift.x[0], radius * self.lift.x[1]],
            [radius * self.lift.y[0], radius * self.lift.y[1]],
            start + fa * sweep,
            (fb - fa) * sweep,
        )
    }
}
