//! Exact evaluation of one transverse profile edge and its derivatives.
use super::super::paths::{dot2, round};
use super::super::*;
use multicad_rust_kernel::{normal_sweep, spur_gear::Edge};

#[derive(Clone, Copy, Debug)]
pub(in crate::drawing) enum Shape {
    Line {
        a: P2,
        b: P2,
    },
    Arc {
        center: P2,
        radius: f64,
        start: f64,
        sweep: f64,
    },
    /// `base * (u - t v)` in the frame `x`, `y` with `t = start + s * span`.
    Involute {
        base: f64,
        start: f64,
        span: f64,
        x: P2,
        y: P2,
    },
}

/// One directed profile edge, material on its left.
#[derive(Clone, Copy, Debug)]
pub(in crate::drawing) struct Flank {
    pub edge: Edge,
    pub shape: Shape,
    /// Least and greatest distance from the axis.
    pub radius: [f64; 2],
    /// Bounds of the position and its first three fraction derivatives.
    pub size: [f64; 4],
    /// Polar angle of the edge's middle and a bound of the edge's angular
    /// distance from it.
    pub azimuth: [f64; 2],
}

pub(super) fn turn([x, y]: P2) -> P2 {
    [-y, x]
}

pub(super) fn rotate([x, y]: P2, angle: f64) -> P2 {
    let (sin, cos) = angle.sin_cos();
    [x * cos - y * sin, x * sin + y * cos]
}

fn polar(angle: f64) -> P2 {
    let (sin, cos) = angle.sin_cos();
    [cos, sin]
}

impl Flank {
    pub(super) fn new(edge: Edge) -> Self {
        let shape = match edge {
            Edge::Planar(planar) => match round(planar) {
                Some((center, radius, start, sweep)) => Shape::Arc {
                    center,
                    radius,
                    start,
                    sweep,
                },
                None => {
                    let normal_sweep::Edge::Line { a, b } = planar else {
                        unreachable!("only lines have no radius")
                    };
                    Shape::Line { a, b }
                }
            },
            Edge::Involute(curve) => {
                let [x, y] = curve.basis();
                let [a, b] = curve.interval();
                Shape::Involute {
                    base: curve.base_radius_mm(),
                    start: a,
                    span: b - a,
                    x,
                    y,
                }
            }
        };
        let length = |p: P2| p[0].hypot(p[1]);
        let (radius, size) = match shape {
            Shape::Line { a, b } => {
                let d = [b[0] - a[0], b[1] - a[1]];
                let l = length(d);
                let foot = (-dot2(a, d) / (l * l)).clamp(0., 1.);
                let near = length([a[0] + foot * d[0], a[1] + foot * d[1]]);
                let far = length(a).max(length(b));
                ([near, far], [far, l, 0., 0.])
            }
            Shape::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                let offset = length(center);
                let mut lo = f64::INFINITY;
                let mut hi = 0_f64;
                let mut visit = |angle: f64| {
                    let [x, y] = polar(angle);
                    let r = length([center[0] + radius * x, center[1] + radius * y]);
                    lo = lo.min(r);
                    hi = hi.max(r);
                };
                visit(start);
                visit(start + sweep);
                if offset > 0. {
                    let toward = center[1].atan2(center[0]);
                    for extreme in [toward, toward + std::f64::consts::PI] {
                        let travel =
                            ((extreme - start) * sweep.signum()).rem_euclid(std::f64::consts::TAU);
                        if travel < sweep.abs() {
                            visit(extreme);
                        }
                    }
                }
                let w = sweep.abs();
                (
                    [lo, hi],
                    [
                        offset + radius,
                        radius * w,
                        radius * w * w,
                        radius * w * w * w,
                    ],
                )
            }
            Shape::Involute {
                base, start, span, ..
            } => {
                let [a, b] = [start, start + span];
                let t = a.abs().max(b.abs());
                let w = span.abs();
                (
                    [base * a.hypot(1.), base * b.hypot(1.)],
                    [
                        base * t.hypot(1.),
                        base * t * w,
                        base * t.hypot(1.) * w * w,
                        base * t.hypot(2.) * w * w * w,
                    ],
                )
            }
        };
        let mut flank = Self {
            edge,
            shape,
            radius: [radius[0].min(radius[1]), radius[0].max(radius[1])],
            size,
            azimuth: [0., std::f64::consts::PI],
        };
        // The edge stays within half its speed bound of its middle point.
        let middle = flank.point(0.5);
        let reach = 0.5 * size[1];
        let distance = length(middle);
        if reach < distance {
            flank.azimuth = [middle[1].atan2(middle[0]), (reach / distance).asin()];
        }
        flank
    }

    /// Position and its first three derivatives by the edge fraction.
    pub(super) fn at(&self, s: f64) -> [P2; 4] {
        match self.shape {
            Shape::Line { a, b } => {
                let d = [b[0] - a[0], b[1] - a[1]];
                [[a[0] + s * d[0], a[1] + s * d[1]], d, [0.; 2], [0.; 2]]
            }
            Shape::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                let e = polar(start + s * sweep);
                let j = turn(e);
                let [r1, r2, r3] = [
                    radius * sweep,
                    radius * sweep * sweep,
                    radius * sweep * sweep * sweep,
                ];
                [
                    [center[0] + radius * e[0], center[1] + radius * e[1]],
                    [r1 * j[0], r1 * j[1]],
                    [-r2 * e[0], -r2 * e[1]],
                    [-r3 * j[0], -r3 * j[1]],
                ]
            }
            Shape::Involute {
                base,
                start,
                span,
                x,
                y,
            } => {
                let t = start + s * span;
                let (sin, cos) = t.sin_cos();
                let map = |[p, q]: P2, scale: f64| {
                    [scale * (p * x[0] + q * y[0]), scale * (p * x[1] + q * y[1])]
                };
                let [u, v] = [[cos, sin], [-sin, cos]];
                [
                    map([u[0] - t * v[0], u[1] - t * v[1]], base),
                    map([t * u[0], t * u[1]], base * span),
                    map([u[0] + t * v[0], u[1] + t * v[1]], base * span * span),
                    map(
                        [2. * v[0] - t * u[0], 2. * v[1] - t * u[1]],
                        base * span * span * span,
                    ),
                ]
            }
        }
    }

    pub(super) fn point(&self, s: f64) -> P2 {
        self.at(s)[0]
    }

    /// Unit tangent, including the one-sided limit at an involute's base.
    pub(super) fn tangent(&self, s: f64) -> P2 {
        match self.shape {
            Shape::Involute {
                start, span, x, y, ..
            } => {
                let (sin, cos) = (start + s * span).sin_cos();
                let sign = span.signum();
                [
                    sign * (cos * x[0] + sin * y[0]),
                    sign * (cos * x[1] + sin * y[1]),
                ]
            }
            _ => {
                let d = self.at(s)[1];
                let l = d[0].hypot(d[1]);
                [d[0] / l, d[1] / l]
            }
        }
    }

    /// Radius of an arc centered on the gear axis.
    pub(super) fn coaxial(&self) -> Option<f64> {
        match self.shape {
            Shape::Arc { center, radius, .. } if center == [0.; 2] => Some(radius),
            _ => None,
        }
    }

    /// Angle of the outward normal at the start and its total turn, for the
    /// shapes whose normal turns uniformly with the fraction.
    pub(super) fn normal_turn(&self) -> Option<[f64; 2]> {
        let [x, y] = self.tangent(0.);
        let start = (-x).atan2(y);
        match self.shape {
            Shape::Line { .. } => None,
            Shape::Arc { sweep, .. } => Some([start, sweep]),
            Shape::Involute { span, x, y, .. } => {
                let mirror = (x[0] * y[1] - x[1] * y[0]).signum();
                Some([start, mirror * span])
            }
        }
    }
}
