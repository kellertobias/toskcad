//! Face support surfaces in their exchange parameterisation.
use super::curve::{PANEL, uniform};
use super::{Curve, NurbsSurface, Point, add, cross, distance, dot, nearest, norm, scale, sub};
use multicad_precision::PortableMath as _;
use std::f64::consts::{FRAC_PI_2, PI, TAU};
mod prepared_inverse;
use prepared_inverse::Candidate;
pub(crate) use prepared_inverse::PreparedInverse;

/// A right-handed orthonormal placement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    pub origin: Point,
    pub x: Point,
    pub y: Point,
    pub z: Point,
}

impl Frame {
    /// Coordinates of `point` along the frame axes.
    pub fn local(&self, point: Point) -> Point {
        let d = sub(point, self.origin);
        [dot(d, self.x), dot(d, self.y), dot(d, self.z)]
    }

    pub(crate) fn mapped(
        &self,
        map: &dyn Fn(Point) -> Point,
        linear: &dyn Fn(Point) -> Point,
    ) -> Self {
        Self {
            origin: map(self.origin),
            x: linear(self.x),
            y: linear(self.y),
            z: linear(self.z),
        }
    }

    fn radial(&self, u: f64) -> (Point, Point) {
        let (sin, cos) = u.portable_sin_cos();
        (
            add(scale(self.x, cos), scale(self.y, sin)),
            add(scale(self.x, -sin), scale(self.y, cos)),
        )
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Surface {
    /// `origin + u * x + v * y`.
    Plane(Frame),
    /// `origin + r * radial(u) + v * z`.
    Cylinder { frame: Frame, radius: f64 },
    /// `origin + (r + v * tan(half_angle)) * radial(u) + v * z`.
    Cone {
        frame: Frame,
        radius: f64,
        half_angle: f64,
    },
    /// `origin + r * cos v * radial(u) + r * sin v * z`.
    Sphere { frame: Frame, radius: f64 },
    /// `origin + (major + minor * cos v) * radial(u) + minor * sin v * z`.
    Torus {
        frame: Frame,
        major: f64,
        minor: f64,
    },
    /// `curve(u) + v * direction`.
    Extrusion { curve: Curve, direction: Point },
    /// `curve(v)` turned by `u` about the axis through `origin`.
    Revolution {
        curve: Curve,
        origin: Point,
        axis: Point,
    },
    Nurbs {
        surface: NurbsSurface,
        closed: [bool; 2],
    },
}

impl Surface {
    pub fn nurbs(surface: NurbsSurface, tolerance: f64) -> Self {
        // Closure must hold along the WHOLE boundary, never at samples.
        let closed = [0, 1].map(|axis| {
            surface
                .boundary_gap(axis)
                .is_some_and(|gap| gap <= tolerance)
        });
        Self::Nurbs { surface, closed }
    }

    /// Parameter periods along u and v.
    pub fn periods(&self) -> [Option<f64>; 2] {
        match self {
            Self::Plane(_) => [None, None],
            Self::Cylinder { .. } | Self::Cone { .. } | Self::Sphere { .. } => [Some(TAU), None],
            Self::Torus { .. } => [Some(TAU), Some(TAU)],
            Self::Extrusion { curve, .. } => [curve.period(), None],
            Self::Revolution { curve, .. } => [Some(TAU), curve.period()],
            Self::Nurbs { surface, closed } => [0, 1].map(|axis| {
                let [a, b] = surface.range(axis);
                closed[axis].then_some(b - a)
            }),
        }
    }

    pub fn is_plane(&self) -> bool {
        matches!(self, Self::Plane(_))
    }

    /// Point with the partial derivatives along u and v.
    pub fn evaluate(&self, [u, v]: [f64; 2]) -> (Point, [Point; 2]) {
        match self {
            Self::Plane(f) => (add(f.origin, add(scale(f.x, u), scale(f.y, v))), [f.x, f.y]),
            Self::Cylinder { frame, radius } => {
                let (r, t) = frame.radial(u);
                (
                    add(frame.origin, add(scale(r, *radius), scale(frame.z, v))),
                    [scale(t, *radius), frame.z],
                )
            }
            Self::Cone {
                frame,
                radius,
                half_angle,
            } => {
                let (r, t) = frame.radial(u);
                let slope = half_angle.portable_tan();
                let rho = radius + v * slope;
                (
                    add(frame.origin, add(scale(r, rho), scale(frame.z, v))),
                    [scale(t, rho), add(scale(r, slope), frame.z)],
                )
            }
            Self::Sphere { frame, radius } => {
                let (r, t) = frame.radial(u);
                let (sin, cos) = v.portable_sin_cos();
                (
                    add(
                        frame.origin,
                        scale(add(scale(r, cos), scale(frame.z, sin)), *radius),
                    ),
                    [
                        scale(t, radius * cos),
                        scale(add(scale(r, -sin), scale(frame.z, cos)), *radius),
                    ],
                )
            }
            Self::Torus {
                frame,
                major,
                minor,
            } => {
                let (r, t) = frame.radial(u);
                let (sin, cos) = v.portable_sin_cos();
                let rho = major + minor * cos;
                (
                    add(
                        frame.origin,
                        add(scale(r, rho), scale(frame.z, minor * sin)),
                    ),
                    [
                        scale(t, rho),
                        scale(add(scale(r, -sin), scale(frame.z, cos)), *minor),
                    ],
                )
            }
            Self::Extrusion { curve, direction } => {
                let (p, d) = curve.evaluate(u);
                (add(p, scale(*direction, v)), [d, *direction])
            }
            Self::Revolution {
                curve,
                origin,
                axis,
            } => {
                let (p, d) = curve.evaluate(v);
                let (sin, cos) = u.portable_sin_cos();
                let turn = |q: Point| {
                    let along = scale(*axis, dot(q, *axis));
                    let radial = sub(q, along);
                    let tangent = cross(*axis, radial);
                    (
                        add(along, add(scale(radial, cos), scale(tangent, sin))),
                        add(scale(radial, -sin), scale(tangent, cos)),
                    )
                };
                let (position, around) = turn(sub(p, *origin));
                (add(*origin, position), [around, turn(d).0])
            }
            Self::Nurbs { surface, closed } => {
                let [u, v] = Self::confine(surface, *closed, [u, v]);
                surface.evaluate(u, v)
            }
        }
    }

    /// Keep spline parameters inside their knot range, wrapping closed ones.
    fn confine(surface: &NurbsSurface, closed: [bool; 2], uv: [f64; 2]) -> [f64; 2] {
        [0, 1].map(|axis| Self::confine_axis(surface, closed, axis, uv[axis]))
    }

    /// [`Self::confine`] along one axis; the axes are independent.
    pub(crate) fn confine_axis(
        surface: &NurbsSurface,
        closed: [bool; 2],
        axis: usize,
        t: f64,
    ) -> f64 {
        let [a, b] = surface.range(axis);
        if closed[axis] {
            a + (t - a).rem_euclid(b - a)
        } else {
            t.clamp(a, b)
        }
    }

    /// The parameter range of each axis carried by an open spline.
    pub(crate) fn bounded(&self) -> [Option<[f64; 2]>; 2] {
        let open = |curve: &Curve| match curve {
            Curve::Nurbs {
                curve,
                closed: false,
            } => Some(curve.range()),
            _ => None,
        };
        match self {
            Self::Nurbs { surface, closed } => {
                [0, 1].map(|axis| (!closed[axis]).then(|| surface.range(axis)))
            }
            Self::Extrusion { curve, .. } => [open(curve), None],
            Self::Revolution { curve, .. } => [None, open(curve)],
            _ => [None, None],
        }
    }

    /// True where the u derivative collapses, as at a sphere pole or apex.
    /// A spline whose parametrization merely stalls at one point, such as an
    /// end segment with a repeated control point, is no pole: a pole
    /// collapses the whole u row at its v, so row points apart by more than
    /// the tolerance prove there is none.
    pub fn is_pole(&self, uv: [f64; 2], tolerance: f64) -> bool {
        let (point, [su, _]) = self.evaluate(uv);
        if norm(su) > tolerance {
            return false;
        }
        let surface = match self {
            Self::Nurbs { surface, .. } => surface,
            // Every u isoline of an extrusion is the whole translated
            // curve: a stalled parametrization, as at a spline end that
            // repeats its control point, never collapses it to a point.
            Self::Extrusion { .. } => return false,
            _ => return true,
        };
        let [a, b] = surface.range(0);
        [a, 0.5 * (a + b), b]
            .into_iter()
            .all(|u| distance(self.evaluate([u, uv[1]]).0, point) <= tolerance)
    }

    /// Parameters of the surface point closest to `point`. A seed selects
    /// the periodic representative and the free coordinate at a pole.
    pub fn invert(&self, point: Point, seed: Option<[f64; 2]>) -> [f64; 2] {
        self.invert_prepared(point, seed, None)
    }

    fn invert_prepared(
        &self,
        point: Point,
        seed: Option<[f64; 2]>,
        prepared: Option<&PreparedInverse<'_>>,
    ) -> [f64; 2] {
        let free = seed.map_or(0., |seed| seed[0]);
        let angle = |frame: &Frame, flip: bool| {
            let l = frame.local(point);
            if l[0].portable_hypot(l[1]) <= 1e-12 * (1. + norm(l)) {
                free
            } else if flip {
                (-l[1]).portable_atan2(-l[0]).rem_euclid(TAU)
            } else {
                l[1].portable_atan2(l[0]).rem_euclid(TAU)
            }
        };
        let uv = match self {
            Self::Plane(frame) => {
                let l = frame.local(point);
                [l[0], l[1]]
            }
            Self::Cylinder { frame, .. } => [angle(frame, false), frame.local(point)[2]],
            Self::Cone {
                frame,
                radius,
                half_angle,
            } => {
                let v = frame.local(point)[2];
                [angle(frame, radius + v * half_angle.portable_tan() < 0.), v]
            }
            Self::Sphere { frame, .. } => {
                let l = frame.local(point);
                [
                    angle(frame, false),
                    l[2].portable_atan2(l[0].portable_hypot(l[1])),
                ]
            }
            Self::Torus { frame, major, .. } => {
                let l = frame.local(point);
                [
                    angle(frame, false),
                    l[2].portable_atan2(l[0].portable_hypot(l[1]) - major)
                        .rem_euclid(TAU),
                ]
            }
            _ => self.closest(point, seed, prepared),
        };
        let periods = self.periods();
        seed.map_or(uv, |seed| {
            [0, 1].map(|axis| nearest(uv[axis], seed[axis], periods[axis]))
        })
    }

    pub(crate) fn candidates(&self, point: Point) -> Vec<[f64; 2]> {
        match self {
            Self::Extrusion { curve, direction } => seeds(curve, point)
                .into_iter()
                .map(|u| {
                    let offset = sub(point, curve.evaluate(u).0);
                    [u, dot(offset, *direction) / dot(*direction, *direction)]
                })
                .collect(),
            Self::Revolution {
                curve,
                origin,
                axis,
            } => {
                let radial = |q: Point| sub(q, scale(*axis, dot(q, *axis)));
                let target = radial(sub(point, *origin));
                seeds(curve, point)
                    .into_iter()
                    .map(|v| {
                        let source = radial(sub(curve.evaluate(v).0, *origin));
                        let u =
                            dot(cross(source, target), *axis).portable_atan2(dot(source, target));
                        [u.rem_euclid(TAU), v]
                    })
                    .collect()
            }
            Self::Nurbs { surface, .. } => {
                let cuts = [0, 1].map(|axis| {
                    let [a, b] = surface.range(axis);
                    surface.breaks(axis, None, a, b)
                });
                let mut grid = Vec::new();
                for &u in &cuts[0] {
                    for &v in &cuts[1] {
                        grid.push([u, v]);
                    }
                }
                grid
            }
            _ => vec![],
        }
    }

    fn closest(
        &self,
        point: Point,
        seed: Option<[f64; 2]>,
        prepared: Option<&PreparedInverse<'_>>,
    ) -> [f64; 2] {
        let residual = |uv: [f64; 2]| distance(self.evaluate(uv).0, point);
        // Newton from the few candidates nearest the point: evaluating every
        // candidate is cheap, refining every one is not.
        let global = || {
            let mut ranked: Vec<(f64, [f64; 2])> = match prepared.and_then(|p| p.candidates()) {
                Some(candidates) => candidates
                    .iter()
                    .map(|Candidate { uv, position }| (distance(*position, point), *uv))
                    .collect(),
                None => self
                    .candidates(point)
                    .into_iter()
                    .map(|uv| (residual(uv), uv))
                    .collect(),
            };
            ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
            let mut best = (f64::INFINITY, [0.; 2]);
            for &(_, candidate) in ranked.iter().take(4) {
                let refined = self.newton(point, candidate);
                let d = residual(refined);
                if d < best.0 {
                    best = (d, refined);
                }
            }
            best
        };
        match seed {
            Some(seed) => {
                let local = self.newton(point, seed);
                let d = residual(local);
                // A seeded search that settles on a foot point near the
                // surface, where the offset meets both tangents squarely,
                // keeps it. One that lands elsewhere fell into another
                // basin, and a distant foot can be a saddle or the far side
                // of the surface: the global search decides.
                let (s, tangents) = self.evaluate(local);
                let offset = sub(point, s);
                let square = tangents
                    .iter()
                    .all(|t| dot(offset, *t).abs() <= 1e-6 * norm(offset) * norm(*t));
                // Both thresholds are geometric lengths plus the rounding of
                // the coordinates, so a translated copy decides alike.
                let rounding = 64. * f64::EPSILON * norm(point);
                let near = 1e-6 + rounding;
                if d <= 1e-9 + rounding || (square && d <= near) {
                    return local;
                }
                let (other, uv) = global();
                if other < d { uv } else { local }
            }
            None => global().1,
        }
    }

    /// The best iterate of a Newton descent. Where a derivative stalls, as
    /// at a spline end, a full step can overshoot along the flat direction;
    /// the search never returns a point farther than one it already found.
    fn newton(&self, point: Point, start: [f64; 2]) -> [f64; 2] {
        let mut uv = start;
        let mut best = (f64::INFINITY, start);
        for _ in 0..60 {
            let (s, [su, sv]) = self.evaluate(uv);
            let r = sub(point, s);
            let d = norm(r);
            if d < best.0 {
                best = (d, uv);
            }
            let (a, b, c) = (dot(su, su), dot(su, sv), dot(sv, sv));
            let (ru, rv) = (dot(r, su), dot(r, sv));
            let determinant = a * c - b * b;
            let step = if determinant > 1e-12 * a * c {
                [
                    (c * ru - b * rv) / determinant,
                    (a * rv - b * ru) / determinant,
                ]
            } else {
                [
                    if a > c * 1e-12 { ru / a } else { 0. },
                    if c > a * 1e-12 { rv / c } else { 0. },
                ]
            };
            // A bounded spline parameter stops at its end: beyond it the
            // spline repeats its end point and would read as a closer match.
            let bounded = self.bounded();
            let next = [0, 1].map(|axis| {
                let moved = uv[axis] + step[axis];
                bounded[axis].map_or(moved, |[a, b]| moved.clamp(a, b))
            });
            let moved = (next[0] - uv[0]).abs() + (next[1] - uv[1]).abs();
            uv = next;
            if moved <= 1e-15 * (1. + uv[0].abs() + uv[1].abs()) {
                break;
            }
        }
        if distance(self.evaluate(uv).0, point) < best.0 {
            uv
        } else {
            best.1
        }
    }

    /// The v parameter where the surface closes to a point above or below
    /// `from`. A horn torus touches its axis at v = pi, once per period.
    pub(crate) fn pole(&self, upper: bool, from: f64) -> Option<f64> {
        match self {
            Self::Sphere { .. } => Some(if upper { FRAC_PI_2 } else { -FRAC_PI_2 }),
            Self::Cone {
                radius, half_angle, ..
            } => {
                let slope = half_angle.portable_tan();
                ((slope < 0.) == upper).then_some(-radius / slope)
            }
            Self::Torus { major, minor, .. } if major == minor => {
                let turns = ((from - PI) / TAU).floor();
                Some(PI + TAU * (turns + f64::from(u8::from(upper))))
            }
            _ => None,
        }
    }

    /// The complete parameter rectangle of a closed surface.
    pub(crate) fn complete(&self) -> Option<[[f64; 2]; 2]> {
        match self {
            Self::Sphere { .. } => Some([[0., -FRAC_PI_2], [TAU, FRAC_PI_2]]),
            Self::Torus { .. } => Some([[0., 0.], [TAU, TAU]]),
            _ => None,
        }
    }

    /// Quadrature and sampling panels along one parameter axis.
    pub(crate) fn breaks(&self, axis: usize, a: f64, b: f64) -> Vec<f64> {
        match (self, axis) {
            (Self::Plane(_), _)
            | (Self::Cylinder { .. } | Self::Cone { .. } | Self::Extrusion { .. }, 1) => {
                vec![a, b]
            }
            (Self::Extrusion { curve, .. }, _) | (Self::Revolution { curve, .. }, 1) => {
                curve.breaks(a, b)
            }
            (Self::Nurbs { surface, .. }, _) => surface.breaks(axis, self.periods()[axis], a, b),
            _ => uniform(a, b, PANEL),
        }
    }

    pub(crate) fn mapped(
        &self,
        map: &dyn Fn(Point) -> Point,
        linear: &dyn Fn(Point) -> Point,
    ) -> Self {
        match self {
            Self::Plane(frame) => Self::Plane(frame.mapped(map, linear)),
            Self::Cylinder { frame, radius } => Self::Cylinder {
                frame: frame.mapped(map, linear),
                radius: *radius,
            },
            Self::Cone {
                frame,
                radius,
                half_angle,
            } => Self::Cone {
                frame: frame.mapped(map, linear),
                radius: *radius,
                half_angle: *half_angle,
            },
            Self::Sphere { frame, radius } => Self::Sphere {
                frame: frame.mapped(map, linear),
                radius: *radius,
            },
            Self::Torus {
                frame,
                major,
                minor,
            } => Self::Torus {
                frame: frame.mapped(map, linear),
                major: *major,
                minor: *minor,
            },
            Self::Extrusion { curve, direction } => Self::Extrusion {
                curve: curve.mapped(map, linear),
                direction: linear(*direction),
            },
            Self::Revolution {
                curve,
                origin,
                axis,
            } => Self::Revolution {
                curve: curve.mapped(map, linear),
                origin: map(*origin),
                axis: linear(*axis),
            },
            Self::Nurbs { surface, closed } => Self::Nurbs {
                surface: surface.mapped(map),
                closed: *closed,
            },
        }
    }
}

/// Starting parameters on a generator curve for a closest-point search.
fn seeds(curve: &Curve, point: Point) -> Vec<f64> {
    match curve {
        Curve::Line { .. } | Curve::Hyperbola { .. } | Curve::Parabola { .. } => {
            vec![curve.invert(point, None)]
        }
        Curve::Circle { .. } | Curve::Ellipse { .. } => uniform(0., TAU, PANEL),
        Curve::Nurbs { curve: spline, .. } => {
            let [a, b] = spline.range();
            curve.breaks(a, b)
        }
    }
}
