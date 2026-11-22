//! A round tube fitted to the outside of a crossing round member: the
//! fish-mouth (saddle) cope of a weldment branch.
//!
//! The tube runs along the unit axis `a` from the point `o` where its axis
//! meets the axis `b` of the receiving envelope cylinder of radius `Rc`.
//! With `e1` the unit component of `b` across `a`, `e2 = a x e1`,
//! `cos = a.b` and `sin = |b - cos a|`, every generator of the tube wall at
//! radius `rho` and angle `theta` leaves the envelope at
//!
//! `t0 = (rho cos cos(theta) + sqrt(Rc^2 - rho^2 sin(theta)^2)) / sin`,
//!
//! the larger root of `dist(o + rho e(theta) + t a, b)^2 = Rc^2`. Beyond it
//! the distance to the envelope axis only grows, so the kept tube is
//! `t0 <= t <= length` on every generator and the plane through the
//! receiving axis that bounds the kept side touches it at most at isolated
//! points. The boundary keeps its actual supports: both tube walls, the
//! square end plane and the envelope cylinder. Where `rho == Rc` the saddle
//! is two exact half ellipses; otherwise it is the quartic above, carried as
//! a C1 piecewise cubic Hermite curve through exact points and tangents
//! whose sampled deviation stays below a sixteenth of the tolerance. The
//! boundary validates and integrates itself and must reproduce the
//! closed-form volume.
use super::*;
use crate::geometry::authored_sweep::{AuthoredNaming, FaceName};
use brep::{Boundary, Edge, Face, Frame, Loop, NurbsCurve, Shell};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

type P = [f64; 3];

/// A round tube coped to a crossing envelope cylinder.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CopedTube {
    /// Where the tube axis meets the envelope axis.
    pub origin_mm: P,
    /// Direction from the cope toward the square far end.
    pub axis: P,
    pub outer_radius_mm: f64,
    pub inner_radius_mm: f64,
    /// Distance along `axis` from `origin_mm` to the square far end.
    pub length_mm: f64,
    /// Direction of the receiving envelope cylinder through `origin_mm`.
    pub envelope_axis: P,
    pub envelope_radius_mm: f64,
}

/// Names of the four coped-tube faces.
#[derive(Clone, Debug, PartialEq)]
pub struct CopedTubeNames {
    pub outer: Origin,
    pub inner: Origin,
    pub end: Origin,
    pub saddle: Origin,
}

/// The material of a body lies outside this infinite cylinder.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ExcludedCylinder {
    pub(crate) point_mm: P,
    pub(crate) axis: P,
    pub(crate) radius_mm: f64,
}

const SAMPLES: usize = 8;
const MAX_SEGMENTS: usize = 1 << 12;
/// Trapezoid nodes for the periodic closed-form volume.
const VOLUME_NODES: usize = 1 << 14;

fn dot(a: P, b: P) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: P, b: P) -> P {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn mix(a: P, s: f64, b: P) -> P {
    std::array::from_fn(|i| a[i] + s * b[i])
}
fn scale(a: P, s: f64) -> P {
    a.map(|x| x * s)
}
fn norm(a: P) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: P) -> Option<P> {
    let n = norm(a);
    (n.is_finite() && n > 0.).then(|| scale(a, 1. / n))
}

fn rejected(message: &str) -> GeometryError {
    error("cope", GeometryErrorKind::Unsupported, message)
}

/// The tube chart and the saddle law.
struct Chart {
    o: P,
    a: P,
    e1: P,
    e2: P,
    b: P,
    cos: f64,
    sin: f64,
    rc: f64,
}

impl Chart {
    fn radial(&self, theta: f64) -> (P, P) {
        let (s, c) = theta.sin_cos();
        (
            mix(scale(self.e1, c), s, self.e2),
            mix(scale(self.e1, -s), c, self.e2),
        )
    }

    /// Saddle point and its theta derivative at wall radius `rho`.
    fn saddle(&self, rho: f64, theta: f64) -> (P, P) {
        let (s, c) = theta.sin_cos();
        let root = (self.rc * self.rc - rho * rho * s * s).max(0.).sqrt();
        let t = (rho * self.cos * c + root) / self.sin;
        let dt = (-rho * self.cos * s - rho * rho * s * c / root) / self.sin;
        let (r, tangent) = self.radial(theta);
        (
            mix(mix(self.o, rho, r), t, self.a),
            mix(scale(tangent, rho), dt, self.a),
        )
    }

    /// Deviation of a point from the tube wall at `rho` and from the
    /// envelope cylinder.
    fn residual(&self, rho: f64, point: P) -> f64 {
        let d: P = std::array::from_fn(|i| point[i] - self.o[i]);
        let along = dot(d, self.a);
        let wall = (norm(mix(d, -along, self.a)) - rho).abs();
        let across = dot(d, self.b);
        let envelope = (norm(mix(d, -across, self.b)) - self.rc).abs();
        wall.max(envelope)
    }

    /// The saddle at radius `rho < Rc` over `[from, from + PI]` as C1 cubic
    /// Bezier pieces through exact points and tangents.
    fn hermite(&self, rho: f64, from: f64, tolerance: f64) -> Result<NurbsCurve, GeometryError> {
        let mut segments = 16;
        loop {
            let h = PI / segments as f64;
            let nodes: Vec<_> = (0..=segments)
                .map(|k| self.saddle(rho, from + h * k as f64))
                .collect();
            let mut points = Vec::with_capacity(3 * segments + 1);
            let mut worst: f64 = 0.;
            for k in 0..segments {
                let [(p0, d0), (p1, d1)] = [nodes[k], nodes[k + 1]];
                let control = [p0, mix(p0, h / 3., d0), mix(p1, -h / 3., d1), p1];
                if k == 0 {
                    points.push(p0);
                }
                points.extend_from_slice(&control[1..]);
                for j in 1..SAMPLES {
                    let u = j as f64 / SAMPLES as f64;
                    let w = [
                        (1. - u).powi(3),
                        3. * u * (1. - u).powi(2),
                        3. * u * u * (1. - u),
                        u.powi(3),
                    ];
                    let bezier: P =
                        std::array::from_fn(|i| (0..4).map(|n| w[n] * control[n][i]).sum());
                    let exact = self.saddle(rho, from + h * (k as f64 + u)).0;
                    let gap = norm(std::array::from_fn(|i| bezier[i] - exact[i]));
                    worst = worst.max(gap).max(self.residual(rho, bezier));
                }
            }
            if worst.is_finite() && worst <= tolerance / 16. {
                let distinct = (0..=segments).map(|k| from + h * k as f64).collect();
                let mut multiplicities = vec![3; segments + 1];
                multiplicities[0] = 4;
                multiplicities[segments] = 4;
                return NurbsCurve::new(3, points, None, distinct, multiplicities).map_err(failure);
            }
            segments *= 2;
            if segments > MAX_SEGMENTS {
                return Err(rejected(
                    "The cope saddle is too close to tangency for an exact Rust boundary",
                ));
            }
        }
    }

    /// Closed-form kept volume `pi T (R^2 - r^2) - (1/sin) int int t0 rho`.
    fn volume(&self, [inner, outer]: [f64; 2], length: f64) -> f64 {
        let rc2 = self.rc * self.rc;
        let band = |nodes: usize| {
            (0..nodes)
                .map(|k| {
                    let s = (TAU * k as f64 / nodes as f64).sin();
                    let a = rc2 - inner * inner * s * s;
                    let b = (rc2 - outer * outer * s * s).max(0.);
                    (outer * outer - inner * inner) * (a * a + a * b + b * b)
                        / (3. * (a.powf(1.5) + b.powf(1.5)))
                })
                .sum::<f64>()
                * TAU
                / nodes as f64
        };
        PI * length * (outer * outer - inner * inner) - band(VOLUME_NODES) / self.sin
    }
}

impl Solid {
    /// The exact retained boundary of a round tube coped to a crossing
    /// envelope cylinder; see the module notes.
    pub fn rust_coped_round_tube(
        tube: &CopedTube,
        names: &CopedTubeNames,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        let tolerance = policy.linear_tolerance_mm;
        let CopedTube {
            origin_mm: o,
            outer_radius_mm: outer,
            inner_radius_mm: inner,
            length_mm: length,
            envelope_radius_mm: rc,
            ..
        } = *tube;
        policy.point(o)?;
        for value in [inner, outer, length, rc, outer - inner] {
            policy.length(value)?;
        }
        if feature.is_nil() || outer > rc {
            return Err(rejected(
                "A Rust cope needs a round branch no wider than the receiving envelope",
            ));
        }
        let (Some(a), Some(b)) = (unit(tube.axis), unit(tube.envelope_axis)) else {
            return Err(input("A cope needs two axis directions"));
        };
        let cos = dot(a, b);
        let across = mix(b, -cos, a);
        let sin = norm(across);
        if sin.is_nan() || sin < 0.1 {
            return Err(rejected("The cope axes are too close to parallel"));
        }
        let e1 = scale(across, 1. / sin);
        let chart = Chart {
            o,
            a,
            e1,
            e2: cross(a, e1),
            b,
            cos,
            sin,
            rc,
        };
        // The saddle may not reach the square far end.
        let reach = (outer * cos.abs() + rc) / sin;
        policy.length(length - reach).map_err(|_| {
            rejected("The coped branch is too short to keep its square end clear of the saddle")
        })?;
        let mut vertices = Vec::with_capacity(6);
        let mut vertex = |p: P| {
            vertices.push(p);
            vertices.len() - 1
        };
        let rings = [outer, inner];
        let side: Vec<[usize; 2]> = rings
            .iter()
            .map(|&rho| [-FRAC_PI_2, FRAC_PI_2].map(|theta| vertex(chart.saddle(rho, theta).0)))
            .collect();
        let end = mix(o, length, a);
        let tops: Vec<usize> = rings.iter().map(|&rho| vertex(mix(end, rho, e1))).collect();
        let mut edges = Vec::with_capacity(6);
        let mut add = |curve: Curve, ends: [usize; 2], vertices: &[P]| {
            let points = ends.map(|v| vertices[v]);
            edges.push(Edge::new(curve, ends, points, tolerance).map_err(failure)?);
            Ok::<_, GeometryError>(edges.len() - 1)
        };
        // Saddle halves: `cos(theta) >= 0` from -90 to 90 degrees, then back.
        let mut saddles = Vec::with_capacity(2);
        for (ring, &rho) in rings.iter().enumerate() {
            let [low, high] = side[ring];
            let halves = if rho == rc {
                // Both cylinders have the same radius: the saddle lies in
                // the two planes t = rho cos(theta) (cos +- 1) / sin.
                [1., -1.].map(|sign| {
                    let k = (cos + sign) / sin;
                    let x = mix(e1, k, a);
                    let major = norm(x);
                    let x = scale(x, 1. / major);
                    Curve::Ellipse {
                        frame: Frame {
                            origin: o,
                            x,
                            y: chart.e2,
                            z: cross(x, chart.e2),
                        },
                        semiaxes: [rho * major, rho],
                    }
                })
            } else {
                [-FRAC_PI_2, FRAC_PI_2]
                    .map(|from| {
                        chart
                            .hermite(rho, from, tolerance)
                            .map(|c| Curve::nurbs(c, tolerance))
                    })
                    .into_iter()
                    .collect::<Result<Vec<_>, _>>()?
                    .try_into()
                    .expect("two halves")
            };
            let [first, second] = halves;
            saddles.push([
                add(first, [low, high], &vertices)?,
                add(second, [high, low], &vertices)?,
            ]);
        }
        let top_frame = Frame {
            origin: end,
            x: e1,
            y: chart.e2,
            z: a,
        };
        let circles: Vec<usize> = rings
            .iter()
            .zip(&tops)
            .map(|(&radius, &v)| {
                add(
                    Curve::Circle {
                        frame: top_frame,
                        radius,
                    },
                    [v, v],
                    &vertices,
                )
            })
            .collect::<Result<_, _>>()?;
        let wall_frame = Frame {
            origin: o,
            x: e1,
            y: chart.e2,
            z: a,
        };
        let w = scale(mix(a, -cos, b), 1. / sin);
        let x = scale(w, -1.);
        let envelope_frame = Frame {
            origin: o,
            x,
            y: cross(b, x),
            z: b,
        };
        let [[o_low, o_high], [i_low, i_high]] = [saddles[0], saddles[1]];
        let faces = vec![
            Face {
                surface: Surface::Cylinder {
                    frame: wall_frame,
                    radius: outer,
                },
                same_sense: true,
                loops: vec![
                    Loop {
                        uses: vec![(o_low, true), (o_high, true)],
                        vertex: None,
                    },
                    Loop {
                        uses: vec![(circles[0], false)],
                        vertex: None,
                    },
                ],
            },
            Face {
                surface: Surface::Cylinder {
                    frame: wall_frame,
                    radius: inner,
                },
                same_sense: false,
                loops: vec![
                    Loop {
                        uses: vec![(i_high, false), (i_low, false)],
                        vertex: None,
                    },
                    Loop {
                        uses: vec![(circles[1], true)],
                        vertex: None,
                    },
                ],
            },
            Face {
                surface: Surface::Plane(top_frame),
                same_sense: true,
                loops: vec![
                    Loop {
                        uses: vec![(circles[0], true)],
                        vertex: None,
                    },
                    Loop {
                        uses: vec![(circles[1], false)],
                        vertex: None,
                    },
                ],
            },
            Face {
                surface: Surface::Cylinder {
                    frame: envelope_frame,
                    radius: rc,
                },
                same_sense: false,
                loops: vec![
                    Loop {
                        uses: vec![(o_high, false), (o_low, false)],
                        vertex: None,
                    },
                    Loop {
                        uses: vec![(i_low, true), (i_high, true)],
                        vertex: None,
                    },
                ],
            },
        ];
        let boundary = Boundary {
            vertices,
            edges,
            shells: vec![Shell {
                faces: (0..faces.len()).collect(),
                reversed: false,
            }],
            faces,
        };
        let solid = brep::Solid::new(boundary, tolerance, None).map_err(failure)?;
        let expected = chart.volume([inner, outer], length);
        if (solid.volume_mm3() - expected).abs() > solid.volume_error_mm3().max(1e-9 * expected) {
            return Err(error(
                "cope",
                GeometryErrorKind::InvalidResult,
                "Rust coped tube boundary does not reproduce its closed-form volume",
            ));
        }
        let plain = |origin: &Origin| FaceName {
            origins: vec![origin.clone()],
            support: None,
            anchor: None,
        };
        let cap = ReferenceFrame {
            origin_mm: end,
            u: e1,
            v: chart.e2,
        };
        let naming = AuthoredNaming {
            faces: vec![
                plain(&names.outer),
                plain(&names.inner),
                FaceName {
                    origins: vec![names.end.clone()],
                    support: Some(cap),
                    anchor: Some(cap),
                },
                plain(&names.saddle),
            ],
        };
        let mut result = Self::from_retained_brep(solid, feature, Rc::new(naming), policy)?;
        if let SolidBody::ImportedBrep(body) = &mut result.body {
            Rc::get_mut(body).expect("a fresh retained body").excluded = Some(ExcludedCylinder {
                point_mm: o,
                axis: b,
                radius_mm: rc,
            });
        }
        Ok(result)
    }

    /// A coped body meets the member it is fitted to only on the envelope
    /// surface: their common part has no volume.
    pub(in crate::geometry) fn cope_contact_intersection(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Intersect
            || policy.boolean_fuzzy_mm != 0.
            || policy != self.policy
            || policy != other.policy
        {
            return None;
        }
        let excluded = |solid: &Self| match &solid.body {
            SolidBody::ImportedBrep(body) => body.excluded,
            _ => None,
        };
        let (cylinder, inside) = match (excluded(self), excluded(other)) {
            (Some(cylinder), _) => (cylinder, other),
            (None, Some(cylinder)) => (cylinder, self),
            (None, None) => return None,
        };
        let (frame, radius) = if let Some((annulus, frame)) = inside.annulus_body() {
            (frame, annulus.geometry.outer_radius_mm())
        } else {
            return None;
        };
        let axis = frame.normal();
        let d: P = std::array::from_fn(|i| frame.origin[i] - cylinder.point_mm[i]);
        let offset = norm(mix(d, -dot(d, cylinder.axis), cylinder.axis));
        let magnitude = frame
            .origin
            .into_iter()
            .chain(cylinder.point_mm)
            .map(f64::abs)
            .fold(radius.max(cylinder.radius_mm), f64::max);
        let rounding = 64. * f64::EPSILON * magnitude;
        let parallel = norm(cross(axis, cylinder.axis)) <= 64. * f64::EPSILON;
        // The member lies within `radius + offset` of the envelope axis.
        (parallel && radius + offset <= cylinder.radius_mm + rounding).then_some(Ok(None))
    }
}

#[cfg(test)]
mod tests;
