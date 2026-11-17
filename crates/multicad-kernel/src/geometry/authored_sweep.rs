//! Exact authored boundaries for curved sweeps that no specialised body owns:
//! oblique and plane-terminated translations of a certified section, and
//! straight rulings between a section and its tapered offset.
//!
//! Every face keeps its source curve or the exact affine image of that curve
//! on its cap plane: lines stay lines, circles become circles or ellipses and
//! splines keep their knots and weights with mapped poles. Side walls are
//! planes, linear extrusions of the source curve, right cones or exact
//! rational rulings. The retained boundary validates its own closure and
//! integrates its own measurements; no display mesh is involved.
use super::*;
use crate::{CurvedEdge2d, CurvedProfile};
use multicad_core::references::ReferenceFrame;
use multicad_rust_kernel::brep::{
    self, Boundary, Curve, Frame, Loop, NurbsCurve, NurbsSurface, Shell, Surface,
};
use multicad_rust_kernel::spline_sweep::Cubic;
use std::f64::consts::{FRAC_1_SQRT_2, TAU};

mod glue;

const OPERATION: &str = "Rust authored sweep";
type P = [f64; 3];
type P2 = [f64; 2];

/// Authored identity of one retained face, in face order.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct FaceName {
    pub(crate) origins: Vec<Origin>,
    /// The planar support with its outward normal.
    pub(crate) support: Option<ReferenceFrame>,
    /// The authored attachment frame, which need not face outward.
    pub(crate) anchor: Option<ReferenceFrame>,
}

/// Authored names for every face of a retained boundary. Edges and vertices
/// take the names of the faces they bound.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AuthoredNaming {
    pub(crate) faces: Vec<FaceName>,
}

impl AuthoredNaming {
    pub(super) fn moved(&self, transform: RigidTransform) -> Self {
        let frame = |frame: ReferenceFrame| {
            let origin = transform.vector(frame.origin_mm);
            ReferenceFrame {
                origin_mm: std::array::from_fn(|i| origin[i] + transform.translation_mm[i]),
                u: transform.vector(frame.u),
                v: transform.vector(frame.v),
            }
        };
        Self {
            faces: self
                .faces
                .iter()
                .map(|face| FaceName {
                    origins: face.origins.clone(),
                    support: face.support.map(frame),
                    anchor: face.anchor.map(frame),
                })
                .collect(),
        }
    }
}

pub(super) fn failure(cause: brep::Error) -> GeometryError {
    let (kind, message) = match cause {
        brep::Error::Invalid => (
            GeometryErrorKind::InvalidResult,
            "Authored sweep boundary does not close into one valid solid",
        ),
        brep::Error::Unsupported => (
            GeometryErrorKind::Unsupported,
            "Authored sweep boundary is outside the exact boundary provider",
        ),
        brep::Error::Limit => (
            GeometryErrorKind::Limit,
            "Authored sweep boundary exceeds the Rust budget",
        ),
        brep::Error::Cancelled => (
            GeometryErrorKind::Cancelled,
            "Authored sweep boundary superseded",
        ),
    };
    error(OPERATION, kind, message)
}

fn add(a: P, b: P) -> P {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: P, b: P) -> P {
    std::array::from_fn(|i| a[i] - b[i])
}
fn scale(a: P, s: f64) -> P {
    a.map(|x| x * s)
}
fn norm(a: P) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: P) -> Result<P, GeometryError> {
    let length = norm(a);
    if !length.is_finite() || length <= 0. {
        return Err(input("Authored sweep has a degenerate direction"));
    }
    Ok(scale(a, 1. / length))
}

/// `point -> m * point + t`.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Affine {
    m: [P; 3],
    t: P,
}

impl Affine {
    fn translation(t: P) -> Self {
        Self {
            m: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            t,
        }
    }

    /// Parallel projection along `direction` onto `plane`.
    fn projection(direction: P, plane: PlaneFrame) -> Self {
        let n = plane.normal();
        let crossing = dot(n, direction);
        let m = std::array::from_fn(|i| {
            std::array::from_fn(|j| f64::from(u8::from(i == j)) - direction[i] * n[j] / crossing)
        });
        Self {
            m,
            t: scale(direction, dot(n, plane.origin) / crossing),
        }
    }

    fn vector(self, v: P) -> P {
        self.m.map(|row| dot(row, v))
    }

    fn point(self, p: P) -> P {
        add(self.vector(p), self.t)
    }
}

/// One profile curve in world coordinates, in its authored direction.
#[derive(Clone, Debug, PartialEq)]
enum Path {
    Line([P; 2]),
    /// `center + axes[0] cos t + axes[1] sin t` for `t` in `range`.
    Conic {
        center: P,
        axes: [P; 2],
        range: [f64; 2],
    },
    Spline {
        degree: usize,
        poles: Vec<P>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<usize>,
    },
}

impl Path {
    fn sketch(frame: PlaneFrame, edge: &CurvedEdge2d) -> Self {
        let point = |p: P2| frame.point(p);
        match edge {
            CurvedEdge2d::Line { a, b } => Self::Line([point(*a), point(*b)]),
            CurvedEdge2d::Circle { center, radius_mm } => Self::Conic {
                center: point(*center),
                axes: [scale(frame.u, *radius_mm), scale(frame.v, *radius_mm)],
                range: [0., TAU],
            },
            CurvedEdge2d::Arc {
                center,
                a,
                b,
                clockwise,
            } => {
                let [x, y] = [a[0] - center[0], a[1] - center[1]];
                let radius = x.hypot(y);
                let start = y.atan2(x);
                let end = (b[1] - center[1]).atan2(b[0] - center[0]);
                let (v, start, sweep) = if *clockwise {
                    (
                        scale(frame.v, -radius),
                        -start,
                        (start - end).rem_euclid(TAU),
                    )
                } else {
                    (scale(frame.v, radius), start, (end - start).rem_euclid(TAU))
                };
                Self::Conic {
                    center: point(*center),
                    axes: [scale(frame.u, radius), v],
                    range: [start, start + sweep],
                }
            }
            CurvedEdge2d::BSpline(spline) => Self::Spline {
                degree: spline.degree as usize,
                poles: spline.poles.iter().map(|&p| point(p)).collect(),
                weights: spline.weights.clone(),
                knots: spline.knots.clone(),
                multiplicities: spline.multiplicities.iter().map(|&m| m as usize).collect(),
            },
        }
    }

    fn mapped(&self, map: Affine) -> Self {
        match self {
            Self::Line(ends) => Self::Line(ends.map(|p| map.point(p))),
            Self::Conic {
                center,
                axes,
                range,
            } => Self::Conic {
                center: map.point(*center),
                axes: axes.map(|axis| map.vector(axis)),
                range: *range,
            },
            Self::Spline {
                degree,
                poles,
                weights,
                knots,
                multiplicities,
            } => Self::Spline {
                degree: *degree,
                poles: poles.iter().map(|&p| map.point(p)).collect(),
                weights: weights.clone(),
                knots: knots.clone(),
                multiplicities: multiplicities.clone(),
            },
        }
    }

    fn ends(&self, tolerance: f64) -> Result<[P; 2], GeometryError> {
        Ok(match self {
            Self::Line(ends) => *ends,
            Self::Conic {
                center,
                axes,
                range,
            } => range.map(|t| {
                let (sin, cos) = t.sin_cos();
                add(*center, add(scale(axes[0], cos), scale(axes[1], sin)))
            }),
            Self::Spline { .. } => {
                let Curve::Nurbs { curve, .. } = self.curve(tolerance)? else {
                    unreachable!("spline path");
                };
                curve.range().map(|t| curve.evaluate(t).0)
            }
        })
    }

    fn curve(&self, tolerance: f64) -> Result<Curve, GeometryError> {
        match self {
            Self::Line([a, b]) => Ok(Curve::Line {
                origin: *a,
                direction: sub(*b, *a),
            }),
            Self::Conic { center, axes, .. } => conic(*center, *axes),
            Self::Spline {
                degree,
                poles,
                weights,
                knots,
                multiplicities,
            } => {
                let curve = NurbsCurve::new(
                    *degree,
                    poles.clone(),
                    weights.clone(),
                    knots.clone(),
                    multiplicities.clone(),
                )
                .map_err(failure)?;
                Ok(Curve::nurbs(curve, tolerance))
            }
        }
    }
}

/// A circle or an ellipse in its principal frame. The parameter keeps the
/// direction of `axes[0]` toward `axes[1]`.
fn conic(center: P, [a, b]: [P; 2]) -> Result<Curve, GeometryError> {
    let (aa, bb, ab) = (dot(a, a), dot(b, b), dot(a, b));
    let z = unit(cross(a, b))?;
    let scale_sq = aa.max(bb);
    if (aa - bb).abs() <= 1e-14 * scale_sq && ab.abs() <= 1e-14 * scale_sq {
        let x = unit(a)?;
        let radius = 0.5 * (aa.sqrt() + bb.sqrt());
        return Ok(Curve::Circle {
            frame: Frame {
                origin: center,
                x,
                y: cross(z, x),
                z,
            },
            radius,
        });
    }
    let phase = 0.5 * (2. * ab).atan2(aa - bb);
    let (sin, cos) = phase.sin_cos();
    let major = add(scale(a, cos), scale(b, sin));
    let minor = add(scale(a, -sin), scale(b, cos));
    let x = unit(major)?;
    Ok(Curve::Ellipse {
        frame: Frame {
            origin: center,
            x,
            y: cross(z, x),
            z,
        },
        semiaxes: [norm(major), norm(minor)],
    })
}

/// Twice the signed area enclosed with the origin, for orientation only.
fn twice_signed_area(edge: &CurvedEdge2d, policy: GeometryPolicy) -> Result<f64, GeometryError> {
    let cross2 = |a: P2, b: P2| a[0] * b[1] - a[1] * b[0];
    Ok(match edge {
        CurvedEdge2d::Line { a, b } => cross2(*a, *b),
        CurvedEdge2d::Circle { radius_mm, .. } => TAU * radius_mm * radius_mm,
        CurvedEdge2d::Arc {
            center,
            a,
            b,
            clockwise,
        } => {
            let [x, y] = [a[0] - center[0], a[1] - center[1]];
            let radius = x.hypot(y);
            let start = y.atan2(x);
            let end = (b[1] - center[1]).atan2(b[0] - center[0]);
            let sweep = if *clockwise {
                -(start - end).rem_euclid(TAU)
            } else {
                (end - start).rem_euclid(TAU)
            };
            let finish = start + sweep;
            radius * center[0] * (finish.sin() - start.sin())
                - radius * center[1] * (finish.cos() - start.cos())
                + radius * radius * sweep
        }
        CurvedEdge2d::BSpline(spline) => {
            // Polynomial spans integrate exactly with eight Gauss nodes; a
            // rational span is resolved finely enough for the sign alone.
            let curve =
                Path::sketch(PlaneFrame::default(), edge).curve(policy.linear_tolerance_mm)?;
            let Curve::Nurbs { curve, .. } = curve else {
                unreachable!("spline curve");
            };
            let mut total = 0.;
            for pair in spline.knots.windows(2) {
                let pieces = 8;
                let width = (pair[1] - pair[0]) / f64::from(pieces);
                for piece in 0..pieces {
                    let middle = pair[0] + width * (f64::from(piece) + 0.5);
                    for (node, weight) in GAUSS {
                        let (point, tangent) = curve.evaluate(middle + 0.5 * width * node);
                        total +=
                            0.5 * width * weight * (point[0] * tangent[1] - point[1] * tangent[0]);
                    }
                }
            }
            total
        }
    })
}

const GAUSS: [(f64, f64); 8] = [
    (-0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    (-0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (-0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (-0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];

/// Least and greatest `g . p` over one closed authored wire, from the true
/// curves rather than their control points.
fn linear_range(
    wire: &[CurvedEdge2d],
    g: P2,
    policy: GeometryPolicy,
) -> Result<[f64; 2], GeometryError> {
    let value = |p: P2| g[0] * p[0] + g[1] * p[1];
    let mut range = [f64::INFINITY, f64::NEG_INFINITY];
    let mut include = |x: f64| {
        range[0] = range[0].min(x);
        range[1] = range[1].max(x);
    };
    for edge in wire {
        match edge {
            CurvedEdge2d::Line { a, b } => {
                include(value(*a));
                include(value(*b));
            }
            CurvedEdge2d::Circle { center, radius_mm } => {
                let reach = radius_mm * g[0].hypot(g[1]);
                include(value(*center) - reach);
                include(value(*center) + reach);
            }
            CurvedEdge2d::Arc { center, a, b, .. } => {
                include(value(*a));
                include(value(*b));
                let Path::Conic {
                    axes, range: span, ..
                } = Path::sketch(PlaneFrame::default(), edge)
                else {
                    unreachable!("arc path");
                };
                let along = [
                    dot([g[0], g[1], 0.], axes[0]),
                    dot([g[0], g[1], 0.], axes[1]),
                ];
                for extremum in [along[1].atan2(along[0]), (-along[1]).atan2(-along[0])] {
                    let travel = (extremum - span[0]).rem_euclid(TAU);
                    if travel <= span[1] - span[0] {
                        let (sin, cos) = extremum.sin_cos();
                        include(value(*center) + along[0] * cos + along[1] * sin);
                    }
                }
            }
            CurvedEdge2d::BSpline(spline) => {
                for patch in
                    crate::spline_edits::regions::polynomial_patches(spline, 0, policy, None)?
                {
                    let cubic =
                        Cubic::new(patch.poles.to_vec(), vec![0., 1.], false).map_err(|cause| {
                            error(
                                OPERATION,
                                GeometryErrorKind::Unsupported,
                                format!("Cannot bound an authored spline span: {cause:?}"),
                            )
                        })?;
                    include(cubic.support(g));
                    include(-cubic.support([-g[0], -g[1]]));
                }
            }
        }
    }
    Ok(range)
}

/// How one side wall rules its start curve to its end curve.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Wall {
    /// Translation along one direction.
    Translation(P),
    /// Straight rulings between corresponding points of two curves.
    Ruling,
}

/// A closed wire in traversal order.
struct Ring {
    first: usize,
    len: usize,
    /// Whether the traversal follows each authored edge direction.
    forward: bool,
    outer: bool,
}

impl Ring {
    fn edge(&self, k: usize) -> usize {
        self.first
            + if self.forward {
                k % self.len
            } else {
                self.len - 1 - k % self.len
            }
    }
}

/// The authored attachment frames of the start and end caps.
struct Caps {
    anchors: [PlaneFrame; 2],
}

struct Sweep<'a> {
    profile: &'a CurvedProfile,
    start: Vec<Path>,
    end: Vec<Path>,
    wall: Wall,
    /// A direction of every generator from start to end.
    up: P,
    caps: Caps,
}

impl Sweep<'_> {
    fn rings(&self, policy: GeometryPolicy) -> Result<Vec<Ring>, GeometryError> {
        let sign = dot(self.up, self.profile.frame.normal()).signum();
        let mut first = 0;
        std::iter::once(&self.profile.outer)
            .chain(&self.profile.holes)
            .enumerate()
            .map(|(index, wire)| {
                let area = wire
                    .iter()
                    .map(|edge| twice_signed_area(edge, policy))
                    .sum::<Result<f64, _>>()?;
                let ring = Ring {
                    first,
                    len: wire.len(),
                    forward: (area * sign > 0.) == (index == 0),
                    outer: index == 0,
                };
                first += wire.len();
                Ok(ring)
            })
            .collect()
    }

    fn build(
        self,
        feature: EntityId,
        origins: &[Origin],
        policy: GeometryPolicy,
    ) -> Result<Solid, GeometryError> {
        let tolerance = policy.linear_tolerance_mm;
        let rings = self.rings(policy)?;
        let count = self.start.len();
        let mut boundary = Boundary::default();
        // Vertex pairs in traversal order, and the authored start vertex of
        // every edge on each cap.
        let mut low = vec![0; count];
        let mut high = vec![0; count];
        let mut ring_vertices = Vec::with_capacity(rings.len());
        for ring in &rings {
            let mut vertices = Vec::with_capacity(ring.len);
            for k in 0..ring.len {
                let edge = ring.edge(k);
                let side = usize::from(!ring.forward);
                let a = self.start[edge].ends(tolerance)?[side];
                let b = self.end[edge].ends(tolerance)?[side];
                boundary.vertices.extend([a, b]);
                let pair = [boundary.vertices.len() - 2, boundary.vertices.len() - 1];
                vertices.push(pair);
                if ring.forward {
                    low[edge] = pair[0];
                    high[edge] = pair[1];
                }
            }
            if !ring.forward {
                // The authored start of edge e is the traversal end.
                for k in 0..ring.len {
                    let edge = ring.edge(k);
                    let [a, b] = vertices[(k + 1) % ring.len];
                    low[edge] = a;
                    high[edge] = b;
                }
            }
            ring_vertices.push(vertices);
        }
        let mut cap_edges = [vec![0; count], vec![0; count]];
        let mut rails = Vec::with_capacity(rings.len());
        for (ring, vertices) in rings.iter().zip(&ring_vertices) {
            for k in 0..ring.len {
                let edge = ring.edge(k);
                let next = vertices[(k + 1) % ring.len];
                let here = vertices[k];
                for (cap, paths) in [&self.start, &self.end].into_iter().enumerate() {
                    let path = &paths[edge];
                    let [from, to] = if ring.forward {
                        [here[cap], next[cap]]
                    } else {
                        [next[cap], here[cap]]
                    };
                    let points = [from, to].map(|v| boundary.vertices[v]);
                    boundary.edges.push(
                        brep::Edge::new(path.curve(tolerance)?, [from, to], points, tolerance)
                            .map_err(failure)?,
                    );
                    cap_edges[cap][edge] = boundary.edges.len() - 1;
                }
            }
            let mut ring_rails = Vec::with_capacity(ring.len);
            for &[a, b] in vertices {
                let [pa, pb] = [a, b].map(|v| boundary.vertices[v]);
                boundary.edges.push(
                    brep::Edge::new(
                        Curve::Line {
                            origin: pa,
                            direction: sub(pb, pa),
                        },
                        [a, b],
                        [pa, pb],
                        tolerance,
                    )
                    .map_err(failure)?,
                );
                ring_rails.push(boundary.edges.len() - 1);
            }
            rails.push(ring_rails);
        }

        let mut names = Vec::with_capacity(count + 2);
        // Caps: outward supports, loops against the traversal at the start.
        for (cap, end) in [false, true].into_iter().enumerate() {
            let anchor = self.caps.anchors[cap];
            let normal = anchor.normal();
            let outward = (dot(normal, self.up) > 0.) == end;
            let support = ReferenceFrame {
                origin_mm: anchor.origin,
                u: anchor.u,
                v: if outward {
                    anchor.v
                } else {
                    anchor.v.map(|x| -x)
                },
            };
            let loops = rings
                .iter()
                .map(|ring| {
                    let uses = (0..ring.len).map(|k| {
                        let edge = ring.edge(k);
                        (cap_edges[cap][edge], ring.forward == end)
                    });
                    Loop {
                        uses: if end {
                            uses.collect()
                        } else {
                            uses.rev().collect()
                        },
                        vertex: None,
                    }
                })
                .collect();
            boundary.faces.push(brep::Face {
                surface: Surface::Plane(Frame {
                    origin: support.origin_mm,
                    x: support.u,
                    y: support.v,
                    z: support.normal(),
                }),
                same_sense: true,
                loops,
            });
            names.push(FaceName {
                origins: vec![Origin {
                    source: feature,
                    role: if end { "cap.end" } else { "cap.start" }.into(),
                }],
                support: Some(support),
                anchor: Some(ReferenceFrame {
                    origin_mm: anchor.origin,
                    u: anchor.u,
                    v: anchor.v,
                }),
            });
        }
        let mut sides = vec![None; count];
        for (ring, (vertices, ring_rails)) in rings.iter().zip(ring_vertices.iter().zip(&rails)) {
            for k in 0..ring.len {
                let edge = ring.edge(k);
                let next = (k + 1) % ring.len;
                let loop_uses = vec![
                    (cap_edges[0][edge], ring.forward),
                    (ring_rails[next], true),
                    (cap_edges[1][edge], !ring.forward),
                    (ring_rails[k], false),
                ];
                let corners = [vertices[k][0], vertices[next][0], vertices[k][1]]
                    .map(|v| boundary.vertices[v]);
                let (surface, same_sense, support, anchor) =
                    self.wall(edge, ring, corners, [low[edge], high[edge]], &boundary)?;
                sides[edge] = Some((
                    brep::Face {
                        surface,
                        same_sense,
                        loops: vec![Loop {
                            uses: loop_uses,
                            vertex: None,
                        }],
                    },
                    FaceName {
                        origins: vec![
                            Origin {
                                source: feature,
                                role: format!(
                                    "side/{}/{}",
                                    origins[edge].source, origins[edge].role
                                ),
                            },
                            origins[edge].clone(),
                        ],
                        support,
                        anchor,
                    },
                ));
            }
        }
        for side in sides {
            let (face, name) = side.expect("every authored edge bounds one wall");
            boundary.faces.push(face);
            names.push(name);
        }
        boundary.shells.push(Shell {
            faces: (0..boundary.faces.len()).collect(),
            reversed: false,
        });
        let solid = brep::Solid::new(boundary, tolerance, None).map_err(failure)?;
        Solid::from_authored_brep(solid, feature, AuthoredNaming { faces: names }, policy)
    }

    /// The support of one side wall, whether its parameterisation faces
    /// outward, and its planar support and authored frames.
    #[allow(clippy::type_complexity)]
    fn wall(
        &self,
        edge: usize,
        ring: &Ring,
        [here, next, above]: [P; 3],
        [start, top]: [usize; 2],
        boundary: &Boundary,
    ) -> Result<
        (
            Surface,
            bool,
            Option<ReferenceFrame>,
            Option<ReferenceFrame>,
        ),
        GeometryError,
    > {
        let path = &self.start[edge];
        if let Path::Line([a, b]) = path {
            let x = unit(sub(next, here))?;
            let z = unit(cross(x, sub(above, here)))?;
            let support = ReferenceFrame {
                origin_mm: here,
                u: x,
                v: cross(z, x),
            };
            let u = unit(sub(*b, *a))?;
            let rail = sub(boundary.vertices[top], boundary.vertices[start]);
            let v = unit(sub(rail, scale(u, dot(rail, u))))?;
            let anchor = ReferenceFrame {
                origin_mm: boundary.vertices[start],
                u,
                v,
            };
            return Ok((
                Surface::Plane(Frame {
                    origin: here,
                    x,
                    y: support.v,
                    z,
                }),
                true,
                Some(support),
                Some(anchor),
            ));
        }
        let curve = path.curve(1e-7)?;
        match (self.wall, path, &self.end[edge]) {
            (Wall::Translation(direction), ..) => Ok((
                Surface::Extrusion { curve, direction },
                ring.forward,
                None,
                None,
            )),
            (
                Wall::Ruling,
                Path::Conic {
                    center: c0,
                    axes: a0,
                    range,
                },
                Path::Conic {
                    center: c1,
                    axes: a1,
                    ..
                },
            ) if range[1] - range[0] == TAU => {
                let axis = sub(*c1, *c0);
                let height = norm(axis);
                let up = scale(axis, 1. / height);
                let [r0, r1] = [norm(a0[0]), norm(a1[0])];
                let coaxial = norm(cross(axis, cross(a0[0], a0[1]))) <= 1e-12 * height * r0 * r0;
                if coaxial {
                    // Keep a positive semi-angle: the frame sits on the
                    // smaller circle and opens toward the larger one.
                    let x = unit(a0[0])?;
                    let (origin, z, radius, growth) = if r1 >= r0 {
                        (*c0, up, r0, r1 - r0)
                    } else {
                        (*c1, up.map(|x| -x), r1, r0 - r1)
                    };
                    return Ok((
                        Surface::Cone {
                            frame: Frame {
                                origin,
                                x,
                                y: cross(z, x),
                                z,
                            },
                            radius,
                            half_angle: (growth / height).atan(),
                        },
                        ring.outer,
                        None,
                        None,
                    ));
                }
                Ok((
                    ruled_circles([*c0, *c1], [*a0, *a1])?,
                    ring.forward,
                    None,
                    None,
                ))
            }
            _ => Err(error(
                OPERATION,
                GeometryErrorKind::Unsupported,
                "Straight rulings currently join only lines and full circles",
            )),
        }
    }
}

/// The exact rational ruling between two circles with corresponding
/// parameterisations: degree two by one, with equal circle weights.
fn ruled_circles(centers: [P; 2], axes: [[P; 2]; 2]) -> Result<Surface, GeometryError> {
    let rim = |center: P, [a, b]: [P; 2]| {
        [
            [1., 0.],
            [1., 1.],
            [0., 1.],
            [-1., 1.],
            [-1., 0.],
            [-1., -1.],
            [0., -1.],
            [1., -1.],
            [1., 0.],
        ]
        .map(|[x, y]| add(center, add(scale(a, x), scale(b, y))))
    };
    let rims = [rim(centers[0], axes[0]), rim(centers[1], axes[1])];
    let mut points = Vec::with_capacity(18);
    let mut weights = Vec::with_capacity(18);
    for i in 0..9 {
        let weight = if i % 2 == 1 { FRAC_1_SQRT_2 } else { 1. };
        for rim in &rims {
            points.push(rim[i]);
            weights.push(weight);
        }
    }
    let surface = NurbsSurface::new(
        [2, 1],
        [9, 2],
        points,
        Some(weights),
        [vec![0., 0.25, 0.5, 0.75, 1.], vec![0., 1.]],
        [vec![3, 2, 2, 2, 3], vec![2, 2]],
    )
    .map_err(failure)?;
    Ok(Surface::nurbs(surface, 1e-9))
}

impl Solid {
    pub(super) fn from_authored_brep(
        solid: brep::Solid,
        feature: EntityId,
        naming: AuthoredNaming,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        Self::from_retained_brep(solid, feature, Rc::new(naming), policy)
    }

    /// Translate a certified curved section by an oblique displacement.
    pub(crate) fn rust_authored_translation(
        profile: &CurvedProfile,
        displacement: P,
        feature: EntityId,
        origins: &[Origin],
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        policy.point(displacement)?;
        policy.length(norm(displacement))?;
        let height = dot(displacement, profile.frame.normal());
        policy.length(height.abs())?;
        certify(profile, feature, origins, policy)?;
        let map = Affine::translation(displacement);
        let start: Vec<_> = edges(profile)
            .map(|e| Path::sketch(profile.frame, e))
            .collect();
        let end = start.iter().map(|path| path.mapped(map)).collect();
        let mut top = profile.frame;
        top.origin = map.point(top.origin);
        Sweep {
            profile,
            start,
            end,
            wall: Wall::Translation(displacement),
            up: displacement,
            caps: Caps {
                anchors: [profile.frame, top],
            },
        }
        .build(feature, origins, policy)
    }

    /// Rule a certified section to its tapered end section.
    pub(crate) fn rust_authored_taper(
        profile: &CurvedProfile,
        end: &CurvedProfile,
        displacement: P,
        feature: EntityId,
        origins: &[Origin],
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        certify(profile, feature, origins, policy)?;
        certify(end, feature, origins, policy)?;
        let start: Vec<_> = edges(profile)
            .map(|e| Path::sketch(profile.frame, e))
            .collect();
        // Both sections carry exactly one curve per authored origin.
        let end_paths: Vec<_> = edges(end).map(|e| Path::sketch(end.frame, e)).collect();
        Sweep {
            profile,
            start,
            end: end_paths,
            wall: Wall::Ruling,
            up: displacement,
            caps: Caps {
                anchors: [profile.frame, end.frame],
            },
        }
        .build(feature, origins, policy)
    }

    /// Draft both halves of a two-sided sweep from the original section and
    /// join them there, so the waist keeps the authored section exactly. The
    /// far backward cap is the logical start cap.
    #[allow(clippy::too_many_arguments)]
    pub fn rust_extrude_tapered_halves_named(
        profile: &CurvedProfile,
        direction: P,
        [forward_mm, backward_mm]: [f64; 2],
        taper_rad: f64,
        feature: EntityId,
        origins: &[Origin],
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        if direction.iter().any(|x| !x.is_finite()) || (norm(direction) - 1.).abs() > 1e-10 {
            return Err(input("Extrusion direction must be finite and unit length"));
        }
        policy.length(forward_mm)?;
        policy.length(backward_mm)?;
        let half = |length: f64| {
            let displacement = scale(direction, length);
            let end = crate::taper::end_profile(profile, displacement, taper_rad, policy)?;
            Self::rust_authored_taper(profile, &end, displacement, feature, origins, policy)
        };
        let forward = half(forward_mm)?;
        let backward =
            half(-backward_mm)?.with_feature_cap_roles(feature, "cap.end", "cap.start")?;
        forward
            .boolean(&backward, BooleanOperation::Union, policy)?
            .ok_or_else(|| {
                error(
                    OPERATION,
                    GeometryErrorKind::InvalidResult,
                    "Two drafted halves must join into one solid",
                )
            })
    }

    /// Translate a certified section along a unit direction between an
    /// optional start plane and an end plane.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn rust_authored_plane_sweep(
        profile: &CurvedProfile,
        direction: P,
        start: Option<PlaneFrame>,
        end: PlaneFrame,
        from_sketch: bool,
        feature: EntityId,
        origins: &[Origin],
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        if direction.iter().any(|x| !x.is_finite()) || (norm(direction) - 1.).abs() > 1e-10 {
            return Err(input("Extrusion direction must be finite and unit length"));
        }
        let crossing = dot(direction, profile.frame.normal());
        if !crossing.is_finite() || crossing.abs() < 1e-10 {
            return Err(input("Extrusion direction must cross the sketch plane"));
        }
        for plane in start.iter().chain([&end]) {
            plane.checked(policy)?;
            let crossing = dot(direction, plane.normal());
            if !crossing.is_finite() || crossing.abs() < 1e-10 {
                return Err(input(
                    "Extrusion direction must cross its termination plane",
                ));
            }
        }
        certify(profile, feature, origins, policy)?;
        // Travel from a sketch point to a plane is affine in sketch
        // coordinates: `offset + gradient . p`.
        let travel = |plane: PlaneFrame| {
            let n = plane.normal();
            let crossing = dot(direction, n);
            let offset = dot(n, sub(plane.origin, profile.frame.origin)) / crossing;
            let gradient = [
                -dot(n, profile.frame.u) / crossing,
                -dot(n, profile.frame.v) / crossing,
            ];
            (offset, gradient)
        };
        let range = |(offset, gradient): (f64, P2)| {
            linear_range(&profile.outer, gradient, policy)
                .map(|[low, high]| [offset + low, offset + high])
        };
        let (end_offset, end_gradient) = travel(end);
        let reach = range((end_offset, end_gradient))?;
        let excessive =
            |value: f64| !value.is_finite() || value.abs() > policy.maximum_coordinate_mm;
        if excessive(reach[0]) || excessive(reach[1]) {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Input,
                "Plane-terminated sweep exceeds the length limit",
            ));
        }
        match start {
            Some(plane) => {
                let (offset, gradient) = travel(plane);
                let behind = range((offset, gradient))?;
                if from_sketch && behind[1] > -policy.minimum_feature_mm {
                    return Err(input(
                        "The start plane must lie behind the whole sketch section",
                    ));
                }
                let material = range((
                    end_offset - offset,
                    [end_gradient[0] - gradient[0], end_gradient[1] - gradient[1]],
                ))?;
                if material[0] < policy.minimum_feature_mm {
                    return Err(input(
                        "Termination planes cross or leave insufficient material",
                    ));
                }
            }
            None if reach[0] < policy.minimum_feature_mm => {
                return Err(input(
                    "The end plane must lie ahead of the whole sketch section",
                ));
            }
            None => {}
        }
        let sketch: Vec<_> = edges(profile)
            .map(|e| Path::sketch(profile.frame, e))
            .collect();
        let first = match start {
            Some(plane) => {
                let map = Affine::projection(direction, plane);
                sketch.iter().map(|path| path.mapped(map)).collect()
            }
            None => sketch.clone(),
        };
        let map = Affine::projection(direction, end);
        let last = sketch.iter().map(|path| path.mapped(map)).collect();
        Sweep {
            profile,
            start: first,
            end: last,
            wall: Wall::Translation(direction),
            up: direction,
            caps: Caps {
                anchors: [start.unwrap_or(profile.frame), end],
            },
        }
        .build(feature, origins, policy)
    }
}

fn edges(profile: &CurvedProfile) -> impl Iterator<Item = &CurvedEdge2d> {
    profile.outer.iter().chain(profile.holes.iter().flatten())
}

/// Admit the section through the certified normal-sweep provider, which
/// rejects open, crossing, overlapping or misnested wires.
fn certify(
    profile: &CurvedProfile,
    feature: EntityId,
    origins: &[Origin],
    policy: GeometryPolicy,
) -> Result<(), GeometryError> {
    profile.checked(policy)?;
    super::construction::validate_origins(feature, origins, edges(profile).count())?;
    Solid::rust_extrude_normal_profile_named(profile, 1., feature, origins, policy).map(drop)
}

#[cfg(test)]
mod tests;
