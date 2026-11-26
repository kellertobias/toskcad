//! The exact retained boundary of a prism whose cap edges are filleted or
//! chamfered.
//!
//! In the local chart the profile ring runs counterclockwise in the xy
//! plane and the prism spans `z0 <= z <= z1`. Every horizontal section of
//! the treated prism is the profile with each treated cap edge offset
//! inward: a constant-radius fillet of size `r` offsets an edge by
//! `r (1 - cos t)` at depth `r (1 - sin t)` below the cap, a chamfer by
//! `d - depth`. The offsets of consecutive edges meet at the corner nearest
//! the original vertex, so
//!
//! - a fillet blend is a quarter cylinder about a line edge and a torus
//!   (a sphere when the spine radius vanishes) about an arc edge, a chamfer
//!   a plane or a cone at 45 degrees;
//! - two blends meet at a convex vertex along the curve where both offsets
//!   agree, and a blend ends against the unchanged side surface of an
//!   untreated neighbour along the curve where its offset meets that side.
//!
//! End curves are exact conics where they are planar: the bisector ellipse
//! between two line blends, the section ellipse of a line blend by a flat
//! neighbour, the meridian circle where a torus meets a plane through its
//! axis, and the cross-section circle across a tangent vertex. Other end
//! curves are carried as C1 piecewise cubic Hermite curves through exact
//! points and tangents of the section law, refined until the sampled
//! deviation from the exact curve and from both surfaces stays below a
//! sixteenth of the tolerance. The boundary validates and integrates itself
//! and must reproduce the section integral of its volume.
use super::super::authored_sweep::{AuthoredNaming, FaceName};
use super::super::imported_brep_body::failure;
use super::super::*;
use super::profile::{self, Carrier, Corner, P2, Seg};
use multicad_rust_kernel::brep::{
    self, Boundary, Curve, Edge, Face, Frame, Loop, NurbsCurve, Shell, Surface,
};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, TAU};

type P3 = [f64; 3];

const SAMPLES: usize = 8;
const MAX_SEGMENTS: usize = 1 << 12;
const VOLUME_NODES: usize = 64;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Kind {
    Fillet,
    Chamfer,
}

/// A cap-treated prism in its local chart.
pub(super) struct Plan {
    pub(super) frame: PlaneFrame,
    pub(super) ring: Vec<Seg>,
    pub(super) span: [f64; 2],
    pub(super) kind: Kind,
    pub(super) size: f64,
    /// Blend face names of the treated cap edges at the bottom and the top.
    pub(super) caps: [Vec<Option<Vec<Origin>>>; 2],
    pub(super) cap_names: [FaceName; 2],
    pub(super) side_names: Vec<FaceName>,
}

pub(super) fn refused(message: impl Into<String>) -> GeometryError {
    error("edge treatment", GeometryErrorKind::Unsupported, message)
}

fn invalid() -> GeometryError {
    error(
        "edge treatment",
        GeometryErrorKind::InvalidResult,
        "Treated prism boundary does not reproduce its section volume",
    )
}

fn sub3(a: P3, b: P3) -> P3 {
    std::array::from_fn(|i| a[i] - b[i])
}
fn mix3(a: P3, s: f64, b: P3) -> P3 {
    std::array::from_fn(|i| a[i] + s * b[i])
}
fn cross3(a: P3, b: P3) -> P3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm3(a: P3) -> f64 {
    dot(a, a).sqrt()
}
fn unit3(a: P3) -> P3 {
    a.map(|x| x / norm3(a))
}
fn flat(p: P2) -> P3 {
    [p[0], p[1], 0.]
}

/// The section law of a treated edge: inset and depth below the cap at
/// curve parameter `t` in `[0, 1]`, with their rates.
fn law(kind: Kind, size: f64, t: f64) -> ([f64; 2], [f64; 2]) {
    match kind {
        Kind::Fillet => {
            let (s, c) = (FRAC_PI_2 * t).portable_sin_cos();
            (
                [size * (1. - c), size * (1. - s)],
                [size * s * FRAC_PI_2, -size * c * FRAC_PI_2],
            )
        }
        Kind::Chamfer => ([size * t, size * (1. - t)], [size, -size]),
    }
}

/// Deviation of a section point `(inset, depth)` from the blend of size
/// `size`.
fn blend_residual(kind: Kind, size: f64, inset: f64, depth: f64) -> f64 {
    match kind {
        Kind::Fillet => ((size - inset).portable_hypot(size - depth) - size).abs(),
        Kind::Chamfer => (inset + depth - size).abs() * FRAC_PI_4.portable_cos(),
    }
}

struct Chart {
    frame: PlaneFrame,
    normal: P3,
}

impl Chart {
    fn point(&self, p: P3) -> P3 {
        std::array::from_fn(|i| {
            self.frame.origin[i]
                + p[0] * self.frame.u[i]
                + p[1] * self.frame.v[i]
                + p[2] * self.normal[i]
        })
    }
    fn vector(&self, d: P3) -> P3 {
        std::array::from_fn(|i| {
            d[0] * self.frame.u[i] + d[1] * self.frame.v[i] + d[2] * self.normal[i]
        })
    }
    fn frame(&self, origin: P3, x: P3, y: P3) -> Frame {
        let (x, y) = (unit3(x), unit3(y));
        Frame {
            origin: self.point(origin),
            x: self.vector(x),
            y: self.vector(y),
            z: self.vector(cross3(x, y)),
        }
    }
}

impl Plan {
    fn flag(&self, level: usize, k: usize) -> bool {
        self.caps[level][k % self.ring.len()].is_some()
    }
    fn seg(&self, k: usize) -> Seg {
        let n = self.ring.len();
        self.ring[(k + n) % n]
    }
    fn previous(&self, k: usize) -> Seg {
        self.seg(k + self.ring.len() - 1)
    }
    /// Height of the cap at `level` and the direction into the material.
    fn level(&self, level: usize) -> (f64, f64) {
        if level == 1 {
            (self.span[1], -1.)
        } else {
            (self.span[0], 1.)
        }
    }
    fn at(&self, level: usize, p: P2, depth: f64) -> P3 {
        let (z, inward) = self.level(level);
        [p[0], p[1], z + inward * depth]
    }
    /// A treated convex arc whose fillet spine collapses to its center.
    fn collapsed(&self, level: usize, k: usize) -> bool {
        self.kind == Kind::Fillet
            && self.flag(level, k)
            && matches!(self.seg(k), Seg::Arc { r, cw: false, .. } if r == self.size)
    }
    fn inset_of(&self, level: usize, k: usize) -> f64 {
        if self.flag(level, k) { self.size } else { 0. }
    }
}

/// One vertex of the ring at one cap: where its side edge ends (`low`) and
/// where the cap face begins (`cap`). They coincide when neither edge at
/// the vertex is treated.
#[derive(Clone, Copy)]
struct Joint {
    low: usize,
    cap: usize,
}

pub(super) fn build(
    plan: &Plan,
    feature: EntityId,
    policy: GeometryPolicy,
) -> Result<Solid, GeometryError> {
    let n = plan.ring.len();
    let size = plan.size;
    let tolerance = policy.linear_tolerance_mm;
    let minimum = policy.minimum_feature_mm;
    let chart = Chart {
        frame: plan.frame,
        normal: plan.frame.normal(),
    };
    let levels = (0..2).filter(|&l| (0..n).any(|k| plan.flag(l, k))).count();
    if plan.span[1] - plan.span[0] - levels as f64 * size < minimum {
        return Err(refused(
            "The edge treatment is too large for the height of the body",
        ));
    }
    let corners: Vec<Corner> = (0..n)
        .map(|k| {
            profile::corner(plan.previous(k), plan.seg(k))
                .ok_or_else(|| refused("The body has a cusp where two treated faces would meet"))
        })
        .collect::<Result<_, _>>()?;
    for level in 0..2 {
        for (k, &corner) in corners.iter().enumerate() {
            let (a, b) = (plan.flag(level, k + n - 1), plan.flag(level, k));
            match corner {
                _ if !(a || b) => {}
                Corner::Concave => {
                    return Err(refused(
                        "Edge treatments cannot yet meet at a concave corner of the cap; select edges that end at convex corners",
                    ));
                }
                Corner::Smooth if a != b => {
                    return Err(refused(
                        "The treated edge continues tangentially into an unselected edge; select the whole tangent chain",
                    ));
                }
                _ => {}
            }
        }
        for k in 0..n {
            if !plan.flag(level, k) {
                continue;
            }
            if let Seg::Arc { r, cw: false, .. } = plan.seg(k)
                && !plan.collapsed(level, k)
                && r - size < minimum
            {
                return Err(refused(
                    "The edge treatment is too large for the radius of the curved edge",
                ));
            }
        }
    }
    // Section corners at the cap of each level.
    let mut caps = [vec![[0.; 2]; n], vec![[0.; 2]; n]];
    for (level, cap) in caps.iter_mut().enumerate() {
        for k in 0..n {
            let deltas = [plan.inset_of(level, k + n - 1), plan.inset_of(level, k)];
            cap[k] = profile::inset_corner(plan.previous(k), plan.seg(k), corners[k], deltas)
                .ok_or_else(|| refused("The edge treatment is too large for the faces it joins"))?;
        }
        for k in 0..n {
            if plan.collapsed(level, k) {
                let Seg::Arc { c, .. } = plan.seg(k) else {
                    unreachable!("collapsed arc")
                };
                if plan.collapsed(level, k + n - 1) || plan.collapsed(level, k + 1) {
                    return Err(refused(
                        "Two fully rounded corners meet; reduce the edge treatment",
                    ));
                }
                cap[k] = c;
                cap[(k + 1) % n] = c;
            }
        }
        // Each cap curve keeps a positive length in its own direction.
        for k in 0..n {
            if plan.collapsed(level, k) {
                continue;
            }
            let seg = plan.seg(k);
            let [from, to] = [cap[k], cap[(k + 1) % n]];
            let carried = match seg.carrier(plan.inset_of(level, k)) {
                Carrier::Line { direction, .. } => profile::dot(profile::sub(to, from), direction),
                Carrier::Circle { center, radius } => {
                    let Seg::Arc { cw, .. } = seg else {
                        unreachable!("circle carrier")
                    };
                    let (p, q) = (profile::sub(from, center), profile::sub(to, center));
                    let angle = profile::cross(p, q).portable_atan2(profile::dot(p, q));
                    let turn = if cw { -angle } else { angle }.rem_euclid(TAU);
                    if turn > seg.sweep() + (TAU - seg.sweep()) / 2. {
                        -1.
                    } else {
                        radius * turn
                    }
                }
            };
            if carried < minimum {
                return Err(refused(
                    "The edge treatment is too large for the faces it joins",
                ));
            }
        }
    }

    // Vertices.
    let mut vertices: Vec<P3> = Vec::new();
    let mut vertex = |p: P3| {
        vertices.push(chart.point(p));
        vertices.len() - 1
    };
    let mut joints: [Vec<Joint>; 2] = [Vec::with_capacity(n), Vec::with_capacity(n)];
    for level in 0..2 {
        for k in 0..n {
            let v = plan.seg(k).start();
            let any = plan.flag(level, k + n - 1) || plan.flag(level, k);
            let low = vertex(plan.at(level, v, if any { size } else { 0. }));
            let cap = if !any {
                low
            } else if k > 0 && plan.collapsed(level, k - 1) {
                joints[level][k - 1].cap
            } else {
                vertex(plan.at(level, caps[level][k], 0.))
            };
            joints[level].push(Joint { low, cap });
        }
        if plan.collapsed(level, n - 1) {
            joints[level][0].cap = joints[level][n - 1].cap;
        }
    }
    let local = |level: usize, k: usize, which: bool| -> P3 {
        let v = plan.seg(k).start();
        let any = plan.flag(level, k + n - 1) || plan.flag(level, k);
        if which {
            plan.at(level, caps[level][k % n], 0.)
        } else {
            plan.at(level, v, if any { size } else { 0. })
        }
    };

    let mut edges: Vec<Edge> = Vec::new();
    let mut add = |curve: Curve, ends: [usize; 2], vertices: &[P3]| {
        let points = ends.map(|v| vertices[v]);
        edges.push(Edge::new(curve, ends, points, tolerance).map_err(failure)?);
        Ok::<_, GeometryError>(edges.len() - 1)
    };
    let line = |a: P3, b: P3| Curve::Line {
        origin: chart.point(a),
        direction: chart.vector(sub3(b, a)),
    };
    // The carrier of edge `k` offset by `delta` at height `z`.
    let carried = |k: usize, delta: f64, z: f64, from: P3, to: P3| -> Curve {
        match plan.seg(k).carrier(delta) {
            Carrier::Line { .. } => line(from, to),
            Carrier::Circle { center, radius } => {
                let Seg::Arc { cw, .. } = plan.seg(k) else {
                    unreachable!("circle carrier")
                };
                let y = if cw { [0., -1., 0.] } else { [0., 1., 0.] };
                Curve::Circle {
                    frame: chart.frame([center[0], center[1], z], [1., 0., 0.], y),
                    radius,
                }
            }
        }
    };

    let vertical: Vec<usize> = (0..n)
        .map(|k| {
            let ends = [joints[0][k].low, joints[1][k].low];
            add(
                line(local(0, k, false), local(1, k, false)),
                ends,
                &vertices,
            )
        })
        .collect::<Result<_, _>>()?;
    let mut ends = [vec![None; n], vec![None; n]];
    let mut rims = [vec![None; n], vec![None; n]];
    let mut tangents = [vec![None; n], vec![None; n]];
    for level in 0..2 {
        let (z, inward) = plan.level(level);
        for k in 0..n {
            let Joint { low, cap } = joints[level][k];
            if low != cap {
                let curve = end_curve(
                    plan,
                    &chart,
                    level,
                    k,
                    corners[k],
                    caps[level][k],
                    tolerance,
                )?;
                ends[level][k] = Some(add(curve, [low, cap], &vertices)?);
            }
            let next = (k + 1) % n;
            if !plan.collapsed(level, k) {
                let (from, to) = (local(level, k, true), local(level, next, true));
                let curve = carried(k, plan.inset_of(level, k), z, from, to);
                rims[level][k] = Some(add(
                    curve,
                    [joints[level][k].cap, joints[level][next].cap],
                    &vertices,
                )?);
            }
            if plan.flag(level, k) {
                let (from, to) = (local(level, k, false), local(level, next, false));
                let curve = carried(k, 0., z + inward * size, from, to);
                tangents[level][k] = Some(add(
                    curve,
                    [joints[level][k].low, joints[level][next].low],
                    &vertices,
                )?);
            }
        }
    }

    let mut faces = Vec::with_capacity(2 + 3 * n);
    let mut names = Vec::with_capacity(faces.capacity());
    let plane = |origin: P3, x: P3, y: P3| Surface::Plane(chart.frame(origin, x, y));
    // Caps: counterclockwise about the outward normal.
    for (level, level_rims) in rims.iter().enumerate() {
        let (z, _) = plan.level(level);
        let mut uses: Vec<(usize, bool)> = level_rims
            .iter()
            .flatten()
            .map(|&edge| (edge, true))
            .collect();
        let y = if level == 1 {
            [0., 1., 0.]
        } else {
            [0., -1., 0.]
        };
        if level == 0 {
            uses.reverse();
            uses.iter_mut().for_each(|use_| use_.1 = false);
        }
        faces.push(Face {
            surface: plane([0., 0., z], [1., 0., 0.], y),
            same_sense: true,
            loops: vec![Loop { uses, vertex: None }],
        });
        names.push(plan.cap_names[level].clone());
    }
    // Sides.
    for k in 0..n {
        let next = (k + 1) % n;
        let previous = (k + n - 1) % n;
        let mut uses = Vec::new();
        if let Some(edge) = tangents[0][k] {
            uses.push((edge, true));
        } else {
            if plan.flag(0, previous) {
                uses.push((ends[0][k].expect("bottom end"), true));
            }
            uses.push((rims[0][k].expect("bottom rim"), true));
            if plan.flag(0, next) {
                uses.push((ends[0][next].expect("bottom end"), false));
            }
        }
        uses.push((vertical[next], true));
        if let Some(edge) = tangents[1][k] {
            uses.push((edge, false));
        } else {
            if plan.flag(1, next) {
                uses.push((ends[1][next].expect("top end"), true));
            }
            uses.push((rims[1][k].expect("top rim"), false));
            if plan.flag(1, previous) {
                uses.push((ends[1][k].expect("top end"), false));
            }
        }
        uses.push((vertical[k], false));
        let surface = match plan.seg(k) {
            Seg::Line { a, b } => {
                let t = flat(profile::unit(profile::sub(b, a)));
                plane([a[0], a[1], plan.span[0]], t, [0., 0., 1.])
            }
            Seg::Arc { c, r, .. } => Surface::Cylinder {
                frame: chart.frame([c[0], c[1], plan.span[0]], [1., 0., 0.], [0., 1., 0.]),
                radius: r,
            },
        };
        let same_sense = !matches!(plan.seg(k), Seg::Arc { cw: true, .. });
        faces.push(Face {
            surface,
            same_sense,
            loops: vec![Loop { uses, vertex: None }],
        });
        names.push(plan.side_names[k].clone());
    }
    // Blends.
    for (level, level_tangents) in tangents.iter().enumerate() {
        let (z, inward) = plan.level(level);
        let up = [0., 0., -inward];
        for k in 0..n {
            let Some(origins) = &plan.caps[level][k] else {
                continue;
            };
            let next = (k + 1) % n;
            let tangent = level_tangents[k].expect("blend tangent");
            let mut uses = Vec::with_capacity(4);
            if level == 1 {
                uses.push((tangent, true));
                uses.push((ends[1][next].expect("end"), true));
                if let Some(rim) = rims[1][k] {
                    uses.push((rim, false));
                }
                uses.push((ends[1][k].expect("end"), false));
            } else {
                if let Some(rim) = rims[0][k] {
                    uses.push((rim, true));
                }
                uses.push((ends[0][next].expect("end"), false));
                uses.push((tangent, false));
                uses.push((ends[0][k].expect("end"), true));
            }
            let seg = plan.seg(k);
            let (surface, same_sense) = match (plan.kind, seg) {
                (Kind::Fillet, Seg::Line { a, b }) => {
                    let t = flat(profile::unit(profile::sub(b, a)));
                    let inward_normal = flat(seg.inward(a));
                    let spine = mix3(flat(a), size, inward_normal);
                    let x = inward_normal.map(|x| -x);
                    (
                        Surface::Cylinder {
                            frame: chart.frame(
                                [spine[0], spine[1], z + inward * size],
                                x,
                                cross3(t, x),
                            ),
                            radius: size,
                        },
                        true,
                    )
                }
                (Kind::Fillet, Seg::Arc { c, r, cw, .. }) => {
                    let center = [c[0], c[1], z + inward * size];
                    if plan.collapsed(level, k) {
                        let m = flat(seg.bisector());
                        let axis = [-m[1], m[0], 0.];
                        (
                            Surface::Sphere {
                                frame: chart.frame(center, up, cross3(axis, up)),
                                radius: size,
                            },
                            true,
                        )
                    } else {
                        (
                            Surface::Torus {
                                frame: chart.frame(center, [1., 0., 0.], [0., 1., 0.]),
                                major: if cw { r + size } else { r - size },
                                minor: size,
                            },
                            true,
                        )
                    }
                }
                (Kind::Chamfer, Seg::Line { a, b }) => {
                    let t = flat(profile::unit(profile::sub(b, a)));
                    let out = flat(seg.inward(a)).map(|x| -x);
                    let normal = unit3(std::array::from_fn(|i| out[i] + up[i]));
                    (
                        plane([a[0], a[1], z + inward * size], t, cross3(normal, t)),
                        true,
                    )
                }
                (Kind::Chamfer, Seg::Arc { c, r, cw, .. }) => {
                    let down = [0., 0., inward];
                    if cw {
                        // Radius grows from `r` at the side to `r + d` at the cap.
                        let frame = chart.frame(
                            [c[0], c[1], z + inward * size],
                            [1., 0., 0.],
                            cross3(up, [1., 0., 0.]),
                        );
                        (
                            Surface::Cone {
                                frame,
                                radius: r,
                                half_angle: FRAC_PI_4,
                            },
                            false,
                        )
                    } else {
                        let frame =
                            chart.frame([c[0], c[1], z], [1., 0., 0.], cross3(down, [1., 0., 0.]));
                        (
                            Surface::Cone {
                                frame,
                                radius: r - size,
                                half_angle: FRAC_PI_4,
                            },
                            true,
                        )
                    }
                }
            };
            faces.push(Face {
                surface,
                same_sense,
                loops: vec![Loop { uses, vertex: None }],
            });
            names.push(FaceName {
                origins: origins.clone(),
                support: None,
                anchor: None,
            });
        }
    }
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
    let expected = section_volume(plan);
    if (solid.volume_mm3() - expected).abs() > solid.volume_error_mm3().max(1e-9 * expected) {
        return Err(invalid());
    }
    Solid::from_retained_brep(
        solid,
        feature,
        Rc::new(AuthoredNaming { faces: names }),
        policy,
    )
}

/// The volume as the integral of exact section areas over the height.
fn section_volume(plan: &Plan) -> f64 {
    let n = plan.ring.len();
    let corners: Vec<Corner> = (0..n)
        .map(|k| profile::corner(plan.previous(k), plan.seg(k)).expect("checked corner"))
        .collect();
    let area = |level: usize, delta: f64| {
        let deltas: Vec<f64> = (0..n)
            .map(|k| if plan.flag(level, k) { delta } else { 0. })
            .collect();
        let points: Vec<P2> = (0..n)
            .map(|k| {
                profile::inset_corner(
                    plan.previous(k),
                    plan.seg(k),
                    corners[k],
                    [deltas[(k + n - 1) % n], deltas[k]],
                )
                .expect("checked corner")
            })
            .collect();
        profile::twice_area(&plan.ring, &deltas, &points) / 2.
    };
    let base = area(0, 0.);
    let nodes = profile::gauss_legendre(VOLUME_NODES);
    let band = |level: usize| {
        if !(0..n).any(|k| plan.flag(level, k)) {
            return base * plan.size;
        }
        nodes
            .iter()
            .map(|&(x, w)| {
                let ([delta, _], [_, rate]) = law(plan.kind, plan.size, (x + 1.) / 2.);
                // Depth falls as t rises; |d depth / dt| / 2 maps [-1, 1].
                w * area(level, delta) * -rate / 2.
            })
            .sum()
    };
    base * (plan.span[1] - plan.span[0] - 2. * plan.size) + band(0) + band(1)
}

/// The curve from the side corner `(V, depth size)` up to the cap corner at
/// vertex `k`.
fn end_curve(
    plan: &Plan,
    chart: &Chart,
    level: usize,
    k: usize,
    kind: Corner,
    cap: P2,
    tolerance: f64,
) -> Result<Curve, GeometryError> {
    let n = plan.ring.len();
    let size = plan.size;
    let (a, b) = (plan.previous(k), plan.seg(k));
    let (sa, sb) = (plan.flag(level, k + n - 1), plan.flag(level, k));
    let v = b.start();
    let (z, inward) = plan.level(level);
    let up = [0., 0., -inward];
    let low = plan.at(level, v, size);
    let top = plan.at(level, cap, 0.);
    let centre = plan.at(level, cap, size);
    let reach = profile::norm(profile::sub(v, cap));
    let toward = flat(profile::unit(profile::sub(v, cap)));
    let conic = |semi: f64| {
        let frame = chart.frame(centre, toward, up);
        if (semi - size).abs() <= 1e-12 * size {
            Curve::Circle {
                frame,
                radius: size,
            }
        } else {
            Curve::Ellipse {
                frame,
                semiaxes: [semi, size],
            }
        }
    };
    let straight = || Curve::Line {
        origin: chart.point(low),
        direction: chart.vector(sub3(top, low)),
    };
    let lines = matches!((a, b), (Seg::Line { .. }, Seg::Line { .. }));
    // A torus meets a plane through its axis along a meridian circle, a
    // cone along a generator.
    let radial = |arc: Seg, other: Seg| match (arc, other) {
        (Seg::Arc { c, .. }, Seg::Line { .. }) => other.inset(c).abs() <= tolerance,
        _ => false,
    };
    let exact = match kind {
        Corner::Smooth => true,
        Corner::Convex if sa && sb => lines,
        Corner::Convex => lines || (sa && radial(a, b)) || (sb && radial(b, a)),
        Corner::Concave => false,
    };
    if exact {
        return Ok(match plan.kind {
            Kind::Fillet => conic(reach),
            Kind::Chamfer => straight(),
        });
    }
    let point = |t: f64| -> Option<(P3, P3)> {
        let ([delta, depth], [delta_rate, depth_rate]) = law(plan.kind, size, t);
        let insets = [if sa { delta } else { 0. }, if sb { delta } else { 0. }];
        let p = profile::inset_corner(a, b, kind, insets)?;
        let rates = [
            if sa { delta_rate } else { 0. },
            if sb { delta_rate } else { 0. },
        ];
        let d = profile::inset_corner_rate(a, b, p, rates)?;
        Some((
            [p[0], p[1], z + inward * depth],
            [d[0], d[1], inward * depth_rate],
        ))
    };
    let residual = |q: P3| {
        let depth = (q[2] - z) * inward;
        let p = [q[0], q[1]];
        [(a, sa), (b, sb)]
            .into_iter()
            .map(|(seg, treated)| {
                let inset = seg.inset(p);
                if treated {
                    blend_residual(plan.kind, size, inset, depth)
                } else {
                    inset.abs()
                }
            })
            .fold(0., f64::max)
    };
    hermite(&point, &residual, tolerance)
        .map(|curve| Curve::nurbs(curve.map_points(|p| chart.point(p)), tolerance))
}

/// Local control polygon of a Hermite curve before it is placed.
struct LocalCurve {
    points: Vec<P3>,
    distinct: Vec<f64>,
    multiplicities: Vec<usize>,
}

impl LocalCurve {
    fn map_points(self, map: impl Fn(P3) -> P3) -> NurbsCurve {
        NurbsCurve::new(
            3,
            self.points.into_iter().map(map).collect(),
            None,
            self.distinct,
            self.multiplicities,
        )
        .expect("valid Hermite knots")
    }
}

/// C1 cubic Bezier pieces through exact points and tangents on `[0, 1]`.
fn hermite(
    exact: &dyn Fn(f64) -> Option<(P3, P3)>,
    residual: &dyn Fn(P3) -> f64,
    tolerance: f64,
) -> Result<LocalCurve, GeometryError> {
    let unresolved =
        || refused("The edge treatment corner is too close to tangency for an exact boundary");
    let mut segments = 8;
    loop {
        let h = 1. / segments as f64;
        let nodes = (0..=segments)
            .map(|k| exact(h * k as f64))
            .collect::<Option<Vec<_>>>()
            .ok_or_else(unresolved)?;
        let mut points = Vec::with_capacity(3 * segments + 1);
        let mut worst: f64 = 0.;
        for k in 0..segments {
            let [(p0, d0), (p1, d1)] = [nodes[k], nodes[k + 1]];
            let control = [p0, mix3(p0, h / 3., d0), mix3(p1, -h / 3., d1), p1];
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
                let bezier: P3 =
                    std::array::from_fn(|i| (0..4).map(|m| w[m] * control[m][i]).sum());
                let (on, _) = exact(h * (k as f64 + u)).ok_or_else(unresolved)?;
                worst = worst.max(norm3(sub3(bezier, on))).max(residual(bezier));
            }
        }
        if worst.is_finite() && worst <= tolerance / 16. {
            let mut multiplicities = vec![3; segments + 1];
            multiplicities[0] = 4;
            multiplicities[segments] = 4;
            return Ok(LocalCurve {
                points,
                distinct: (0..=segments).map(|k| h * k as f64).collect(),
                multiplicities,
            });
        }
        segments *= 2;
        if segments > MAX_SEGMENTS {
            return Err(unresolved());
        }
    }
}
