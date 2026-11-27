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
//!   untreated neighbour along the curve where its offset meets that side;
//! - where a treated edge continues tangentially into an untreated one,
//!   the blend stops at the end of the selected edge. A flat step face
//!   through the vertex, perpendicular to both edges, closes the blend's
//!   cross-section there, and the untreated neighbour keeps its sharp edge
//!   (CAD-421). In section, the offset edge ends on the common normal and
//!   a straight chord leads back to the vertex.
//!
//! Hole rings run clockwise so that every ring keeps the material on its
//! left; each cap is one face with a loop per ring.
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
    /// The edges of every profile ring, outer ring first. Each ring keeps
    /// the material on its left: the outer ring runs counterclockwise and
    /// hole rings clockwise.
    pub(super) ring: Vec<Seg>,
    /// The edge after and before each edge in its own ring.
    pub(super) next: Vec<usize>,
    pub(super) prev: Vec<usize>,
    /// The edge range of each ring.
    pub(super) loops: Vec<std::ops::Range<usize>>,
    pub(super) span: [f64; 2],
    pub(super) kind: Kind,
    pub(super) size: f64,
    /// Blend face names of the treated cap edges at the bottom and the top.
    pub(super) caps: [Vec<Option<Vec<Origin>>>; 2],
    pub(super) cap_names: [FaceName; 2],
    pub(super) side_names: Vec<FaceName>,
}

/// Successor and predecessor indices of rings stored one after another.
pub(super) fn ring_links(
    lengths: &[usize],
) -> (Vec<usize>, Vec<usize>, Vec<std::ops::Range<usize>>) {
    let (mut next, mut prev, mut loops) = (vec![], vec![], vec![]);
    let mut first = 0;
    for &len in lengths {
        for i in 0..len {
            next.push(first + (i + 1) % len);
            prev.push(first + (i + len - 1) % len);
        }
        loops.push(first..first + len);
        first += len;
    }
    (next, prev, loops)
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
        self.caps[level][k].is_some()
    }
    fn seg(&self, k: usize) -> Seg {
        self.ring[k]
    }
    fn previous(&self, k: usize) -> Seg {
        self.ring[self.prev[k]]
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
    /// Either edge at vertex `k` is treated at `level`.
    fn touched(&self, level: usize, k: usize) -> bool {
        self.flag(level, self.prev[k]) || self.flag(level, k)
    }
}

/// One vertex of the ring at one cap: where its side edge ends (`low`),
/// where the cap edge before it ends and where the cap edge after it
/// begins (`rim`). They coincide when neither edge at the vertex is
/// treated. The two rim vertices differ only at a step: a tangent vertex
/// with one treated edge, where a flat face perpendicular to the edges
/// closes the end of the blend.
#[derive(Clone, Copy)]
struct Joint {
    low: usize,
    rim: [usize; 2],
}

impl Joint {
    fn step(self) -> bool {
        self.rim[0] != self.rim[1]
    }
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
            let (a, b) = (plan.flag(level, plan.prev[k]), plan.flag(level, k));
            match corner {
                _ if !(a || b) => {}
                Corner::Concave => {
                    return Err(refused(
                        "Edge treatments cannot yet meet at a concave corner of the cap; select edges that end at convex corners",
                    ));
                }
                Corner::Smooth
                    if a != b
                        && (plan.collapsed(level, plan.prev[k]) || plan.collapsed(level, k)) =>
                {
                    return Err(refused(
                        "A fully rounded corner cannot end beside an unselected tangent edge; select the whole tangent chain or reduce the edge treatment",
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
    // Section corners at the cap of each level: where the edge before each
    // vertex ends and where the edge after it begins.
    let mut caps = [vec![[[0.; 2]; 2]; n], vec![[[0.; 2]; 2]; n]];
    for (level, cap) in caps.iter_mut().enumerate() {
        for k in 0..n {
            let deltas = [plan.inset_of(level, plan.prev[k]), plan.inset_of(level, k)];
            cap[k] = profile::joint_ends(plan.previous(k), plan.seg(k), corners[k], deltas)
                .ok_or_else(|| refused("The edge treatment is too large for the faces it joins"))?;
        }
        for k in 0..n {
            if plan.collapsed(level, k) {
                let Seg::Arc { c, .. } = plan.seg(k) else {
                    unreachable!("collapsed arc")
                };
                if plan.collapsed(level, plan.prev[k]) || plan.collapsed(level, plan.next[k]) {
                    return Err(refused(
                        "Two fully rounded corners meet; reduce the edge treatment",
                    ));
                }
                cap[k] = [c, c];
                cap[plan.next[k]] = [c, c];
            }
        }
        // Each cap curve keeps a positive length in its own direction.
        for k in 0..n {
            if plan.collapsed(level, k) {
                continue;
            }
            let seg = plan.seg(k);
            let [from, to] = [cap[k][1], cap[plan.next[k]][0]];
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

    // Vertices. The ends of a fully rounded corner share the vertex at its
    // center.
    let mut vertices: Vec<P3> = Vec::new();
    let mut vertex = |p: P3| {
        vertices.push(chart.point(p));
        vertices.len() - 1
    };
    let mut joints: [Vec<Joint>; 2] = [Vec::with_capacity(n), Vec::with_capacity(n)];
    for level in 0..2 {
        for (k, &[p, q]) in caps[level].iter().enumerate() {
            let v = plan.seg(k).start();
            let any = plan.touched(level, k);
            let low = vertex(plan.at(level, v, if any { size } else { 0. }));
            let rim = if !any || plan.collapsed(level, plan.prev[k]) {
                [low; 2]
            } else if p == q {
                [vertex(plan.at(level, p, 0.)); 2]
            } else {
                [vertex(plan.at(level, p, 0.)), vertex(plan.at(level, q, 0.))]
            };
            joints[level].push(Joint { low, rim });
        }
        for k in 0..n {
            if plan.collapsed(level, plan.prev[k]) {
                joints[level][k].rim = joints[level][plan.prev[k]].rim;
            }
        }
    }
    let rim_point = |level: usize, k: usize, side: usize| plan.at(level, caps[level][k][side], 0.);
    let low_point = |level: usize, k: usize| {
        let depth = if plan.touched(level, k) { size } else { 0. };
        plan.at(level, plan.seg(k).start(), depth)
    };
    // The rim side of a treated vertex where the blend ends: the treated
    // edge's side of a step.
    let blend_side = |level: usize, k: usize| usize::from(!plan.flag(level, plan.prev[k]));

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
            add(line(low_point(0, k), low_point(1, k)), ends, &vertices)
        })
        .collect::<Result<_, _>>()?;
    let mut ends = [vec![None; n], vec![None; n]];
    // At a step: the side edge up to the cap (`risers`) and the cap edge
    // across the end of the blend (`steps`).
    let mut risers = [vec![None; n], vec![None; n]];
    let mut steps = [vec![None; n], vec![None; n]];
    let mut rims = [vec![None; n], vec![None; n]];
    let mut tangents = [vec![None; n], vec![None; n]];
    for level in 0..2 {
        let (z, inward) = plan.level(level);
        for k in 0..n {
            let joint = joints[level][k];
            if plan.touched(level, k) {
                let side = blend_side(level, k);
                let curve = end_curve(
                    plan,
                    &chart,
                    level,
                    k,
                    corners[k],
                    caps[level][k][side],
                    tolerance,
                )?;
                ends[level][k] = Some(add(curve, [joint.low, joint.rim[side]], &vertices)?);
                if joint.step() {
                    let (cap, corner) = (joint.rim[side], joint.rim[1 - side]);
                    let corner_point = rim_point(level, k, 1 - side);
                    risers[level][k] = Some(add(
                        line(low_point(level, k), corner_point),
                        [joint.low, corner],
                        &vertices,
                    )?);
                    steps[level][k] = Some(add(
                        line(rim_point(level, k, side), corner_point),
                        [cap, corner],
                        &vertices,
                    )?);
                }
            }
            let next = plan.next[k];
            if !plan.collapsed(level, k) {
                let (from, to) = (rim_point(level, k, 1), rim_point(level, next, 0));
                let curve = carried(k, plan.inset_of(level, k), z, from, to);
                rims[level][k] = Some(add(
                    curve,
                    [joint.rim[1], joints[level][next].rim[0]],
                    &vertices,
                )?);
            }
            if plan.flag(level, k) {
                let (from, to) = (low_point(level, k), low_point(level, next));
                let curve = carried(k, 0., z + inward * size, from, to);
                tangents[level][k] =
                    Some(add(curve, [joint.low, joints[level][next].low], &vertices)?);
            }
        }
    }
    // The edge an untreated side face takes up to the cap beside a treated
    // neighbour: the riser of a step, otherwise the blend's end curve.
    let side_end = |level: usize, k: usize| -> usize {
        risers[level][k]
            .or(ends[level][k])
            .expect("side end at a treated vertex")
    };

    let mut faces = Vec::with_capacity(2 + 4 * n);
    let mut names = Vec::with_capacity(faces.capacity());
    let plane = |origin: P3, x: P3, y: P3| Surface::Plane(chart.frame(origin, x, y));
    // Caps: each ring keeps the material on its left, so its rims run
    // counterclockwise about the outward normal of the top.
    for level in 0..2 {
        let (z, _) = plan.level(level);
        let loops = plan
            .loops
            .iter()
            .map(|range| {
                let mut uses = Vec::new();
                for k in range.clone() {
                    if let Some(step) = steps[level][k] {
                        uses.push((step, plan.flag(level, plan.prev[k])));
                    }
                    if let Some(rim) = rims[level][k] {
                        uses.push((rim, true));
                    }
                }
                if level == 0 {
                    uses.reverse();
                    uses.iter_mut().for_each(|use_| use_.1 = !use_.1);
                }
                Loop { uses, vertex: None }
            })
            .collect();
        let y = if level == 1 {
            [0., 1., 0.]
        } else {
            [0., -1., 0.]
        };
        faces.push(Face {
            surface: plane([0., 0., z], [1., 0., 0.], y),
            same_sense: true,
            loops,
        });
        names.push(plan.cap_names[level].clone());
    }
    // Sides.
    for k in 0..n {
        let next = plan.next[k];
        let previous = plan.prev[k];
        let mut uses = Vec::new();
        if let Some(edge) = tangents[0][k] {
            uses.push((edge, true));
        } else {
            if plan.flag(0, previous) {
                uses.push((side_end(0, k), true));
            }
            uses.push((rims[0][k].expect("bottom rim"), true));
            if plan.flag(0, next) {
                uses.push((side_end(0, next), false));
            }
        }
        uses.push((vertical[next], true));
        if let Some(edge) = tangents[1][k] {
            uses.push((edge, false));
        } else {
            if plan.flag(1, next) {
                uses.push((side_end(1, next), true));
            }
            uses.push((rims[1][k].expect("top rim"), false));
            if plan.flag(1, previous) {
                uses.push((side_end(1, k), false));
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
    // Steps: the flat end of a blend at a tangent vertex, perpendicular to
    // the edges there. Its outward normal points along the treated edge.
    for level in 0..2 {
        for k in 0..n {
            let (Some(end), Some(riser), Some(step)) =
                (ends[level][k], risers[level][k], steps[level][k])
            else {
                continue;
            };
            let v = plan.seg(k).start();
            let inward = flat(plan.seg(k).inward(v));
            let ahead = plan.flag(level, k);
            let uses = if ahead == (level == 1) {
                vec![(end, true), (step, true), (riser, false)]
            } else {
                vec![(riser, true), (step, false), (end, false)]
            };
            faces.push(Face {
                surface: plane(plan.at(level, v, 0.), inward, [0., 0., 1.]),
                same_sense: ahead,
                loops: vec![Loop { uses, vertex: None }],
            });
            let (treated, kept) = if ahead {
                (k, plan.prev[k])
            } else {
                (plan.prev[k], k)
            };
            let mut origins = plan.caps[level][treated].clone().expect("treated edge");
            origins.extend(plan.side_names[kept].origins.iter().cloned());
            names.push(FaceName {
                origins,
                support: None,
                anchor: None,
            });
        }
    }
    // Blends.
    for (level, level_tangents) in tangents.iter().enumerate() {
        let (z, inward) = plan.level(level);
        let up = [0., 0., -inward];
        for k in 0..n {
            let Some(origins) = &plan.caps[level][k] else {
                continue;
            };
            let next = plan.next[k];
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
        let joints: Vec<[P2; 2]> = (0..n)
            .map(|k| {
                profile::joint_ends(
                    plan.previous(k),
                    plan.seg(k),
                    corners[k],
                    [deltas[plan.prev[k]], deltas[k]],
                )
                .expect("checked corner")
            })
            .collect();
        let starts: Vec<P2> = joints.iter().map(|joint| joint[1]).collect();
        let ends: Vec<P2> = (0..n).map(|k| joints[plan.next[k]][0]).collect();
        profile::twice_area(&plan.ring, &deltas, &starts, &ends, &plan.next) / 2.
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
    let size = plan.size;
    let (a, b) = (plan.previous(k), plan.seg(k));
    let (sa, sb) = (plan.flag(level, plan.prev[k]), plan.flag(level, k));
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
