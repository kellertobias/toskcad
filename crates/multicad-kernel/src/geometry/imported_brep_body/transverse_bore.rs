//! A profile prism cut by a coaxial stepped tool whose axis crosses the
//! extrusion at a right angle: drilled, counterbored and double-entry access
//! bores across interrupted profile stock, and cross-bored round barrels.
//!
//! In the chart of the cut, `x` runs along the extrusion over `[x0, x1]`,
//! the tool axis is the line `x = y = 0` along `z`, and the profile is a set
//! of line and arc rings in the `(y, z)` plane. The tool removes
//! `x^2 + y^2 < r(z)^2`, where the radius `r` is piecewise constant over
//! axial bands; a band that ends inside the profile leaves a planar shoulder
//! or floor. Every point `(y, z)` of the profile inside a band and with
//! `|y| < r` loses the chord `|x| < h = sqrt(r^2 - y^2)`, so the whole
//! boundary is the lift of a planar arrangement:
//!
//! - profile edges inside the tool lift to the window curves of their side
//!   faces on both branches `x = +-h`: exact lines, circles and ellipses on
//!   planar sides, and on round sides the quartic of two perpendicular
//!   cylinders, carried like a weldment cope as a C1 cubic Hermite curve
//!   through exact points and tangents whose sampled deviation from both
//!   surfaces stays below a sixteenth of the tolerance;
//! - band ends lift to circular arcs on their shoulder planes, and profile
//!   crossings of a band end to straight shoulder segments along `x`;
//! - the silhouette `|y| = r` of a band is a smooth fold of its tool wall
//!   and no edge; a profile vertex exactly on it pinches its side face.
//!
//! The faces are traced from the edges each surface carries, oriented by the
//! exact material indicator, and validated and integrated by the retained
//! boundary. The integrated volume must reproduce an independent quadrature
//! of the prism less the removed chords. Tangency of a profile edge to a
//! silhouette or band end, a profile vertex on a band end inside the tool
//! and a tool reaching the prism ends are refused.
use super::*;
use crate::geometry::authored_sweep::{AuthoredNaming, FaceName};
use brep::{Boundary, Frame, Loop, NurbsCurve, Shell};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

type P = [f64; 3];
type Q2 = [f64; 2];

const OPERATION: &str = "transverse bore";
/// Relative margin every admitted classification must clear.
const RESOLUTION: f64 = 1. / (1_u64 << 30) as f64;
const SAMPLES: usize = 8;
const MAX_SEGMENTS: usize = 1 << 12;
const MAX_EDGES: usize = 20_000;

fn refused(message: impl Into<String>) -> GeometryError {
    error(OPERATION, GeometryErrorKind::Unsupported, message)
}

fn invalid(message: impl Into<String>) -> GeometryError {
    error(OPERATION, GeometryErrorKind::Input, message)
}

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
fn sub(a: P, b: P) -> P {
    std::array::from_fn(|i| a[i] - b[i])
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

/// One axial span of the tool and its radius, measured along the tool axis
/// from the tool origin. Overlapping steps remove their union.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoreStep {
    pub start_mm: f64,
    pub end_mm: f64,
    pub radius_mm: f64,
}

/// A coaxial stepped cutter and the authored names of its faces.
#[derive(Clone, Debug, PartialEq)]
pub struct TransverseBoreTool {
    pub origin_mm: P,
    pub axis: P,
    pub steps: Vec<BoreStep>,
    /// Name of every retained tool wall.
    pub wall: Origin,
    /// Name of every retained shoulder or floor.
    pub shoulder: Origin,
}

/// The original prism and the tool of a transverse bore, kept with its
/// retained boundary for later proofs about the cut material.
pub(crate) struct TransverseRecord {
    pub(crate) source: Solid,
    pub(crate) origin_mm: P,
    pub(crate) axis: P,
    /// `(start, end, radius)` along the axis from the origin.
    pub(crate) bands: Vec<[f64; 3]>,
}

impl TransverseRecord {
    pub(super) fn transformed(&self, transform: RigidTransform) -> Result<Self, GeometryError> {
        let turned = transform.vector(self.origin_mm);
        Ok(Self {
            source: self.source.transformed(transform)?,
            origin_mm: std::array::from_fn(|i| turned[i] + transform.translation_mm[i]),
            axis: transform.vector(self.axis),
            bands: self.bands.clone(),
        })
    }
}

/// The prism to cut: a profile in a placed plane, swept along the plane
/// normal over `span`, with the names of its caps and sides.
pub(crate) struct PrismSource {
    pub(crate) frame: PlaneFrame,
    pub(crate) span: [f64; 2],
    /// Outer ring counterclockwise about the normal, holes clockwise.
    pub(crate) rings: Vec<Vec<multicad_rust_kernel::normal_sweep::Edge>>,
    pub(crate) caps: [Vec<Origin>; 2],
    /// Names of each side, in ring order.
    pub(crate) sides: Vec<Vec<Origin>>,
    pub(crate) feature: EntityId,
}

/// A profile edge in the `(y, z)` chart plane.
#[derive(Clone, Copy, Debug)]
enum Seg {
    Line {
        a: Q2,
        b: Q2,
    },
    /// `c + r (cos(start + f sweep), sin(start + f sweep))`; a full circle
    /// has `|sweep| = TAU` and `a == b`.
    Arc {
        c: Q2,
        r: f64,
        a: Q2,
        b: Q2,
        start: f64,
        sweep: f64,
    },
}

impl Seg {
    fn ends(&self) -> [Q2; 2] {
        match *self {
            Seg::Line { a, b } | Seg::Arc { a, b, .. } => [a, b],
        }
    }

    fn full(&self) -> bool {
        matches!(*self, Seg::Arc { sweep, .. } if sweep.abs() == TAU)
    }

    fn angle(&self, f: f64) -> f64 {
        match *self {
            Seg::Arc { start, sweep, .. } => start + f * sweep,
            Seg::Line { .. } => 0.,
        }
    }

    fn point(&self, f: f64) -> Q2 {
        let [a, b] = self.ends();
        if f == 0. {
            return a;
        }
        if f == 1. && !self.full() {
            return b;
        }
        match *self {
            Seg::Line { a, b } => [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])],
            Seg::Arc { c, r, .. } => {
                let (s, co) = self.angle(f).sin_cos();
                [c[0] + r * co, c[1] + r * s]
            }
        }
    }

    /// Fraction of an arc angle, if the angle lies on the arc.
    fn fraction(&self, angle: f64) -> Option<f64> {
        let Seg::Arc { start, sweep, .. } = *self else {
            return None;
        };
        let t = ((angle - start) * sweep.signum()).rem_euclid(TAU) / sweep.abs();
        (t <= 1.).then_some(t)
    }

    fn signed_area(&self) -> f64 {
        match *self {
            Seg::Line { a, b } => 0.5 * (a[0] * b[1] - b[0] * a[1]),
            Seg::Arc {
                c, r, a, b, sweep, ..
            } => 0.5 * (r * r * sweep + c[0] * (b[1] - a[1]) - c[1] * (b[0] - a[0])),
        }
    }

    fn reversed(self) -> Self {
        match self {
            Seg::Line { a, b } => Seg::Line { a: b, b: a },
            Seg::Arc {
                c,
                r,
                a,
                b,
                start,
                sweep,
            } => Seg::Arc {
                c,
                r,
                a: b,
                b: a,
                start: start + sweep,
                sweep: -sweep,
            },
        }
    }

    /// Crossings of the vertical ray from `(y, z)` toward `+z`.
    fn ray_crossings(&self, y: f64, z: f64) -> usize {
        match *self {
            Seg::Line { a, b } => {
                if (a[0] > y) == (b[0] > y) {
                    return 0;
                }
                let zc = a[1] + (y - a[0]) * (b[1] - a[1]) / (b[0] - a[0]);
                usize::from(zc > z)
            }
            Seg::Arc { c, r, .. } => {
                let dy = y - c[0];
                if dy.abs() >= r {
                    return 0;
                }
                let root = (r * r - dy * dy).sqrt();
                [root, -root]
                    .into_iter()
                    .filter(|&dz| {
                        c[1] + dz > z
                            && self
                                .fraction(dz.atan2(dy))
                                .is_some_and(|t| t > 0. && t < 1.)
                    })
                    .count()
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Band {
    z: [f64; 2],
    r: f64,
}

#[derive(Clone, Copy, Debug)]
struct Level {
    z: f64,
    below: Option<usize>,
    above: Option<usize>,
}

impl Level {
    fn bands(&self) -> impl Iterator<Item = usize> + '_ {
        self.below.into_iter().chain(self.above)
    }
}

/// A point of the planar arrangement.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Pid {
    Vertex(usize, usize),
    Split(usize, usize, usize),
    /// `(level, band, upper side)`: where a band end meets its silhouette.
    Corner(usize, usize, bool),
}

/// Where a lifted vertex lies along `x`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Xk {
    Low,
    High,
    Zero,
    Pos(usize),
    Neg(usize),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Surf {
    Cap(bool),
    Side(usize),
    Wall(usize),
    Shoulder(usize),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Family {
    Frame,
    Window,
}

struct EdgeRec {
    edge: brep::Edge,
    surfaces: [Surf; 2],
    family: Family,
}

#[derive(Clone, Copy, Debug)]
enum Kind {
    /// A profile vertex or split point outside every band.
    Free,
    Inside(usize),
    Silhouette(usize),
    Level(usize),
}

#[derive(Clone, Copy, Debug)]
struct Point2 {
    yz: Q2,
    kind: Kind,
}

/// The planar arrangement and its lift.
struct Cut {
    rings: Vec<Vec<Seg>>,
    /// Flattened side index of each ring edge.
    side: Vec<Vec<usize>>,
    x: [f64; 2],
    bands: Vec<Band>,
    levels: Vec<Level>,
    margin: f64,
    snap: f64,
    tolerance: f64,
    points: BTreeMap<Pid, Point2>,
    /// Split fractions of each ring edge, sorted, without the ends.
    splits: Vec<Vec<Vec<f64>>>,
    vertices: Vec<P>,
    lifted: BTreeMap<(Pid, Xk), usize>,
    edges: Vec<EdgeRec>,
}

fn classify(value: f64, snap: f64, margin: f64) -> Result<std::cmp::Ordering, GeometryError> {
    if value.abs() <= snap {
        return Ok(std::cmp::Ordering::Equal);
    }
    if value.abs() < margin {
        return Err(refused(
            "Transverse bore classification is unresolved within the certified margin",
        ));
    }
    Ok(value.total_cmp(&0.))
}

impl Cut {
    fn band_at(&self, z: f64) -> Option<usize> {
        self.bands.iter().position(|b| b.z[0] < z && z < b.z[1])
    }

    fn radius_at(&self, z: f64) -> f64 {
        self.band_at(z).map_or(0., |k| self.bands[k].r)
    }

    fn in_profile(&self, y: f64, z: f64) -> bool {
        // A generic offset far below the margin keeps the ray off profile
        // vertices that lie on the chart axes.
        let y = y + 0.3713 * self.margin;
        self.rings
            .iter()
            .flatten()
            .map(|seg| seg.ray_crossings(y, z))
            .sum::<usize>()
            % 2
            == 1
    }

    /// The exact material indicator of the cut prism.
    fn inside(&self, [x, y, z]: P) -> bool {
        let r = self.radius_at(z);
        x > self.x[0] && x < self.x[1] && self.in_profile(y, z) && x * x + y * y >= r * r
    }

    fn level_radius(&self, level: &Level) -> [f64; 2] {
        let radii: Vec<f64> = level.bands().map(|k| self.bands[k].r).collect();
        [
            radii.iter().copied().fold(f64::INFINITY, f64::min),
            radii.iter().copied().fold(0., f64::max),
        ]
    }

    /// Classify a profile vertex against every band and band end.
    fn vertex_kind(&self, [y, z]: Q2) -> Result<Kind, GeometryError> {
        for level in &self.levels {
            if classify(z - level.z, self.snap, self.margin)? == std::cmp::Ordering::Equal {
                let [_, outer] = self.level_radius(level);
                if classify(y.abs() - outer, self.snap, self.margin)? == std::cmp::Ordering::Less {
                    return Err(refused(
                        "A profile vertex lies on a tool shoulder inside the tool",
                    ));
                }
                return Ok(Kind::Free);
            }
        }
        let Some(k) = self.band_at(z) else {
            return Ok(Kind::Free);
        };
        Ok(
            match classify(y.abs() - self.bands[k].r, self.snap, self.margin)? {
                std::cmp::Ordering::Less => Kind::Inside(k),
                std::cmp::Ordering::Equal => Kind::Silhouette(k),
                std::cmp::Ordering::Greater => Kind::Free,
            },
        )
    }

    /// Interior fractions of `seg` where coordinate `axis` equals `value`,
    /// with refusals for tangency that `relevant` accepts.
    fn crossings(
        &self,
        seg: &Seg,
        axis: usize,
        value: f64,
        relevant: &dyn Fn(Q2) -> Result<bool, GeometryError>,
    ) -> Result<Vec<f64>, GeometryError> {
        let [a, b] = seg.ends();
        let mut found = Vec::new();
        match *seg {
            Seg::Line { .. } => {
                let p = classify(a[axis] - value, self.snap, self.margin)?;
                let q = classify(b[axis] - value, self.snap, self.margin)?;
                use std::cmp::Ordering::Equal;
                if p == Equal && q == Equal {
                    // Lying on the line: the relevant part must be empty.
                    for f in [0.25, 0.5, 0.75] {
                        if relevant(seg.point(f))? {
                            return Err(refused(
                                "A profile edge lies on a tool silhouette or shoulder",
                            ));
                        }
                    }
                } else if p != Equal && q != Equal && p != q {
                    let f = (a[axis] - value) / (a[axis] - b[axis]);
                    if relevant(seg.point(f))? {
                        found.push(f);
                    }
                }
            }
            Seg::Arc { c, r, .. } => {
                let offset = value - c[axis];
                let gap = classify(offset.abs() - r, self.snap, self.margin)?;
                if gap == std::cmp::Ordering::Equal {
                    let angle = if axis == 0 {
                        if offset > 0. { 0. } else { PI }
                    } else if offset > 0. {
                        FRAC_PI_2
                    } else {
                        -FRAC_PI_2
                    };
                    let mut point = c;
                    point[axis] = value;
                    if seg.fraction(angle).is_some() && relevant(point)? {
                        return Err(refused(
                            "A profile arc is tangent to a tool silhouette or shoulder",
                        ));
                    }
                } else if gap == std::cmp::Ordering::Less {
                    let angles = if axis == 0 {
                        let base = (offset / r).clamp(-1., 1.).acos();
                        [base, -base]
                    } else {
                        let base = (offset / r).clamp(-1., 1.).asin();
                        [base, PI - base]
                    };
                    for angle in angles {
                        let Some(f) = seg.fraction(angle) else {
                            continue;
                        };
                        let mut point = seg.point(f);
                        point[axis] = value;
                        let at_end = [a, b].iter().any(|end| {
                            (end[0] - point[0]).abs() <= self.margin
                                && (end[1] - point[1]).abs() <= self.margin
                        });
                        if at_end || f <= 0. || f >= 1. {
                            continue;
                        }
                        if relevant(point)? {
                            found.push(f);
                        }
                    }
                }
            }
        }
        Ok(found)
    }

    /// Split every profile edge at silhouettes, shoulders and the `y`
    /// extremes of arcs inside the tool, recording each arrangement point.
    fn split(&mut self) -> Result<(), GeometryError> {
        let mut splits = Vec::new();
        for (ring, segs) in self.rings.clone().iter().enumerate() {
            let mut ring_splits = Vec::new();
            for (index, seg) in segs.iter().enumerate() {
                if !seg.full() {
                    let kind = self.vertex_kind(seg.ends()[0])?;
                    self.points.insert(
                        Pid::Vertex(ring, index),
                        Point2 {
                            yz: seg.ends()[0],
                            kind,
                        },
                    );
                }
                let mut found: Vec<(f64, Kind)> = Vec::new();
                for (k, band) in self.bands.iter().enumerate() {
                    for sign in [1., -1.] {
                        let relevant = |p: Q2| -> Result<bool, GeometryError> {
                            let low = classify(p[1] - band.z[0], self.snap, self.margin)?;
                            let high = classify(p[1] - band.z[1], self.snap, self.margin)?;
                            use std::cmp::Ordering::*;
                            if low == Equal || high == Equal {
                                return Err(refused(
                                    "A profile edge crosses a tool silhouette at a shoulder",
                                ));
                            }
                            Ok(low == Greater && high == Less)
                        };
                        for f in self.crossings(seg, 0, sign * band.r, &relevant)? {
                            found.push((f, Kind::Silhouette(k)));
                        }
                    }
                }
                for (s, level) in self.levels.iter().enumerate() {
                    let [inner, outer] = self.level_radius(level);
                    let relevant = |p: Q2| -> Result<bool, GeometryError> {
                        let y = p[0].abs();
                        let to_outer = classify(y - outer, self.snap, self.margin)?;
                        if to_outer == std::cmp::Ordering::Equal
                            || (level.below.is_some()
                                && level.above.is_some()
                                && classify(y - inner, self.snap, self.margin)?
                                    == std::cmp::Ordering::Equal)
                        {
                            return Err(refused(
                                "A profile edge crosses a tool shoulder at a silhouette",
                            ));
                        }
                        Ok(to_outer == std::cmp::Ordering::Less)
                    };
                    for f in self.crossings(seg, 1, level.z, &relevant)? {
                        found.push((f, Kind::Level(s)));
                    }
                }
                if let Seg::Arc { c, r, .. } = *seg {
                    for angle in [0., PI] {
                        let Some(f) = seg.fraction(angle) else {
                            continue;
                        };
                        if f <= 0. || f >= 1. {
                            continue;
                        }
                        let p = [c[0] + r * angle.cos(), c[1]];
                        if let Some(k) = self.band_at(p[1])
                            && p[0].abs() < self.bands[k].r
                            && found.iter().all(|(g, _)| (g - f).abs() > 1e-9)
                        {
                            found.push((f, Kind::Inside(k)));
                        }
                    }
                }
                found.sort_by(|a, b| a.0.total_cmp(&b.0));
                if found.windows(2).any(|w| w[1].0 - w[0].0 <= 1e-12) {
                    return Err(refused(
                        "Tool silhouettes and shoulders meet a profile edge at one point",
                    ));
                }
                for (i, &(f, kind)) in found.iter().enumerate() {
                    let mut yz = seg.point(f);
                    match kind {
                        Kind::Silhouette(k) => yz[0] = self.bands[k].r.copysign(yz[0]),
                        Kind::Level(s) => yz[1] = self.levels[s].z,
                        _ => {}
                    }
                    self.points
                        .insert(Pid::Split(ring, index, i), Point2 { yz, kind });
                }
                ring_splits.push(found.iter().map(|(f, _)| *f).collect());
            }
            splits.push(ring_splits);
        }
        self.splits = splits;
        Ok(())
    }

    /// Points of edge `(ring, index)` in order, with their fractions.
    fn edge_points(&self, ring: usize, index: usize) -> Vec<(Pid, f64)> {
        let seg = &self.rings[ring][index];
        let count = self.rings[ring].len();
        let mut points = Vec::new();
        if !seg.full() {
            points.push((Pid::Vertex(ring, index), 0.));
        }
        for (i, &f) in self.splits[ring][index].iter().enumerate() {
            points.push((Pid::Split(ring, index, i), f));
        }
        if seg.full() {
            if points.is_empty() {
                // A closed ring without splits starts and ends at its start.
                return vec![];
            }
            let first = points[0];
            points.push((first.0, first.1 + 1.));
        } else {
            points.push((Pid::Vertex(ring, (index + 1) % count), 1.));
        }
        points
    }

    /// The band of the arrangement piece between two fractions, if the
    /// tool removes its chord.
    fn piece_band(&self, seg: &Seg, f0: f64, f1: f64) -> Option<usize> {
        let [y, z] = seg.point(0.5 * (f0 + f1));
        let k = self.band_at(z)?;
        (y.abs() < self.bands[k].r).then_some(k)
    }

    fn half(&self, k: usize, y: f64) -> f64 {
        let r = self.bands[k].r;
        ((r - y) * (r + y)).max(0.).sqrt()
    }

    fn on_silhouette(&self, pid: Pid, k: usize) -> bool {
        match pid {
            Pid::Corner(_, band, _) => band == k,
            _ => matches!(self.points[&pid].kind, Kind::Silhouette(b) if b == k),
        }
    }

    fn yz(&self, pid: Pid) -> Q2 {
        match pid {
            Pid::Corner(level, band, upper) => [
                if upper {
                    self.bands[band].r
                } else {
                    -self.bands[band].r
                },
                self.levels[level].z,
            ],
            _ => self.points[&pid].yz,
        }
    }

    fn vertex(&mut self, pid: Pid, xk: Xk) -> usize {
        let xk = match xk {
            Xk::Pos(k) | Xk::Neg(k) if self.on_silhouette(pid, k) => Xk::Zero,
            other => other,
        };
        if let Some(&v) = self.lifted.get(&(pid, xk)) {
            return v;
        }
        let [y, z] = self.yz(pid);
        let x = match xk {
            Xk::Low => self.x[0],
            Xk::High => self.x[1],
            Xk::Zero => 0.,
            Xk::Pos(k) => self.half(k, y),
            Xk::Neg(k) => -self.half(k, y),
        };
        self.vertices.push([x, y, z]);
        self.lifted.insert((pid, xk), self.vertices.len() - 1);
        self.vertices.len() - 1
    }

    fn add(
        &mut self,
        curve: Curve,
        ends: [usize; 2],
        surfaces: [Surf; 2],
        family: Family,
    ) -> Result<(), GeometryError> {
        if self.edges.len() >= MAX_EDGES {
            return Err(failure(brep::Error::Limit));
        }
        let points = ends.map(|v| self.vertices[v]);
        let edge = brep::Edge::new(curve, ends, points, self.tolerance).map_err(failure)?;
        self.edges.push(EdgeRec {
            edge,
            surfaces,
            family,
        });
        Ok(())
    }

    /// A line between two lifted vertices.
    fn line(&mut self, ends: [usize; 2], surfaces: [Surf; 2]) -> Result<(), GeometryError> {
        let [a, b] = ends.map(|v| self.vertices[v]);
        self.add(
            Curve::Line {
                origin: a,
                direction: sub(b, a),
            },
            ends,
            surfaces,
            Family::Frame,
        )
    }

    /// Order the ends of a periodic `theta` curve so that its forward range
    /// passes through `middle`.
    fn periodic_ends(&self, ends: [usize; 2], middle: P) -> [usize; 2] {
        let angle = |p: P| p[1].atan2(p[0]);
        let [a, b] = ends.map(|v| angle(self.vertices[v]));
        let m = angle(middle);
        if (m - a).rem_euclid(TAU) < (b - a).rem_euclid(TAU) {
            ends
        } else {
            [ends[1], ends[0]]
        }
    }

    fn cap_edges(&mut self) -> Result<(), GeometryError> {
        for ring in 0..self.rings.len() {
            for index in 0..self.rings[ring].len() {
                let seg = self.rings[ring][index];
                let side = Surf::Side(self.side[ring][index]);
                for (end, xk) in [(false, Xk::Low), (true, Xk::High)] {
                    let x = self.x[usize::from(end)];
                    let start = if seg.full() {
                        // A closed ring keeps one vertex at its start.
                        self.points
                            .entry(Pid::Vertex(ring, index))
                            .or_insert(Point2 {
                                yz: seg.ends()[0],
                                kind: Kind::Free,
                            });
                        self.vertex(Pid::Vertex(ring, index), xk)
                    } else {
                        self.vertex(Pid::Vertex(ring, index), xk)
                    };
                    let finish = if seg.full() {
                        start
                    } else {
                        let next = (index + 1) % self.rings[ring].len();
                        self.vertex(Pid::Vertex(ring, next), xk)
                    };
                    let surfaces = [Surf::Cap(end), side];
                    match seg {
                        Seg::Line { .. } => self.line([start, finish], surfaces)?,
                        Seg::Arc { c, r, sweep, .. } => {
                            let s = sweep.signum();
                            self.add(
                                Curve::Circle {
                                    frame: Frame {
                                        origin: [x, c[0], c[1]],
                                        x: [0., 1., 0.],
                                        y: [0., 0., s],
                                        z: [s, 0., 0.],
                                    },
                                    radius: r,
                                },
                                [start, finish],
                                surfaces,
                                Family::Frame,
                            )?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn long_edges(&mut self) -> Result<(), GeometryError> {
        for ring in 0..self.rings.len() {
            let count = self.rings[ring].len();
            if count == 1 && self.rings[ring][0].full() {
                continue;
            }
            for index in 0..count {
                let previous = Surf::Side(self.side[ring][(index + count - 1) % count]);
                let next = Surf::Side(self.side[ring][index]);
                let pid = Pid::Vertex(ring, index);
                let low = self.vertex(pid, Xk::Low);
                let high = self.vertex(pid, Xk::High);
                match self.points[&pid].kind {
                    Kind::Inside(k) | Kind::Silhouette(k) => {
                        let neg = self.vertex(pid, Xk::Neg(k));
                        let pos = self.vertex(pid, Xk::Pos(k));
                        self.line([low, neg], [previous, next])?;
                        self.line([pos, high], [previous, next])?;
                    }
                    _ => self.line([low, high], [previous, next])?,
                }
            }
        }
        Ok(())
    }

    /// The exact point and parameter derivative of a window on a round side.
    fn round_window(
        &self,
        seg: &Seg,
        k: usize,
        branch: f64,
        by_tool: Option<f64>,
        s: f64,
    ) -> (P, P) {
        let Seg::Arc { c, r: radius, .. } = *seg else {
            unreachable!("round window on an arc");
        };
        let r = self.bands[k].r;
        if let Some(lobe) = by_tool {
            // Tool angle: x = r cos, y = r sin, z on the profile circle.
            let (sin, cos) = s.sin_cos();
            let y = r * sin;
            let dy = y - c[0];
            let root = (radius * radius - dy * dy).max(0.).sqrt();
            let z = c[1] + lobe * root;
            let dz = -lobe * dy * r * cos / root;
            ([r * cos, y, z], [-r * sin, r * cos, dz])
        } else {
            // Profile angle: y and z on the circle, x on the tool wall.
            let (sin, cos) = s.sin_cos();
            let y = c[0] + radius * cos;
            let z = c[1] + radius * sin;
            let dy = -radius * sin;
            let root = (r * r - y * y).max(0.).sqrt();
            (
                [branch * root, y, z],
                [-branch * y * dy / root, dy, radius * cos],
            )
        }
    }

    /// Residual of a point from the tool wall and the round side.
    fn round_residual(&self, seg: &Seg, k: usize, p: P) -> f64 {
        let Seg::Arc { c, r: radius, .. } = *seg else {
            return f64::INFINITY;
        };
        let r = self.bands[k].r;
        let tool = (p[0].hypot(p[1]) - r).abs();
        let side = ((p[1] - c[0]).hypot(p[2] - c[1]) - radius).abs();
        tool.max(side)
    }

    #[allow(clippy::too_many_arguments)]
    fn hermite(
        &self,
        seg: &Seg,
        k: usize,
        branch: f64,
        by_tool: Option<f64>,
        range: [f64; 2],
        first: P,
        last: P,
    ) -> Result<NurbsCurve, GeometryError> {
        let mut segments = 4;
        loop {
            let h = (range[1] - range[0]) / segments as f64;
            let nodes: Vec<(P, P)> = (0..=segments)
                .map(|i| self.round_window(seg, k, branch, by_tool, range[0] + h * i as f64))
                .collect();
            let mut points = Vec::with_capacity(3 * segments + 1);
            let mut worst: f64 = 0.;
            for i in 0..segments {
                let p0 = if i == 0 { first } else { nodes[i].0 };
                let p1 = if i + 1 == segments {
                    last
                } else {
                    nodes[i + 1].0
                };
                let (d0, d1) = (nodes[i].1, nodes[i + 1].1);
                let control = [p0, mix(p0, h / 3., d0), mix(p1, -h / 3., d1), p1];
                if i == 0 {
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
                        std::array::from_fn(|n| (0..4).map(|m| w[m] * control[m][n]).sum());
                    let exact = self
                        .round_window(seg, k, branch, by_tool, range[0] + h * (i as f64 + u))
                        .0;
                    worst = worst
                        .max(norm(sub(bezier, exact)))
                        .max(self.round_residual(seg, k, bezier));
                }
            }
            if worst.is_finite() && worst <= self.tolerance / 16. {
                let distinct = (0..=segments).map(|i| range[0] + h * i as f64).collect();
                let mut multiplicities = vec![3; segments + 1];
                multiplicities[0] = 4;
                multiplicities[segments] = 4;
                return NurbsCurve::new(3, points, None, distinct, multiplicities).map_err(failure);
            }
            segments *= 2;
            if segments > MAX_SEGMENTS {
                return Err(refused(
                    "A transverse bore window is too close to tangency for an exact Rust boundary",
                ));
            }
        }
    }

    /// The window curve of one inside piece on one branch.
    #[allow(clippy::too_many_arguments)]
    fn window(
        &mut self,
        seg: Seg,
        k: usize,
        branch: f64,
        f: [f64; 2],
        pids: [Pid; 2],
        side: usize,
    ) -> Result<(), GeometryError> {
        let xk = if branch > 0. { Xk::Pos(k) } else { Xk::Neg(k) };
        let ends = pids.map(|pid| self.vertex(pid, xk));
        let surfaces = [Surf::Side(side), Surf::Wall(k)];
        let r = self.bands[k].r;
        let middle2 = seg.point(0.5 * (f[0] + f[1]));
        let middle = [branch * self.half(k, middle2[0]), middle2[0], middle2[1]];
        match seg {
            Seg::Line { a, b } => {
                let (dy, dz) = (b[0] - a[0], b[1] - a[1]);
                if dy.abs() <= self.snap {
                    self.add(
                        Curve::Line {
                            origin: self.vertices[ends[0]],
                            direction: sub(self.vertices[ends[1]], self.vertices[ends[0]]),
                        },
                        ends,
                        surfaces,
                        Family::Window,
                    )
                } else if dz.abs() <= self.snap {
                    let ordered = self.periodic_ends(ends, middle);
                    self.add(
                        Curve::Circle {
                            frame: Frame {
                                origin: [0., 0., a[1]],
                                x: [1., 0., 0.],
                                y: [0., 1., 0.],
                                z: [0., 0., 1.],
                            },
                            radius: r,
                        },
                        ordered,
                        surfaces,
                        Family::Window,
                    )
                } else {
                    let m = dz / dy;
                    let stretch = m.hypot(1.);
                    let y = [0., 1. / stretch, m / stretch];
                    let ordered = self.periodic_ends(ends, middle);
                    self.add(
                        Curve::Ellipse {
                            frame: Frame {
                                origin: [0., 0., a[1] - a[0] * m],
                                x: [1., 0., 0.],
                                y,
                                z: cross([1., 0., 0.], y),
                            },
                            semiaxes: [r, r * stretch],
                        },
                        ordered,
                        surfaces,
                        Family::Window,
                    )
                }
            }
            Seg::Arc { c, .. } => {
                let touches = pids.iter().any(|&pid| self.on_silhouette(pid, k));
                let (by_tool, range) = if touches {
                    // The quartic is a graph over the tool angle on a piece
                    // without a y extreme of the profile circle.
                    let lobe = (middle2[1] - c[1]).signum();
                    let angle = |p: P| {
                        let t = p[1].atan2(p[0]);
                        if branch < 0. { t.rem_euclid(TAU) } else { t }
                    };
                    (Some(lobe), ends.map(|v| angle(self.vertices[v])))
                } else {
                    (None, f.map(|g| seg.angle(g)))
                };
                let mut range = range;
                let mut points = ends.map(|v| self.vertices[v]);
                let mut ordered = ends;
                if range[1] < range[0] {
                    range.swap(0, 1);
                    points.swap(0, 1);
                    ordered.swap(0, 1);
                }
                let curve = self.hermite(&seg, k, branch, by_tool, range, points[0], points[1])?;
                self.add(
                    Curve::nurbs(curve, self.tolerance),
                    ordered,
                    surfaces,
                    Family::Window,
                )
            }
        }
    }

    fn window_edges(&mut self) -> Result<(), GeometryError> {
        for ring in 0..self.rings.len() {
            for index in 0..self.rings[ring].len() {
                let seg = self.rings[ring][index];
                let side = self.side[ring][index];
                let points = self.edge_points(ring, index);
                if points.is_empty() {
                    // A closed ring without arrangement points.
                    let [y, z] = seg.point(0.);
                    if let Some(k) = self.band_at(z)
                        && y.abs() < self.bands[k].r
                    {
                        let pid = Pid::Vertex(ring, index);
                        for branch in [1., -1.] {
                            self.window(seg, k, branch, [0., 1.], [pid, pid], side)?;
                        }
                    }
                    continue;
                }
                for pair in points.windows(2) {
                    let [(p0, f0), (p1, f1)] = [pair[0], pair[1]];
                    let Some(k) = self.piece_band(&seg, f0, f1) else {
                        continue;
                    };
                    for branch in [1., -1.] {
                        self.window(seg, k, branch, [f0, f1], [p0, p1], side)?;
                    }
                }
            }
        }
        Ok(())
    }

    fn step_edges(&mut self) -> Result<(), GeometryError> {
        let pids: Vec<(Pid, usize, usize)> = self
            .points
            .iter()
            .filter_map(|(&pid, point)| match (pid, point.kind) {
                (Pid::Split(ring, index, _), Kind::Level(s)) => {
                    Some((pid, s, self.side[ring][index]))
                }
                _ => None,
            })
            .collect();
        for (pid, s, side) in pids {
            let level = self.levels[s];
            let y = self.yz(pid)[0].abs();
            let mut bands: Vec<usize> = level.bands().collect();
            bands.sort_by(|a, b| self.bands[*a].r.total_cmp(&self.bands[*b].r));
            let outer = *bands.last().expect("a level bounds a band");
            let inner = (bands.len() == 2 && y < self.bands[bands[0]].r).then_some(bands[0]);
            let surfaces = [Surf::Side(side), Surf::Shoulder(s)];
            match inner {
                Some(inner) => {
                    let ends = [
                        self.vertex(pid, Xk::Pos(inner)),
                        self.vertex(pid, Xk::Pos(outer)),
                    ];
                    self.line(ends, surfaces)?;
                    let ends = [
                        self.vertex(pid, Xk::Neg(outer)),
                        self.vertex(pid, Xk::Neg(inner)),
                    ];
                    self.line(ends, surfaces)?;
                }
                None => {
                    let ends = [
                        self.vertex(pid, Xk::Neg(outer)),
                        self.vertex(pid, Xk::Pos(outer)),
                    ];
                    self.line(ends, surfaces)?;
                }
            }
        }
        Ok(())
    }

    fn shoulder_edges(&mut self) -> Result<(), GeometryError> {
        for s in 0..self.levels.len() {
            let level = self.levels[s];
            let crossings: Vec<Pid> = self
                .points
                .iter()
                .filter(|(_, p)| matches!(p.kind, Kind::Level(l) if l == s))
                .map(|(&pid, _)| pid)
                .collect();
            for k in level.bands().collect::<Vec<_>>() {
                let r = self.bands[k].r;
                let mut stops: Vec<(f64, Pid)> = crossings
                    .iter()
                    .map(|&pid| (self.yz(pid)[0], pid))
                    .filter(|(y, _)| y.abs() < r)
                    .collect();
                stops.push((-r, Pid::Corner(s, k, false)));
                stops.push((r, Pid::Corner(s, k, true)));
                stops.sort_by(|a, b| a.0.total_cmp(&b.0));
                for pair in stops.windows(2) {
                    let [(ya, a), (yb, b)] = [pair[0], pair[1]];
                    if !self.in_profile(0.5 * (ya + yb), level.z) {
                        continue;
                    }
                    let circle = Curve::Circle {
                        frame: Frame {
                            origin: [0., 0., level.z],
                            x: [1., 0., 0.],
                            y: [0., 1., 0.],
                            z: [0., 0., 1.],
                        },
                        radius: r,
                    };
                    let surfaces = [Surf::Shoulder(s), Surf::Wall(k)];
                    let front = [self.vertex(a, Xk::Pos(k)), self.vertex(b, Xk::Pos(k))];
                    self.add(circle.clone(), front, surfaces, Family::Frame)?;
                    let back = [self.vertex(b, Xk::Neg(k)), self.vertex(a, Xk::Neg(k))];
                    self.add(circle, back, surfaces, Family::Frame)?;
                }
            }
        }
        Ok(())
    }

    fn surface(&self, surf: Surf) -> Surface {
        match surf {
            Surf::Cap(end) => Surface::Plane(Frame {
                origin: [self.x[usize::from(end)], 0., 0.],
                x: [0., 1., 0.],
                y: [0., 0., 1.],
                z: [1., 0., 0.],
            }),
            Surf::Side(side) => {
                let seg = self.seg_of(side);
                match seg {
                    Seg::Line { a, b } => {
                        let d = [0., b[0] - a[0], b[1] - a[1]];
                        let d = unit(d).expect("validated profile edge");
                        Surface::Plane(Frame {
                            origin: [0., a[0], a[1]],
                            x: d,
                            y: [1., 0., 0.],
                            z: cross(d, [1., 0., 0.]),
                        })
                    }
                    Seg::Arc { c, r, .. } => Surface::Cylinder {
                        frame: Frame {
                            origin: [0., c[0], c[1]],
                            x: [0., 1., 0.],
                            y: [0., 0., 1.],
                            z: [1., 0., 0.],
                        },
                        radius: r,
                    },
                }
            }
            Surf::Wall(k) => Surface::Cylinder {
                frame: Frame {
                    origin: [0.; 3],
                    x: [1., 0., 0.],
                    y: [0., 1., 0.],
                    z: [0., 0., 1.],
                },
                radius: self.bands[k].r,
            },
            Surf::Shoulder(s) => Surface::Plane(Frame {
                origin: [0., 0., self.levels[s].z],
                x: [1., 0., 0.],
                y: [0., 1., 0.],
                z: [0., 0., 1.],
            }),
        }
    }

    /// Whether `p` projects into the part of a support surface that can
    /// carry faces of this cut: a side over its profile edge, a wall over
    /// its band and a shoulder between its radii. Tangent profile edges
    /// continue each other's supports beyond their spans.
    fn on_patch(&self, surf: Surf, p: P) -> bool {
        match surf {
            Surf::Cap(_) => true,
            Surf::Side(side) => match self.seg_of(side) {
                Seg::Line { a, b } => {
                    let d = [b[0] - a[0], b[1] - a[1]];
                    let t =
                        ((p[1] - a[0]) * d[0] + (p[2] - a[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1]);
                    t > 0. && t < 1.
                }
                seg @ Seg::Arc { c, .. } => {
                    seg.full()
                        || seg
                            .fraction((p[2] - c[1]).atan2(p[1] - c[0]))
                            .is_some_and(|t| t > 0. && t < 1.)
                }
            },
            Surf::Wall(k) => p[2] > self.bands[k].z[0] && p[2] < self.bands[k].z[1],
            Surf::Shoulder(s) => {
                let level = &self.levels[s];
                let [inner, outer] = self.level_radius(level);
                let inner = if level.below.is_some() && level.above.is_some() {
                    inner
                } else {
                    0.
                };
                let rho = p[0].hypot(p[1]);
                rho > inner && rho < outer
            }
        }
    }

    fn seg_of(&self, side: usize) -> Seg {
        for (ring, sides) in self.side.iter().enumerate() {
            if let Some(index) = sides.iter().position(|&s| s == side) {
                return self.rings[ring][index];
            }
        }
        unreachable!("every side belongs to a ring")
    }

    /// Geometric normal and chart coordinates `(a, b)` on a surface; the
    /// chart basis is right-handed about the normal.
    fn chart(&self, surface: &Surface, p: P) -> (P, Q2) {
        match surface {
            Surface::Plane(frame) => {
                let local = frame.local(p);
                (frame.z, [local[0], local[1]])
            }
            Surface::Cylinder { frame, .. } => {
                let local = frame.local(p);
                let normal =
                    unit(mix(scale(frame.x, local[0]), local[1], frame.y)).unwrap_or(frame.x);
                (normal, [local[1].atan2(local[0]), local[2]])
            }
            _ => unreachable!("transverse bores carry planes and cylinders"),
        }
    }

    /// Faces of one surface from the edges it carries.
    fn faces_of(&self, surf: Surf, step: f64) -> Result<Vec<(brep::Face, Surf)>, GeometryError> {
        let surface = self.surface(surf);
        let periodic = matches!(surface, Surface::Cylinder { .. });
        let mut uses = Vec::new();
        let mut sense = None;
        for (index, record) in self.edges.iter().enumerate() {
            if !record.surfaces.contains(&surf) {
                continue;
            }
            let edge = &record.edge;
            let t = 0.5 * (edge.range[0] + edge.range[1]);
            let (m, tangent) = edge.curve.evaluate(t);
            let (normal, _) = self.chart(&surface, m);
            let left = unit(cross(normal, tangent)).ok_or_else(|| refused("Degenerate edge"))?;
            let mut found = None;
            for side in [1., -1.] {
                let q = mix(m, side * step, left);
                let below = self.inside(mix(q, -step / 8., normal));
                let above = self.inside(mix(q, step / 8., normal));
                if below != above && self.on_patch(surf, q) {
                    if found.is_some() {
                        return Err(refused(format!(
                            "A transverse bore edge bounds faces on both sides of one surface ({surf:?} at {m:?})"
                        )));
                    }
                    found = Some((side, if below { 1. } else { -1. }));
                }
            }
            let (side, sigma) = found.ok_or_else(|| {
                refused(format!(
                    "A transverse bore edge bounds no face of its surface ({surf:?} at {m:?})"
                ))
            })?;
            if sense.is_some_and(|s| s != sigma) {
                return Err(refused(format!(
                    "A transverse bore surface changes its material side ({surf:?} at {m:?})"
                )));
            }
            sense = Some(sigma);
            uses.push((index, side == sigma));
        }
        let Some(sigma) = sense else {
            return Ok(vec![]);
        };
        let ends = |&(edge, forward): &(usize, bool)| {
            let v = self.edges[edge].edge.vertices;
            if forward { [v[0], v[1]] } else { [v[1], v[0]] }
        };
        let mut outgoing: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (i, u) in uses.iter().enumerate() {
            outgoing.entry(ends(u)[0]).or_default().push(i);
        }
        let mut used = vec![false; uses.len()];
        let mut loops: Vec<Vec<(usize, bool)>> = Vec::new();
        while let Some(first) = (0..uses.len())
            .filter(|&i| !used[i])
            .min_by_key(|&i| outgoing[&ends(&uses[i])[0]].len())
        {
            let origin = ends(&uses[first])[0];
            let mut current = first;
            let mut trace = Vec::new();
            loop {
                used[current] = true;
                trace.push(uses[current]);
                let end = ends(&uses[current])[1];
                if end == origin {
                    break;
                }
                let candidates: Vec<usize> = outgoing
                    .get(&end)
                    .map(|c| c.iter().copied().filter(|&i| !used[i]).collect())
                    .unwrap_or_default();
                let next = match candidates.as_slice() {
                    [] => return Err(refused("A transverse bore face boundary is open")),
                    [only] => *only,
                    many => {
                        let family = self.edges[uses[current].0].family;
                        let other: Vec<usize> = many
                            .iter()
                            .copied()
                            .filter(|&i| self.edges[uses[i].0].family != family)
                            .collect();
                        match other.as_slice() {
                            [only] => *only,
                            _ => {
                                return Err(refused(
                                    "A transverse bore face meets itself ambiguously",
                                ));
                            }
                        }
                    }
                };
                current = next;
            }
            loops.push(trace);
        }
        // Classify the loops in the surface chart.
        let period = if periodic { Some(TAU) } else { None };
        struct Traced {
            uses: Vec<(usize, bool)>,
            points: Vec<Q2>,
            winding: i64,
            area: f64,
        }
        let mut traced = Vec::new();
        for trace in loops {
            let mut points: Vec<Q2> = Vec::new();
            for &(edge, forward) in &trace {
                let mut polyline = self.edges[edge]
                    .edge
                    .polyline(self.tolerance.max(1e-4), 0.2)
                    .map_err(failure)?;
                if !forward {
                    polyline.reverse();
                }
                for p in polyline {
                    let (_, mut ab) = self.chart(&surface, p);
                    if let (Some(period), Some(last)) = (period, points.last()) {
                        ab[0] += ((last[0] - ab[0]) / period).round() * period;
                    }
                    points.push(ab);
                }
            }
            let winding = period.map_or(0, |period| {
                ((points.last().unwrap()[0] - points[0][0]) / period).round() as i64
            });
            let area = 0.5
                * points
                    .windows(2)
                    .map(|w| w[0][0] * w[1][1] - w[1][0] * w[0][1])
                    .sum::<f64>();
            traced.push(Traced {
                uses: trace,
                points,
                winding,
                area,
            });
        }
        let contains = |outer: &Traced, point: Q2| {
            let mut point = point;
            if let Some(period) = period {
                let low = outer
                    .points
                    .iter()
                    .map(|p| p[0])
                    .fold(f64::INFINITY, f64::min);
                point[0] += ((low - point[0]) / period).ceil() * period;
            }
            let mut inside = false;
            for w in outer.points.windows(2) {
                let [a, b] = [w[0], w[1]];
                if (a[1] > point[1]) != (b[1] > point[1]) {
                    let x = a[0] + (point[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
                    if x > point[0] {
                        inside = !inside;
                    }
                }
            }
            inside
        };
        let mut faces: Vec<Vec<usize>> = Vec::new();
        let outers: Vec<usize> = (0..traced.len())
            .filter(|&i| traced[i].winding == 0 && traced[i].area * sigma > 0.)
            .collect();
        for &o in &outers {
            faces.push(vec![o]);
        }
        // Winding loops pair from below: material above opens a band.
        let mut winding: Vec<usize> = (0..traced.len())
            .filter(|&i| traced[i].winding != 0)
            .collect();
        let mean = |i: usize| {
            traced[i].points.iter().map(|p| p[1]).sum::<f64>() / traced[i].points.len() as f64
        };
        winding.sort_by(|&a, &b| mean(a).total_cmp(&mean(b)));
        let mut bands: Vec<(usize, Vec<usize>, [f64; 2])> = Vec::new();
        let mut open: Option<usize> = None;
        for i in winding {
            let above = traced[i].winding as f64 * sigma > 0.;
            match (open, above) {
                (None, true) => open = Some(i),
                (Some(low), false) => {
                    bands.push((faces.len(), vec![low, i], [mean(low), mean(i)]));
                    faces.push(vec![low, i]);
                    open = None;
                }
                _ => return Err(refused("Transverse bore band loops do not pair")),
            }
        }
        if open.is_some() {
            return Err(refused("Transverse bore band loops do not pair"));
        }
        for i in 0..traced.len() {
            if traced[i].winding != 0 || traced[i].area * sigma > 0. {
                continue;
            }
            let probe = traced[i].points[0];
            let owner = outers
                .iter()
                .enumerate()
                .filter(|&(_, &o)| contains(&traced[o], probe))
                .min_by(|a, b| traced[*a.1].area.abs().total_cmp(&traced[*b.1].area.abs()))
                .map(|(face, _)| face)
                .or_else(|| {
                    bands
                        .iter()
                        .find(|(_, _, [low, high])| probe[1] > *low && probe[1] < *high)
                        .map(|(face, _, _)| *face)
                })
                .ok_or_else(|| refused("A transverse bore hole has no surrounding face"))?;
            faces[owner].push(i);
        }
        Ok(faces
            .into_iter()
            .map(|members| {
                (
                    brep::Face {
                        surface: surface.clone(),
                        same_sense: sigma > 0.,
                        loops: members
                            .into_iter()
                            .map(|i| Loop {
                                uses: traced[i].uses.clone(),
                                vertex: None,
                            })
                            .collect(),
                    },
                    surf,
                )
            })
            .collect())
    }

    /// Removed volume by quadrature of the chord areas over each band.
    fn removed_volume(&self) -> f64 {
        let mut total = 0.;
        for band in &self.bands {
            let r = band.r;
            let mut cuts = vec![band.z[0], band.z[1]];
            for seg in self.rings.iter().flatten() {
                let [a, b] = seg.ends();
                cuts.push(a[1]);
                cuts.push(b[1]);
                if let Seg::Arc { c, r: radius, .. } = *seg {
                    cuts.push(c[1] + radius);
                    cuts.push(c[1] - radius);
                    for y in [r, -r] {
                        let dy = y - c[0];
                        if dy.abs() < radius {
                            let root = (radius * radius - dy * dy).sqrt();
                            cuts.push(c[1] + root);
                            cuts.push(c[1] - root);
                        }
                    }
                } else if a[0] != b[0] {
                    for y in [r, -r] {
                        let f = (y - a[0]) / (b[0] - a[0]);
                        if (0. ..=1.).contains(&f) {
                            cuts.push(a[1] + f * (b[1] - a[1]));
                        }
                    }
                }
            }
            cuts.retain(|z| *z >= band.z[0] && *z <= band.z[1]);
            cuts.sort_by(f64::total_cmp);
            cuts.dedup();
            let area = |z: f64| self.chord_area(z, r);
            for w in cuts.windows(2) {
                total += tanh_sinh(&area, w[0], w[1]);
            }
        }
        total
    }

    /// Area of `x^2 + y^2 < r^2` over the profile slice at `z`.
    fn chord_area(&self, z: f64, r: f64) -> f64 {
        let mut ys = Vec::new();
        for seg in self.rings.iter().flatten() {
            match *seg {
                Seg::Line { a, b } => {
                    if (a[1] > z) != (b[1] > z) {
                        ys.push(a[0] + (z - a[1]) * (b[0] - a[0]) / (b[1] - a[1]));
                    }
                }
                Seg::Arc { c, r: radius, .. } => {
                    let dz = z - c[1];
                    if dz.abs() < radius {
                        let root = (radius * radius - dz * dz).sqrt();
                        for dy in [root, -root] {
                            if seg.fraction(dz.atan2(dy)).is_some_and(|t| t > 0. && t < 1.) {
                                ys.push(c[0] + dy);
                            }
                        }
                    }
                }
            }
        }
        ys.sort_by(f64::total_cmp);
        let g = |y: f64| {
            let y = y.clamp(-r, r);
            y * ((r - y) * (r + y)).max(0.).sqrt() + r * r * (y / r).asin()
        };
        ys.chunks(2)
            .filter(|pair| pair.len() == 2)
            .map(|pair| g(pair[1]) - g(pair[0]))
            .sum()
    }

    fn profile_area(&self) -> f64 {
        self.rings.iter().flatten().map(Seg::signed_area).sum()
    }
}

/// Double-exponential quadrature, robust to endpoint root singularities.
fn tanh_sinh(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    if b <= a {
        return 0.;
    }
    let half = 0.5 * (b - a);
    let center = 0.5 * (a + b);
    let mut previous = f64::NAN;
    let mut h = 0.5;
    for _ in 0..8 {
        let mut sum = 0.;
        let mut k = 0i64;
        loop {
            let t = k as f64 * h;
            let u = FRAC_PI_2 * t.sinh();
            let x = u.tanh();
            let w = FRAC_PI_2 * t.cosh() / u.cosh().powi(2);
            if w < 1e-300 || x >= 1. {
                break;
            }
            let terms = if k == 0 {
                f(center)
            } else {
                f(center + half * x) + f(center - half * x)
            };
            sum += w * terms;
            k += 1;
            if k > 100_000 {
                break;
            }
        }
        let value = h * half * sum;
        if (value - previous).abs() <= 1e-14 * value.abs().max(1e-300) {
            return value;
        }
        previous = value;
        h *= 0.5;
    }
    previous
}

/// The union of the steps as disjoint bands of positive radius, and the
/// band ends.
fn bands(steps: &[BoreStep]) -> (Vec<Band>, Vec<Level>) {
    let mut stops: Vec<f64> = steps.iter().flat_map(|s| [s.start_mm, s.end_mm]).collect();
    stops.sort_by(f64::total_cmp);
    stops.dedup();
    let mut bands: Vec<Band> = Vec::new();
    for w in stops.windows(2) {
        let middle = 0.5 * (w[0] + w[1]);
        let r = steps
            .iter()
            .filter(|s| s.start_mm < middle && middle < s.end_mm)
            .map(|s| s.radius_mm)
            .fold(0., f64::max);
        if r <= 0. {
            continue;
        }
        match bands.last_mut() {
            Some(last) if last.z[1] == w[0] && last.r == r => last.z[1] = w[1],
            _ => bands.push(Band { z: [w[0], w[1]], r }),
        }
    }
    let mut levels: Vec<Level> = Vec::new();
    for (k, band) in bands.iter().enumerate() {
        match levels.last_mut() {
            Some(last) if last.z == band.z[0] => last.above = Some(k),
            _ => levels.push(Level {
                z: band.z[0],
                below: None,
                above: Some(k),
            }),
        }
        levels.push(Level {
            z: band.z[1],
            below: Some(k),
            above: None,
        });
    }
    (bands, levels)
}

impl Solid {
    /// The exact retained boundary of a prism less a coaxial stepped tool
    /// crossing its extrusion at a right angle; see the module notes.
    pub(crate) fn rust_transverse_bore_from(
        source: &PrismSource,
        tool: &TransverseBoreTool,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        let tolerance = policy.linear_tolerance_mm;
        policy.point(tool.origin_mm)?;
        if tool.steps.is_empty()
            || tool.steps.iter().any(|s| {
                !(s.start_mm.is_finite() && s.end_mm.is_finite() && s.end_mm > s.start_mm)
                    || policy.length(s.radius_mm).is_err()
            })
        {
            return Err(invalid(
                "A transverse bore needs finite steps with positive radii",
            ));
        }
        source.frame.checked(policy)?;
        let frame = source.frame;
        let ex = unit(frame.normal()).ok_or_else(|| invalid("Degenerate prism frame"))?;
        let ez = unit(tool.axis).ok_or_else(|| invalid("A transverse bore needs an axis"))?;
        let magnitude = tool
            .origin_mm
            .iter()
            .chain(&frame.origin)
            .map(|v| v.abs())
            .chain(source.span.iter().map(|v| v.abs()))
            .chain(
                tool.steps
                    .iter()
                    .flat_map(|s| [s.start_mm.abs(), s.end_mm.abs()]),
            )
            .fold(1., f64::max);
        if dot(ex, ez).abs() > 32. * f64::EPSILON {
            return Err(refused(
                "The transverse bore axis is not perpendicular to the extrusion",
            ));
        }
        let ez = unit(mix(ez, -dot(ex, ez), ex)).expect("perpendicular axis");
        let ey = cross(ez, ex);
        let o = tool.origin_mm;
        let to_chart = |p: P| {
            let d = sub(p, o);
            [dot(d, ex), dot(d, ey), dot(d, ez)]
        };
        let profile = |q: Q2| {
            let world = mix(mix(frame.origin, q[0], frame.u), q[1], frame.v);
            let c = to_chart(world);
            [c[1], c[2]]
        };
        let base = dot(sub(frame.origin, o), ex);
        let x = [base + source.span[0], base + source.span[1]];
        let mut rings = Vec::new();
        let mut scale_mm: f64 = magnitude;
        for ring in &source.rings {
            let mut segs = Vec::new();
            for edge in ring {
                use multicad_rust_kernel::normal_sweep::Edge as E;
                let seg = match *edge {
                    E::Line { a, b } => Seg::Line {
                        a: profile(a),
                        b: profile(b),
                    },
                    E::Arc { center, a, b, .. } => {
                        let c = profile(center);
                        let [pa, pb] = [profile(a), profile(b)];
                        let r = edge.radius_mm().ok_or_else(|| invalid("Degenerate arc"))?;
                        Seg::Arc {
                            c,
                            r,
                            a: pa,
                            b: pb,
                            start: (pa[1] - c[1]).atan2(pa[0] - c[0]),
                            sweep: edge
                                .signed_sweep_rad()
                                .ok_or_else(|| invalid("Degenerate arc"))?,
                        }
                    }
                    E::Circle {
                        center,
                        radius_mm,
                        clockwise,
                    } => {
                        let c = profile(center);
                        let start = profile([center[0] + radius_mm, center[1]]);
                        Seg::Arc {
                            c,
                            r: radius_mm,
                            a: start,
                            b: start,
                            start: (start[1] - c[1]).atan2(start[0] - c[0]),
                            sweep: if clockwise { -TAU } else { TAU },
                        }
                    }
                };
                for p in seg.ends() {
                    scale_mm = scale_mm.max(p[0].abs()).max(p[1].abs());
                }
                segs.push(seg);
            }
            rings.push(segs);
        }
        let mut grouped = Vec::new();
        let mut i = 0;
        for ring in &rings {
            grouped.push((i..i + ring.len()).collect::<Vec<_>>());
            i += ring.len();
        }
        // The chart keeps the profile orientation: y, z is right-handed
        // about x, and the outer ring stays counterclockwise.
        for (index, ring) in rings.iter_mut().enumerate() {
            let area: f64 = ring.iter().map(Seg::signed_area).sum();
            if (area > 0.) != (index == 0) {
                ring.reverse();
                for seg in ring.iter_mut() {
                    *seg = seg.reversed();
                }
                grouped[index].reverse();
            }
        }
        let (bands, levels) = bands(&tool.steps);
        let reach = bands.iter().map(|b| b.r).fold(0., f64::max);
        let margin = RESOLUTION * scale_mm;
        if policy.length(-x[0] - reach - margin).is_err()
            || policy.length(x[1] - reach - margin).is_err()
        {
            return Err(refused(format!(
                "A transverse bore must stay clear of the prism ends (extent {x:?}, tool radius {reach})"
            )));
        }
        let mut cut = Cut {
            rings,
            side: grouped,
            x,
            bands,
            levels,
            margin,
            snap: 64. * f64::EPSILON * scale_mm,
            tolerance,
            points: BTreeMap::new(),
            splits: vec![],
            vertices: vec![],
            lifted: BTreeMap::new(),
            edges: vec![],
        };
        cut.split()?;
        cut.cap_edges()?;
        cut.long_edges()?;
        cut.window_edges()?;
        cut.step_edges()?;
        cut.shoulder_edges()?;
        let mut surfaces: BTreeSet<Surf> = BTreeSet::new();
        for record in &cut.edges {
            surfaces.extend(record.surfaces);
        }
        let step = (1e-6_f64 * scale_mm).max(64. * tolerance);
        let mut faces = Vec::new();
        for surf in surfaces {
            if interrupt.is_some_and(|stop| stop()) {
                return Err(error(OPERATION, GeometryErrorKind::Cancelled, "Cancelled"));
            }
            faces.extend(cut.faces_of(surf, step)?);
        }
        let removed = cut.removed_volume();
        let expected = cut.profile_area() * (x[1] - x[0]) - removed;
        let names: Vec<FaceName> = faces
            .iter()
            .map(|(_, surf)| FaceName {
                origins: match *surf {
                    Surf::Cap(end) => source.caps[usize::from(end)].clone(),
                    Surf::Side(side) => source.sides[side].clone(),
                    Surf::Wall(_) => vec![tool.wall.clone()],
                    Surf::Shoulder(_) => vec![tool.shoulder.clone()],
                },
                support: None,
                anchor: None,
            })
            .collect();
        let boundary = Boundary {
            vertices: cut.vertices.clone(),
            edges: cut.edges.iter().map(|r| r.edge.clone()).collect(),
            shells: vec![Shell {
                faces: (0..faces.len()).collect(),
                reversed: false,
            }],
            faces: faces.into_iter().map(|(face, _)| face).collect(),
        };
        let solid = brep::Solid::new(boundary, tolerance, interrupt).map_err(failure)?;
        if (solid.volume_mm3() - expected).abs()
            > solid.volume_error_mm3().max(1e-9 * expected.abs())
        {
            return Err(error(
                OPERATION,
                GeometryErrorKind::InvalidResult,
                format!(
                    "Rust transverse bore boundary does not reproduce its quadrature volume: {} != {}",
                    solid.volume_mm3(),
                    expected
                ),
            ));
        }
        let world = |p: P| mix(mix(mix(o, p[0], ex), p[1], ey), p[2], ez);
        let linear = |v: P| mix(mix(scale(ex, v[0]), v[1], ey), v[2], ez);
        let solid = solid.moved(&world, &linear).map_err(failure)?;
        Self::from_retained_brep(
            solid,
            source.feature,
            Rc::new(AuthoredNaming { faces: names }),
            policy,
        )
    }
}

impl Solid {
    /// Cut an unpierced normal sweep or a plain cylinder with a coaxial
    /// stepped tool whose axis crosses the extrusion at a right angle.
    pub fn rust_transverse_bore(
        &self,
        tool: &TransverseBoreTool,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if policy != self.policy || policy.boolean_fuzzy_mm != 0. {
            return Err(invalid(
                "A transverse bore requires the source policy without Boolean fuzz",
            ));
        }
        let source = self.prism_source()?.ok_or_else(|| {
            refused("A transverse bore needs an unpierced normal sweep or a plain cylinder")
        })?;
        let mut result = Self::rust_transverse_bore_from(&source, tool, policy, None)?;
        let (bands, _) = bands(&tool.steps);
        let record = TransverseRecord {
            source: self.share()?,
            origin_mm: tool.origin_mm,
            axis: unit(tool.axis).expect("validated tool axis"),
            bands: bands.iter().map(|b| [b.z[0], b.z[1], b.r]).collect(),
        };
        if let SolidBody::ImportedBrep(body) = &mut result.body {
            Rc::get_mut(body).expect("a fresh retained body").transverse = Some(Rc::new(record));
        }
        Ok(result)
    }

    /// The original prism and tool of a transverse bore.
    pub(crate) fn transverse_record(&self) -> Option<&TransverseRecord> {
        match &self.body {
            SolidBody::ImportedBrep(body) => body.transverse.as_deref(),
            _ => None,
        }
    }

    /// The profile prism of an unpierced normal sweep or plain cylinder.
    pub(crate) fn prism_source(&self) -> Result<Option<PrismSource>, GeometryError> {
        let owner = EntityId::new();
        let faces =
            |snapshot: &TopologySnapshot| -> Vec<(u32, Vec<Origin>, Option<AnalyticGeometry>)> {
                let mut faces: Vec<_> = snapshot
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == TopologyKind::Face)
                    .map(|e| {
                        (
                            e.key.ordinal,
                            e.origins.iter().cloned().collect(),
                            e.geometry,
                        )
                    })
                    .collect();
                faces.sort_by_key(|f| f.0);
                faces
            };
        if let Some((body, frame)) = self.normal_sweep_body() {
            let Some((rings, span)) = body.prism_profile() else {
                return Ok(None);
            };
            let snapshot = self.topology_snapshot(owner, 0)?;
            let mut faces = faces(&snapshot).into_iter().map(|f| f.1);
            let count: usize = rings.iter().map(Vec::len).sum();
            let caps = [faces.next(), faces.next()];
            let sides: Vec<Vec<Origin>> = faces.take(count).collect();
            let [Some(start), Some(end)] = caps else {
                return Ok(None);
            };
            if sides.len() != count {
                return Ok(None);
            }
            return Ok(Some(PrismSource {
                frame,
                span,
                rings,
                caps: [start, end],
                sides,
                feature: body.feature,
            }));
        }
        // A placed Rust cylinder starts at its base rim centre and runs
        // along its frame normal.
        let (Some((cylinder, frame)), Some((_, _, radius, height))) =
            (self.axial_cylinder_frame(), self.rust_cylinder_axis())
        else {
            return Ok(None);
        };
        let span = [0., height];
        let snapshot = self.topology_snapshot(owner, 0)?;
        let normal = frame.normal();
        let mut caps: [Vec<Origin>; 2] = [vec![], vec![]];
        let mut sides = vec![];
        for (_, origins, geometry) in faces(&snapshot) {
            match geometry {
                Some(AnalyticGeometry::Plane { frame: plane }) => {
                    let level = dot(sub(plane.origin_mm, frame.origin), normal);
                    let end = (level - span[1]).abs() < (level - span[0]).abs();
                    caps[usize::from(end)] = origins;
                }
                _ => sides = origins,
            }
        }
        Ok(Some(PrismSource {
            frame,
            span,
            rings: vec![vec![multicad_rust_kernel::normal_sweep::Edge::Circle {
                center: [0., 0.],
                radius_mm: radius,
                clockwise: false,
            }]],
            caps,
            sides: vec![sides],
            feature: EntityId::from_uuid_bits(cylinder.feature()),
        }))
    }
}

#[cfg(test)]
mod tests;
