//! Parallel planar faces whose boundaries include circular arcs.
//!
//! Every pair of points on two parallel planes is at least their separation
//! apart, so a point strictly inside one face whose normal projection lies
//! strictly inside the other proves the separation is the exact minimum.
//! Membership uses exact distances to the line and arc edges: a point
//! farther from every edge than the arc sagitta of a fixed chord subdivision
//! has the same inclusion in the face as in its chord polygon, which an
//! even-odd crossing count decides. Faces whose projections do not overlap
//! that way are refused, not approximated.
use super::*;
use multicad_core::references::{AnalyticGeometry, TopologyKind};
use multicad_precision::PortableMath as _;
use std::collections::BTreeMap;
use std::f64::consts::TAU;

type Point = [f64; 3];
type Flat = [f64; 2];
/// Chords per full turn of an arc for the crossing count.
const CHORDS: f64 = 128.;
const GRID: usize = 16;

#[derive(Clone, Copy)]
enum Edge {
    Line([Flat; 2]),
    Arc {
        center: Flat,
        radius: f64,
        start: f64,
        sweep: f64,
    },
}

struct Plane {
    origin: Point,
    normal: Point,
    axes: [Point; 2],
}

fn dot(a: Point, b: Point) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn norm(a: Point) -> f64 {
    a[0].portable_hypot(a[1]).portable_hypot(a[2])
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn unit(a: Point) -> Option<Point> {
    let n = norm(a);
    (n.is_finite() && n > 0.).then(|| a.map(|x| x / n))
}

impl Plane {
    fn new(origin: Point, normal: Point) -> Option<Self> {
        let normal = unit(normal)?;
        let seed = if normal[0].abs() < 0.9 {
            [1., 0., 0.]
        } else {
            [0., 1., 0.]
        };
        let u = unit(cross(normal, seed))?;
        let v = cross(normal, u);
        Some(Self {
            origin,
            normal,
            axes: [u, v],
        })
    }

    fn flat(&self, p: Point) -> Flat {
        let d = sub(p, self.origin);
        self.axes.map(|axis| dot(d, axis))
    }

    fn point(&self, [x, y]: Flat) -> Point {
        std::array::from_fn(|i| self.origin[i] + x * self.axes[0][i] + y * self.axes[1][i])
    }
}

fn segment_distance(p: Flat, [a, b]: [Flat; 2]) -> f64 {
    let v = [b[0] - a[0], b[1] - a[1]];
    let w = [p[0] - a[0], p[1] - a[1]];
    let t = ((w[0] * v[0] + w[1] * v[1]) / (v[0] * v[0] + v[1] * v[1])).clamp(0., 1.);
    (w[0] - t * v[0]).portable_hypot(w[1] - t * v[1])
}

impl Edge {
    fn distance(&self, p: Flat) -> f64 {
        match *self {
            Self::Line(segment) => segment_distance(p, segment),
            Self::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                let at = |angle: f64| {
                    let (s, c) = angle.portable_sin_cos();
                    [center[0] + radius * c, center[1] + radius * s]
                };
                let d = [p[0] - center[0], p[1] - center[1]];
                let offset = (d[1].portable_atan2(d[0]) - start).rem_euclid(TAU);
                let within = if sweep >= 0. {
                    offset <= sweep
                } else {
                    TAU - offset <= -sweep || offset == 0.
                };
                if within || sweep.abs() >= TAU {
                    (d[0].portable_hypot(d[1]) - radius).abs()
                } else {
                    let [a, b] = [at(start), at(start + sweep)];
                    (p[0] - a[0])
                        .portable_hypot(p[1] - a[1])
                        .min((p[0] - b[0]).portable_hypot(p[1] - b[1]))
                }
            }
        }
    }

    /// Chord subdivision and the largest distance of a chord from its arc.
    fn chords(&self) -> (Vec<[Flat; 2]>, f64) {
        match *self {
            Self::Line(segment) => (vec![segment], 0.),
            Self::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                let count = (sweep.abs() / TAU * CHORDS).ceil().max(1.) as usize;
                let step = sweep / count as f64;
                let at = |i: usize| {
                    let (s, c) = (start + step * i as f64).portable_sin_cos();
                    [center[0] + radius * c, center[1] + radius * s]
                };
                let sagitta = radius * (1. - (step.abs() / 2.).portable_cos());
                (
                    (0..count).map(|i| [at(i), at(i + 1)]).collect(),
                    sagitta * (1. + 1e-9) + 1e-12 * radius,
                )
            }
        }
    }
}

struct Region {
    edges: Vec<Edge>,
    chords: Vec<[Flat; 2]>,
    sagitta: f64,
}

impl Region {
    fn new(edges: Vec<Edge>) -> Self {
        let mut chords = vec![];
        let mut sagitta = 0_f64;
        for edge in &edges {
            let (c, s) = edge.chords();
            chords.extend(c);
            sagitta = sagitta.max(s);
        }
        Self {
            edges,
            chords,
            sagitta,
        }
    }

    /// Strict interior with a clearance from every edge.
    fn inside(&self, p: Flat, clearance: f64) -> bool {
        if self
            .edges
            .iter()
            .any(|edge| edge.distance(p) <= self.sagitta + clearance)
        {
            return false;
        }
        let mut inside = false;
        for &[a, b] in &self.chords {
            if (a[1] > p[1]) != (b[1] > p[1]) {
                let x = a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
                inside ^= x > p[0];
            }
        }
        inside
    }
}

/// Face plane, its boundary in the plane chart of `chart`, and whether any
/// boundary edge is curved. `None` for faces this provider does not cover.
fn region(solid: &Solid, face: u32, chart: Option<&Plane>) -> Option<(Plane, Region, bool)> {
    let (origin, normal) = solid
        .faces()
        .ok()?
        .into_iter()
        .find(|f| f.index == face)?
        .plane?;
    let own = Plane::new(origin, normal)?;
    let chart = chart.unwrap_or(&own);
    let topology = solid.rust_topology()?.ok()?;
    let geometry: BTreeMap<_, _> = topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Edge)
        .map(|e| (e.key.ordinal, e.geometry))
        .collect();
    let mut edges = vec![];
    let mut curved = false;
    for (index, _) in crate::rust_display::planar_loops(solid, face)
        .into_iter()
        .flatten()
    {
        edges.push(match (*geometry.get(&index)?)? {
            AnalyticGeometry::Line { a_mm, b_mm } => {
                Edge::Line([chart.flat(a_mm), chart.flat(b_mm)])
            }
            AnalyticGeometry::CircularArc {
                frame,
                radius_mm,
                start_rad,
                sweep_rad,
            } => {
                curved = true;
                let n = unit(frame.normal())?;
                let u = unit(frame.u)?;
                // The arc chart must lie in the face plane, with its own
                // rotation sense carried by the sign of the sweep.
                if norm(cross(n, chart.normal)) > 1e-12 {
                    return None;
                }
                let sense = dot(n, chart.normal).signum();
                let [ux, uy] = chart.axes.map(|axis| dot(u, axis));
                let phase = uy.portable_atan2(ux);
                Edge::Arc {
                    center: chart.flat(frame.origin_mm),
                    radius: radius_mm,
                    start: phase + sense * start_rad,
                    sweep: sense * sweep_rad,
                }
            }
            _ => return None,
        });
    }
    (!edges.is_empty()).then(|| (own, Region::new(edges), curved))
}

pub(super) fn measure(
    a: &Solid,
    sa: BRepSelection,
    b: &Solid,
    sb: BRepSelection,
    tolerance: f64,
) -> Option<Result<MinimumDistance, GeometryError>> {
    let (BRepSelection::Face(fa), BRepSelection::Face(fb)) = (sa, sb) else {
        return None;
    };
    let (plane, left, curved_a) = region(a, fa, None)?;
    let (other, right, curved_b) = region(b, fb, Some(&plane))?;
    let tilt = norm(cross(plane.normal, other.normal));
    if !(curved_a || curved_b) || tilt > 1e-12 {
        return None;
    }
    let separation = dot(sub(other.origin, plane.origin), plane.normal);
    // Chart and projection arithmetic over the face extents.
    let extent = left
        .chords
        .iter()
        .chain(&right.chords)
        .flatten()
        .flat_map(|p| p.map(f64::abs))
        .chain(plane.origin.map(f64::abs))
        .chain(other.origin.map(f64::abs))
        .fold(separation.abs().max(1.), f64::max);
    let clearance = tolerance * 4. + (tilt + 64. * f64::EPSILON) * extent;
    Some((|| {
        let [lo, hi] = left.chords.iter().flatten().fold(
            [[f64::INFINITY; 2], [f64::NEG_INFINITY; 2]],
            |[lo, hi], p| {
                [
                    [lo[0].min(p[0]), lo[1].min(p[1])],
                    [hi[0].max(p[0]), hi[1].max(p[1])],
                ]
            },
        );
        let candidates = (0..=GRID).flat_map(|i| {
            (0..=GRID).map(move |j| {
                // Interleave so central points come first.
                let f = |k: usize| 0.5 + (k as f64 / GRID as f64 - 0.5) * 0.999;
                [
                    lo[0] + (hi[0] - lo[0]) * f(i),
                    lo[1] + (hi[1] - lo[1]) * f(j),
                ]
            })
        });
        let found = candidates
            .into_iter()
            .find(|&p| left.inside(p, clearance) && right.inside(p, clearance))
            .ok_or_else(|| {
                error(
                    GeometryErrorKind::Unsupported,
                    "Parallel planar faces with circular boundaries need overlapping normal projections for an exact distance",
                )
            })?;
        let pa = plane.point(found);
        let pb: Point = std::array::from_fn(|i| pa[i] + separation * plane.normal[i]);
        let distance = norm(sub(pb, pa));
        if !distance.is_finite() || (distance - separation.abs()).abs() > tolerance / 4. {
            return Err(error(
                GeometryErrorKind::Unsupported,
                "Parallel face distance witness exceeds the numerical budget",
            ));
        }
        Ok(MinimumDistance {
            distance_mm: separation.abs(),
            witnesses: vec![DistanceWitness { a_mm: pa, b_mm: pb }],
            witnesses_complete: false,
            total_solutions: None,
            inner_solution: false,
        })
    })())
}
