//! The common region of two placed exact gear sections.
//!
//! Each section is a set of directed exact loops (lines, arcs and involutes)
//! with material on the left. Every crossing of two placed edges lies in a
//! pair of parameter boxes whose chord tubes meet, because each tube encloses
//! its exact edge by the curvature bound; boxes whose tubes are apart are
//! discarded and the rest are halved until a Newton step on the exact edges
//! settles the crossing. Edges are split at the crossings; a piece is kept
//! when it lies inside the other section. A piece that starts or ends at a
//! transversal crossing is inside exactly when it turns into the other
//! section's material there; the disc proof decides pieces without
//! crossings and cross-checks every piece whose midpoint it resolves.
//! Kept pieces must link into closed loops through shared crossings and edge
//! vertices. Tangential contacts, crossings at vertices, unresolved or
//! conflicting pieces and broken links leave the region unresolved
//! (`Ok(None)`): nothing is
//! sampled, approximated or reconstructed.
use super::{Error, Placement, checkpoint, curvature_bound, witness};
use crate::{
    normal_sweep,
    revolution::length,
    spur_gear::{Edge, Point, Profile, internal::Ring, sub},
};
use std::collections::BTreeMap;

const MAX_BOXES: usize = 4_000_000;
const MAX_CROSSINGS: usize = 65_536;

/// One exact gear section in its own transverse chart.
#[derive(Clone, Copy)]
pub enum Section<'a> {
    /// An external tooth outline minus its optional bore.
    Spur(&'a Profile),
    /// The material between an internal ring's rim and its tooth cavity.
    Ring(&'a Ring),
}

/// A trimmed exact edge of one section, in that section's own chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Piece {
    pub section: usize,
    pub edge: Edge,
}

/// One connected common region: its outer loop and any holes, each loop
/// directed with the common material on its left.
#[derive(Clone, Debug, PartialEq)]
pub struct Region {
    pub outer: Vec<Piece>,
    pub holes: Vec<Vec<Piece>>,
}

/// The Green integral of one placed directed edge, one half of the
/// integral of p x dp, in closed form.
fn signed_area(piece: Piece, placement: Placement) -> f64 {
    let cross = |a: Point, b: Point| a[0] * b[1] - a[1] * b[0];
    let edge = piece.edge;
    let local = match edge {
        Edge::Planar(normal_sweep::Edge::Line { a, b }) => cross(a, b) / 2.,
        Edge::Planar(arc @ normal_sweep::Edge::Arc { center, a, b, .. }) => {
            let radius = arc.radius_mm().unwrap_or(0.);
            let sweep = arc.signed_sweep_rad().unwrap_or(0.);
            (cross(center, sub(b, a)) + radius * radius * sweep) / 2.
        }
        Edge::Planar(normal_sweep::Edge::Circle {
            radius_mm,
            clockwise,
            ..
        }) => {
            let sign = if clockwise { -1. } else { 1. };
            sign * std::f64::consts::PI * radius_mm * radius_mm
        }
        Edge::Involute(flank) => flank.signed_area_mm2(),
    };
    // A rigid motion keeps cross products; the translation adds c x (b - a).
    let [u, v] = placement.axes;
    let chord = sub(edge.point(1.), edge.point(0.));
    let moved = [
        u[0] * chord[0] + v[0] * chord[1],
        u[1] * chord[0] + v[1] * chord[1],
    ];
    local + cross(placement.center_mm, moved) / 2.
}

impl Region {
    /// The exact area of this region in the common chart.
    pub fn area_mm2(&self, placements: [Placement; 2]) -> f64 {
        self.outer
            .iter()
            .chain(self.holes.iter().flatten())
            .map(|piece| signed_area(*piece, placements[piece.section]))
            .sum()
    }
}

fn reversed(edge: Edge) -> Edge {
    match edge {
        Edge::Planar(edge) => Edge::Planar(edge.reversed()),
        Edge::Involute(edge) => Edge::Involute(edge.reversed()),
    }
}

/// The edge between two fractions, directed as the original.
fn trimmed(edge: Edge, [from, to]: [f64; 2]) -> Edge {
    if from == 0. && to == 1. {
        return edge;
    }
    match edge {
        Edge::Planar(normal_sweep::Edge::Line { .. }) => Edge::Planar(normal_sweep::Edge::Line {
            a: edge.point(from),
            b: edge.point(to),
        }),
        Edge::Planar(normal_sweep::Edge::Arc {
            center, clockwise, ..
        }) => Edge::Planar(normal_sweep::Edge::Arc {
            center,
            a: edge.point(from),
            b: edge.point(to),
            clockwise,
        }),
        Edge::Planar(normal_sweep::Edge::Circle { .. }) => {
            unreachable!("circles are split into arcs")
        }
        Edge::Involute(flank) => Edge::Involute(flank.trimmed([from, to])),
    }
}

/// A full circle as two directed half arcs.
fn circle_arcs(edge: normal_sweep::Edge) -> Result<Vec<Edge>, Error> {
    let normal_sweep::Edge::Circle {
        center,
        radius_mm,
        clockwise,
    } = edge
    else {
        return Ok(vec![Edge::Planar(edge)]);
    };
    let a = [center[0] + radius_mm, center[1]];
    let b = [center[0] - radius_mm, center[1]];
    if radius_mm.is_nan() || radius_mm <= 0. {
        return Err(Error::Invalid);
    }
    Ok([(a, b), (b, a)]
        .into_iter()
        .map(|(a, b)| {
            Edge::Planar(normal_sweep::Edge::Arc {
                center,
                a,
                b,
                clockwise,
            })
        })
        .collect())
}

/// Shoelace area of a loop from its edge samples.
fn sampled_area(edges: &[Edge]) -> f64 {
    let mut points = Vec::new();
    for edge in edges {
        let count = 32;
        points.extend((0..count).map(|i| edge.point(i as f64 / count as f64)));
    }
    points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .map(|(a, b)| a[0] * b[1] - a[1] * b[0])
        .sum::<f64>()
        / 2.
}

/// Directed loops of a section: positive outer loops, negative holes.
pub(super) fn loops(section: Section<'_>) -> Result<Vec<Vec<Edge>>, Error> {
    let (outer, hole): (Vec<Edge>, Option<Vec<Edge>>) = match section {
        Section::Spur(profile) => (
            profile.outline().iter().map(|item| item.edge).collect(),
            profile.bore().map(circle_arcs).transpose()?,
        ),
        Section::Ring(ring) => (
            circle_arcs(ring.rim())?,
            Some(ring.cavity().iter().map(|item| item.edge).collect()),
        ),
    };
    let directed = |edges: Vec<Edge>, positive: bool| {
        if (sampled_area(&edges) > 0.) == positive {
            edges
        } else {
            edges.into_iter().rev().map(reversed).collect()
        }
    };
    Ok(std::iter::once(directed(outer, true))
        .chain(hole.map(|edges| directed(edges, false)))
        .collect())
}

/// The local point of a chart point under a rigid placement.
fn local(placement: Placement, point: Point) -> Point {
    let d = sub(point, placement.center_mm);
    let [u, v] = placement.axes;
    [d[0] * u[0] + d[1] * u[1], d[0] * v[0] + d[1] * v[1]]
}

fn strictly_inside(
    section: Section<'_>,
    point: Point,
    tolerance: f64,
) -> Result<Option<bool>, Error> {
    use super::DiscSide::*;
    let side = match section {
        Section::Spur(profile) => super::classify_spur_disc(profile, point, 0., tolerance)?,
        Section::Ring(ring) => super::classify_ring_disc(ring, point, 0., tolerance)?,
    };
    Ok(match side {
        Inside => Some(true),
        Outside => Some(false),
        Unresolved => None,
    })
}

/// One placed exact edge with its chord-error bound.
#[derive(Clone, Copy)]
struct Placed {
    edge: Edge,
    placement: Placement,
    second: f64,
}

impl Placed {
    fn at(&self, t: f64) -> Point {
        self.placement.point(self.edge.point(t))
    }

    /// Chord over [t0, t1], the sagitta bound of the exact edge from it,
    /// and the radius of a tube enclosing the placed edge.
    fn chord(&self, [t0, t1]: [f64; 2], guard: f64) -> ([Point; 2], f64, f64) {
        let span = t1 - t0;
        let sagitta = self.second * span * span / 8.;
        (
            [self.at(t0), self.at(t1)],
            sagitta,
            sagitta + self.placement.uncertainty_mm + guard,
        )
    }

    fn derivative(&self, t: f64) -> Point {
        let h = 1e-7;
        let (a, b) = ((t - h).max(0.), (t + h).min(1.));
        let [p, q] = [self.at(a), self.at(b)];
        [(q[0] - p[0]) / (b - a), (q[1] - p[1]) / (b - a)]
    }
}

enum Found {
    Crossings(Vec<[f64; 2]>),
    Degenerate,
}

/// Every transversal crossing of two placed edges, as fraction pairs.
fn crossings(
    pair: [Placed; 2],
    guard: f64,
    budget: &mut usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Found, Error> {
    let mut stack = vec![[[0., 1.], [0., 1.]]];
    let mut roots: Vec<[f64; 2]> = Vec::new();
    let mut seen = 0usize;
    while let Some(boxes) = stack.pop() {
        seen += 1;
        if seen.is_multiple_of(4096) {
            checkpoint(interrupt)?;
        }
        if *budget == 0 {
            return Err(Error::Limit);
        }
        *budget -= 1;
        let [(a, sa, ta), (b, sb, tb)] = [0, 1].map(|i| pair[i].chord(boxes[i], guard));
        if witness::segment_distance(a, b) > ta + tb {
            continue;
        }
        if sa.max(sb) > 1e-9 {
            let split = usize::from(sb > sa);
            let [lo, hi] = boxes[split];
            let mid = lo.midpoint(hi);
            for half in [[lo, mid], [mid, hi]] {
                let mut next = boxes;
                next[split] = half;
                stack.push(next);
            }
            continue;
        }
        // Newton on the exact edges from the leaf's middle.
        let mut t = boxes.map(|[lo, hi]| lo.midpoint(hi));
        for _ in 0..20 {
            let f = sub(pair[0].at(t[0]), pair[1].at(t[1]));
            let [p, q] = [pair[0].derivative(t[0]), pair[1].derivative(t[1])];
            let det = -p[0] * q[1] + p[1] * q[0];
            if det == 0. {
                break;
            }
            // Solve [p, -q] [dt, ds] = -f.
            let dt = (-f[0] * -q[1] - -q[0] * -f[1]) / det;
            let ds = (p[0] * -f[1] - p[1] * -f[0]) / det;
            t = [t[0] + dt, t[1] + ds];
            if dt.abs().max(ds.abs()) < 1e-15 {
                break;
            }
        }
        if !t
            .iter()
            .all(|v| v.is_finite() && (-1e-9..=1. + 1e-9).contains(v))
        {
            continue;
        }
        let gap = length(sub(pair[0].at(t[0]), pair[1].at(t[1])));
        if gap > 1e-9 + guard {
            continue;
        }
        if t.iter().any(|&v| v < 1e-9 || v > 1. - 1e-9) {
            return Ok(Found::Degenerate);
        }
        let [p, q] = [pair[0].derivative(t[0]), pair[1].derivative(t[1])];
        let sine = (p[0] * q[1] - p[1] * q[0]) / (length(p) * length(q));
        if sine.abs() < 1e-6 {
            return Ok(Found::Degenerate);
        }
        let point = pair[0].at(t[0]);
        if roots
            .iter()
            .all(|r| length(sub(pair[0].at(r[0]), point)) > 1e-7)
        {
            roots.push(t);
        }
    }
    Ok(Found::Crossings(roots))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Key {
    Crossing(usize),
    /// The start vertex of an edge: (section, loop, edge).
    Vertex(usize, usize, usize),
}

struct Kept {
    piece: Piece,
    end: Key,
}

/// The common region of two placed sections, `Ok(None)` when a contact,
/// vertex crossing or boundary point leaves it unresolved.
pub fn intersect_sections(
    sections: [Section<'_>; 2],
    placements: [Placement; 2],
    tolerance_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Option<Vec<Region>>, Error> {
    if !(1e-7..=1e-3).contains(&tolerance_mm)
        || placements.into_iter().any(|p| !p.valid(tolerance_mm))
    {
        return Err(Error::Invalid);
    }
    let loops = [loops(sections[0])?, loops(sections[1])?];
    let scale = placements
        .iter()
        .flat_map(|p| p.center_mm)
        .map(f64::abs)
        .fold(1., f64::max)
        * 4.;
    let guard = 1024. * f64::EPSILON * scale;
    let placed = |side: usize, edge: Edge| -> Result<Placed, Error> {
        Ok(Placed {
            edge,
            placement: placements[side],
            second: curvature_bound(edge)?,
        })
    };
    // Fractions of every crossing on each edge, by (section, loop, edge).
    let mut cuts: BTreeMap<(usize, usize, usize), Vec<(f64, usize)>> = BTreeMap::new();
    let mut count = 0usize;
    // The turn `p × q` from the first edge's direction to the second's at
    // every crossing, in the common chart.
    let mut turns: Vec<f64> = Vec::new();
    let mut budget = MAX_BOXES;
    let bound = |side: usize, edge: Edge| {
        // Enclosing box of the placed edge from its chord tube.
        let p = placed(side, edge)?;
        let (chord, _, tube) = p.chord([0., 1.], guard);
        let r = tube + edge.length_mm();
        let mid = [
            chord[0][0].midpoint(chord[1][0]),
            chord[0][1].midpoint(chord[1][1]),
        ];
        Ok::<_, Error>((p, mid, r))
    };
    let second: Vec<Vec<_>> = loops[1]
        .iter()
        .map(|edges| {
            edges
                .iter()
                .map(|&e| bound(1, e))
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<_, _>>()?;
    for (la, edges) in loops[0].iter().enumerate() {
        for (ea, &edge) in edges.iter().enumerate() {
            let (pa, ma, ra) = bound(0, edge)?;
            for (lb, others) in second.iter().enumerate() {
                for (eb, &(pb, mb, rb)) in others.iter().enumerate() {
                    if length(sub(ma, mb)) > ra + rb {
                        continue;
                    }
                    match crossings([pa, pb], guard, &mut budget, interrupt)? {
                        Found::Degenerate => return Ok(None),
                        Found::Crossings(found) => {
                            for t in found {
                                count += 1;
                                if count > MAX_CROSSINGS {
                                    return Err(Error::Limit);
                                }
                                let id = count - 1;
                                let [p, q] = [pa.derivative(t[0]), pb.derivative(t[1])];
                                turns.push(p[0] * q[1] - p[1] * q[0]);
                                cuts.entry((0, la, ea)).or_default().push((t[0], id));
                                cuts.entry((1, lb, eb)).or_default().push((t[1], id));
                            }
                        }
                    }
                }
            }
        }
    }
    // Split, classify and index the kept pieces by their start.
    let mut kept: Vec<Kept> = Vec::new();
    let mut by_start: BTreeMap<Key, usize> = BTreeMap::new();
    for side in 0..2 {
        let other = 1 - side;
        for (l, edges) in loops[side].iter().enumerate() {
            for (e, &edge) in edges.iter().enumerate() {
                checkpoint(interrupt)?;
                let mut marks = cuts.get(&(side, l, e)).cloned().unwrap_or_default();
                marks.sort_by(|a, b| a.0.total_cmp(&b.0));
                let next = (side, l, (e + 1) % edges.len());
                let mut stops = vec![(0., Key::Vertex(side, l, e))];
                stops.extend(marks.iter().map(|&(t, id)| (t, Key::Crossing(id))));
                stops.push((1., Key::Vertex(next.0, next.1, next.2)));
                // Material lies left of every directed loop, so at a
                // transversal crossing an edge enters the other section
                // when it turns left of the other edge's direction.
                let enters = |key: Key| match key {
                    Key::Crossing(id) => Some((turns[id] < 0.) == (side == 0)),
                    Key::Vertex(..) => None,
                };
                for pair in stops.windows(2) {
                    let [(t0, start), (t1, end)] = [pair[0], pair[1]];
                    let middle = placements[side].point(edge.point(t0.midpoint(t1)));
                    let point = local(placements[other], middle);
                    // A piece in a sliver thinner than the disc proof's
                    // tube still follows from its crossings; every answer
                    // that resolves must agree.
                    let answers = [
                        enters(start),
                        enters(end).map(|entering| !entering),
                        strictly_inside(sections[other], point, tolerance_mm)?,
                    ];
                    let mut resolved = answers.into_iter().flatten();
                    let Some(inside) = resolved.next() else {
                        return Ok(None);
                    };
                    if resolved.any(|answer| answer != inside) {
                        return Ok(None);
                    }
                    match inside {
                        false => {}
                        true => {
                            if by_start.insert(start, kept.len()).is_some() {
                                return Ok(None);
                            }
                            kept.push(Kept {
                                piece: Piece {
                                    section: side,
                                    edge: trimmed(edge, [t0, t1]),
                                },
                                end,
                            });
                        }
                    }
                }
            }
        }
    }
    // Link the kept pieces into closed loops.
    let mut used = vec![false; kept.len()];
    let mut closed: Vec<Vec<Piece>> = Vec::new();
    for first in 0..kept.len() {
        if used[first] {
            continue;
        }
        let mut current = first;
        let mut path = Vec::new();
        loop {
            if used[current] {
                if current != first {
                    return Ok(None);
                }
                break;
            }
            used[current] = true;
            path.push(kept[current].piece);
            let Some(&next) = by_start.get(&kept[current].end) else {
                return Ok(None);
            };
            current = next;
        }
        closed.push(path);
    }
    // Positive loops bound regions; each hole joins the outer loop that
    // encloses its first point.
    let chart = |piece: &Piece| {
        let placement = placements[piece.section];
        let Piece { edge, .. } = *piece;
        (0..32)
            .map(move |i| placement.point(edge.point(i as f64 / 32.)))
            .collect::<Vec<_>>()
    };
    let polygon = |path: &[Piece]| path.iter().flat_map(chart).collect::<Vec<_>>();
    let area = |points: &[Point]| {
        points
            .iter()
            .zip(points.iter().cycle().skip(1))
            .map(|(a, b)| a[0] * b[1] - a[1] * b[0])
            .sum::<f64>()
            / 2.
    };
    let mut regions = Vec::new();
    let mut holes = Vec::new();
    for path in closed {
        let points = polygon(&path);
        if area(&points) > 0. {
            regions.push((
                points,
                Region {
                    outer: path,
                    holes: vec![],
                },
            ));
        } else {
            holes.push((points[0], path));
        }
    }
    for (point, path) in holes {
        let owner = regions.iter_mut().find(|(outline, _)| {
            let mut inside = false;
            for (a, b) in outline.iter().zip(outline.iter().cycle().skip(1)) {
                if (a[1] > point[1]) != (b[1] > point[1])
                    && point[0] < a[0] + (b[0] - a[0]) * (point[1] - a[1]) / (b[1] - a[1])
                {
                    inside = !inside;
                }
            }
            inside
        });
        match owner {
            Some((_, region)) => region.holes.push(path),
            None => return Ok(None),
        }
    }
    Ok(Some(
        regions.into_iter().map(|(_, region)| region).collect(),
    ))
}

#[cfg(test)]
mod tests;
