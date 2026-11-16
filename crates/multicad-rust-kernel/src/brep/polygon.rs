//! Triangulation of a face's parameter region from its boundary points.
//!
//! The region's loops are merged into one ring: a winding pair is joined
//! across the band it bounds, a lone winding loop is closed along its pole,
//! and holes are bridged in. The ring is then clipped ear by ear, shortest
//! diagonal first. No point is added, so neighbouring faces see exactly the
//! shared edge samples.
use super::domain::{Domain, Key, LoopPoint, UvLoop};
use super::{Error, Result, Surface, ensure, norm};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Corner {
    pub uv: [f64; 2],
    pub key: Key,
    pub pole: bool,
    /// Distinguishes the copies of a pole along its parameter row.
    pub copy: usize,
}

pub(crate) type Triangle = [Corner; 3];

/// Sides of the synthetic ring around a complete closed surface.
const SIDE: usize = 12;

fn corner(point: &LoopPoint, uv: [f64; 2]) -> Corner {
    Corner {
        uv,
        key: point.key,
        pole: point.pole,
        copy: 0,
    }
}

fn orient(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// True when the segments meet anywhere but at a shared end point.
fn crosses(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    if a == c || a == d || b == c || b == d {
        return false;
    }
    let (o1, o2, o3, o4) = (
        orient(a, b, c),
        orient(a, b, d),
        orient(c, d, a),
        orient(c, d, b),
    );
    let overlap = |axis: usize| {
        a[axis].min(b[axis]) <= c[axis].max(d[axis]) && c[axis].min(d[axis]) <= a[axis].max(b[axis])
    };
    o1 * o2 <= 0. && o3 * o4 <= 0. && overlap(0) && overlap(1)
}

/// True when `p` lies on the segment from `a` to `b`, strictly between its
/// ends, within `tolerance` of the line in parameter units.
fn between(a: [f64; 2], b: [f64; 2], p: [f64; 2], tolerance: f64) -> bool {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let along = (p[0] - a[0]) * dx + (p[1] - a[1]) * dy;
    p != a
        && p != b
        && orient(a, b, p).abs() <= tolerance * dx.hypot(dy)
        && along > 0.
        && along < dx * dx + dy * dy
}

/// The widest parameter distance of the ring from its first point.
fn spread(ring: &[Corner]) -> f64 {
    ring.iter().fold(0., |e: f64, c| {
        e.max((c.uv[0] - ring[0].uv[0]).hypot(c.uv[1] - ring[0].uv[1]))
    })
}

fn inside_ring(ring: &[Corner], p: [f64; 2]) -> bool {
    let mut inside = false;
    for (i, a) in ring.iter().enumerate() {
        let b = ring[(i + 1) % ring.len()];
        if (a.uv[1] > p[1]) != (b.uv[1] > p[1]) {
            let x = a.uv[0] + (p[1] - a.uv[1]) / (b.uv[1] - a.uv[1]) * (b.uv[0] - a.uv[0]);
            inside ^= x > p[0];
        }
    }
    inside
}

struct Metric([f64; 2]);

impl Metric {
    fn distance(&self, a: [f64; 2], b: [f64; 2]) -> f64 {
        ((a[0] - b[0]) * self.0[0]).hypot((a[1] - b[1]) * self.0[1])
    }
}

/// A structured grid over a complete sphere or torus. Opposite sides
/// carry the same keys, and a sphere's lower and upper rows are its poles.
fn complete(face: usize, [low, high]: [[f64; 2]; 2], torus: bool) -> Vec<Triangle> {
    let at = |i: usize, j: usize| {
        let uv = [
            low[0] + (high[0] - low[0]) * i as f64 / SIDE as f64,
            low[1] + (high[1] - low[1]) * j as f64 / SIDE as f64,
        ];
        let pole = !torus && j.is_multiple_of(SIDE);
        let key = if pole {
            1 + j / SIDE
        } else {
            3 + (i % SIDE) * SIDE + j % SIDE
        };
        Corner {
            uv,
            key: Key::Extra(face, key),
            pole,
            copy: i % SIDE,
        }
    };
    let mut triangles = Vec::new();
    for i in 0..SIDE {
        for j in 0..SIDE {
            triangles.push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
            triangles.push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    triangles
}

/// One traversal of a loop from point `start`, closing point included.
fn unrolled(
    trim: &UvLoop,
    start: usize,
    periods: [Option<f64>; 2],
    shift: [f64; 2],
) -> Vec<Corner> {
    let n = trim.points.len();
    (start..=start + n)
        .map(|i| {
            let uv = trim.at(i, periods);
            corner(&trim.points[i % n], [uv[0] + shift[0], uv[1] + shift[1]])
        })
        .collect()
}

/// Join a winding pair across the band between them at their closest points.
fn band(a: &UvLoop, b: &UvLoop, periods: [Option<f64>; 2], metric: &Metric) -> Vec<Corner> {
    let axis = usize::from(a.winding[0] == 0);
    let period = periods[axis].unwrap_or(0.);
    let mut best = (f64::INFINITY, 0, 0, 0.);
    for i in 0..a.points.len() {
        let mut end = a.points[i].uv;
        end[axis] += period;
        for (j, point) in b.points.iter().enumerate() {
            let turns = ((end[axis] - point.uv[axis]) / period).round();
            let mut near = point.uv;
            near[axis] += turns * period;
            let d = metric.distance(end, near);
            if d < best.0 {
                best = (d, i, j, turns * period);
            }
        }
    }
    let mut shift = [0.; 2];
    shift[axis] = best.3;
    let mut ring = unrolled(a, best.1, periods, [0.; 2]);
    ring.extend(unrolled(b, best.2, periods, shift));
    ring
}

/// Move a hole by whole periods into the ring.
fn place(ring: &[Corner], hole: &UvLoop, periods: [Option<f64>; 2]) -> Result<Vec<Corner>> {
    let turns = |axis: usize| {
        if periods[axis].is_some() {
            -2..=2
        } else {
            0..=0
        }
    };
    // A hole may touch the ring at one point; its centre lies inside.
    let count = hole.points.len() as f64;
    let centre = [0, 1].map(|axis| hole.points.iter().map(|p| p.uv[axis]).sum::<f64>() / count);
    let mut placed = None;
    for ku in turns(0) {
        for kv in turns(1) {
            let shift = [
                f64::from(ku) * periods[0].unwrap_or(0.),
                f64::from(kv) * periods[1].unwrap_or(0.),
            ];
            if placed.is_none() && inside_ring(ring, [centre[0] + shift[0], centre[1] + shift[1]]) {
                placed = Some(shift);
            }
        }
    }
    let shift = placed.ok_or(Error::Invalid)?;
    Ok(hole
        .points
        .iter()
        .map(|p| corner(p, [p.uv[0] + shift[0], p.uv[1] + shift[1]]))
        .collect())
}

fn rightmost(points: &[Corner]) -> usize {
    let mut h = 0;
    for (i, point) in points.iter().enumerate() {
        if point.uv[0] > points[h].uv[0] {
            h = i;
        }
    }
    h
}

/// Join a placed hole to the ring by the shortest unobstructed bridge,
/// preferring the hole's rightmost point. `holes[index..]` are still apart.
fn bridge(
    ring: &mut Vec<Corner>,
    holes: &[Vec<Corner>],
    index: usize,
    metric: &Metric,
) -> Result<()> {
    let points = &holes[index];
    let first = rightmost(points);
    let line = 1e-9 * spread(ring);
    for offset in 0..points.len() {
        let h = (first + offset) % points.len();
        let from = points[h].uv;
        let mut order: Vec<usize> = (0..ring.len()).collect();
        order.sort_by(|&a, &b| {
            metric
                .distance(ring[a].uv, from)
                .total_cmp(&metric.distance(ring[b].uv, from))
        });
        // A bridge through a boundary point would pass that point by.
        let blocked = |target: [f64; 2], path: &[Corner]| {
            (0..path.len()).any(|k| {
                crosses(from, target, path[k].uv, path[(k + 1) % path.len()].uv)
                    || between(from, target, path[k].uv, line)
            })
        };
        let found = order.into_iter().find(|&r| {
            let target = ring[r].uv;
            // Among copies of one ring point, the bridge leaves from the
            // copy whose corner opens towards the hole.
            let before = ring[(r + ring.len() - 1) % ring.len()].uv;
            let after = ring[(r + 1) % ring.len()].uv;
            let (left, right) = (orient(before, target, from), orient(target, after, from));
            let opens = if orient(before, target, after) >= 0. {
                left > 0. && right > 0.
            } else {
                left >= 0. || right >= 0.
            };
            opens
                && !ring[r].pole
                && !blocked(target, ring)
                && !holes[index..].iter().any(|hole| blocked(target, hole))
        });
        if let Some(r) = found {
            let mut merged = ring[..=r].to_vec();
            merged.extend_from_slice(&points[h..]);
            merged.extend_from_slice(&points[..=h]);
            merged.extend_from_slice(&ring[r..]);
            *ring = merged;
            return Ok(());
        }
    }
    Err(Error::Unsupported)
}

fn ring(
    loops: &[UvLoop],
    pole: &Option<(f64, Key)>,
    surface: &Surface,
    metric: &Metric,
) -> Result<Vec<Corner>> {
    let periods = surface.periods();
    let (mut ring, holes) = if let Some((v, key)) = pole {
        let mut ring = unrolled(&loops[0], 0, periods, [0.; 2]);
        let (start, end) = (ring[0].uv[0], ring[ring.len() - 1].uv[0]);
        ring.extend([end, start].map(|u| Corner {
            uv: [u, *v],
            key: *key,
            pole: true,
            copy: 0,
        }));
        (ring, &loops[1..])
    } else if loops[0].winding != [0, 0] {
        (band(&loops[0], &loops[1], periods, metric), &loops[2..])
    } else {
        let outer = loops[0].points.iter().map(|p| corner(p, p.uv)).collect();
        (outer, &loops[1..])
    };
    let mut placed = holes
        .iter()
        .map(|hole| place(&ring, hole, periods))
        .collect::<Result<Vec<_>>>()?;
    // Rightmost first: a later hole can then reach one already joined.
    placed.sort_by(|a, b| b[rightmost(b)].uv[0].total_cmp(&a[rightmost(a)].uv[0]));
    for index in 0..placed.len() {
        bridge(&mut ring, &placed, index, metric)?;
    }
    Ok(ring)
}

/// Cut every pole row into short pieces. All of them are the same point
/// of the surface, and each then fans out to the ring points near its own
/// parameter instead of one copy reaching around the pole.
fn rows(ring: Vec<Corner>, period: Option<f64>) -> Vec<Corner> {
    let mut cut = Vec::with_capacity(ring.len());
    for (i, &corner) in ring.iter().enumerate() {
        cut.push(corner);
        let next = ring[(i + 1) % ring.len()];
        if corner.pole && next.pole && corner.key == next.key {
            let span = next.uv[0] - corner.uv[0];
            let pieces = (span.abs() / period.map_or(span.abs() / 8., |p| p / 24.))
                .ceil()
                .max(1.) as usize;
            for k in 1..pieces {
                let u = corner.uv[0] + span * k as f64 / pieces as f64;
                cut.push(Corner {
                    uv: [u, corner.uv[1]],
                    ..corner
                });
            }
        }
    }
    cut
}

struct Clip<'a> {
    ring: &'a [Corner],
    next: Vec<usize>,
    previous: Vec<usize>,
    active: Vec<bool>,
    /// End vertices of every edge, on a curved face.
    edges: Option<&'a [[usize; 2]]>,
    /// Parameters of each ring side and laid diagonal, by ordered key pair.
    sides: BTreeMap<(Key, Key), Vec<[[f64; 2]; 2]>>,
    metric: &'a Metric,
    epsilon: f64,
    /// How far from a diagonal a ring point still lies on it.
    line: f64,
}

/// The ordered key pair of a side with its parameters in the same order.
fn side(a: Corner, b: Corner) -> ((Key, Key), [[f64; 2]; 2]) {
    if a.key <= b.key {
        ((a.key, b.key), [a.uv, b.uv])
    } else {
        ((b.key, a.key), [b.uv, a.uv])
    }
}

impl Clip<'_> {
    /// Whether all three corners lie on one edge, with the middle one
    /// inside it. On a curved face such a triangle is a chord fan of the
    /// edge: it lies off the surface and coincides with the neighbouring
    /// face's triangle where two faces meet at a thin wedge, so it is taken
    /// only when no other ear remains. A plane contains its chord fans.
    fn chord(&self, a: Key, b: Key, c: Key) -> bool {
        let (Some(ends), Key::Sample(edge, _)) = (self.edges, b) else {
            return false;
        };
        let on = |key: Key| match key {
            Key::Sample(other, _) => other == edge,
            Key::Vertex(vertex) => ends[edge].contains(&vertex),
            Key::Extra(..) => false,
        };
        on(a) && on(c)
    }

    /// Whether a side or diagonal joins the same two boundary points
    /// elsewhere in the parameter plane, as across a seam. Another diagonal
    /// there would lay further triangles on the same chord.
    fn copies(&self, a: Corner, c: Corner) -> bool {
        let (keys, uv) = side(a, c);
        self.sides
            .get(&keys)
            .is_some_and(|sides| sides.iter().any(|&other| other != uv))
    }

    /// The length of the diagonal that clipping `b` would add, or `None`
    /// when `b` is not an admissible ear.
    fn ear(&self, b: usize) -> Option<(bool, f64)> {
        let (a, c) = (self.previous[b], self.next[b]);
        let [pa, pb, pc] = [a, b, c].map(|i| self.ring[i]);
        let twins = |p: Corner, q: Corner| p.key == q.key && !(p.pole && q.pole);
        let row = |p: Corner, q: Corner| p.key == q.key && p.pole && q.pole;
        if orient(pa.uv, pb.uv, pc.uv) <= self.epsilon
            || twins(pa, pb)
            || twins(pb, pc)
            || twins(pa, pc)
            || (self.copies(pa, pc) && !(pa.pole && pc.pole) && !row(pa, pb) && !row(pb, pc))
        {
            return None;
        }
        let mut i = c;
        while i != a {
            // A side through the diagonal shows a ring pinched at a shared
            // point, where the triangle would leave the region.
            if crosses(pa.uv, pc.uv, self.ring[i].uv, self.ring[self.next[i]].uv) {
                return None;
            }
            let p = self.ring[i].uv;
            // A ring point on the diagonal would be a corner of the
            // neighbouring face but not of this triangle's side.
            let row = pa.pole && pc.pole;
            if !row && !self.ring[i].pole && between(pa.uv, pc.uv, p, self.line) {
                return None;
            }
            let outside = p == pa.uv
                || p == pb.uv
                || p == pc.uv
                || orient(pa.uv, pb.uv, p) < 0.
                || orient(pb.uv, pc.uv, p) < 0.
                || orient(pc.uv, pa.uv, p) < 0.;
            if !outside {
                return None;
            }
            i = self.next[i];
        }
        let chord = self.chord(pa.key, pb.key, pc.key);
        Some((chord, self.metric.distance(pa.uv, pc.uv)))
    }
}

/// Triangles covering the face's parameter region, counter-clockwise.
/// `edges` holds the end vertices of every edge of the boundary.
pub(crate) fn triangulate(
    domain: &Domain,
    surface: &Surface,
    face: usize,
    edges: &[[usize; 2]],
) -> Result<Vec<Triangle>> {
    let (loops, pole) = match domain {
        Domain::Complete(range) => {
            return Ok(complete(face, *range, surface.periods()[1].is_some()));
        }
        Domain::Trimmed { loops, pole } => (loops, pole),
    };
    let metric = Metric(scale(loops, surface));
    let mut ring = rows(ring(loops, pole, surface, &metric)?, surface.periods()[0]);
    // Copies of a pole whole periods apart along u are one point of its
    // row: where a lone loop closes along its pole without a seam edge, the
    // seam sides at both ends of the row then share their corners, as the
    // seam sides of a complete surface do.
    let period = surface.periods()[0];
    let copies: Vec<usize> = ring
        .iter()
        .enumerate()
        .map(|(index, corner)| {
            let same = |other: &Corner| {
                other.key == corner.key
                    && other.pole
                    && period.is_some_and(|period| {
                        let turns = (other.uv[0] - corner.uv[0]) / period;
                        (turns - turns.round()).abs() <= 1e-9
                    })
            };
            if corner.pole {
                ring.iter().position(same).unwrap_or(index)
            } else {
                0
            }
        })
        .collect();
    for (corner, copy) in ring.iter_mut().zip(copies) {
        corner.copy = copy;
    }
    let n = ring.len();
    ensure(n >= 3, Error::Invalid)?;
    let extent = ring
        .iter()
        .fold(0., |e: f64, c| e.max(metric.distance(c.uv, ring[0].uv)));
    let mut clip = Clip {
        ring: &ring,
        edges: (!surface.is_plane()).then_some(edges),
        next: (0..n).map(|i| (i + 1) % n).collect(),
        previous: (0..n).map(|i| (i + n - 1) % n).collect(),
        active: vec![true; n],
        sides: BTreeMap::new(),
        metric: &metric,
        epsilon: 1e-13 * extent * extent / (metric.0[0] * metric.0[1]),
        line: 1e-9 * spread(&ring),
    };
    for i in 0..n {
        let (keys, uv) = side(ring[i], ring[(i + 1) % n]);
        clip.sides.entry(keys).or_default().push(uv);
    }
    let mut ears: Vec<Option<(bool, f64)>> = (0..n).map(|b| clip.ear(b)).collect();
    let mut triangles = Vec::with_capacity(n - 2);
    // Twice the area still to cover. A ring pinched at a shared boundary
    // point leaves slivers without area, whose sides are already covered.
    let mut left: f64 = (0..n)
        .map(|i| orient(ring[0].uv, ring[i].uv, ring[(i + 1) % n].uv))
        .sum();
    let enough = clip.epsilon.max(1e-9 * left);
    let mut count = n;
    while count > 3 && left > enough {
        let mut best: Option<((bool, f64), usize)> = None;
        for (b, ear) in ears.iter().enumerate() {
            if let (true, Some(cost)) = (clip.active[b], *ear)
                && best.is_none_or(|(least, _)| {
                    (cost.0, cost.1).partial_cmp(&(least.0, least.1))
                        == Some(std::cmp::Ordering::Less)
                })
            {
                best = Some((cost, b));
            }
        }
        let (_, b) = best.ok_or(Error::Unsupported)?;
        // A diagonal laid since this ear was judged may now forbid it.
        ears[b] = clip.ear(b);
        if ears[b].is_none() {
            continue;
        }
        let (a, c) = (clip.previous[b], clip.next[b]);
        let (keys, uv) = side(ring[a], ring[c]);
        clip.sides.entry(keys).or_default().push(uv);
        triangles.push([ring[a], ring[b], ring[c]]);
        left -= orient(ring[a].uv, ring[b].uv, ring[c].uv);
        count -= 1;
        clip.active[b] = false;
        clip.next[a] = c;
        clip.previous[c] = a;
        ears[a] = clip.ear(a);
        ears[c] = clip.ear(c);
    }
    if left > enough {
        let a = clip.active.iter().position(|&active| active).unwrap_or(0);
        let [b, c] = [clip.next[a], clip.next[clip.next[a]]];
        let twins = |p: Corner, q: Corner| p.key == q.key && !(p.pole && q.pole);
        ensure(
            !twins(ring[a], ring[b]) && !twins(ring[b], ring[c]) && !twins(ring[c], ring[a]),
            Error::Unsupported,
        )?;
        triangles.push([ring[a], ring[b], ring[c]]);
    }
    Ok(triangles)
}

/// Lengths of the parameter directions near the region, for comparing
/// distances across u and v.
fn scale(loops: &[UvLoop], surface: &Surface) -> [f64; 2] {
    let uv = loops[0]
        .points
        .iter()
        .find(|point| !point.pole)
        .map_or([0.; 2], |point| point.uv);
    let (_, [su, sv]) = surface.evaluate(uv);
    [norm(su).max(1e-9), norm(sv).max(1e-9)]
}

/// Whether the parameter point, or a periodic copy, lies in the region.
pub(crate) fn contains(triangles: &[Triangle], uv: [f64; 2], periods: [Option<f64>; 2]) -> bool {
    let turns = |axis: usize| {
        if periods[axis].is_some() {
            -2..=2
        } else {
            0..=0
        }
    };
    turns(0).any(|ku| {
        turns(1).any(|kv| {
            let p = [
                uv[0] + f64::from(ku) * periods[0].unwrap_or(0.),
                uv[1] + f64::from(kv) * periods[1].unwrap_or(0.),
            ];
            triangles.iter().any(|t| {
                orient(t[0].uv, t[1].uv, p) >= 0.
                    && orient(t[1].uv, t[2].uv, p) >= 0.
                    && orient(t[2].uv, t[0].uv, p) >= 0.
            })
        })
    })
}
