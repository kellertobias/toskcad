//! Trim loops unrolled into the parameter plane of their surface.
//!
//! Loops are stored counter-clockwise about the surface normal, whatever
//! the face sense. Periodic parameters are continuous along a loop; a loop
//! that closes one period away from its start winds about the surface.
use super::{Boundary, Error, Point, Result, Surface, distance, ensure, nearest};

/// Identity of a boundary point shared between the faces that meet there.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Key {
    Vertex(usize),
    /// Interior sample `1` of edge `0`.
    Sample(usize, usize),
    /// A point private to face `0`.
    Extra(usize, usize),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct LoopPoint {
    pub uv: [f64; 2],
    /// The u parameter is free here: the surface collapses to a point.
    pub pole: bool,
    pub key: Key,
}

/// The piece of an edge between two consecutive loop points.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Segment {
    pub edge: usize,
    pub t: [f64; 2],
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct UvLoop {
    pub points: Vec<LoopPoint>,
    /// `segments[i]` joins point `i` to point `i + 1`; `None` is the
    /// zero-length row across a pole.
    pub segments: Vec<Option<Segment>>,
    /// Whole periods gained along u and v by one traversal.
    pub winding: [i32; 2],
}

impl UvLoop {
    /// Parameters of point `index`, continued past the end of the loop.
    pub fn at(&self, index: usize, periods: [Option<f64>; 2]) -> [f64; 2] {
        let turns = (index / self.points.len()) as f64;
        let uv = self.points[index % self.points.len()].uv;
        [0, 1].map(|axis| {
            uv[axis] + turns * f64::from(self.winding[axis]) * periods[axis].unwrap_or(0.)
        })
    }

    pub fn area(&self) -> f64 {
        let n = self.points.len();
        (0..n)
            .map(|i| {
                let (a, b) = (self.points[i].uv, self.points[(i + 1) % n].uv);
                0.5 * (a[0] * b[1] - b[0] * a[1])
            })
            .sum()
    }

    fn shift(&mut self, axis: usize, by: f64) {
        for point in &mut self.points {
            point.uv[axis] += by;
        }
    }

    fn mean(&self, axis: usize) -> f64 {
        self.points.iter().map(|p| p.uv[axis]).sum::<f64>() / self.points.len() as f64
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Domain {
    /// The whole parameter rectangle of a closed surface.
    Complete([[f64; 2]; 2]),
    Trimmed {
        /// The outer loop, or the winding pair lower/right loop first;
        /// hole loops follow.
        loops: Vec<UvLoop>,
        /// A lone winding loop closes at this v, optionally at a vertex.
        pole: Option<(f64, Key)>,
    },
}

pub(crate) struct Raw {
    pub position: Point,
    key: Key,
    segment: Segment,
}

/// The loop's boundary points in counter-clockwise order.
pub(crate) fn raw(
    boundary: &Boundary,
    uses: &[(usize, bool)],
    samples: &[Vec<f64>],
    reverse: bool,
) -> Vec<Raw> {
    let mut points = Vec::new();
    for &(index, forward) in uses {
        let edge = &boundary.edges[index];
        let parameters = &samples[index];
        let last = parameters.len() - 1;
        for step in 0..last {
            let (k, next) = if forward {
                (step, step + 1)
            } else {
                (last - step, last - step - 1)
            };
            let key = match k {
                0 => Key::Vertex(edge.vertices[0]),
                k if k == last => Key::Vertex(edge.vertices[1]),
                k => Key::Sample(index, k),
            };
            let position = match key {
                Key::Vertex(vertex) => boundary.vertices[vertex],
                _ => edge.curve.evaluate(parameters[k]).0,
            };
            points.push(Raw {
                position,
                key,
                segment: Segment {
                    edge: index,
                    t: [parameters[k], parameters[next]],
                },
            });
        }
    }
    if reverse {
        // Point i keeps its identity; its outgoing segment becomes the
        // reversed segment that used to arrive at it.
        let n = points.len();
        let segments: Vec<_> = points.iter().map(|p| p.segment).collect();
        points.reverse();
        for (i, point) in points.iter_mut().enumerate() {
            let old = segments[(2 * n - 2 - i) % n];
            point.segment = Segment {
                edge: old.edge,
                t: [old.t[1], old.t[0]],
            };
        }
    }
    points
}

/// The axes along which `uv` sits on a seam of the surface.
fn seams(surface: &Surface, uv: [f64; 2]) -> [bool; 2] {
    let periods = surface.periods();
    [0, 1].map(|axis| {
        periods[axis].is_some_and(|period| {
            let turn = (uv[axis] / period).rem_euclid(1.);
            turn.min(1. - turn) <= 1e-9
        })
    })
}

/// Unroll one loop. `None` means every point lies on a seam or a pole, so
/// the loop carries no information beyond the closed surface itself.
pub(crate) fn unroll(
    surface: &Surface,
    mut points: Vec<Raw>,
    tolerance: f64,
) -> Result<Option<UvLoop>> {
    let periods = surface.periods();
    let mut free = Vec::with_capacity(points.len());
    for point in &points {
        let uv = surface.invert(point.position, None);
        ensure(
            distance(surface.evaluate(uv).0, point.position) <= tolerance,
            Error::Invalid,
        )?;
        let pole = surface.is_pole(uv, tolerance);
        let seam = seams(surface, relative(surface, uv));
        free.push((uv, pole, pole || seam[0] || seam[1], seam));
    }
    // A loop that runs along one seam alone, such as a meridian closing a
    // band of a torus, is anchored at a seam point: the band's partner loop
    // places it. A loop on both seams or only at poles, or one that winds
    // nowhere, describes nothing beyond the closed surface itself.
    let along = |axis: usize| {
        free.iter().all(|state| state.1 || state.3[axis])
            && free.iter().any(|state| !state.1 && !state.3[1 - axis])
    };
    let anchor = free.iter().position(|state| !state.2).or_else(|| {
        (along(0) || along(1))
            .then(|| free.iter().position(|state| !state.1))
            .flatten()
    });
    let Some(anchor) = anchor else {
        return Ok(None);
    };
    let on_seams = free.iter().all(|state| state.2);
    points.rotate_left(anchor);
    free.rotate_left(anchor);
    let mut unrolled = UvLoop {
        points: vec![],
        segments: vec![],
        winding: [0; 2],
    };
    // Start indices of the chains that follow a pole.
    let mut chains = vec![];
    let mut previous = free[0].0;
    for (i, point) in points.iter().enumerate() {
        let (uv, pole, _, _) = free[i];
        let after_pole = unrolled.points.last().is_some_and(|last| last.pole);
        let uv = if i == 0 || after_pole {
            uv
        } else {
            surface.invert(point.position, Some(previous))
        };
        if pole {
            let uv = [previous[0], uv[1]];
            unrolled.points.push(LoopPoint {
                uv,
                pole,
                key: point.key,
            });
            unrolled.segments.push(None);
            chains.push(unrolled.points.len());
            unrolled.points.push(LoopPoint {
                uv,
                pole,
                key: point.key,
            });
        } else {
            unrolled.points.push(LoopPoint {
                uv,
                pole,
                key: point.key,
            });
        }
        unrolled.segments.push(Some(point.segment));
        previous = uv;
    }
    ensure(chains.len() <= 2, Error::Unsupported)?;
    let start = unrolled.points[0].uv;
    let closing = surface.invert(points[0].position, Some(previous));
    let turns = [0, 1]
        .map(|axis| periods[axis].map_or(0., |p| ((closing[axis] - start[axis]) / p).round()));
    // A seam loop that returns to its start encloses the whole surface,
    // like a meridian seam from pole to pole.
    if on_seams && turns == [0., 0.] {
        return Ok(None);
    }
    let Some(&last) = chains.last() else {
        unrolled.winding = turns.map(|t| t as i32);
        return Ok(Some(unrolled));
    };
    ensure(turns[1] == 0., Error::Unsupported)?;
    let period = periods[0].unwrap_or(0.);
    // The chain after the last pole must close onto the anchor.
    for point in &mut unrolled.points[last..] {
        point.uv[0] -= turns[0] * period;
    }
    let count = unrolled.points.len();
    let follow = |points: &[LoopPoint], pole: usize| points[(pole + 1) % count].uv[0];
    unrolled.points[last].uv[0] = if last + 1 == count {
        start[0]
    } else {
        follow(&unrolled.points, last)
    };
    if chains.len() == 2 {
        // A chain between two poles is the right side of the face when it
        // climbs and the left side when it descends.
        let (first, end) = (chains[0], last - 1);
        let reference = unrolled.points[first - 1].uv[0];
        let climbs = unrolled.points[end].uv[1] > unrolled.points[first].uv[1];
        let lead = unrolled.points[first + 1].uv[0];
        let offset = (lead - reference).rem_euclid(period);
        // Both meridians at one u would bound a lune without width; a face
        // bounded by one meridian alone is its seam and covers the surface.
        ensure(offset > 1e-9 * period, Error::Invalid)?;
        let offset = if climbs { offset } else { offset - period };
        let by = reference + offset - lead;
        for point in &mut unrolled.points[first + 1..=end] {
            point.uv[0] += by;
        }
        unrolled.points[first].uv[0] = lead + by;
    }
    Ok(Some(unrolled))
}

/// Parameters measured from the surface's seam.
fn relative(surface: &Surface, uv: [f64; 2]) -> [f64; 2] {
    match surface {
        Surface::Nurbs { surface, .. } => [0, 1].map(|axis| uv[axis] - surface.range(axis)[0]),
        Surface::Extrusion {
            curve: super::Curve::Nurbs { curve, .. },
            ..
        } => [uv[0] - curve.range()[0], uv[1]],
        Surface::Revolution {
            curve: super::Curve::Nurbs { curve, .. },
            ..
        } => [uv[0], uv[1] - curve.range()[0]],
        _ => uv,
    }
}

/// The parameter rectangle of a closed surface whose face is bounded only
/// by its own seams: edges the face uses in both directions.
fn whole(face: &super::Face) -> Option<[[f64; 2]; 2]> {
    let seams = face.loops.iter().all(|trim| {
        trim.uses
            .iter()
            .all(|&(edge, forward)| trim.uses.contains(&(edge, !forward)))
    });
    face.surface
        .complete()
        .filter(|_| seams && !face.loops.is_empty())
}

/// Unroll every loop of a face and arrange them as one region.
pub(crate) fn build(
    boundary: &Boundary,
    face_index: usize,
    samples: &[Vec<f64>],
    tolerance: f64,
) -> Result<Domain> {
    let face = &boundary.faces[face_index];
    let surface = &face.surface;
    let periods = surface.periods();
    if let Some(range) = whole(face) {
        return Ok(Domain::Complete(range));
    }
    let mut loops = Vec::new();
    let mut unanchored = 0;
    for trim in face.loops.iter().filter(|trim| !trim.uses.is_empty()) {
        let points = raw(boundary, &trim.uses, samples, !face.same_sense);
        match unroll(surface, points, tolerance)? {
            Some(unrolled) => loops.push(unrolled),
            None => unanchored += 1,
        }
    }
    if loops.is_empty() {
        ensure(unanchored > 0 || !face.loops.is_empty(), Error::Invalid)?;
        return surface
            .complete()
            .map(Domain::Complete)
            .ok_or(Error::Invalid);
    }
    ensure(unanchored == 0, Error::Unsupported)?;
    ensure(
        loops
            .iter()
            .all(|l| l.winding.iter().all(|w| w.abs() <= 1) && l.winding.contains(&0)),
        Error::Unsupported,
    )?;
    let (mut winding, mut flat): (Vec<_>, Vec<_>) =
        loops.into_iter().partition(|l| l.winding != [0, 0]);
    let (mut outer, holes): (Vec<_>, Vec<_>) = flat.drain(..).partition(|l| l.area() > 0.);
    let mut pole = None;
    match (winding.len(), outer.len()) {
        (0, 1) => {}
        (1, 0) => {
            let turn = winding[0].winding;
            ensure(turn[1] == 0, Error::Unsupported)?;
            let v = surface.pole(turn[0] > 0).ok_or(Error::Invalid)?;
            let apex = surface.evaluate([0., v]).0;
            let vertex = face
                .loops
                .iter()
                .filter_map(|trim| trim.vertex)
                .find(|&vertex| distance(boundary.vertices[vertex], apex) <= tolerance);
            pole = Some((v, vertex.map_or(Key::Extra(face_index, 0), Key::Vertex)));
        }
        (2, 0) => {
            let axis = usize::from(winding[0].winding[0] == 0);
            let turns = [winding[0].winding[axis], winding[1].winding[axis]];
            ensure(
                turns[0] + turns[1] == 0 && winding[1].winding[1 - axis] == 0,
                Error::Invalid,
            )?;
            if turns[0] < 0 {
                winding.swap(0, 1);
            }
            // Along u the forward loop lies below its partner; along v it
            // lies to the right of it.
            let across = 1 - axis;
            let (lead, partner) = if axis == 0 { (0, 1) } else { (1, 0) };
            let gap = winding[partner].mean(across) - winding[lead].mean(across);
            let by = nearest(gap, periods[across].unwrap_or(0.) * 0.5, periods[across]) - gap;
            ensure(gap + by > 0., Error::Invalid)?;
            winding[partner].shift(across, by);
            if axis == 0 {
                // Both starts share one period window for the mass integral.
                let from = winding[0].points[0].uv[0];
                let other = winding[1].points[0].uv[0];
                let period = periods[0].unwrap_or(0.);
                winding[1].shift(0, from + (other - from).rem_euclid(period) - other);
            }
        }
        _ => return Err(Error::Invalid),
    }
    winding.append(&mut outer);
    winding.extend(holes);
    Ok(Domain::Trimmed {
        loops: winding,
        pole,
    })
}

/// Whether the face's loops, as stored, run counter-clockwise about the
/// surface normal. Writers disagree on whether the loops of a reversed face
/// run about the face normal or the surface normal, so the direction is
/// read from the geometry: the outer loop encloses a finite region. `None`
/// means the loops alone do not decide it.
pub(crate) fn orientation(
    boundary: &Boundary,
    face_index: usize,
    samples: &[Vec<f64>],
    tolerance: f64,
) -> Result<Option<bool>> {
    let face = &boundary.faces[face_index];
    if whole(face).is_some() {
        return Ok(None);
    }
    let mut outer = 0.;
    let mut winding = Vec::new();
    for trim in face.loops.iter().filter(|trim| !trim.uses.is_empty()) {
        let points = raw(boundary, &trim.uses, samples, false);
        let Some(unrolled) = unroll(&face.surface, points, tolerance)? else {
            continue;
        };
        if unrolled.winding == [0, 0] {
            let area = unrolled.area();
            if area.abs() > f64::abs(outer) {
                outer = area;
            }
        } else {
            winding.push(unrolled);
        }
    }
    if outer != 0. {
        return Ok(Some(outer > 0.));
    }
    let above = |low: &UvLoop, high: &UvLoop| low.winding[0] > 0 && low.mean(1) < high.mean(1);
    Ok(match winding.as_slice() {
        [a, b] if a.winding[0] + b.winding[0] == 0 && face.surface.periods()[1].is_none() => {
            Some(above(a, b) || above(b, a))
        }
        _ => None,
    })
}
