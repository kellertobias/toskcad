//! Plane sections of the retained boundary.
//!
//! The section follows the triangles of a tessellation only to know which
//! boundary pieces the plane crosses. Every published point is then solved
//! on the exact edge curve or surface, so neighbouring faces meet in the
//! same point and the loops close.
use super::mesh::{self, Mesh};
use super::{
    Interrupt, Point, Result, Solid, Surface, add, cross, distance, dot, norm, scale, sub,
};
use crate::trim_curve::{Ellipse, EllipseArc};
use crate::warped_union::{Curve as SectionCurve, SectionEdge};
use std::collections::BTreeMap;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

/// The parameter in `[0, 1]` where `height` changes sign, by bisection.
/// `above` is the sign class at zero, taken from the mesh vertex so that
/// every triangle agrees on which sides the plane crosses.
fn crossing(above: bool, height: impl Fn(f64) -> f64) -> f64 {
    let (mut low, mut high) = (0., 1.);
    for _ in 0..60 {
        let middle = 0.5 * (low + high);
        if (height(middle) >= 0.) == above {
            low = middle;
        } else {
            high = middle;
        }
    }
    0.5 * (low + high)
}

/// The point where side `i` of triangle `index` meets the plane.
fn point(
    solid: &Solid,
    mesh: &Mesh,
    index: usize,
    i: usize,
    above: bool,
    plane: (Point, Point),
) -> Point {
    let height = |p: Point| dot(sub(p, plane.0), plane.1);
    match mesh.sides[index][i] {
        Some(piece) => {
            let curve = &solid.boundary().edges[piece.edge].curve;
            let at = |s: f64| curve.evaluate(piece.t[0] + (piece.t[1] - piece.t[0]) * s).0;
            at(crossing(above, |s| height(at(s))))
        }
        None => {
            let surface = &solid.boundary().faces[mesh.faces[index] as usize].surface;
            let at = |s: f64| surface.evaluate(mesh.cells[index].along(i, s)).0;
            at(crossing(above, |s| height(at(s))))
        }
    }
}

/// Section points closer than this fraction of the deflection are one.
const MERGE: f64 = 1e-3;

/// Closed loops where the plane through `origin` with unit `normal` cuts
/// the solid. Seen against the normal, material lies to the left of each
/// loop. Chords between consecutive points stay within the tessellation
/// deflection of the exact section curve, plus the merge distance.
pub fn section(
    solid: &Solid,
    origin: Point,
    normal: Point,
    options: mesh::Options,
    interrupt: Interrupt<'_>,
) -> Result<Vec<Vec<Point>>> {
    Ok(crossings(solid, origin, normal, options, interrupt)?
        .into_iter()
        .map(|ring| {
            simplified(
                ring.into_iter().map(|(p, _)| p).collect(),
                MERGE * options.deflection_mm,
            )
        })
        .filter(|ring| ring.len() >= 3)
        .collect())
}

/// Each loop's crossing points with the face whose triangle the chord
/// leaving the point runs through.
fn crossings(
    solid: &Solid,
    origin: Point,
    normal: Point,
    options: mesh::Options,
    interrupt: Interrupt<'_>,
) -> Result<Vec<Vec<(Point, u32)>>> {
    let mesh = mesh::tessellate(solid, options, interrupt)?;
    let above: Vec<bool> = mesh
        .positions
        .iter()
        .map(|&p| dot(sub(p, origin), normal) >= 0.)
        .collect();
    let mut points = BTreeMap::<(u32, u32), Point>::new();
    // Each crossed triangle joins the side it is entered by to the side it
    // is left by.
    let mut next = BTreeMap::<(u32, u32), ((u32, u32), u32)>::new();
    for (index, triangle) in mesh.triangles.iter().enumerate() {
        let mut ends = [None; 2];
        for i in 0..3 {
            let (a, b) = (triangle[i], triangle[(i + 1) % 3]);
            if above[a as usize] != above[b as usize] {
                let key = (a.min(b), a.max(b));
                let found = point(solid, &mesh, index, i, above[a as usize], (origin, normal));
                points.entry(key).or_insert(found);
                // Material is left of the travel direction `normal x outward`:
                // the loop leaves through the side that climbs above the plane.
                ends[usize::from(above[b as usize])] = Some(key);
            }
        }
        if let [Some(from), Some(to)] = ends {
            next.insert(from, (to, mesh.faces[index]));
        }
    }
    let mut loops = Vec::new();
    while let Some((&start, _)) = next.iter().next() {
        let mut ring = Vec::new();
        let mut key = start;
        while let Some((following, face)) = next.remove(&key) {
            ring.push((points[&key], face));
            key = following;
        }
        loops.push(ring);
    }
    Ok(loops)
}

/// Section loops as exact curves where the plane meets a face in a line or
/// an ellipse: a plane, a cylinder, a cone or a sphere cut so. Elsewhere,
/// as on a cone cut in a parabola or hyperbola, a torus or a spline, the
/// loop keeps its chords, which stay within the tessellation deflection
/// plus the merge distance of the exact section curve. Consecutive edges
/// share their end points, and material lies to the left seen against the
/// normal.
pub fn section_edges(
    solid: &Solid,
    origin: Point,
    normal: Point,
    options: mesh::Options,
    interrupt: Interrupt<'_>,
) -> Result<Vec<Vec<SectionEdge>>> {
    let unit = scale(normal, 1. / norm(normal));
    let mut loops = Vec::new();
    for ring in crossings(solid, origin, normal, options, interrupt)? {
        let ring = merged(ring, MERGE * options.deflection_mm);
        if ring.len() < 3 {
            continue;
        }
        let n = ring.len();
        // Start at a change of face, so that every run lies on one face.
        let start = (0..n)
            .find(|&i| ring[(i + n - 1) % n].1 != ring[i].1)
            .unwrap_or(0);
        let closed = (0..n).all(|i| ring[i].1 == ring[0].1);
        let mut edges = Vec::new();
        let mut i = 0;
        while i < n {
            let face = ring[(start + i) % n].1;
            let mut run = vec![ring[(start + i) % n].0];
            while i < n && ring[(start + i) % n].1 == face {
                i += 1;
                run.push(ring[(start + i) % n].0);
            }
            let surface = &solid.boundary().faces[face as usize].surface;
            edges.extend(exact(
                surface,
                &run,
                closed,
                (origin, unit),
                solid.tolerance(),
            ));
        }
        loops.push(edges);
    }
    Ok(loops)
}

/// Drop points within `merge` of the point kept before them; the kept
/// point keeps its own face.
fn merged(ring: Vec<(Point, u32)>, merge: f64) -> Vec<(Point, u32)> {
    let mut kept: Vec<(Point, u32)> = Vec::with_capacity(ring.len());
    for (point, face) in ring {
        if kept
            .last()
            .is_none_or(|&(last, _)| norm(sub(point, last)) > merge)
        {
            kept.push((point, face));
        }
    }
    while kept.len() > 1 && norm(sub(kept[kept.len() - 1].0, kept[0].0)) <= merge {
        kept.pop();
    }
    kept
}

fn lines(run: &[Point]) -> Vec<SectionEdge> {
    run.windows(2)
        .map(|pair| SectionEdge {
            ends_mm: [pair[0], pair[1]],
            curve: SectionCurve::Line,
        })
        .collect()
}

/// The exact curve through the points of one run on `surface`, or its
/// chords when the plane meets the surface in another curve.
fn exact(
    surface: &Surface,
    run: &[Point],
    closed: bool,
    plane: (Point, Point),
    tolerance: f64,
) -> Vec<SectionEdge> {
    let ends = [run[0], run[run.len() - 1]];
    let near = |p: Point| tolerance + 1e-12 * (1. + norm(p));
    if !closed && distance(ends[0], ends[1]) > near(ends[0]) {
        // A straight run: every point on the chord joining its ends.
        let along = sub(ends[1], ends[0]);
        let straight = run.iter().all(|&p| {
            let t = dot(sub(p, ends[0]), along) / dot(along, along);
            (-1e-12..=1. + 1e-12).contains(&t)
                && distance(p, add(ends[0], scale(along, t))) <= near(p)
        });
        if straight {
            return vec![SectionEdge {
                ends_mm: ends,
                curve: SectionCurve::Line,
            }];
        }
    }
    match conic(surface, plane).and_then(|ellipse| arc(ellipse, run, closed, near)) {
        Some(arc) => vec![SectionEdge {
            ends_mm: ends,
            curve: SectionCurve::Ellipse(arc),
        }],
        None => lines(run),
    }
}

/// The ellipse where the plane meets a sphere, a cylinder or one nappe of
/// a cone, from the surface's implicit quadric `(x - c)' M (x - c) = k`
/// restricted to the plane. `None` for any other section.
fn conic(surface: &Surface, (origin, normal): (Point, Point)) -> Option<Ellipse> {
    let outer = |z: Point, by: f64| -> [Point; 3] {
        std::array::from_fn(|i| {
            std::array::from_fn(|j| f64::from(u8::from(i == j)) - by * z[i] * z[j])
        })
    };
    let (centre, m, k) = match surface {
        Surface::Sphere { frame, radius } => (frame.origin, outer(frame.z, 0.), radius * radius),
        Surface::Cylinder { frame, radius } => (frame.origin, outer(frame.z, 1.), radius * radius),
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => {
            let slope = half_angle.tan();
            (
                sub(frame.origin, scale(frame.z, radius / slope)),
                outer(frame.z, 1. + slope * slope),
                0.,
            )
        }
        _ => return None,
    };
    let apply = |x: Point| -> Point { std::array::from_fn(|i| dot(m[i], x)) };
    let helper = if normal[0].abs() < 0.9 {
        [1., 0., 0.]
    } else {
        [0., 1., 0.]
    };
    let e1 = {
        let e = cross(normal, helper);
        scale(e, 1. / norm(e))
    };
    let e2 = cross(normal, e1);
    let w = sub(origin, centre);
    let (a, b, c) = (dot(e1, apply(e1)), dot(e1, apply(e2)), dot(e2, apply(e2)));
    let linear = [dot(e1, apply(w)), dot(e2, apply(w))];
    let constant = dot(w, apply(w)) - k;
    let determinant = a * c - b * b;
    if determinant <= 1e-12 * (a * a + c * c) {
        return None;
    }
    // Centre y0 = -A^-1 l; the conic is (y - y0)' A (y - y0) = h.
    let y0 = [
        -(c * linear[0] - b * linear[1]) / determinant,
        -(a * linear[1] - b * linear[0]) / determinant,
    ];
    let h = -(linear[0] * y0[0] + linear[1] * y0[1]) - constant;
    let mean = 0.5 * (a + c);
    let spread = (0.25 * (a - c) * (a - c) + b * b).sqrt();
    let (low, high) = (mean - spread, mean + spread);
    if low <= 0. || h <= 0. {
        return None;
    }
    // The major axis is the eigenvector of the smaller eigenvalue.
    let direction = if spread <= 1e-12 * mean {
        [1., 0.]
    } else if (a - low).abs() >= (c - low).abs() {
        [b, low - a]
    } else {
        [low - c, b]
    };
    let major = add(scale(e1, direction[0]), scale(e2, direction[1]));
    let centre = add(origin, add(scale(e1, y0[0]), scale(e2, y0[1])));
    Ellipse::new(centre, normal, major, [(h / low).sqrt(), (h / high).sqrt()]).ok()
}

/// The arc of `ellipse` through every point of the run, in run order.
fn arc(
    ellipse: Ellipse,
    run: &[Point],
    closed: bool,
    near: impl Fn(Point) -> f64,
) -> Option<EllipseArc> {
    let [major, minor] = ellipse.semiaxes_mm();
    let mut angles: Vec<f64> = Vec::with_capacity(run.len());
    for &p in run {
        let d = sub(p, ellipse.center_mm());
        let mut t =
            (dot(d, ellipse.minor_axis()) / minor).atan2(dot(d, ellipse.major_axis()) / major);
        if distance(ellipse.point_at(t), p) > near(p) {
            return None;
        }
        if let Some(&previous) = angles.last() {
            t = previous + (t - previous + PI).rem_euclid(TAU) - PI;
            // Consecutive crossings are a chord apart, far less than a
            // half turn: a larger step would leave the run ambiguous.
            if (t - previous).abs() >= FRAC_PI_2 {
                return None;
            }
        }
        angles.push(t);
    }
    let sweep = angles[angles.len() - 1] - angles[0];
    let sweep = if closed {
        if (sweep.abs() - TAU).abs() > 1e-6 {
            return None;
        }
        TAU.copysign(sweep)
    } else {
        sweep
    };
    // Every crossing moves the same way round the ellipse.
    let monotonic = angles
        .windows(2)
        .all(|pair| (pair[1] - pair[0]) * sweep > 0.);
    (monotonic && sweep.abs() <= TAU)
        .then(|| EllipseArc::new(ellipse, angles[0], sweep).ok())
        .flatten()
}

/// Drop points within `merge` of the point kept before them, and points
/// inside a straight run, around the loop. Merging moves the loop by at
/// most `merge`.
fn simplified(ring: Vec<Point>, merge: f64) -> Vec<Point> {
    let mut kept: Vec<Point> = Vec::with_capacity(ring.len());
    for point in ring {
        if kept
            .last()
            .is_none_or(|&last| norm(sub(point, last)) > merge)
        {
            kept.push(point);
        }
    }
    while kept.len() > 1 && norm(sub(kept[kept.len() - 1], kept[0])) <= merge {
        kept.pop();
    }
    let straight = |a: Point, b: Point, c: Point| {
        let (ab, bc) = (sub(b, a), sub(c, b));
        norm(cross(ab, bc)) <= 1e-12 * norm(ab) * norm(bc) && dot(ab, bc) > 0.
    };
    let mut changed = true;
    while changed && kept.len() > 3 {
        changed = false;
        let n = kept.len();
        if let Some(i) =
            (0..n).find(|&i| straight(kept[(i + n - 1) % n], kept[i], kept[(i + 1) % n]))
        {
            kept.remove(i);
            changed = true;
        }
    }
    kept
}

/// Length of a closed loop.
pub fn length(ring: &[Point]) -> f64 {
    (0..ring.len())
        .map(|i| norm(sub(ring[(i + 1) % ring.len()], ring[i])))
        .sum()
}
