//! Plane sections of the retained boundary.
//!
//! The section follows the triangles of a tessellation only to know which
//! boundary pieces the plane crosses. Every published point is then solved
//! on the exact edge curve or surface, so neighbouring faces meet in the
//! same point and the loops close.
use super::mesh::{self, Mesh};
use super::{Interrupt, Point, Result, Solid, cross, dot, norm, sub};
use std::collections::BTreeMap;

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
    let mesh = mesh::tessellate(solid, options, interrupt)?;
    let above: Vec<bool> = mesh
        .positions
        .iter()
        .map(|&p| dot(sub(p, origin), normal) >= 0.)
        .collect();
    let mut points = BTreeMap::<(u32, u32), Point>::new();
    // Each crossed triangle joins the side it is entered by to the side it
    // is left by.
    let mut next = BTreeMap::<(u32, u32), (u32, u32)>::new();
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
            next.insert(from, to);
        }
    }
    let mut loops = Vec::new();
    while let Some((&start, _)) = next.iter().next() {
        let mut ring: Vec<Point> = Vec::new();
        let mut key = start;
        while let Some(following) = next.remove(&key) {
            ring.push(points[&key]);
            key = following;
        }
        let ring = simplified(ring, MERGE * options.deflection_mm);
        if ring.len() >= 3 {
            loops.push(ring);
        }
    }
    Ok(loops)
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
