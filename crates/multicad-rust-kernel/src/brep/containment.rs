//! Whole-box material classification against a retained boundary.
//!
//! A box is settled in two steps. First, every face's support surface is
//! proved apart from the box. Analytic supports use their exact
//! 1-Lipschitz distance functions at the box centre; spline supports use
//! the convex hulls of their Bézier pieces, halved until the hull boxes
//! clear the probe. The least of these lower bounds is the box's certified
//! clearance from the whole boundary, so the box lies wholly on one side.
//! Second, that side is the ray parity of the box centre against the
//! watertight tessellation, admitted only when the clearance and the centre's
//! distance from every triangle both exceed twice the tessellation's
//! deflection, the standard of the shell-nesting proof. Anything closer is
//! `Unresolved`.
use super::nurbs::Patch;
use super::{
    Curve, Interrupt, Point, Result, Solid, Surface, checkpoint, cross, crosses_ray, dot, mesh,
    norm, point_triangle, sub,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoxSide {
    Inside,
    Outside,
    Unresolved,
}

type Homogeneous = [f64; 4];

/// The most hull boxes one classification may visit.
const MAX_HULLS: usize = 1 << 20;
const MAX_DEPTH: usize = 40;

struct Probe {
    lo: Point,
    hi: Point,
    center: Point,
    /// Half the box diagonal.
    reach: f64,
}

fn gap(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[1]).max(b[0] - a[1]).max(0.)
}

impl Probe {
    /// Euclidean gap between the probe and an axis-aligned box.
    fn box_gap(&self, lo: Point, hi: Point) -> f64 {
        (0..3)
            .map(|i| gap([self.lo[i], self.hi[i]], [lo[i], hi[i]]))
            .fold(0., f64::hypot)
    }
}

fn projected(h: Homogeneous) -> Option<Point> {
    (h[3] > 0. && h.iter().all(|v| v.is_finite())).then(|| [h[0] / h[3], h[1] / h[3], h[2] / h[3]])
}

/// The axis-aligned box of a hull's projected controls, with a rounding
/// guard. Positive weights keep the curve inside that hull.
fn hull_box(control: &[Homogeneous]) -> Option<(Point, Point)> {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &h in control {
        let p = projected(h)?;
        for i in 0..3 {
            lo[i] = lo[i].min(p[i]);
            hi[i] = hi[i].max(p[i]);
        }
    }
    let scale = lo.iter().chain(&hi).map(|v| v.abs()).fold(1., f64::max);
    let guard = 64. * f64::EPSILON * scale;
    Some((lo.map(|v| v - guard), hi.map(|v| v + guard)))
}

/// de Casteljau halves of one homogeneous Bézier sequence.
fn halves(control: &[Homogeneous]) -> [Vec<Homogeneous>; 2] {
    let mut work = control.to_vec();
    let n = work.len();
    let mut left = Vec::with_capacity(n);
    let mut right = vec![[0.; 4]; n];
    left.push(work[0]);
    right[n - 1] = work[n - 1];
    for level in 1..n {
        for i in 0..n - level {
            work[i] = std::array::from_fn(|k| 0.5 * (work[i][k] + work[i + 1][k]));
        }
        left.push(work[0]);
        right[n - 1 - level] = work[n - 1 - level];
    }
    [left, right]
}

fn split_patch(patch: &Patch, along_u: bool) -> [Patch; 2] {
    let (rows, columns) = (patch.rows, patch.columns);
    let mut parts = [patch.clone(), patch.clone()];
    if along_u {
        for j in 0..columns {
            let line: Vec<_> = (0..rows).map(|i| patch.control[i * columns + j]).collect();
            let [a, b] = halves(&line);
            for i in 0..rows {
                parts[0].control[i * columns + j] = a[i];
                parts[1].control[i * columns + j] = b[i];
            }
        }
    } else {
        for i in 0..rows {
            let line = &patch.control[i * columns..(i + 1) * columns];
            let [a, b] = halves(line);
            parts[0].control[i * columns..(i + 1) * columns].copy_from_slice(&a);
            parts[1].control[i * columns..(i + 1) * columns].copy_from_slice(&b);
        }
    }
    parts
}

/// Spread of the projected net along u rows versus v columns.
fn longer_along_u(patch: &Patch) -> bool {
    let at =
        |i: usize, j: usize| projected(patch.control[i * patch.columns + j]).unwrap_or([0.; 3]);
    let mut du: f64 = 0.;
    let mut dv: f64 = 0.;
    for i in 0..patch.rows {
        for j in 0..patch.columns {
            if i + 1 < patch.rows {
                du = du.max(norm(sub(at(i + 1, j), at(i, j))));
            }
            if j + 1 < patch.columns {
                dv = dv.max(norm(sub(at(i, j + 1), at(i, j))));
            }
        }
    }
    du * (patch.rows as f64) >= dv * (patch.columns as f64)
}

/// A lower bound of the gap between the probe and a spline surface, or
/// `None` when the hulls still meet the probe within the margin.
fn patch_clearance(
    patches: Vec<Patch>,
    probe: &Probe,
    margin: f64,
    budget: &mut usize,
) -> Option<f64> {
    let mut least = f64::INFINITY;
    let mut stack: Vec<(Patch, usize)> = patches.into_iter().map(|p| (p, 0)).collect();
    while let Some((patch, depth)) = stack.pop() {
        *budget = budget.checked_sub(1)?;
        let (lo, hi) = hull_box(&patch.control)?;
        let distance = probe.box_gap(lo, hi);
        if distance > margin {
            least = least.min(distance);
            continue;
        }
        if depth >= MAX_DEPTH {
            return None;
        }
        let along_u = longer_along_u(&patch);
        for part in split_patch(&patch, along_u) {
            stack.push((part, depth + 1));
        }
    }
    Some(least)
}

/// An orthonormal pair across `direction`.
fn across(direction: Point) -> Option<[Point; 2]> {
    let length = norm(direction);
    if length.is_nan() || length <= 0. {
        return None;
    }
    let d = direction.map(|v| v / length);
    let seed = if d[0].abs() < 0.6 {
        [1., 0., 0.]
    } else {
        [0., 1., 0.]
    };
    let a = cross(d, seed);
    let a = a.map(|v| v / norm(a));
    Some([a, cross(d, a)])
}

/// A lower bound of the gap between the probe and the extrusion of a spline
/// curve along `direction`, compared in the plane across the direction,
/// onto which projection is 1-Lipschitz.
fn extrusion_clearance(
    segments: Vec<Vec<Homogeneous>>,
    direction: Point,
    probe: &Probe,
    margin: f64,
    budget: &mut usize,
) -> Option<f64> {
    let [e1, e2] = across(direction)?;
    let flat = |p: Point| [dot(p, e1), dot(p, e2)];
    let corners = (0..8).map(|k| {
        std::array::from_fn(|i| {
            if k >> i & 1 == 1 {
                probe.hi[i]
            } else {
                probe.lo[i]
            }
        })
    });
    let mut probe2 = [[f64::INFINITY, f64::NEG_INFINITY]; 2];
    for corner in corners {
        let f = flat(corner);
        for i in 0..2 {
            probe2[i] = [probe2[i][0].min(f[i]), probe2[i][1].max(f[i])];
        }
    }
    let mut least = f64::INFINITY;
    let mut stack: Vec<(Vec<Homogeneous>, usize)> = segments.into_iter().map(|s| (s, 0)).collect();
    while let Some((segment, depth)) = stack.pop() {
        *budget = budget.checked_sub(1)?;
        let mut range = [[f64::INFINITY, f64::NEG_INFINITY]; 2];
        let mut scale: f64 = 1.;
        for &h in &segment {
            let f = flat(projected(h)?);
            for i in 0..2 {
                range[i] = [range[i][0].min(f[i]), range[i][1].max(f[i])];
                scale = scale.max(f[i].abs());
            }
        }
        let guard = 64. * f64::EPSILON * scale;
        let distance = (0..2)
            .map(|i| gap(probe2[i], [range[i][0] - guard, range[i][1] + guard]))
            .fold(0., f64::hypot);
        if distance > margin {
            least = least.min(distance);
            continue;
        }
        if depth >= MAX_DEPTH {
            return None;
        }
        for half in halves(&segment) {
            stack.push((half, depth + 1));
        }
    }
    Some(least)
}

/// A lower bound of the gap between the probe and one support surface.
fn surface_clearance(
    surface: &Surface,
    probe: &Probe,
    margin: f64,
    budget: &mut usize,
) -> Option<f64> {
    let c = probe.center;
    let h = probe.reach;
    let lipschitz = |distance: f64| Some(distance - h);
    match surface {
        Surface::Plane(frame) => lipschitz(frame.local(c)[2].abs()),
        Surface::Cylinder { frame, radius } => {
            let l = frame.local(c);
            lipschitz((l[0].hypot(l[1]) - radius).abs())
        }
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => {
            // Both generating lines of the double cone in the meridian
            // half-plane; the surface lies on their union.
            let l = frame.local(c);
            let (rho, z) = (l[0].hypot(l[1]), l[2]);
            let t = half_angle.tan();
            let scale = t.hypot(1.);
            let first = (rho - radius - t * z).abs() / scale;
            let second = (rho + radius + t * z).abs() / scale;
            lipschitz(first.min(second))
        }
        Surface::Sphere { frame, radius } => lipschitz((norm(frame.local(c)) - radius).abs()),
        Surface::Torus {
            frame,
            major,
            minor,
        } => {
            let l = frame.local(c);
            lipschitz(((l[0].hypot(l[1]) - major).hypot(l[2]) - minor).abs())
        }
        Surface::Extrusion { curve, direction } => match curve {
            Curve::Line {
                origin,
                direction: along,
            } => {
                let n = cross(*along, *direction);
                let length = norm(n);
                if length.is_nan() || length <= 0. {
                    return None;
                }
                lipschitz((dot(sub(c, *origin), n) / length).abs())
            }
            Curve::Circle { frame, radius } => {
                let n = norm(*direction);
                if n.is_nan() || n <= 0. || (dot(*direction, frame.z).abs() / n - 1.).abs() > 1e-12
                {
                    return None;
                }
                let l = frame.local(c);
                lipschitz((l[0].hypot(l[1]) - radius).abs())
            }
            Curve::Nurbs { curve, .. } => {
                extrusion_clearance(curve.bezier()?, *direction, probe, margin, budget)
            }
            _ => None,
        },
        Surface::Nurbs { surface, .. } => patch_clearance(surface.bezier()?, probe, margin, budget),
        Surface::Revolution { .. } => None,
    }
}

impl Solid {
    /// A certified lower bound of the distance between the axis-aligned box
    /// `[lo, hi]` and every face support of the boundary, or `None` when a
    /// support comes within `margin` or cannot be bounded.
    pub fn box_clearance(
        &self,
        [lo, hi]: [Point; 2],
        margin: f64,
        interrupt: Interrupt<'_>,
    ) -> Result<Option<f64>> {
        if (0..3).any(|i| !(lo[i].is_finite() && hi[i].is_finite() && lo[i] <= hi[i]))
            || !(margin.is_finite() && margin >= 0.)
        {
            return Err(super::Error::Invalid);
        }
        let center = std::array::from_fn(|i| 0.5 * (lo[i] + hi[i]));
        let half: Point = std::array::from_fn(|i| 0.5 * (hi[i] - lo[i]));
        let scale = center
            .iter()
            .chain(&half)
            .map(|v| v.abs())
            .fold(1., f64::max);
        let probe = Probe {
            lo,
            hi,
            center,
            reach: norm(half) * (1. + 1e-12) + 64. * f64::EPSILON * scale,
        };
        let mut budget = MAX_HULLS;
        let mut least = f64::INFINITY;
        for face in &self.boundary().faces {
            checkpoint(interrupt)?;
            match surface_clearance(&face.surface, &probe, margin, &mut budget) {
                Some(distance) if distance > margin => least = least.min(distance),
                _ => {
                    if budget == 0 {
                        return Err(super::Error::Limit);
                    }
                    return Ok(None);
                }
            }
        }
        Ok(Some(least))
    }

    /// The side of a whole axis-aligned box whose certified `clearance`
    /// from every face support came from `box_clearance`. `tessellation`
    /// must tessellate this solid within `deflection_mm`; its ray parity is
    /// admitted only beyond twice that deflection, as in the shell-nesting
    /// proof.
    pub fn box_parity(
        &self,
        [lo, hi]: [Point; 2],
        clearance: f64,
        tessellation: &mesh::Mesh,
        deflection_mm: f64,
    ) -> BoxSide {
        let floor = 2. * deflection_mm.max(tessellation.maximum_deviation_mm);
        if clearance.is_nan() || clearance <= floor {
            return BoxSide::Unresolved;
        }
        let center = std::array::from_fn(|i| 0.5 * (lo[i] + hi[i]));
        let mut crossings = 0usize;
        for triangle in &tessellation.triangles {
            let t = triangle.map(|i| tessellation.positions[i as usize]);
            if point_triangle(center, t) <= floor {
                return BoxSide::Unresolved;
            }
            if crosses_ray(center, &t) {
                crossings += 1;
            }
        }
        if crossings % 2 == 1 {
            BoxSide::Inside
        } else {
            BoxSide::Outside
        }
    }
}

#[cfg(test)]
mod tests;
