//! Whole-probe relations against any retained warped shell, chart or not.
//!
//! A box probe whose enclosure clears every face by more than the
//! construction tolerance cannot meet the boundary, so all of it lies on one
//! side. A face clears it when its whole support (plane or full cylinder)
//! does, or, for a planar face, when the box of its outer loop does. Its centre is then classified by
//! the winding number of the shell's exact plane section through it. Section
//! curves lie on the supports, so they keep at least that clearance from the
//! centre; an elliptic arc is followed in steps whose sagitta stays below
//! half of it, and each step therefore sweeps the same angle as its chord.
//! A proved relation shares an unchanged operand; nothing is cut.
use super::*;
use crate::geometry::helical_body::contained_boolean::{
    chart_bounds, probe_bounds, whole_operand_result,
};
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use multicad_rust_kernel::helical_thread::BoxClassification;
use multicad_rust_kernel::warped_union::SectionEdge;
use std::f64::consts::TAU;

const SECTION_SEGMENTS: usize = 1 << 14;
const ARC_STEPS: usize = 1 << 16;

/// A conservative lower bound of the distance from the box to a cylinder
/// support, or `None` when the support may touch it.
fn cylinder_clearance(support: Support, local: [B; 3], margin: f64) -> Option<f64> {
    match support {
        Support::Plane { .. } => None,
        Support::Cylinder {
            axis_origin_mm,
            axis_direction,
            radius_mm,
            ..
        } => {
            let length = dot(axis_direction, axis_direction).sqrt();
            if !length.is_finite() || length <= 0. {
                return None;
            }
            let unit = axis_direction.map(|a| a / length);
            let helper = if unit[0].abs() < 0.9 {
                [1., 0., 0.]
            } else {
                [0., 1., 0.]
            };
            let first = {
                let c = cross(unit, helper);
                let n = dot(c, c).sqrt();
                c.map(|x| x / n)
            };
            let second = cross(unit, first);
            let v: [B; 3] = std::array::from_fn(|i| local[i] - B::point(axis_origin_mm[i]));
            let along = |direction: [f64; 3]| {
                (0..3).fold(B::ZERO, |sum, i| sum + v[i] * B::point(direction[i]))
            };
            let radial = along(first).square() + along(second).square();
            if !radial.valid() {
                return None;
            }
            let magnitude = local
                .iter()
                .flat_map(|b| [b.lo, b.hi])
                .chain(axis_origin_mm)
                .map(f64::abs)
                .fold(radius_mm.max(1.), f64::max);
            let slack = margin + 256. * f64::EPSILON * magnitude;
            let outer = radius_mm + slack;
            let inner = radius_mm - slack;
            if radial.lo > outer * outer {
                Some(radial.lo.sqrt() - outer)
            } else if inner > 0. && radial.hi < inner * inner {
                Some(inner - radial.hi.sqrt())
            } else {
                None
            }
        }
    }
}

fn plane_clearance(normal: Point, point: Point, local: [B; 3], margin: f64) -> Option<f64> {
    let distance = (0..3).fold(B::ZERO, |sum, i| {
        sum + B::point(normal[i]) * (local[i] - B::point(point[i]))
    });
    if !distance.valid() {
        return None;
    }
    let magnitude = local
        .iter()
        .flat_map(|b| [b.lo, b.hi])
        .chain(point)
        .map(f64::abs)
        .fold(1., f64::max);
    let slack = margin + 256. * f64::EPSILON * magnitude;
    if distance.lo > slack {
        Some(distance.lo - slack)
    } else if distance.hi < -slack {
        Some(-distance.hi - slack)
    } else {
        None
    }
}

/// The gap between the box and a planar face's own extent. A planar region
/// lies within the box of its outer loop, so a face whose infinite support
/// crosses the probe can still clear it.
fn face_box_gap(
    shell: &Shell,
    face: &warped_union::Face,
    local: [B; 3],
    margin: f64,
) -> Option<f64> {
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for used in &face.boundary {
        let edge = shell.edges().get(used.edge)?;
        for &vertex in &edge.vertices {
            let p = *shell.vertices().get(vertex)?;
            for i in 0..3 {
                low[i] = low[i].min(p[i]);
                high[i] = high[i].max(p[i]);
            }
        }
        if let Curve::Ellipse(arc) = edge.curve {
            for i in 0..3 {
                let [lo, hi] = arc.range_along(std::array::from_fn(|j| f64::from(i == j)));
                low[i] = low[i].min(lo);
                high[i] = high[i].max(hi);
            }
        }
    }
    let magnitude = low
        .into_iter()
        .chain(high)
        .chain(local.iter().flat_map(|b| [b.lo, b.hi]))
        .map(f64::abs)
        .fold(1., f64::max);
    let slack = margin + 256. * f64::EPSILON * magnitude;
    let gap = (0..3)
        .map(|i| (local[i].lo - high[i]).max(low[i] - local[i].hi))
        .fold(f64::NEG_INFINITY, f64::max);
    (gap.is_finite() && gap > slack).then_some(gap - slack)
}

/// The gap along the axis between the box and a cylindrical face. The axial
/// coordinate is linear on the unrolled patch, so its range over the patch
/// is its range over the patch boundary.
fn cylinder_axial_gap(
    shell: &Shell,
    face: &warped_union::Face,
    local: [B; 3],
    margin: f64,
) -> Option<f64> {
    let Support::Cylinder {
        axis_origin_mm,
        axis_direction,
        ..
    } = face.support
    else {
        return None;
    };
    let length = dot(axis_direction, axis_direction).sqrt();
    let axis = axis_direction.map(|a| a / length);
    let along = (0..3).fold(B::ZERO, |sum, i| {
        sum + (local[i] - B::point(axis_origin_mm[i])) * B::point(axis[i])
    });
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    for used in face.boundary.iter().chain(face.holes.iter().flatten()) {
        let edge = shell.edges().get(used.edge)?;
        for &vertex in &edge.vertices {
            let t = dot(
                subtract(*shell.vertices().get(vertex)?, axis_origin_mm),
                axis,
            );
            low = low.min(t);
            high = high.max(t);
        }
        if let Curve::Ellipse(arc) = edge.curve {
            let offset = dot(axis_origin_mm, axis);
            let [lo, hi] = arc.range_along(axis);
            low = low.min(lo - offset);
            high = high.max(hi - offset);
        }
    }
    if !along.valid() || !low.is_finite() || !high.is_finite() {
        return None;
    }
    let magnitude = [along.lo, along.hi, low, high]
        .into_iter()
        .chain(axis_origin_mm)
        .map(f64::abs)
        .fold(1., f64::max);
    let slack = margin + 256. * f64::EPSILON * magnitude;
    let gap = (along.lo - high).max(low - along.hi);
    (gap > slack).then_some(gap - slack)
}

/// The signed angle `a` to `b` subtends at `center` in the section plane.
fn turn(center: Point, a: Point, b: Point) -> f64 {
    let [ax, ay] = [a[0] - center[0], a[1] - center[1]];
    let [bx, by] = [b[0] - center[0], b[1] - center[1]];
    (ax * by - ay * bx).atan2(ax * bx + ay * by)
}

/// Total angle swept by one oriented section edge, or `None` if the arc
/// needs more steps than the budget allows.
fn edge_turn(edge: &SectionEdge, center: Point, clearance: f64, budget: &mut usize) -> Option<f64> {
    let [start, end] = edge.ends_mm;
    let Curve::Ellipse(arc) = edge.curve else {
        return Some(turn(center, start, end));
    };
    let ellipse = arc.ellipse();
    let [major, _] = ellipse.semiaxes_mm();
    let [first, last] = arc.endpoints();
    let gap = |p: Point, q: Point| norm(subtract(p, q));
    // The edge may run against the arc's parameter direction.
    let reversed = gap(first, start) > gap(last, start);
    // Sagitta of a parameter step h is at most major * h^2 / 8.
    let step = (4. * clearance / major).sqrt();
    let sweep = arc.sweep_rad();
    let steps = (sweep.abs() / step).ceil().max(1.);
    if !steps.is_finite() || steps as usize > *budget {
        return None;
    }
    let steps = steps as usize;
    *budget -= steps;
    let mut previous = start;
    let mut total = 0.;
    for k in 1..=steps {
        let fraction = k as f64 / steps as f64;
        let point = if k == steps {
            end
        } else {
            let fraction = if reversed { 1. - fraction } else { fraction };
            ellipse.point_at(arc.start_rad() + sweep * fraction)
        };
        total += turn(center, previous, point);
        previous = point;
    }
    Some(total)
}

type Planar = [f64; 2];

fn sub2(a: Planar, b: Planar) -> Planar {
    [a[0] - b[0], a[1] - b[1]]
}

fn cross2(a: Planar, b: Planar) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn dot2(a: Planar, b: Planar) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

/// Counter-clockwise convex hull of a few points.
fn hull(mut points: Vec<Planar>) -> Vec<Planar> {
    points.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
    points.dedup();
    if points.len() < 3 {
        return points;
    }
    let mut lower: Vec<Planar> = Vec::new();
    for &p in &points {
        while lower.len() >= 2
            && cross2(
                sub2(lower[lower.len() - 1], lower[lower.len() - 2]),
                sub2(p, lower[lower.len() - 2]),
            ) <= 0.
        {
            lower.pop();
        }
        lower.push(p);
    }
    let mut upper: Vec<Planar> = Vec::new();
    for &p in points.iter().rev() {
        while upper.len() >= 2
            && cross2(
                sub2(upper[upper.len() - 1], upper[upper.len() - 2]),
                sub2(p, upper[upper.len() - 2]),
            ) <= 0.
        {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

fn point_segment(p: Planar, a: Planar, b: Planar) -> f64 {
    let d = sub2(b, a);
    let length = dot2(d, d);
    let t = if length > 0. {
        (dot2(sub2(p, a), d) / length).clamp(0., 1.)
    } else {
        0.
    };
    let q = [a[0] + t * d[0], a[1] + t * d[1]];
    let r = sub2(p, q);
    r[0].hypot(r[1])
}

fn inside_convex(p: Planar, polygon: &[Planar]) -> bool {
    polygon.len() >= 3
        && (0..polygon.len()).all(|i| {
            let [a, b] = [polygon[i], polygon[(i + 1) % polygon.len()]];
            cross2(sub2(b, a), sub2(p, a)) >= 0.
        })
}

/// Distance between a segment and a convex polygon, zero when they meet.
fn segment_polygon(a: Planar, b: Planar, polygon: &[Planar]) -> f64 {
    if inside_convex(a, polygon) || inside_convex(b, polygon) {
        return 0.;
    }
    let mut best = f64::INFINITY;
    for i in 0..polygon.len() {
        let [c, d] = [polygon[i], polygon[(i + 1) % polygon.len()]];
        let crosses = {
            let [d1, d2] = [
                cross2(sub2(b, a), sub2(c, a)),
                cross2(sub2(b, a), sub2(d, a)),
            ];
            let [d3, d4] = [
                cross2(sub2(d, c), sub2(a, c)),
                cross2(sub2(d, c), sub2(b, c)),
            ];
            d1 * d2 < 0. && d3 * d4 < 0.
        };
        if crosses {
            return 0.;
        }
        best = best
            .min(point_segment(a, c, d))
            .min(point_segment(b, c, d))
            .min(point_segment(c, a, b))
            .min(point_segment(d, a, b));
    }
    best
}

/// A loop of oriented shell edges as a polyline within `sagitta` of it.
fn loop_polyline(
    shell: &Shell,
    ring: &[warped_union::OrientedEdge],
    sagitta: f64,
    budget: &mut usize,
) -> Option<Vec<Point>> {
    let mut points = Vec::new();
    for used in ring {
        let edge = shell.edges().get(used.edge)?;
        let [a, b] = edge.vertices.map(|v| shell.vertices().get(v).copied());
        let (start, end) = if used.forward { (a?, b?) } else { (b?, a?) };
        points.push(start);
        if let Curve::Ellipse(arc) = edge.curve {
            let ellipse = arc.ellipse();
            let [major, _] = ellipse.semiaxes_mm();
            let [first, last] = arc.endpoints();
            let reversed = norm(subtract(first, start)) > norm(subtract(last, start));
            let step = (8. * sagitta / major).sqrt();
            let steps = (arc.sweep_rad().abs() / step).ceil().max(1.);
            if !steps.is_finite() || steps as usize > *budget {
                return None;
            }
            let steps = steps as usize;
            *budget -= steps;
            for k in 1..steps {
                let fraction = k as f64 / steps as f64;
                let fraction = if reversed { 1. - fraction } else { fraction };
                points.push(ellipse.point_at(arc.start_rad() + arc.sweep_rad() * fraction));
            }
        }
        let _ = end;
    }
    Some(points)
}

/// The in-plane clearance of a planar face from the box when the box's
/// projection lies strictly inside one hole loop or strictly outside the
/// outer loop. The polyline stays within its sagitta of the loop, so a
/// projection that clears the polyline by more than the sagitta clears the
/// loop and has the same winding number about it.
fn loop_clearance(
    shell: &Shell,
    face: &warped_union::Face,
    normal: Point,
    local: [B; 3],
    margin: f64,
) -> Option<f64> {
    let helper = if normal[0].abs() < 0.9 {
        [1., 0., 0.]
    } else {
        [0., 1., 0.]
    };
    let u = {
        let c = cross(normal, helper);
        let n = dot(c, c).sqrt();
        c.map(|x| x / n)
    };
    let v = cross(normal, u);
    let flat = |p: Point| [dot(p, u), dot(p, v)];
    let mut corners = Vec::with_capacity(8);
    for i in 0..8 {
        let corner: Point = std::array::from_fn(|k| {
            if i >> k & 1 == 0 {
                local[k].lo
            } else {
                local[k].hi
            }
        });
        corners.push(flat(corner));
    }
    let shadow = hull(corners);
    if shadow.len() < 3 {
        return None;
    }
    let magnitude = shadow.iter().flatten().map(|x| x.abs()).fold(1., f64::max);
    let sagitta = 1e-4_f64.max(margin);
    let slack = margin + sagitta + 1024. * f64::EPSILON * magnitude;
    let mut budget = ARC_STEPS;
    let candidates = face
        .holes
        .iter()
        .map(|ring| (ring, true))
        .chain([(&face.boundary, false)]);
    for (ring, hole) in candidates {
        let polyline: Vec<Planar> = loop_polyline(shell, ring, sagitta, &mut budget)?
            .into_iter()
            .map(flat)
            .collect();
        let gap = (0..polyline.len())
            .map(|i| segment_polygon(polyline[i], polyline[(i + 1) % polyline.len()], &shadow))
            .fold(f64::INFINITY, f64::min);
        // A NaN gap is never clear.
        if gap.is_nan() || gap <= slack {
            continue;
        }
        let wound = shadow.iter().all(|&p| {
            let total: f64 = (0..polyline.len())
                .map(|i| {
                    let [a, b] = [
                        sub2(polyline[i], p),
                        sub2(polyline[(i + 1) % polyline.len()], p),
                    ];
                    cross2(a, b).atan2(dot2(a, b))
                })
                .sum();
            let winding = (total / TAU).round();
            if hole {
                winding.abs() == 1.
            } else {
                winding == 0.
            }
        });
        if wound {
            return Some(gap - slack);
        }
    }
    None
}

/// Whole-box relation of one local box to the shell.
fn classify(shell: &Shell, local: [B; 3], margin: f64) -> Option<BoxClassification> {
    if !local.iter().all(|b| b.valid()) {
        return None;
    }
    let mut nearest = f64::INFINITY;
    for face in shell.faces() {
        let gap = match face.support {
            Support::Plane { normal } => plane_clearance(normal, face.point_mm, local, margin)
                .or_else(|| face_box_gap(shell, face, local, margin))
                .or_else(|| loop_clearance(shell, face, normal, local, margin)),
            support => cylinder_clearance(support, local, margin)
                .or_else(|| cylinder_axial_gap(shell, face, local, margin)),
        };
        let gap = gap?;
        nearest = nearest.min(gap);
    }
    if nearest.is_nan() || nearest <= margin {
        return None;
    }
    // Every point of the box has the same relation; try a few section
    // levels in case one meets a horizontal face exactly.
    let center = local.map(|b| b.lo.midpoint(b.hi));
    let rise = 0.25 * (local[2].hi - local[2].lo);
    for offset in [0., rise, -rise] {
        let point = [center[0], center[1], center[2] + offset];
        let Ok(loops) = shell.section(point, [0., 0., 1.], margin.max(1e-7), SECTION_SEGMENTS)
        else {
            continue;
        };
        let mut budget = ARC_STEPS;
        let mut total = 0.;
        for ring in &loops {
            for edge in ring {
                total += edge_turn(edge, point, nearest, &mut budget)?;
            }
        }
        let winding = total / TAU;
        let rounded = winding.round();
        if (winding - rounded).abs() > 0.25 || rounded.abs() > 1. {
            return None;
        }
        return Some(if rounded == 0. {
            BoxClassification::Outside
        } else {
            BoxClassification::Inside
        });
    }
    None
}

impl Solid {
    pub(in crate::geometry) fn warped_shell_probe_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        let (probe, source_is_left, (body, frame)) = match self.warped_union_body() {
            Some(found) => (other, true, found),
            None => (self, false, other.warped_union_body()?),
        };
        let margin = policy.linear_tolerance_mm;
        // A rectilinear probe is the union of its cells; any other probe is
        // enclosed by its outer box.
        let pieces: Vec<[B; 3]> = match probe.orthogonal_cells() {
            Ok(Some(cells)) => cells
                .iter()
                .map(|cell| std::array::from_fn(|i| B::new(cell.min[i], cell.max[i])))
                .collect(),
            _ => match probe_bounds(probe) {
                Some(bounds) => vec![bounds],
                // Any other Rust solid is enclosed by its bounds, widened by
                // the construction tolerance.
                None if probe.uses_rust_geometry() => {
                    let bounds = probe.bounds().ok()?;
                    vec![std::array::from_fn(|i| {
                        B::new(bounds.min[i] - margin, bounds.max[i] + margin)
                    })]
                }
                None => return None,
            },
        };
        let mut relation = None;
        for piece in pieces {
            let local = chart_bounds(piece, frame.origin, [frame.u, frame.v, frame.normal()])
                .map(|[lo, hi]| B::new(lo, hi));
            let next = classify(&body.geometry, local, margin)?;
            if relation.is_some_and(|known| known != next) {
                return None;
            }
            relation = Some(next);
        }
        whole_operand_result(self, other, source_is_left, relation?, operation)
    }
}

#[cfg(test)]
mod tests;
