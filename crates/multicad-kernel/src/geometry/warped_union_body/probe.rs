//! Whole-probe relations against any retained warped shell, chart or not.
//!
//! A box probe whose enclosure clears every face support (plane or full
//! cylinder) by more than the construction tolerance cannot meet the
//! boundary, so all of it lies on one side. Its centre is then classified by
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
        let local = chart_bounds(
            probe_bounds(probe)?,
            frame.origin,
            [frame.u, frame.v, frame.normal()],
        )
        .map(|[lo, hi]| B::new(lo, hi));
        if !local.iter().all(|b| b.valid()) {
            return None;
        }
        let shell = &body.geometry;
        let mut nearest = f64::INFINITY;
        for face in shell.faces() {
            let gap = match face.support {
                Support::Plane { normal } => plane_clearance(normal, face.point_mm, local, margin),
                support => cylinder_clearance(support, local, margin),
            }?;
            nearest = nearest.min(gap);
        }
        if !(nearest > margin) {
            return None;
        }
        let center = local.map(|b| b.lo.midpoint(b.hi));
        let loops = shell
            .section(center, [0., 0., 1.], margin.max(1e-7), SECTION_SEGMENTS)
            .ok()?;
        let mut budget = ARC_STEPS;
        let mut total = 0.;
        for ring in &loops {
            for edge in ring {
                total += edge_turn(edge, center, nearest, &mut budget)?;
            }
        }
        let winding = total / TAU;
        let rounded = winding.round();
        if (winding - rounded).abs() > 0.25 || rounded.abs() > 1. {
            return None;
        }
        let relation = if rounded == 0. {
            BoxClassification::Outside
        } else {
            BoxClassification::Inside
        };
        whole_operand_result(self, other, source_is_left, relation, operation)
    }
}

#[cfg(test)]
mod tests;
