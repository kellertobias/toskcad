//! A cylinder in certified holes and clear of the stock meets no material.
//!
//! A through bore removes the whole cylindrical column between its two
//! original stock faces, and a blind pocket the column from its entry face
//! to its floor. Along a parallel cylinder whose disk lies in a hole's disk,
//! that hole's axial interval is void. Every remaining interval of the
//! cylinder must lie in a box that is separated from, or only touches,
//! every stock cell; cuts only remove stock, so the cells enclose the
//! material. Then the cylinder shares no interior with the stock: the
//! intersection is empty and cutting it removes nothing. Ends may meet
//! within coordinate rounding, which is far below the stock's minimum
//! feature.
//! A rectangular prism with circular through holes is read as the identical
//! orthogonal stock with those bores.
use super::*;

impl Solid {
    pub(in crate::geometry) fn through_bore_clearance(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation == BooleanOperation::Union
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
        {
            return None;
        }
        let bored = |solid: &Self| {
            solid.orthogonal_body().is_some() || solid.drawing_bored_equivalent().is_some()
        };
        let (stock, tool) = if bored(self) && other.axial_cylinder_frame().is_some() {
            (self, other)
        } else if operation == BooleanOperation::Intersect
            && bored(other)
            && self.axial_cylinder_frame().is_some()
        {
            (other, self)
        } else {
            return None;
        };
        // A rectangular prism with circular through holes is the identical
        // point set as orthogonal stock with those through bores.
        let equivalent;
        let stock = if stock.orthogonal_body().is_some() {
            stock
        } else {
            equivalent = stock.drawing_bored_equivalent()?.ok()?;
            &equivalent
        };
        let (body, frame) = stock.orthogonal_body()?;
        let (cylinder, tool_frame) = tool.axial_cylinder_frame()?;
        let axes = [frame.u, frame.v, frame.normal()];
        let local = |point| axes.map(|axis| dot(subtract(point, frame.origin), axis));
        let start = local(tool_frame.origin);
        let end = local(transform_point(
            tool_frame,
            [0., 0., cylinder.height().as_f64()],
        ));
        let radius = cylinder.circle().radius().as_f64();
        let bounds = body.geometry.bounds();
        let numeric = 64.
            * f64::EPSILON
            * bounds
                .min
                .into_iter()
                .chain(bounds.max)
                .chain(start)
                .chain(end)
                .map(f64::abs)
                .fold(1., f64::max);
        let delta = subtract(end, start);
        let axis = (0..3).max_by(|&a, &b| delta[a].abs().total_cmp(&delta[b].abs()))?;
        if (0..3).any(|i| i != axis && delta[i].abs() > numeric) {
            return None;
        }
        let [i, j] = [(axis + 1) % 3, (axis + 2) % 3];
        let inside = |center: Point, hole_radius: f64| {
            (start[i] - center[i]).hypot(start[j] - center[j]) + radius + numeric <= hole_radius
        };
        // Axial intervals of the tool's own column that a hole leaves void.
        let mut voids: Vec<[f64; 2]> = Vec::new();
        for bore in body.through_bores.iter() {
            let [a, b] = bore.centers_mm;
            if (0..3).all(|k| k == axis || a[k] == b[k])
                && a[axis] != b[axis]
                && inside(a, bore.radius_mm)
            {
                voids.push([a[axis].min(b[axis]), a[axis].max(b[axis])]);
            }
        }
        for pocket in body.pockets.iter() {
            let normal = body.boundary.faces.get(pocket.face)?.normal;
            if normal[axis].abs() == 1. && inside(pocket.center_mm, pocket.radius_mm) {
                let floor = pocket.floor_center_mm(normal)[axis];
                let entry = pocket.center_mm[axis];
                voids.push([entry.min(floor), entry.max(floor)]);
            }
        }
        // What remains of the tool's interval outside every void must clear
        // every stock cell, contact allowed.
        voids.sort_by(|a, b| a[0].total_cmp(&b[0]));
        let mut remaining = Vec::new();
        let mut low = start[axis].min(end[axis]);
        let high = start[axis].max(end[axis]);
        for [a, b] in voids {
            if b < low - numeric || a > high + numeric {
                continue;
            }
            if a > low + numeric {
                remaining.push([low, a.min(high)]);
            }
            low = low.max(b);
        }
        if low < high - numeric {
            remaining.push([low, high]);
        }
        let inside = remaining.iter().all(|&[a, b]| {
            let mut min = [0.; 3];
            let mut max = [0.; 3];
            min[axis] = a;
            max[axis] = b;
            for k in [i, j] {
                min[k] = start[k] - radius;
                max[k] = start[k] + radius;
            }
            body.geometry.cells().iter().all(|cell| {
                (0..3).any(|k| cell.min[k] >= max[k] - numeric || cell.max[k] <= min[k] + numeric)
            })
        });
        inside.then(|| match operation {
            BooleanOperation::Intersect => Ok(None),
            _ => self.share().map(Some),
        })
    }
}
