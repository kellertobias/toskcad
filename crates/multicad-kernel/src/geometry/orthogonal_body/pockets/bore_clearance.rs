//! A cylinder inside a certified through bore meets no stock material.
//!
//! A through bore removes the whole cylindrical column between its two
//! original stock faces, so a parallel cylinder whose disk lies in the bore
//! disk and whose axial interval lies between those faces is entirely void.
//! Its intersection with the stock is empty and cutting it removes nothing.
//! The interval may overshoot a face only by coordinate rounding, which is
//! far below the stock's minimum feature, so no material lies there either.
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
        let inside = body.through_bores.iter().any(|bore| {
            let [a, b] = bore.centers_mm;
            let bore_axis =
                (0..3).max_by(|&i, &j| (b[i] - a[i]).abs().total_cmp(&(b[j] - a[j]).abs()));
            if bore_axis != Some(axis) || (0..3).any(|i| i != axis && a[i] != b[i]) {
                return false;
            }
            let [i, j] = [(axis + 1) % 3, (axis + 2) % 3];
            let offset = (start[i] - a[i]).hypot(start[j] - a[j]);
            let [low, high] = [a[axis].min(b[axis]), a[axis].max(b[axis])];
            offset + radius + numeric <= bore.radius_mm
                && start[axis].min(end[axis]) >= low - numeric
                && start[axis].max(end[axis]) <= high + numeric
        });
        inside.then(|| match operation {
            BooleanOperation::Intersect => Ok(None),
            _ => self.share().map(Some),
        })
    }
}
