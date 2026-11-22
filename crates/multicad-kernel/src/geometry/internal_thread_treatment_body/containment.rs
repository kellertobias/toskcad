//! Reuse an unchanged actual box only after whole-source membership is earned.
use super::*;
use crate::geometry::helical_body::contained_boolean::{
    chart_bounds, probe_bounds, whole_operand_result,
};
use interval::Interval as B;
impl Solid {
    pub(in crate::geometry) fn contained_treatment_box_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Intersect
            || policy != self.policy
            || policy != other.policy
            || policy.boolean_fuzzy_mm != 0.
        {
            return None;
        }
        let (source, probe, left) = if matches!(
            self.body,
            SolidBody::EntryChamfer(..) | SolidBody::InternalThreadTreatment(..)
        ) {
            (self, other, true)
        } else if matches!(
            other.body,
            SolidBody::EntryChamfer(..) | SolidBody::InternalThreadTreatment(..)
        ) {
            (other, self, false)
        } else {
            return None;
        };
        // Only actual finite rectangular-prism source geometry supplies the
        // probe; its full outward source hull is transformed, never corners
        // sampled from a display mesh or an assumed point-volume surrogate.
        let (box_source, _) = probe.rust_prism_frame()?;
        if box_source.rounded_end().is_some() {
            return None;
        }
        let world = probe_bounds(probe)?;
        let frame = match &source.body {
            SolidBody::EntryChamfer(body, frame) if body.external_material().is_some() => *frame,
            SolidBody::InternalThreadTreatment(_, frame) => *frame,
            _ => return None,
        };
        // Existing interval cofactors invert the COMPLETE literal source map.
        // Any invalid/inexact inverse hull remains unresolved in the predicate.
        let bounds = chart_bounds(world, frame.origin, [frame.u, frame.v, frame.normal()]);
        if bounds.iter().any(|[lo, hi]| !B::new(*lo, *hi).valid()) {
            return None;
        }
        let relation = match &source.body {
            SolidBody::EntryChamfer(body, _) => body
                .external_material()?
                .classify_treatment_box(bounds)
                .ok()?,
            SolidBody::InternalThreadTreatment(body, _) => {
                body.material.classify_treatment_box(bounds).ok()?
            }
            _ => return None,
        };
        whole_operand_result(self, other, left, relation, operation)
    }
}

#[cfg(test)]
mod tests;
