//! Whole-probe material relations against the retained paired insert source.
use super::*;
use crate::geometry::helical_body::contained_boolean::{
    chart_bounds, probe_bounds, whole_operand_result,
};
use multicad_rust_kernel::{helical_thread::BoxClassification, wood_insert::Classification};

impl Solid {
    pub(in crate::geometry) fn contained_insert_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        let (probe, source_is_left, body, frame) =
            if let Some((body, frame)) = self.wood_insert_body() {
                (other, true, body, frame)
            } else {
                let (body, frame) = other.wood_insert_body()?;
                (self, false, body, frame)
            };
        let bounds = chart_bounds(
            probe_bounds(probe)?,
            frame.origin,
            [frame.u, frame.v, frame.normal()],
        );
        let relation = match body.source.classify_box(bounds, None).ok()? {
            Classification::Inside => BoxClassification::Inside,
            Classification::Outside => BoxClassification::Outside,
            Classification::UnresolvedBoundary => BoxClassification::UnresolvedBoundary,
        };
        whole_operand_result(self, other, source_is_left, relation, operation)
    }
}

#[cfg(test)]
mod tests;
