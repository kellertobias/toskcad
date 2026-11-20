//! Conservative whole-prism relations against retained clipped NPT material.
//!
//! The world probe box is mapped through the true affine inverse of the
//! occurrence pose (with its rounding enclosure) and classified by the
//! original source's certified predicate. A proved relation shares one
//! unchanged operand; undetermined probes are left to other providers.
use super::*;
use crate::geometry::helical_body::contained_boolean::{probe_bounds, whole_operand_result};
use multicad_rust_kernel::helical_thread::BoxClassification;

impl Solid {
    pub(in crate::geometry) fn contained_clipped_npt_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        let (probe, source_is_left, body) = match self.clipped_npt_body() {
            Some(found) => (other, true, found),
            None => (self, false, other.clipped_npt_body()?),
        };
        let world = probe_bounds(probe)?.map(|b| [b.lo, b.hi]);
        let local = if body.pose == RigidTransform::default() {
            world
        } else {
            body.calibration().inverse_point_bounds(world).ok()?
        };
        let relation = match body
            .source
            .classify_box(local, Default::default(), None)
            .ok()?
            .classification()
        {
            StockClipClassification::Inside => BoxClassification::Inside,
            StockClipClassification::Outside => BoxClassification::Outside,
            StockClipClassification::Undetermined => return None,
        };
        whole_operand_result(self, other, source_is_left, relation, operation)
    }
}
