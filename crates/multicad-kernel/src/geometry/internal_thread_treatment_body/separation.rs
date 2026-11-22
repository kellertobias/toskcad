//! Source-bound whole-pair separation, never a generic Boolean approximation.
use super::*;
use multicad_rust_kernel::helical_thread::ProfileContainmentOptions;
impl Solid {
    /// Earn separation of THIS actual external treated thread and ONE actual
    /// stock-minus-authored-MAX-cavity member in their exact common chart.
    /// Receiving-cutter Chamfer bodies and all added/composite geometry refuse.
    /// Returned gaps remain source-local radial quantities, not world distances.
    pub fn rust_chamfer_receiver_separation<'a>(
        &'a self,
        receiver: &'a Self,
        options: ProfileContainmentOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<source::ChamferSeparation<'a>, GeometryError> {
        checkpoint(interrupt)?;
        self.policy.validate()?;
        if self.policy != receiver.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(unsupported(
                "Whole-pair separation requires the same actual policy without Boolean fuzz",
            ));
        }
        let (male, male_frame) = self.entry_chamfer_body().ok_or_else(|| {
            unsupported("Whole-pair separation requires actual external Chamfer material")
        })?;
        let external = male.external_material().ok_or_else(|| {
            unsupported("A receiving-cutter Chamfer is not external occupied material")
        })?;
        let (receiving, frame) = receiver.internal_thread_treatment_body().ok_or_else(|| {
            unsupported("Whole-pair separation requires actual stock minus an authored MAX cavity")
        })?;
        // These body variants retain their COMPLETE actual affine placement in
        // this frame; unlike imported/mapped bodies there is no hidden operator.
        // Constructor/transform metric guards earned invertibility. The SAME
        // literal affine map preserves the proved local disjointness.
        if male_frame != frame {
            return Err(unsupported(
                "Whole-pair separation requires the exact common retained source frame",
            ));
        }
        receiving
            .material
            .prove_coaxial_chamfer_separation(external, options, interrupt)
            .map_err(failure)
    }
    pub(in crate::geometry) fn separated_treatment_intersection(
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
        let (male, female) = match (&self.body, &other.body) {
            (SolidBody::EntryChamfer(..), SolidBody::InternalThreadTreatment(..)) => (self, other),
            (SolidBody::InternalThreadTreatment(..), SolidBody::EntryChamfer(..)) => (other, self),
            _ => return None,
        };
        // An unearned result supplies no Boolean shortcut. The existing
        // unsupported path remains responsible for every other case.
        male.rust_chamfer_receiver_separation(female, Default::default(), None)
            .ok()?;
        Some(Ok(None))
    }
}

#[cfg(test)]
mod tests;
