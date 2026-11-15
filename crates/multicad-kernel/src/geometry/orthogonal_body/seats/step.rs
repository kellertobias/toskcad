//! Cached exact-source qualification and fresh actual-file seated verification.
use super::*;
use multicad_rust_kernel::{
    step::{self, orthogonal_seats as writer},
    step_read::owned_insert as read,
};
#[cfg(test)]
mod tests;
fn failure(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error("qualified perpendicular-seat STEP", kind, message)
}
fn proof_error(e: writer::QualificationError) -> GeometryError {
    use writer::QualificationError as E;
    failure(
        match e {
            E::Invalid => GeometryErrorKind::Input,
            E::Unsupported => GeometryErrorKind::Unsupported,
            E::Limit => GeometryErrorKind::Limit,
            E::Cancelled => GeometryErrorKind::Cancelled,
            E::Encoding(step::StepError::Limit) => GeometryErrorKind::Limit,
            E::Encoding(step::StepError::Unsupported) => GeometryErrorKind::Unsupported,
            E::Encoding(_) => GeometryErrorKind::InvalidResult,
        },
        format!("Finite source/file qualification failed: {e:?}"),
    )
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(failure(
            GeometryErrorKind::Cancelled,
            "Seated STEP verification superseded",
        ))
    } else {
        Ok(())
    }
}
fn options(policy: GeometryPolicy) -> writer::QualificationOptions {
    let mut options = writer::QualificationOptions::default();
    options.material.clearance_mm = policy.linear_tolerance_mm;
    options.material.minimum_feature_mm = policy.minimum_feature_mm;
    options.material.transfer.maximum_position_error_mm = policy.linear_tolerance_mm;
    options
        .material
        .transfer
        .maximum_normalized_derivative_error_mm = policy.linear_tolerance_mm;
    options.material.transfer.maximum_sewing_error_mm = policy.linear_tolerance_mm;
    options
}
impl OrthogonalSeatBody {
    pub(super) fn qualified_step(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Rc<writer::QualifiedOrthogonalSeatsStep>, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        if policy.boolean_fuzzy_mm != 0. {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Fuzzy source rewriting is not a seated export proof",
            ));
        }
        let placement = transport::placement(frame)?;
        let cached = self.qualified_step.borrow().as_ref().cloned();
        if let Some((quality, pose, file)) = cached {
            if quality == policy && pose == frame && file.matches_source(&self.source, placement) {
                checkpoint(interrupt)?;
                return Ok(file);
            }
        }
        self.validate_frame(frame, policy)?;
        let file = Rc::new(
            writer::qualify(&self.source, placement, options(policy), &|| {
                interrupt.is_some_and(|f| f())
            })
            .map_err(proof_error)?,
        );
        if !file.matches_source(&self.source, placement)
            || file.source_to_healed_position_error_mm() > policy.linear_tolerance_mm
            || file.source_to_healed_normalized_derivative_error_mm() > policy.linear_tolerance_mm
            || file
                .raw_to_healed_correspondence()
                .maximum_sewing_error_mm()
                > policy.linear_tolerance_mm
        {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Earned seated transfer does not match the current source, pose and policy",
            ));
        }
        checkpoint(interrupt)?;
        *self.qualified_step.borrow_mut() = Some((policy, frame, Rc::clone(&file)));
        Ok(file)
    }
    pub(in crate::geometry) fn supports_step(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> bool {
        self.qualified_step(frame, policy, None).is_ok()
    }
    pub(in crate::geometry) fn append_step(
        &self,
        target: &mut step::Writer,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<(), step::StepError> {
        let convert = |e: GeometryError| match e.kind {
            GeometryErrorKind::Limit => step::StepError::Limit,
            GeometryErrorKind::Unsupported | GeometryErrorKind::Cancelled => {
                step::StepError::Unsupported
            }
            _ => step::StepError::InvalidGeometry,
        };
        let file = self.qualified_step(frame, policy, None).map_err(convert)?;
        target
            .append_qualified_orthogonal_seats(&file, Default::default(), &|| false)
            .map(|_| ())
            .map_err(|e| convert(proof_error(e)))
    }
    pub(in crate::geometry) fn verify_step_member(
        &self,
        actual: &read::CompoundDocument<'_>,
        member: read::Member,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HealedSupportStepVerification, GeometryError> {
        checkpoint(interrupt)?;
        let qualified = self.qualified_step(frame, policy, interrupt)?;
        let proof = qualified
            .verify_actual_member(
                actual,
                member.item_entity(),
                writer::VerificationOptions {
                    material: options(policy).material,
                    ..Default::default()
                },
                &|| interrupt.is_some_and(|f| f()),
            )
            .map_err(|cause| {
                let mut error = proof_error(cause);
                if !matches!(
                    error.kind,
                    GeometryErrorKind::Limit | GeometryErrorKind::Cancelled
                ) {
                    error.kind = GeometryErrorKind::InvalidResult;
                }
                error
            })?;
        let source_volume = self.world_volume_bounds(frame)?;
        let verified = proof.original_source_volume_bounds_mm3();
        if source_volume[0] > verified[1] || verified[0] > source_volume[1] {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Current affine source and sealed source mass enclosures disagree",
            ));
        }
        let bounds = self.world_bounds(frame)?;
        let residual = proof.source_to_healed_position_error_mm();
        let source = proof.source_to_raw_correspondence();
        let healed = proof.raw_to_healed_correspondence();
        if residual > policy.linear_tolerance_mm
            || proof.source_to_healed_normalized_derivative_error_mm() > policy.linear_tolerance_mm
            || healed.maximum_sewing_error_mm() > policy.linear_tolerance_mm
        {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Actual source/file correspondence exceeds the current policy",
            ));
        }
        checkpoint(interrupt)?;
        Ok(HealedSupportStepVerification {
            healed_volume_bounds_mm3: proof.healed_volume_bounds_mm3(),
            source_volume_bounds_mm3: source_volume,
            correspondence_enclosure_mm: [
                bounds.min.map(|v| (v - residual).next_down()),
                bounds.max.map(|v| (v + residual).next_up()),
            ],
            source_to_raw_position_error_mm: source.maximum_position_error_mm(),
            sewing_correspondence_bound_mm: healed.maximum_sewing_error_mm(),
            source_to_healed_position_error_mm: residual,
            source_support_derivative_error_mm: [source.maximum_normalized_derivative_error_mm();
                2],
            original_coedge_derivative_error_mm: source.maximum_normalized_derivative_error_mm(),
            verified_tolerance_mm: policy.linear_tolerance_mm,
            matched_entities: proof.matched_entities(),
            literal_binding_work: proof.literal_binding_work_used(),
        })
    }
}
