//! Qualified local hanger export with explicit support-healing provenance.
use super::*;
use multicad_rust_kernel::{
    step::{self, hanger as writer},
    step_read::{literal_member_binding as binding, owned_insert as read},
};
#[cfg(test)]
mod tests;
#[derive(Clone, Debug, PartialEq)]
pub struct HealedHangerStepVerification {
    pub healed_volume_bounds_mm3: [f64; 2],
    pub source_volume_bounds_mm3: [f64; 2],
    /// Source envelope expanded by measured correspondence, not raw file bounds.
    pub correspondence_enclosure_mm: [[f64; 3]; 2],
    pub source_to_raw_position_error_mm: f64,
    pub sewing_correspondence_bound_mm: f64,
    pub source_to_healed_position_error_mm: f64,
    pub source_support_derivative_error_mm: [f64; 2],
    /// Original source to transmitted coedges; this is not healed-cap D1.
    pub original_coedge_derivative_error_mm: f64,
    pub verified_tolerance_mm: f64,
    pub matched_entities: usize,
    pub literal_binding_work: usize,
}
fn failure(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error("qualified hanger STEP", kind, message)
}
fn writer_error(e: writer::Error) -> GeometryError {
    failure(
        match e {
            writer::Error::Invalid => GeometryErrorKind::Input,
            writer::Error::Limit => GeometryErrorKind::Limit,
            writer::Error::Cancelled => GeometryErrorKind::Cancelled,
            writer::Error::Unresolved => GeometryErrorKind::Unsupported,
        },
        format!("Original hanger export qualification failed: {e:?}"),
    )
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(failure(
            GeometryErrorKind::Cancelled,
            "Hanger export superseded",
        ))
    } else {
        Ok(())
    }
}
impl HangerBody {
    fn qualified_step(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Rc<writer::QualifiedHangerStep>, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        if frame != PlaneFrame::default() {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Placed hanger STEP correspondence is not yet qualified",
            ));
        }
        let cached = self.qualified_step.borrow().as_ref().cloned();
        let qualified = match cached {
            Some(file) => file,
            None => {
                let file = Rc::new(
                    writer::export_qualified_hanger(&self.source, Default::default(), interrupt)
                        .map_err(writer_error)?,
                );
                *self.qualified_step.borrow_mut() = Some(Rc::clone(&file));
                file
            }
        };
        if qualified.source_to_healed_position_error_mm() > policy.linear_tolerance_mm {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Measured hanger source/healing correspondence exceeds the requested tolerance",
            ));
        }
        checkpoint(interrupt)?;
        Ok(qualified)
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
            .append_qualified_hanger(&file, Default::default(), None)
            .map(|_| ())
            .map_err(|e| convert(writer_error(e)))
    }
    pub(in crate::geometry) fn verify_step_member(
        &self,
        actual: &read::CompoundDocument<'_>,
        member: read::Member,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HealedHangerStepVerification, GeometryError> {
        if member.mapping_link().is_some() {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Mapped hanger STEP needs a separate placement proof",
            ));
        }
        let qualified = self.qualified_step(frame, policy, interrupt)?;
        // The independently decoded candidate already earned finite lifts,
        // healing and volume. Exact local graph binding transfers those proofs
        // to the actual written member, preserving every literal leaf/owner.
        let expected =
            read::CompoundDocument::parse(qualified.as_str(), Default::default(), interrupt)
                .map_err(|e| {
                    failure(
                        match e {
                            read::Error::Limit => GeometryErrorKind::Limit,
                            read::Error::Cancelled => GeometryErrorKind::Cancelled,
                            _ => GeometryErrorKind::InvalidResult,
                        },
                        format!("Qualified candidate representation failed: {e:?}"),
                    )
                })?;
        let matched = binding::bind_local_members(
            &expected,
            expected.represented_members()[0].item_entity(),
            actual,
            member.item_entity(),
            binding::Options {
                // Large literal tensors contain many exact decimal digits.
                // This is the explicit local-record binding phase; original
                // encoder/decoder/material query frontiers remain unchanged.
                maximum_work: 256_000_000,
                ..Default::default()
            },
            interrupt,
        )
        .map_err(|e| {
            failure(
                match e {
                    binding::Error::Limit
                    | binding::Error::WorkLimit
                    | binding::Error::StorageLimit
                    | binding::Error::EntityLimit => GeometryErrorKind::Limit,
                    binding::Error::Cancelled => GeometryErrorKind::Cancelled,
                    binding::Error::Invalid | binding::Error::Different => {
                        GeometryErrorKind::InvalidResult
                    }
                },
                format!("Actual hanger member differs from its qualified literal geometry: {e:?}"),
            )
        })?;
        let extent = self
            .source
            .bounds_mm()
            .into_iter()
            .flatten()
            .map(f64::abs)
            .fold(1., f64::max);
        let source_volume = self.material(frame, policy)?.volume_bounds(
            QueryOptions {
                volume_error_mm3: (policy.linear_tolerance_mm * extent * extent).max(1e-6),
                maximum_cells: self.source.wood().bands().len() + 6,
            },
            interrupt,
        )?;
        let bounds = self.source.bounds_mm();
        let residual = qualified.source_to_healed_position_error_mm();
        checkpoint(interrupt)?;
        Ok(HealedHangerStepVerification {
            healed_volume_bounds_mm3: qualified.healed_volume_bounds_mm3(),
            source_volume_bounds_mm3: source_volume,
            correspondence_enclosure_mm: [
                bounds.map(|v| (v[0] - residual).next_down()),
                bounds.map(|v| (v[1] + residual).next_up()),
            ],
            source_to_raw_position_error_mm: qualified.source_to_raw_position_error_mm(),
            sewing_correspondence_bound_mm: qualified.sewing_correspondence_bound_mm(),
            source_to_healed_position_error_mm: residual,
            source_support_derivative_error_mm: qualified.source_support_derivative_error_mm(),
            original_coedge_derivative_error_mm: qualified.original_coedge_derivative_error_mm(),
            verified_tolerance_mm: policy.linear_tolerance_mm,
            matched_entities: matched.matched_entities(),
            literal_binding_work: matched.work_used(),
        })
    }
}
