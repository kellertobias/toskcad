//! Actual retained original BRep replay, separately proved output operator.
use super::*;
use multicad_rust_kernel::step::{self, literal_npt as replay};
fn failure(cause: read::Error) -> step::StepError {
    match cause {
        read::Error::Invalid => step::StepError::InvalidGeometry,
        read::Error::Limit => step::StepError::Limit,
        read::Error::Unsupported | read::Error::Cancelled => step::StepError::Unsupported,
    }
}
impl LiteralNptBody {
    fn step_options(&self) -> Result<replay::Options, step::StepError> {
        let extra = self.plan.retained_budget_bytes() - self.source.retained_budget_bytes();
        let maximum_bytes = (BYTES - OVERHEAD)
            .checked_sub(extra)
            .and_then(|v| v.checked_sub(self.cached_bytes()))
            .ok_or(step::StepError::Limit)?;
        Ok(replay::Options {
            maximum_bytes,
            position_error_mm: self.policy.linear_tolerance_mm,
            ..Default::default()
        })
    }
    fn candidate(
        &self,
        options: replay::Options,
    ) -> Result<replay::VerifiedReplay<'_>, step::StepError> {
        match &self.occurrence {
            Some(p) => replay::verify_placed_replay(&self.source, p, options, None),
            None => replay::verify_exact_replay(&self.source, options, None),
        }
        .map_err(failure)
    }
    pub(in crate::geometry) fn supports_step(&self) -> bool {
        let Ok(options) = self.step_options() else {
            return false;
        };
        if let Some((bytes, ready)) = self
            .step_ready
            .get()
            .filter(|(bytes, _)| *bytes == options.maximum_bytes)
        {
            let _ = bytes;
            return ready;
        }
        let ready = self.candidate(options).is_ok();
        self.step_ready.set(Some((options.maximum_bytes, ready)));
        ready
    }
    pub(in crate::geometry) fn append_step(
        &self,
        writer: &mut step::Writer,
    ) -> Result<(), step::StepError> {
        let options = self.step_options()?;
        let candidate = self.candidate(options)?;
        let receipt = writer
            .append_verified_literal_npt(&candidate, options, None)
            .map_err(failure)?;
        self.last_step.set(Some(receipt));
        Ok(())
    }
}
impl Solid {
    /// Fresh actual-file verification of exactly one represented member.
    pub fn verify_literal_npt_step(
        &self,
        actual_text: &str,
        tolerance: f64,
    ) -> Result<replay::ReadbackReceipt, GeometryError> {
        let result = self.verify_literal_npt_step_member(actual_text, 0, tolerance)?;
        if result.represented_member_count() != 1 {
            return Err(super::failure(read::Error::Unsupported));
        }
        Ok(result)
    }
    /// Fresh proof for one explicitly ordered actual member. This does not
    /// validate the other compound members or their expected cardinality.
    pub fn verify_literal_npt_step_member(
        &self,
        actual_text: &str,
        member_index: usize,
        tolerance: f64,
    ) -> Result<replay::ReadbackReceipt, GeometryError> {
        self.verify_literal_npt_step_member_with_resident_bytes(
            actual_text,
            member_index,
            tolerance,
            0,
        )
    }
    /// As above, with other operation-owned retained state explicitly charged.
    /// Ordinary per-member receipts cover this subject and the actual whole
    /// input document, rather than all scene or application allocations.
    pub fn verify_literal_npt_step_member_with_resident_bytes(
        &self,
        actual_text: &str,
        member_index: usize,
        tolerance: f64,
        other_resident_bytes: usize,
    ) -> Result<replay::ReadbackReceipt, GeometryError> {
        let body = self
            .literal_npt_body()
            .ok_or_else(|| super::failure(read::Error::Unsupported))?;
        let already_charged = body.source.retained_budget_bytes()
            + body
                .occurrence
                .as_ref()
                .map_or(0, |p| p.retained_budget_bytes());
        let additional_resident_bytes = body
            .resident_bytes()
            .checked_sub(already_charged)
            .and_then(|n| n.checked_add(OVERHEAD))
            .and_then(|n| n.checked_add(other_resident_bytes))
            .ok_or_else(|| super::failure(read::Error::Limit))?;
        replay::verify_readback(
            &body.source,
            body.occurrence.as_ref(),
            actual_text,
            member_index,
            replay::ReadbackOptions {
                additional_resident_bytes,
                position_error_mm: tolerance,
                ..Default::default()
            },
            None,
        )
        .map_err(super::failure)
    }
    /// Last successful verified in-memory append. Filesystem publication is a
    /// separate final IO action and is not asserted by this geometry receipt.
    pub fn rust_literal_npt_step_evidence(&self) -> Option<replay::Receipt> {
        self.literal_npt_body()?.last_step.get()
    }
}
