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
    /// Last successful verified in-memory append. Filesystem publication is a
    /// separate final IO action and is not asserted by this geometry receipt.
    pub fn rust_literal_npt_step_evidence(&self) -> Option<replay::Receipt> {
        self.literal_npt_body()?.last_step.get()
    }
}
