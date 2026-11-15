//! Native STEP admission requires independent literal material and measured
//! correspondence to this unchanged source. Placement is a separate gate.
use super::*;
use multicad_rust_kernel::step::{self, npt_stock as npt};

fn step_error(cause: npt::Error) -> step::StepError {
    match cause {
        npt::Error::Limit => step::StepError::Limit,
        npt::Error::Invalid => step::StepError::InvalidGeometry,
        npt::Error::Cancelled | npt::Error::Unresolved => step::StepError::Unsupported,
    }
}

impl ClippedNptBody {
    fn with_step_plan<T>(
        &self,
        policy: GeometryPolicy,
        f: impl FnOnce(&npt::SurfacePlan<'_>) -> Result<T, npt::Error>,
    ) -> Result<T, step::StepError> {
        // Native's literal near-rigid map and STEP's normalized operator are
        // distinct mathematical maps. Do not silently identify them here.
        if self.pose != RigidTransform::default() {
            return Err(step::StepError::Unsupported);
        }
        policy
            .validate()
            .map_err(|_| step::StepError::InvalidGeometry)?;
        let options = npt::CurveOptions {
            position_error_mm: (policy.linear_tolerance_mm / 4.).min(1e-7),
            ..Default::default()
        };
        let curves = npt::CurvePlan::new(&self.source, options, None).map_err(step_error)?;
        let plan = npt::SurfacePlan::new(curves, options.into(), None).map_err(step_error)?;
        f(&plan).map_err(step_error)
    }

    pub(in crate::geometry) fn supports_step(&self, policy: GeometryPolicy) -> bool {
        if self.pose != RigidTransform::default() {
            return false;
        }
        let key = policy.linear_tolerance_mm.to_bits();
        if let Some((previous, supported)) = self.step_capability.get()
            && previous == key
        {
            return supported;
        }
        let supported = self
            .with_step_plan(policy, |plan| {
                plan.verify_candidate(
                    npt::VerifiedOptions {
                        position_error_mm: policy.linear_tolerance_mm,
                        ..Default::default()
                    },
                    None,
                )
                .map(|_| ())
            })
            .is_ok();
        self.step_capability.set(Some((key, supported)));
        supported
    }

    pub(in crate::geometry) fn append_step(
        &self,
        writer: &mut step::Writer,
        policy: GeometryPolicy,
    ) -> Result<(), step::StepError> {
        self.with_step_plan(policy, |plan| {
            let options = npt::VerifiedOptions {
                position_error_mm: policy.linear_tolerance_mm,
                ..Default::default()
            };
            let candidate = plan.verify_candidate(options, None)?;
            writer
                .append_verified_source_npt_candidate(&candidate, options, None)
                .map(|_| ())
        })
    }
}
