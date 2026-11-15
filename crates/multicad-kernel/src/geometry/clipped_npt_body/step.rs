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
    /// The shared document has its own bounded global arena. This member query
    /// retains its unchanged 16M/64M phase limits and a shared 64MiB query ledger.
    pub(in crate::geometry) fn verify_step_member(
        &self,
        document: &multicad_rust_kernel::step_read::owned_insert::CompoundDocument<'_>,
        member: multicad_rust_kernel::step_read::owned_insert::Member,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<WoodInsertStepVerification, GeometryError> {
        use multicad_rust_kernel::step_read::npt_stock as read;
        let fail = |kind, message| error("clipped NPT STEP readback", kind, message);
        let source_error = |cause: npt::Error| {
            fail(
                match cause {
                    npt::Error::Invalid => GeometryErrorKind::Input,
                    npt::Error::Limit => GeometryErrorKind::Limit,
                    npt::Error::Cancelled => GeometryErrorKind::Cancelled,
                    npt::Error::Unresolved => GeometryErrorKind::Unsupported,
                },
                format!("Original-source correspondence failed: {cause:?}"),
            )
        };
        let read_error = |cause: read::Error| {
            fail(
                match cause {
                    read::Error::Invalid => GeometryErrorKind::InvalidResult,
                    read::Error::Limit => GeometryErrorKind::Limit,
                    read::Error::Cancelled => GeometryErrorKind::Cancelled,
                    read::Error::Unsupported => GeometryErrorKind::Unsupported,
                },
                format!("Actual member material failed: {cause:?}"),
            )
        };
        if self.pose != RigidTransform::default() || member.mapping_link().is_some() {
            return Err(fail(
                GeometryErrorKind::Unsupported,
                "Placed NPT file correspondence is not yet qualified".into(),
            ));
        }
        policy.validate()?;
        let options = npt::CurveOptions {
            position_error_mm: (policy.linear_tolerance_mm / 4.).min(1e-7),
            ..Default::default()
        };
        let curves = npt::CurvePlan::new(&self.source, options, interrupt).map_err(source_error)?;
        let plan =
            npt::SurfacePlan::new(curves, options.into(), interrupt).map_err(source_error)?;
        let available = (64_usize * 1024 * 1024)
            .checked_sub(plan.retained_budget_bytes() + 16 * 1024)
            .filter(|n| *n > 0)
            .ok_or_else(|| fail(GeometryErrorKind::Limit, "Member storage limit".into()))?;
        let selected = read::Selected::parse(
            document,
            member.item_entity(),
            read::Options {
                maximum_bytes: available,
                ..Default::default()
            },
            interrupt,
        )
        .map_err(read_error)?;
        let material = selected
            .retain(
                read::MaterialOptions {
                    maximum_bytes: available,
                    ..Default::default()
                },
                interrupt,
            )
            .map_err(read_error)?;
        let proof = plan
            .correspond_literal(
                &material,
                npt::CorrespondenceOptions {
                    position_error_mm: policy.linear_tolerance_mm,
                    ..Default::default()
                },
                interrupt,
            )
            .map_err(source_error)?;
        let bounds = material.local_enclosure_mm();
        let derivative = proof.original_derivative_error_mm();
        // The shared receipt retains its historical row/ruling field names.
        // Here they bound the original atlas u/v parameters, respectively;
        // both also conservatively include the finite-edge parameter bound.
        Ok(WoodInsertStepVerification {
            decoded_bounds_mm: [bounds.map(|b| b[0]), bounds.map(|b| b[1])],
            decoded_volume_bounds_mm3: material.local_volume_mm3(),
            source_volume_bounds_mm3: self.volume,
            maximum_position_error_mm: proof.position_error_mm(),
            maximum_row_derivative_error_mm: derivative[0].max(proof.edge_derivative_error_mm()),
            maximum_ruling_derivative_error_mm: derivative[1].max(proof.edge_derivative_error_mm()),
            verified_tolerance_mm: policy.linear_tolerance_mm,
        })
    }

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
