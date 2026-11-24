//! Native STEP admission requires independent literal material and measured
//! correspondence to this unchanged source. A placed occurrence is written
//! as one map of the local chart, and readback bounds that decoded map
//! against the occurrence before admitting it.
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
        // A placed occurrence must be written as one explicit map of the
        // unchanged local chart; an unplaced one must not carry any map.
        if (self.pose != RigidTransform::default()) != member.mapping_link().is_some() {
            return Err(fail(
                GeometryErrorKind::InvalidResult,
                "Actual member mapping does not match the occurrence placement".into(),
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
        let derivative = proof.original_derivative_error_mm();
        let edge = proof.edge_derivative_error_mm();
        let local = [
            proof.position_error_mm(),
            derivative[0].max(edge),
            derivative[1].max(edge),
        ];
        let [position, row, ruling] = if self.pose == RigidTransform::default() {
            local
        } else {
            self.transport(&material, local)
                .map_err(|message| fail(GeometryErrorKind::InvalidResult, message.into()))?
        };
        if position > policy.linear_tolerance_mm {
            return Err(fail(
                GeometryErrorKind::InvalidResult,
                "Actual STEP placement exceeds the source position tolerance".into(),
            ));
        }
        // World values follow the file's own decoded operator; the identity
        // member's world and local charts coincide.
        let bounds = material.world_enclosure_mm();
        // The shared receipt retains its historical row/ruling field names.
        // Here they bound the original atlas u/v parameters, respectively;
        // both also conservatively include the finite-edge parameter bound.
        Ok(WoodInsertStepVerification {
            decoded_bounds_mm: [bounds.map(|b| b[0]), bounds.map(|b| b[1])],
            decoded_volume_bounds_mm3: material.world_volume_mm3(),
            source_volume_bounds_mm3: self.placed_volume(),
            maximum_position_error_mm: position,
            maximum_row_derivative_error_mm: row,
            maximum_ruling_derivative_error_mm: ruling,
            verified_tolerance_mm: policy.linear_tolerance_mm,
        })
    }

    /// Bound the placed file against the placed source from the local proof.
    ///
    /// With the decoded file map `x -> A x + o` and the occurrence
    /// `x -> R x + t`, a file point `A F + o` differs from the source point
    /// `R S + t` by `R (F - S) + (A - R) F + (o - t)`. `|R y| <=
    /// sqrt(1 + m) |y|` for the occurrence's metric error `m`, and `F` lies
    /// in the decoded local enclosure, so the position error is at most
    /// `sqrt(1 + m) e + |(o - t) + (A - R) box|`. Derivatives transport the
    /// same way without the offset. Every decoded surface cell's local
    /// derivative hull, scaled by its uniform span count to the global
    /// parameter, bounds `F'`; edges are surface boundaries in those same
    /// parameters, and the planar caps' unit tangents add the X/Y columns.
    fn transport(
        &self,
        material: &multicad_rust_kernel::step_read::npt_stock::OwnedMaterial,
        [position, row, ruling]: [f64; 3],
    ) -> Result<[f64; 3], &'static str> {
        use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
        let interval = |[lo, hi]: [f64; 2]| I::new(lo, hi);
        let norm = |v: [I; 3]| -> Result<f64, &'static str> {
            let sum = v.into_iter().fold(I::point(0.), |s, x| s + x.square());
            let sum = sum.checked().ok_or("Placed STEP enclosure overflow")?;
            I::new(sum.lo.max(0.), sum.hi)
                .sqrt()
                .checked()
                .map(|r| r.hi)
                .ok_or("Placed STEP enclosure overflow")
        };
        let actual = material.linear_columns_bounds().map(|c| c.map(interval));
        let rotation = self.pose.rotation;
        let delta: [[I; 3]; 3] = std::array::from_fn(|j| {
            std::array::from_fn(|i| actual[j][i] - I::point(rotation[i][j]))
        });
        let apply = |d: [I; 3]| -> [I; 3] {
            std::array::from_fn(|i| (0..3).fold(I::point(0.), |s, j| s + delta[j][i] * d[j]))
        };
        let origin = material.origin_mm().map(interval);
        let local = material.local_enclosure_mm().map(interval);
        let moved = apply(local);
        let offset = norm(std::array::from_fn(|i| {
            origin[i] - I::point(self.pose.translation_mm[i]) + moved[i]
        }))?;
        let stretch = (I::point(1.) + I::point(self.calibration().metric_error()))
            .sqrt()
            .checked()
            .ok_or("Placed STEP enclosure overflow")?;
        let mut derivative = norm(apply([I::point(1.), I::point(0.), I::point(0.)]))?
            .max(norm(apply([I::point(0.), I::point(1.), I::point(0.)]))?);
        for (t, tensor) in material.tensors().iter().enumerate() {
            let spans = tensor.spans();
            for u in 0..spans[0] {
                for v in 0..spans[1] {
                    let unresolved = "Decoded tensor derivative bounds are unresolved";
                    let cell = material.cell(t, [u, v]).map_err(|_| unresolved)?;
                    let jet = cell.first_derivative_bounds().map_err(|_| unresolved)?;
                    for (axis, d) in [jet.du, jet.dv].into_iter().enumerate() {
                        let scale = I::point(spans[axis] as f64);
                        let d = d.map(|b| interval(b) * scale);
                        derivative = derivative.max(norm(apply(d))?);
                    }
                }
            }
        }
        let sum = |a: f64, b: I| (I::point(a) * stretch + b).checked().map(|r| r.hi);
        let derivative = I::point(derivative);
        match (
            sum(position, I::point(offset)),
            sum(row, derivative),
            sum(ruling, derivative),
        ) {
            (Some(position), Some(row), Some(ruling)) => Ok([position, row, ruling]),
            _ => Err("Placed STEP enclosure overflow"),
        }
    }

    /// The occurrence map written around the local chart.
    fn occurrence_map(&self) -> npt::OccurrenceMap {
        let r = self.pose.rotation;
        npt::OccurrenceMap {
            columns: std::array::from_fn(|j| std::array::from_fn(|i| r[i][j])),
            origin_mm: self.pose.translation_mm,
        }
    }

    fn with_step_plan<T>(
        &self,
        policy: GeometryPolicy,
        f: impl FnOnce(&npt::SurfacePlan<'_>) -> Result<T, npt::Error>,
    ) -> Result<T, step::StepError> {
        // The plan is the unchanged LOCAL chart. Native's literal near-rigid
        // map and STEP's normalized operator are distinct maps: a placed
        // member is admitted only after readback bounds their difference.
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
        if self.pose != RigidTransform::default() && !self.occurrence_map().admitted() {
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
            if self.pose == RigidTransform::default() {
                writer.append_verified_source_npt_candidate(&candidate, options, None)
            } else {
                writer.append_placed_verified_source_npt_candidate(
                    &candidate,
                    options,
                    self.occurrence_map(),
                    None,
                )
            }
            .map(|_| ())
        })
    }
}
