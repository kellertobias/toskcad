//! Explicit staged ACTUAL-FILE qualification. Existing stage proof/byte/entity
//! caps are unchanged. The caller separately reserves the SUM of delegated
//! work caps before any reader runs; reused prefixes are conservatively paid
//! again in that reservation. This is not the old16M reader contract and does
//! not claim a new aggregate process-memory bound. Each original source-stage
//! memory contract remains in force. Final geometry/mass are distinctly repaired,
//! not exact occupied material of the unmodified literal STEP surfaces.
use super::*;
const STAGES: usize = 17;
const HARD: [usize; STAGES] = [
    1_000_000, 10_000_000, 1_000_000, 10_000_000, 16_000_000, 10_000_000, 16_000_000, 16_000_000,
    16_000_000, 16_000_000, 16_000_000, 16_000_000, 16_000_000, 16_000_000, 16_000_000, 16_000_000,
    16_000_000,
];
const TOTAL: usize = 256_000_000;
const INPUT_BYTES: usize = 64 * 1024 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(usize)]
pub enum PeriodicReadStage {
    Body,
    Outer,
    Endpoints,
    Analytic,
    Lifts,
    Cubics,
    Sewn,
    Representative,
    ClosedCavity,
    CommonMouth,
    Head,
    Periodic,
    Finite,
    Shoulder,
    Material,
    Volume,
    Owned,
}
#[derive(Clone, Copy)]
pub struct PeriodicReadOptions {
    pub maximum_reserved_work: usize,
    pub maximum_input_bytes: usize,
    /// Ordered by PeriodicReadStage. Limits can only lower the fixed stage caps.
    pub stage_work_limits: [usize; STAGES],
}
impl Default for PeriodicReadOptions {
    fn default() -> Self {
        Self {
            maximum_reserved_work: TOTAL,
            maximum_input_bytes: INPUT_BYTES,
            stage_work_limits: HARD,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PeriodicReadError {
    pub stage: Option<PeriodicReadStage>,
    pub cause: PhillipsReadError,
}
pub struct PeriodicMaterialRead {
    material: OwnedPeriodicMaterial,
    reserved: usize,
}
impl PeriodicMaterialRead {
    pub fn material(&self) -> &OwnedPeriodicMaterial {
        &self.material
    }
    pub fn into_material(self) -> OwnedPeriodicMaterial {
        self.material
    }
    pub fn reserved_delegated_work(&self) -> usize {
        self.reserved
    }
    /// Existing child receipt, with its original aggregate64M proof cap. The
    /// independent body-reader reservation is included in reserved_delegated_work.
    pub fn material_proof_work_used(&self) -> usize {
        self.material.total_work_used()
    }
}
fn fail(cause: PhillipsReadError) -> PeriodicReadError {
    PeriodicReadError { stage: None, cause }
}
struct Context<'a> {
    interrupt: Option<&'a dyn Fn() -> bool>,
    progress: Option<&'a dyn Fn(PeriodicReadStage)>,
}
impl Context<'_> {
    fn stage<T>(
        &self,
        s: PeriodicReadStage,
        f: impl FnOnce() -> std::result::Result<T, PhillipsReadError>,
    ) -> std::result::Result<T, PeriodicReadError> {
        if self.interrupt.is_some_and(|c| c()) {
            return Err(PeriodicReadError {
                stage: Some(s),
                cause: PhillipsReadError::Cancelled,
            });
        }
        if let Some(progress) = self.progress {
            progress(s)
        }
        if self.interrupt.is_some_and(|c| c()) {
            return Err(PeriodicReadError {
                stage: Some(s),
                cause: PhillipsReadError::Cancelled,
            });
        }
        f().map_err(|cause| PeriodicReadError {
            stage: Some(s),
            cause,
        })
    }
}
pub fn read_periodic_phillips_material(
    input: &str,
    options: PeriodicReadOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> std::result::Result<PeriodicMaterialRead, PeriodicReadError> {
    read_periodic_phillips_material_with_progress(input, options, interrupt, None)
}
/// Reports the actual stage being started. Completion is the returned owned
/// receipt; observing the last stage alone is not a success/admission signal.
pub fn read_periodic_phillips_material_with_progress(
    input: &str,
    options: PeriodicReadOptions,
    interrupt: Option<&dyn Fn() -> bool>,
    progress: Option<&dyn Fn(PeriodicReadStage)>,
) -> std::result::Result<PeriodicMaterialRead, PeriodicReadError> {
    if interrupt.is_some_and(|f| f()) {
        return Err(fail(PhillipsReadError::Cancelled));
    }
    let context = Context {
        interrupt,
        progress,
    };
    if options.maximum_reserved_work > TOTAL
        || options.maximum_input_bytes > INPUT_BYTES
        || options
            .stage_work_limits
            .iter()
            .zip(HARD)
            .any(|(&n, hard)| n == 0 || n > hard)
    {
        return Err(fail(PhillipsReadError::Invalid));
    }
    //17 fixed validated limits sum to at most224M, below the32-bit usize bound.
    let reserved: usize = options.stage_work_limits.iter().sum();
    if reserved > options.maximum_reserved_work || input.len() > options.maximum_input_bytes {
        return Err(fail(PhillipsReadError::Limit));
    }
    let w = options.stage_work_limits;
    let body = context.stage(PeriodicReadStage::Body, || {
        read_phillips_body_with_interrupt(input, 1e-7, w[0], interrupt)
    })?;
    let outer = context.stage(PeriodicReadStage::Outer, || {
        body.outer_carrier_bounds(
            OuterCarrierOptions {
                maximum_distance_mm: 1e-7,
                work_limit: w[1],
            },
            interrupt,
        )
    })?;
    let cavity = body.cavity();
    let ends = context.stage(PeriodicReadStage::Endpoints, || {
        cavity.endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: w[2],
            },
            interrupt,
        )
    })?;
    let analytic = context.stage(PeriodicReadStage::Analytic, || {
        cavity.analytic_support_bounds(
            &ends,
            AnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: w[3],
            },
            interrupt,
        )
    })?;
    let lifts = context.stage(PeriodicReadStage::Lifts, || {
        analytic.finite_edge_lifts(
            ends,
            AnalyticEdgeLiftOptions {
                maximum_work: w[4],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let cubics = context.stage(PeriodicReadStage::Cubics, || {
        cavity.cubic_support_bounds(
            CubicSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: w[5],
            },
            interrupt,
        )
    })?;
    let sewn = context.stage(PeriodicReadStage::Sewn, || {
        lifts.sewn_edge_sources(
            &body,
            &cubics,
            SewnEdgeOptions {
                maximum_work: w[6],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let rep = context.stage(PeriodicReadStage::Representative, || {
        body.derived_cavity_representative(
            &sewn,
            DerivedCavityOptions {
                maximum_work: w[7],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    // The representative retains Body and newly earned source/proof data;
    // it does not borrow the prerequisite buffers. Release those immediately.
    drop(sewn);
    drop(lifts);
    drop(cubics);
    let closed = context.stage(PeriodicReadStage::ClosedCavity, || {
        rep.closed_boundary_certificate(
            crate::phillips_drive::ClosedBoundaryOptions {
                maximum_work: w[8],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let common = context.stage(PeriodicReadStage::CommonMouth, || {
        closed.common_mouth_head_envelope(
            &outer,
            CommonMouthOptions {
                maximum_work: w[9],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let head = context.stage(PeriodicReadStage::Head, || {
        common.closed_repaired_head(
            ClosedHeadOptions {
                maximum_work: w[10],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let field = context.stage(PeriodicReadStage::Periodic, || {
        head.common_chart_periodic_thread_repair(
            AnalyticThreadOptions {
                maximum_work: w[11],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let finite = context.stage(PeriodicReadStage::Finite, || {
        head.periodic_thread_finite_curves(
            &field,
            PeriodicFiniteOptions {
                maximum_work: w[12],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let join = context.stage(PeriodicReadStage::Shoulder, || {
        head.periodic_shared_crest_shoulder(
            &finite,
            PeriodicShoulderOptions {
                maximum_work: w[13],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let material = context.stage(PeriodicReadStage::Material, || {
        head.periodic_closed_material(
            &join,
            PeriodicMaterialOptions {
                maximum_work: w[14],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let volume = context.stage(PeriodicReadStage::Volume, || {
        material.volume_bounds(
            PeriodicVolumeOptions {
                maximum_work: w[15],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let owned = context.stage(PeriodicReadStage::Owned, || {
        volume.owned_source(
            PeriodicOwnedOptions {
                maximum_work: w[16],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    if interrupt.is_some_and(|f| f()) {
        return Err(PeriodicReadError {
            stage: Some(PeriodicReadStage::Owned),
            cause: PhillipsReadError::Cancelled,
        });
    }
    Ok(PeriodicMaterialRead {
        material: owned,
        reserved,
    })
}

#[cfg(test)]
mod tests;
