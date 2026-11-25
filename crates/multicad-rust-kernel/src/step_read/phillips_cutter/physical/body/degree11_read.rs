//! Independent actual-file qualification of distinctly repaired D11 material.
//! All original child geometry/byte/entity caps remain unchanged. Independent
//! reader limits are conservatively reserved; inherited D11 graph prefixes
//! are paid once through their nested aggregate. This is neither the old16M
//! reader contract nor a new whole-process memory/RSS guarantee. The returned
//! actual source is not a nominal recipe or exact unmodified raw STEP shell.
use super::*;
const STAGES: usize = 16;
const HARD: [usize; STAGES] = [
    1_000_000,
    10_000_000,
    1_000_000,
    10_000_000,
    16_000_000,
    10_000_000,
    16_000_000,
    16_000_000,
    16_000_000,
    16_000_000,
    16_000_000,
    16_000_000,
    64_000_000,
    64_000_000,
    512_000_000,
    512_000_000,
];
const FRONTIER_WORK: usize = 16_384;
const TOTAL: usize = 624_000_000 + FRONTIER_WORK;
const INPUT_BYTES: usize = 64 * 1024 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(usize)]
pub enum Degree11ReadStage {
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
    Carrier,
    Occupied,
    Material,
    Volume,
    Owned,
    /// Fresh supplied-author correspondence of the same decoded volume. Only
    /// the author-matching entry point starts it; its work is reported apart.
    AuthorMatch,
}
#[derive(Clone, Copy)]
pub struct Degree11ReadOptions {
    pub maximum_reserved_work: usize,
    pub maximum_input_bytes: usize,
    /// Ordered by Degree11ReadStage. Limits can only lower the fixed stage caps.
    pub stage_work_limits: [usize; STAGES],
}
impl Default for Degree11ReadOptions {
    fn default() -> Self {
        Self {
            maximum_reserved_work: TOTAL,
            maximum_input_bytes: INPUT_BYTES,
            stage_work_limits: HARD,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Degree11ReadError {
    pub stage: Option<Degree11ReadStage>,
    pub cause: PhillipsReadError,
}
pub struct Degree11MaterialRead {
    material: OwnedDegree11Material,
    reserved: usize,
}
impl Degree11MaterialRead {
    pub fn material(&self) -> &OwnedDegree11Material {
        &self.material
    }
    pub fn into_material(self) -> OwnedDegree11Material {
        self.material
    }
    pub fn reserved_delegated_work(&self) -> usize {
        self.reserved
    }
    /// Actual final nested proof/copy aggregate, with its unchanged512M cap.
    /// Head/carrier/occupied/material/volume prefixes are included once.
    pub fn material_proof_work_used(&self) -> usize {
        self.material.work_used()
    }
}
fn fail(cause: PhillipsReadError) -> Degree11ReadError {
    Degree11ReadError { stage: None, cause }
}
struct Context<'a> {
    interrupt: Option<&'a dyn Fn() -> bool>,
    progress: Option<&'a dyn Fn(Degree11ReadStage)>,
}
impl Context<'_> {
    fn stage<T>(
        &self,
        s: Degree11ReadStage,
        f: impl FnOnce() -> std::result::Result<T, PhillipsReadError>,
    ) -> std::result::Result<T, Degree11ReadError> {
        if self.interrupt.is_some_and(|c| c()) {
            return Err(Degree11ReadError {
                stage: Some(s),
                cause: PhillipsReadError::Cancelled,
            });
        }
        if let Some(progress) = self.progress {
            progress(s)
        }
        if self.interrupt.is_some_and(|c| c()) {
            return Err(Degree11ReadError {
                stage: Some(s),
                cause: PhillipsReadError::Cancelled,
            });
        }
        f().map_err(|cause| Degree11ReadError {
            stage: Some(s),
            cause,
        })
    }
}
pub fn read_degree11_phillips_material(
    input: &str,
    options: Degree11ReadOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> std::result::Result<Degree11MaterialRead, Degree11ReadError> {
    read_degree11_phillips_material_with_progress(input, options, interrupt, None)
}
/// Reports the actual stage being started. Completion is the returned owned
/// receipt; observing the last stage alone is not a success/admission signal.
pub fn read_degree11_phillips_material_with_progress(
    input: &str,
    options: Degree11ReadOptions,
    interrupt: Option<&dyn Fn() -> bool>,
    progress: Option<&dyn Fn(Degree11ReadStage)>,
) -> std::result::Result<Degree11MaterialRead, Degree11ReadError> {
    read_with(input, options, interrupt, progress, |_, _| Ok(())).map(|(read, ())| read)
}
/// Complete supplied-author/decoded-material correspondence of these bytes.
/// The match is earned on the SAME decoded volume the owned snapshot copies,
/// before that snapshot exists. Its fresh work is NOT inside the reader
/// reservation; callers reserve it separately through `matching.maximum_work`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Degree11AuthorMatchReceipt {
    pub maximum_world_distance_mm: f64,
    pub maximum_world_chart_derivative_mm: f64,
    /// Author minus independently decoded material.
    pub independent_volume_difference_mm3: [f64; 2],
    pub fresh_work: usize,
    pub total_work: usize,
    pub accounted_peak_bytes: usize,
    pub cavity_cells: usize,
}
pub fn read_degree11_phillips_material_for_author(
    input: &str,
    options: Degree11ReadOptions,
    author: &crate::phillips_drive::PhillipsScrew,
    frame: [Point; 3],
    matching: SuppliedDegree11MaterialOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> std::result::Result<(Degree11MaterialRead, Degree11AuthorMatchReceipt), Degree11ReadError> {
    read_with(input, options, interrupt, None, |volume, context| {
        context.stage(Degree11ReadStage::AuthorMatch, || {
            let matched = volume.match_supplied_author(author, frame, matching, interrupt)?;
            Ok(Degree11AuthorMatchReceipt {
                maximum_world_distance_mm: matched.maximum_world_distance_mm(),
                maximum_world_chart_derivative_mm: matched.maximum_world_chart_derivative_mm(),
                independent_volume_difference_mm3: matched.independent_volume_difference_mm3(),
                fresh_work: matched.new_work_used(),
                total_work: matched.total_work_used(),
                accounted_peak_bytes: matched.accounted_peak_bytes(),
                cavity_cells: matched.cavity_cells_used(),
            })
        })
    })
}
fn read_with<R>(
    input: &str,
    options: Degree11ReadOptions,
    interrupt: Option<&dyn Fn() -> bool>,
    progress: Option<&dyn Fn(Degree11ReadStage)>,
    check: impl FnOnce(
        &Degree11MaterialVolume<'_, '_, '_, '_>,
        &Context<'_>,
    ) -> std::result::Result<R, Degree11ReadError>,
) -> std::result::Result<(Degree11MaterialRead, R), Degree11ReadError> {
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
    // The first ten independent source readers retain their original limits.
    // Head + NEW carrier, then occupied/material/volume/copy form ONE nested
    // aggregate: reserve its maximum, not the sum of repeatedly inherited
    // prefixes. Frontier allowance covers fixed stage-entry refusal work.
    // Every term is below625M, including on a32-bit usize platform.
    let w = options.stage_work_limits;
    let aggregate = [w[10] + w[11], w[12], w[13], w[14], w[15]]
        .into_iter()
        .max()
        .unwrap();
    let reserved = w[..10].iter().sum::<usize>() + aggregate + FRONTIER_WORK;
    if reserved > options.maximum_reserved_work || input.len() > options.maximum_input_bytes {
        return Err(fail(PhillipsReadError::Limit));
    }
    let w = options.stage_work_limits;
    let body = context.stage(Degree11ReadStage::Body, || {
        read_phillips_body_with_interrupt(input, 1e-7, w[0], interrupt)
    })?;
    let outer = context.stage(Degree11ReadStage::Outer, || {
        body.outer_carrier_bounds(
            OuterCarrierOptions {
                maximum_distance_mm: 1e-7,
                work_limit: w[1],
            },
            interrupt,
        )
    })?;
    let cavity = body.cavity();
    let ends = context.stage(Degree11ReadStage::Endpoints, || {
        cavity.endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: w[2],
            },
            interrupt,
        )
    })?;
    let analytic = context.stage(Degree11ReadStage::Analytic, || {
        cavity.analytic_support_bounds(
            &ends,
            AnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: w[3],
            },
            interrupt,
        )
    })?;
    let lifts = context.stage(Degree11ReadStage::Lifts, || {
        analytic.finite_edge_lifts(
            ends,
            AnalyticEdgeLiftOptions {
                maximum_work: w[4],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let cubics = context.stage(Degree11ReadStage::Cubics, || {
        cavity.cubic_support_bounds(
            CubicSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: w[5],
            },
            interrupt,
        )
    })?;
    let sewn = context.stage(Degree11ReadStage::Sewn, || {
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
    let rep = context.stage(Degree11ReadStage::Representative, || {
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
    let closed = context.stage(Degree11ReadStage::ClosedCavity, || {
        rep.closed_boundary_certificate(
            crate::phillips_drive::ClosedBoundaryOptions {
                maximum_work: w[8],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let common = context.stage(Degree11ReadStage::CommonMouth, || {
        closed.common_mouth_head_envelope(
            &outer,
            CommonMouthOptions {
                maximum_work: w[9],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let head = context.stage(Degree11ReadStage::Head, || {
        common.closed_repaired_head(
            ClosedHeadOptions {
                maximum_work: w[10],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let carrier = context.stage(Degree11ReadStage::Carrier, || {
        head.degree11_thread_carrier(
            Degree11ThreadOptions {
                maximum_work: w[11],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let occupied = context.stage(Degree11ReadStage::Occupied, || {
        carrier.certify_occupied(
            Degree11OccupiedOptions {
                maximum_work: w[12],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let material = context.stage(Degree11ReadStage::Material, || {
        head.closed_degree11_material(
            &occupied,
            Degree11MaterialOptions {
                maximum_work: w[13],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let volume = context.stage(Degree11ReadStage::Volume, || {
        material.volume_bounds(
            Degree11MaterialVolumeOptions {
                maximum_work: w[14],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    let checked = check(&volume, &context)?;
    let owned = context.stage(Degree11ReadStage::Owned, || {
        volume.owned_snapshot(
            Degree11OwnedMaterialOptions {
                maximum_work: w[15],
                ..Default::default()
            },
            interrupt,
        )
    })?;
    if interrupt.is_some_and(|f| f()) {
        return Err(Degree11ReadError {
            stage: Some(Degree11ReadStage::Owned),
            cause: PhillipsReadError::Cancelled,
        });
    }
    Ok((
        Degree11MaterialRead {
            material: owned,
            reserved,
        },
        checked,
    ))
}

#[cfg(test)]
mod tests;
