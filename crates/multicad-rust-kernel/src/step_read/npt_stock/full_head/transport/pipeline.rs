//! One cumulative source-known mapped export gate; independent public import
//! remains the ordinary cold reader. No literal whole-head mesh is admitted.
use super::*;
use crate::{step::npt_stock::SurfacePlan, step_read::npt_stock::LiteralOccurrencePlacement};
#[cfg(test)]
mod tests;
const ORCHESTRATION_WORK: usize = 128;
#[derive(Clone, Copy, Debug)]
pub struct VerifiedMappedReceipt {
    pub encoding: encoding::EncodingReceipt,
    pub source: MappedSourceReceipt,
    pub orchestration_work: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
}
pub struct VerifiedMappedFile<'plan, 'source, 'pose> {
    actual: MappedActualHead,
    source: &'plan SurfacePlan<'source>,
    placement: &'pose LiteralOccurrencePlacement,
    head: HexHeadIntent,
    receipt: VerifiedMappedReceipt,
}
impl<'plan, 'source, 'pose> VerifiedMappedFile<'plan, 'source, 'pose> {
    pub fn actual(&self) -> &MappedActualHead {
        &self.actual
    }
    pub fn source_plan(&self) -> &'plan SurfacePlan<'source> {
        self.source
    }
    pub fn original_placement(&self) -> &'pose LiteralOccurrencePlacement {
        self.placement
    }
    pub fn head_definition(&self) -> HexHeadIntent {
        self.head
    }
    pub fn receipt(&self) -> VerifiedMappedReceipt {
        self.receipt
    }
    pub fn into_actual(self) -> MappedActualHead {
        self.actual
    }
}
/// Original local material, source plan and exact factor proof are inherited
/// ONCE. Encoder, actual-output binding/operator and raw/source/frame work all
/// consume ONE unchanged cumulative272M/64MiB envelope. Additional residents
/// exclude these explicit owners. No recipe or trace supplies output material.
pub fn encode_verified_mapped<'plan, 'source, 'pose>(
    local: Arc<OwnedActualHead>,
    source: &'plan SurfacePlan<'source>,
    head: HexHeadIntent,
    placement: &'pose LiteralOccurrencePlacement,
    options: SourceOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<VerifiedMappedFile<'plan, 'source, 'pose>> {
    super::super::source::validate_source_options(head, options)?;
    if interrupt.is_some_and(|f| f()) {
        return Err(Error::Cancelled);
    }
    let r = options.resources;
    // Fixed validation, stage-limit arithmetic, Arc publication and result
    // assembly are reserved BEFORE entering any child or allocating output.
    let available = r
        .maximum_work
        .checked_sub(ORCHESTRATION_WORK)
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)?;
    let outside = r
        .additional_resident_bytes
        .checked_add(size_of::<VerifiedMappedFile<'_, '_, '_>>())
        .ok_or(Error::Limit)?;
    let candidate = encoding::propose_mapped(
        &local,
        placement,
        encoding::EncodingOptions {
            maximum_work: available
                .checked_sub(source.work_used())
                .filter(|n| *n > 0)
                .ok_or(Error::Limit)?,
            maximum_bytes: r.maximum_bytes,
            additional_resident_bytes: outside
                .checked_add(source.retained_budget_bytes())
                .ok_or(Error::Limit)?,
        },
        interrupt,
    )?;
    let encoding = candidate.receipt();
    let bytes = candidate.actual_bytes();
    drop(candidate);
    let actual = MappedActualHead::read_reusing_local(
        bytes,
        local,
        OwnedHeadOptions {
            resources: Options {
                maximum_work: available
                    .checked_sub(source.work_used())
                    .and_then(|n| n.checked_sub(placement.work_used()))
                    .and_then(|n| n.checked_sub(encoding.fresh_encoding_work))
                    .filter(|n| *n > 0)
                    .ok_or(Error::Limit)?,
                maximum_bytes: r.maximum_bytes,
                additional_resident_bytes: outside
                    .checked_add(source.retained_budget_bytes())
                    .and_then(|n| n.checked_add(placement.retained_budget_bytes()))
                    .ok_or(Error::Limit)?,
            },
            material_tolerance_mm: options.material_tolerance_mm,
            minimum_wall_mm: options.minimum_wall_mm,
        },
        interrupt,
    )?;
    // Borrowed qualification is scoped here; the owning result below retains
    // actual bytes/material plus explicit ORIGINAL source/frame associations.
    let matched = {
        let qualified = actual.qualify_original_source(
            source,
            head,
            placement,
            SourceOptions {
                resources: Options {
                    maximum_work: available
                        .checked_sub(encoding.fresh_encoding_work)
                        .filter(|n| *n > 0)
                        .ok_or(Error::Limit)?,
                    maximum_bytes: r.maximum_bytes,
                    additional_resident_bytes: outside,
                },
                ..options
            },
            interrupt,
        )?;
        qualified.receipt()
    };
    let aggregate = matched
        .aggregate_work
        .checked_add(encoding.fresh_encoding_work)
        .and_then(|n| n.checked_add(ORCHESTRATION_WORK))
        .ok_or(Error::Limit)?;
    let peak = encoding
        .accounted_peak_bytes
        .max(matched.accounted_peak_bytes);
    if aggregate > r.maximum_work || peak > r.maximum_bytes {
        return Err(Error::Limit);
    }
    if interrupt.is_some_and(|f| f()) {
        return Err(Error::Cancelled);
    }
    Ok(VerifiedMappedFile {
        actual,
        source,
        placement,
        head,
        receipt: VerifiedMappedReceipt {
            encoding,
            source: matched,
            orchestration_work: ORCHESTRATION_WORK,
            aggregate_work: aggregate,
            accounted_peak_bytes: peak,
        },
    })
}
