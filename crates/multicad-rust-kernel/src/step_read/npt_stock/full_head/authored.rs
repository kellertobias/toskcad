//! Authored clipped-NPT hex-head screw export. The authored stock clip plan is
//! written, its own physical head is earned and encoded inside one mapped item
//! carrying the authored occurrence chain, and the cumulative gate qualifies
//! the actual mapped bytes against the same authored source, head intent and
//! chain. A separate entry point re-proves any actual file the same way.
//! Volumes are the sewn reference, not raw occupied STEP mass.
use super::transport::{MappedActualHead, MappedSourceReceipt, VerifiedMappedReceipt};
use super::*;
use crate::{
    axial_pipe_thread::authored::StockClipPlan,
    step::npt_stock::{CurveOptions, CurvePlan, SurfacePlan},
    step_read::npt_stock::{LiteralOccurrencePlacement, head},
};
use std::sync::Arc;
/// One rigid factor: row-major matrix and translation (world = M p + t).
pub type OccurrenceFactor = ([[f64; 3]; 3], [f64; 3]);
#[derive(Clone, Copy, Debug)]
pub struct AuthoredHeadReceipt {
    pub source: MappedSourceReceipt,
    /// Sewn reference of the actual mapped file, never raw occupied mass.
    pub world_sewn_volume_mm3: [f64; 2],
    pub world_sewn_enclosure_mm: [[f64; 2]; 3],
}
pub struct AuthoredHeadExport {
    text: String,
    receipt: VerifiedMappedReceipt,
    volume: [f64; 2],
    enclosure: [[f64; 2]; 3],
}
impl AuthoredHeadExport {
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn into_string(self) -> String {
        self.text
    }
    pub fn receipt(&self) -> VerifiedMappedReceipt {
        self.receipt
    }
    pub fn world_sewn_volume_mm3(&self) -> [f64; 2] {
        self.volume
    }
    pub fn world_sewn_enclosure_mm(&self) -> [[f64; 2]; 3] {
        self.enclosure
    }
}
fn npt(e: crate::step::npt_stock::Error) -> Error {
    source::source_error(e)
}
fn with_plan<R>(
    clip: &StockClipPlan,
    interrupt: Option<&dyn Fn() -> bool>,
    f: impl FnOnce(&SurfacePlan<'_>) -> Result<R>,
) -> Result<R> {
    let options = CurveOptions::default();
    let curves = CurvePlan::new(clip, options, interrupt).map_err(npt)?;
    let plan = SurfacePlan::new(curves, options.into(), interrupt).map_err(npt)?;
    f(&plan)
}
fn placement(
    factors: &[OccurrenceFactor],
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<LiteralOccurrencePlacement> {
    let (first, rest) = factors.split_first().ok_or(Error::Invalid)?;
    let mut p = LiteralOccurrencePlacement::new(first.0, first.1, Default::default(), interrupt)?;
    for (m, t) in rest {
        p = p.transformed(*m, *t, Default::default(), interrupt)?;
    }
    Ok(p)
}
/// Earn the local physical head from the written authored stock.
fn local_head(
    plan: &SurfacePlan<'_>,
    head: HexHeadIntent,
    options: SourceOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<OwnedActualHead> {
    let stock = crate::step::npt_stock::write_source_candidate(
        plan,
        CurveOptions::default().maximum_work,
        CurveOptions::default().maximum_bytes,
        interrupt,
    )
    .map_err(npt)?
    .into_string();
    let doc = CompoundDocument::parse(&stock, DocumentOptions::default(), interrupt)
        .map_err(neutral_error)?;
    let root = doc
        .represented_members()
        .first()
        .ok_or(Error::Invalid)?
        .item_entity();
    let parsed = Parsed::parse(&doc, root, Default::default(), interrupt)?;
    // The auxiliary verification cap is the highest planar face: the neck.
    let cap = parsed
        .faces
        .iter()
        .filter_map(|f| match f.support {
            Support::Plane { height } => Some((parsed.points[height].coordinates[2].value(), f.id)),
            _ => None,
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .ok_or(Error::Unsupported)?
        .1;
    let input = head::HeadInput::regular_hex(&parsed, cap, head.across_flats_mm, head.height_mm)?;
    let proof = head::certify(
        &parsed,
        input,
        options.material_tolerance_mm,
        options.minimum_wall_mm,
        Default::default(),
        interrupt,
    )?;
    let local = crate::step::npt_stock::full_head::propose(&proof, Default::default(), interrupt)?;
    OwnedActualHead::read(
        Arc::from(local.as_str()),
        OwnedHeadOptions {
            material_tolerance_mm: options.material_tolerance_mm,
            minimum_wall_mm: options.minimum_wall_mm,
            ..Default::default()
        },
        interrupt,
    )
}
/// Write and qualify one authored clipped hex-head screw in its occurrence.
pub fn export_authored_hex_head(
    clip: &StockClipPlan,
    head: HexHeadIntent,
    factors: &[OccurrenceFactor],
    options: SourceOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<AuthoredHeadExport> {
    source::validate_source_options(head, options)?;
    with_plan(clip, interrupt, |plan| {
        let local = Arc::new(local_head(plan, head, options, interrupt)?);
        let pose = placement(factors, interrupt)?;
        let verified =
            transport::encode_verified_mapped(local, plan, head, &pose, options, interrupt)?;
        let actual = verified.actual();
        Ok(AuthoredHeadExport {
            text: actual.actual_file().to_owned(),
            receipt: verified.receipt(),
            volume: actual.world_sewn_volume_mm3(),
            enclosure: actual.world_sewn_enclosure_mm(),
        })
    })
}
/// Re-prove ACTUAL mapped bytes against the authored source, head and chain.
pub fn qualify_authored_hex_head_file(
    actual: &str,
    clip: &StockClipPlan,
    head: HexHeadIntent,
    factors: &[OccurrenceFactor],
    options: SourceOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<AuthoredHeadReceipt> {
    source::validate_source_options(head, options)?;
    with_plan(clip, interrupt, |plan| {
        let pose = placement(factors, interrupt)?;
        let read = MappedActualHead::read(
            Arc::from(actual),
            OwnedHeadOptions {
                material_tolerance_mm: options.material_tolerance_mm,
                minimum_wall_mm: options.minimum_wall_mm,
                ..Default::default()
            },
            interrupt,
        )?;
        let qualified = read.qualify_original_source(plan, head, &pose, options, interrupt)?;
        Ok(AuthoredHeadReceipt {
            source: qualified.receipt(),
            world_sewn_volume_mm3: read.world_sewn_volume_mm3(),
            world_sewn_enclosure_mm: read.world_sewn_enclosure_mm(),
        })
    })
}
#[cfg(test)]
#[path = "authored_tests.rs"]
mod tests;
