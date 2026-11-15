//! Atomic source-qualified bytes with independently decoded healed material.
//! Original affine-source mass and healed mass are deliberately separate.
use super::*;
use crate::step_read::{
    owned_insert::{CompoundDocument, DocumentOptions},
    owned_orthogonal_seats as read,
};
use std::collections::BTreeSet;
mod replay;
mod verify;
pub use verify::{VerificationOptions, VerifiedMember};
mod snapshot;
pub use replay::{AppendOptions, QualifiedAppendReceipt};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy)]
pub struct QualificationOptions {
    /// Finite whole-operation budget, including the one original source
    /// constructor prerequisite. Individual unchanged query caps still apply.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub material: read::MaterialOptions,
}
impl Default for QualificationOptions {
    fn default() -> Self {
        Self {
            maximum_work: 64_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            material: read::MaterialOptions::default(),
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QualificationError {
    Invalid,
    Unsupported,
    Limit,
    Cancelled,
    Encoding(StepError),
}
impl From<StepError> for QualificationError {
    fn from(e: StepError) -> Self {
        Self::Encoding(e)
    }
}
impl From<read::Error> for QualificationError {
    fn from(e: read::Error) -> Self {
        match e {
            read::Error::Invalid => Self::Invalid,
            read::Error::Unsupported => Self::Unsupported,
            read::Error::Limit => Self::Limit,
            read::Error::Cancelled => Self::Cancelled,
        }
    }
}
impl From<super::Error> for QualificationError {
    fn from(e: super::Error) -> Self {
        match e {
            super::Error::Invalid => Self::Invalid,
            super::Error::Limit => Self::Limit,
            super::Error::Cancelled => Self::Cancelled,
            super::Error::Encoding(e) => Self::Encoding(e),
        }
    }
}
fn document_error(e: crate::step_read::owned_insert::Error) -> QualificationError {
    use crate::step_read::owned_insert::Error as E;
    match e {
        E::Invalid => QualificationError::Invalid,
        E::Unsupported => QualificationError::Unsupported,
        E::Limit => QualificationError::Limit,
        E::Cancelled => QualificationError::Cancelled,
    }
}
struct Ledger<'a> {
    used: usize,
    left: usize,
    maximum_bytes: usize,
    peak: usize,
    interrupt: &'a dyn Fn() -> bool,
}
impl Ledger<'_> {
    fn charge(&mut self, n: usize) -> Result<(), QualificationError> {
        if (self.interrupt)() {
            return Err(QualificationError::Cancelled);
        }
        self.left = self.left.checked_sub(n).ok_or(QualificationError::Limit)?;
        self.used = self.used.checked_add(n).ok_or(QualificationError::Limit)?;
        Ok(())
    }
    fn bytes(&mut self, n: usize) -> Result<(), QualificationError> {
        if n > self.maximum_bytes {
            return Err(QualificationError::Limit);
        }
        self.peak = self.peak.max(n);
        Ok(())
    }
    fn remaining_bytes(&self, live: usize) -> Result<usize, QualificationError> {
        self.maximum_bytes
            .checked_sub(live)
            .filter(|n| *n > 0)
            .ok_or(QualificationError::Limit)
    }
    fn remaining_work(&self, cap: usize) -> Result<usize, QualificationError> {
        let n = self.left.min(cap);
        if n == 0 {
            Err(QualificationError::Limit)
        } else {
            Ok(n)
        }
    }
}
pub struct QualifiedOrthogonalSeatsStep {
    text: String,
    source: Vec<u64>,
    placement: StockClipPlacement,
    quality: read::MaterialOptions,
    original_volume: [f64; 2],
    healed_volume: [f64; 2],
    original: read::Correspondence,
    healed: read::Correspondence,
    composed_position: f64,
    composed_derivative: f64,
    work: usize,
    encoding_work: usize,
    peak: usize,
}
impl QualifiedOrthogonalSeatsStep {
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn into_string(self) -> String {
        self.text
    }
    /// Exact finite source inputs, local owner order and actual stored pose;
    /// this does not replace checking the requested material quality policy.
    pub fn matches_source(&self, source: &Seats, placement: StockClipPlacement) -> bool {
        snapshot::words(source, placement).is_ok_and(|words| self.source == words)
    }
    pub fn placement(&self) -> StockClipPlacement {
        self.placement
    }
    pub fn qualification_material_options(&self) -> read::MaterialOptions {
        self.quality
    }
    pub fn original_source_volume_bounds_mm3(&self) -> [f64; 2] {
        self.original_volume
    }
    pub fn healed_volume_bounds_mm3(&self) -> [f64; 2] {
        self.healed_volume
    }
    pub fn source_to_raw_correspondence(&self) -> read::Correspondence {
        self.original
    }
    pub fn raw_to_healed_correspondence(&self) -> read::Correspondence {
        self.healed
    }
    pub fn source_to_healed_position_error_mm(&self) -> f64 {
        self.composed_position
    }
    pub fn source_to_healed_normalized_derivative_error_mm(&self) -> f64 {
        self.composed_derivative
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn encoding_work_used(&self) -> usize {
        self.encoding_work
    }
    /// Conservative phase-wide allocation accounting; not process RSS.
    pub fn peak_accounted_bytes(&self) -> usize {
        self.peak
    }
}
fn source_scope(source: &Seats, w: &mut Ledger<'_>) -> Result<(), QualificationError> {
    if source.boundary().faces.len() > 8 {
        return Err(QualificationError::Unsupported);
    }
    let mut points = BTreeSet::new();
    let mut uses = 0usize;
    for p in source
        .boundary()
        .faces
        .iter()
        .flat_map(|f| &f.rings)
        .flatten()
    {
        w.charge(64)?;
        uses += 1;
        points.insert(point_key(*p));
        if points.len() > 16 || uses > 48 {
            return Err(QualificationError::Unsupported);
        }
    }
    Ok(())
}
/// No bytes or receipt escape on cancellation, budget refusal, unsupported
/// source correspondence, or independent material failure.
pub fn qualify(
    source: &Seats,
    placement: StockClipPlacement,
    options: QualificationOptions,
    interrupt: &dyn Fn() -> bool,
) -> Result<QualifiedOrthogonalSeatsStep, QualificationError> {
    let mut w = Ledger {
        used: 0,
        left: options.maximum_work,
        maximum_bytes: options.maximum_bytes,
        peak: 0,
        interrupt,
    };
    w.charge(0)?;
    if options.maximum_work == 0
        || options.maximum_work > 256_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
    {
        return Err(QualificationError::Invalid);
    }
    options.material.validate()?;
    w.charge(source.work_used)?;
    w.bytes(source.retained_bytes())?;
    source_scope(source, &mut w)?;
    // All admitted source arrays are already constructor-bounded. The exact
    // words contain <64KiB; reserve their pinned Vec growth before capture.
    w.bytes(
        source
            .retained_bytes()
            .checked_add(128 * 1024)
            .ok_or(QualificationError::Limit)?,
    )?;
    let snapshot = snapshot::words(source, placement)?;
    w.charge(snapshot.len() * 4)?;
    let snapshot_bytes = snapshot.capacity() * std::mem::size_of::<u64>();
    if snapshot_bytes > 128 * 1024 {
        return Err(QualificationError::Limit);
    }
    let candidate = super::candidate(
        source,
        placement,
        Options {
            maximum_work: w.remaining_work(16_000_000)?,
            maximum_bytes: w.remaining_bytes(snapshot_bytes)?,
        },
        interrupt,
    )?;
    w.charge(candidate.encoding_work())?;
    w.bytes(snapshot_bytes + candidate.accounted_peak_bytes())?;
    let live = source.retained_bytes()
        + snapshot_bytes
        + candidate.text.capacity()
        + candidate.faces.capacity() * std::mem::size_of::<u32>();
    w.bytes(live)?;
    // Generated local schema is bounded above by8 stock faces,16 vertices,
    // 24 stock edges,16 bores and32 seats. 128 bytes per input byte bounds the
    // parser's live Value/component/BTree growth before any arena is allocated.
    let parse = live
        .checked_add(
            candidate
                .as_str()
                .len()
                .checked_mul(128)
                .ok_or(QualificationError::Limit)?,
        )
        .ok_or(QualificationError::Limit)?;
    w.bytes(parse)?;
    let document = CompoundDocument::parse(
        candidate.as_str(),
        DocumentOptions {
            maximum_work: w.remaining_work(512_000_000)?,
            retained_budget_bytes: w
                .remaining_bytes(live)?
                .checked_add(candidate.as_str().len())
                .ok_or(QualificationError::Limit)?,
        },
        Some(interrupt),
    )
    .map_err(document_error)?;
    w.charge(document.work_used())?;
    // The borrowed candidate text is already live, so its bytes are not charged
    // twice through the document's retained-input accounting.
    let doc_extra = document
        .retained_budget_bytes()
        .checked_sub(candidate.as_str().len())
        .ok_or(QualificationError::Invalid)?;
    let live = live
        .checked_add(doc_extra)
        .ok_or(QualificationError::Limit)?;
    w.bytes(live)?;
    let [member] = document.represented_members() else {
        return Err(QualificationError::Invalid);
    };
    let boundary = read::decode(
        &document,
        member.item_entity(),
        read::Options {
            maximum_work: w.remaining_work(16_000_000)?,
            maximum_bytes: w.remaining_bytes(live)?,
        },
        interrupt,
    )?;
    w.charge(boundary.work_used())?;
    let live = live
        .checked_add(boundary.query_budget_bytes())
        .ok_or(QualificationError::Limit)?;
    w.bytes(live)?;
    let mut material_options = options.material;
    material_options.transfer.maximum_work = material_options
        .transfer
        .maximum_work
        .min(w.remaining_work(16_000_000)?);
    material_options.transfer.maximum_bytes = material_options
        .transfer
        .maximum_bytes
        .min(w.remaining_bytes(live)?);
    let material = boundary.heal_material(material_options, interrupt)?;
    w.charge(material.work_used())?;
    let live = live
        .checked_add(material.query_budget_bytes())
        .ok_or(QualificationError::Limit)?;
    w.bytes(live)?;
    let mut transfer = options.material.transfer;
    transfer.maximum_work = transfer.maximum_work.min(w.remaining_work(16_000_000)?);
    transfer.maximum_bytes = transfer.maximum_bytes.min(w.remaining_bytes(live)?);
    let matched = material.check_original(
        source,
        placement,
        candidate.face_entities(),
        transfer,
        interrupt,
    )?;
    w.charge(matched.work_used())?;
    w.bytes(
        live.checked_add(matched.query_budget_bytes())
            .ok_or(QualificationError::Limit)?,
    )?;
    let original_volume = matched.original_source_volume_bounds_mm3();
    let healed_volume = material.healed_volume_bounds_mm3();
    let original = matched.source_to_raw_correspondence();
    let healed = material.raw_to_healed_correspondence();
    let composed_position = matched.source_to_healed_position_error_mm();
    let composed_derivative = matched.source_to_healed_normalized_derivative_error_mm();
    w.charge(0)?;
    let encoding_work = candidate.encoding_work();
    let text = candidate.text;
    Ok(QualifiedOrthogonalSeatsStep {
        text,
        source: snapshot,
        placement,
        quality: options.material,
        original_volume,
        healed_volume,
        original,
        healed,
        composed_position,
        composed_derivative,
        work: w.used,
        encoding_work,
        peak: w.peak,
    })
}
