//! Exact retained-geometry STEP replay, independently decoded before append.
//! No authored thread recipe or display mesh participates in this path.
use super::retained_replay::{self, Node, OUTPUT, Trace, View, digits};
use super::*;
use crate::step_read::{
    literal_member_binding::{self as binding, RetainedRecords, Token, token},
    npt_stock::{self as read, MaterialOptions, OwnedMaterial, Selected},
    owned_insert::{CompoundDocument, DocumentOptions},
};
use crate::tapered_thread::interval::I;
pub use read::Error;
use std::mem::size_of;
mod placed;
mod readback;
pub use readback::{ReadbackOptions, ReadbackReceipt, verify_readback};
type Result<T> = std::result::Result<T, Error>;
const WORK: usize = 96_000_000;
const BYTES: usize = 64 * 1024 * 1024;
const FIXED: usize = 128 * 1024;
#[derive(Clone, Copy)]
pub struct Options {
    /// One NEW export query: encode, independent readback, binding and append.
    /// Completed input-document/source qualification is reported separately.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub position_error_mm: f64,
    pub original_derivative_error_mm: f64,
    pub normal_chord_error: f64,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: WORK,
            maximum_bytes: BYTES,
            position_error_mm: 1e-7,
            original_derivative_error_mm: 1e-6,
            normal_chord_error: 1e-6,
        }
    }
}
impl Options {
    fn check(self) -> Result<()> {
        if self.maximum_work == 0
            || self.maximum_work > WORK
            || self.maximum_bytes == 0
            || self.maximum_bytes > BYTES
            || [
                self.position_error_mm,
                self.original_derivative_error_mm,
                self.normal_chord_error,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0.)
        {
            return Err(Error::Invalid);
        }
        Ok(())
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Receipt {
    inherited_source_work: usize,
    source_document_work: usize,
    graph_preparation_work: usize,
    encoding_work: usize,
    document_work: usize,
    material_work: usize,
    correspondence_work: usize,
    binding_work: usize,
    append_work: usize,
    work: usize,
    peak: usize,
    matched_entities: usize,
    position: f64,
    derivative: f64,
    normal: f64,
    literal_volume: [f64; 2],
}
impl Receipt {
    /// Completed input import/qualification scope; NOT charged again to
    /// this new query. All actual encode/readback/binding phases share work.
    pub fn source_qualification_work_used(self) -> usize {
        self.inherited_source_work
    }
    pub fn source_document_work_used(self) -> usize {
        self.source_document_work
    }
    /// One input-document scope, one selected-source qualification and this
    /// export. Do not sum over members sharing the same original document.
    pub fn end_to_end_work_used(self) -> usize {
        self.source_document_work + self.inherited_source_work + self.work
    }
    pub fn graph_preparation_work_used(self) -> usize {
        self.graph_preparation_work
    }
    pub fn encoding_work_used(self) -> usize {
        self.encoding_work
    }
    pub fn document_work_used(self) -> usize {
        self.document_work
    }
    pub fn material_work_used(self) -> usize {
        self.material_work
    }
    pub fn correspondence_work_used(self) -> usize {
        self.correspondence_work
    }
    pub fn binding_work_used(self) -> usize {
        self.binding_work
    }
    pub fn append_work_used(self) -> usize {
        self.append_work
    }
    /// Aggregate of the NEW export phases, bounded by maximum_work.
    pub fn aggregate_work_used(self) -> usize {
        self.work
    }
    pub fn query_work_used(self) -> usize {
        self.work
    }
    pub fn accounted_peak_bytes(self) -> usize {
        self.peak
    }
    pub fn matched_entities(self) -> usize {
        self.matched_entities
    }
    /// Whole-source world correspondence. Exact replay binds the complete
    /// original operator; Native replay measures the independently decoded map.
    pub fn position_error_mm(self) -> f64 {
        self.position
    }
    pub fn original_derivative_error_mm(self) -> f64 {
        self.derivative
    }
    pub fn normal_chord_error(self) -> f64 {
        self.normal
    }
    pub fn literal_volume_mm3(self) -> [f64; 2] {
        self.literal_volume
    }
}
struct SourceView<'a>(&'a OwnedMaterial);
impl RetainedRecords for SourceView<'_> {
    fn record_count(&self) -> usize {
        self.0.record_count()
    }
    fn record_index(&self, id: u32) -> Option<usize> {
        self.0.record_index(id)
    }
    fn component(&self, record: usize, index: usize) -> Option<&str> {
        self.0.component(record, index)
    }
}
/// Sealed immutable actual packet and its verified record spans. Original
/// components are retained by the source; this owns only its actual encoding.
pub struct VerifiedReplay<'source> {
    source: &'source OwnedMaterial,
    occurrence: Option<&'source read::LiteralOccurrencePlacement>,
    packet: String,
    trace: Trace,
    receipt: Receipt,
}
impl VerifiedReplay<'_> {
    pub fn source(&self) -> &OwnedMaterial {
        self.source
    }
    pub fn occurrence(&self) -> Option<&read::LiteralOccurrencePlacement> {
        self.occurrence
    }
    pub fn receipt(&self) -> Receipt {
        self.receipt
    }
    pub fn standalone_step(&self) -> &str {
        &self.packet
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.source.retained_budget_bytes()
            + self.occurrence.map_or(0, |p| p.retained_budget_bytes())
            + self.packet.capacity()
            + self.trace.nodes.capacity() * size_of::<Node>()
            + size_of::<Self>()
    }
}
struct Ledger<'a> {
    used: usize,
    maximum: usize,
    peak: usize,
    bytes: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Ledger<'_> {
    fn step(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        if self.used > self.maximum {
            return Err(Error::Limit);
        }
        Ok(())
    }
    fn observe(&mut self, n: usize) -> Result<()> {
        if n > self.bytes {
            return Err(Error::Limit);
        }
        self.peak = self.peak.max(n);
        Ok(())
    }
    fn available(&self, cap: usize) -> Result<usize> {
        self.maximum
            .checked_sub(self.used)
            .filter(|n| *n > 0)
            .map(|n| n.min(cap))
            .ok_or(Error::Limit)
    }
}
impl retained_replay::Budget for Ledger<'_> {
    type Error = Error;
    fn step(&mut self, n: usize) -> Result<()> {
        self.step(n)
    }
    fn observe(&mut self, n: usize) -> Result<()> {
        self.observe(n)
    }
    fn invalid() -> Error {
        Error::Invalid
    }
    fn limit() -> Error {
        Error::Limit
    }
}
fn encode(
    view: &impl RetainedRecords,
    root: u32,
    file: &mut StepFile,
    base: usize,
    work: &mut Ledger<'_>,
) -> Result<Trace> {
    retained_replay::encode(view, root, file, base, work, 16_000_000)
}
fn room(total: usize, used: usize) -> Result<usize> {
    total
        .checked_sub(used)
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)
}
fn bind_error(e: binding::Error) -> Error {
    match e {
        binding::Error::Cancelled => Error::Cancelled,
        binding::Error::Invalid => Error::Invalid,
        binding::Error::Different => Error::Unsupported,
        _ => Error::Limit,
    }
}
fn document_error(e: crate::step_read::owned_insert::Error) -> Error {
    match e {
        crate::step_read::owned_insert::Error::Cancelled => Error::Cancelled,
        crate::step_read::owned_insert::Error::Invalid => Error::Invalid,
        crate::step_read::owned_insert::Error::Unsupported => Error::Unsupported,
        _ => Error::Limit,
    }
}
impl From<StepError> for Error {
    fn from(e: StepError) -> Self {
        match e {
            StepError::Limit => Self::Limit,
            _ => Self::Unsupported,
        }
    }
}
/// Verify the ACTUAL standalone output, including occupied local geometry and
/// its proper original mapping. Every numeric leaf and reference binding is
/// compared to retained originals after the neutral material phase.
pub fn verify_exact_replay<'source>(
    source: &'source OwnedMaterial,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<VerifiedReplay<'source>> {
    verify(source, None, options, interrupt)
}
pub fn verify_placed_replay<'source>(
    source: &'source OwnedMaterial,
    occurrence: &'source read::LiteralOccurrencePlacement,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<VerifiedReplay<'source>> {
    verify(source, Some(occurrence), options, interrupt)
}
fn verify<'source>(
    source: &'source OwnedMaterial,
    occurrence: Option<&'source read::LiteralOccurrencePlacement>,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<VerifiedReplay<'source>> {
    options.check()?;
    let mut work = Ledger {
        used: 0,
        maximum: options.maximum_work,
        peak: 0,
        bytes: options.maximum_bytes,
        interrupt,
    };
    let encoding_start = work.used;
    work.step(1)?;
    let mut writer = Writer::new();
    let mut trace = None;
    let mut graph_preparation_work = 0;
    writer.append_typed(|file| {
        let encoded = if let Some(placement) = occurrence {
            let (encoded, preparation) = placed::encode(source, placement, file, &mut work)?;
            graph_preparation_work = preparation;
            encoded
        } else {
            encode(
                &SourceView(source),
                source.member().item_entity(),
                file,
                source.retained_budget_bytes(),
                &mut work,
            )?
        };
        // The emitted local BRep stays literal. Only a separately decoded
        // output operator may approximate a NEW Native occurrence.
        let member = encoded.member;
        trace = Some(encoded);
        Ok::<_, Error>(Entity(member as usize))
    })?;
    let trace = trace.expect("successful standalone replay encoding");
    let packet = writer.finish().map_err(Error::from)?;
    let encoding_work = work.used - encoding_start - graph_preparation_work;
    let trace_bytes = trace.nodes.capacity() * size_of::<Node>();
    let reserved = source.retained_budget_bytes()
        + occurrence.map_or(0, |p| p.retained_budget_bytes())
        + trace_bytes
        + (packet.capacity() - packet.len())
        + FIXED;
    let document = CompoundDocument::parse(
        &packet,
        DocumentOptions {
            maximum_work: work.available(512_000_000)?,
            retained_budget_bytes: room(options.maximum_bytes, reserved)?,
        },
        interrupt,
    )
    .map_err(document_error)?;
    let document_work = document.work_used();
    work.step(document_work)?;
    let document_peak = reserved + document.retained_budget_bytes();
    work.observe(document_peak)?;
    assert_eq!(
        document.represented_members().len(),
        1,
        "one typed retained member encoding"
    );
    let q = readback::qualify(
        source,
        occurrence,
        &document,
        readback::SelectionScope {
            member_index: 0,
            document_peak,
        },
        options,
        &mut work,
        interrupt,
    )?;
    let receipt = Receipt {
        inherited_source_work: source.work_used(),
        source_document_work: source.original_document_work_used(),
        graph_preparation_work,
        encoding_work,
        document_work,
        material_work: q.material_work,
        binding_work: q.binding_work,
        correspondence_work: q.correspondence_work,
        append_work: 0,
        work: work.used,
        peak: work.peak,
        matched_entities: q.matched_entities,
        position: q.position,
        derivative: q.derivative,
        normal: q.normal,
        literal_volume: q.volume,
    };
    drop(document);
    Ok(VerifiedReplay {
        source,
        occurrence,
        packet,
        trace,
        receipt,
    })
}
impl Writer {
    /// Append only a sealed independently read-back packet. Actual rebasing is
    /// bound again to its verified token graph before the atomic append commits.
    pub fn append_verified_literal_npt(
        &mut self,
        candidate: &VerifiedReplay<'_>,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Receipt> {
        options.check()?;
        let mut work = Ledger {
            used: 0,
            maximum: options.maximum_work,
            peak: candidate.receipt.peak,
            bytes: options.maximum_bytes,
            interrupt,
        };
        work.step(candidate.receipt.work)?;
        work.observe(candidate.receipt.peak)?;
        let mut receipt = None;
        self.append_typed(|file| {
            work.observe(candidate.retained_budget_bytes() + file.data.capacity() + FIXED)?;
            let view = View {
                data: &candidate.packet,
                trace: &candidate.trace,
            };
            let before = work.used;
            let trace = encode(
                &view,
                candidate.trace.member,
                file,
                candidate.retained_budget_bytes(),
                &mut work,
            )?;
            let binding = binding::bind_retained(
                &view,
                candidate.trace.member,
                &View {
                    data: &file.data,
                    trace: &trace,
                },
                trace.member,
                binding::Options {
                    maximum_work: work.available(512_000_000)?,
                    maximum_bytes: room(
                        options.maximum_bytes,
                        candidate.retained_budget_bytes()
                            + file.data.capacity()
                            + trace.nodes.capacity() * size_of::<Node>()
                            + FIXED,
                    )?,
                    maximum_entities: 1_048_576,
                },
                interrupt,
            )
            .map_err(bind_error)?;
            work.step(binding.work_used())?;
            work.observe(
                candidate.retained_budget_bytes()
                    + file.data.capacity()
                    + trace.nodes.capacity() * size_of::<Node>()
                    + FIXED
                    + binding.query_peak_bytes(),
            )?;
            receipt = Some(Receipt {
                append_work: work.used - before,
                work: work.used,
                peak: work.peak,
                ..candidate.receipt
            });
            Ok::<_, Error>(Entity(trace.member as usize))
        })?;
        Ok(receipt.expect("successful atomic append receipt"))
    }
}
#[cfg(test)]
mod tests;
