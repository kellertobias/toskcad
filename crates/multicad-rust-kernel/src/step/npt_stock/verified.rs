//! Explicit bounded encode/independent-readback/source-correspondence pipeline.
use super::writer::OccurrenceMap;
use super::*;
use crate::step_read::{
    npt_stock::{MaterialOptions, OwnedMaterial, Selected},
    owned_insert::{CompoundDocument, DocumentOptions},
};
const PIPELINE_SCRATCH: usize = 128 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct VerifiedOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub position_error_mm: f64,
    pub original_derivative_error_mm: f64,
}
impl Default for VerifiedOptions {
    fn default() -> Self {
        Self {
            maximum_work: 96_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            position_error_mm: 1e-7,
            original_derivative_error_mm: 1e-6,
        }
    }
}
impl VerifiedOptions {
    fn validate(self) -> Result<(), Error> {
        if self.maximum_work == 0
            || self.maximum_work > 96_000_000
            || self.maximum_bytes == 0
            || self.maximum_bytes > 64 * 1024 * 1024
            || !self.position_error_mm.is_finite()
            || self.position_error_mm <= 0.
            || !self.original_derivative_error_mm.is_finite()
            || self.original_derivative_error_mm <= 0.
        {
            return Err(Error::Invalid);
        }
        Ok(())
    }
}
#[derive(Clone, Copy, Debug)]
pub struct VerifiedReceipt {
    encoding_work: usize,
    document_work: usize,
    material_work: usize,
    correspondence_work: usize,
    orchestration_work: usize,
    aggregate_work: usize,
    encoding_peak: usize,
    document_peak: usize,
    material_peak: usize,
    correspondence_peak: usize,
    peak: usize,
    literal_bytes: usize,
    append_work: usize,
    correspondence: CorrespondenceReceipt,
}
/// The same accounting fields after a successful atomic verified replay.
pub type VerifiedAppendReceipt = VerifiedReceipt;
impl VerifiedReceipt {
    pub fn encoding_work_used(&self) -> usize {
        self.encoding_work
    }
    pub fn document_work_used(&self) -> usize {
        self.document_work
    }
    pub fn material_work_used(&self) -> usize {
        self.material_work
    }
    pub fn correspondence_work_used(&self) -> usize {
        self.correspondence_work
    }
    pub fn orchestration_work_used(&self) -> usize {
        self.orchestration_work
    }
    pub fn aggregate_work_used(&self) -> usize {
        self.aggregate_work
    }
    pub fn encoding_peak_bytes(&self) -> usize {
        self.encoding_peak
    }
    pub fn document_peak_bytes(&self) -> usize {
        self.document_peak
    }
    pub fn material_peak_bytes(&self) -> usize {
        self.material_peak
    }
    pub fn correspondence_peak_bytes(&self) -> usize {
        self.correspondence_peak
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
    pub fn literal_retained_bytes(&self) -> usize {
        self.literal_bytes
    }
    pub fn append_work_used(&self) -> usize {
        self.append_work
    }
    pub fn position_error_mm(&self) -> f64 {
        self.correspondence.position_error_mm()
    }
    pub fn original_derivative_error_mm(&self) -> [f64; 2] {
        self.correspondence.original_derivative_error_mm()
    }
    pub fn edge_derivative_error_mm(&self) -> f64 {
        self.correspondence.edge_derivative_error_mm()
    }
    pub fn local_correspondence(&self) -> CorrespondenceReceipt {
        self.correspondence
    }
}
pub struct VerifiedCandidate<'plan, 'source> {
    plan: &'plan SurfacePlan<'source>,
    material: OwnedMaterial,
    receipt: VerifiedReceipt,
}
impl<'plan, 'source> VerifiedCandidate<'plan, 'source> {
    pub fn source_plan(&self) -> &'plan SurfacePlan<'source> {
        self.plan
    }
    pub fn literal_material(&self) -> &OwnedMaterial {
        &self.material
    }
    pub fn receipt(&self) -> VerifiedReceipt {
        self.receipt
    }
}
struct Pipeline<'a> {
    used: usize,
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Pipeline<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        self.used += n;
        Ok(())
    }
    fn available(&self, cap: usize) -> Result<usize, Error> {
        if self.remaining == 0 {
            return Err(Error::Limit);
        }
        Ok(self.remaining.min(cap))
    }
}
// These private seams consume ONLY a fresh sealed typed-writer packet. Its
// fixed millimetre representation, dimensions2/3, unique IDs, degree11 uniform
// knots,20edges/12vertices/eight tensor supports, exact UV rows and cap bindings
// establish the decoder schema/topology. Encoder limits bound Part21 entity,
// value, token and nesting limits; all phase options are positive remaining
// validated balances. Structural Invalid is an internal producer invariant;
// independent public malformed-file readers retain their complete error APIs.
fn decoder_error(error: crate::step_read::npt_stock::Error) -> Error {
    assert!(
        !matches!(error, crate::step_read::npt_stock::Error::Invalid),
        "sealed typed NPT packet schema/topology"
    );
    if matches!(error, crate::step_read::npt_stock::Error::Unsupported) {
        Error::Unresolved
    } else if matches!(error, crate::step_read::npt_stock::Error::Limit) {
        Error::Limit
    } else {
        Error::Cancelled
    }
}
fn document_error(error: crate::step_read::owned_insert::Error) -> Error {
    // The neutral representation is fixed and valid independently of occupied
    // geometry, whose genuine Unsupported remains handled by decoder_error.
    if matches!(error, crate::step_read::owned_insert::Error::Limit) {
        Error::Limit
    } else {
        assert_eq!(
            error,
            crate::step_read::owned_insert::Error::Cancelled,
            "sealed typed NPT neutral document"
        );
        Error::Cancelled
    }
}
fn room(maximum: usize, retained: usize) -> Result<usize, Error> {
    maximum
        .checked_sub(retained)
        .filter(|n| *n > 0)
        .ok_or(Error::Limit)
}
impl<'source> SurfacePlan<'source> {
    /// Verify ACTUAL emitted numeric leaves and occupied geometry in a fresh
    /// independent reader. This local receipt never proves Native-world pose.
    pub fn verify_candidate<'plan>(
        &'plan self,
        options: VerifiedOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<VerifiedCandidate<'plan, 'source>, Error> {
        options.validate()?;
        let mut pipeline = Pipeline {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        pipeline.charge(2048)?;
        let phase_bytes = room(options.maximum_bytes, PIPELINE_SCRATCH)?;
        let file = write_source_candidate(
            self,
            pipeline.available(16_000_000)?,
            phase_bytes,
            interrupt,
        )?;
        let encoded = file.receipt();
        pipeline.charge(encoded.work_used())?;
        let encoding_peak = encoded.accounted_peak_bytes() + PIPELINE_SCRATCH;
        let (material, document_work, document_peak, material_peak) = {
            let text = file.into_string();
            let slack = text.capacity() - text.len();
            let reserved = self.retained_budget_bytes() + slack + PIPELINE_SCRATCH;
            let document = CompoundDocument::parse(
                &text,
                DocumentOptions {
                    maximum_work: pipeline.available(512_000_000)?,
                    retained_budget_bytes: room(options.maximum_bytes, reserved)?,
                },
                interrupt,
            )
            .map_err(document_error)?;
            let document_work = document.work_used();
            pipeline.charge(document_work)?;
            let document_peak = reserved + document.retained_budget_bytes();
            let material_bytes = room(options.maximum_bytes, document_peak)?;
            assert_eq!(
                document.represented_members().len(),
                1,
                "one typed candidate append"
            );
            let selected = Selected::parse(
                &document,
                document.represented_members()[0].item_entity(),
                crate::step_read::npt_stock::Options {
                    maximum_work: pipeline.available(16_000_000)?,
                    maximum_bytes: material_bytes,
                },
                interrupt,
            )
            .map_err(decoder_error)?;
            let material = selected
                .retain(
                    MaterialOptions {
                        maximum_work: pipeline.available(64_000_000)?,
                        maximum_bytes: material_bytes,
                        ..Default::default()
                    },
                    interrupt,
                )
                .map_err(decoder_error)?;
            pipeline.charge(material.work_used())?;
            let material_peak = document_peak + material.construction_query_peak_bytes();
            (material, document_work, document_peak, material_peak)
        };
        // Reader input/arena and temporary packet have been destroyed. The
        // following proof reads only immutable retained ORIGINAL components.
        let inherited = self.work_used() + material.work_used();
        let correspondence = self.correspond_literal(
            &material,
            CorrespondenceOptions {
                position_error_mm: options.position_error_mm,
                original_derivative_error_mm: options.original_derivative_error_mm,
                maximum_work: (inherited + pipeline.remaining).min(96_000_000),
                maximum_bytes: phase_bytes,
            },
            interrupt,
        )?;
        pipeline.charge(correspondence.query_work_used())?;
        let correspondence_peak = correspondence.accounted_peak_bytes() + PIPELINE_SCRATCH;
        let peak = encoding_peak
            .max(document_peak)
            .max(material_peak)
            .max(correspondence_peak);
        let receipt = VerifiedReceipt {
            encoding_work: encoded.work_used(),
            document_work,
            material_work: material.work_used(),
            correspondence_work: correspondence.query_work_used(),
            orchestration_work: 2048,
            aggregate_work: pipeline.used,
            encoding_peak,
            document_peak,
            material_peak,
            correspondence_peak,
            peak,
            literal_bytes: material.retained_budget_bytes(),
            append_work: 0,
            correspondence,
        };
        Ok(VerifiedCandidate {
            plan: self,
            material,
            receipt,
        })
    }
}
impl super::super::Writer {
    /// Atomic replay of the SAME immutable plan and deterministic numeric
    /// encoder that earned independent readback. IDs may be rebased; no numeric
    /// leaves, patch/edge owners, knots or orientation decisions are changed.
    /// This proves that appended member, not unrelated existing compound data.
    pub fn append_verified_source_npt_candidate(
        &mut self,
        candidate: &VerifiedCandidate<'_, '_>,
        options: VerifiedOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<VerifiedAppendReceipt, Error> {
        self.append_verified_npt_with(candidate, options, None, interrupt)
    }
    /// The same atomic replay of the verified LOCAL chart inside one mapped
    /// item carrying the occurrence's normalized operator. The local proof
    /// never certifies that operator: a consumer must decode the written map
    /// and bound its difference from the occurrence before admitting it.
    pub fn append_placed_verified_source_npt_candidate(
        &mut self,
        candidate: &VerifiedCandidate<'_, '_>,
        options: VerifiedOptions,
        placement: OccurrenceMap,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<VerifiedAppendReceipt, Error> {
        self.append_verified_npt_with(candidate, options, Some(placement), interrupt)
    }
    fn append_verified_npt_with(
        &mut self,
        candidate: &VerifiedCandidate<'_, '_>,
        options: VerifiedOptions,
        placement: Option<OccurrenceMap>,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<VerifiedAppendReceipt, Error> {
        options.validate()?;
        let mut pipeline = Pipeline {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        pipeline.charge(candidate.receipt.aggregate_work)?;
        pipeline.charge(256)?;
        if candidate.receipt.peak > options.maximum_bytes {
            return Err(Error::Limit);
        }
        if candidate.receipt.position_error_mm() > options.position_error_mm
            || candidate
                .receipt
                .original_derivative_error_mm()
                .iter()
                .any(|e| *e > options.original_derivative_error_mm)
            || candidate.receipt.edge_derivative_error_mm() > options.original_derivative_error_mm
        {
            return Err(Error::Unresolved);
        }
        let maximum_work = pipeline.available(16_000_000)?;
        // The inherited verified construction peak already includes these
        // immutable literal bytes and pipeline scratch; the earlier peak gate
        // dominates this subtraction. Actual writer growth remains fallible.
        let maximum_bytes =
            options.maximum_bytes - candidate.material.retained_budget_bytes() - PIPELINE_SCRATCH;
        let encoded = writer::append_verified_replay(
            self,
            candidate,
            placement,
            maximum_work,
            maximum_bytes,
            interrupt,
        )?;
        let additional = encoded.work_used();
        // Encoder was bounded by the remaining aggregate BEFORE its atomic
        // commit. Accounting after commit cannot introduce a new cancellation
        // or error that would leave a failed operation appended to the writer.
        assert!(
            additional <= pipeline.remaining,
            "sealed replay work allowance"
        );
        pipeline.remaining -= additional;
        pipeline.used += additional;
        let mut receipt = candidate.receipt;
        receipt.encoding_work += additional;
        receipt.append_work += additional;
        receipt.orchestration_work += 256;
        receipt.aggregate_work = pipeline.used;
        receipt.encoding_peak = receipt.encoding_peak.max(
            encoded.accounted_peak_bytes()
                + candidate.material.retained_budget_bytes()
                + PIPELINE_SCRATCH,
        );
        receipt.peak = receipt.peak.max(receipt.encoding_peak);
        Ok(receipt)
    }
}
#[cfg(test)]
mod tests;
