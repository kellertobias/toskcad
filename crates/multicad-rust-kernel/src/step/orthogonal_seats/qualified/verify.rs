//! Selected ACTUAL-file binding plus independent healed material readback.
use super::*;
use crate::{step_read::literal_member_binding as binding, tapered_thread::interval::I};
#[derive(Clone, Copy)]
pub struct VerificationOptions {
    /// Additional whole verification work. The immutable qualification and
    /// actual shared document are separately reported prerequisites.
    pub maximum_work: usize,
    /// Additional query allocations, excluding the actual shared document.
    pub maximum_bytes: usize,
    pub material: read::MaterialOptions,
}
impl Default for VerificationOptions {
    fn default() -> Self {
        Self {
            maximum_work: 64_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            material: read::MaterialOptions::default(),
        }
    }
}
pub struct VerifiedMember {
    item: u32,
    solid: u32,
    matched_entities: usize,
    original_volume: [f64; 2],
    healed_volume: [f64; 2],
    original: read::Correspondence,
    healed: read::Correspondence,
    position: f64,
    derivative: f64,
    work: usize,
    known: usize,
    peak: usize,
    document_bytes: usize,
}
impl VerifiedMember {
    pub fn item_entity(&self) -> u32 {
        self.item
    }
    pub fn solid_entity(&self) -> u32 {
        self.solid
    }
    pub fn matched_entities(&self) -> usize {
        self.matched_entities
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
        self.position
    }
    pub fn source_to_healed_normalized_derivative_error_mm(&self) -> f64 {
        self.derivative
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn aggregate_known_work(&self) -> usize {
        self.known
    }
    pub fn query_peak_accounted_bytes(&self) -> usize {
        self.peak
    }
    pub fn actual_document_retained_bytes(&self) -> usize {
        self.document_bytes
    }
}
fn binding_error(e: binding::Error) -> QualificationError {
    match e {
        binding::Error::Invalid => QualificationError::Invalid,
        binding::Error::Different => QualificationError::Unsupported,
        binding::Error::Cancelled => QualificationError::Cancelled,
        binding::Error::Limit
        | binding::Error::WorkLimit
        | binding::Error::StorageLimit
        | binding::Error::EntityLimit => QualificationError::Limit,
    }
}
impl QualifiedOrthogonalSeatsStep {
    /// The full reachable local graph and EXACT numeric/string tokens must
    /// match the sealed producer, permitting entity-ID rebasing only. Actual
    /// selected material is then independently reconstructed and measured.
    /// The receipt remains healed support-derived mass, never raw STEP mass.
    pub fn verify_actual_member(
        &self,
        actual: &CompoundDocument<'_>,
        item: u32,
        options: VerificationOptions,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<VerifiedMember, QualificationError> {
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
        let held =
            self.text.capacity() + self.source.capacity() * std::mem::size_of::<u64>() + 64 * 1024;
        w.bytes(held)?;
        w.bytes(
            held.checked_add(
                self.text
                    .len()
                    .checked_mul(128)
                    .ok_or(QualificationError::Limit)?,
            )
            .ok_or(QualificationError::Limit)?,
        )?;
        let expected = CompoundDocument::parse(
            self.as_str(),
            DocumentOptions {
                maximum_work: w.remaining_work(512_000_000)?,
                retained_budget_bytes: w
                    .remaining_bytes(held)?
                    .checked_add(self.text.len())
                    .ok_or(QualificationError::Limit)?,
            },
            Some(interrupt),
        )
        .map_err(document_error)?;
        w.charge(expected.work_used())?;
        let doc_extra = expected
            .retained_budget_bytes()
            .checked_sub(self.text.len())
            .ok_or(QualificationError::Invalid)?;
        let live = held
            .checked_add(doc_extra)
            .ok_or(QualificationError::Limit)?;
        w.bytes(live)?;
        let [expected_member] = expected.represented_members() else {
            return Err(QualificationError::Invalid);
        };
        let bound = binding::bind_local_members(
            &expected,
            expected_member.item_entity(),
            actual,
            item,
            binding::Options {
                maximum_work: w.remaining_work(16_000_000)?,
                maximum_bytes: w.remaining_bytes(live)?,
                maximum_entities: 4096,
            },
            Some(interrupt),
        )
        .map_err(binding_error)?;
        w.charge(bound.work_used())?;
        w.bytes(
            live.checked_add(bound.query_peak_bytes())
                .ok_or(QualificationError::Limit)?,
        )?;
        let boundary = read::decode(
            actual,
            item,
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
        w.bytes(
            live.checked_add(material.query_budget_bytes())
                .ok_or(QualificationError::Limit)?,
        )?;
        let healed = material.raw_to_healed_correspondence();
        let position = (I::point(self.original.maximum_position_error_mm())
            + I::point(healed.maximum_position_error_mm()))
        .hi;
        let derivative = (I::point(self.original.maximum_normalized_derivative_error_mm())
            + I::point(healed.maximum_normalized_derivative_error_mm()))
        .hi;
        if position > options.material.transfer.maximum_position_error_mm
            || derivative
                > options
                    .material
                    .transfer
                    .maximum_normalized_derivative_error_mm
            || self.original.maximum_sewing_error_mm()
                > options.material.transfer.maximum_sewing_error_mm
        {
            return Err(QualificationError::Unsupported);
        }
        let known = self
            .work
            .checked_add(actual.work_used())
            .and_then(|n| n.checked_add(w.used))
            .ok_or(QualificationError::Limit)?;
        w.charge(0)?;
        Ok(VerifiedMember {
            item,
            solid: boundary.solid_entity(),
            matched_entities: bound.matched_entities(),
            original_volume: self.original_volume,
            healed_volume: material.healed_volume_bounds_mm3(),
            original: self.original,
            healed,
            position,
            derivative,
            work: w.used,
            known,
            peak: w.peak,
            document_bytes: actual.retained_budget_bytes(),
        })
    }
}
