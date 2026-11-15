//! Fresh actual-file qualification shared by emitted replay and Model readback.
use super::*;
/// Fresh actual-file query. Its aggregate128MiB live scope is separate from
/// replay Options; each selected geometry/material/binding child keeps64MiB.
#[derive(Clone, Copy)]
pub struct ReadbackOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    /// Other live caller geometry/cache allocations, excluding source and
    /// occurrence which this query already charges once.
    pub additional_resident_bytes: usize,
    pub position_error_mm: f64,
    pub original_derivative_error_mm: f64,
    pub normal_chord_error: f64,
}
impl Default for ReadbackOptions {
    fn default() -> Self {
        Self {
            maximum_work: WORK,
            maximum_bytes: 128 * 1024 * 1024,
            additional_resident_bytes: 0,
            position_error_mm: 1e-7,
            original_derivative_error_mm: 1e-6,
            normal_chord_error: 1e-6,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct ReadbackReceipt {
    proof: Receipt,
    bounds: [[f64; 2]; 3],
    count: usize,
    index: usize,
}
impl ReadbackReceipt {
    pub fn receipt(self) -> Receipt {
        self.proof
    }
    pub fn represented_member_count(self) -> usize {
        self.count
    }
    pub fn member_index(self) -> usize {
        self.index
    }
    pub fn actual_volume_mm3(self) -> [f64; 2] {
        self.proof.literal_volume_mm3()
    }
    pub fn actual_bounds_mm(self) -> [[f64; 2]; 3] {
        self.bounds
    }
}
/// Independently qualify one explicitly ordered represented member of ACTUAL
/// output bytes. This neither generates expected STEP nor validates other
/// compound members; callers must separately enforce count/order/ownership.
pub fn verify_readback(
    source: &OwnedMaterial,
    occurrence: Option<&read::LiteralOccurrencePlacement>,
    actual_text: &str,
    member_index: usize,
    options: ReadbackOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<ReadbackReceipt> {
    let quality = Options {
        maximum_work: options.maximum_work,
        maximum_bytes: BYTES,
        position_error_mm: options.position_error_mm,
        original_derivative_error_mm: options.original_derivative_error_mm,
        normal_chord_error: options.normal_chord_error,
    };
    quality.check()?;
    if options.maximum_bytes == 0 || options.maximum_bytes > 128 * 1024 * 1024 {
        return Err(Error::Invalid);
    }
    let reserved = source
        .retained_budget_bytes()
        .checked_add(occurrence.map_or(0, |p| p.retained_budget_bytes()))
        .and_then(|n| n.checked_add(options.additional_resident_bytes))
        .and_then(|n| n.checked_add(FIXED))
        .ok_or(Error::Limit)?;
    let mut work = Ledger {
        used: 0,
        maximum: options.maximum_work,
        peak: 0,
        bytes: options.maximum_bytes,
        interrupt,
    };
    work.step(1)?;
    work.observe(reserved)?;
    let document = CompoundDocument::parse(
        actual_text,
        DocumentOptions {
            maximum_work: work.available(512_000_000)?,
            retained_budget_bytes: room(options.maximum_bytes, reserved)?,
        },
        interrupt,
    )
    .map_err(document_error)?;
    let document_work = document.work_used() + 1;
    work.step(document.work_used())?;
    let document_peak = reserved + document.retained_budget_bytes();
    work.observe(document_peak)?;
    let count = document.represented_members().len();
    let q = qualify(
        source,
        occurrence,
        &document,
        SelectionScope {
            member_index,
            document_peak,
        },
        Options {
            maximum_bytes: options.maximum_bytes,
            ..quality
        },
        &mut work,
        interrupt,
    )?;
    Ok(ReadbackReceipt {
        proof: Receipt {
            inherited_source_work: source.work_used(),
            source_document_work: source.original_document_work_used(),
            graph_preparation_work: 0,
            encoding_work: 0,
            document_work,
            material_work: q.material_work,
            correspondence_work: q.correspondence_work,
            binding_work: q.binding_work,
            append_work: 0,
            work: work.used,
            peak: work.peak,
            matched_entities: q.matched_entities,
            position: q.position,
            derivative: q.derivative,
            normal: q.normal,
            literal_volume: q.volume,
        },
        bounds: q.bounds,
        count,
        index: member_index,
    })
}
pub(super) struct SelectionScope {
    pub(super) member_index: usize,
    pub(super) document_peak: usize,
}
pub(super) struct Qualified {
    pub(super) material_work: usize,
    pub(super) correspondence_work: usize,
    pub(super) binding_work: usize,
    pub(super) matched_entities: usize,
    pub(super) position: f64,
    pub(super) derivative: f64,
    pub(super) normal: f64,
    pub(super) volume: [f64; 2],
    pub(super) bounds: [[f64; 2]; 3],
}
pub(super) fn qualify(
    source: &OwnedMaterial,
    occurrence: Option<&read::LiteralOccurrencePlacement>,
    document: &CompoundDocument<'_>,
    scope: SelectionScope,
    options: Options,
    work: &mut Ledger<'_>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Qualified> {
    let SelectionScope {
        member_index,
        document_peak,
    } = scope;
    let members = document.represented_members();
    let member = *members.get(member_index).ok_or(Error::Invalid)?;
    let selected = Selected::parse(
        document,
        member.item_entity(),
        read::Options {
            maximum_work: work.available(16_000_000)?,
            maximum_bytes: room(options.maximum_bytes, document_peak)?.min(BYTES),
        },
        interrupt,
    )?;
    let pose = selected.placement(
        read::Options {
            maximum_work: work.available(16_000_000)?,
            maximum_bytes: room(options.maximum_bytes, document_peak)?.min(BYTES),
        },
        interrupt,
    )?;
    let pose_work = pose.work_used() - selected.parsed().work_used();
    let remaining = room(work.available(64_000_000)?, pose_work)?;
    let material_options = MaterialOptions {
        maximum_work: remaining,
        maximum_bytes: room(options.maximum_bytes, document_peak + size_of_val(&pose))?.min(BYTES),
        ..Default::default()
    };
    let occupied = selected
        .parsed()
        .projection(material_options, interrupt)?
        .occupied(material_options, interrupt)?;
    let d = pose.determinant_bounds();
    let literal_volume =
        I::new(occupied.volume_mm3()[0], occupied.volume_mm3()[1]) * I::new(d[0], d[1]);
    if d[0] <= 0. || !literal_volume.valid() || literal_volume.lo <= 0. {
        return Err(Error::Unsupported);
    }
    let material_work = occupied.work_used() + pose_work;
    work.step(material_work)?;
    let correspondence_start = work.used;
    let (position, derivative, normal) = if let Some(placement) = occurrence {
        placed::correspond(source, placement, &pose, options, work)?
    } else {
        (0., 0., 0.)
    };
    let correspondence_work = work.used - correspondence_start;
    work.observe(document_peak + occupied.query_budget_bytes() + size_of_val(&pose))?;
    assert!(
        occupied.volume_mm3()[0] > 0.,
        "earned independent occupied receipt"
    );
    let binding_start = work.used;
    let mut local = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
    for tensor in selected.parsed().tensors() {
        for index in tensor.point_indices() {
            work.step(32)?;
            let point = selected.parsed().points()[*index].coordinates();
            for axis in 0..3 {
                let value = point[axis].bounds();
                local[axis][0] = local[axis][0].min(value[0]);
                local[axis][1] = local[axis][1].max(value[1]);
            }
        }
    }
    let bounds = pose.map_box_mm(local)?;
    let mut ids = Vec::new();
    work.observe(
        document_peak + selected.parsed().query_budget_bytes() + selected.entity_ids().len() * 4,
    )?;
    ids.try_reserve_exact(selected.entity_ids().len())
        .map_err(|_| Error::Limit)?;
    ids.extend_from_slice(selected.entity_ids());
    work.step(ids.len())?;
    work.observe(document_peak + selected.parsed().query_budget_bytes() + ids.capacity() * 4)?;
    drop(selected);
    let binding = binding::bind_owned_document(
        source,
        if occurrence.is_some() {
            source.member().local_solid_entity()
        } else {
            source.member().item_entity()
        },
        document,
        &ids,
        if occurrence.is_some() {
            member.local_solid_entity()
        } else {
            member.item_entity()
        },
        binding::Options {
            maximum_work: work.available(512_000_000)?,
            maximum_bytes: room(options.maximum_bytes, document_peak + ids.capacity() * 4)?
                .min(BYTES),
            maximum_entities: 1_048_576,
        },
        interrupt,
    )
    .map_err(bind_error)?;
    work.step(binding.work_used())?;
    work.observe(document_peak + ids.capacity() * 4 + binding.query_peak_bytes())?;
    Ok(Qualified {
        material_work,
        correspondence_work,
        binding_work: work.used - binding_start,
        matched_entities: binding.matched_entities(),
        position,
        derivative,
        normal,
        volume: [literal_volume.lo, literal_volume.hi],
        bounds,
    })
}
