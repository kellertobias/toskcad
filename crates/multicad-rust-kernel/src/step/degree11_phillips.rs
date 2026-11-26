//! Original retained singleton proposal, independently decoded repaired material.
//! No nominal reconstruction, periodic cast or exact raw-material claim. Extra
//! occurrence placement/compound insertion requires a separate selected-member
//! mapping proof and remains outside this flat-singleton exchange contract.
use crate::step_read::{
    self, PhillipsDegree11ReadOptions as ReaderOptions,
    PhillipsDegree11ReplayMatchOptions as MatchOptions, PhillipsOwnedDegree11Material as Material,
    PhillipsReadError,
};
const FILE: usize = 64 * 1024 * 1024;
const KNOWN_BYTES: usize = 256 * 1024 * 1024;
const PROPOSAL: usize = 16_000_000;
const READER: usize = 624_016_384;
const MATCH: usize = 512_000_000;
const TOTAL: usize = PROPOSAL + READER + MATCH;
const FIXED: usize = 16 * 1024;
const IDENTITY: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
#[derive(Clone, Copy)]
pub struct Options {
    /// Conservative aggregate: original source prefix, encoding, independent
    /// reader reservation and FRESH match work. Reader's owning proof prefix
    /// is not counted again through match.total_work_used().
    pub maximum_aggregate_work: usize,
    pub maximum_proposal_work: usize,
    pub maximum_output_bytes: usize,
    /// Source/candidate/result/receipt numerical residency. Parser/reader child
    /// arenas retain their separately enforced caps; this is not process RSS
    /// or a claim that the complete reader fits one64MiB allocation envelope.
    pub maximum_known_resident_bytes: usize,
    pub reader: ReaderOptions,
    pub correspondence: MatchOptions,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_aggregate_work: TOTAL,
            maximum_proposal_work: PROPOSAL,
            maximum_output_bytes: FILE,
            maximum_known_resident_bytes: KNOWN_BYTES,
            reader: Default::default(),
            correspondence: Default::default(),
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Limit,
    Cancelled,
    Unsupported,
    Reader(step_read::PhillipsDegree11ReadError),
}
fn source_error(e: PhillipsReadError) -> Error {
    match e {
        PhillipsReadError::Invalid => Error::Invalid,
        PhillipsReadError::Limit => Error::Limit,
        PhillipsReadError::Cancelled => Error::Cancelled,
        PhillipsReadError::Unsupported => Error::Unsupported,
    }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn valid(o: Options) -> bool {
    o.maximum_aggregate_work > 0
        && o.maximum_aggregate_work <= TOTAL
        && o.maximum_proposal_work > 0
        && o.maximum_proposal_work <= PROPOSAL
        && o.maximum_output_bytes > 0
        && o.maximum_output_bytes <= FILE
        && o.maximum_known_resident_bytes > 0
        && o.maximum_known_resident_bytes <= KNOWN_BYTES
        && o.correspondence.maximum_work <= MATCH
        && o.correspondence.maximum_bytes <= FILE
        && [
            o.correspondence.maximum_distance_mm,
            o.correspondence.maximum_derivative_mm,
        ]
        .into_iter()
        .all(|x| x.is_finite() && x > 0. && x <= 1e-7)
        && o.correspondence.maximum_volume_discrepancy_mm3.is_finite()
        && o.correspondence.maximum_volume_discrepancy_mm3 > 0.
        && o.correspondence.maximum_volume_discrepancy_mm3 <= 1e-5
        && o.reader.maximum_reserved_work <= READER
        && o.reader.maximum_input_bytes <= FILE
        && o.reader
            .stage_work_limits
            .into_iter()
            .zip(ReaderOptions::default().stage_work_limits)
            .all(|(actual, hard)| actual > 0 && actual <= hard)
}
#[derive(Clone, Copy, Debug)]
pub struct Verification {
    pub maximum_world_distance_mm: f64,
    /// ORIGINAL normalized root-cell/source-atlas chart derivatives.
    pub maximum_world_chart_derivative_mm: f64,
    pub independent_repaired_volume_discrepancy_mm3: f64,
    pub decoded_repaired_volume_bounds_mm3: [f64; 2],
    pub original_source_prefix_work: usize,
    pub reader_reserved_work: usize,
    pub actual_decoded_material_proof_work: usize,
    pub fresh_correspondence_work: usize,
    pub proposal_work: usize,
    pub aggregate_reserved_work: usize,
    /// Logical source/candidate/result bytes plus bounded wrapper workspace.
    /// Borrowed actual-input spans are counted by length; unknown caller
    /// backing-buffer capacity and delegated child arenas are separate.
    pub known_resident_peak_bytes: usize,
}
pub struct QualifiedMaterial {
    material: Material,
    verification: Verification,
}
impl QualifiedMaterial {
    pub fn material(&self) -> &Material {
        &self.material
    }
    pub fn verification(&self) -> Verification {
        self.verification
    }
    pub fn into_material(self) -> Material {
        self.material
    }
}
pub struct QualifiedStep {
    text: String,
    qualification: QualifiedMaterial,
}
impl QualifiedStep {
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn qualification(&self) -> &QualifiedMaterial {
        &self.qualification
    }
}
const IDENTITY_OPERATOR: [[[f64; 2]; 3]; 4] = [
    [[0., 0.], [0., 0.], [0., 0.]],
    [[1., 1.], [0., 0.], [0., 0.]],
    [[0., 0.], [1., 1.], [0., 0.]],
    [[0., 0.], [0., 0.], [1., 1.]],
];
fn qualification(
    source: &Material,
    text: &str,
    candidate_capacity: usize,
    proposal_work: usize,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedMaterial, Error> {
    qualification_through(
        source,
        text,
        candidate_capacity,
        proposal_work,
        IDENTITY,
        IDENTITY_OPERATOR,
        options,
        interrupt,
    )
}
/// `text` is the decoded file's LOCAL member; `operator` encloses its file
/// map (identity for a flat singleton) and `occurrence` is the source's.
#[allow(clippy::too_many_arguments)]
fn qualification_through(
    source: &Material,
    text: &str,
    candidate_capacity: usize,
    proposal_work: usize,
    occurrence: [[f64; 3]; 4],
    operator: [[[f64; 2]; 3]; 4],
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedMaterial, Error> {
    checkpoint(interrupt)?;
    if !valid(options) {
        return Err(Error::Invalid);
    }
    if text.len() > options.maximum_output_bytes {
        return Err(Error::Limit);
    }
    let before = source
        .work_used()
        .checked_add(proposal_work)
        .ok_or(Error::Limit)?;
    // Reserve the complete independently delegated reader scope BEFORE reading.
    // It includes the eventual decoded owner proof; fresh matching below adds
    // only its genuinely new work, never the decoded owning prefix again.
    let reservation = options.reader.maximum_reserved_work;
    if reservation > READER {
        return Err(Error::Invalid);
    }
    let used = before
        .checked_add(reservation)
        .filter(|x| *x <= options.maximum_aggregate_work)
        .ok_or(Error::Limit)?;
    let resident = source
        .retained_budget_bytes()
        .checked_add(candidate_capacity)
        .and_then(|x| x.checked_add(FIXED + std::mem::size_of::<QualifiedStep>()))
        .filter(|x| *x <= options.maximum_known_resident_bytes)
        .ok_or(Error::Limit)?;
    let read = step_read::read_degree11_phillips_material(text, options.reader, interrupt)
        .map_err(Error::Reader)?;
    let reader_reserved_work = read.reserved_delegated_work();
    let actual_decoded_material_proof_work = read.material_proof_work_used();
    // The upfront caller reservation was conservative. The sealed reader now
    // reports its ACTUAL delegated reservation after any lowered stage caps.
    // Replace only that unused reservation; never remove actual child work.
    let used = before
        .checked_add(reader_reserved_work)
        .filter(|x| *x <= used)
        .ok_or(Error::Limit)?;
    let inherited = source
        .work_used()
        .checked_add(read.material().work_used())
        .ok_or(Error::Limit)?;
    let remaining = options
        .maximum_aggregate_work
        .checked_sub(used)
        .ok_or(Error::Limit)?;
    let maximum_work = options
        .correspondence
        .maximum_work
        .min(inherited.checked_add(remaining).ok_or(Error::Limit)?);
    let maximum_bytes = options.correspondence.maximum_bytes.min(
        options
            .maximum_known_resident_bytes
            .checked_sub(candidate_capacity + FIXED + std::mem::size_of::<QualifiedStep>())
            .ok_or(Error::Limit)?,
    );
    let matched = source
        .match_replayed_material_through(
            read.material(),
            occurrence,
            operator,
            MatchOptions {
                maximum_work,
                maximum_bytes,
                ..options.correspondence
            },
            interrupt,
        )
        .map_err(source_error)?;
    let fresh = matched.new_work_used();
    let aggregate = used
        .checked_add(fresh)
        .filter(|x| *x <= options.maximum_aggregate_work)
        .ok_or(Error::Limit)?;
    let peak = resident.max(
        matched
            .accounted_peak_bytes()
            .checked_add(candidate_capacity + FIXED + std::mem::size_of::<QualifiedStep>())
            .filter(|x| *x <= options.maximum_known_resident_bytes)
            .ok_or(Error::Limit)?,
    );
    let verification = Verification {
        maximum_world_distance_mm: matched.maximum_world_distance_mm(),
        maximum_world_chart_derivative_mm: matched.maximum_world_chart_derivative_mm(),
        independent_repaired_volume_discrepancy_mm3: matched.independent_volume_discrepancy_mm3(),
        decoded_repaired_volume_bounds_mm3: read.material().world_volume_mm3(),
        original_source_prefix_work: source.work_used(),
        reader_reserved_work,
        actual_decoded_material_proof_work,
        fresh_correspondence_work: fresh,
        proposal_work,
        aggregate_reserved_work: aggregate,
        known_resident_peak_bytes: peak,
    };
    checkpoint(interrupt)?;
    Ok(QualifiedMaterial {
        material: read.into_material(),
        verification,
    })
}
/// Fresh qualification of the ACTUAL singleton bytes. Valid but stale/other
/// geometry is refused by the complete independently earned material match.
pub fn qualify_actual_file(
    source: &Material,
    text: &str,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedMaterial, Error> {
    qualification(source, text, text.len(), 0, options, interrupt)
}
/// Preserve original retained tokens as an UNTRUSTED proposal. The accepted
/// interpretation is the independently decoded DISTINCTLY repaired material,
/// not exact physical mass/closure of the original literal STEP surfaces.
pub fn export_retained(
    source: &Material,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedStep, Error> {
    checkpoint(interrupt)?;
    if !valid(options) {
        return Err(Error::Invalid);
    }
    let input = source.original_input();
    if input.len() > options.maximum_output_bytes {
        return Err(Error::Limit);
    }
    // Two linear byte actions cover original reads and destination writes;
    // fixed metadata/checkpoint/try-reserve/publication are prepaid separately.
    let proposal = input
        .len()
        .checked_mul(2)
        .and_then(|x| x.checked_add(FIXED))
        .filter(|x| *x <= options.maximum_proposal_work)
        .ok_or(Error::Limit)?;
    let preflight = source
        .retained_budget_bytes()
        .checked_add(input.len().checked_mul(2).ok_or(Error::Limit)?)
        .and_then(|x| x.checked_add(FIXED + std::mem::size_of::<QualifiedStep>()))
        .filter(|x| *x <= options.maximum_known_resident_bytes)
        .ok_or(Error::Limit)?;
    let mut text = String::new();
    text.try_reserve_exact(input.len())
        .map_err(|_| Error::Limit)?;
    if source
        .retained_budget_bytes()
        .checked_add(text.capacity())
        .and_then(|x| x.checked_add(FIXED))
        .is_none_or(|x| x > options.maximum_known_resident_bytes)
    {
        return Err(Error::Limit);
    }
    text.push_str(input);
    checkpoint(interrupt)?;
    let mut qualification =
        qualification(source, &text, text.capacity(), proposal, options, interrupt)?;
    qualification.verification.known_resident_peak_bytes = qualification
        .verification
        .known_resident_peak_bytes
        .max(preflight);
    checkpoint(interrupt)?;
    Ok(QualifiedStep {
        text,
        qualification,
    })
}
mod placed;
pub use placed::{export_retained_placed, qualify_placed_actual_file};
mod authored;
pub use authored::{
    AuthoredAppendReceipt, AuthoredOptions, AuthoredVerification, QualifiedAuthoredMaterial,
    QualifiedAuthoredStep, export_authored, qualify_authored_file, supports_author,
};
