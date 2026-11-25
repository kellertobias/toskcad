//! Authored Phillips source -> STEP proposal -> independently decoded repaired
//! material, matched back to the SAME author in its occurrence frame. The
//! accepted material is the decoded repaired volume; the writer's surfaces are
//! an untrusted proposal and no exact raw-file mass claim is made.
use super::*;
use crate::{
    phillips_drive::PhillipsScrew,
    step::{self, periodic_phillips},
    step_read::{
        PhillipsDegree11AuthorMatchReceipt as MatchReceipt,
        PhillipsSuppliedDegree11MaterialOptions as SuppliedOptions,
    },
};
const AUTHORED_TOTAL: usize = PROPOSAL + READER + MATCH;
#[derive(Clone, Copy)]
pub struct AuthoredOptions {
    /// Proposal encoding plus the reserved reader scope plus the complete
    /// author-match allowance, each charged once.
    pub maximum_aggregate_work: usize,
    pub maximum_proposal_work: usize,
    pub maximum_output_bytes: usize,
    pub reader: ReaderOptions,
    pub matching: SuppliedOptions,
}
impl Default for AuthoredOptions {
    fn default() -> Self {
        Self {
            maximum_aggregate_work: AUTHORED_TOTAL,
            maximum_proposal_work: PROPOSAL,
            maximum_output_bytes: FILE,
            reader: Default::default(),
            matching: Default::default(),
        }
    }
}
fn valid_authored(o: AuthoredOptions) -> bool {
    o.maximum_aggregate_work > 0
        && o.maximum_aggregate_work <= AUTHORED_TOTAL
        && o.maximum_proposal_work > 0
        && o.maximum_proposal_work <= PROPOSAL
        && o.maximum_output_bytes > 0
        && o.maximum_output_bytes <= FILE
        && o.matching.maximum_work > 0
        && o.matching.maximum_work <= MATCH
        && o.reader.maximum_reserved_work <= READER
        && o.reader.maximum_input_bytes <= FILE
        && o.reader
            .stage_work_limits
            .into_iter()
            .zip(ReaderOptions::default().stage_work_limits)
            .all(|(actual, hard)| actual > 0 && actual <= hard)
}
fn valid_frame([origin, u, v]: [[f64; 3]; 3]) -> bool {
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    origin.iter().chain(&u).chain(&v).all(|x| x.is_finite())
        && (dot(u, u) - 1.).abs() <= 1e-12
        && (dot(v, v) - 1.).abs() <= 1e-12
        && dot(u, v).abs() <= 1e-12
}
/// Nonadmitting routing check, never a geometry certificate: one start, a Pan
/// or Countersunk head, and a countersunk recess whose deepest tip leaves the
/// required wall above the head bottom (the periodic family's own rule).
/// Qualification of the actual written bytes remains the only admission.
pub fn supports_author(author: &PhillipsScrew) -> bool {
    use crate::head_profile::Head;
    let spec = author.spec();
    if author.thread().spec().starts != 1 {
        return false;
    }
    match spec.head {
        Some(Head::Crowned { .. }) => true,
        Some(Head::Countersunk { .. }) => {
            let tip = author.cavity().source().envelope().tip_depth_mm_bounds();
            let profile = author.profile();
            // Directed lower bound of (end - top) - tip.
            let floor = ((profile.end()[1] - profile.top()[1]).next_down() - tip[1]).next_down();
            spec.minimum_wall_mm.is_finite()
                && spec.minimum_wall_mm > 0.
                && floor.is_finite()
                && floor > spec.minimum_wall_mm
        }
        _ => false,
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AuthoredVerification {
    pub maximum_world_distance_mm: f64,
    /// Whole-source paired normalized chart derivatives, not unit normals.
    pub maximum_world_chart_derivative_mm: f64,
    /// Author minus independently decoded repaired material.
    pub independent_volume_difference_mm3: [f64; 2],
    pub decoded_repaired_volume_bounds_mm3: [f64; 2],
    pub reader_reserved_work: usize,
    pub actual_decoded_material_proof_work: usize,
    pub fresh_correspondence_work: usize,
    pub correspondence_cavity_cells: usize,
    pub correspondence_peak_bytes: usize,
    pub proposal_work: usize,
    pub aggregate_reserved_work: usize,
}
pub struct QualifiedAuthoredMaterial {
    material: Material,
    verification: AuthoredVerification,
}
impl QualifiedAuthoredMaterial {
    pub fn material(&self) -> &Material {
        &self.material
    }
    pub fn verification(&self) -> AuthoredVerification {
        self.verification
    }
}
pub struct QualifiedAuthoredStep {
    text: String,
    qualification: QualifiedAuthoredMaterial,
}
impl QualifiedAuthoredStep {
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn qualification(&self) -> &QualifiedAuthoredMaterial {
        &self.qualification
    }
}
fn qualify(
    author: &PhillipsScrew,
    frame: [[f64; 3]; 3],
    text: &str,
    proposal_work: usize,
    options: AuthoredOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedAuthoredMaterial, Error> {
    checkpoint(interrupt)?;
    if !valid_authored(options) || !valid_frame(frame) {
        return Err(Error::Invalid);
    }
    if text.len() > options.maximum_output_bytes {
        return Err(Error::Limit);
    }
    // Reserve the complete reader scope and the whole match allowance BEFORE
    // reading; the actual reader reservation then replaces its upper bound.
    let reserved = proposal_work
        .checked_add(options.reader.maximum_reserved_work)
        .and_then(|x| x.checked_add(options.matching.maximum_work))
        .filter(|x| *x <= options.maximum_aggregate_work)
        .ok_or(Error::Limit)?;
    let (read, matched): (_, MatchReceipt) = step_read::read_degree11_phillips_material_for_author(
        text,
        options.reader,
        author,
        frame,
        options.matching,
        interrupt,
    )
    .map_err(Error::Reader)?;
    let aggregate = proposal_work
        .checked_add(read.reserved_delegated_work())
        .and_then(|x| x.checked_add(matched.fresh_work))
        .filter(|x| *x <= reserved)
        .ok_or(Error::Limit)?;
    let verification = AuthoredVerification {
        maximum_world_distance_mm: matched.maximum_world_distance_mm,
        maximum_world_chart_derivative_mm: matched.maximum_world_chart_derivative_mm,
        independent_volume_difference_mm3: matched.independent_volume_difference_mm3,
        decoded_repaired_volume_bounds_mm3: read.material().world_volume_mm3(),
        reader_reserved_work: read.reserved_delegated_work(),
        actual_decoded_material_proof_work: read.material_proof_work_used(),
        fresh_correspondence_work: matched.fresh_work,
        correspondence_cavity_cells: matched.cavity_cells,
        correspondence_peak_bytes: matched.accounted_peak_bytes,
        proposal_work,
        aggregate_reserved_work: aggregate,
    };
    checkpoint(interrupt)?;
    Ok(QualifiedAuthoredMaterial {
        material: read.into_material(),
        verification,
    })
}
/// Fresh qualification of these ACTUAL bytes against the authored source in
/// `frame` = [origin, u, v]. Valid but different geometry is refused by the
/// complete author/decoded-material correspondence, never by a marker.
pub fn qualify_authored_file(
    author: &PhillipsScrew,
    frame: [[f64; 3]; 3],
    text: &str,
    options: AuthoredOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedAuthoredMaterial, Error> {
    qualify(author, frame, text, 0, options, interrupt)
}
/// Write the authored proposal in its world frame, then qualify exactly those
/// bytes. Nothing is returned unless the decoded material matches the author.
pub fn export_authored(
    author: &PhillipsScrew,
    frame: [[f64; 3]; 3],
    tolerance_mm: f64,
    options: AuthoredOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<QualifiedAuthoredStep, Error> {
    checkpoint(interrupt)?;
    if !valid_authored(options)
        || !valid_frame(frame)
        || !tolerance_mm.is_finite()
        || tolerance_mm <= 0.
        || tolerance_mm > options.matching.maximum_distance_mm
    {
        return Err(Error::Invalid);
    }
    let [origin, u, v] = frame;
    let text =
        step::write_phillips_screw_with_interrupt(author, origin, u, v, tolerance_mm, interrupt)
            .map_err(|e| match e {
                step::phillips_cutter::Error::Limit => Error::Limit,
                step::phillips_cutter::Error::Unsupported => Error::Unsupported,
                step::phillips_cutter::Error::Cancelled => Error::Cancelled,
                step::phillips_cutter::Error::InvalidGeometry
                | step::phillips_cutter::Error::InvalidPlacement => Error::Invalid,
            })?;
    if text.len() > options.maximum_output_bytes {
        return Err(Error::Limit);
    }
    // Two linear byte actions cover encoding and the actual-byte handoff;
    // fixed metadata is prepaid. Writer geometry stays under its own caps.
    let proposal = text
        .len()
        .checked_mul(2)
        .and_then(|x| x.checked_add(FIXED))
        .filter(|x| *x <= options.maximum_proposal_work)
        .ok_or(Error::Limit)?;
    let qualification = qualify(author, frame, &text, proposal, options, interrupt)?;
    checkpoint(interrupt)?;
    Ok(QualifiedAuthoredStep {
        text,
        qualification,
    })
}
#[derive(Clone, Copy, Debug)]
pub struct AuthoredAppendReceipt {
    pub matched_entities: usize,
    pub document_work: usize,
    pub selection_work: usize,
    pub encoding_work: usize,
    pub binding_work: usize,
    pub work_used: usize,
    pub accounted_peak_bytes: usize,
    pub candidate: AuthoredVerification,
}
impl step::Writer {
    /// Atomically append a sealed qualified authored candidate by exact
    /// rebased replay of its own tokens. Material authority is the candidate's
    /// completed qualification; the replay only binds the emitted graph.
    pub fn append_qualified_authored_degree11_phillips(
        &mut self,
        candidate: &QualifiedAuthoredStep,
        options: periodic_phillips::AppendOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<AuthoredAppendReceipt, Error> {
        let bytes = candidate.text.capacity()
            + candidate.qualification.material.retained_budget_bytes()
            + std::mem::size_of::<QualifiedAuthoredStep>();
        let core = self
            .append_qualified_singleton(candidate.as_str(), bytes, options, None, interrupt)
            .map_err(|e| match e {
                periodic_phillips::Error::Limit => Error::Limit,
                periodic_phillips::Error::Cancelled => Error::Cancelled,
                periodic_phillips::Error::Unsupported => Error::Unsupported,
                _ => Error::Invalid,
            })?;
        Ok(AuthoredAppendReceipt {
            matched_entities: core.matched_entities,
            document_work: core.document_work,
            selection_work: core.selection_work,
            encoding_work: core.encoding_work,
            binding_work: core.binding_work,
            work_used: core.work_used,
            accounted_peak_bytes: core.accounted_peak_bytes,
            candidate: candidate.qualification.verification,
        })
    }
}
