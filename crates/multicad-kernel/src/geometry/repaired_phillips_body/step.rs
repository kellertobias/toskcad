//! Physical exchange via an independently decoded, distinctly repaired source.
//! Raw-file exact material and unit-normal correspondence are not asserted.
use super::*;
use multicad_rust_kernel::{step::periodic_phillips as exchange, step_read::owned_insert as read};
use std::{
    fs::{self, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};
const MAX_FILE: usize = 64 * 1024 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PeriodicPhillipsStepCandidateResources {
    pub proposal_work: usize,
    pub reader_reserved_work: usize,
    pub decoded_material_proof_work: usize,
    pub correspondence_fixed_work: usize,
    pub correspondence_cavity_cells: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DecodedPeriodicPhillipsStepVerification {
    pub represented_member_count: usize,
    pub member_index: usize,
    pub decoded_bounds_mm: [[f64; 3]; 2],
    pub decoded_volume_bounds_mm3: [f64; 2],
    pub source_volume_bounds_mm3: [f64; 2],
    pub source_occurrence: [[f64; 3]; 4],
    pub maximum_world_position_error_mm: f64,
    /// Complete paired normalized chart derivatives, not unit-normal angles.
    pub maximum_world_chart_derivative_error_mm: f64,
    pub independent_volume_discrepancy_mm3: f64,
    pub verified_tolerance_mm: f64,
    pub reader_reserved_work: usize,
    pub decoded_material_proof_work: usize,
    pub correspondence_fixed_work: usize,
    pub correspondence_cavity_cells: usize,
    /// Additional neutral document scope, separate from reader-owned scopes.
    pub document_work: usize,
    pub document_retained_bytes: usize,
    pub selection_work: usize,
    pub selection_query_peak_bytes: usize,
    pub binding_work: usize,
    pub binding_query_peak_bytes: usize,
    pub matched_entities: usize,
    /// Present only when this operation also built/qualified a new export
    /// candidate before the fresh ACTUAL-file scope reported above. The two
    /// reservations and adaptive-cell counts must not be silently conflated.
    pub export_candidate: Option<super::PeriodicPhillipsStepCandidateResources>,
}
fn fail(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error("periodic Phillips STEP", kind, message)
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(fail(GeometryErrorKind::Cancelled, "Exchange superseded"))
    } else {
        Ok(())
    }
}
fn exchange_error(e: exchange::Error) -> GeometryError {
    use multicad_rust_kernel::step_read::PhillipsReadError as R;
    let kind = match e {
        exchange::Error::Invalid => GeometryErrorKind::Input,
        exchange::Error::Limit => GeometryErrorKind::Limit,
        exchange::Error::Cancelled => GeometryErrorKind::Cancelled,
        exchange::Error::Unsupported => GeometryErrorKind::Unsupported,
        exchange::Error::Reader(e) => match e.cause {
            R::Invalid => GeometryErrorKind::InvalidResult,
            R::Limit => GeometryErrorKind::Limit,
            R::Cancelled => GeometryErrorKind::Cancelled,
            R::Unsupported => GeometryErrorKind::Unsupported,
        },
    };
    fail(
        kind,
        format!("Independent decoded periodic repair qualification failed:{e:?}"),
    )
}
fn document_error(e: read::Error) -> GeometryError {
    fail(
        match e {
            read::Error::Invalid => GeometryErrorKind::InvalidResult,
            read::Error::Limit => GeometryErrorKind::Limit,
            read::Error::Cancelled => GeometryErrorKind::Cancelled,
            read::Error::Unsupported => GeometryErrorKind::Unsupported,
        },
        format!("Actual represented STEP document failed:{e:?}"),
    )
}
fn options(policy: GeometryPolicy, tolerance: f64) -> Result<exchange::Options, GeometryError> {
    policy.validate()?;
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(fail(
            GeometryErrorKind::Input,
            "Positive finite exchange tolerance required",
        ));
    }
    let mut options = exchange::Options::default();
    let tol = tolerance.min(policy.linear_tolerance_mm).min(1e-7);
    options.correspondence.maximum_distance_mm = tol;
    options.correspondence.maximum_derivative_mm = tol;
    Ok(options)
}
fn occurrence(frame: PlaneFrame) -> [[f64; 3]; 4] {
    [frame.origin, frame.u, frame.v, frame.normal()]
}
impl RepairedPhillipsBody {
    fn verify_singleton_step(
        &self,
        text: &str,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<DecodedPeriodicPhillipsStepVerification, GeometryError> {
        checkpoint(interrupt)?;
        let options = options(policy, tolerance)?;
        let source = self
            .source
            .periodic()
            .ok_or_else(|| unsupported("No periodic source"))?;
        let document = read::CompoundDocument::parse(text, Default::default(), interrupt)
            .map_err(document_error)?;
        if document.represented_members().len() != 1
            || document.represented_members()[0].mapping_link().is_some()
        {
            return Err(fail(
                GeometryErrorKind::Unsupported,
                "Exactly one ordinary represented member required",
            ));
        }
        let verified =
            exchange::qualify_actual_file(source, text, occurrence(frame), options, interrupt)
                .map_err(exchange_error)?;
        let material = verified.material();
        let b = periodic_enclosure_bounds(material);
        let bounds = Placement::new(material.chart(), PlaneFrame::default())?.bounds(Bounds {
            min: b.map(|x| x[0]),
            max: b.map(|x| x[1]),
        })?;
        let proof = verified.verification();
        checkpoint(interrupt)?;
        Ok(DecodedPeriodicPhillipsStepVerification {
            represented_member_count: 1,
            member_index: 0,
            decoded_bounds_mm: [bounds.min, bounds.max],
            decoded_volume_bounds_mm3: proof.decoded_volume_bounds_mm3,
            source_volume_bounds_mm3: self.volume(frame)?,
            source_occurrence: occurrence(frame),
            maximum_world_position_error_mm: proof.maximum_world_distance_mm,
            maximum_world_chart_derivative_error_mm: proof.maximum_world_chart_derivative_mm,
            independent_volume_discrepancy_mm3: proof.independent_volume_discrepancy_mm3,
            verified_tolerance_mm: options.correspondence.maximum_distance_mm,
            reader_reserved_work: proof.reserved_reader_work,
            decoded_material_proof_work: proof.decoded_material_proof_work,
            correspondence_fixed_work: proof.correspondence_fixed_work,
            correspondence_cavity_cells: proof.correspondence_cavity_cells,
            document_work: document.work_used(),
            document_retained_bytes: document.retained_budget_bytes(),
            selection_work: 0,
            selection_query_peak_bytes: 0,
            binding_work: 0,
            binding_query_peak_bytes: 0,
            matched_entities: 0,
            export_candidate: None,
        })
    }
}
struct Temporary(PathBuf);
impl Drop for Temporary {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}
fn temporary(path: &Path) -> Result<(Temporary, fs::File), GeometryError> {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let name = path
        .file_name()
        .ok_or_else(|| fail(GeometryErrorKind::Input, "Output filename required"))?;
    for _ in 0..32 {
        let mut tmp = std::ffi::OsString::from(".");
        tmp.push(name);
        tmp.push(format!(
            ".{}.{}.step-tmp",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let actual = parent.join(tmp);
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&actual)
        {
            Ok(file) => return Ok((Temporary(actual), file)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(e) => return Err(fail(GeometryErrorKind::InvalidResult, e.to_string())),
        }
    }
    Err(fail(
        GeometryErrorKind::Limit,
        "Temporary output namespace exhausted",
    ))
}
impl Solid {
    pub fn verify_periodic_phillips_step_source(
        &self,
        text: &str,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<super::DecodedPeriodicPhillipsStepVerification, GeometryError>> {
        let (body, frame) = self.repaired_phillips_body()?;
        body.source.periodic()?;
        Some(body.verify_singleton_step(text, frame, self.policy, tolerance, interrupt))
    }
    /// Atomically replace a destination only after independently qualifying the
    /// actual synchronized temporary-file bytes against the retained source.
    pub fn write_periodic_phillips_step_with_verification(
        &self,
        path: &Path,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<super::DecodedPeriodicPhillipsStepVerification, GeometryError>> {
        let (body, frame) = self.repaired_phillips_body()?;
        let source = body.source.periodic()?;
        Some((|| {
            checkpoint(interrupt)?;
            let candidate = exchange::export(
                source,
                occurrence(frame),
                options(self.policy, tolerance)?,
                interrupt,
            )
            .map_err(exchange_error)?;
            let (tmp, mut file) = temporary(path)?;
            file.write_all(candidate.as_str().as_bytes())
                .and_then(|()| file.sync_all())
                .map_err(|e| fail(GeometryErrorKind::InvalidResult, e.to_string()))?;
            drop(file);
            checkpoint(interrupt)?;
            let mut actual = String::new();
            fs::File::open(&tmp.0)
                .and_then(|f| f.take((MAX_FILE + 1) as u64).read_to_string(&mut actual))
                .map_err(|e| fail(GeometryErrorKind::InvalidResult, e.to_string()))?;
            if actual.len() > MAX_FILE {
                return Err(fail(
                    GeometryErrorKind::Limit,
                    "Actual output exceeds file quota",
                ));
            }
            let mut receipt =
                body.verify_singleton_step(&actual, frame, self.policy, tolerance, interrupt)?;
            let candidate_proof = candidate.qualification().verification();
            receipt.export_candidate = Some(PeriodicPhillipsStepCandidateResources {
                proposal_work: candidate.proposal_work_used(),
                reader_reserved_work: candidate_proof.reserved_reader_work,
                decoded_material_proof_work: candidate_proof.decoded_material_proof_work,
                correspondence_fixed_work: candidate_proof.correspondence_fixed_work,
                correspondence_cavity_cells: candidate_proof.correspondence_cavity_cells,
            });
            checkpoint(interrupt)?;
            fs::rename(&tmp.0, path)
                .map_err(|e| fail(GeometryErrorKind::InvalidResult, e.to_string()))?;
            Ok(receipt)
        })())
    }
}
mod degree11;
pub use degree11::DecodedDegree11PhillipsStepVerification;
mod members;
#[cfg(test)]
pub(in crate::geometry) mod tests;
