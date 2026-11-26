//! Re-export of imported degree-eleven Phillips material. The original
//! retained tokens are only a proposal: the ACTUAL written bytes are decoded
//! independently and their distinctly repaired material must match the
//! retained source within the unchanged correspondence bounds.
//! The imported (identity) occurrence writes the original flat singleton. A
//! moved or instanced occurrence replays the same local records inside one
//! mapped item; the decoded file operator is enclosed and every bound is
//! transported through it to this occurrence.
use super::*;
use multicad_rust_kernel::{
    step::degree11_phillips as d11,
    step_read::PhillipsOwnedDegree11Material as OwnedDegree11Material,
};
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecodedDegree11PhillipsStepVerification {
    /// Independently decoded repaired material, not exact raw STEP mass.
    pub decoded_volume_bounds_mm3: [f64; 2],
    pub source_volume_bounds_mm3: [f64; 2],
    pub maximum_world_position_error_mm: f64,
    /// Original normalized root-cell/source-atlas chart derivatives.
    pub maximum_world_chart_derivative_error_mm: f64,
    pub independent_volume_discrepancy_mm3: f64,
    pub verified_tolerance_mm: f64,
    pub reader_reserved_work: usize,
    pub decoded_material_proof_work: usize,
    pub correspondence_work: usize,
    pub proposal_work: usize,
    pub aggregate_reserved_work: usize,
    pub known_resident_peak_bytes: usize,
}
fn d11_error(e: d11::Error) -> GeometryError {
    use multicad_rust_kernel::step_read::PhillipsReadError as R;
    let kind = match e {
        d11::Error::Invalid => GeometryErrorKind::Input,
        d11::Error::Limit => GeometryErrorKind::Limit,
        d11::Error::Cancelled => GeometryErrorKind::Cancelled,
        d11::Error::Unsupported => GeometryErrorKind::Unsupported,
        d11::Error::Reader(e) => match e.cause {
            R::Invalid => GeometryErrorKind::InvalidResult,
            R::Limit => GeometryErrorKind::Limit,
            R::Cancelled => GeometryErrorKind::Cancelled,
            R::Unsupported => GeometryErrorKind::Unsupported,
        },
    };
    error(
        "degree-eleven Phillips STEP",
        kind,
        format!("Independent decoded degree-eleven qualification failed:{e:?}"),
    )
}
fn d11_options(policy: GeometryPolicy, tolerance: f64) -> Result<d11::Options, GeometryError> {
    policy.validate()?;
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(fail(
            GeometryErrorKind::Input,
            "Positive finite exchange tolerance required",
        ));
    }
    let mut options = d11::Options::default();
    let tol = tolerance.min(policy.linear_tolerance_mm).min(1e-7);
    options.correspondence.maximum_distance_mm = tol;
    options.correspondence.maximum_derivative_mm = tol;
    Ok(options)
}
impl RepairedPhillipsBody {
    fn degree11_receipt(
        &self,
        proof: d11::Verification,
        frame: PlaneFrame,
        options: d11::Options,
    ) -> Result<DecodedDegree11PhillipsStepVerification, GeometryError> {
        Ok(DecodedDegree11PhillipsStepVerification {
            decoded_volume_bounds_mm3: proof.decoded_repaired_volume_bounds_mm3,
            source_volume_bounds_mm3: self.volume(frame)?,
            maximum_world_position_error_mm: proof.maximum_world_distance_mm,
            maximum_world_chart_derivative_error_mm: proof.maximum_world_chart_derivative_mm,
            independent_volume_discrepancy_mm3: proof.independent_repaired_volume_discrepancy_mm3,
            verified_tolerance_mm: options.correspondence.maximum_distance_mm,
            reader_reserved_work: proof.reader_reserved_work,
            decoded_material_proof_work: proof.actual_decoded_material_proof_work,
            correspondence_work: proof.fresh_correspondence_work,
            proposal_work: proof.proposal_work,
            aggregate_reserved_work: proof.aggregate_reserved_work,
            known_resident_peak_bytes: proof.known_resident_peak_bytes,
        })
    }
}
fn flat(frame: PlaneFrame) -> bool {
    frame == PlaneFrame::default()
}
fn qualify(
    source: &OwnedDegree11Material,
    text: &str,
    frame: PlaneFrame,
    options: d11::Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<d11::QualifiedMaterial, GeometryError> {
    if flat(frame) {
        d11::qualify_actual_file(source, text, options, interrupt)
    } else {
        d11::qualify_placed_actual_file(source, text, occurrence(frame), options, interrupt)
    }
    .map_err(d11_error)
}
impl Solid {
    /// Whether this solid re-exports its imported degree-eleven material.
    pub fn has_degree11_phillips_step_source(&self) -> bool {
        self.repaired_phillips_body()
            .is_some_and(|(body, _)| body.source.degree11().is_some())
    }
    /// Fresh qualification of ACTUAL bytes against the retained material.
    pub fn verify_degree11_phillips_step_source(
        &self,
        text: &str,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<DecodedDegree11PhillipsStepVerification, GeometryError>> {
        let (body, frame) = self.repaired_phillips_body()?;
        let source = body.source.degree11()?;
        Some((|| {
            checkpoint(interrupt)?;
            let options = d11_options(self.policy, tolerance)?;
            let verified = qualify(source, text, frame, options, interrupt)?;
            checkpoint(interrupt)?;
            body.degree11_receipt(verified.verification(), frame, options)
        })())
    }
    /// Atomically replace a destination only after independently qualifying
    /// the actual synchronized temporary-file bytes against the source.
    pub fn write_degree11_phillips_step_with_verification(
        &self,
        path: &Path,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<DecodedDegree11PhillipsStepVerification, GeometryError>> {
        let (body, frame) = self.repaired_phillips_body()?;
        let source = body.source.degree11()?;
        Some((|| {
            checkpoint(interrupt)?;
            let options = d11_options(self.policy, tolerance)?;
            let candidate = if flat(frame) {
                d11::export_retained(source, options, interrupt)
            } else {
                d11::export_retained_placed(source, occurrence(frame), options, interrupt)
            }
            .map_err(d11_error)?;
            let (tmp, mut file) = temporary(path)?;
            file.write_all(candidate.as_str().as_bytes())
                .and_then(|()| file.sync_all())
                .map_err(|e| fail(GeometryErrorKind::InvalidResult, e.to_string()))?;
            drop(file);
            drop(candidate);
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
            let verified = qualify(source, &actual, frame, options, interrupt)?;
            let receipt = body.degree11_receipt(verified.verification(), frame, options)?;
            checkpoint(interrupt)?;
            fs::rename(&tmp.0, path)
                .map_err(|e| fail(GeometryErrorKind::InvalidResult, e.to_string()))?;
            Ok(receipt)
        })())
    }
}
