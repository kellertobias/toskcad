//! Physical STEP of an authored clipped-NPT hex-head screw: one mapped item
//! carrying this body's exact occurrence chain, written only after the actual
//! bytes requalify against the authored clip plan, head intent and chain.
//! Standalone only; volumes are the sewn reference, not raw occupied mass.
use super::*;
use multicad_rust_kernel::{
    axial_pipe_thread::authored::{AuthoredSection, StockClipPlan},
    step_read::npt_stock::{
        Error as HeadError,
        full_head::{
            self as head, AuthoredHeadReceipt, HexHeadIntent, OccurrenceFactor, SourceOptions,
        },
    },
};
const STEP_OPERATION: &str = "pipe hardware STEP";
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PipeHardwareStepVerification {
    pub world_sewn_volume_mm3: [f64; 2],
    pub world_sewn_enclosure_mm: [[f64; 2]; 3],
    pub world_sewn_reference_position_error_mm: f64,
    pub head_author_corner_error_mm: f64,
    pub factor_chain_discrepancy_mm: f64,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
}
fn fail(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error(STEP_OPERATION, kind, message)
}
fn head_error(e: HeadError) -> GeometryError {
    fail(
        match e {
            HeadError::Invalid => GeometryErrorKind::InvalidResult,
            HeadError::Unsupported => GeometryErrorKind::Unsupported,
            HeadError::Limit => GeometryErrorKind::Limit,
            HeadError::Cancelled => GeometryErrorKind::Cancelled,
        },
        format!("Authored pipe hardware STEP qualification failed:{e:?}"),
    )
}
fn receipt(r: AuthoredHeadReceipt) -> PipeHardwareStepVerification {
    PipeHardwareStepVerification {
        world_sewn_volume_mm3: r.world_sewn_volume_mm3,
        world_sewn_enclosure_mm: r.world_sewn_enclosure_mm,
        world_sewn_reference_position_error_mm: r.source.world_sewn_reference_position_error_mm,
        head_author_corner_error_mm: r.source.head_author_corner_error_mm,
        factor_chain_discrepancy_mm: r.source.factor_chain_discrepancy_mm,
        aggregate_work: r.source.aggregate_work,
        accounted_peak_bytes: r.source.accounted_peak_bytes,
    }
}
impl PipeHardwareBody {
    /// Clipped NPT-flat external hex head only; nuts, BSPT and unclipped
    /// sources are outside this exchange family.
    fn step_source(&self) -> Option<(&StockClipPlan, HexHeadIntent)> {
        let PipeHardwareKind::HexHead {
            across_flats_mm,
            height_mm,
        } = self.source.intent.kind()
        else {
            return None;
        };
        let clip = self.source.intent.clipped_thread()?;
        matches!(
            self.source.intent.thread().spec().section,
            AuthoredSection::NptFlat { .. }
        )
        .then_some((
            clip,
            HexHeadIntent {
                across_flats_mm,
                height_mm,
            },
        ))
    }
    fn options(
        &self,
        tolerance: f64,
        policy: GeometryPolicy,
    ) -> Result<SourceOptions, GeometryError> {
        if !tolerance.is_finite() || tolerance <= 0. {
            return Err(fail(
                GeometryErrorKind::Input,
                "Positive finite STEP tolerance required",
            ));
        }
        let tol = tolerance.min(policy.linear_tolerance_mm).min(1e-7);
        Ok(SourceOptions {
            position_error_mm: tol,
            material_tolerance_mm: tol,
            ..Default::default()
        })
    }
    pub(in crate::geometry) fn supports_step(&self) -> bool {
        self.step_source().is_some()
    }
    fn factors(&self) -> Vec<OccurrenceFactor> {
        self.placement.factor_rows()
    }
}
impl Solid {
    /// Whether this standalone solid writes authored pipe hardware STEP.
    pub fn has_pipe_hardware_step_source(&self) -> bool {
        self.pipe_hardware_body()
            .is_some_and(PipeHardwareBody::supports_step)
    }
    /// Fresh proof of ACTUAL bytes against this body's authored source.
    pub fn verify_pipe_hardware_step_source(
        &self,
        text: &str,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<PipeHardwareStepVerification, GeometryError>> {
        let body = self.pipe_hardware_body()?;
        let (clip, intent) = body.step_source()?;
        Some((|| {
            let options = body.options(tolerance, self.policy)?;
            head::qualify_authored_hex_head_file(
                text,
                clip,
                intent,
                &body.factors(),
                options,
                interrupt,
            )
            .map(receipt)
            .map_err(head_error)
        })())
    }
    /// Write via a synchronized sibling, requalify the ACTUAL read-back bytes,
    /// then rename. An existing destination is untouched on any failure.
    pub fn write_pipe_hardware_step_with_verification(
        &self,
        path: &std::path::Path,
        tolerance: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<PipeHardwareStepVerification, GeometryError>> {
        use std::io::{Read as _, Write as _};
        let body = self.pipe_hardware_body()?;
        let (clip, intent) = body.step_source()?;
        Some((|| {
            let io = |e: std::io::Error| fail(GeometryErrorKind::InvalidResult, e.to_string());
            let options = body.options(tolerance, self.policy)?;
            let candidate =
                head::export_authored_hex_head(clip, intent, &body.factors(), options, interrupt)
                    .map_err(head_error)?;
            let name = path
                .file_name()
                .ok_or_else(|| fail(GeometryErrorKind::Input, "Output filename required"))?;
            let mut staged = std::ffi::OsString::from(".");
            staged.push(name);
            staged.push(format!(".{}.pipe-step-tmp", std::process::id()));
            let staged = path.with_file_name(staged);
            struct Remove(std::path::PathBuf);
            impl Drop for Remove {
                fn drop(&mut self) {
                    let _ = std::fs::remove_file(&self.0);
                }
            }
            let guard = Remove(staged.clone());
            let mut file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&staged)
                .map_err(io)?;
            file.write_all(candidate.as_str().as_bytes())
                .and_then(|()| file.sync_all())
                .map_err(io)?;
            drop(file);
            drop(candidate);
            let mut actual = String::new();
            std::fs::File::open(&staged)
                .and_then(|f| f.take(64 * 1024 * 1024 + 1).read_to_string(&mut actual))
                .map_err(io)?;
            if actual.len() > 64 * 1024 * 1024 {
                return Err(fail(
                    GeometryErrorKind::Limit,
                    "Actual output exceeds file quota",
                ));
            }
            let proof = head::qualify_authored_hex_head_file(
                &actual,
                clip,
                intent,
                &body.factors(),
                options,
                interrupt,
            )
            .map_err(head_error)?;
            std::fs::rename(&staged, path).map_err(io)?;
            std::mem::forget(guard);
            Ok(receipt(proof))
        })())
    }
}
