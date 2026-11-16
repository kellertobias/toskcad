//! Standalone actual-import re-export. Physical source records remain exact;
//! the declared sewn reference is compared under original B then optional A.
use super::*;
use read::full_head::transport::replay::{self, ReplayOptions};
use std::{
    fs::{self, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};
const NATIVE_WORK: usize = 1024;
const IO_HEADER: usize = 8192;
#[derive(Clone, Copy, Debug)]
pub struct ActualHeadStepVerification {
    pub original_imported_sewn_volume_bounds_mm3: [f64; 2],
    pub output_sewn_volume_bounds_mm3: [f64; 2],
    pub output_sewn_bounds_mm: [[f64; 2]; 3],
    pub maximum_world_position_error_mm: f64,
    pub source_work: usize,
    pub occurrence_work: usize,
    /// Actual candidate byte encoding plus its separately counted operator proposal.
    pub candidate_preparation_work: usize,
    pub operator_comparison_work: usize,
    pub actual_output_binding_work: usize,
    pub native_and_io_work: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
}
fn checkpoint(stop: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if stop.is_some_and(|f| f()) {
        Err(failure(read::Error::Cancelled))
    } else {
        Ok(())
    }
}
fn add(a: usize, b: usize) -> Result<usize, GeometryError> {
    a.checked_add(b).ok_or_else(|| failure(read::Error::Limit))
}
fn subtract(a: usize, b: usize) -> Result<usize, GeometryError> {
    a.checked_sub(b)
        .filter(|n| *n > 0)
        .ok_or_else(|| failure(read::Error::Limit))
}
impl ActualClippedHeadBody {
    fn step_options(&self, mut options: ReplayOptions) -> Result<ReplayOptions, GeometryError> {
        if options.resources.maximum_work == 0
            || options.resources.maximum_work > 272_000_000
            || options.resources.maximum_bytes == 0
            || options.resources.maximum_bytes > BYTES
            || !options.maximum_position_error_mm.is_finite()
            || options.maximum_position_error_mm <= 0.
            || options.maximum_position_error_mm > 1e-4
        {
            return Err(failure(read::Error::Invalid));
        }
        options.maximum_position_error_mm = options
            .maximum_position_error_mm
            .min(self.policy.linear_tolerance_mm);
        options.resources.maximum_work = subtract(
            options.resources.maximum_work,
            add(self.preparation, NATIVE_WORK)?,
        )?;
        options.resources.additional_resident_bytes = add(
            options.resources.additional_resident_bytes,
            self.external_bytes(),
        )?;
        Ok(options)
    }
    fn step_receipt(
        &self,
        proof: &replay::VerifiedReplay<'_, '_>,
        native: usize,
        encoding: usize,
        peak: usize,
    ) -> Result<ActualHeadStepVerification, GeometryError> {
        let p = proof.receipt();
        let v = self.source.world_sewn_volume_mm3();
        let v = self.occurrence.as_ref().map_or(I::new(v[0], v[1]), |a| {
            let d = a.determinant_bounds();
            I::new(v[0], v[1]) * I::new(d[0], d[1])
        });
        Ok(ActualHeadStepVerification {
            original_imported_sewn_volume_bounds_mm3: [v.lo, v.hi],
            output_sewn_volume_bounds_mm3: proof.actual().world_sewn_volume_mm3(),
            output_sewn_bounds_mm: proof.actual().world_sewn_enclosure_mm(),
            maximum_world_position_error_mm: p.maximum_world_position_error_mm,
            source_work: p.inherited_source_work,
            occurrence_work: p.inherited_occurrence_work,
            candidate_preparation_work: encoding,
            operator_comparison_work: p.operator_comparison_work,
            actual_output_binding_work: p.fresh_output_binding_work,
            native_and_io_work: native,
            aggregate_work: add(add(p.aggregate_work, native)?, encoding)?,
            accounted_peak_bytes: p.accounted_peak_bytes.max(peak),
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
        .ok_or_else(|| input("STEP output filename required"))?;
    for _ in 0..32 {
        let mut name = std::ffi::OsString::from(name);
        name.push(format!(
            ".{}.{}.actual-head-tmp",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let p = parent.join(name);
        match OpenOptions::new().create_new(true).write(true).open(&p) {
            Ok(file) => return Ok((Temporary(p), file)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(e) => {
                return Err(error(
                    "actual head STEP",
                    GeometryErrorKind::InvalidResult,
                    e.to_string(),
                ));
            }
        }
    }
    Err(failure(read::Error::Limit))
}
impl Solid {
    /// Actual supplied output reader and full exact local source/operator bind.
    /// Additional resident bytes exclude this body and the explicit input Arc.
    pub fn verify_actual_clipped_head_step(
        &self,
        input: Arc<str>,
        options: ReplayOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<ActualHeadStepVerification, GeometryError>> {
        let b = self.actual_clipped_head_body()?;
        Some((|| {
            checkpoint(interrupt)?;
            let o = b.step_options(options)?;
            let proof = replay::qualify(&b.source, b.occurrence.as_ref(), input, o, interrupt)
                .map_err(failure)?;
            b.step_receipt(&proof, add(b.preparation, NATIVE_WORK)?, 0, 0)
        })())
    }
    /// Stage, sync, independently decode ACTUAL read-back bytes, then rename.
    /// Only one output material bind executes; source/A prefixes are never reset.
    pub fn write_actual_clipped_head_step_with_verification(
        &self,
        path: &Path,
        options: ReplayOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<ActualHeadStepVerification, GeometryError>> {
        let b = self.actual_clipped_head_body()?;
        Some((|| {
            checkpoint(interrupt)?;
            crate::step_path(path).map_err(|e| input(e.to_string()))?;
            let mut o = b.step_options(options)?;
            let candidate =
                replay::propose(&b.source, b.occurrence.as_ref(), o, interrupt).map_err(failure)?;
            let preparation = candidate.receipt();
            let bytes = candidate.into_bytes();
            let len = bytes.len();
            // Five real full-byte actions: write, zero-initialize read buffer,
            // read, UTF8 validate, Arc copy. Constant covers filesystem calls,
            // EOF probe/metadata/receipt setup and publication bookkeeping.
            let io = add(
                len.checked_mul(5)
                    .ok_or_else(|| failure(read::Error::Limit))?,
                IO_HEADER,
            )?;
            let encoded = add(preparation.fresh_encoding_work, preparation.proposal_work)?;
            o.resources.maximum_work = subtract(o.resources.maximum_work, add(io, encoded)?)?;
            // Candidate Arc is dropped before allocating read storage. At Arc
            // conversion only Vec+newArc overlap, plus the complete body/caches.
            let mut io_peak = add(
                add(
                    add(
                        b.resident_bytes(),
                        options.resources.additional_resident_bytes,
                    )?,
                    len.checked_mul(2)
                        .ok_or_else(|| failure(read::Error::Limit))?,
                )?,
                IO_HEADER + 4 * std::mem::size_of::<usize>(),
            )?;
            if io_peak > options.resources.maximum_bytes {
                return Err(failure(read::Error::Limit));
            }
            // Reserve the remaining bind before ANY file side effect.
            let prefix = add(
                b.source.receipt().aggregate_work,
                b.occurrence
                    .as_ref()
                    .map_or(0, read::LiteralOccurrencePlacement::work_used),
            )?;
            if o.resources.maximum_work <= prefix {
                return Err(failure(read::Error::Limit));
            }
            checkpoint(interrupt)?;
            let (tmp, mut file) = temporary(path)?;
            file.write_all(bytes.as_bytes())
                .and_then(|()| file.sync_all())
                .map_err(|e| {
                    error(
                        "actual head STEP",
                        GeometryErrorKind::InvalidResult,
                        e.to_string(),
                    )
                })?;
            drop(file);
            drop(bytes);
            checkpoint(interrupt)?;
            let mut file = fs::File::open(&tmp.0).map_err(|e| {
                error(
                    "actual head STEP",
                    GeometryErrorKind::InvalidResult,
                    e.to_string(),
                )
            })?;
            if file.metadata().map_err(|e| input(e.to_string()))?.len() != len as u64 {
                return Err(failure(read::Error::Unsupported));
            }
            let mut actual = Vec::new();
            actual
                .try_reserve_exact(len)
                .map_err(|_| failure(read::Error::Limit))?;
            let actual_io_peak = add(
                add(
                    b.resident_bytes(),
                    options.resources.additional_resident_bytes,
                )?,
                add(
                    add(actual.capacity(), len)?,
                    IO_HEADER + 4 * std::mem::size_of::<usize>(),
                )?,
            )?;
            if actual_io_peak > options.resources.maximum_bytes {
                return Err(failure(read::Error::Limit));
            }
            io_peak = io_peak.max(actual_io_peak);
            actual.resize(len, 0);
            file.read_exact(&mut actual).map_err(|e| {
                error(
                    "actual head STEP",
                    GeometryErrorKind::InvalidResult,
                    e.to_string(),
                )
            })?;
            let mut end = [0u8; 1];
            if file.read(&mut end).map_err(|e| input(e.to_string()))? != 0 {
                return Err(failure(read::Error::Unsupported));
            }
            let text = std::str::from_utf8(&actual).map_err(|_| failure(read::Error::Invalid))?;
            let actual_text: Arc<str> = Arc::from(text);
            drop(actual);
            drop(file);
            checkpoint(interrupt)?;
            let proof =
                replay::qualify(&b.source, b.occurrence.as_ref(), actual_text, o, interrupt)
                    .map_err(failure)?;
            let native = add(add(b.preparation, NATIVE_WORK)?, io)?;
            let receipt = b.step_receipt(
                &proof,
                native,
                encoded,
                io_peak.max(preparation.accounted_peak_bytes),
            )?;
            if receipt.aggregate_work > options.resources.maximum_work
                || receipt.accounted_peak_bytes > options.resources.maximum_bytes
            {
                return Err(failure(read::Error::Limit));
            }
            checkpoint(interrupt)?;
            fs::rename(&tmp.0, path).map_err(|e| {
                error(
                    "actual head STEP",
                    GeometryErrorKind::InvalidResult,
                    e.to_string(),
                )
            })?;
            Ok(receipt)
        })())
    }
}
