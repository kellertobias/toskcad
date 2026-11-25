//! Physical STEP of authored pipe hardware, written only after the actual
//! bytes requalify against the authored source:
//! - a clipped-NPT hex head is one mapped item carrying this body's exact
//!   occurrence chain, requalified against the clip plan, head intent and
//!   chain (standalone only; volumes are the sewn reference);
//! - any other hex head or hex nut keeps its retained radial thread tables
//!   unchanged beside exact hex stock, placed by the composed occurrence
//!   whose discrepancy from the chain is bounded and reported.
use super::*;
use multicad_rust_kernel::{
    axial_pipe_thread::authored::{AuthoredSection, StockClipPlan},
    radial_tensor::RadialBoundary,
    step_read::{
        HexStock, HexStockKind,
        npt_stock::{
            Error as HeadError,
            full_head::{
                self as head, AuthoredHeadReceipt, HexHeadIntent, OccurrenceFactor, SourceOptions,
            },
        },
    },
};
const STEP_OPERATION: &str = "pipe hardware STEP";
/// Which proof the receipt carries; its fields keep their units per method.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PipeHardwareStepMethod {
    /// Clipped-NPT hex head: one mapped item requalified from bytes.
    MappedClippedHead,
    /// Retained radial thread or bore tables decoded bit-identical beside
    /// exact hex stock, placed by the bounded composed occurrence.
    RetainedHexStock,
}
/// Volumes: the sewn reference (mapped clipped head) or the source net
/// interval at this occurrence (retained hex stock). Position error: the
/// sewn reference's (mapped) or the retained boundary's own error plus the
/// occurrence discrepancy (retained). Corner error: decoded head corners
/// (mapped) or written stock corners, against the authored hex.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PipeHardwareStepVerification {
    pub method: PipeHardwareStepMethod,
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
fn exchange_error(e: multicad_rust_kernel::step::StepError) -> GeometryError {
    fail(
        match e {
            multicad_rust_kernel::step::StepError::Limit => GeometryErrorKind::Limit,
            multicad_rust_kernel::step::StepError::Unsupported => GeometryErrorKind::Unsupported,
            _ => GeometryErrorKind::InvalidResult,
        },
        format!("Retained pipe hardware STEP readback failed:{e:?}"),
    )
}
fn receipt(r: AuthoredHeadReceipt) -> PipeHardwareStepVerification {
    PipeHardwareStepVerification {
        method: PipeHardwareStepMethod::MappedClippedHead,
        world_sewn_volume_mm3: r.world_sewn_volume_mm3,
        world_sewn_enclosure_mm: r.world_sewn_enclosure_mm,
        world_sewn_reference_position_error_mm: r.source.world_sewn_reference_position_error_mm,
        head_author_corner_error_mm: r.source.head_author_corner_error_mm,
        factor_chain_discrepancy_mm: r.source.factor_chain_discrepancy_mm,
        aggregate_work: r.source.aggregate_work,
        accounted_peak_bytes: r.source.accounted_peak_bytes,
    }
}
enum Exchange<'a> {
    Clipped(&'a StockClipPlan, HexHeadIntent),
    Retained {
        local: &'a RadialBoundary,
        kind: HexStockKind,
        across_flats_mm: f64,
        span_mm: [f64; 2],
    },
}
/// The placed retained boundary, its stock and the earned transport bounds.
struct Retained {
    boundary: Rc<RadialBoundary>,
    stock: HexStock,
    proof: PipeHardwareStepVerification,
}
fn upper(v: I) -> Result<f64, GeometryError> {
    v.checked()
        .map(|v| v.hi)
        .filter(|h| h.is_finite())
        .ok_or_else(|| {
            fail(
                GeometryErrorKind::Unsupported,
                "Pipe hardware STEP bound overflow",
            )
        })
}
impl PipeHardwareBody {
    fn exchange(&self) -> Option<Exchange<'_>> {
        let kind = self.source.intent.kind();
        if let Some(clip) = self.source.intent.clipped_thread() {
            // Clipped NPT-flat external hex head only.
            let PipeHardwareKind::HexHead {
                across_flats_mm,
                height_mm,
            } = kind
            else {
                return None;
            };
            return matches!(
                self.source.intent.thread().spec().section,
                AuthoredSection::NptFlat { .. }
            )
            .then_some(Exchange::Clipped(
                clip,
                HexHeadIntent {
                    across_flats_mm,
                    height_mm,
                },
            ));
        }
        let (local, pose) = self.source.thread.rust_radial_boundary()?;
        if pose != PlaneFrame::default() {
            return None;
        }
        let length = self.source.intent.thread().spec().length_mm;
        let (kind, across_flats_mm, span_mm) = match kind {
            PipeHardwareKind::HexHead {
                across_flats_mm,
                height_mm,
            } => (
                HexStockKind::Head,
                across_flats_mm,
                [length, length + height_mm],
            ),
            PipeHardwareKind::HexNut { across_flats_mm } => {
                (HexStockKind::Nut, across_flats_mm, [0., length])
            }
        };
        Some(Exchange::Retained {
            local,
            kind,
            across_flats_mm,
            span_mm,
        })
    }
    fn options(
        &self,
        tolerance: f64,
        policy: GeometryPolicy,
    ) -> Result<SourceOptions, GeometryError> {
        let tol = Self::tolerance(tolerance, policy)?;
        Ok(SourceOptions {
            position_error_mm: tol,
            material_tolerance_mm: tol,
            ..Default::default()
        })
    }
    fn tolerance(tolerance: f64, policy: GeometryPolicy) -> Result<f64, GeometryError> {
        if !tolerance.is_finite() || tolerance <= 0. {
            return Err(fail(
                GeometryErrorKind::Input,
                "Positive finite STEP tolerance required",
            ));
        }
        Ok(tolerance.min(policy.linear_tolerance_mm).min(1e-7))
    }
    pub(in crate::geometry) fn supports_step(&self) -> bool {
        self.exchange().is_some()
    }
    fn factors(&self) -> Vec<OccurrenceFactor> {
        self.placement.factor_rows()
    }
    /// Place the retained boundary and stock by the composed occurrence.
    /// The exact chain is y = A x + t (enclosed); the written frame is
    /// y = M x + o with M = [u v u*v] from the enclosure midpoints. For every
    /// local point |x| <= R the written image differs from the chain image by
    /// at most |t - o| + |A - M|_F R plus the f64 rounding of the transform.
    fn retained(
        &self,
        local: &RadialBoundary,
        kind: HexStockKind,
        across_flats_mm: f64,
        span_mm: [f64; 2],
        tolerance: f64,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Retained, GeometryError> {
        checkpoint(interrupt)?;
        let tol = Self::tolerance(tolerance, policy)?;
        let (columns, translation) = self.placement.composed();
        let o = midpoint(translation);
        let u = midpoint(columns[0]);
        let v = midpoint(columns[1]);
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let boundary = if self.placement.is_identity() {
            Rc::new(local.clone())
        } else {
            Rc::new(
                local
                    .transformed_with_interrupt(
                        [o, u, v],
                        policy.linear_tolerance_mm,
                        policy.minimum_feature_mm,
                        interrupt,
                    )
                    .map_err(radial_failure)?,
            )
        };
        checkpoint(interrupt)?;
        // Radius of every local source point: thread box and hex corners.
        let b = self.source.thread.bounds()?;
        let mut radius = I::ZERO;
        let mut corner = |p: [I; 3]| {
            let r = norm(p.map(|c| I::point(c.lo.abs().max(c.hi.abs()))));
            radius = I::new(radius.lo.max(r.lo), radius.hi.max(r.hi));
        };
        for c in 0..8 {
            corner(std::array::from_fn(|k| {
                I::point(if c >> k & 1 == 0 { b.min[k] } else { b.max[k] })
            }));
        }
        for p in self.source.hex_vertices {
            corner(p);
        }
        let written = [u, v, n];
        let mut frobenius = I::ZERO;
        for j in 0..3 {
            for k in 0..3 {
                frobenius = frobenius + (columns[j][k] - I::point(written[j][k])).square();
            }
        }
        let shift = norm(std::array::from_fn(|k| translation[k] - I::point(o[k])));
        let scale = o.iter().map(|x| x.abs()).sum::<f64>() + 3. * upper(radius)?;
        let rounding = I::point(8. * f64::EPSILON) * I::point(scale);
        // A sum of squares is nonnegative; outward subnormal rounding may
        // leave its enclosure's lower end just below zero.
        let frobenius = I::new(frobenius.lo.max(0.), frobenius.hi).sqrt();
        let discrepancy = upper(shift + frobenius * radius + rounding)?;
        if discrepancy > tol {
            return Err(fail(
                GeometryErrorKind::Unsupported,
                format!("The occurrence discrepancy {discrepancy:e} mm exceeds the STEP tolerance"),
            ));
        }
        let stock = HexStock {
            kind,
            origin: o,
            normal: n,
            reference: u,
            across_flats_mm,
            span_mm,
        };
        let rings = stock.corners().ok_or_else(|| {
            fail(
                GeometryErrorKind::Unsupported,
                "The placed hex stock is invalid",
            )
        })?;
        // Written corners against the authored hex at this occurrence.
        let mut corner_error = 0f64;
        for (i, p) in rings.iter().flatten().enumerate() {
            let exact = self.placement.point_bounds(self.source.hex_vertices[i])?;
            let d = norm(std::array::from_fn(|k| {
                I::point((exact[k].lo - p[k]).abs().max((exact[k].hi - p[k]).abs()))
            }));
            corner_error = corner_error.max(upper(d)?);
        }
        if corner_error > tol {
            return Err(fail(
                GeometryErrorKind::Unsupported,
                format!(
                    "The hex stock corner error {corner_error:e} mm exceeds the STEP tolerance"
                ),
            ));
        }
        let volume = self.source.volume * self.placement.determinant;
        let bounds = self.bounds()?;
        if !volume.valid() || volume.lo <= 0. {
            return Err(fail(
                GeometryErrorKind::Unsupported,
                "The pipe hardware volume is not resolved",
            ));
        }
        let proof = PipeHardwareStepVerification {
            method: PipeHardwareStepMethod::RetainedHexStock,
            world_sewn_volume_mm3: [volume.lo, volume.hi],
            world_sewn_enclosure_mm: std::array::from_fn(|k| [bounds.min[k], bounds.max[k]]),
            world_sewn_reference_position_error_mm: upper(
                I::point(local.boundary_error_mm()) + I::point(discrepancy),
            )?,
            head_author_corner_error_mm: corner_error,
            factor_chain_discrepancy_mm: discrepancy,
            aggregate_work: 0,
            accounted_peak_bytes: boundary.retained_budget_bytes(),
        };
        Ok(Retained {
            boundary,
            stock,
            proof,
        })
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
        let exchange = body.exchange()?;
        Some((|| match exchange {
            Exchange::Clipped(clip, intent) => {
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
            }
            Exchange::Retained {
                local,
                kind,
                across_flats_mm,
                span_mm,
            } => {
                let r = body.retained(
                    local,
                    kind,
                    across_flats_mm,
                    span_mm,
                    tolerance,
                    self.policy,
                    interrupt,
                )?;
                let tol = PipeHardwareBody::tolerance(tolerance, self.policy)?;
                multicad_rust_kernel::step::verify_pipe_hex_hardware(
                    text,
                    &r.boundary,
                    r.stock,
                    tol,
                )
                .map_err(exchange_error)?;
                checkpoint(interrupt)?;
                Ok(r.proof)
            }
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
        let exchange = body.exchange()?;
        Some((|| {
            let io = |e: std::io::Error| fail(GeometryErrorKind::InvalidResult, e.to_string());
            let candidate = match &exchange {
                Exchange::Clipped(clip, intent) => {
                    let options = body.options(tolerance, self.policy)?;
                    head::export_authored_hex_head(
                        clip,
                        *intent,
                        &body.factors(),
                        options,
                        interrupt,
                    )
                    .map_err(head_error)?
                    .into_string()
                }
                Exchange::Retained {
                    local,
                    kind,
                    across_flats_mm,
                    span_mm,
                } => {
                    let r = body.retained(
                        local,
                        *kind,
                        *across_flats_mm,
                        *span_mm,
                        tolerance,
                        self.policy,
                        interrupt,
                    )?;
                    let tol = PipeHardwareBody::tolerance(tolerance, self.policy)?;
                    multicad_rust_kernel::step::write_pipe_hex_hardware(&r.boundary, r.stock, tol)
                        .map_err(exchange_error)?
                }
            };
            checkpoint(interrupt)?;
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
            file.write_all(candidate.as_bytes())
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
            let proof = self
                .verify_pipe_hardware_step_source(&actual, tolerance, interrupt)
                .expect("typed authored pipe hardware")?;
            std::fs::rename(&staged, path).map_err(io)?;
            std::mem::forget(guard);
            Ok(proof)
        })())
    }
}
