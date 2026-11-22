//! A trimmed-tensor radial-graph solid read back from STEP text alone.
//!
//! `step_read::owned_hanger` decodes the literal file (a product item or one
//! declared rigid `MAPPED_ITEM`), recovers its exact pole grid, polar
//! supports and whole finite edge lifts, heals the cap sewing and integrates
//! the occupied material. That is the same independent proof the qualified
//! hanger export runs on its own bytes before writing them; no source recipe,
//! header or label is read. The admitted solid keeps the healed volume
//! interval and a conservative radial envelope. Its boundary is not rebuilt,
//! so meshes, faces, topology and re-export are typed refusals.
use super::*;
use multicad_rust_kernel::step_read::owned_hanger as read;

const OPERATION: &str = "STEP import";
/// Same widths and stage budgets as the qualified hanger export's own readback.
const MAXIMUM_VOLUME_WIDTH_MM3: f64 = 1e-6;

pub(crate) struct HealedImportBody {
    pub(super) feature: EntityId,
    volume: [f64; 2],
    radius_mm: f64,
    heights_mm: [f64; 2],
    sewing_mm: f64,
}

fn failure(stage: &str, cause: read::Error) -> GeometryError {
    error(
        OPERATION,
        match cause {
            read::Error::Invalid => GeometryErrorKind::InvalidResult,
            read::Error::Unsupported => GeometryErrorKind::Unsupported,
            read::Error::Limit => GeometryErrorKind::Limit,
            read::Error::Cancelled => GeometryErrorKind::Cancelled,
        },
        format!("Healed trimmed-tensor STEP material failed at {stage}: {cause:?}"),
    )
}

impl HealedImportBody {
    pub(super) fn volume_mm3(&self) -> f64 {
        self.volume[0].midpoint(self.volume[1])
    }

    pub(super) fn volume_bounds_mm3(&self) -> [f64; 2] {
        self.volume
    }

    pub(super) fn sewing_error_bound_mm(&self) -> f64 {
        self.sewing_mm
    }

    /// Box of the local cylinder `r <= R`, `z0 <= z <= z1` under the frame.
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Bounds {
        let n = frame.normal();
        let [z0, z1] = self.heights_mm;
        let ends =
            [z0, z1].map(|z| std::array::from_fn::<f64, 3, _>(|i| frame.origin[i] + n[i] * z));
        // A few rounding units for the declared operator's projected basis.
        let slack = 1e-11 * (self.radius_mm + z0.abs().max(z1.abs()) + 1.);
        let half = |i: usize| self.radius_mm * (1. - n[i] * n[i]).max(0.).sqrt() + slack;
        Bounds {
            min: std::array::from_fn(|i| ends[0][i].min(ends[1][i]) - half(i)),
            max: std::array::from_fn(|i| ends[0][i].max(ends[1][i]) + half(i)),
        }
    }

    pub(crate) fn unsupported(what: &'static str) -> GeometryError {
        error(
            what,
            GeometryErrorKind::Unsupported,
            "The imported healed trimmed-tensor solid keeps its proved volume but has no rebuilt boundary",
        )
    }
}

impl Solid {
    pub(crate) fn healed_import_body(&self) -> Option<(&HealedImportBody, PlaneFrame)> {
        match &self.body {
            SolidBody::HealedImport(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }

    /// The proved healed volume interval and cap sewing bound of an
    /// imported trimmed-tensor solid.
    pub fn healed_import_proof(&self) -> Option<([f64; 2], f64)> {
        self.healed_import_body()
            .map(|(body, _)| (body.volume_bounds_mm3(), body.sewing_error_bound_mm()))
    }

    /// `Ok(None)` when the text is not a single literal trimmed-tensor
    /// radial-graph solid; the caller then keeps its own refusal.
    pub(super) fn import_healed_tensor_step(
        input: &str,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Option<Self>, GeometryError> {
        let decoded = match read::decode(
            input,
            read::ReadOptions {
                maximum_work: 128_000_000,
                retained_budget_bytes: 64 * 1024 * 1024,
            },
            None,
        ) {
            Ok(decoded) => decoded,
            Err(read::Error::Unsupported | read::Error::Invalid) => return Ok(None),
            Err(cause) => return Err(failure("decode", cause)),
        };
        let bytes = 64 * 1024 * 1024;
        let grid = read::ExactPoleGrid::prepare(
            &decoded,
            read::GridOptions {
                maximum_preparation_work: 224_000_000,
                maximum_bytes: bytes,
            },
            None,
        )
        .map_err(|e| failure("pole grid", e))?;
        let polar = read::PolarSupports::prepare(
            &grid,
            read::PolarOptions {
                maximum_work: 64_000_000,
                maximum_bytes: bytes,
            },
            None,
        )
        .map_err(|e| failure("polar supports", e))?;
        let lifts = decoded
            .surface_lifts(
                read::LiftOptions {
                    maximum_work: 160_000_000,
                    maximum_bytes: bytes,
                    ..read::LiftOptions::default()
                },
                None,
            )
            .map_err(|e| failure("edge lifts", e))?;
        let healed = read::HealedMaterial::prepare(
            &polar,
            &lifts,
            read::HealingOptions {
                maximum_work: 64_000_000,
                maximum_bytes: bytes,
                ..read::HealingOptions::default()
            },
            None,
        )
        .map_err(|e| failure("healing", e))?;
        let volume = healed
            .healed_volume(
                read::VolumeOptions {
                    maximum_work: 16_000_000,
                    maximum_bytes: bytes,
                    maximum_width_mm3: MAXIMUM_VOLUME_WIDTH_MM3,
                },
                None,
            )
            .map_err(|e| failure("volume", e))?
            .bounds_mm3();
        let sewing = healed.sewing_error_bound_mm();
        if sewing > policy.linear_tolerance_mm {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Unsupported,
                "Healed STEP cap sewing exceeds the requested tolerance",
            ));
        }
        let radius = healed
            .component_radius_bounds_mm()
            .into_iter()
            .map(|r| r[1])
            .fold(0., f64::max);
        let heights = polar
            .plane_heights_mm()
            .map_err(|e| failure("cap planes", e))?;
        let low = heights.iter().map(|h| h[0]).fold(f64::INFINITY, f64::min);
        let high = heights
            .iter()
            .map(|h| h[1])
            .fold(f64::NEG_INFINITY, f64::max);
        let positive = |x: f64| x.is_finite() && x > 0.;
        if !positive(volume[0])
            || !positive(volume[1])
            || !positive(radius)
            || !positive(high - low)
            || !high.is_finite()
            || !low.is_finite()
        {
            return Err(error(
                OPERATION,
                GeometryErrorKind::InvalidResult,
                "Healed STEP material has no finite positive extent",
            ));
        }
        let frame = match decoded.placement() {
            None => PlaneFrame::default(),
            Some(p) => {
                let [u, v, _] = p.axes();
                PlaneFrame {
                    origin: p.origin_mm(),
                    u,
                    v,
                }
            }
        };
        frame.checked(policy)?;
        Ok(Some(Self {
            body: SolidBody::HealedImport(
                Rc::new(HealedImportBody {
                    feature,
                    volume,
                    radius_mm: radius,
                    heights_mm: [low, high],
                    sewing_mm: sewing,
                }),
                frame,
            ),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        }))
    }
}
