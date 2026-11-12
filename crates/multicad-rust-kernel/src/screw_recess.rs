//! Typed physical recesses: shape and floor remain geometry, never size labels.
use crate::{
    hex_drive::HexDrive,
    hexalobular_drive::{Cone, Drive, Error},
    tapered_thread::interval::I,
};
use std::sync::Arc;
#[cfg(test)]
mod tests;

#[derive(Clone, Debug, PartialEq)]
pub enum Recess {
    FlatHex(FlatHex),
    ConicalHexalobular(ConicalHexalobular),
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlatHex {
    contour: HexDrive,
    depth_mm: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct ConicalHexalobular {
    source: Arc<Drive>,
    volume_bounds_mm3: [f64; 2],
}
impl FlatHex {
    pub fn geometry(&self) -> (&HexDrive, f64) {
        (&self.contour, self.depth_mm)
    }
}
impl ConicalHexalobular {
    pub fn source(&self) -> &Drive {
        &self.source
    }
}
impl Recess {
    pub fn canonical_flat_hex(across_flats_mm: f64, depth_mm: f64) -> Result<Self, Error> {
        let contour = HexDrive::canonical(
            across_flats_mm,
            64. * f64::EPSILON * across_flats_mm.abs().max(1.),
        )
        .map_err(|_| Error::Invalid)?;
        Self::flat_hex(contour, depth_mm)
    }
    pub fn flat_hex(contour: HexDrive, depth_mm: f64) -> Result<Self, Error> {
        if !depth_mm.is_finite() || depth_mm <= 0. {
            return Err(Error::Invalid);
        }
        if depth_mm > 1e7 {
            return Err(Error::Limit);
        }
        Ok(Self::FlatHex(FlatHex { contour, depth_mm }))
    }
    pub fn conical_hexalobular(
        source: Drive,
        width_mm3: f64,
        work_limit: usize,
    ) -> Result<Self, Error> {
        Self::conical_hexalobular_with_cancel(source, width_mm3, work_limit, &|| false)
    }
    pub fn conical_hexalobular_with_cancel(
        source: Drive,
        width_mm3: f64,
        work_limit: usize,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, Error> {
        let volume_bounds_mm3 =
            source.volume_bounds_with_cancel(width_mm3, work_limit, cancelled)?;
        Ok(Self::ConicalHexalobular(ConicalHexalobular {
            source: Arc::new(source),
            volume_bounds_mm3,
        }))
    }
    pub fn flat_hex_geometry(&self) -> Option<(&HexDrive, f64)> {
        match self {
            Self::FlatHex(h) => Some(h.geometry()),
            Self::ConicalHexalobular(_) => None,
        }
    }
    pub fn hexalobular_geometry(&self) -> Option<&Drive> {
        match self {
            Self::ConicalHexalobular(h) => Some(h.source()),
            Self::FlatHex(_) => None,
        }
    }
    pub fn source_error_mm(&self) -> f64 {
        match self {
            Self::FlatHex(_) => 0.,
            Self::ConicalHexalobular(h) => h.source.source_error_mm(),
        }
    }
    pub fn maximum_radius_bound_mm(&self) -> f64 {
        match self {
            Self::FlatHex(h) => h.contour.maximum_radius_mm(),
            Self::ConicalHexalobular(h) => h.source.maximum_radius_bound_mm(),
        }
    }
    pub fn maximum_depth_mm(&self) -> f64 {
        match self {
            Self::FlatHex(h) => h.depth_mm,
            Self::ConicalHexalobular(h) => h.source.cone().tip_depth_mm,
        }
    }
    pub fn wall_count(&self) -> u32 {
        match self {
            Self::FlatHex(_) => 6,
            Self::ConicalHexalobular(_) => 12,
        }
    }
    pub fn radius_at_angle_mm(&self, angle: f64) -> Option<f64> {
        if !angle.is_finite() {
            return None;
        }
        match self {
            Self::FlatHex(h) => Some(h.contour.radius_at_angle_mm(angle)),
            Self::ConicalHexalobular(h) => h.source.radius_at_angle_mm(angle),
        }
    }
    /// Nominal evaluation only. A conical curved recess uses the midpoint of
    /// its cached certificate; callers needing precision use the bounds.
    pub fn nominal_volume_mm3(&self) -> f64 {
        match self {
            Self::FlatHex(h) => h.contour.area_mm2() * h.depth_mm,
            Self::ConicalHexalobular(h) => {
                I::new(h.volume_bounds_mm3[0], h.volume_bounds_mm3[1]).midpoint()
            }
        }
    }
    pub fn volume_bounds_mm3(&self) -> [f64; 2] {
        match self {
            Self::FlatHex(h) => {
                let [lo, hi] = h
                    .contour
                    .area_bounds_mm2(0.)
                    .expect("validated finite hex geometry");
                let bounds = I::new(lo, hi) * I::point(h.depth_mm);
                [bounds.lo, bounds.hi]
            }
            Self::ConicalHexalobular(h) => h.volume_bounds_mm3,
        }
    }
    /// Classify the cut in recess coordinates, extending above its open mouth.
    /// Its volume certificate covers only nonnegative depths. Clearance shrinks the cut;
    /// `None` preserves a numerically unresolved circular ray intersection.
    pub fn contains_void(&self, point: [f64; 3], clearance_mm: f64) -> Option<bool> {
        if !point.into_iter().all(f64::is_finite) || !clearance_mm.is_finite() || clearance_mm < 0.
        {
            return None;
        }
        if point[2] >= self.maximum_depth_mm() - clearance_mm {
            return Some(false);
        }
        let [x, y, z] = point;
        let radius = x.hypot(y);
        let angle = y.atan2(x);
        let (rim, depth) = match self {
            Self::FlatHex(h) => (
                h.contour.radius_at_angle_mm(angle),
                h.depth_mm - clearance_mm,
            ),
            Self::ConicalHexalobular(h) => {
                let rim = h.source.radius_at_angle_mm(angle)?;
                let Cone {
                    tip_depth_mm,
                    depth_per_radius,
                } = h.source.cone();
                (
                    rim,
                    tip_depth_mm
                        - radius * depth_per_radius
                        - clearance_mm * (1. + depth_per_radius),
                )
            }
        };
        Some(radius < rim - clearance_mm && z < depth)
    }
}
