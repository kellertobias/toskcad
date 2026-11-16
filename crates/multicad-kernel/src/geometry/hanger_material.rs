//! Native placement of the distinct continuous authored Hanger source.
//! Solid topology and normal Model publication require separate admission.
use super::internal_thread_treatment_body::interval::Interval as I;
use super::*;
use multicad_rust_kernel::authored_hanger_bolt::{self as source, Hanger, QueryOptions};
mod faces;
mod mesh;
pub use faces::HangerFaceProperties;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HangerMaterial {
    source: Hanger,
    frame: PlaneFrame,
}

pub(super) fn failure(cause: source::Error) -> GeometryError {
    error(
        "Rust authored hanger",
        match cause {
            source::Error::Invalid => GeometryErrorKind::Input,
            source::Error::Unsupported | source::Error::Unresolved => {
                GeometryErrorKind::Unsupported
            }
            source::Error::Limit => GeometryErrorKind::Limit,
            source::Error::Cancelled => GeometryErrorKind::Cancelled,
        },
        format!("Actual connected hanger source proof failed: {cause:?}"),
    )
}

fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(failure(source::Error::Cancelled))
    } else {
        Ok(())
    }
}

impl HangerMaterial {
    pub fn new(
        source: Hanger,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        let spec = source.spec();
        policy.length(spec.plain_length_mm)?;
        policy.length(spec.plain_radius_mm)?;
        policy.length(source.metric().minimum_radius_mm()[0])?;
        policy.length(source.wood().minimum_radius_mm())?;
        let metric_shoulder =
            (I::point(spec.plain_radius_mm) - I::point(source.metric().maximum_radius_mm())).lo;
        let wood_shoulder =
            (I::point(source.wood().minimum_radius_mm()) - I::point(spec.plain_radius_mm)).lo;
        policy.length(metric_shoulder)?;
        policy.length(wood_shoulder)?;
        let local = source.bounds_mm();
        frame.checked_metric_extent(local.map(|b| b[0].abs().max(b[1].abs())), policy)?;
        let result = Self { source, frame };
        let bounds = result.bounds();
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        checkpoint(interrupt)?;
        Ok(result)
    }
    pub fn source(&self) -> &Hanger {
        &self.source
    }
    pub fn frame(&self) -> PlaneFrame {
        self.frame
    }
    /// Infallible after admission: source dimensions and frame origin are
    /// bounded by 1e7mm and admitted axes are near-unit. Each outward affine
    /// operation stays finite. World policy limits are checked by `new`.
    pub fn bounds(&self) -> Bounds {
        let local = self.source.bounds_mm().map(|[lo, hi]| I::new(lo, hi));
        let axes = [self.frame.u, self.frame.v, self.frame.normal()];
        let world: [I; 3] = std::array::from_fn(|i| {
            (0..3).fold(I::point(self.frame.origin[i]), |s, j| {
                s + I::point(axes[j][i]) * local[j]
            })
        });
        Bounds {
            min: world.map(|v| v.lo),
            max: world.map(|v| v.hi),
        }
    }
    /// Actual source volume under the literal affine frame. This does not
    /// substitute the manufacturing mesh's separately bounded signed volume.
    pub fn volume_bounds(
        &self,
        options: QueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<[f64; 2], GeometryError> {
        let v = self
            .source
            .volume(options, interrupt)
            .map_err(failure)?
            .material_mm3;
        let [u, vaxis, n] =
            [self.frame.u, self.frame.v, self.frame.normal()].map(|a| a.map(I::point));
        let c = [(1, 2), (2, 0), (0, 1)].map(|(j, k)| u[j] * vaxis[k] - u[k] * vaxis[j]);
        let determinant = (0..3).fold(I::ZERO, |s, i| s + c[i] * n[i]);
        let world = I::new(v[0], v[1]) * determinant;
        if !world.valid() || world.lo <= 0. || world.width() > options.volume_error_mm3 {
            return Err(failure(source::Error::Unresolved));
        }
        checkpoint(interrupt)?;
        Ok([world.lo, world.hi])
    }
}

#[cfg(test)]
mod face_tests;
#[cfg(test)]
mod mesh_tests;
#[cfg(test)]
pub(in crate::geometry) mod source_tests;
