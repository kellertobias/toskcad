//! Connected canonical wood/plain/continuous rounded metric material.
//! This source is distinct from the legacy literal hanger; it does not imply
//! a mesh, Native topology or exchange representation.
use crate::{
    authored_flat_thread::AuthoredThread, authored_metric_thread::AuthoredMetricThread,
    tapered_thread::interval::I,
};
mod boundary;
mod classify;
pub use boundary::{
    Boundary, BoundaryOptions, Coedge, Component, Curve, CurveBounds, CurveKey, End, Knot,
    VertexKey,
};
mod faces;
pub use faces::FaceKey;
mod properties;
#[cfg(test)]
pub(crate) mod tests;
mod volume;
pub use properties::{FaceProperties, FacePropertyOptions};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Unresolved,
    Limit,
    Cancelled,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spec {
    pub wood: AuthoredThread,
    pub metric: AuthoredMetricThread,
    pub plain_length_mm: f64,
    pub plain_radius_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SectionKind {
    Metric,
    Plain,
    Wood,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxialSection {
    pub kind: SectionKind,
    pub start_mm: [f64; 2],
    pub end_mm: [f64; 2],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hanger {
    spec: Spec,
    metric_bottom: [f64; 2],
    core_radius_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Classification {
    Inside,
    Outside,
    UnresolvedBoundary,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QueryOptions {
    pub volume_error_mm3: f64,
    /// Five actual bands per thread and the real plain cylinder consume11cells.
    pub maximum_cells: usize,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Volume {
    pub material_mm3: [f64; 2],
    pub wood_mm3: [f64; 2],
    pub plain_mm3: [f64; 2],
    pub metric_mm3: [f64; 2],
}
impl Hanger {
    pub fn new(spec: Spec, interrupt: Option<&dyn Fn() -> bool>) -> Result<Self, Error> {
        checkpoint(interrupt)?;
        if [spec.plain_length_mm, spec.plain_radius_mm]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0. || v > 1e7)
        {
            return Err(Error::Invalid);
        }
        // These are exact comparisons of admitted finite binary values. No
        // tolerance adjusts a radius or turns a touching shoulder into a land.
        if spec.plain_radius_mm <= spec.metric.maximum_radius_mm()
            || spec.plain_radius_mm >= spec.wood.minimum_radius_mm()
        {
            return Err(Error::Unsupported);
        }
        let bottom = -I::point(spec.plain_length_mm) - I::point(spec.metric.full_length_mm());
        // Each dimension is <=1e7: the directed sum is finite without a second
        // independently selected physical station or a renamed local origin.
        let result = Self {
            metric_bottom: pair(bottom),
            core_radius_mm: spec.metric.minimum_radius_mm()[0]
                .min(spec.plain_radius_mm)
                .min(spec.wood.minimum_radius_mm()),
            spec,
        };
        checkpoint(interrupt)?;
        Ok(result)
    }
    pub fn spec(&self) -> Spec {
        self.spec
    }
    pub fn wood(&self) -> &AuthoredThread {
        &self.spec.wood
    }
    pub fn metric(&self) -> &AuthoredMetricThread {
        &self.spec.metric
    }
    pub fn connection_core_radius_mm(&self) -> f64 {
        self.core_radius_mm
    }
    pub fn sections(&self) -> [AxialSection; 3] {
        let p = -self.spec.plain_length_mm;
        [
            AxialSection {
                kind: SectionKind::Metric,
                start_mm: self.metric_bottom,
                end_mm: [p, p],
            },
            AxialSection {
                kind: SectionKind::Plain,
                start_mm: [p, p],
                end_mm: [0., 0.],
            },
            AxialSection {
                kind: SectionKind::Wood,
                start_mm: [0., 0.],
                end_mm: [self.wood().full_length_mm(); 2],
            },
        ]
    }
    /// Complete source support, for separation rather than material membership.
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        let r = self.wood().maximum_radius_mm();
        [
            [-r, r],
            [-r, r],
            [self.metric_bottom[0], self.wood().full_length_mm()],
        ]
    }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn pair(x: I) -> [f64; 2] {
    [x.lo, x.hi]
}
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}

// Detached closed candidate; Native manufacturing must still qualify the actual
// placed vertices and global literal embedding before publication.
pub use crate::authored_internal_treatment::HangerCandidate as Candidate;
impl Hanger {
    pub fn tessellate_candidate(
        &self,
        options: crate::display::Tessellation,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Candidate, Error> {
        crate::authored_internal_treatment::hanger_candidate(self, options, maximum_work, interrupt)
    }
}
