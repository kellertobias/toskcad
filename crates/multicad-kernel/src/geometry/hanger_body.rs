//! Retained connected Hanger geometry. Finite reference admission is separate.
use super::*;
use multicad_rust_kernel::authored_hanger_bolt::{Hanger, QueryOptions};
mod contained_boolean;
mod display;
mod step;
mod topology;
use multicad_rust_kernel::authored_hanger_bolt as source;
pub use step::HealedHangerStepVerification;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HangerData {
    pub source: Hanger,
    pub frame: PlaneFrame,
}
pub(crate) struct HangerBody {
    source: Hanger,
    pub(super) feature: EntityId,
    boundary: source::Boundary,
    qualified_step:
        std::cell::RefCell<Option<Rc<multicad_rust_kernel::step::hanger::QualifiedHangerStep>>>,
}
impl Solid {
    pub fn rust_hanger(
        source: Hanger,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        if feature.is_nil() {
            return Err(input("Hanger needs its original feature identity"));
        }
        HangerMaterial::new(source, frame, policy, interrupt)?;
        let boundary = source
            .boundary(source::BoundaryOptions::default(), interrupt)
            .map_err(super::hanger_material::failure)?;
        Ok(Self {
            body: SolidBody::Hanger(
                Rc::new(HangerBody {
                    source,
                    feature,
                    boundary,
                    qualified_step: std::cell::RefCell::new(None),
                }),
                frame,
            ),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn hanger_body(&self) -> Option<(&HangerBody, PlaneFrame)> {
        match &self.body {
            SolidBody::Hanger(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
    pub fn rust_hanger_data(&self) -> Option<HangerData> {
        self.hanger_body().map(|(body, frame)| HangerData {
            source: body.source,
            frame,
        })
    }
}
impl HangerBody {
    pub(super) fn source(&self) -> &Hanger {
        &self.source
    }
    pub(super) fn material(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<HangerMaterial, GeometryError> {
        HangerMaterial::new(self.source, frame, policy, None)
    }
    pub(super) fn volume(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<f64, GeometryError> {
        let scale = self
            .source
            .bounds_mm()
            .into_iter()
            .flatten()
            .map(f64::abs)
            .fold(1., f64::max);
        let bounds = self.material(frame, policy)?.volume_bounds(
            QueryOptions {
                volume_error_mm3: (policy.linear_tolerance_mm * scale * scale).max(1e-7),
                // One cell per actual wood band, one plain section and
                // five sections of the rounded metric profile.
                maximum_cells: self.source.wood().bands().len() + 6,
            },
            None,
        )?;
        Ok(bounds[0].midpoint(bounds[1]))
    }
}
