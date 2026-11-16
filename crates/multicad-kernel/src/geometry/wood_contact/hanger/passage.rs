//! Bind empty-space passage proofs to the actual retained Native materials.
use super::*;
use source::{HangerPassage, PassageRequirements, ThroughStock};
#[cfg(test)]
mod tests;
impl Solid {
    pub fn rust_hanger_through_clearance(
        &self,
        hanger: &Self,
        minimum_wall_mm: f64,
        requirements: PassageRequirements,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HangerPassage, GeometryError> {
        self.policy.validate()?;
        if self.policy != hanger.policy
            || self.policy.boolean_fuzzy_mm != 0.
            || requirements.linear_tolerance_mm != self.policy.linear_tolerance_mm
            || requirements.minimum_feature_mm != self.policy.minimum_feature_mm
        {
            return Err(error(
                "Hanger clamp clearance",
                GeometryErrorKind::Input,
                "Clearance requires the actual source geometry policies",
            ));
        }
        let unsupported = || {
            error(
                "Hanger clamp clearance",
                GeometryErrorKind::Unsupported,
                "Clearance requires one actual retained through bore and a connected authored hanger",
            )
        };
        let (body, placed) = self.orthogonal_body().ok_or_else(unsupported)?;
        let (material, bore) = body.hanger_through_source().ok_or_else(unsupported)?;
        let hanger = hanger.rust_hanger_data().ok_or_else(unsupported)?;
        placed.checked(self.policy)?;
        hanger.frame.checked(self.policy)?;
        source::certify_hanger_through_stock(
            ThroughStock {
                material,
                bore,
                frame: frame(placed),
            },
            &hanger.source,
            frame(hanger.frame),
            minimum_wall_mm,
            requirements,
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
    pub fn rust_hanger_annulus_clearance(
        &self,
        hanger: &Self,
        requirements: PassageRequirements,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HangerPassage, GeometryError> {
        self.policy.validate()?;
        if self.policy != hanger.policy
            || self.policy.boolean_fuzzy_mm != 0.
            || requirements.linear_tolerance_mm != self.policy.linear_tolerance_mm
            || requirements.minimum_feature_mm != self.policy.minimum_feature_mm
        {
            return Err(error(
                "Hanger washer clearance",
                GeometryErrorKind::Input,
                "Clearance requires the actual source geometry policies",
            ));
        }
        let unsupported = || {
            error(
                "Hanger washer clearance",
                GeometryErrorKind::Unsupported,
                "Clearance requires an actual retained annulus and a connected authored hanger",
            )
        };
        let (ring, placed) = self.rust_annulus_data().ok_or_else(unsupported)?;
        let hanger = hanger.rust_hanger_data().ok_or_else(unsupported)?;
        placed.checked(self.policy)?;
        hanger.frame.checked(self.policy)?;
        source::certify_hanger_annulus(
            &ring,
            frame(placed),
            &hanger.source,
            frame(hanger.frame),
            requirements,
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
}
