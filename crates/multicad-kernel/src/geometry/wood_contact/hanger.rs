//! Actual connected hanger source queries. These independent certificates do
//! not waive clamp, washer, stock-interface or other component collisions.
use super::*;
use source::{HangerContact, HangerMating, HangerMatingRequirements};
#[cfg(test)]
mod tests;
impl Solid {
    /// Intentional wood bite against the actual current blind-pilot stock.
    pub fn rust_hanger_wood_contact(
        &self,
        hanger: &Self,
        requirements: Requirements,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HangerContact, GeometryError> {
        self.policy.validate()?;
        if self.policy != hanger.policy
            || self.policy.boolean_fuzzy_mm != 0.
            || requirements.linear_tolerance_mm != self.policy.linear_tolerance_mm
            || requirements.minimum_feature_mm != self.policy.minimum_feature_mm
        {
            return Err(error(
                "Hanger wood contact",
                GeometryErrorKind::Input,
                "Contact requires the actual source geometry policies",
            ));
        }
        let unsupported = || {
            error(
                "Hanger wood contact",
                GeometryErrorKind::Unsupported,
                "Contact requires a connected authored hanger and one retained orthogonal blind pilot",
            )
        };
        let (stock, stock_frame) = self.orthogonal_body().ok_or_else(unsupported)?;
        let (material, pilot) = stock.wood_contact_source().ok_or_else(unsupported)?;
        let hanger = hanger.rust_hanger_data().ok_or_else(unsupported)?;
        stock_frame.checked(self.policy)?;
        hanger.frame.checked(self.policy)?;
        source::certify_hanger_wood(
            Stock {
                material,
                pilot,
                frame: frame(stock_frame),
            },
            &hanger.source,
            frame(hanger.frame),
            requirements,
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
    /// Whole-phase clearance of the actual continuous metric source against
    /// an actual finite retained hex nut, including arbitrary relative roll.
    pub fn rust_hanger_nut_mating(
        &self,
        nut: &Self,
        requirements: HangerMatingRequirements,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HangerMating, GeometryError> {
        self.policy.validate()?;
        if self.policy != nut.policy
            || self.policy.boolean_fuzzy_mm != 0.
            || requirements.linear_tolerance_mm != self.policy.linear_tolerance_mm
        {
            return Err(error(
                "Hanger nut mating",
                GeometryErrorKind::Input,
                "Mating requires the actual source geometry policies",
            ));
        }
        let unsupported = || {
            error(
                "Hanger nut mating",
                GeometryErrorKind::Unsupported,
                "Mating requires a connected authored hanger and an actual retained hex nut",
            )
        };
        let hanger = self.rust_hanger_data().ok_or_else(unsupported)?;
        let (nut, nut_frame) = nut.rust_hex_nut_geometry().ok_or_else(unsupported)?;
        hanger.frame.checked(self.policy)?;
        nut_frame.checked(self.policy)?;
        source::certify_hanger_nut(
            &hanger.source,
            frame(hanger.frame),
            nut,
            frame(nut_frame),
            requirements,
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
}
