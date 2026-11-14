//! Actual whole screw inside the drilled pilot or free entry halfspace.
use super::*;
use source::{MatingExternal, ScrewClearance, ScrewClearanceRequirements};
impl Solid {
    pub fn rust_wood_insert_screw_clearance(
        &self,
        screw: &Self,
        minimum_tip_clearance_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ScrewClearance, GeometryError> {
        let operation = "Wood insert screw blind-bottom clearance";
        self.policy.validate()?;
        if self.policy != screw.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(error(
                operation,
                GeometryErrorKind::Input,
                "Clearance requires equal actual source policies without Boolean fuzz",
            ));
        }
        let unsupported = || {
            error(
                operation,
                GeometryErrorKind::Unsupported,
                "Clearance requires actual retained screw and blind pilot material",
            )
        };
        let (body, placed) = self.orthogonal_body().ok_or_else(unsupported)?;
        let (material, pilot) = body.wood_contact_source().ok_or_else(unsupported)?;
        let (source_screw, screw_frame) = mating::actual_screw(screw).ok_or_else(unsupported)?;
        placed.checked(self.policy)?;
        screw_frame.checked(self.policy)?;
        source::certify_insert_screw_clearance(
            Stock {
                material,
                pilot,
                frame: frame(placed),
            },
            MatingExternal {
                source: source_screw,
                frame: frame(screw_frame),
            },
            ScrewClearanceRequirements {
                minimum_tip_clearance_mm,
                linear_tolerance_mm: self.policy.linear_tolerance_mm,
                minimum_feature_mm: self.policy.minimum_feature_mm,
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(|cause| {
            let mut e = source_error(cause);
            e.operation = operation;
            e
        })
    }
}
#[cfg(test)]
mod tests;
