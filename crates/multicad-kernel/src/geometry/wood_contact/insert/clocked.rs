//! Separate earned raw-clock providers; existing quarter-clock queries stay first.
use super::*;
impl Solid {
    pub fn rust_clocked_wood_insert_mating(
        &self,
        insert: &Self,
        minimum: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<source::ClockedInsertMating, GeometryError> {
        self.policy.validate()?;
        if self.policy != insert.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(input(
                "Clocked mating requires equal actual source policies",
            ));
        }
        let unsupported = || {
            error(
                "Clocked insert mating",
                GeometryErrorKind::Unsupported,
                "Clocked mating requires actual regular screw and original paired insert",
            )
        };
        let (screw, placed) = self.rust_machine_screw_geometry().ok_or_else(unsupported)?;
        let insert = insert.rust_wood_insert_data().ok_or_else(unsupported)?;
        source::certify_clocked_insert_mating(
            &insert.source,
            frame(insert.frame),
            screw,
            frame(placed),
            source::ClockedMatingOptions {
                requirements: source::MatingRequirements {
                    minimum_engagement_mm: minimum,
                    minimum_clearance_mm: 4. * self.policy.linear_tolerance_mm,
                },
                maximum_cells: 4096,
                linear_tolerance_mm: self.policy.linear_tolerance_mm,
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
    pub fn rust_clocked_wood_insert_screw_clearance(
        &self,
        screw: &Self,
        minimum: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<source::ScrewClearance, GeometryError> {
        self.policy.validate()?;
        if self.policy != screw.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(input(
                "Clocked clearance requires equal actual source policies",
            ));
        }
        let unsupported = || {
            error(
                "Clocked insert floor clearance",
                GeometryErrorKind::Unsupported,
                "Clocked clearance requires original regular screw and actual blind pilot",
            )
        };
        let (stock, placed) = self.orthogonal_body().ok_or_else(unsupported)?;
        let (material, pilot) = stock.wood_contact_source().ok_or_else(unsupported)?;
        let (screw, screw_frame) = screw
            .rust_machine_screw_geometry()
            .ok_or_else(unsupported)?;
        source::certify_clocked_insert_screw_clearance(
            Stock {
                material,
                pilot,
                frame: frame(placed),
            },
            screw,
            frame(screw_frame),
            source::ScrewClearanceRequirements {
                minimum_tip_clearance_mm: minimum,
                linear_tolerance_mm: self.policy.linear_tolerance_mm,
                minimum_feature_mm: self.policy.minimum_feature_mm,
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
}
#[cfg(test)]
mod diagnostic;
