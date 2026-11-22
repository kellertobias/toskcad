//! Whole actual screw material in the retained cavity or free entry halfspace.
use super::*;
use source::{InsertMating, MatingExternal, MatingRequirements, MatingScrew};
impl Solid {
    /// Actual head/plain/neck stays above entry; every overlapping part of the
    /// unchanged occupied thread clears the actual finite insert cavity.
    /// Thread continuation beyond the insert requires the separate actual
    /// stock.rust_wood_insert_screw_clearance proof; mating never waives stock.
    pub fn rust_wood_insert_mating(
        &self,
        insert: &Self,
        minimum_engagement_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<InsertMating, GeometryError> {
        self.policy.validate()?;
        let fail = |kind, message| error("Wood insert mating", kind, message);
        if self.policy != insert.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(fail(
                GeometryErrorKind::Input,
                "Mating requires equal actual source policies without Boolean fuzz",
            ));
        }
        let unsupported = || {
            fail(
                GeometryErrorKind::Unsupported,
                "Mating requires an actual retained screw and paired authored insert",
            )
        };
        let insert = insert.rust_wood_insert_data().ok_or_else(unsupported)?;
        let (screw, placed) = actual_screw(self).ok_or_else(unsupported)?;
        insert.frame.checked(self.policy)?;
        placed.checked(self.policy)?;
        source::certify_insert_mating(
            &insert.source,
            frame(insert.frame),
            MatingExternal {
                source: screw,
                frame: frame(placed),
            },
            MatingRequirements {
                minimum_engagement_mm,
                minimum_clearance_mm: 4. * self.policy.linear_tolerance_mm,
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
}
pub(super) fn actual_screw(solid: &Solid) -> Option<(MatingScrew<'_>, PlaneFrame)> {
    if let Some((s, f)) = solid.rust_machine_screw_geometry() {
        Some((MatingScrew::Regular(s), f))
    } else if let Some((s, f)) = solid.rust_profiled_screw_geometry() {
        Some((MatingScrew::Profiled(s), f))
    } else if let Some((s, f)) = solid.rust_periodic_repaired_phillips_source() {
        Some((MatingScrew::Periodic(s), f))
    } else {
        solid
            .rust_phillips_screw_geometry()
            .map(|(s, f)| (MatingScrew::Phillips(s), f))
    }
}
#[cfg(test)]
pub(super) mod tests;
