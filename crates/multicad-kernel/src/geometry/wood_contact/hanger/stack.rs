//! Bind the complete related-material proof to current Native source owners.
use super::*;
use source::{HangerStack, PlacedAnnulus, StackRequirements, StackSources, ThroughStock};
#[cfg(test)]
mod tests;
impl Solid {
    /// Prove intentional wood bite, every related material separation, actual
    /// cap bearing lands and whole-phase nut mating. Scene components outside
    /// this stack still require their own source-owned collision proof.
    pub fn rust_hanger_stack(
        &self,
        wood: &Self,
        clamp: &Self,
        nut: &Self,
        washers: &[&Self],
        requirements: StackRequirements,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HangerStack, GeometryError> {
        self.policy.validate()?;
        if self.policy.boolean_fuzzy_mm != 0.
            || [wood, clamp, nut]
                .into_iter()
                .chain(washers.iter().copied())
                .any(|solid| solid.policy != self.policy)
            || requirements.wood.linear_tolerance_mm != self.policy.linear_tolerance_mm
            || requirements.wood.minimum_feature_mm != self.policy.minimum_feature_mm
            || requirements.passage.linear_tolerance_mm != self.policy.linear_tolerance_mm
            || requirements.passage.minimum_feature_mm != self.policy.minimum_feature_mm
            || requirements.mating.linear_tolerance_mm != self.policy.linear_tolerance_mm
        {
            return Err(error(
                "Hanger material stack",
                GeometryErrorKind::Input,
                "Stack verification requires the actual source geometry policies",
            ));
        }
        if washers.len() > 62 {
            return Err(error(
                "Hanger material stack",
                GeometryErrorKind::Limit,
                "Hanger stack exceeds the 64 hardware member budget",
            ));
        }
        let unsupported = || {
            error(
                "Hanger material stack",
                GeometryErrorKind::Unsupported,
                "Stack verification requires an actual authored hanger, blind pilot, through bore, retained nut and annular washers",
            )
        };
        let hanger = self.rust_hanger_data().ok_or_else(unsupported)?;
        let (wood, wood_frame) = wood.orthogonal_body().ok_or_else(unsupported)?;
        let (material, pilot) = wood.wood_contact_source().ok_or_else(unsupported)?;
        let (clamp, clamp_frame) = clamp.orthogonal_body().ok_or_else(unsupported)?;
        let (clamp_material, bore) = clamp.hanger_through_source().ok_or_else(unsupported)?;
        let (nut, nut_frame) = nut.rust_hex_nut_geometry().ok_or_else(unsupported)?;
        for placed in [hanger.frame, wood_frame, clamp_frame, nut_frame] {
            placed.checked(self.policy)?;
        }
        let mut rings = Vec::with_capacity(washers.len());
        for washer in washers {
            let (ring, placed) = washer.rust_annulus_data().ok_or_else(unsupported)?;
            placed.checked(self.policy)?;
            rings.push((ring, placed));
        }
        let rings: Vec<_> = rings
            .iter()
            .map(|(ring, placed)| PlacedAnnulus {
                source: ring,
                frame: frame(*placed),
            })
            .collect();
        source::certify_hanger_stack(
            StackSources {
                wood: Stock {
                    material,
                    pilot,
                    frame: frame(wood_frame),
                },
                clamp: ThroughStock {
                    material: clamp_material,
                    bore,
                    frame: frame(clamp_frame),
                },
                hanger: &hanger.source,
                hanger_frame: frame(hanger.frame),
                washers: &rings,
                nut,
                nut_frame: frame(nut_frame),
            },
            requirements,
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
}
