//! Immutable source owners for strict component envelope separation.
use super::*;
use multicad_rust_kernel::{
    annulus::Annulus, authored_hanger_bolt::Hanger, hex_nut::Nut, orthogonal::OrthogonalSolid,
    prism::Prism,
};
use source::{HangerComponent, PlacedHangerComponent};
#[cfg(test)]
mod tests;
enum Actual<'a> {
    Hanger(&'a Hanger, PlaneFrame),
    Nut(&'a Nut, PlaneFrame),
    Ring(Annulus, PlaneFrame),
    Stock(&'a OrthogonalSolid, PlaneFrame),
    Prism(&'a Prism, PlaneFrame),
}
impl<'a> Actual<'a> {
    fn read(solid: &'a Solid) -> Option<Self> {
        if let Some((body, frame)) = solid.hanger_body() {
            Some(Self::Hanger(body.source(), frame))
        } else if let Some((source, frame)) = solid.rust_hex_nut_geometry() {
            Some(Self::Nut(source, frame))
        } else if let Some((source, frame)) = solid.rust_annulus_data() {
            Some(Self::Ring(source, frame))
        } else if let Some((source, frame)) = solid.orthogonal_body() {
            source
                .wood_envelope_source()
                .map(|source| Self::Stock(source, frame))
        } else {
            solid
                .rust_prism_frame()
                .map(|(source, frame)| Self::Prism(source, frame))
        }
    }
    fn view(&self) -> PlacedHangerComponent<'_> {
        let (source, placed) = match self {
            Self::Hanger(source, placed) => (HangerComponent::Hanger(source), placed),
            Self::Nut(source, placed) => (HangerComponent::Nut(source), placed),
            Self::Ring(source, placed) => (HangerComponent::Annulus(source), placed),
            Self::Stock(source, placed) => (HangerComponent::Stock(source), placed),
            Self::Prism(source, placed) => (HangerComponent::Prism(source), placed),
        };
        PlacedHangerComponent {
            source,
            frame: frame(*placed),
        }
    }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(source_error(source::Error::Cancelled))
    } else {
        Ok(())
    }
}
impl Solid {
    /// Only Some(true) proves actual source material separation. Other results
    /// are unresolved; no nominal Solid bounds or Boolean fallback is used.
    pub fn rust_hanger_component_separation(
        &self,
        other: &Self,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<bool>, GeometryError> {
        checkpoint(interrupt)?;
        self.policy.validate()?;
        if self.policy != other.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(error(
                "Hanger scene clearance",
                GeometryErrorKind::Input,
                "Clearance requires the actual source geometry policies",
            ));
        }
        let result = match (Actual::read(self), Actual::read(other)) {
            (Some(a), Some(b)) => source::hanger_component_separation(a.view(), b.view(), &|| {
                interrupt.is_some_and(|f| f())
            })
            .map_err(source_error)?,
            _ => None,
        };
        checkpoint(interrupt)?;
        Ok(result)
    }
}
