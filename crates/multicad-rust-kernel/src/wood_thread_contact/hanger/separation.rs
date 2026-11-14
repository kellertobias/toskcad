//! Outer enclosures of actual retained component sources. Overlap means only
//! unresolved; these envelopes never become material or a Boolean exemption.
use super::super::separation as envelopes;
use super::*;
use crate::{annulus::Annulus, hex_nut::Nut, prism::Prism};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy)]
pub enum HangerComponent<'a> {
    Hanger(&'a Hanger),
    Nut(&'a Nut),
    Annulus(&'a Annulus),
    Stock(&'a OrthogonalSolid),
    Prism(&'a Prism),
}
#[derive(Clone, Copy)]
pub struct PlacedHangerComponent<'a> {
    pub source: HangerComponent<'a>,
    pub frame: Frame,
}
fn bounds(component: PlacedHangerComponent<'_>) -> R<Option<[I; 3]>> {
    let local = match component.source {
        HangerComponent::Hanger(source) => source.bounds_mm().map(|b| I::new(b[0], b[1])),
        HangerComponent::Nut(source) => {
            // AF contains both the actual finite hex vertices and the complete
            // regular hex support. The retained bore only removes material.
            let r = source.across_flats_mm();
            [
                I::new(-r, r),
                I::new(-r, r),
                I::new(0., source.thread().spec().length_mm),
            ]
        }
        HangerComponent::Annulus(source) => {
            let r = source.outer_radius_mm();
            [I::new(-r, r), I::new(-r, r), I::new(0., source.height_mm())]
        }
        HangerComponent::Stock(source) => {
            let b = source.bounds();
            std::array::from_fn(|i| I::new(b.min[i], b.max[i]))
        }
        HangerComponent::Prism(source) => {
            if source.rounded_end().is_some() {
                return Ok(None);
            }
            return envelopes::rectangular_box(source, component.frame).map(Some);
        }
    };
    envelopes::world_box(local, component.frame).map(Some)
}
/// Some(true) earns strict separation of the actual material. Some(false) and
/// None are both unresolved, never a positive collision or disjointness claim.
pub fn hanger_component_separation(
    a: PlacedHangerComponent<'_>,
    b: PlacedHangerComponent<'_>,
    interrupt: &dyn Fn() -> bool,
) -> R<Option<bool>> {
    checkpoint(interrupt)?;
    let a = bounds(a)?;
    let b = bounds(b)?;
    let result = match (a, b) {
        (Some(a), Some(b)) => Some(envelopes::separated(a, b)),
        _ => None,
    };
    checkpoint(interrupt)?;
    Ok(result)
}
