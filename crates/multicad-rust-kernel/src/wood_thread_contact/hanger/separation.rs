//! Outer enclosures of actual retained component sources. Overlap means only
//! unresolved; these envelopes never become material or a Boolean exemption.
use super::super::separation as envelopes;
use super::*;
use crate::{
    annulus::Annulus, cylinder::AxialCylinder, hex_nut::Nut, prism::Prism, sphere::Sphere,
};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy)]
pub enum HangerComponent<'a> {
    Hanger(&'a Hanger),
    Nut(&'a Nut),
    Annulus(&'a Annulus),
    Stock(&'a OrthogonalSolid),
    Prism(&'a Prism),
    /// Placement is the first rim centre, matching the native sweep datum.
    Cylinder(&'a AxialCylinder),
    Sphere(&'a Sphere),
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
        HangerComponent::Cylinder(source) => {
            let exact = |r: crate::Rational| {
                let (n, d) = r.ratio();
                enclosure(&Q::new(n.into(), d.into()))
            };
            let radius = positive(exact(source.circle().radius())?)?.hi;
            let height = positive(exact(source.height())?)?.hi;
            [
                I::new(-radius, radius),
                I::new(-radius, radius),
                I::new(0., height),
            ]
        }
        HangerComponent::Sphere(source) => {
            let radius = source.radius_mm();
            let local = I::new(-radius, radius);
            let mut world = envelopes::world_box([local; 3], component.frame)?;
            // The analytic sphere is rotation invariant. Also enclose the raw
            // frame image used by retained surface/mesh consumers, without
            // assuming that admitted floating axes are exactly orthonormal.
            for (i, coordinate) in world.iter_mut().enumerate() {
                let analytic = finite(ip(component.frame.origin_mm[i]) + local)?;
                coordinate.lo = coordinate.lo.min(analytic.lo);
                coordinate.hi = coordinate.hi.max(analytic.hi);
            }
            return Ok(Some(world));
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
