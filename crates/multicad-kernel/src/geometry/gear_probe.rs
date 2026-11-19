//! Whole-prism relations against actual Rust involute gear solids.
//!
//! A rectangular or polygon prism wholly inside or wholly outside the gear
//! material reuses one unchanged operand, exactly as for the other
//! containment providers. Each gear is split into axial slabs whose section
//! is constant (spur bands, hubs, internal rings, gap cores) or rigidly
//! rotated by a Lipschitz phase (helical and chevron bands). The probe's
//! enclosing chart box becomes one transverse disc per slab, settled against
//! the exact section by the chord-tube proof of `classify_*_disc`. A probe
//! reaching any cap, shoulder, flank, bore or rim stays unresolved.
use super::*;
use crate::geometry::helical_body::contained_boolean::{
    chart_bounds, probe_bounds, whole_operand_result,
};
use multicad_rust_kernel::{
    helical_thread::BoxClassification,
    spur_gear::{
        Profile,
        clearance::{self, DiscSide},
        herringbone::Herringbone,
        internal::Ring,
    },
};

#[derive(Clone, Copy)]
enum Section<'a> {
    Spur(&'a Profile),
    Twisted(&'a Herringbone),
    Annulus([f64; 2]),
    Ring(&'a Ring),
}

struct Slab<'a> {
    z: [f64; 2],
    section: Section<'a>,
}

fn bore_radius(profile: &Profile) -> f64 {
    profile
        .bore()
        .and_then(|edge| edge.radius_mm())
        .unwrap_or(0.)
}

fn hub_slabs<'a>(
    extensions: [Option<multicad_rust_kernel::spur_gear::HubExtension>; 2],
    bore: f64,
    slabs: &mut Vec<Slab<'a>>,
) {
    for hub in extensions.into_iter().flatten() {
        let (lo, hi) = (
            hub.gear_z_mm.min(hub.end_z_mm),
            hub.gear_z_mm.max(hub.end_z_mm),
        );
        slabs.push(Slab {
            z: [lo, hi],
            section: Section::Annulus([bore, hub.radius_mm]),
        });
    }
}

fn slabs(solid: &Solid) -> Option<(Vec<Slab<'_>>, PlaneFrame)> {
    let mut slabs = Vec::new();
    let frame = if let Some((body, frame)) = solid.spur_gear_body() {
        let profile = &body.profile;
        slabs.push(Slab {
            z: [0., profile.height_mm()],
            section: Section::Spur(profile),
        });
        hub_slabs(profile.hub_extensions(), bore_radius(profile), &mut slabs);
        frame
    } else if let Some((body, frame)) = solid.herringbone_body() {
        let gear = &body.geometry;
        let profile = gear.profile();
        let height = profile.height_mm();
        let bore = bore_radius(profile);
        match gear.shoulder_heights_mm() {
            Some([low, high]) => {
                slabs.push(Slab {
                    z: [0., low],
                    section: Section::Twisted(gear),
                });
                slabs.push(Slab {
                    z: [low, high],
                    section: Section::Annulus([bore, gear.root_radius_mm()]),
                });
                slabs.push(Slab {
                    z: [high, height],
                    section: Section::Twisted(gear),
                });
            }
            None => slabs.push(Slab {
                z: [0., height],
                section: Section::Twisted(gear),
            }),
        }
        hub_slabs(gear.hub_extensions(), bore, &mut slabs);
        frame
    } else {
        let (body, frame) = solid.internal_ring_body()?;
        slabs.push(Slab {
            z: [0., body.ring.height_mm()],
            section: Section::Ring(&body.ring),
        });
        frame
    };
    Some((slabs, frame))
}

/// The local probe box `[x, y]` as a disc about its center.
fn disc([x, y]: [[f64; 2]; 2]) -> ([f64; 2], f64) {
    let center = [(x[0] + x[1]) / 2., (y[0] + y[1]) / 2.];
    let radius = ((x[1] - x[0]) / 2.).hypot((y[1] - y[0]) / 2.);
    // Enclose the rounding of the midpoint and the half diagonal.
    let guard = 8. * f64::EPSILON * (center[0].abs() + center[1].abs() + radius).max(1.);
    (center, radius + guard)
}

fn section_side(
    section: Section<'_>,
    [x, y, z]: [[f64; 2]; 3],
    tolerance: f64,
) -> Result<DiscSide, clearance::Error> {
    let (center, radius) = disc([x, y]);
    match section {
        Section::Spur(profile) => clearance::classify_spur_disc(profile, center, radius, tolerance),
        Section::Ring(ring) => clearance::classify_ring_disc(ring, center, radius, tolerance),
        Section::Annulus(radii) => clearance::classify_annulus_disc(radii, center, radius),
        Section::Twisted(gear) => {
            // Material at height z is the profile rotated by phase(z), and
            // the phase is |twist|-Lipschitz in z (continuous at a chevron
            // apex). Unrotating by the mid phase moves any probe point by at
            // most its radius times half the phase range.
            let middle = (z[0] + z[1]) / 2.;
            let swing = gear.twist_rad_per_mm().abs() * (z[1] - z[0]) / 2.;
            let (sin, cos) = (-gear.phase_at(middle)).sin_cos();
            let rotated = [
                cos * center[0] - sin * center[1],
                sin * center[0] + cos * center[1],
            ];
            let reach = center[0].hypot(center[1]) + radius;
            let radius = radius * (1. + 4. * f64::EPSILON) + reach * swing * (1. + 1e-12);
            clearance::classify_spur_disc(gear.profile(), rotated, radius, tolerance)
        }
    }
}

/// Classify the local probe chart box against every axial slab.
fn classify(
    slabs: &[Slab<'_>],
    bounds: [[f64; 2]; 3],
    tolerance: f64,
) -> Result<BoxClassification, clearance::Error> {
    use BoxClassification::*;
    if bounds.iter().flatten().any(|v| !v.is_finite()) || bounds.iter().any(|[lo, hi]| lo > hi) {
        return Ok(UnresolvedBoundary);
    }
    let z = bounds[2];
    let scale = slabs
        .iter()
        .flat_map(|slab| slab.z)
        .chain(z)
        .map(f64::abs)
        .fold(1., f64::max);
    let guard = 64. * f64::EPSILON * scale;
    let mut relation = Outside;
    for slab in slabs {
        let [lo, hi] = slab.z;
        if z[1] < lo - guard || z[0] > hi + guard {
            continue;
        }
        if z[0] <= lo + guard || z[1] >= hi - guard {
            // The probe reaches this slab's cap or shoulder plane. Only a
            // probe missing the material on both sides is still settled.
            let clipped = [z[0].max(lo), z[1].min(hi)];
            if section_side(slab.section, [bounds[0], bounds[1], clipped], tolerance)?
                != DiscSide::Outside
            {
                return Ok(UnresolvedBoundary);
            }
            continue;
        }
        match section_side(slab.section, bounds, tolerance)? {
            DiscSide::Inside => relation = Inside,
            DiscSide::Outside => {}
            DiscSide::Unresolved => return Ok(UnresolvedBoundary),
        }
    }
    Ok(relation)
}

impl Solid {
    /// A finite rectangular or polygon prism wholly inside or wholly outside
    /// an actual Rust spur, helical, chevron or internal-ring gear reuses one
    /// unchanged operand.
    pub(in crate::geometry) fn contained_gear_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        let (probe, source_is_left, (slabs, frame)) = if let Some(gear) = slabs(self) {
            (other, true, gear)
        } else {
            (self, false, slabs(other)?)
        };
        let bounds = chart_bounds(
            probe_bounds(probe)?,
            frame.origin,
            [frame.u, frame.v, frame.normal()],
        );
        let relation = match classify(&slabs, bounds, policy.linear_tolerance_mm) {
            Ok(relation) => relation,
            Err(clearance::Error::Cancelled) => return None,
            Err(clearance::Error::Invalid) => return None,
            Err(clearance::Error::Limit) => {
                return Some(Err(error(
                    "Boolean",
                    GeometryErrorKind::Limit,
                    "Rust gear containment exceeded its chord budget",
                )));
            }
        };
        whole_operand_result(self, other, source_is_left, relation, operation)
    }
}

#[cfg(test)]
mod tests;
