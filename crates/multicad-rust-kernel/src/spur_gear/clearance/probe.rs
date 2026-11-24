//! Whole-disc material classification against one exact gear section.
//!
//! A disc is settled only when no chord tube of the exact outline comes
//! within the disc radius of its center. Every exact boundary point then
//! lies farther than the radius from the center, so the whole disc lies on
//! the side the center's ray parity reports. Circles (bores, hubs, rims and
//! gap cores) are compared exactly with a rounding guard. Anything closer is
//! `Unresolved`; nothing is approximated or reconstructed.
use super::{Error, MAX_REFINEMENTS, Placement, Rim};
use crate::{
    revolution::length,
    spur_gear::{Edge, Point, Profile, internal::Ring},
};
use std::sync::{Arc, Mutex};

/// The chord-tube rims of one fixed outline already built for disc
/// classification, by the exact bits of their deflection.
///
/// A rim is a pure function of the outline and its deflection, and each
/// classification refines through the same deflections, so a remembered
/// rim is the one a fresh build returns. It is derived data: a copy starts
/// empty and it never affects equality.
#[derive(Default)]
pub(crate) struct DiscRims(Mutex<Vec<(u64, Arc<Rim>)>>);

/// Rims kept per outline: every refinement of a few tolerances.
const DISC_RIMS: usize = 4 * MAX_REFINEMENTS;

impl Clone for DiscRims {
    fn clone(&self) -> Self {
        Self::default()
    }
}

impl PartialEq for DiscRims {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

impl std::fmt::Debug for DiscRims {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("DiscRims")
    }
}

impl DiscRims {
    /// The rim at `deflection`, built by `build` unless remembered. A
    /// poisoned lock only loses the memo.
    fn get(
        &self,
        deflection: f64,
        build: impl FnOnce() -> Result<Rim, Error>,
    ) -> Result<Arc<Rim>, Error> {
        let key = deflection.to_bits();
        let remembered = self.0.lock().ok().and_then(|rims| {
            rims.iter()
                .find(|(bits, _)| *bits == key)
                .map(|(_, rim)| Arc::clone(rim))
        });
        if let Some(rim) = remembered {
            return Ok(rim);
        }
        let rim = Arc::new(build()?);
        if let Ok(mut rims) = self.0.lock()
            && rims.len() < DISC_RIMS
        {
            rims.push((key, Arc::clone(&rim)));
        }
        Ok(rim)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DiscSide {
    /// Every disc point lies strictly in material.
    Inside,
    /// No disc point lies in material.
    Outside,
    Unresolved,
}

fn valid(center: Point, radius_mm: f64, tolerance_mm: f64) -> bool {
    center
        .into_iter()
        .all(|v| v.is_finite() && v.abs() <= super::super::MAX_COORDINATE_MM)
        && radius_mm.is_finite()
        && radius_mm >= 0.
        && (1e-7..=1e-3).contains(&tolerance_mm)
}

/// Rounding guard for exact circle comparisons at this scale.
fn guard(center: Point, radius_mm: f64, circle_mm: f64) -> f64 {
    1024. * f64::EPSILON * (length(center) + radius_mm + circle_mm).max(1.)
}

/// Whether the disc lies strictly inside (`Some(true)`) or strictly outside
/// (`Some(false)`) the closed exact outline, refining its chord tubes.
fn enclosed_by(
    outline: impl Fn() -> Vec<Edge>,
    rims: &DiscRims,
    center: Point,
    radius_mm: f64,
    tolerance_mm: f64,
    size_mm: f64,
) -> Result<Option<bool>, Error> {
    let mut deflection = (size_mm * 2e-4).clamp(tolerance_mm, 0.01);
    for _ in 0..MAX_REFINEMENTS {
        let rim = rims.get(deflection, || {
            Rim::from_outline(outline(), 0., Placement::default(), deflection, None)
        })?;
        // The tube encloses the exact edge, so a chord beyond this margin
        // proves the exact boundary lies beyond the disc.
        let margin = rim.tube_mm + rim.guard_mm + radius_mm;
        if let Some(inside) = rim.winding(center, margin) {
            return Ok(Some(inside));
        }
        deflection *= 0.25;
    }
    Ok(None)
}

/// The disc against the exact circle of `circle_mm` about the origin:
/// `Some(true)` strictly within, `Some(false)` strictly beyond.
fn within_circle(center: Point, radius_mm: f64, circle_mm: f64) -> Option<bool> {
    let distance = length(center);
    let guard = guard(center, radius_mm, circle_mm);
    if distance + radius_mm < circle_mm - guard {
        Some(true)
    } else if distance - radius_mm > circle_mm + guard {
        Some(false)
    } else {
        None
    }
}

/// Classify a disc in the transverse section of an external spur band:
/// the exact tooth outline minus its optional bore.
pub fn classify_spur_disc(
    profile: &Profile,
    center: Point,
    radius_mm: f64,
    tolerance_mm: f64,
) -> Result<DiscSide, Error> {
    if !valid(center, radius_mm, tolerance_mm) {
        return Err(Error::Invalid);
    }
    let bore = profile
        .bore_radius_mm()
        .map_or(Some(false), |bore| within_circle(center, radius_mm, bore));
    if bore == Some(true) {
        return Ok(DiscSide::Outside);
    }
    let outline = || profile.outline().iter().map(|item| item.edge).collect();
    Ok(
        match enclosed_by(
            outline,
            &profile.disc_rims,
            center,
            radius_mm,
            tolerance_mm,
            profile.root_radius_mm,
        )? {
            Some(false) => DiscSide::Outside,
            Some(true) if bore == Some(false) => DiscSide::Inside,
            _ => DiscSide::Unresolved,
        },
    )
}

/// Classify a disc in the section of an internal ring: inside the outside
/// rim and outside the exact inward-tooth cavity.
pub fn classify_ring_disc(
    ring: &Ring,
    center: Point,
    radius_mm: f64,
    tolerance_mm: f64,
) -> Result<DiscSide, Error> {
    if !valid(center, radius_mm, tolerance_mm) {
        return Err(Error::Invalid);
    }
    let rim = within_circle(center, radius_mm, ring.outer_radius_mm());
    if rim == Some(false) {
        return Ok(DiscSide::Outside);
    }
    let cavity = || ring.cavity().iter().map(|item| item.edge).collect();
    Ok(
        match enclosed_by(
            cavity,
            &ring.disc_rims,
            center,
            radius_mm,
            tolerance_mm,
            ring.tip_radius_mm(),
        )? {
            Some(true) => DiscSide::Outside,
            Some(false) if rim == Some(true) => DiscSide::Inside,
            _ => DiscSide::Unresolved,
        },
    )
}

/// Classify a disc against the exact annulus `inner < r < outer` (a hub or
/// a gap core around its bore; `inner` zero for none).
pub fn classify_annulus_disc(
    [inner, outer]: [f64; 2],
    center: Point,
    radius_mm: f64,
) -> Result<DiscSide, Error> {
    if !valid(center, radius_mm, 1e-7) || !(0. <= inner && inner < outer && outer.is_finite()) {
        return Err(Error::Invalid);
    }
    let beyond_inner = inner == 0. || within_circle(center, radius_mm, inner) == Some(false);
    let in_inner = inner > 0. && within_circle(center, radius_mm, inner) == Some(true);
    Ok(match within_circle(center, radius_mm, outer) {
        Some(false) => DiscSide::Outside,
        _ if in_inner => DiscSide::Outside,
        Some(true) if beyond_inner => DiscSide::Inside,
        _ => DiscSide::Unresolved,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spur_gear::Spec;
    use multicad_precision::PortableMath as _;
    use std::f64::consts::PI;

    fn spur(bore: f64) -> Profile {
        let pressure = 20_f64.to_radians();
        let pitch = 20.;
        Profile::new(
            Spec {
                teeth: 20,
                pitch_radius_mm: pitch,
                base_radius_mm: pitch * pressure.portable_cos(),
                tip_radius_mm: pitch + 2.,
                root_radius_mm: pitch - 2.5,
                root_fillet_radius_mm: 0.4,
                half_base_angle_rad: PI / 40. + pressure.portable_tan() - pressure,
                face_width_mm: 8.,
                bore_radius_mm: bore,
            },
            1e-6,
            1e7,
        )
        .unwrap()
    }

    fn polar(r: f64, a: f64) -> Point {
        [r * a.portable_cos(), r * a.portable_sin()]
    }

    #[test]
    fn tooth_space_bore_and_root_discs_are_settled_and_flank_discs_are_not() {
        let gear = spur(3.);
        let side = |c, r| classify_spur_disc(&gear, c, r, 1e-7).unwrap();
        // The tooth centered on +x is material at the pitch circle; the
        // space half a pitch away is not.
        assert_eq!(side(polar(20., 0.), 0.01), DiscSide::Inside);
        assert_eq!(side(polar(20., PI / 20.), 0.01), DiscSide::Outside);
        assert_eq!(side(polar(10., 0.3), 0.01), DiscSide::Inside);
        assert_eq!(side([0., 0.], 1.), DiscSide::Outside);
        assert_eq!(side([0., 0.], 3.5), DiscSide::Unresolved);
        assert_eq!(side(polar(3., 1.), 0.01), DiscSide::Unresolved);
        assert_eq!(side(polar(30., 0.), 1.), DiscSide::Outside);
        // A disc spanning a flank is never reported on either side.
        assert_eq!(side(polar(20., PI / 40.), 0.2), DiscSide::Unresolved);
        assert_eq!(side(polar(20., 0.), 2.), DiscSide::Unresolved);
    }

    #[test]
    fn ring_material_lies_between_the_cavity_and_the_outside_rim() {
        let a = 20_f64.to_radians();
        let pitch = 24.;
        let ring = Ring::new(
            crate::spur_gear::internal::RingSpec {
                teeth: 48,
                pitch_radius_mm: pitch,
                base_radius_mm: pitch * a.portable_cos(),
                tip_radius_mm: pitch - 1.,
                root_radius_mm: pitch + 1.25,
                root_fillet_radius_mm: 0.2,
                half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) - a.portable_tan() + a,
                face_width_mm: 2.,
                outer_radius_mm: pitch + 4.,
            },
            1e-5,
            1e7,
            1e-7,
        )
        .unwrap();
        let side = |c, r| classify_ring_disc(&ring, c, r, 1e-7).unwrap();
        // An inward tooth is centered on +x; its space is half a pitch away.
        assert_eq!(side(polar(pitch, 0.), 0.01), DiscSide::Inside);
        assert_eq!(side(polar(pitch, PI / 48.), 0.01), DiscSide::Outside);
        assert_eq!(side(polar(pitch + 2.5, 0.7), 0.01), DiscSide::Inside);
        assert_eq!(side([0., 0.], 5.), DiscSide::Outside);
        assert_eq!(side(polar(pitch + 4., 0.2), 0.01), DiscSide::Unresolved);
        assert_eq!(side(polar(pitch + 6., 0.2), 1.), DiscSide::Outside);
    }

    #[test]
    fn remembered_rims_classify_as_fresh_ones_and_stay_bounded() {
        let gear = spur(3.);
        // Discs on and near the flanks need every refinement; tolerances
        // above the start deflection add their own levels.
        let discs: Vec<(Point, f64)> = (0..240)
            .map(|k| {
                let k = f64::from(k);
                (
                    polar(17. + 0.0291 * k, 0.0137 * k),
                    [0.3, 0.01, 1e-4][k as usize % 3],
                )
            })
            .collect();
        for tolerance in [1e-7, 1e-4, 1e-3, 1e-7] {
            for &(center, radius) in &discs {
                // A copy starts with no rims and builds every one afresh.
                let fresh = classify_spur_disc(&gear.clone(), center, radius, tolerance);
                assert_eq!(classify_spur_disc(&gear, center, radius, tolerance), fresh);
            }
        }
        let kept = gear.disc_rims.0.lock().unwrap().len();
        assert!((1..=DISC_RIMS).contains(&kept), "{kept}");
        assert!(gear.clone().disc_rims.0.lock().unwrap().is_empty());
        assert_eq!(gear.clone(), gear);
    }

    #[test]
    fn annuli_compare_exact_circles_and_invalid_requests_fail_closed() {
        let side = |c, r| classify_annulus_disc([2., 5.], c, r).unwrap();
        assert_eq!(side([3.5, 0.], 1.), DiscSide::Inside);
        assert_eq!(side([0., 0.], 1.), DiscSide::Outside);
        assert_eq!(side([7., 0.], 1.), DiscSide::Outside);
        assert_eq!(side([5., 0.], 0.1), DiscSide::Unresolved);
        assert_eq!(side([2., 0.], 0.1), DiscSide::Unresolved);
        assert_eq!(
            classify_annulus_disc([0., 5.], [0., 0.], 1.).unwrap(),
            DiscSide::Inside
        );
        assert_eq!(
            classify_annulus_disc([5., 2.], [0., 0.], 1.),
            Err(Error::Invalid)
        );
        assert_eq!(
            classify_spur_disc(&spur(0.), [f64::NAN, 0.], 1., 1e-7),
            Err(Error::Invalid)
        );
        assert_eq!(
            classify_spur_disc(&spur(0.), [0., 0.], -1., 1e-7),
            Err(Error::Invalid)
        );
        assert_eq!(
            classify_spur_disc(&spur(0.), [0., 0.], 1., 1.),
            Err(Error::Invalid)
        );
    }
}
