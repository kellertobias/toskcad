//! Positive clearance or a strict common-interior point for an external
//! spur pinion against one exact, finite straight-sided outline such as a
//! rack. Like the gear certificate, a result applies to one fixed pose.
use super::{
    Error, MAX_PAIR_TESTS, MAX_REFINEMENTS, Placement, ResultKind, Rim, boundaries_clear,
    checkpoint, tip_radius,
};
use crate::{
    normal_sweep,
    revolution::dot,
    spur_gear::{Edge, Profile, add, scale, sub},
};

/// Certify a pinion (`placements[0]`) against a closed counterclockwise
/// outline (`placements[1]`) whose local +Y points toward the pinion and
/// whose local Y=0 is its pitch line. The direction only steers the
/// common-interior candidate; Clear is proved for any outline position.
pub fn certify_rack(
    pinion: &Profile,
    rack: &[normal_sweep::Edge],
    placements: [Placement; 2],
    tolerance_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<ResultKind, Error> {
    if !(1e-7..=1e-3).contains(&tolerance_mm)
        || placements.into_iter().any(|p| !p.valid(1e-3))
        || rack.len() < 3
        || rack
            .iter()
            .any(|edge| matches!(edge, normal_sweep::Edge::Circle { .. }))
    {
        return Err(Error::Invalid);
    }
    // Only a candidate, midway between the rack tips and the pinion tips
    // on the perpendicular from the pinion center to the pitch line. The
    // winding tests certify membership in both outlines.
    let [along, toward] = placements[1].axes;
    let origin = placements[1].center_mm;
    let center = sub(placements[0].center_mm, origin);
    let rack_tip = rack
        .iter()
        .map(|edge| edge.support([0., 1.]))
        .fold(f64::NEG_INFINITY, f64::max);
    let height = (rack_tip + dot(center, toward) - tip_radius(pinion)) * 0.5;
    let witness = add(
        origin,
        add(scale(along, dot(center, along)), scale(toward, height)),
    );
    let mut deflection = (pinion.root_radius_mm * 2e-4).clamp(tolerance_mm, 0.01);
    for _ in 0..MAX_REFINEMENTS {
        checkpoint(interrupt)?;
        let teeth = Rim::build(pinion, placements[0], deflection, interrupt)?;
        let outline = Rim::from_outline(
            rack.iter().map(|&edge| Edge::Planar(edge)),
            0.,
            placements[1],
            deflection,
            interrupt,
        )?;
        if boundaries_clear(&teeth, &outline, MAX_PAIR_TESTS, interrupt)? {
            return Ok(ResultKind::Clear);
        }
        if teeth.contains_strictly(witness)
            && outline.winding(witness, outline.tube_mm + outline.guard_mm) == Some(true)
        {
            return Ok(ResultKind::CommonInterior);
        }
        deflection *= 0.25;
    }
    Ok(ResultKind::Unresolved)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        spur_gear::{Spec, clearance::nearest_boundary},
        straight_rack,
    };
    use std::f64::consts::{FRAC_PI_2, PI};

    const TOL: f64 = 1e-7;

    /// A module-1, 20° pinion with the same thinning as the rack.
    fn pinion(teeth: u32, thinning: f64) -> Profile {
        let pitch = f64::from(teeth) / 2.;
        let pressure = 20_f64.to_radians();
        let base = pitch * pressure.cos();
        let thickness = PI / 2. - thinning;
        Profile::new(
            Spec {
                teeth,
                pitch_radius_mm: pitch,
                base_radius_mm: base,
                tip_radius_mm: pitch + 1.,
                root_radius_mm: pitch - 1.25,
                root_fillet_radius_mm: 0.2,
                half_base_angle_rad: thickness / (2. * pitch) + pressure.tan() - pressure,
                face_width_mm: 2.,
                bore_radius_mm: 0.,
            },
            1e-6,
            1e7,
        )
        .unwrap()
    }

    fn rack(thinning: f64) -> Vec<normal_sweep::Edge> {
        straight_rack::outline(
            straight_rack::Spec {
                teeth: 20,
                pitch_mm: PI,
                pitch_thickness_mm: PI / 2. - thinning,
                pressure_angle_rad: 20_f64.to_radians(),
                addendum_mm: 1.,
                dedendum_mm: 1.25,
                root_fillet_radius_mm: 0.2,
                backing_thickness_mm: 3.,
                face_width_mm: 2.,
            },
            1e-6,
            1e7,
        )
        .unwrap()
        .0
    }

    /// The pinion turns by `drive` while the rack translates by the
    /// pitch-line arc. With 20 cells the rack centre is a tooth space, so
    /// a pinion tooth centred on the -Y axis is the nominal mesh.
    fn pose(teeth: u32, drive: f64) -> [Placement; 2] {
        let pitch = f64::from(teeth) / 2.;
        let rotation = -FRAC_PI_2 + drive;
        let (sin, cos) = rotation.sin_cos();
        [
            Placement {
                center_mm: [0., pitch],
                axes: [[cos, sin], [-sin, cos]],
                uncertainty_mm: 0.,
            },
            Placement {
                center_mm: [pitch * drive, 0.],
                ..Placement::default()
            },
        ]
    }

    #[test]
    fn meshed_rack_samples_clear_with_actual_flank_proximity_and_wrong_phase_shares_material() {
        let teeth = 24;
        let gear = pinion(teeth, 0.05);
        let bar = rack(0.05);
        let edges: Vec<_> = gear
            .outline()
            .iter()
            .enumerate()
            .map(|(i, item)| (i, item.edge))
            .collect();
        let rack_edges: Vec<_> = bar
            .iter()
            .enumerate()
            .map(|(i, &edge)| (i, Edge::Planar(edge)))
            .collect();
        let half_period = PI / f64::from(teeth);
        for i in 0..9 {
            let drive = half_period * (f64::from(i) / 4. - 1.);
            let placements = pose(teeth, drive);
            assert_eq!(
                certify_rack(&gear, &bar, placements, TOL, None),
                Ok(ResultKind::Clear),
                "{i}"
            );
            let witness = nearest_boundary([&edges, &rack_edges], placements, 0.11, None)
                .unwrap()
                .expect("flank proximity within total thinning");
            assert!(
                witness.distance_mm > 0. && witness.distance_mm <= 0.11,
                "{witness:?}"
            );
        }
        // Half a pitch out of phase puts a pinion tooth onto a rack tooth.
        let wrong = pose(teeth, 0.);
        let shifted = [
            wrong[0],
            Placement {
                center_mm: [PI / 2., 0.],
                ..wrong[1]
            },
        ];
        assert_eq!(
            certify_rack(&gear, &bar, shifted, TOL, None),
            Ok(ResultKind::CommonInterior)
        );
        // Withdrawn far enough, the pair still clears but has no witness.
        let mut apart = pose(teeth, 0.);
        apart[0].center_mm[1] += 1.;
        assert_eq!(
            certify_rack(&gear, &bar, apart, TOL, None),
            Ok(ResultKind::Clear)
        );
        assert_eq!(
            nearest_boundary([&edges, &rack_edges], apart, 0.11, None),
            Ok(None)
        );
    }

    #[test]
    fn invalid_inputs_cancellation_and_unresolved_contact_fail_closed() {
        let gear = pinion(24, 0.05);
        let bar = rack(0.05);
        let placements = pose(24, 0.);
        assert_eq!(
            certify_rack(&gear, &bar, placements, 0.1, None),
            Err(Error::Invalid)
        );
        assert_eq!(
            certify_rack(&gear, &bar[..2], placements, TOL, None),
            Err(Error::Invalid)
        );
        let circle = [normal_sweep::Edge::Circle {
            center: [0.; 2],
            radius_mm: 1.,
            clockwise: false,
        }; 3];
        assert_eq!(
            certify_rack(&gear, &circle, placements, TOL, None),
            Err(Error::Invalid)
        );
        let mut skewed = placements;
        skewed[1].axes[0] = [2., 0.];
        assert_eq!(
            certify_rack(&gear, &bar, skewed, TOL, None),
            Err(Error::Invalid)
        );
        assert_eq!(
            certify_rack(&gear, &bar, placements, TOL, Some(&|| true)),
            Err(Error::Cancelled)
        );
        // Unthinned teeth touch at the pitch point: neither a clearance
        // nor an interior witness can be certified.
        let tight = pinion(24, 0.);
        assert_eq!(
            certify_rack(&tight, &rack(0.), placements, TOL, None),
            Ok(ResultKind::Unresolved)
        );
    }
}
