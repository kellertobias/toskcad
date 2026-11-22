//! Closed axial interval separation for a fixed continuous-chevron pose.
use super::*;
use crate::spur_gear::herringbone::Herringbone;

const MAX_INTERVALS: usize = 4096;
const MAX_TWIST_UNCERTAINTY_MM: f64 = 0.04;

/// Certify an entire fixed herringbone pose over its aligned face width.
/// Callers supply an orthonormal common transverse chart and a bound on
/// any projection/placement roundoff. The two apex planes must coincide.
pub fn certify_herringbone(
    gears: [&Herringbone; 2],
    placements: [Placement; 2],
    tolerance_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<ResultKind, Error> {
    let profiles = gears.map(Herringbone::profile);
    // Added hubs stay strictly inside the root disk, which belongs to
    // every transverse outer profile. A Clear result excludes intersecting
    // outer Jordan regions and nesting, so these coaxial disks remain
    // disjoint even beyond the aligned tooth bands. Interior witnesses
    // continue to refer to real tooth material with the bore excluded.
    // Gapped chevrons with the same gap repeat their lower band in the
    // upper band; between the shoulders each keeps only its root disk,
    // which a Clear tooth band already separates. Other gap mixtures
    // would need the core-only section and are refused.
    let height = profiles[0].height_mm();
    if height != profiles[1].height_mm()
        || gears[0].gap_mm() != gears[1].gap_mm()
        || placements.into_iter().any(|p| !p.valid(1e-3))
        || !(1e-7..=1e-3).contains(&tolerance_mm)
    {
        return Err(Error::Invalid);
    }

    // Two chevrons repeat their lower sections in exact reverse order.
    // A single-hand sweep requires the full width. The absolute twist
    // rate bounds displacement even across a chevron's middle plane.
    let width = if gears.iter().any(|gear| gear.is_helical()) {
        height
    } else {
        gears[0].band_length_mm()
    };
    let mut pending = vec![[0., width]];
    let mut remaining = MAX_INTERVALS;
    while let Some([start, end]) = pending.pop() {
        checkpoint(interrupt)?;
        remaining = remaining.checked_sub(1).ok_or(Error::Limit)?;
        let half = (end - start) / 2.;
        let midpoint = start + half;
        let placed = placed_at(gears, placements, midpoint, half);
        if placed
            .iter()
            .all(|p| p.uncertainty_mm <= MAX_TWIST_UNCERTAINTY_MM)
        {
            match certify_with_uncertainty(
                profiles,
                placed,
                tolerance_mm,
                MAX_TWIST_UNCERTAINTY_MM,
                interrupt,
            )? {
                ResultKind::Clear => continue,
                ResultKind::CommonInterior => return Ok(ResultKind::CommonInterior),
                ResultKind::Unresolved => {}
            }
        }
        if half <= tolerance_mm {
            return Ok(ResultKind::Unresolved);
        }
        pending.extend([[start, midpoint], [midpoint, end]]);
    }
    Ok(ResultKind::Clear)
}

/// Both transverse sections at height `z_mm`, turned by their own phases.
/// Each uncertainty also encloses every section within `half_mm` of
/// `z_mm`: a circle displacement is no greater than its arc length.
pub(super) fn placed_at(
    gears: [&Herringbone; 2],
    placements: [Placement; 2],
    z_mm: f64,
    half_mm: f64,
) -> [Placement; 2] {
    std::array::from_fn(|i| {
        let base = placements[i];
        let [u, v] = base.axes;
        let (sin, cos) = gears[i].phase_at(z_mm).sin_cos();
        let radial_speed = gears[i].tip_radius_mm() * gears[i].twist_rad_per_mm().abs();
        Placement {
            center_mm: base.center_mm,
            axes: [
                add(scale(u, cos), scale(v, sin)),
                add(scale(u, -sin), scale(v, cos)),
            ],
            uncertainty_mm: base.uncertainty_mm + radial_speed * half_mm,
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spur_gear::herringbone::tests::fixture;

    #[test]
    fn coaxial_hubs_preserve_clearance_outside_the_aligned_tooth_bands() {
        use crate::spur_gear::Hub;
        let phase = std::f64::consts::PI * (1. + 1. / 36.);
        let (sin, cos) = phase.sin_cos();
        let center = 60. / 15_f64.to_radians().cos();
        let placements = [
            Placement::default(),
            Placement {
                center_mm: [center, 0.],
                axes: [[cos, sin], [-sin, cos]],
                uncertainty_mm: 0.,
            },
        ];
        for helical in [false, true] {
            let gears = [24, 36].map(|teeth| {
                let (spec, twist) = fixture(teeth);
                let gear = if helical {
                    Herringbone::new_helical(spec, twist, 1e-5, 1e7)
                } else {
                    Herringbone::new(spec, twist, 1e-5, 1e7)
                }
                .unwrap();
                gear.with_hub(
                    Hub {
                        radius_mm: 8.,
                        axial_start_mm: -10.,
                        axial_end_mm: 20.,
                    },
                    1e-5,
                    1e7,
                )
                .unwrap()
            });
            assert_eq!(
                certify_herringbone([&gears[0], &gears[1]], placements, 1e-7, None),
                Ok(ResultKind::Clear)
            );
            assert!(gears[0].axial_bounds()[0] < 0. && gears[1].axial_bounds()[1] > 4.);
            assert!(center > 2. * 8.);
        }
    }

    #[test]
    fn single_hand_and_mixed_sweeps_check_material_past_the_middle_plane() {
        let (spec_a, twist_a) = fixture(24);
        let (spec_b, twist_b) = fixture(36);
        let center = spec_a.pitch_radius_mm + spec_b.pitch_radius_mm;
        let nominal = std::f64::consts::PI * (1. + 1. / 36.);
        let at_phase = |phase: f64| {
            let (sin, cos) = phase.sin_cos();
            [
                Placement::default(),
                Placement {
                    center_mm: [center, 0.],
                    axes: [[cos, sin], [-sin, cos]],
                    uncertainty_mm: 0.,
                },
            ]
        };
        let tolerance = 1e-7;

        let linear_a = Herringbone::new_helical(spec_a, twist_a, 1e-5, 1e7).unwrap();
        let linear_b = Herringbone::new_helical(spec_b, twist_b, 1e-5, 1e7).unwrap();
        let placement = at_phase(nominal);
        assert_eq!(
            certify_herringbone([&linear_a, &linear_b], placement, tolerance, None),
            Ok(ResultKind::Clear)
        );
        assert_eq!(
            certify_herringbone(
                [&linear_a, &linear_b],
                at_phase(std::f64::consts::PI),
                tolerance,
                None
            ),
            Ok(ResultKind::CommonInterior)
        );

        let chevron_b = Herringbone::new(spec_b, twist_b, 1e-5, 1e7).unwrap();
        let reflected_a = Herringbone::new(spec_a, twist_a, 1e-5, 1e7).unwrap();
        assert_eq!(
            certify_herringbone([&reflected_a, &chevron_b], placement, tolerance, None),
            Ok(ResultKind::Clear)
        );
        // The pitch witness lies in B's gap, but the pair collides: a tooth
        // tip lies strictly in both real outlines. The independent point
        // below confirms the collision at the far cap.
        assert_eq!(
            certify_herringbone([&linear_a, &chevron_b], placement, tolerance, None),
            Ok(ResultKind::CommonInterior)
        );
        let (sin, cos) = linear_a.phase_at(spec_a.face_width_mm).sin_cos();
        let end_a = Placement {
            axes: [[cos, sin], [-sin, cos]],
            ..placement[0]
        };
        let witness = [spec_a.pitch_radius_mm, 1.8];
        assert!(
            Rim::build(linear_a.profile(), end_a, 1e-3, None)
                .unwrap()
                .contains_strictly(witness)
        );
        assert!(
            Rim::build(chevron_b.profile(), placement[1], 1e-3, None)
                .unwrap()
                .contains_strictly(witness)
        );
    }

    #[test]
    fn a_wide_helical_pair_turned_half_a_tooth_has_common_interior_off_the_centre_line() {
        // KHK normal-system 20/30 pair, mn = 2, beta = 20 deg, 8 mm wide:
        // the drive phase advances by almost half a pitch over the face,
        // so at many heights the pitch witness lies in both gaps while
        // neighbouring tips still collide.
        let beta = 20_f64.to_radians();
        let pressure = (20_f64.to_radians().tan() / beta.cos()).atan();
        let gear = |teeth: u32, hand: f64| {
            let pitch = 2. * f64::from(teeth) / (2. * beta.cos());
            let spec = crate::spur_gear::Spec {
                teeth,
                pitch_radius_mm: pitch,
                base_radius_mm: pitch * pressure.cos(),
                tip_radius_mm: pitch + 2.,
                root_radius_mm: pitch - 2.5,
                root_fillet_radius_mm: 0.4,
                half_base_angle_rad: std::f64::consts::PI / (2. * f64::from(teeth))
                    + pressure.tan()
                    - pressure
                    - 0.05 / (2. * pitch),
                face_width_mm: 8.,
                bore_radius_mm: 2.,
            };
            let twist = hand * beta.tan() / pitch;
            (
                Herringbone::new_helical(spec, twist, 1e-5, 1e7).unwrap(),
                pitch,
            )
        };
        let ((a, ra), (b, rb)) = (gear(20, 1.), gear(30, -1.));
        let at_phase = |phase: f64| {
            let (sin, cos) = phase.sin_cos();
            [
                Placement::default(),
                Placement {
                    center_mm: [ra + rb, 0.],
                    axes: [[cos, sin], [-sin, cos]],
                    uncertainty_mm: 0.,
                },
            ]
        };
        use std::f64::consts::PI;
        assert_eq!(
            certify_herringbone([&a, &b], at_phase(PI + PI / 30.), 1e-7, None),
            Ok(ResultKind::Clear)
        );
        assert_eq!(
            certify_herringbone([&a, &b], at_phase(PI), 1e-7, None),
            Ok(ResultKind::CommonInterior)
        );
    }

    #[test]
    fn full_width_clearance_separates_contact_from_material_overlap() {
        let [a, b] = [24, 36].map(|teeth| {
            let (spec, twist) = fixture(teeth);
            Herringbone::new(spec, twist, 1e-5, 1e7).unwrap()
        });
        let pair = [&a, &b];
        let placement = |distance: f64, phase: f64| {
            let (sin, cos) = phase.sin_cos();
            [
                Placement::default(),
                Placement {
                    center_mm: [distance, 0.],
                    axes: [[cos, sin], [-sin, cos]],
                    uncertainty_mm: 0.,
                },
            ]
        };
        let center = 60. / 15_f64.to_radians().cos();
        let tolerance = 1e-7;
        let clear = placement(center, std::f64::consts::PI * (1. + 1. / 36.));
        let overlaps = placement(center, std::f64::consts::PI);
        // At the end cap the two tips are exactly tangent. The certified
        // clearance proof must not promote this to a strict separation.
        let touch = placement(a.tip_radius_mm() + b.tip_radius_mm(), std::f64::consts::PI);
        assert_eq!(
            certify_herringbone(pair, clear, tolerance, None),
            Ok(ResultKind::Clear)
        );
        assert_eq!(
            certify_herringbone(pair, overlaps, tolerance, None),
            Ok(ResultKind::CommonInterior)
        );
        assert_eq!(
            certify_herringbone(pair, touch, tolerance, None),
            Ok(ResultKind::Unresolved)
        );
        assert_eq!(
            certify_herringbone(pair, clear, tolerance, Some(&|| true)),
            Err(Error::Cancelled)
        );
        assert_eq!(
            certify_herringbone(pair, [Placement::default(); 2], tolerance, None),
            Err(Error::Invalid)
        );
        assert_eq!(
            certify_herringbone(pair, clear, f64::NAN, None),
            Err(Error::Invalid)
        );
    }

    #[test]
    fn equal_gap_double_helical_pairs_certify_their_bands_and_refuse_mixtures() {
        let [(spec_a, twist_a), (spec_b, twist_b)] = [24, 36].map(|teeth| {
            let (spec, twist) = fixture(teeth);
            (
                crate::spur_gear::Spec {
                    face_width_mm: 10.,
                    ..spec
                },
                twist,
            )
        });
        let center = spec_a.pitch_radius_mm + spec_b.pitch_radius_mm;
        let at_phase = |phase: f64| {
            let (sin, cos) = phase.sin_cos();
            [
                Placement::default(),
                Placement {
                    center_mm: [center, 0.],
                    axes: [[cos, sin], [-sin, cos]],
                    uncertainty_mm: 0.,
                },
            ]
        };
        let nominal = at_phase(std::f64::consts::PI * (1. + 1. / 36.));
        let gapped =
            |spec, twist, gap| Herringbone::new_gapped(spec, twist, gap, 1e-5, 1e7).unwrap();
        let a = gapped(spec_a, twist_a, 3.);
        let b = gapped(spec_b, twist_b, 3.);
        assert_eq!(
            certify_herringbone([&a, &b], nominal, 1e-7, None),
            Ok(ResultKind::Clear)
        );
        // Turning B by half a tooth puts its teeth into A's teeth.
        assert_eq!(
            certify_herringbone([&a, &b], at_phase(std::f64::consts::PI), 1e-7, None),
            Ok(ResultKind::CommonInterior)
        );
        let narrower = gapped(spec_b, twist_b, 2.);
        let continuous = Herringbone::new(spec_b, twist_b, 1e-5, 1e7).unwrap();
        let helical = Herringbone::new_helical(spec_b, twist_b, 1e-5, 1e7).unwrap();
        for other in [&narrower, &continuous, &helical] {
            assert_eq!(
                certify_herringbone([&a, other], nominal, 1e-7, None),
                Err(Error::Invalid)
            );
        }
    }
}
