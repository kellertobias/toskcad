//! Rust fixed-pose certificates and actual flank witnesses for hubbed and
//! profile-shifted spur pairs, finite racks with a spur pinion and twisted
//! bands. Every positive case also has a wrong-phase or separated control.
use multicad_core::document::EntityId;
use multicad_kernel::{
    GearHub, GearSweep, GeometryErrorKind, GeometryPolicy, InvoluteGear, RigidTransform, Solid,
    StraightRack, rust_kernel::spur_gear::clearance::ResultKind,
};
use multicad_precision::PortableMath as _;
use std::f64::consts::{FRAC_PI_2, PI};

const THINNING: f64 = 0.05;

/// A module-2, 20° external gear with `shift` profile-shift coefficient
/// and the tooth thinned by THINNING on the pitch circle.
fn spur(teeth: u32, shift: f64, hub: Option<GearHub>) -> InvoluteGear {
    let module = 2.;
    let pressure = 20_f64.to_radians();
    let pitch = module * f64::from(teeth) / 2.;
    let thickness = module * (PI / 2. + 2. * shift * pressure.portable_tan()) - THINNING;
    InvoluteGear {
        teeth,
        pitch_radius_mm: pitch,
        base_radius_mm: pitch * pressure.portable_cos(),
        tip_radius_mm: pitch + module * (1. + shift),
        root_radius_mm: pitch - module * (1.25 - shift),
        root_fillet_radius_mm: 0.4,
        half_base_angle_rad: thickness / (2. * pitch) + pressure.portable_tan() - pressure,
        face_width_mm: 4.,
        bore_radius_mm: 3.,
        hub,
        twist_rad_per_mm: 0.,
        sweep: GearSweep::Spur,
    }
}

fn rotated(solid: &Solid, phase: f64, translation_mm: [f64; 3]) -> Solid {
    let (s, c) = phase.portable_sin_cos();
    solid
        .transformed(RigidTransform {
            rotation: [[c, -s, 0.], [s, c, 0.], [0., 0., 1.]],
            translation_mm,
        })
        .unwrap()
}

/// Working pressure angle and centre distance of a shifted external pair.
fn working_center(teeth: [u32; 2], shifts: [f64; 2]) -> f64 {
    let pressure = 20_f64.to_radians();
    let inv = |a: f64| a.portable_tan() - a;
    let target = inv(pressure)
        + 2. * pressure.portable_tan() * (shifts[0] + shifts[1]) / f64::from(teeth[0] + teeth[1]);
    let mut working = pressure;
    for _ in 0..60 {
        working -= (inv(working) - target) / working.portable_tan().powi(2);
    }
    f64::from(teeth[0] + teeth[1]) * pressure.portable_cos() / working.portable_cos()
}

#[test]
fn hubbed_shifted_spurs_clear_with_named_flank_witnesses_and_wrong_phase_overlaps() {
    let policy = GeometryPolicy::default();
    let hub = |radius_mm| {
        Some(GearHub {
            radius_mm,
            axial_start_mm: -3.,
            axial_end_mm: 9.,
        })
    };
    let shifts = [0.4, -0.2];
    let a =
        Solid::rust_spur_gear_named(spur(12, shifts[0], hub(6.)), EntityId::new(), policy).unwrap();
    let b = Solid::rust_spur_gear_named(spur(30, shifts[1], hub(12.)), EntityId::new(), policy)
        .unwrap();
    // Working backlash shrinks the centre distance slightly below the
    // zero-backlash value, never into the thinned teeth.
    let center = working_center([12, 30], shifts);
    let nominal = PI + PI / 30.;
    let reach = 2. * THINNING + 10. * policy.linear_tolerance_mm;
    for i in 0..9 {
        let drive = 2. * PI / 12. * f64::from(i) / 8.;
        let first = rotated(&a, drive, [0.; 3]);
        let second = rotated(&b, nominal - drive * 12. / 30., [center, 0., 0.]);
        assert_eq!(
            first.certify_rust_spur_pair(&second, policy, None).unwrap(),
            ResultKind::Clear,
            "{i}"
        );
        let bands = first
            .rust_gear_flank_witnesses(&second, policy, reach, None)
            .unwrap();
        assert_eq!(bands.len(), 1);
        let witness = bands[0].expect("actual flank proximity");
        assert!(witness.distance_mm > 0. && witness.distance_mm <= reach);
        let gap = (0..3)
            .map(|k| (witness.points_mm[0][k] - witness.points_mm[1][k]).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!((gap - witness.distance_mm).abs() < 1e-12);
        // Both points lie in the common tooth band, on involute walls.
        assert!(witness.points_mm.iter().all(|p| (0. ..=4.).contains(&p[2])));
        let faces = witness.faces.unwrap();
        for (solid, face) in [(&first, faces[0]), (&second, faces[1])] {
            assert!(face >= 2);
            let named = solid
                .faces()
                .unwrap()
                .into_iter()
                .find(|f| f.index == face)
                .unwrap();
            assert!(named.plane.is_none(), "an involute wall is not planar");
        }
    }
    let opposed = rotated(&b, PI, [center, 0., 0.]);
    assert_eq!(
        a.certify_rust_spur_pair(&opposed, policy, None).unwrap(),
        ResultKind::CommonInterior
    );
    // Withdrawn by a millimetre, the pair is clear but has no engagement.
    let apart = rotated(&b, nominal, [center + 1., 0., 0.]);
    assert_eq!(
        a.certify_rust_spur_pair(&apart, policy, None).unwrap(),
        ResultKind::Clear
    );
    assert_eq!(
        a.rust_gear_flank_witnesses(&apart, policy, reach, None)
            .unwrap(),
        vec![None]
    );
}

#[test]
fn hub_extents_bound_axial_drift_and_witness_inputs_fail_closed() {
    let policy = GeometryPolicy::default();
    let long = Some(GearHub {
        radius_mm: 6.,
        axial_start_mm: -400.,
        axial_end_mm: 400.,
    });
    let a = Solid::rust_spur_gear_named(spur(12, 0., long), EntityId::new(), policy).unwrap();
    let b = Solid::rust_spur_gear_named(spur(30, 0., None), EntityId::new(), policy).unwrap();
    // A tilt the 4 mm tooth band alone could absorb exceeds policy over
    // the 800 mm hub: 1e-9 rad drifts 4e-9 mm across the band, 4e-7 mm here.
    let tilt = 1e-9_f64;
    let leaning = rotated(&b, PI + PI / 30., [0.; 3])
        .transformed(RigidTransform {
            rotation: [[1., 0., tilt], [0., 1., 0.], [-tilt, 0., 1.]],
            translation_mm: [42., 0., 0.],
        })
        .unwrap();
    let short = Solid::rust_spur_gear_named(spur(12, 0., None), EntityId::new(), policy).unwrap();
    assert_eq!(
        short
            .certify_rust_spur_pair(&leaning, policy, None)
            .unwrap(),
        ResultKind::Clear
    );
    let refused = leaning
        .certify_rust_spur_pair(&a, policy, None)
        .unwrap_err();
    assert_eq!(refused.kind, GeometryErrorKind::Unsupported);
    assert!(refused.message.contains("axial placement"), "{refused:?}");

    let placed = rotated(&b, PI + PI / 30., [42., 0., 0.]);
    for reach in [0., -1., f64::INFINITY, f64::NAN] {
        assert_eq!(
            short
                .rust_gear_flank_witnesses(&placed, policy, reach, None)
                .unwrap_err()
                .kind,
            GeometryErrorKind::Input
        );
    }
    let coarse = GeometryPolicy {
        linear_tolerance_mm: 1e-6,
        ..policy
    };
    assert_eq!(
        short
            .rust_gear_flank_witnesses(&placed, coarse, 0.1, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        short
            .rust_gear_flank_witnesses(&placed, policy, 0.1, Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let block = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [1, 1, 1].map(|n| multicad_kernel::rust_kernel::Rational::new(n, 1).unwrap()),
        Default::default(),
        policy,
    )
    .unwrap();
    assert_eq!(
        short
            .rust_gear_flank_witnesses(&block, policy, 0.1, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
}

fn rack(cells: u32) -> StraightRack {
    let module = 2.;
    StraightRack {
        teeth: cells,
        pitch_mm: PI * module,
        pitch_thickness_mm: PI * module / 2. - THINNING,
        pressure_angle_rad: 20_f64.to_radians(),
        addendum_mm: module,
        dedendum_mm: 1.25 * module,
        root_fillet_radius_mm: 0.4,
        backing_thickness_mm: 3.,
        face_width_mm: 4.,
    }
}

#[test]
fn finite_rack_and_spur_pinion_clear_in_either_order_with_actual_flank_points() {
    let policy = GeometryPolicy::default();
    let bar = Solid::rust_straight_rack_named(rack(20), EntityId::new(), policy).unwrap();
    let hub = Some(GearHub {
        radius_mm: 6.,
        axial_start_mm: 0.,
        axial_end_mm: 10.,
    });
    let pinion = Solid::rust_spur_gear_named(spur(24, 0., hub), EntityId::new(), policy).unwrap();
    let pitch = 24.;
    let reach = 2. * THINNING + 10. * policy.linear_tolerance_mm;
    let half_period = PI / 24.;
    for i in 0..9 {
        let drive = half_period * (f64::from(i) / 4. - 1.);
        // A pinion tooth points at the rack's central tooth space.
        let gear = rotated(&pinion, -FRAC_PI_2 + drive, [0., pitch, 0.]);
        let moved = rotated(&bar, 0., [pitch * drive, 0., 0.]);
        for (first, second, pinion_first) in [(&gear, &moved, true), (&moved, &gear, false)] {
            assert_eq!(
                first.certify_rust_rack_pair(second, policy, None).unwrap(),
                ResultKind::Clear,
                "{i}"
            );
            let bands = first
                .rust_gear_flank_witnesses(second, policy, reach, None)
                .unwrap();
            let witness = bands[0].expect("rack flank proximity");
            assert!(witness.distance_mm > 0. && witness.distance_mm <= reach);
            assert_eq!(witness.faces, None);
            // The rack point lies on its tooth flanks, below the pitch
            // circle's top, and the pinion point is farther from the rack.
            let [p, q] = witness.points_mm;
            let (pinion_point, rack_point) = if pinion_first { (p, q) } else { (q, p) };
            assert!(rack_point[1] > -2.5 && rack_point[1] < 2.);
            assert!(pinion_point[1] > rack_point[1] - reach);
        }
    }
    let wrong = rotated(&bar, 0., [PI, 0., 0.]);
    let gear = rotated(&pinion, -FRAC_PI_2, [0., pitch, 0.]);
    assert_eq!(
        gear.certify_rust_rack_pair(&wrong, policy, None).unwrap(),
        ResultKind::CommonInterior
    );
    let lifted = rotated(&pinion, -FRAC_PI_2, [0., pitch + 1., 0.]);
    assert_eq!(
        lifted.certify_rust_rack_pair(&bar, policy, None).unwrap(),
        ResultKind::Clear
    );
    assert_eq!(
        lifted
            .rust_gear_flank_witnesses(&bar, policy, reach, None)
            .unwrap(),
        vec![None]
    );
    // Two racks, a rack against nothing like a pinion and a mixed policy
    // are explicit refusals rather than results.
    assert_eq!(
        bar.certify_rust_rack_pair(&bar, policy, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let fuzzy = GeometryPolicy {
        boolean_fuzzy_mm: 1e-7,
        ..policy
    };
    assert_eq!(
        gear.certify_rust_rack_pair(&bar, fuzzy, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        gear.certify_rust_rack_pair(&bar, policy, Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let crossed = bar
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]],
            translation_mm: [0., 0., 0.],
        })
        .unwrap();
    assert_eq!(
        gear.certify_rust_rack_pair(&crossed, policy, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
}

fn twisted(teeth: u32, sweep: GearSweep) -> InvoluteGear {
    let beta = 15_f64.to_radians();
    let at = (20_f64.to_radians().portable_tan() / beta.portable_cos()).portable_atan();
    let transverse_module = 2. / beta.portable_cos();
    let pitch = transverse_module * f64::from(teeth) / 2.;
    InvoluteGear {
        teeth,
        pitch_radius_mm: pitch,
        base_radius_mm: pitch * at.portable_cos(),
        tip_radius_mm: pitch + 2.,
        root_radius_mm: pitch - 2.5,
        root_fillet_radius_mm: 0.4,
        half_base_angle_rad: (PI * transverse_module / 2. - THINNING) / (2. * pitch)
            + at.portable_tan()
            - at,
        face_width_mm: 4.,
        bore_radius_mm: 4.,
        hub: None,
        twist_rad_per_mm: if teeth == 24 { 1. } else { -1. } * beta.portable_tan() / pitch,
        sweep,
    }
}

#[test]
fn twisted_pairs_report_one_band_or_both_chevron_bands() {
    let policy = GeometryPolicy::default();
    let reach = 2. * THINNING + 10. * policy.linear_tolerance_mm;
    for (sweep, bands) in [(GearSweep::Helical, 1), (GearSweep::Chevron, 2)] {
        let build = |teeth| {
            let spec = twisted(teeth, sweep);
            if sweep == GearSweep::Helical {
                Solid::rust_helical_gear_named(spec, EntityId::new(), policy)
            } else {
                Solid::rust_herringbone_named(spec, EntityId::new(), policy)
            }
            .unwrap()
        };
        let a = build(24);
        let b = build(36);
        let center = twisted(24, sweep).pitch_radius_mm + twisted(36, sweep).pitch_radius_mm;
        let meshed = rotated(&b, PI + PI / 36., [center, 0., 0.]);
        assert_eq!(
            a.certify_rust_herringbone_pair(&meshed, policy, None)
                .unwrap(),
            ResultKind::Clear
        );
        let found = a
            .rust_gear_flank_witnesses(&meshed, policy, reach, None)
            .unwrap();
        assert_eq!(found.len(), bands, "{sweep:?}");
        for (band, witness) in found.into_iter().enumerate() {
            let witness = witness.expect("flank proximity in each band");
            assert!(witness.distance_mm > 0. && witness.distance_mm <= reach);
            let z = witness.points_mm[0][2];
            assert!(
                if band == 0 { z < 2. } else { z > 2. },
                "{sweep:?} band {band} at {z}"
            );
            assert!(witness.faces.unwrap().iter().all(|&face| face >= 2));
        }
        let apart = rotated(&b, PI + PI / 36., [center + 1., 0., 0.]);
        assert!(
            a.rust_gear_flank_witnesses(&apart, policy, reach, None)
                .unwrap()
                .iter()
                .all(Option::is_none)
        );
    }
}
