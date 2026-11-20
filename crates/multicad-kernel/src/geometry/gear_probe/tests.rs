use super::*;
use crate::{GearHub, GearSweep, InvoluteGear};
use std::f64::consts::PI;

fn p() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn spec(beta: f64, sweep: GearSweep) -> InvoluteGear {
    let module = 2.;
    let teeth = 20;
    let pressure = (20_f64.to_radians().tan() / beta.cos()).atan();
    let pitch = module * f64::from(teeth) / (2. * beta.cos());
    InvoluteGear {
        teeth,
        pitch_radius_mm: pitch,
        base_radius_mm: pitch * pressure.cos(),
        tip_radius_mm: pitch + module,
        root_radius_mm: pitch - 1.25 * module,
        root_fillet_radius_mm: 0.2 * module,
        half_base_angle_rad: PI / (2. * f64::from(teeth)) + pressure.tan() - pressure,
        face_width_mm: 8.,
        bore_radius_mm: 2.,
        hub: None,
        twist_rad_per_mm: beta.tan() / pitch,
        sweep,
    }
}

fn cube(center: [f64; 3], side: f64) -> Solid {
    Solid::box_solid([side; 3], p())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: center.map(|v| v - side / 2.),
            ..Default::default()
        })
        .unwrap()
}

fn polar(r: f64, a: f64, z: f64) -> [f64; 3] {
    [r * a.cos(), r * a.sin(), z]
}

fn intersect(a: &Solid, b: &Solid) -> Option<Result<Option<Solid>, GeometryError>> {
    a.contained_gear_boolean(b, BooleanOperation::Intersect, p())
}

#[test]
fn spur_probes_reuse_an_operand_only_when_wholly_inside_or_outside() {
    let gear =
        Solid::rust_spur_gear_named(spec(0., GearSweep::Spur), EntityId::new(), p()).unwrap();
    let inside = cube(polar(20., 0., 4.), 0.01);
    let kept = intersect(&gear, &inside).unwrap().unwrap().unwrap();
    assert!((kept.volume_mm3().unwrap() - 1e-6).abs() < 1e-15);
    // Either operand order, and a difference keeps the unchanged gear.
    assert!(intersect(&inside, &gear).unwrap().unwrap().is_some());
    let space = cube(polar(20., PI / 20., 4.), 0.01);
    assert!(intersect(&gear, &space).unwrap().unwrap().is_none());
    let difference = gear
        .contained_gear_boolean(&space, BooleanOperation::Subtract, p())
        .unwrap()
        .unwrap()
        .unwrap();
    assert_eq!(difference.volume_mm3().unwrap(), gear.volume_mm3().unwrap());
    let bore = cube([0., 0., 4.], 0.5);
    assert!(intersect(&gear, &bore).unwrap().unwrap().is_none());
    // A probe across a flank, a cap or the bore wall is never settled.
    for straddling in [
        cube(polar(20., PI / 40., 4.), 0.5),
        cube(polar(20., 0., 8.), 0.01),
        cube([2., 0., 4.], 0.01),
    ] {
        assert!(intersect(&gear, &straddling).is_none());
    }
    // Wholly above the cap is outside even over the tooth.
    assert!(
        intersect(&gear, &cube(polar(20., 0., 9.), 0.01))
            .unwrap()
            .unwrap()
            .is_none()
    );
    // Fuzzy policies are left to other providers.
    let fuzzy = GeometryPolicy {
        boolean_fuzzy_mm: 1e-6,
        ..p()
    };
    assert!(
        gear.contained_gear_boolean(&inside, BooleanOperation::Intersect, fuzzy)
            .is_none()
    );
}

#[test]
fn hubs_twisted_bands_and_gap_cores_use_their_own_sections() {
    let hubbed = Solid::rust_spur_gear_named(
        InvoluteGear {
            hub: Some(GearHub {
                radius_mm: 10.,
                axial_start_mm: -4.,
                axial_end_mm: 8.,
            }),
            ..spec(0., GearSweep::Spur)
        },
        EntityId::new(),
        p(),
    )
    .unwrap();
    assert!(
        intersect(&hubbed, &cube(polar(9., 1., -2.), 0.01))
            .unwrap()
            .unwrap()
            .is_some()
    );
    assert!(
        intersect(&hubbed, &cube(polar(11., 1., -2.), 0.01))
            .unwrap()
            .unwrap()
            .is_none()
    );
    // The tooth band and hub meet at z = 0 inside the root: still unsettled.
    assert!(intersect(&hubbed, &cube(polar(9., 1., 0.), 0.01)).is_none());

    let beta = 20_f64.to_radians();
    let s = spec(beta, GearSweep::Helical);
    let helical = Solid::rust_helical_gear_named(s, EntityId::new(), p()).unwrap();
    let z = 6.;
    let phase = s.twist_rad_per_mm * z;
    let tooth = cube(polar(s.pitch_radius_mm, phase, z), 0.01);
    assert!(intersect(&helical, &tooth).unwrap().unwrap().is_some());
    // The unrotated tooth position is a space at this height.
    let space = cube(polar(s.pitch_radius_mm, phase + PI / 20., z), 0.01);
    assert!(intersect(&helical, &space).unwrap().unwrap().is_none());
    // A tall probe sweeps through the twisted flanks.
    assert!(intersect(&helical, &cube(polar(s.pitch_radius_mm, phase, z), 3.)).is_none());

    let s = spec(beta, GearSweep::GappedChevron { gap_mm: 2. });
    let gapped = Solid::rust_double_helical_named(s, EntityId::new(), p()).unwrap();
    assert!(
        intersect(&gapped, &cube(polar(s.pitch_radius_mm, 0., 4.), 0.01))
            .unwrap()
            .unwrap()
            .is_none()
    );
    assert!(
        intersect(&gapped, &cube(polar(s.root_radius_mm - 0.5, 0.3, 4.), 0.01))
            .unwrap()
            .unwrap()
            .is_some()
    );
}

#[test]
fn internal_ring_material_is_between_cavity_and_rim() {
    let a = 20_f64.to_radians();
    let pitch = 24.;
    let ring = Solid::rust_internal_ring_named(
        InvoluteGear {
            teeth: 48,
            pitch_radius_mm: pitch,
            base_radius_mm: pitch * a.cos(),
            tip_radius_mm: pitch - 1.,
            root_radius_mm: pitch + 1.25,
            root_fillet_radius_mm: 0.2,
            half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) - a.tan() + a,
            face_width_mm: 2.,
            bore_radius_mm: 0.,
            hub: None,
            twist_rad_per_mm: 0.,
            sweep: GearSweep::InternalSpur {
                outer_radius_mm: pitch + 4.,
            },
        },
        EntityId::new(),
        p(),
    )
    .unwrap();
    let probe = |r, a| intersect(&ring, &cube(polar(r, a, 1.), 0.01));
    assert!(probe(pitch, 0.).unwrap().unwrap().is_some());
    assert!(probe(pitch, PI / 48.).unwrap().unwrap().is_none());
    assert!(probe(pitch + 3., 0.5).unwrap().unwrap().is_some());
    assert!(probe(5., 0.).unwrap().unwrap().is_none());
    assert!(probe(pitch + 4., 0.5).is_none());
}

#[test]
fn certified_gear_pairs_settle_disjoint_poses_and_refuse_interference() {
    let a = spec(0., GearSweep::Spur);
    let b = InvoluteGear {
        teeth: 30,
        pitch_radius_mm: 30.,
        base_radius_mm: 30. * 20_f64.to_radians().cos(),
        tip_radius_mm: 32.,
        root_radius_mm: 27.5,
        half_base_angle_rad: PI / 60. + 20_f64.to_radians().tan()
            - 20_f64.to_radians()
            - 0.05 / 60.,
        ..a
    };
    let a = InvoluteGear {
        half_base_angle_rad: a.half_base_angle_rad - 0.05 / 40.,
        ..a
    };
    let one = Solid::rust_spur_gear_named(a, EntityId::new(), p()).unwrap();
    let two = Solid::rust_spur_gear_named(b, EntityId::new(), p()).unwrap();
    let place = |angle: f64| {
        two.transformed(RigidTransform {
            rotation: [
                [angle.cos(), -angle.sin(), 0.],
                [angle.sin(), angle.cos(), 0.],
                [0., 0., 1.],
            ],
            translation_mm: [50., 0., 0.],
        })
        .unwrap()
    };
    let meshing = place(PI + PI / 30.);
    assert!(
        one.boolean(&meshing, BooleanOperation::Intersect, p())
            .unwrap()
            .is_none()
    );
    let kept = one
        .boolean(&meshing, BooleanOperation::Subtract, p())
        .unwrap()
        .unwrap();
    assert_eq!(kept.volume_mm3().unwrap(), one.volume_mm3().unwrap());
    let error = one
        .boolean(&place(PI), BooleanOperation::Intersect, p())
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert!(error.message.contains("interfere"), "{}", error.message);
    // The serial and native-parallel entry points agree for Rust solids.
    assert!(one.intersect_parallel(&meshing, p()).unwrap().is_none());
    // The clear pose has a certified gap below the total thinning; its
    // witness points are actual material points.
    let gap = one
        .minimum_distance(
            crate::BRepSelection::Whole,
            &meshing,
            crate::BRepSelection::Whole,
        )
        .unwrap();
    assert!(
        gap.distance_mm > 0. && gap.distance_mm < 0.2,
        "{}",
        gap.distance_mm
    );
    let w = gap.witnesses[0];
    let span = (0..3)
        .map(|i| (w.a_mm[i] - w.b_mm[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!((span - gap.distance_mm).abs() < 1e-12);
    // An interfering pose has no distance certificate.
    assert!(
        one.minimum_distance(
            crate::BRepSelection::Whole,
            &place(PI),
            crate::BRepSelection::Whole
        )
        .is_err()
    );
}

#[test]
fn worm_probes_follow_the_placed_thread() {
    let worm = Solid::rust_za_worm_named(
        crate::ZaWorm {
            axial_module_mm: 1.,
            threaded_length_mm: 16.,
            bore_radius_mm: 1.,
            axial_thinning_mm: 0.05,
        },
        EntityId::new(),
        p(),
    )
    .unwrap()
    .transformed(RigidTransform {
        translation_mm: [10., 0., 0.],
        ..Default::default()
    })
    .unwrap();
    let at = |x, y, z| intersect(&worm, &cube([10. + x, y, z], 0.01));
    assert!(
        at(2., 1., 3.)
            .unwrap()
            .unwrap()
            .unwrap()
            .volume_mm3()
            .unwrap()
            > 0.
    );
    assert!(at(0., 0., 0.).unwrap().unwrap().is_none());
    assert!(at(7., 0., 0.).unwrap().unwrap().is_none());
    assert!(at(1., 0., 0.).is_none());
}

#[test]
fn generated_wheel_probes_use_the_certified_distance_field() {
    let wheel = Solid::rust_generated_za_worm_wheel_named(
        crate::GeneratedZaWormWheel {
            teeth: 32,
            axial_module_mm: 1.,
            face_width_mm: 2.,
            bore_radius_mm: 2.,
            pitch_thinning_mm: 0.05,
        },
        EntityId::new(),
        p(),
    )
    .unwrap();
    let tooth = PI / 32.;
    assert!(
        intersect(&wheel, &cube(polar(16., tooth, 0.), 0.005))
            .unwrap()
            .unwrap()
            .is_some()
    );
    assert!(
        intersect(&wheel, &cube(polar(16., 0., 0.), 0.005))
            .unwrap()
            .unwrap()
            .is_none()
    );
    assert!(intersect(&wheel, &cube(polar(16., tooth / 2., 0.), 0.5)).is_none());
}
