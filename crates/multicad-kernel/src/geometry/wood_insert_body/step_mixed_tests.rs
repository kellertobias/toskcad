//! Real mixed-source export regressions; every leaf retains its own proof authority.
use super::*;
use crate::{
    CurvedEdge2d, GearSweep, HelicalThread, InvoluteGear, NormalReliefMotif, RigidTransform,
};
use multicad_rust_kernel::{Rational, normal_relief::Kind};

fn pose() -> PlaneFrame {
    PlaneFrame {
        origin: [31., -17., 9.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    }
}
fn insert() -> Solid {
    Solid::rust_wood_insert(
        crate::geometry::wood_insert_world_tests::source(true),
        EntityId::new(),
        pose(),
        Default::default(),
        None,
    )
    .unwrap()
}
fn text(solid: &Solid) -> String {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("actual-mixed.step");
    solid.write_step(&path).unwrap();
    std::fs::read_to_string(path).unwrap()
}
fn compound(a: &Solid, b: &Solid) -> Solid {
    Solid::compound_named(&[a, b], EntityId::new(), Default::default()).unwrap()
}
fn gear(sweep: GearSweep, teeth: u32) -> Solid {
    let beta = 15_f64.to_radians();
    let pressure = (20_f64.to_radians().tan() / beta.cos()).atan();
    let module = 2. / beta.cos();
    let pitch = module * f64::from(teeth) / 2.;
    let spec = InvoluteGear {
        teeth,
        pitch_radius_mm: pitch,
        base_radius_mm: pitch * pressure.cos(),
        tip_radius_mm: pitch + 2.,
        root_radius_mm: pitch - 2.5,
        root_fillet_radius_mm: 0.4,
        half_base_angle_rad: (std::f64::consts::PI * module / 2. - 0.05) / (2. * pitch)
            + pressure.tan()
            - pressure,
        face_width_mm: 4.,
        bore_radius_mm: 4.,
        hub: None,
        twist_rad_per_mm: beta.tan() / pitch,
        sweep,
    };
    let build = if sweep == GearSweep::Chevron {
        Solid::rust_herringbone_named
    } else {
        Solid::rust_helical_gear_named
    };
    build(spec, EntityId::new(), Default::default())
        .unwrap()
        .transformed(RigidTransform {
            rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
            translation_mm: [-40., 11., 7.],
        })
        .unwrap()
}
#[test]
fn actual_gear_and_placed_insert_keep_independent_member_proofs() {
    for sweep in [GearSweep::Helical, GearSweep::Chevron] {
        let insert = insert();
        let ordinary = gear(sweep, 24);
        let solid = compound(&ordinary, &insert);
        let encoded = text(&solid);
        let proof = solid
            .verify_wood_insert_step_compound(&encoded, 1e-7, None)
            .unwrap()
            .unwrap();
        assert_eq!(proof.members.solid_count, 2);
        assert_eq!(proof.members.retained_members.len(), 1);
        assert_eq!(proof.members.retained_members[0].member_index, 1);
        assert_eq!(proof.members.decoded_scalar_members.len(), 1);
        let changed = compound(&gear(sweep, 36), &insert);
        assert!(
            changed
                .verify_wood_insert_step_compound(&encoded, 1e-7, None)
                .unwrap()
                .is_err()
        );
    }
}
fn relief(kind: Kind, height: f64) -> Solid {
    let policy = GeometryPolicy::default();
    let frame = pose();
    let blank = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [12, 10, 3].map(|n| Rational::new(n, 1).unwrap()),
        frame,
        policy,
    )
    .unwrap();
    let cap = PlaneFrame {
        origin: std::array::from_fn(|i| frame.origin[i] + 3. * frame.normal()[i]),
        ..frame
    };
    let face = blank
        .faces()
        .unwrap()
        .into_iter()
        .find(|f| {
            f.plane.is_some_and(|(_, n)| {
                n.into_iter()
                    .zip(frame.normal())
                    .map(|(a, b)| a * b)
                    .sum::<f64>()
                    > 0.999
            })
        })
        .unwrap();
    let feature = EntityId::new();
    blank
        .rust_normal_relief_on_cap(
            face.index,
            cap,
            &[NormalReliefMotif {
                boundary: vec![CurvedEdge2d::Circle {
                    center: [3., 3.],
                    radius_mm: 1.,
                }],
                height_mm: height,
                origin: Origin {
                    source: feature,
                    role: "mixed/real-relief".into(),
                },
            }],
            kind,
            feature,
            policy,
            None,
        )
        .unwrap()
}
#[test]
fn actual_relief_and_placed_insert_reject_changed_ordinary_source() {
    for kind in [Kind::Raised, Kind::Recessed] {
        let insert = insert();
        let ordinary = relief(kind, 0.5);
        let solid = compound(&insert, &ordinary);
        let encoded = text(&solid);
        let proof = solid
            .verify_wood_insert_step_compound(&encoded, 1e-7, None)
            .unwrap()
            .unwrap();
        assert_eq!(proof.members.solid_count, 2);
        assert_eq!(proof.members.retained_members[0].member_index, 0);
        assert_eq!(proof.members.decoded_scalar_members.len(), 1);
        assert_eq!(proof.members.decoded_scalar_members[0].0, 1);
        let changed = compound(&insert, &relief(kind, 0.75));
        assert!(
            changed
                .verify_wood_insert_step_compound(&encoded, 1e-7, None)
                .unwrap()
                .is_err()
        );
    }
}
fn entry(depth: f64) -> Solid {
    let root = 1. / (4. * 3_f64.sqrt());
    let thread = HelicalThread {
        major_radius_mm: 3.,
        pitch_mm: 1.,
        length_mm: 6.6,
        depth_mm: 0.875 * 3_f64.sqrt() / 2. - root,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand: true,
    };
    Solid::rust_helical_thread_entry_chamfer(
        thread,
        Some(root),
        depth,
        None,
        EntityId::new(),
        pose(),
        Default::default(),
    )
    .unwrap()
}
#[test]
fn actual_entry_chamfer_with_insert_preserves_source_known_mass_authority() {
    let insert = insert();
    let chamfer = entry(0.3);
    let solid = compound(&chamfer, &insert);
    let encoded = text(&solid);
    let proof = solid
        .verify_wood_insert_step_compound(&encoded, 1e-7, None)
        .unwrap()
        .unwrap();
    assert_eq!(proof.members.solid_count, 2);
    assert_eq!(proof.members.retained_members.len(), 1);
    assert_eq!(proof.members.retained_members[0].member_index, 1);
    assert!(proof.members.decoded_scalar_members.is_empty());
    assert_eq!(proof.source_known_members.len(), 1);
    assert_eq!(proof.source_known_members[0].0, 0);
    assert!(proof.source_known_members[0].1.source_volume_bounds_mm3[0] > 0.);
    let changed = compound(&entry(0.4), &insert);
    assert!(
        changed
            .verify_wood_insert_step_compound(&encoded, 1e-7, None)
            .unwrap()
            .is_err()
    );
}
