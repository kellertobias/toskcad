//! Gapped double-helical gears: two opposite-hand tooth bands joined by
//! the root core through the authored central gap.
use multicad_core::{
    document::EntityId,
    references::{AnalyticGeometry, TopologyKind},
};
use multicad_kernel::{
    DrawingProjection, DrawingSolid, GearSweep, GeometryErrorKind, GeometryPolicy, InvoluteGear,
    MeshOptions, PlaneFrame, RigidTransform, Solid, project_drawing,
    rust_kernel::spur_gear::clearance::ResultKind,
};
use multicad_precision::PortableMath as _;
use std::f64::consts::PI;

const GAP: f64 = 3.;
const WIDTH: f64 = 10.;

fn specification(teeth: u32) -> InvoluteGear {
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
        half_base_angle_rad: (PI * transverse_module / 2. - 0.05) / (2. * pitch)
            + at.portable_tan()
            - at,
        face_width_mm: WIDTH,
        bore_radius_mm: 4.,
        hub: None,
        twist_rad_per_mm: if teeth == 24 { 1. } else { -1. } * beta.portable_tan() / pitch,
        sweep: GearSweep::GappedChevron { gap_mm: GAP },
    }
}

fn transform() -> RigidTransform {
    RigidTransform {
        rotation: [[0.8, 0., 0.6], [0., 1., 0.], [-0.6, 0., 0.8]],
        translation_mm: [10., -20., 30.],
    }
}

#[test]
fn named_world_topology_metrics_mesh_and_verified_step() {
    let id = EntityId::new();
    let spec = specification(24);
    let solid = Solid::rust_double_helical_named(spec, id, Default::default()).unwrap();
    let world = solid.transformed(transform()).unwrap();
    let snapshot = solid.topology_snapshot(id, 3).unwrap();
    snapshot.validate().unwrap();
    let display = world.topology_with_display(id, 3, Some(0.005)).unwrap();
    display.topology.validate().unwrap();
    let roles = |role: &str| {
        snapshot
            .elements
            .iter()
            .filter(|element| {
                element.key.kind == TopologyKind::Face
                    && element.origins.iter().any(|origin| origin.role == role)
            })
            .count()
    };
    assert_eq!(roles("gear/root.core"), 1);
    assert_eq!(roles("gear/gap.lower"), 24);
    assert_eq!(roles("gear/gap.upper"), 24);
    assert_eq!(roles("gear/bore"), 1);
    assert_eq!(roles("gear/lower/tooth/0/involute.left"), 1);
    assert_eq!(roles("gear/upper/tooth/23/root"), 1);
    // Every root arc and the core seam carry exact analytic geometry.
    let circular = snapshot
        .elements
        .iter()
        .filter(|element| {
            element.key.kind == TopologyKind::Edge
                && matches!(element.geometry, Some(AnalyticGeometry::CircularArc { .. }))
        })
        .count();
    assert!(circular >= 48, "{circular}");

    let band = (WIDTH - GAP) / 2.;
    let root = spec.root_radius_mm;
    let bore = spec.bore_radius_mm;
    let faces = world.faces().unwrap();
    let local = solid.faces().unwrap();
    let core = local
        .iter()
        .find(|face| (face.area_mm2 - 2. * PI * root * GAP).abs() < 1e-9)
        .expect("root core area");
    assert!((core.center_mm[2] - WIDTH / 2.).abs() < 1e-12);
    // Both shoulder rows cover the profile outside the root circle.
    let shoulders = local
        .iter()
        .filter(|face| {
            face.plane
                .is_some_and(|(origin, _)| (origin[2] - band).abs() < 1e-12)
                && face.center_mm[2] == band
        })
        .collect::<Vec<_>>();
    assert_eq!(shoulders.len(), 24);
    let section = shoulders.iter().map(|face| face.area_mm2).sum::<f64>();
    for face in &shoulders {
        let radial = face.center_mm[0].portable_hypot(face.center_mm[1]);
        assert!(radial > root && radial < spec.tip_radius_mm);
    }
    let volume = solid.volume_mm3().unwrap();
    let profile_area = section + PI * (root * root - bore * bore);
    assert!(
        (volume - (2. * band * profile_area + GAP * PI * (root * root - bore * bore))).abs() < 1e-6,
        "{volume}"
    );
    assert_eq!(world.volume_mm3().unwrap(), volume);
    assert_eq!(faces.len(), local.len());

    let result = world
        .manufacturing_mesh_with_receipt(MeshOptions {
            linear_deflection_mm: 0.005,
            ..Default::default()
        })
        .unwrap();
    assert!(
        result
            .mesh
            .face_indices
            .iter()
            .all(|&face| (face as usize) < faces.len())
    );

    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("double-helical.step");
    assert!(world.supports_step_export());
    world.write_step(&path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let verified = world.verify_herringbone_step(&text).unwrap().unwrap();
    assert!((verified - volume).abs() < 1e-5);
    assert!(solid.verify_herringbone_step(&text).unwrap().is_err());
    let compound =
        Solid::compound_named(&[&solid, &world], EntityId::new(), Default::default()).unwrap();
    let compound_path = directory.path().join("double-helical-compound.step");
    compound.write_step(&compound_path).unwrap();
    let text = std::fs::read_to_string(&compound_path).unwrap();
    let (count, total, _) = compound
        .verify_twisted_step_compound(&text, 1e-7)
        .unwrap()
        .unwrap();
    assert_eq!(count, 2);
    assert!((total - 2. * volume).abs() < 1e-5);
    if multicad_kernel::OCCT_LINKED {
        let read = multicad_kernel::Shape::read_step(&path).unwrap();
        assert!(read.valid().unwrap());
        assert_eq!(Solid::step_solids(&path).unwrap().len(), 1);
        assert!((read.volume().unwrap() - verified).abs() < 2e-4);
    }
}

#[test]
fn native_double_helical_matches_volume_area_and_bounds() {
    if !multicad_kernel::OCCT_LINKED {
        return;
    }
    for (teeth, hub_span, bore) in [
        (24, None, 4.),
        (36, None, 0.),
        (24, Some([-2., 2.]), 4.),
        (24, Some([8., 12.]), 0.),
        (24, Some([-2., 12.]), 4.),
        (24, Some([4., 6.]), 4.),
    ] {
        let mut spec = specification(teeth);
        spec.bore_radius_mm = bore;
        spec.hub = hub_span.map(|[start, end]| multicad_kernel::GearHub {
            radius_mm: 8.,
            axial_start_mm: start,
            axial_end_mm: end,
        });
        let feature = EntityId::new();
        let rust = Solid::rust_double_helical_named(spec, feature, Default::default())
            .unwrap()
            .transformed(transform())
            .unwrap();
        let native = Solid::involute_gear_named(spec, feature, Default::default())
            .unwrap()
            .transformed(transform())
            .unwrap();
        let label = format!("{teeth} {hub_span:?} {bore}");
        let volume = rust.volume_mm3().unwrap();
        assert!(
            (native.volume_mm3().unwrap() - volume).abs() < 2e-4,
            "{label}: {} vs {volume}",
            native.volume_mm3().unwrap()
        );
        let area = |solid: &Solid| {
            solid
                .faces()
                .unwrap()
                .iter()
                .map(|face| face.area_mm2)
                .sum::<f64>()
        };
        assert!(
            (area(&rust) - area(&native)).abs() < 2e-4,
            "{label}: area {} vs {}",
            area(&rust),
            area(&native)
        );
        for (a, b) in [rust.bounds().unwrap().min, rust.bounds().unwrap().max]
            .into_iter()
            .zip([native.bounds().unwrap().min, native.bounds().unwrap().max])
        {
            assert!(
                a.into_iter().zip(b).all(|(a, b)| (a - b).abs() < 1e-5),
                "{label}: {a:?} {b:?}"
            );
        }
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("rust.step");
        rust.write_step(&path).unwrap();
        let read = multicad_kernel::Shape::read_step(&path).unwrap();
        assert!(read.valid().unwrap(), "{label}");
        assert!((read.volume().unwrap() - volume).abs() < 2e-4, "{label}");
    }
}

#[test]
fn invalid_gaps_and_mismatched_sweeps_are_refused() {
    let policy = GeometryPolicy::default();
    let mut chevron = specification(24);
    chevron.sweep = GearSweep::Chevron;
    assert_eq!(
        Solid::rust_double_helical_named(chevron, EntityId::new(), policy)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    for gap in [0., -1., WIDTH, WIDTH + 1., f64::NAN] {
        let mut spec = specification(24);
        spec.sweep = GearSweep::GappedChevron { gap_mm: gap };
        assert!(
            Solid::rust_double_helical_named(spec, EntityId::new(), policy).is_err(),
            "{gap}"
        );
    }
    let mut gapped = specification(24);
    gapped.sweep = GearSweep::GappedChevron { gap_mm: GAP };
    assert!(Solid::rust_herringbone_named(gapped, EntityId::new(), policy).is_err());
}

#[test]
fn pair_certificate_clears_the_meshing_pose_and_rejects_wrong_ones() {
    let policy = GeometryPolicy::default();
    let a = Solid::rust_double_helical_named(specification(24), EntityId::new(), policy).unwrap();
    let b = Solid::rust_double_helical_named(specification(36), EntityId::new(), policy).unwrap();
    let center = specification(24).pitch_radius_mm + specification(36).pitch_radius_mm;
    let placed = |solid: &Solid, phase: f64| {
        let (s, c) = phase.portable_sin_cos();
        solid
            .transformed(RigidTransform {
                rotation: [[c, -s, 0.], [s, c, 0.], [0., 0., 1.]],
                translation_mm: [center, 0., 0.],
            })
            .unwrap()
    };
    assert_eq!(
        a.certify_rust_herringbone_pair(&placed(&b, PI + PI / 36.), policy, None)
            .unwrap(),
        ResultKind::Clear
    );
    assert_eq!(
        a.certify_rust_herringbone_pair(&placed(&b, PI), policy, None)
            .unwrap(),
        ResultKind::CommonInterior
    );
    // A different gap needs the core-only section and is not certified.
    let mut narrower = specification(36);
    narrower.sweep = GearSweep::GappedChevron { gap_mm: 2. };
    let narrower = Solid::rust_double_helical_named(narrower, EntityId::new(), policy).unwrap();
    assert!(
        a.certify_rust_herringbone_pair(&placed(&narrower, PI + PI / 36.), policy, None)
            .is_err()
    );
}

/// Gapped gears draw through the exact-curve scene engine. Expectations
/// come from the specification and the exact B-rep face areas: the axial
/// view spans the bore and tip circles, a cut through a tooth band encloses
/// the band's transverse profile and a cut through the gap the root core.
#[test]
fn drawings_of_gapped_gears_follow_the_bands_and_the_core() {
    use multicad_core::{
        drawings::{VectorCurve, VectorRole},
        materials::HatchPattern,
    };
    use multicad_kernel::DrawingHatch;
    let spec = specification(24);
    let solid =
        Solid::rust_double_helical_named(spec, EntityId::new(), Default::default()).unwrap();
    let bounds = solid.bounds().unwrap();
    assert!(bounds.min[2].abs() < 1e-9 && (bounds.max[2] - WIDTH).abs() < 1e-9);
    let (root, bore, tip) = (spec.root_radius_mm, spec.bore_radius_mm, spec.tip_radius_mm);
    let band = (WIDTH - GAP) / 2.;
    let shoulders: f64 = solid
        .faces()
        .unwrap()
        .iter()
        .filter(|face| {
            face.plane
                .is_some_and(|(origin, _)| (origin[2] - band).abs() < 1e-12)
                && face.center_mm[2] == band
        })
        .map(|face| face.area_mm2)
        .sum();
    let core = PI * (root * root - bore * bore);
    let radius = |p: &[f64; 2]| p[0].portable_hypot(p[1]);

    let view = project_drawing(
        &[DrawingSolid {
            solid: &solid,
            hatch: None,
        }],
        DrawingProjection::default(),
    )
    .unwrap();
    assert!(view.iter().all(|v| v.role == VectorRole::Visible));
    let outer = view
        .iter()
        .flat_map(|v| v.preview_mm.iter().map(radius))
        .fold(0., f64::max);
    assert!((outer - tip).abs() < 1e-6, "{outer}");
    assert!(
        view.iter()
            .all(|v| v.preview_mm.iter().all(|p| radius(p) >= bore - 1e-9))
    );
    assert!(view.iter().any(|v| matches!(
        v.curve,
        VectorCurve::Ellipse { major_mm, minor_mm, .. }
            if (major_mm - bore).abs() < 1e-9 && (minor_mm - bore).abs() < 1e-9
    )));

    let cut = |z: f64, section_only: bool| {
        project_drawing(
            &[DrawingSolid {
                solid: &solid,
                hatch: Some(DrawingHatch {
                    angle_rad: 0.,
                    spacing_mm: 1.,
                    line_width_mm: 0.18,
                    pattern: HatchPattern::Solid,
                }),
            }],
            DrawingProjection {
                frame: PlaneFrame {
                    origin: [0., 0., z],
                    ..Default::default()
                },
                section_offset_mm: Some(0.),
                section_only,
                ..Default::default()
            },
        )
        .unwrap()
    };
    let area = |vectors: &[multicad_kernel::ProjectedVector]| -> f64 {
        vectors
            .iter()
            .filter_map(|v| match v.curve {
                VectorCurve::HatchTriangle { points: [a, b, c] } => {
                    Some(((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.)
                }
                _ => None,
            })
            .sum()
    };
    for z in [1., WIDTH - 1.] {
        let section = cut(z, true);
        let profile = shoulders + core;
        println!(
            "{{\"case\":\"gapped cut at {z}\",\"hatch_mm2\":{},\"profile_mm2\":{profile},\"vectors\":{}}}",
            area(&section),
            section.len()
        );
        assert!(
            (area(&section) - profile).abs() <= 1e-6 * profile,
            "{z}: {} against {profile}",
            area(&section)
        );
        for vector in section.iter().filter(|v| v.role == VectorRole::Section) {
            assert!(vector.preview_mm.iter().all(|p| {
                let r = radius(p);
                (r - bore).abs() < 1e-6 || (r >= root - 1e-6 && r <= tip + 1e-6)
            }));
        }
    }
    let gap = cut(WIDTH / 2., true);
    println!(
        "{{\"case\":\"gapped cut in the gap\",\"hatch_mm2\":{},\"core_mm2\":{core},\"view_vectors\":{},\"outer_mm\":{outer}}}",
        area(&gap),
        view.len()
    );
    assert!((area(&gap) - core).abs() <= 1e-6 * core, "{}", area(&gap));
    let mut circles: Vec<f64> = gap
        .iter()
        .filter(|v| v.role == VectorRole::Section)
        .map(|v| match v.curve {
            VectorCurve::Ellipse {
                major_mm, minor_mm, ..
            } if (major_mm - minor_mm).abs() < 1e-9 => major_mm,
            ref other => panic!("{other:?}"),
        })
        .collect();
    circles.sort_by(f64::total_cmp);
    circles.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    assert_eq!(circles.len(), 2, "{circles:?}");
    assert!((circles[0] - bore).abs() < 1e-9 && (circles[1] - root).abs() < 1e-9);
    // The retained view behind the cut keeps the far outline.
    let retained = cut(1., false);
    assert!(retained.iter().any(|v| v.role == VectorRole::Section));
    assert!(retained.iter().any(|v| v.role == VectorRole::Visible));
}
