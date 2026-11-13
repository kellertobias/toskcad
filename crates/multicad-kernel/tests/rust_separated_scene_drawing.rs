//! Mixed Rust scenes with single-body drawing families.
//!
//! When the members' projected bounds are separated, no line of sight meets
//! two members and each section cut is independent, so the scene drawing is
//! exactly the union of the member drawings. Overlapping scenes keep their
//! explicit refusal.
use multicad_core::document::EntityId;
use multicad_kernel::{
    DrawingProjection, DrawingSolid, GearSweep, GeometryErrorKind, GeometryPolicy, InvoluteGear,
    PlaneFrame, ProjectedVector, RigidTransform, Solid, project_drawing,
};
use std::f64::consts::PI;

#[cfg(feature = "occt")]
mod native_drawing;

fn spur() -> Solid {
    let a = 20_f64.to_radians();
    let pitch = 12.;
    Solid::rust_spur_gear_named(
        InvoluteGear {
            hub: None,
            teeth: 24,
            pitch_radius_mm: pitch,
            base_radius_mm: pitch * a.cos(),
            tip_radius_mm: pitch + 1.,
            root_radius_mm: pitch - 1.25,
            root_fillet_radius_mm: 0.2,
            half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) + a.tan() - a,
            face_width_mm: 4.,
            bore_radius_mm: 3.,
            twist_rad_per_mm: 0.,
            sweep: GearSweep::Spur,
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn stock(x: f64) -> Solid {
    use multicad_rust_kernel::Rational;
    Solid::rust_box(
        EntityId::new().uuid_bits(),
        [20, 12, 6].map(|value| Rational::new(value, 1).unwrap()),
        PlaneFrame {
            origin: [x, -6., 0.],
            ..Default::default()
        },
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn shifted(solid: &Solid, x: f64) -> Solid {
    solid
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: [x, 0., 0.],
        })
        .unwrap()
}

fn draw(
    solids: &[&Solid],
    options: DrawingProjection,
) -> Result<Vec<ProjectedVector>, multicad_kernel::GeometryError> {
    let inputs: Vec<_> = solids
        .iter()
        .map(|&solid| DrawingSolid { solid, hatch: None })
        .collect();
    project_drawing(&inputs, options)
}

fn options() -> Vec<(&'static str, DrawingProjection)> {
    let axial = PlaneFrame::default();
    let cut = PlaneFrame {
        origin: [0., 0., 2.],
        ..Default::default()
    };
    let section = DrawingProjection {
        frame: cut,
        section_offset_mm: Some(0.),
        ..Default::default()
    };
    vec![
        (
            "axial",
            DrawingProjection {
                frame: axial,
                ..Default::default()
            },
        ),
        (
            "axial hidden",
            DrawingProjection {
                frame: axial,
                hidden_lines: true,
                ..Default::default()
            },
        ),
        ("section", section),
        (
            "section only",
            DrawingProjection {
                section_only: true,
                ..section
            },
        ),
        (
            "finite depth",
            DrawingProjection {
                section_depth_mm: Some(1.),
                ..section
            },
        ),
    ]
}

#[test]
fn separated_members_draw_as_the_union_of_their_exact_drawings() {
    let gear = spur();
    let other = shifted(&gear, 40.);
    let plate = stock(-50.);
    for (label, options) in options() {
        for scene in [vec![&gear, &plate], vec![&plate, &gear, &other]] {
            let composed = draw(&scene, options).unwrap();
            let mut expected = Vec::new();
            for (source, member) in scene.iter().enumerate() {
                expected.extend(
                    draw(&[member], options)
                        .unwrap()
                        .into_iter()
                        .map(|mut vector| {
                            vector.source_index = source;
                            vector
                        }),
                );
            }
            assert_eq!(composed, expected, "{label}");
            assert!(
                composed.iter().any(|v| v.source_index == scene.len() - 1),
                "{label}"
            );
        }
    }
}

#[test]
fn overlapping_members_and_exhausted_budgets_refuse_explicitly() {
    let gear = spur();
    // The plate lies under the gear in the axial view.
    let plate = stock(-5.);
    for (label, options) in options() {
        let error = draw(&[&gear, &plate], options).unwrap_err();
        assert_eq!(
            error.kind,
            GeometryErrorKind::Unsupported,
            "{label}: {error:?}"
        );
    }
    let apart = stock(-50.);
    let error = draw(
        &[&gear, &apart],
        DrawingProjection {
            maximum_curves: 120,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert_eq!(error.kind, GeometryErrorKind::Limit);
}

#[cfg(feature = "occt")]
#[test]
fn separated_scenes_match_the_native_projection_of_both_bodies() {
    let gear = spur();
    let plate = stock(-50.);
    let axial = DrawingProjection::default();
    native_drawing::assert_matches_native(&[&gear, &plate], axial, "separated axial");
    native_drawing::assert_matches_native(
        &[&gear, &plate],
        DrawingProjection {
            hidden_lines: true,
            ..axial
        },
        "separated axial hidden",
    );
    // OpenCascade bridges gaps in its gear flank sections; its axial view of
    // the same import is the reference for the gear's section rims.
    native_drawing::assert_matches_native_with(
        &[&gear, &plate],
        DrawingProjection {
            frame: PlaneFrame {
                origin: [0., 0., 2.],
                ..Default::default()
            },
            section_offset_mm: Some(0.),
            ..Default::default()
        },
        &[axial],
        false,
        "separated section",
    );
}
