//! Sections of straight Rust spur gears and internal rings across their axis.
//!
//! Every cut through a straight tooth band repeats the exact cap profile. The
//! section rims must equal the exact axial-view rims, and the hatch must cover
//! the exact cap area of material without entering the bore or the cavity.
use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
};
use multicad_kernel::{
    DrawingHatch, DrawingProjection, DrawingSolid, GearHub, GearSweep, GeometryErrorKind,
    GeometryPolicy, InvoluteGear, PlaneFrame, ProjectedVector, RigidTransform, Solid,
    project_drawing,
};
use std::f64::consts::PI;

#[cfg(feature = "occt")]
mod native_drawing;

const SPACING: f64 = 0.25;
/// Hatch boundaries sample exact edges with chords no farther than a quarter
/// of the default minimum feature from the edge, as for herringbone sections.
const SAMPLED: f64 = 2.5e-6 + 1e-9;

fn spur(bore: f64, hub: Option<GearHub>) -> InvoluteGear {
    let a = 20_f64.to_radians();
    let pitch = 12.;
    InvoluteGear {
        hub,
        teeth: 24,
        pitch_radius_mm: pitch,
        base_radius_mm: pitch * a.cos(),
        tip_radius_mm: pitch + 1.,
        root_radius_mm: pitch - 1.25,
        root_fillet_radius_mm: 0.2,
        half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) + a.tan() - a,
        face_width_mm: 4.,
        bore_radius_mm: bore,
        twist_rad_per_mm: 0.,
        sweep: GearSweep::Spur,
    }
}

fn ring() -> InvoluteGear {
    let a = 20_f64.to_radians();
    let pitch = 24.;
    InvoluteGear {
        hub: None,
        teeth: 48,
        pitch_radius_mm: pitch,
        base_radius_mm: pitch * a.cos(),
        tip_radius_mm: pitch - 1.,
        root_radius_mm: pitch + 1.25,
        root_fillet_radius_mm: 0.2,
        half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) - a.tan() + a,
        face_width_mm: 2.,
        bore_radius_mm: 0.,
        twist_rad_per_mm: 0.,
        sweep: GearSweep::InternalSpur {
            outer_radius_mm: pitch + 4.,
        },
    }
}

fn hatch() -> Option<DrawingHatch> {
    hatch_spaced(SPACING)
}

fn hatch_spaced(spacing_mm: f64) -> Option<DrawingHatch> {
    Some(DrawingHatch {
        angle_rad: 0.37,
        spacing_mm,
        line_width_mm: 0.18,
        pattern: multicad_core::materials::HatchPattern::Lines,
    })
}

/// A tilted placement exercises proper rotations and translation.
fn placed(solid: Solid) -> (Solid, PlaneFrame) {
    let n = 14_f64.sqrt();
    let [x, y, z] = [1. / n, 2. / n, 3. / n];
    let (s, c) = 0.7_f64.sin_cos();
    let t = 1. - c;
    let rotation = [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ];
    let translation_mm = [30., -20., 10.];
    let solid = solid
        .transformed(RigidTransform {
            rotation,
            translation_mm,
        })
        .unwrap();
    let column = |j: usize| [rotation[0][j], rotation[1][j], rotation[2][j]];
    (
        solid,
        PlaneFrame {
            origin: translation_mm,
            u: column(0),
            v: column(1),
        },
    )
}

fn section_frame(placement: PlaneFrame, z: f64, toward_cap: bool) -> PlaneFrame {
    let n = placement.normal();
    PlaneFrame {
        origin: std::array::from_fn(|i| placement.origin[i] + n[i] * z),
        u: placement.u,
        v: if toward_cap {
            placement.v
        } else {
            placement.v.map(|x| -x)
        },
    }
}

fn draw(solid: &Solid, options: DrawingProjection) -> Vec<ProjectedVector> {
    project_drawing(
        &[DrawingSolid {
            solid,
            hatch: hatch(),
        }],
        options,
    )
    .unwrap()
}

fn curves(vectors: &[ProjectedVector], role: VectorRole) -> Vec<&VectorCurve> {
    vectors
        .iter()
        .filter(|v| v.role == role)
        .map(|v| &v.curve)
        .collect()
}

/// Exact curve data equal within rounding of the translated section frame.
fn same_curves(actual: &[&VectorCurve], expected: &[&VectorCurve], label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label}");
    let close = |a: &[f64], b: &[f64]| {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-9)
    };
    let flat = |curve: &VectorCurve| -> Vec<f64> {
        match curve {
            VectorCurve::Line { a, b } => [a, b].into_iter().flatten().copied().collect(),
            VectorCurve::Ellipse {
                center,
                u,
                v,
                major_mm,
                minor_mm,
                start_rad,
                end_rad,
            } => [center, u, v]
                .into_iter()
                .flatten()
                .copied()
                .chain([*major_mm, *minor_mm, *start_rad, *end_rad])
                .collect(),
            VectorCurve::BSpline {
                poles,
                weights,
                knots,
                first,
                last,
                ..
            } => poles
                .iter()
                .flatten()
                .chain(weights)
                .chain(knots)
                .copied()
                .chain([*first, *last])
                .collect(),
            VectorCurve::HatchTriangle { points } => points.iter().flatten().copied().collect(),
        }
    };
    for (a, b) in actual.iter().zip(expected) {
        assert_eq!(
            std::mem::discriminant(*a),
            std::mem::discriminant(*b),
            "{label}"
        );
        assert!(close(&flat(a), &flat(b)), "{label}: {a:?} vs {b:?}");
    }
}

fn hatch_area(vectors: &[ProjectedVector], spacing: f64) -> f64 {
    vectors
        .iter()
        .filter(|v| v.role == VectorRole::Hatch)
        .map(|v| match v.curve {
            VectorCurve::Line { a, b } => (a[0] - b[0]).hypot(a[1] - b[1]) * spacing,
            _ => panic!("hatch strokes are lines"),
        })
        .sum()
}

fn cap_area(solid: &Solid) -> f64 {
    solid
        .faces()
        .unwrap()
        .iter()
        .filter(|face| face.plane.is_some())
        .map(|face| face.area_mm2)
        .fold(0., f64::max)
}

fn hatch_radii(vectors: &[ProjectedVector]) -> [f64; 2] {
    let mut range = [f64::INFINITY, 0_f64];
    for vector in vectors.iter().filter(|v| v.role == VectorRole::Hatch) {
        let VectorCurve::Line { a, b } = vector.curve else {
            panic!("hatch strokes are lines");
        };
        for q in [0., 0.5, 1.] {
            let r = (a[0] + (b[0] - a[0]) * q).hypot(a[1] + (b[1] - a[1]) * q);
            range = [range[0].min(r), range[1].max(r)];
        }
    }
    range
}

#[test]
fn spur_sections_repeat_the_exact_cap_rims_and_hatch_only_material() {
    let policy = GeometryPolicy::default();
    let (solid, placement) =
        placed(Solid::rust_spur_gear_named(spur(3., None), EntityId::new(), policy).unwrap());
    let area = cap_area(&solid);
    let annulus = PI * (13_f64.powi(2) - 3_f64.powi(2));
    assert!(area > 0.7 * annulus && area < annulus, "{area}");
    for toward_cap in [true, false] {
        for z in [0.5, 3.7] {
            let frame = section_frame(placement, z, toward_cap);
            let axial = draw(
                &solid,
                DrawingProjection {
                    frame,
                    ..Default::default()
                },
            );
            let base = DrawingProjection {
                frame,
                section_offset_mm: Some(0.),
                ..Default::default()
            };
            for (label, options) in [
                (
                    "section only",
                    DrawingProjection {
                        section_only: true,
                        ..base
                    },
                ),
                ("through", base),
                (
                    "finite depth",
                    DrawingProjection {
                        section_depth_mm: Some(0.3),
                        ..base
                    },
                ),
            ] {
                let vectors = draw(&solid, options);
                let label = format!("{label} z={z} toward cap {toward_cap}");
                same_curves(
                    &curves(&vectors, VectorRole::Section),
                    &curves(&axial, VectorRole::Visible),
                    &label,
                );
                let behind = curves(&vectors, VectorRole::Visible);
                if label.starts_with("section only") {
                    assert!(behind.is_empty());
                } else {
                    same_curves(&behind, &curves(&axial, VectorRole::Visible), &label);
                }
                let hatched = hatch_area(&vectors, SPACING);
                assert!(
                    (hatched - area).abs() < 0.01 * area,
                    "{label} z={z}: hatch {hatched} vs exact cap {area}"
                );
                let [inner, outer] = hatch_radii(&vectors);
                assert!(
                    inner >= 3. - SAMPLED && outer <= 13. + SAMPLED,
                    "{inner} {outer}"
                );
            }
        }
    }
}

#[test]
fn ring_sections_hatch_between_the_cavity_and_the_outside_rim() {
    let policy = GeometryPolicy::default();
    let (solid, placement) =
        placed(Solid::rust_internal_ring_named(ring(), EntityId::new(), policy).unwrap());
    let area = cap_area(&solid);
    let frame = section_frame(placement, 1., true);
    let axial = draw(
        &solid,
        DrawingProjection {
            frame,
            ..Default::default()
        },
    );
    // A finer spacing over the 48-tooth cavity exhausts the hatch work
    // budget and refuses with Limit rather than publishing a partial hatch.
    let options = DrawingProjection {
        frame,
        section_offset_mm: Some(0.),
        ..Default::default()
    };
    let error = project_drawing(
        &[DrawingSolid {
            solid: &solid,
            hatch: hatch_spaced(0.1),
        }],
        options,
    )
    .unwrap_err();
    assert_eq!(error.kind, GeometryErrorKind::Limit);
    let vectors = project_drawing(
        &[DrawingSolid {
            solid: &solid,
            hatch: hatch_spaced(0.5),
        }],
        options,
    )
    .unwrap();
    same_curves(
        &curves(&vectors, VectorRole::Section),
        &curves(&axial, VectorRole::Visible),
        "ring",
    );
    let hatched = hatch_area(&vectors, 0.5);
    assert!((hatched - area).abs() < 0.02 * area, "{hatched} vs {area}");
    let [inner, outer] = hatch_radii(&vectors);
    // Material lies between the tip circle (23) and the outside rim (28).
    assert!(
        inner >= 23. - SAMPLED && outer <= 28. + SAMPLED,
        "{inner} {outer}"
    );
}

#[test]
fn planes_beside_the_band_retain_the_whole_gear_or_nothing() {
    let policy = GeometryPolicy::default();
    let (solid, placement) =
        placed(Solid::rust_spur_gear_named(spur(3., None), EntityId::new(), policy).unwrap());
    let above = section_frame(placement, 6., true);
    let axial = draw(
        &solid,
        DrawingProjection {
            frame: above,
            ..Default::default()
        },
    );
    let behind = draw(
        &solid,
        DrawingProjection {
            frame: above,
            section_offset_mm: Some(0.),
            ..Default::default()
        },
    );
    assert_eq!(behind, axial, "the whole gear lies behind the plane");
    for options in [
        DrawingProjection {
            frame: above,
            section_offset_mm: Some(0.),
            section_only: true,
            ..Default::default()
        },
        DrawingProjection {
            frame: above,
            section_offset_mm: Some(0.),
            section_depth_mm: Some(1.),
            ..Default::default()
        },
        DrawingProjection {
            frame: section_frame(placement, 6., false),
            section_offset_mm: Some(0.),
            ..Default::default()
        },
    ] {
        assert!(draw(&solid, options).is_empty(), "{options:?}");
    }
}

#[test]
fn lateral_windows_clip_the_sampled_section_inside_the_window() {
    let policy = GeometryPolicy::default();
    let (solid, placement) =
        placed(Solid::rust_spur_gear_named(spur(3., None), EntityId::new(), policy).unwrap());
    let vectors = draw(
        &solid,
        DrawingProjection {
            frame: section_frame(placement, 2., true),
            section_offset_mm: Some(0.),
            section_u_mm: Some([2., 20.]),
            section_v_mm: Some([-5., 5.]),
            ..Default::default()
        },
    );
    assert!(vectors.iter().any(|v| v.role == VectorRole::Section));
    for vector in &vectors {
        for p in &vector.preview_mm {
            assert!((2. - 1e-6..=20. + 1e-6).contains(&p[0]), "{p:?}");
            assert!((-5. - 1e-6..=5. + 1e-6).contains(&p[1]), "{p:?}");
        }
    }
}

#[test]
fn hub_sections_and_inclined_planes_remain_explicit_findings() {
    let policy = GeometryPolicy::default();
    let hubbed = Solid::rust_spur_gear_named(
        spur(
            3.,
            Some(GearHub {
                radius_mm: 6.,
                axial_start_mm: -3.,
                axial_end_mm: 4.,
            }),
        ),
        EntityId::new(),
        policy,
    )
    .unwrap();
    let plain = Solid::rust_spur_gear_named(spur(3., None), EntityId::new(), policy).unwrap();
    let s = 0.5_f64.sqrt();
    for (solid, frame) in [
        (
            &hubbed,
            PlaneFrame {
                origin: [0., 0., -1.],
                ..Default::default()
            },
        ),
        (
            &plain,
            PlaneFrame {
                origin: [0., 0., 2.],
                u: [1., 0., 0.],
                v: [0., s, s],
            },
        ),
    ] {
        let error = project_drawing(
            &[DrawingSolid {
                solid,
                hatch: hatch(),
            }],
            DrawingProjection {
                frame,
                section_offset_mm: Some(0.),
                ..Default::default()
            },
        )
        .unwrap_err();
        assert_eq!(error.kind, GeometryErrorKind::Unsupported, "{error:?}");
    }
}

#[cfg(feature = "occt")]
#[test]
fn gear_sections_match_native_open_cascade() {
    let policy = GeometryPolicy::default();
    for (label, solid) in [
        (
            "spur",
            placed(Solid::rust_spur_gear_named(spur(3., None), EntityId::new(), policy).unwrap()),
        ),
        (
            "ring",
            placed(Solid::rust_internal_ring_named(ring(), EntityId::new(), policy).unwrap()),
        ),
    ] {
        let (solid, placement) = solid;
        for toward_cap in [true, false] {
            let frame = section_frame(placement, 1., toward_cap);
            let base = DrawingProjection {
                frame,
                section_offset_mm: Some(0.),
                ..Default::default()
            };
            let axial = DrawingProjection {
                frame,
                ..Default::default()
            };
            native_drawing::assert_matches_native(
                &[&solid],
                axial,
                &format!("{label} axial view toward cap {toward_cap}"),
            );
            // OpenCascade's section of the imported flank B-spline surfaces
            // misses short stretches of the boundary and bridges them with
            // chords whose midpoints lie up to 0.102 mm off the profile, in
            // both section-only and retained results. Its axial view of the
            // same import stays within 0.006 mm of the Rust rims. For these
            // straight bands that axial view is the section profile, so it is
            // the independent reference; the native section's own vertices
            // are also accepted, but its gap chords are not required to lie
            // on the Rust rims.
            for (view, options) in [
                ("section", base),
                (
                    "section only",
                    DrawingProjection {
                        section_only: true,
                        ..base
                    },
                ),
                (
                    "finite depth",
                    DrawingProjection {
                        section_depth_mm: Some(0.5),
                        ..base
                    },
                ),
            ] {
                native_drawing::assert_matches_native_with(
                    &[&solid],
                    options,
                    &[axial],
                    false,
                    &format!("{label} {view} toward cap {toward_cap}"),
                );
            }
        }
    }
}
