//! Non-axial, cropped and sectioned drawings of ZA worms and generated worm
//! wheels through the exact-curve scene engine (CAD-346).
use multicad_core::{document::EntityId, drawings::VectorRole};
use multicad_kernel::{
    DrawingProjection, DrawingSolid, GeneratedZaWormWheel, GeometryPolicy, PlaneFrame,
    ProjectedVector, RigidTransform, Solid, ZaWorm, project_drawing, rust_kernel::Rational,
};

fn worm(bore: f64) -> Solid {
    Solid::rust_za_worm_named(
        ZaWorm {
            axial_module_mm: 1.,
            threaded_length_mm: 20.,
            bore_radius_mm: bore,
            axial_thinning_mm: 0.02,
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn wheel(bore: f64) -> Solid {
    Solid::rust_generated_za_worm_wheel_named(
        GeneratedZaWormWheel {
            teeth: 32,
            axial_module_mm: 1.,
            face_width_mm: 1.6,
            bore_radius_mm: bore,
            pitch_thinning_mm: 0.03,
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn stock() -> Solid {
    let r = |v: i64| Rational::new(v, 1).unwrap();
    Solid::rust_box(
        EntityId::new().uuid_bits(),
        [r(8), r(6), r(5)],
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// Turn a solid about its axis so that no face seam falls on a silhouette
/// of the axis-aligned views, and move it.
fn turned(solid: Solid, translation_mm: [f64; 3]) -> Solid {
    let (sin, cos) = 0.1_f64.sin_cos();
    solid
        .transformed(RigidTransform {
            rotation: [[cos, -sin, 0.], [sin, cos, 0.], [0., 0., 1.]],
            translation_mm,
        })
        .unwrap()
}

fn frames(origin: [f64; 3]) -> [(&'static str, PlaneFrame); 3] {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    [
        ("axial", [[1., 0., 0.], [0., 1., 0.]]),
        ("side", [[1., 0., 0.], [0., 0., 1.]]),
        ("oblique", [[s, -s, 0.], [t, t, -2. * t]]),
    ]
    .map(|(name, [u, v])| (name, PlaneFrame { origin, u, v }))
}

/// Plain and hidden-line views, then sections through `origin`: whole,
/// section-only, finite depth, lateral window and detail.
fn cases(origin: [f64; 3]) -> Vec<(String, DrawingProjection)> {
    let mut cases = Vec::new();
    for (name, frame) in frames(origin) {
        let view = DrawingProjection {
            frame,
            ..Default::default()
        };
        let section = DrawingProjection {
            section_offset_mm: Some(0.),
            ..view
        };
        cases.extend([
            (name.to_string(), view),
            (
                format!("{name} hidden"),
                DrawingProjection {
                    hidden_lines: true,
                    ..view
                },
            ),
            (
                format!("{name} detail"),
                DrawingProjection {
                    detail: Some(([0.5, 0.25], 3.)),
                    hidden_lines: true,
                    ..view
                },
            ),
            (format!("{name} section"), section),
            (
                format!("{name} section hidden"),
                DrawingProjection {
                    hidden_lines: true,
                    section_offset_mm: Some(0.3),
                    ..section
                },
            ),
            (
                format!("{name} section-only"),
                DrawingProjection {
                    section_only: true,
                    ..section
                },
            ),
            (
                format!("{name} finite depth"),
                DrawingProjection {
                    section_depth_mm: Some(1.),
                    ..section
                },
            ),
            (
                format!("{name} window"),
                DrawingProjection {
                    section_u_mm: Some([-4., 4.5]),
                    section_v_mm: Some([-1.5, 3.]),
                    section_depth_mm: Some(2.5),
                    ..section
                },
            ),
            (
                format!("{name} section detail"),
                DrawingProjection {
                    detail: Some(([1., 0.5], 4.)),
                    ..section
                },
            ),
        ]);
    }
    cases
}

fn project(solids: &[&Solid], options: DrawingProjection) -> Vec<ProjectedVector> {
    let inputs: Vec<_> = solids
        .iter()
        .map(|&solid| DrawingSolid { solid, hatch: None })
        .collect();
    project_drawing(&inputs, options).unwrap_or_else(|error| panic!("{error:?}"))
}

fn count(vectors: &[ProjectedVector], role: VectorRole) -> usize {
    vectors.iter().filter(|v| v.role == role).count()
}

#[test]
#[ignore = "timing report; run explicitly with --ignored --nocapture"]
fn report_worm_scene_drawings() {
    let only = std::env::var("SCENE_CASE").ok();
    let bodies = [("worm", worm(2.), [0.; 3]), ("wheel", wheel(5.), [0.; 3])];
    for (label, solid, origin) in &bodies {
        for (name, options) in cases(*origin) {
            let label = format!("{label} {name}");
            if only
                .as_ref()
                .is_some_and(|only| !label.contains(only.as_str()))
            {
                continue;
            }
            let start = std::time::Instant::now();
            let inputs = [DrawingSolid { solid, hatch: None }];
            let text = match project_drawing(&inputs, options) {
                Ok(v) => format!(
                    "ok v{} h{} s{}",
                    count(&v, VectorRole::Visible),
                    count(&v, VectorRole::Hidden),
                    count(&v, VectorRole::Section)
                ),
                Err(e) => format!("ERR {:?} {}", e.kind, e.message),
            };
            println!("{label}: {text} in {:?}", start.elapsed());
        }
    }
}

/// Independent OpenCascade 7.9.3 reference: the same solids exchanged
/// through STEP and drawn by native hidden-line removal and sectioning.
#[cfg(feature = "occt")]
mod native {
    use super::*;

    /// The largest distance from a preview point or chord midpoint of `from`
    /// to the previews of `to` with the same source and, if required, role.
    fn deviation(
        from: &[ProjectedVector],
        to: &[ProjectedVector],
        same_role: bool,
    ) -> (f64, [f64; 2], VectorRole, usize) {
        let mut worst = (0., [0.; 2], VectorRole::Visible, 0);
        for vector in from {
            let spans: Vec<_> = to
                .iter()
                .filter(|other| {
                    other.source_index == vector.source_index
                        && (!same_role || other.role == vector.role)
                })
                .flat_map(|other| other.preview_mm.windows(2).map(|w| [w[0], w[1]]))
                .collect();
            for pair in vector.preview_mm.windows(2) {
                let mid = [
                    pair[0][0].midpoint(pair[1][0]),
                    pair[0][1].midpoint(pair[1][1]),
                ];
                for [x, y] in [pair[0], mid, pair[1]] {
                    let distance = spans
                        .iter()
                        .map(|&[[a, b], [c, d]]| {
                            let [u, v] = [c - a, d - b];
                            let along = (((x - a) * u + (y - b) * v) / (u * u + v * v).max(1e-300))
                                .clamp(0., 1.);
                            (x - a - along * u).hypot(y - b - along * v)
                        })
                        .fold(f64::INFINITY, f64::min);
                    if distance > worst.0 {
                        worst = (distance, [x, y], vector.role, vector.source_index);
                    }
                }
            }
        }
        worst
    }

    /// Rust curves lie on native curves of the same role and source, and
    /// native curves on Rust curves of the same source, within `bound`.
    pub fn check(label: &str, solids: &[&Solid], origin: [f64; 3], bound: f64, skip: &[&str]) {
        let only = std::env::var("SCENE_CASE").ok();
        let directory = tempfile::tempdir().unwrap();
        let native: Vec<Solid> = solids
            .iter()
            .enumerate()
            .map(|(index, solid)| {
                let path = directory.path().join(format!("solid-{index}.step"));
                solid.write_step(&path).unwrap();
                Solid::read_step(&path, solid.geometry_policy()).unwrap()
            })
            .collect();
        let native: Vec<&Solid> = native.iter().collect();
        let mut failures = Vec::new();
        for (name, options) in cases(origin) {
            let label = format!("{label} {name}");
            if only
                .as_ref()
                .is_some_and(|only| !label.contains(only.as_str()))
                || skip.iter().any(|skipped| name == *skipped)
            {
                continue;
            }
            let rust = project(solids, options);
            let reference = project(&native, options);
            let complete = project(
                solids,
                DrawingProjection {
                    hidden_lines: true,
                    tangent_edges: true,
                    ..options
                },
            );
            let forward = deviation(&rust, &reference, true);
            let backward = deviation(&reference, &complete, false);
            println!(
                "{{\"case\":\"{label}\",\"rust_to_native_mm\":{},\"native_to_rust_mm\":{},\"rust\":{},\"native\":{}}}",
                forward.0,
                backward.0,
                rust.len(),
                reference.len()
            );
            if forward.0 > bound {
                failures.push(format!("{label}: Rust {forward:?}"));
            }
            if backward.0 > bound {
                failures.push(format!("{label}: OpenCascade {backward:?}"));
            }
        }
        assert!(failures.is_empty(), "{failures:#?}");
    }

    #[test]
    fn worm_matches_native_open_cascade() {
        let solid = turned(worm(2.), [0.; 3]);
        check("worm", &[&solid], [0.3, -0.2, 0.4], 0.03, &[]);
    }

    #[test]
    fn worm_wheel_matches_native_open_cascade() {
        let solid = turned(wheel(5.), [0.; 3]);
        check("wheel", &[&solid], [0.3, -0.2, 0.1], 0.03, &[]);
    }

    #[test]
    fn worm_scenes_match_native_open_cascade() {
        let block = stock();
        let screw = turned(worm(0.), [3., 2., 1.]);
        check("box worm", &[&block, &screw], [1., 1., 1.], 0.03, &[]);
        let gear = turned(wheel(5.), [3., 2., 1.]);
        check("box wheel", &[&block, &gear], [1., 1., 1.], 0.03, &[]);
    }
}
