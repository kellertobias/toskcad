//! Non-axial, cropped and sectioned drawings of ZA worms and generated worm
//! wheels through the exact-curve scene engine (CAD-346).
use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
};
use multicad_kernel::{
    DrawingProjection, DrawingSolid, GeneratedZaWormWheel, GeometryErrorKind, GeometryPolicy,
    PlaneFrame, ProjectedVector, RigidTransform, Solid, ZaWorm, project_drawing,
    rust_kernel::{
        Rational,
        worm_gear::{Containment, Wheel, WheelSpec, Worm, WormSpec},
    },
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

fn stock(size: [i64; 3]) -> Solid {
    let r = |v: i64| Rational::new(v, 1).unwrap();
    Solid::rust_box(
        EntityId::new().uuid_bits(),
        size.map(r),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// Turn a solid about its axis so that no face seam falls on a silhouette
/// of the axis-aligned views, and move it.
#[cfg(feature = "occt")]
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

/// The provider solids of `worm` and `wheel`, built independently.
fn sources(bore: [f64; 2]) -> (Worm, Wheel) {
    (
        Worm::new(
            WormSpec {
                axial_module_mm: 1.,
                threaded_length_mm: 20.,
                bore_radius_mm: bore[0],
                axial_thinning_mm: 0.02,
            },
            1e-3,
        )
        .unwrap(),
        Wheel::new(
            WheelSpec {
                teeth: 32,
                axial_module_mm: 1.,
                face_width_mm: 1.6,
                bore_radius_mm: bore[1],
                pitch_thinning_mm: 0.03,
            },
            1e-3,
            1e7,
        )
        .unwrap(),
    )
}

fn section_planes() -> Vec<(&'static str, PlaneFrame, f64)> {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let frame = |u, v| PlaneFrame {
        origin: [0.2, -0.1, 0.15],
        u,
        v,
    };
    vec![
        ("side", frame([1., 0., 0.], [0., 0., 1.]), 0.),
        ("side offset", frame([1., 0., 0.], [0., 0., 1.]), 1.3),
        ("oblique", frame([s, -s, 0.], [t, t, -2. * t]), 0.4),
        ("tilted", frame([1., 0., 0.], [0., 0.96, 0.28]), 0.2),
    ]
}

/// The provider's exact membership predicate of one solid.
type Classifier = Box<dyn Fn([f64; 3]) -> Containment>;

/// The smallest probe radius of a doubling ladder at which the provider's
/// exact membership predicate finds both material and void around `p`,
/// probing 26 directions: an upper bound, for a locally smooth boundary, on
/// the distance from `p` to the boundary.
fn gap(classify: &dyn Fn([f64; 3]) -> Containment, p: [f64; 3]) -> f64 {
    if classify(p) == Containment::Boundary {
        return 0.;
    }
    let mut radius = 1e-9;
    while radius < 1. {
        let mut seen = [false; 2];
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let d = [dx, dy, dz].map(f64::from);
                    let length = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                    if length == 0. {
                        continue;
                    }
                    let q: [f64; 3] = std::array::from_fn(|i| p[i] + radius * d[i] / length);
                    match classify(q) {
                        Containment::Inside => seen[0] = true,
                        Containment::Outside => seen[1] = true,
                        Containment::Boundary => seen = [true; 2],
                    }
                }
            }
        }
        if seen == [true; 2] {
            return radius;
        }
        radius *= 2.;
    }
    f64::INFINITY
}

/// Every vertex of a published cut lies in the cutting plane and on the
/// boundary of the provider's solid, and the cut is not empty. The engine
/// joins collinear published lines of one source onto a common line within
/// its weld distance of four linear tolerances, so published vertices are
/// held to twice that distance; the members alone are exact to 1e-8 mm
/// (`scene/worm/tests.rs`).
#[test]
fn published_sections_lie_on_the_provider_boundary_in_the_cutting_plane() {
    let (screw, gear) = sources([2., 5.]);
    let weld = 4. * GeometryPolicy::default().linear_tolerance_mm;
    let bodies: [(&str, Solid, Classifier); 2] = [
        ("worm", worm(2.), Box::new(move |p| screw.classify(p, 0.))),
        ("wheel", wheel(5.), Box::new(move |p| gear.classify(p, 0.))),
    ];
    let mut worst = 0_f64;
    for (label, solid, classify) in &bodies {
        for (name, frame, offset) in section_planes() {
            let vectors = project(
                &[solid],
                DrawingProjection {
                    frame,
                    section_offset_mm: Some(offset),
                    section_only: true,
                    ..Default::default()
                },
            );
            let u = frame.u;
            let n = [
                u[1] * frame.v[2] - u[2] * frame.v[1],
                u[2] * frame.v[0] - u[0] * frame.v[2],
                u[0] * frame.v[1] - u[1] * frame.v[0],
            ];
            let v = [
                n[1] * u[2] - n[2] * u[1],
                n[2] * u[0] - n[0] * u[2],
                n[0] * u[1] - n[1] * u[0],
            ];
            let mut vertices = 0;
            for vector in &vectors {
                assert_eq!(vector.role, VectorRole::Section, "{label} {name}");
                // Lines carry every vertex; conics carry their two ends.
                let ends = [vector.preview_mm[0], *vector.preview_mm.last().unwrap()];
                let exact = matches!(vector.curve, VectorCurve::Line { .. });
                let points: &[[f64; 2]] = if exact { &vector.preview_mm } else { &ends };
                for [x, y] in points {
                    let world: [f64; 3] = std::array::from_fn(|i| {
                        frame.origin[i] + x * u[i] + y * v[i] + offset * n[i]
                    });
                    let distance = gap(classify.as_ref(), world);
                    worst = worst.max(distance);
                    assert!(
                        distance <= 2. * weld,
                        "{label} {name} {world:?} {distance} {:?}",
                        vector.curve
                    );
                    vertices += 1;
                }
            }
            assert!(vertices >= 8, "{label} {name} {vertices}");
        }
    }
    println!("worst published section vertex gap {worst} mm");
}

fn centred(solid: Solid, at: [f64; 3]) -> Solid {
    let b = solid.bounds().unwrap();
    solid
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: std::array::from_fn(|i| at[i] - (b.min[i] + b.max[i]) / 2.),
        })
        .unwrap()
}

/// A plate between the viewer and the gear hides every curve of the gear,
/// a plate behind it none, and crops bound what is published.
#[test]
fn a_plate_in_front_hides_the_gear_and_crops_bound_the_drawing() {
    let side = PlaneFrame {
        origin: [0.; 3],
        u: [1., 0., 0.],
        v: [0., 0., 1.],
    };
    for (label, solid, reach) in [("worm", worm(0.), 9.), ("wheel", wheel(5.), 20.)] {
        // The side view looks along +Y: the viewer is at -Y.
        let size = [2 * reach as i64 + 8, 2, 2 * reach as i64 + 8];
        let front = centred(stock(size), [0., -reach - 2., 0.]);
        let behind = centred(stock(size), [0., reach + 2., 0.]);
        let view = DrawingProjection {
            frame: side,
            hidden_lines: true,
            ..Default::default()
        };
        let alone = project(&[&solid], view);
        // Published length of the gear's curves of one role.
        let length = |vectors: &[ProjectedVector], source: usize, role: Option<VectorRole>| {
            vectors
                .iter()
                .filter(|v| v.source_index == source && role.is_none_or(|role| v.role == role))
                .flat_map(|v| v.preview_mm.windows(2))
                .map(|pair| (pair[1][0] - pair[0][0]).hypot(pair[1][1] - pair[0][1]))
                .sum::<f64>()
        };
        let same = |a: f64, b: f64| (a - b).abs() <= 1e-6 * a.max(b).max(1.);
        let visible = length(&alone, 0, Some(VectorRole::Visible));
        let hidden = length(&alone, 0, Some(VectorRole::Hidden));
        assert!(
            visible > reach && hidden > reach,
            "{label} {visible} {hidden}"
        );
        let covered = project(&[&front, &solid], view);
        assert_eq!(
            length(&covered, 1, Some(VectorRole::Visible)),
            0.,
            "{label}"
        );
        assert!(
            same(
                length(&covered, 1, Some(VectorRole::Hidden)),
                visible + hidden
            ),
            "{label}: every curve of the gear is hidden behind the plate"
        );
        let open = project(&[&behind, &solid], view);
        assert!(
            same(length(&open, 1, Some(VectorRole::Visible)), visible),
            "{label}"
        );
        assert!(
            same(length(&open, 1, Some(VectorRole::Hidden)), hidden),
            "{label}"
        );
        // A detail circle and a lateral window keep only what lies inside.
        let centre = [0.4 * reach, 0.2];
        let detail = project(
            &[&solid],
            DrawingProjection {
                detail: Some((centre, 3.)),
                ..view
            },
        );
        assert!(!detail.is_empty(), "{label}");
        for vector in &detail {
            for p in &vector.preview_mm {
                assert!(
                    (p[0] - centre[0]).hypot(p[1] - centre[1]) <= 3. + 1e-5,
                    "{label}"
                );
            }
        }
        let window = project(
            &[&solid],
            DrawingProjection {
                section_offset_mm: Some(0.),
                section_u_mm: Some([0.3 * reach, 0.3 * reach + 3.]),
                section_v_mm: Some([-0.5, 0.6]),
                ..view
            },
        );
        assert!(count(&window, VectorRole::Section) > 0, "{label}");
        for vector in &window {
            for p in &vector.preview_mm {
                assert!(
                    p[0] >= 0.3 * reach - 1e-5
                        && p[0] <= 0.3 * reach + 3. + 1e-5
                        && p[1] >= -0.5 - 1e-5
                        && p[1] <= 0.6 + 1e-5,
                    "{label} {p:?}"
                );
            }
        }
        // An exhausted budget is a typed failure without partial output.
        let inputs = [DrawingSolid {
            solid: &solid,
            hatch: None,
        }];
        let error = project_drawing(
            &inputs,
            DrawingProjection {
                maximum_curves: 5,
                ..view
            },
        )
        .unwrap_err();
        assert_eq!(error.kind, GeometryErrorKind::Limit, "{label}");
    }
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

    /// Diagnostic overlay: native curves wide and pale under the Rust curves.
    fn svg(directory: &str, label: &str, rust: &[ProjectedVector], native: &[ProjectedVector]) {
        let all = || rust.iter().chain(native).flat_map(|v| v.preview_mm.iter());
        let lo = [0, 1].map(|k| all().map(|p| p[k]).fold(f64::INFINITY, f64::min) - 0.5);
        let hi = [0, 1].map(|k| all().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max) + 0.5);
        let mut text = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{} {} {} {}\"><rect x=\"{}\" y=\"{}\" width=\"100%\" height=\"100%\" fill=\"white\"/>",
            lo[0],
            -hi[1],
            hi[0] - lo[0],
            hi[1] - lo[1],
            lo[0],
            -hi[1]
        );
        for (vectors, width, native) in [(native, 0.06, true), (rust, 0.015, false)] {
            for vector in vectors {
                let colour = match (native, vector.role) {
                    (true, VectorRole::Visible) => "#9cf",
                    (true, VectorRole::Hidden) => "#fc9",
                    (true, _) => "#9e9",
                    (false, VectorRole::Visible) => "#00c",
                    (false, VectorRole::Hidden) => "#c00",
                    (false, _) => "#070",
                };
                let points: Vec<String> = vector
                    .preview_mm
                    .iter()
                    .map(|p| format!("{:.4},{:.4}", p[0], -p[1]))
                    .collect();
                text.push_str(&format!(
                    "<polyline fill=\"none\" stroke=\"{colour}\" stroke-width=\"{width}\" points=\"{}\"/>",
                    points.join(" ")
                ));
            }
        }
        text.push_str("</svg>");
        let name = label.replace(' ', "-");
        std::fs::write(format!("{directory}/{name}.svg"), text).unwrap();
    }

    /// Every native curve lies on a Rust curve of the same source (any role,
    /// tangent seams and hidden lines included) within `bound`, except in
    /// the `uncovered` cases; in every section-only case the Rust and native
    /// cuts also lie on each other by role. Every other Rust-to-native
    /// deviation is printed, not asserted: on the STEP-exchanged worm the
    /// native hidden-line removal omits the helicoid flank silhouettes and
    /// draws the bore generators as visible through the core, which the
    /// member tests in `scene/worm/tests.rs` certify against the provider's
    /// exact surfaces instead.
    pub fn check(label: &str, solids: &[&Solid], origin: [f64; 3], bound: f64, uncovered: &[&str]) {
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
            if let Ok(directory) = std::env::var("SCENE_SVG") {
                svg(&directory, &label, &rust, &reference);
            }
            let forward = deviation(&rust, &reference, true);
            let backward = deviation(&reference, &complete, false);
            println!(
                "{{\"case\":\"{label}\",\"rust_to_native_mm\":{},\"native_to_rust_mm\":{},\"rust\":{},\"native\":{}}}",
                forward.0,
                backward.0,
                rust.len(),
                reference.len()
            );
            if name.ends_with("section-only") && forward.0 > bound {
                failures.push(format!("{label}: Rust {forward:?}"));
            }
            if !uncovered.contains(&name.as_str()) && backward.0 > bound {
                failures.push(format!("{label}: OpenCascade {backward:?}"));
            }
        }
        assert!(failures.is_empty(), "{failures:#?}");
    }

    /// Measured: the native axial detail publishes the hidden root circle
    /// of the worm beyond the 3 mm detail circle (1.76 mm off the Rust
    /// curves alone, 6.73 mm in the box scene).
    const DETAIL: &[&str] = &["axial detail"];

    #[test]
    fn worm_matches_native_open_cascade() {
        let solid = turned(worm(2.), [0.; 3]);
        check("worm", &[&solid], [0.3, -0.2, 0.4], 0.03, DETAIL);
    }

    #[test]
    fn worm_wheel_matches_native_open_cascade() {
        let solid = turned(wheel(5.), [0.; 3]);
        check("wheel", &[&solid], [0.3, -0.2, 0.1], 0.03, &[]);
    }

    #[test]
    fn box_and_worm_scenes_match_native_open_cascade() {
        let block = stock([8, 6, 5]);
        let screw = turned(worm(0.), [3., 2., 1.]);
        check("box worm", &[&block, &screw], [1., 1., 1.], 0.03, DETAIL);
    }

    #[test]
    fn box_and_wheel_scenes_match_native_open_cascade() {
        let block = stock([8, 6, 5]);
        let gear = turned(wheel(5.), [3., 2., 1.]);
        check("box wheel", &[&block, &gear], [1., 1., 1.], 0.03, &[]);
    }
}
