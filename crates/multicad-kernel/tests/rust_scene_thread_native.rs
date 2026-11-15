//! OpenCascade 7.9.3 references for thread drawings through the scene engine
//! (CAD-346): a terminal thread relief and a rounded thread, imported from
//! their own STEP export.
//!
//! Only cut faces are compared. OpenCascade's hidden-line removal of these
//! imported helicoidal B-spline faces is not a reference: measured in the
//! side view of the relief, it publishes the rear crest helices as visible and
//! no contour of the flanks at all, and in the rounded thread's side view it
//! draws the edge-on entry cap as hidden. Visibility of both families is
//! certified instead by marching the source solid's own containment test in
//! `drawing::scene::tests`.
#![cfg(feature = "occt")]
use multicad_core::{document::EntityId, drawings::VectorRole};
use multicad_kernel::{
    DrawingProjection, DrawingSolid, GeometryPolicy, HelicalThread, PlaneFrame, ProjectedVector,
    Solid, project_drawing, rust_kernel::rounded_thread::Spec,
};

fn relief() -> Solid {
    let root = 1. / (4. * 3_f64.sqrt());
    Solid::rust_helical_thread_end_relief(
        HelicalThread {
            major_radius_mm: 3.,
            pitch_mm: 1.,
            length_mm: 6.6,
            depth_mm: (1. - 0.125) / (2. * 30_f64.to_radians().tan()) - root,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: false,
        },
        Some(root),
        0.6,
        0.1,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn rounded(left_hand: bool) -> Solid {
    let p = 1.27;
    let a = 27.5_f64.to_radians();
    let r = p * a.cos() / (12. * (1. - a.sin()));
    Solid::rust_rounded_thread(
        Spec {
            major_radius_mm: 3.175,
            pitch_mm: p,
            length_mm: 3.7 * p,
            depth_mm: p / (3. * a.tan()),
            crest_radius_mm: r,
            root_radius_mm: r,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand,
        },
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn project(solid: &Solid, options: DrawingProjection) -> Vec<ProjectedVector> {
    project_drawing(&[DrawingSolid { solid, hatch: None }], options)
        .unwrap_or_else(|error| panic!("{error:?}"))
}

/// The largest distance from a preview point or chord midpoint of `from` to
/// the previews of `to`, restricted to the same role if required.
fn deviation(
    from: &[ProjectedVector],
    to: &[ProjectedVector],
    same_role: bool,
) -> (f64, [f64; 2], VectorRole) {
    let mut worst = (0., [0.; 2], VectorRole::Visible);
    for vector in from {
        let spans: Vec<_> = to
            .iter()
            .filter(|other| !same_role || other.role == vector.role)
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
                    worst = (distance, [x, y], vector.role);
                }
            }
        }
    }
    worst
}

fn cases(origin: [f64; 3]) -> Vec<(String, DrawingProjection)> {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let mut cases = Vec::new();
    for (name, [u, v]) in [
        ("axial", [[1., 0., 0.], [0., 1., 0.]]),
        ("side", [[1., 0., 0.], [0., 0., 1.]]),
        ("oblique", [[s, -s, 0.], [t, t, -2. * t]]),
    ] {
        for (suffix, offset) in [("", 0.), (" offset", 0.73)] {
            cases.push((
                format!("{name}{suffix} section-only"),
                DrawingProjection {
                    frame: PlaneFrame { origin, u, v },
                    section_offset_mm: Some(offset),
                    section_only: true,
                    ..Default::default()
                },
            ));
        }
    }
    cases
}

fn check(label: &str, solid: &Solid, origin: [f64; 3]) {
    let only = std::env::var("THREAD_CASE").ok();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("thread.step");
    solid.write_step(&path).unwrap();
    let native = Solid::read_step(&path, solid.geometry_policy()).unwrap();
    let mut failures = Vec::new();
    for (name, options) in cases(origin) {
        let label = format!("{label} {name}");
        if only
            .as_ref()
            .is_some_and(|only| !label.contains(only.as_str()))
        {
            continue;
        }
        let rust = project(solid, options);
        let reference = project(&native, options);
        let tangent = project(
            solid,
            DrawingProjection {
                tangent_edges: true,
                hidden_lines: true,
                ..options
            },
        );
        if let Ok(directory) = std::env::var("THREAD_DUMP") {
            for (who, list) in [("rust", &rust), ("native", &reference)] {
                let mut text = String::new();
                for vector in list.iter() {
                    let role = format!("{:?}", vector.role);
                    let points: Vec<String> = vector
                        .preview_mm
                        .iter()
                        .map(|p| format!("{} {}", p[0], p[1]))
                        .collect();
                    text.push_str(&format!("{role}|{}\n", points.join(",")));
                }
                std::fs::write(
                    format!("{directory}/{}-{who}.txt", label.replace(' ', "_")),
                    text,
                )
                .unwrap();
            }
        }
        let forward = deviation(&rust, &reference, true);
        let backward = deviation(&reference, &tangent, false);
        println!(
            "{{\"case\":\"{label}\",\"rust_to_native_mm\":{},\"native_to_rust_mm\":{},\"rust\":{},\"native\":{}}}",
            forward.0,
            backward.0,
            rust.len(),
            reference.len()
        );
        let bound = 3. * options.chord_tolerance_mm;
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
fn thread_end_relief_matches_native_open_cascade() {
    check("relief", &relief(), [0.2, 0.1, 3.3]);
}

#[test]
fn rounded_thread_matches_native_open_cascade() {
    check("rounded", &rounded(false), [0.2, 0.1, 2.35]);
    check("rounded left", &rounded(true), [0.2, 0.1, 2.35]);
}
