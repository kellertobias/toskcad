//! OpenCascade 7.9.3 references for thread drawings through the scene engine
//! (CAD-346): threads imported from their own STEP export, and a treated
//! internal thread bore rebuilt by an OpenCascade Boolean.
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

/// The solid's own STEP export, imported by OpenCascade.
fn imported(solid: &Solid) -> Solid {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("thread.step");
    solid.write_step(&path).unwrap();
    Solid::read_step(&path, solid.geometry_policy()).unwrap()
}

fn check(label: &str, solid: &Solid, origin: [f64; 3]) -> Vec<String> {
    compare(label, solid, &imported(solid), origin)
}

fn compare(label: &str, solid: &Solid, native: &Solid, origin: [f64; 3]) -> Vec<String> {
    let only = std::env::var("THREAD_CASE").ok();
    let mut failures = Vec::new();
    let mut unavailable = Vec::new();
    for (name, options) in cases(origin) {
        let label = format!("{label} {name}");
        if only
            .as_ref()
            .is_some_and(|only| !label.contains(only.as_str()))
        {
            continue;
        }
        let rust = project(solid, options);
        // OpenCascade cannot publish every curve family it computes; such a
        // case has no native reference and is reported, never skipped silently.
        let reference = match project_drawing(
            &[DrawingSolid {
                solid: native,
                hatch: None,
            }],
            options,
        ) {
            Ok(reference) => reference,
            Err(error) => {
                println!(
                    "{{\"case\":\"{label}\",\"native_unavailable\":\"{}\"}}",
                    error.message
                );
                unavailable.push(label.clone());
                continue;
            }
        };
        if reference.is_empty() && !rust.is_empty() {
            // A native cut that publishes nothing is no reference either.
            println!("{{\"case\":\"{label}\",\"native_unavailable\":\"no cut face published\"}}");
            unavailable.push(label.clone());
            continue;
        }
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
    unavailable
}

#[test]
fn thread_end_relief_matches_native_open_cascade() {
    assert!(check("relief", &relief(), [0.2, 0.1, 3.3]).is_empty());
}

fn rounded_bore() -> Solid {
    let p = 1.27;
    let a = 27.5_f64.to_radians();
    let r = p * a.cos() / (12. * (1. - a.sin()));
    Solid::rust_rounded_thread_bore(
        Spec {
            major_radius_mm: 3.175,
            pitch_mm: p,
            length_mm: 3.7 * p,
            depth_mm: p / (3. * a.tan()),
            crest_radius_mm: r,
            root_radius_mm: r,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand: false,
        },
        5.,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn chamfered(terminal: bool) -> Solid {
    let root = 1. / (4. * 3_f64.sqrt());
    Solid::rust_helical_thread_entry_chamfer(
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
        0.3,
        terminal.then_some(multicad_kernel::rust_kernel::thread_end_relief::Intent {
            width_mm: 0.6,
            radial_depth_mm: 0.1,
        }),
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

#[test]
fn thread_entry_chamfer_matches_native_open_cascade() {
    // Cuts through the cone zone and through the whole thread.
    let mut unavailable = check("chamfer", &chamfered(false), [0.2, 0.1, 0.15]);
    unavailable.extend(check("chamfer relief", &chamfered(true), [0.2, 0.1, 3.3]));
    // Measured: OpenCascade 7.9.3 cannot publish any side or oblique cut that
    // crosses the cone zone of the imported chamfer ("Projected curve family
    // is unsupported"): 7 of the 12 cuts. Every axial cut and the oblique cut
    // beside the cone are compared; the others are certified against the
    // source in drawing::scene::tests.
    assert!(
        unavailable.len() == 7
            && unavailable
                .iter()
                .all(|label| label.contains("side") || label.contains("oblique")),
        "{unavailable:?}"
    );
}

fn screw(socket: bool) -> Solid {
    use multicad_kernel::rust_kernel::headed_screw::{ScrewHead, ScrewSpec};
    let root = 1. / (4. * 3_f64.sqrt());
    Solid::rust_machine_screw(
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
        if socket {
            ScrewSpec {
                head: ScrewHead::SocketCap {
                    radius_mm: 4.5,
                    height_mm: 3.,
                    socket_across_flats_mm: 3.,
                    socket_depth_mm: 1.5,
                },
                shank_length_mm: 2.,
                shank_radius_mm: 3.,
            }
        } else {
            ScrewSpec {
                head: ScrewHead::ExternalHex {
                    across_flats_mm: 10.,
                    height_mm: 4.,
                },
                shank_length_mm: 0.,
                shank_radius_mm: 3.,
            }
        },
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

#[test]
fn machine_screws_match_native_open_cascade() {
    for (label, socket, origin) in [
        ("socket screw", true, [0.2, 0.1, 1.2]),
        ("hex screw", false, [0.2, 0.1, -0.8]),
    ] {
        let unavailable = check(label, &screw(socket), origin);
        assert!(unavailable.is_empty(), "{label}: {unavailable:?}");
    }
}

#[test]
fn rounded_thread_bore_matches_native_open_cascade() {
    assert!(check("bore", &rounded_bore(), [0.2, 0.1, 2.35]).is_empty());
}

#[test]
fn rounded_thread_matches_native_open_cascade() {
    assert!(check("rounded", &rounded(false), [0.2, 0.1, 2.35]).is_empty());
    assert!(check("rounded left", &rounded(true), [0.2, 0.1, 2.35]).is_empty());
}

/// A radius-4.5 authored stock less an ISO-like internal thread cavity with
/// an entry cone.
fn treated() -> Solid {
    treated_with(None)
}

/// As `treated`, with a terminal groove 1.2 mm wide and 0.4 mm deep.
fn grooved() -> Solid {
    treated_with(Some(
        multicad_kernel::rust_kernel::authored_internal_treatment::Terminal {
            width_mm: 1.2,
            radial_depth_mm: 0.4,
        },
    ))
}

fn treated_with(
    terminal: Option<multicad_kernel::rust_kernel::authored_internal_treatment::Terminal>,
) -> Solid {
    use multicad_kernel::rust_kernel::{Circle, Rational, cylinder::AxialCylinder};
    let policy = GeometryPolicy::default();
    let zero = Rational::new(0, 1).unwrap();
    let stock = AxialCylinder::new(
        7,
        Circle::new([zero; 2], Rational::new(9, 2).unwrap()).unwrap(),
        zero,
        Rational::new(13, 2).unwrap(),
    )
    .unwrap();
    Solid::from_axial_cylinder(stock, PlaneFrame::default(), policy)
        .unwrap()
        .rust_cut_authored_internal_thread_treatments_with_interrupt(
            multicad_kernel::InternalThreadTreatmentCut {
                profile: HelicalThread {
                    major_radius_mm: 3.,
                    pitch_mm: 1.,
                    length_mm: 6.5,
                    depth_mm: 5. * 3_f64.sqrt() / 16.,
                    crest_width_mm: 0.125,
                    lower_flank_deg: 30.,
                    upper_flank_deg: 30.,
                    radial_taper: 0.,
                    starts: 1,
                    left_hand: false,
                },
                intent: multicad_kernel::rust_kernel::authored_internal_treatment::Intent {
                    entry_depth_mm: Some(0.8),
                    terminal,
                },
                minimum_wall_mm: 0.5,
                frame: PlaneFrame::default(),
            },
            EntityId::new(),
            policy,
            None,
        )
        .unwrap()
}

/// The treated bore's OpenCascade reference, built without the treatment's
/// own export (which has none): the revolved stock less the entry cone and
/// a core inside the minor radius, less the imported external thread of the
/// same profile, one pitch longer at each end so no cap is coplanar.
fn treated_reference() -> Solid {
    reference(false)
}

/// The grooved bore's reference: the revolved stock also less the tail
/// cylinder `r < 3.4` past the shoulder at 5.3 mm; the same plug cuts the
/// thread, which the tail void already contains beyond the shoulder.
fn grooved_reference() -> Solid {
    reference(true)
}

fn reference(groove: bool) -> Solid {
    use multicad_kernel::{BooleanOperation, CurvedEdge2d, CurvedProfile, RevolutionAxis};
    let policy = GeometryPolicy::default();
    let minor = 3. - 5. * 3_f64.sqrt() / 16.;
    let entry = minor + 0.8;
    // The cone r = entry - z, continued inside the minor radius to the core.
    let core = minor - 0.2;
    let points: Vec<[f64; 2]> = if groove {
        vec![
            [entry, 0.],
            [4.5, 0.],
            [4.5, 6.5],
            [3.4, 6.5],
            [3.4, 5.3],
            [core, 5.3],
            [core, entry - core],
        ]
    } else {
        vec![
            [entry, 0.],
            [4.5, 0.],
            [4.5, 6.5],
            [core, 6.5],
            [core, entry - core],
        ]
    };
    let stock = Solid::revolve_curved_profile(
        &CurvedProfile {
            frame: PlaneFrame {
                origin: [0.; 3],
                u: [1., 0., 0.],
                v: [0., 0., 1.],
            },
            outer: (0..points.len())
                .map(|i| CurvedEdge2d::Line {
                    a: points[i],
                    b: points[(i + 1) % points.len()],
                })
                .collect(),
            holes: vec![],
        },
        RevolutionAxis {
            origin_mm: [0.; 3],
            direction: [0., 0., 1.],
        },
        std::f64::consts::TAU,
        policy,
    )
    .unwrap();
    let plug = Solid::rust_helical_thread(
        HelicalThread {
            major_radius_mm: 3.,
            pitch_mm: 1.,
            length_mm: 8.5,
            depth_mm: 5. * 3_f64.sqrt() / 16.,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: false,
        },
        None,
        EntityId::new(),
        PlaneFrame {
            origin: [0., 0., -1.],
            ..PlaneFrame::default()
        },
        policy,
    )
    .unwrap();
    stock
        .boolean(&imported(&plug), BooleanOperation::Subtract, policy)
        .unwrap()
        .unwrap()
}

#[test]
fn internal_thread_treatment_matches_native_open_cascade() {
    // Cuts through the entry cone zone and through the whole bore.
    let (solid, native) = (treated(), treated_reference());
    let mut unavailable = compare("treatment cone", &solid, &native, [0.2, 0.1, 0.4]);
    unavailable.extend(compare("treatment", &solid, &native, [0.2, 0.1, 3.25]));
    // Measured: as for the external chamfer, OpenCascade 7.9.3 cannot publish
    // a side or oblique cut that crosses the cone zone of its own Boolean
    // result ("Projected curve family is unsupported"): 7 of the 12 cuts.
    // Every axial cut, including one through the cone zone, and the oblique
    // cut beside the cone are compared; the others are certified against the
    // source in drawing::scene::tests.
    assert!(
        unavailable.len() == 7
            && unavailable
                .iter()
                .all(|label| label.contains("side") || label.contains("oblique")),
        "{unavailable:?}"
    );
}

#[test]
fn internal_thread_treatment_groove_matches_native_open_cascade() {
    // Cuts through the threaded bore and through the groove and shoulder.
    let (solid, native) = (grooved(), grooved_reference());
    let mut unavailable = compare("groove", &solid, &native, [0.2, 0.1, 3.25]);
    unavailable.extend(compare("groove tail", &solid, &native, [0.2, 0.1, 5.6]));
    println!("{{\"groove_native_unavailable\":{unavailable:?}}}");
    // Side and oblique cuts cross the entry cone zone, which OpenCascade
    // 7.9.3 cannot publish for its own Boolean result (as measured for the
    // treatment without a groove); every axial cut must be compared.
    assert!(
        unavailable
            .iter()
            .all(|label| label.contains("side") || label.contains("oblique")),
        "{unavailable:?}"
    );
}

/// The inventory's BSPT-like rounded pipe thread and an NPT-like flat one,
/// each imported from its own STEP export.
fn pipe(rounded: bool) -> Solid {
    use multicad_kernel::rust_kernel::axial_pipe_thread::authored::{
        AuthoredSection, AuthoredSpec,
    };
    let spec = if rounded {
        let p = 25.4 / 28.;
        AuthoredSpec {
            major_radius_mm: 9.728 / 2. - 4. / 32.,
            pitch_mm: p,
            length_mm: 6.5,
            radial_taper: 1. / 32.,
            section: AuthoredSection::BsptRounded {
                depth_mm: p / (3. * 27.5_f64.to_radians().tan()),
            },
            crest_radius_limit_mm: None,
        }
    } else {
        let p = 25.4 / 18.;
        AuthoredSpec {
            major_radius_mm: 6.8,
            pitch_mm: p,
            length_mm: 5.5,
            radial_taper: 1. / 32.,
            section: AuthoredSection::NptFlat {
                depth_mm: 0.8 * p,
                flat_width_mm: 0.038 * p,
            },
            crest_radius_limit_mm: None,
        }
    };
    Solid::rust_authored_axial_pipe_thread_with_interrupt(
        spec,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap()
}

#[test]
fn axial_pipe_threads_match_native_open_cascade() {
    let mut unavailable = check("bspt pipe", &pipe(true), [0.2, 0.1, 3.]);
    unavailable.extend(check("npt pipe", &pipe(false), [0.2, 0.1, 2.75]));
    // Measured: OpenCascade 7.9.3 publishes no cut face at all for the
    // imported NPT thread in the side plane (both offsets), while it does
    // for the axial and oblique cuts; those two are certified against the
    // source's radial bounds in drawing::scene::pipe_tests.
    assert_eq!(
        unavailable,
        [
            "npt pipe side section-only",
            "npt pipe side offset section-only"
        ],
    );
}

fn iso(length_mm: f64, left_hand: bool) -> HelicalThread {
    HelicalThread {
        major_radius_mm: 3.,
        pitch_mm: 1.,
        length_mm,
        depth_mm: 5. * 3_f64.sqrt() / 16.,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand,
    }
}

/// A plain helical thread with a rounded root.
fn plain(left_hand: bool) -> Solid {
    let root = 1. / (4. * 3_f64.sqrt());
    Solid::rust_helical_thread(
        HelicalThread {
            depth_mm: (1. - 0.125) / (2. * 30_f64.to_radians().tan()) - root,
            ..iso(6.6, left_hand)
        },
        Some(root),
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// A radius-5, 4 mm authored stock cylinder with an internal thread of
/// `length`: an overrunning tool, or a partial-depth cut with a floor.
fn internal(length: f64) -> Solid {
    use multicad_kernel::rust_kernel::{Circle, Rational, cylinder::AxialCylinder};
    let policy = GeometryPolicy::default();
    let zero = Rational::new(0, 1).unwrap();
    let stock = AxialCylinder::new(
        7,
        Circle::new([zero; 2], Rational::new(5, 1).unwrap()).unwrap(),
        zero,
        Rational::new(4, 1).unwrap(),
    )
    .unwrap();
    Solid::from_axial_cylinder(stock, PlaneFrame::default(), policy)
        .unwrap()
        .rust_cut_internal_thread(
            iso(length, false),
            0.2,
            EntityId::new(),
            PlaneFrame::default(),
            policy,
        )
        .unwrap()
}

fn nut(left_hand: bool) -> Solid {
    Solid::rust_hex_nut(
        HelicalThread {
            major_radius_mm: 1.55,
            pitch_mm: 0.5,
            length_mm: 2.4,
            depth_mm: 0.28,
            crest_width_mm: 0.08,
            ..iso(2.4, left_hand)
        },
        5.5,
        0.2,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

#[test]
fn plain_helical_threads_match_native_open_cascade() {
    assert!(check("plain", &plain(false), [0.2, 0.1, 3.3]).is_empty());
    assert!(check("plain left", &plain(true), [0.2, 0.1, 2.1]).is_empty());
}

#[test]
fn internal_threads_match_native_open_cascade() {
    for (label, length, origin) in [
        ("overrun internal", 5.2, [0.2, 0.1, 2.]),
        ("floor internal", 2.5, [0.2, 0.1, 1.2]),
        ("floor internal past void", 2.5, [0.2, 0.1, 3.2]),
    ] {
        let unavailable = check(label, &internal(length), origin);
        assert!(unavailable.is_empty(), "{label}: {unavailable:?}");
    }
}

#[test]
fn hex_nuts_match_native_open_cascade() {
    assert!(check("nut", &nut(false), [0.2, 0.1, 1.2]).is_empty());
    assert!(check("nut left", &nut(true), [0.2, 0.1, 0.9]).is_empty());
}
