//! Hidden-line drawings of a rounded-corner plate with axial through-holes,
//! the shape of an imported circuit board, among other bodies.
use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
    materials::HatchPattern,
    references::Origin,
};
use multicad_kernel::{
    CurvedEdge2d, CurvedProfile, DrawingHatch, DrawingProjection, DrawingSolid, GeometryErrorKind,
    GeometryPolicy, PlaneFrame, ProjectedVector, Solid, project_drawing, rust_kernel::Rational,
};
use std::f64::consts::FRAC_1_SQRT_2;

const HOLES: [[f64; 2]; 4] = [[8., 8.], [152., 8.], [8., 92.], [152., 92.]];
const TOL: f64 = 1e-6;

/// 160 x 100 x 1.6 with R3 corners and four r1.6 holes. `skew` moves each
/// corner centre away from its corner along the diagonal, as a three-point
/// arc with a rounded midpoint does.
fn plate(skew: f64) -> Solid {
    let line = |a, b| CurvedEdge2d::Line { a, b };
    let arc = |[x, y]: [f64; 2], [sx, sy]: [f64; 2], a, b| CurvedEdge2d::Arc {
        center: [x - sx * skew, y - sy * skew],
        a,
        b,
        clockwise: false,
    };
    let outer = vec![
        line([3., 0.], [157., 0.]),
        arc([157., 3.], [1., -1.], [157., 0.], [160., 3.]),
        line([160., 3.], [160., 97.]),
        arc([157., 97.], [1., 1.], [160., 97.], [157., 100.]),
        line([157., 100.], [3., 100.]),
        arc([3., 97.], [-1., 1.], [3., 100.], [0., 97.]),
        line([0., 97.], [0., 3.]),
        arc([3., 3.], [-1., -1.], [0., 3.], [3., 0.]),
    ];
    let holes: Vec<_> = HOLES
        .iter()
        .map(|&center| {
            vec![CurvedEdge2d::Circle {
                center,
                radius_mm: 1.6,
            }]
        })
        .collect();
    let origins: Vec<_> = (0..outer.len() + holes.len())
        .map(|_| Origin {
            source: EntityId::new(),
            role: "plate edge".into(),
        })
        .collect();
    let profile = CurvedProfile {
        frame: PlaneFrame::default(),
        outer,
        holes,
    };
    Solid::rust_extrude_normal_profile_named(
        &profile,
        1.6,
        EntityId::new(),
        &origins,
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn cube(id: u128, origin: [f64; 3], size: Rational) -> Solid {
    Solid::rust_box(
        id,
        [size; 3],
        PlaneFrame {
            origin,
            ..Default::default()
        },
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn top() -> DrawingProjection {
    DrawingProjection {
        hidden_lines: true,
        ..Default::default()
    }
}

fn front() -> DrawingProjection {
    DrawingProjection {
        frame: PlaneFrame {
            origin: [0.; 3],
            u: [1., 0., 0.],
            v: [0., 0., 1.],
        },
        hidden_lines: true,
        ..Default::default()
    }
}

fn project(solids: &[&Solid], options: DrawingProjection) -> Vec<ProjectedVector> {
    let scene: Vec<_> = solids
        .iter()
        .map(|&solid| DrawingSolid { solid, hatch: None })
        .collect();
    project_drawing(&scene, options).unwrap()
}

fn line_length(v: &ProjectedVector) -> f64 {
    match v.curve {
        VectorCurve::Line { a, b } => (b[0] - a[0]).hypot(b[1] - a[1]),
        _ => 0.,
    }
}

fn visible_length(vectors: &[ProjectedVector], source: usize) -> f64 {
    vectors
        .iter()
        .filter(|v| v.source_index == source && v.role == VectorRole::Visible)
        .map(line_length)
        .sum()
}

#[test]
fn bodies_under_a_hole_or_beyond_a_corner_round_stay_visible() {
    let half = Rational::new(1, 2).unwrap();
    let two = Rational::new(2, 1).unwrap();
    for skew in [0., 1.487e-7] {
        let plate = plate(skew);
        // Under hole (8, 8): every corner of the 2 mm cube is within 1.6 mm.
        let through_hole = cube(1, [7., 7., -10.], two);
        let under_material = cube(2, [50., 50., -10.], two);
        // Inside the corner square but outside the R3 round.
        let beyond_round = cube(3, [0.1, 0.1, -10.], half);
        // Inside the corner square and under the round.
        let under_round = cube(4, [1.5, 1.5, -10.], half);
        let vectors = project(
            &[
                &plate,
                &through_hole,
                &under_material,
                &beyond_round,
                &under_round,
            ],
            top(),
        );
        assert!(
            (visible_length(&vectors, 1) - 8.).abs() < TOL,
            "skew {skew}"
        );
        assert!(
            (visible_length(&vectors, 3) - 2.).abs() < TOL,
            "skew {skew}"
        );
        for source in [2, 4] {
            assert!(
                vectors.iter().any(|v| v.source_index == source),
                "skew {skew}"
            );
            assert_eq!(visible_length(&vectors, source), 0., "skew {skew}");
        }

        // The plate: four visible straight sides, visible corner arcs, a
        // visible near rim and a hidden far rim for every hole.
        let own: Vec<_> = vectors.iter().filter(|v| v.source_index == 0).collect();
        let circles = |role| {
            own.iter()
                .filter(|v| {
                    v.role == role
                        && matches!(v.curve, VectorCurve::Ellipse { major_mm, minor_mm, .. }
                            if (major_mm - 1.6).abs() < TOL && (minor_mm - 1.6).abs() < TOL)
                })
                .count()
        };
        assert_eq!(circles(VectorRole::Visible), 4);
        assert_eq!(circles(VectorRole::Hidden), 4);
        let straight: f64 = own
            .iter()
            .filter(|v| v.role == VectorRole::Visible)
            .map(|v| line_length(v))
            .sum();
        assert!((straight - 2. * (154. + 94.)).abs() < TOL, "{straight}");
    }
}

#[test]
fn transverse_view_shows_the_plate_band_and_hides_the_hole_walls() {
    let plate = plate(1.487e-7);
    // A separate body makes this a mixed scene rather than a lone panel.
    let beside = cube(1, [200., 50., 0.], Rational::new(2, 1).unwrap());
    let vectors: Vec<_> = project(&[&plate, &beside], front())
        .into_iter()
        .filter(|v| v.source_index == 0)
        .collect();
    let visible: Vec<_> = vectors
        .iter()
        .filter(|v| v.role == VectorRole::Visible)
        .collect();
    // The band outline is 160 x 1.6.
    let horizontal: f64 = visible
        .iter()
        .filter_map(|v| match v.curve {
            VectorCurve::Line { a, b } if (a[1] - b[1]).abs() < TOL => Some((b[0] - a[0]).abs()),
            _ => None,
        })
        .sum();
    assert!(horizontal >= 2. * 160. - TOL, "{horizontal}");
    for x in [0., 160.] {
        assert!(visible.iter().any(|v| matches!(v.curve,
            VectorCurve::Line { a, b } if (a[0] - x).abs() < TOL && (b[0] - x).abs() < TOL
                && (a[1] - b[1]).abs() > 1.6 - TOL)));
    }
    // Every hole wall silhouette is hidden at x = 8 ± 1.6 and 152 ± 1.6.
    for x in [6.4, 9.6, 150.4, 153.6] {
        let wall = |v: &&ProjectedVector| {
            matches!(v.curve, VectorCurve::Line { a, b }
                if (a[0] - x).abs() < TOL && (b[0] - x).abs() < TOL)
        };
        assert!(
            vectors
                .iter()
                .filter(wall)
                .any(|v| v.role == VectorRole::Hidden),
            "x = {x}"
        );
        assert!(
            !vectors
                .iter()
                .filter(wall)
                .any(|v| v.role == VectorRole::Visible)
        );
    }
}

#[test]
fn oblique_views_use_the_sweep_shell_and_unbounded_corner_errors_are_refused() {
    let board = plate(1.487e-7);
    let box_behind = cube(1, [50., 50., -10.], Rational::new(2, 1).unwrap());
    // An oblique view never reaches the axial/transverse decomposition: the
    // exact sweep-shell elevation draws it.
    let oblique = DrawingProjection {
        frame: PlaneFrame {
            origin: [0.; 3],
            u: [1., 0., 0.],
            v: [0., FRAC_1_SQRT_2, FRAC_1_SQRT_2],
        },
        hidden_lines: true,
        ..Default::default()
    };
    let vectors = project(&[&board, &box_behind], oblique);
    assert!(vectors.iter().any(|v| v.source_index == 0
        && v.role == VectorRole::Visible
        && matches!(v.curve, VectorCurve::Line { a, b }
            if a[1].abs() < TOL && b[1].abs() < TOL && (a[0] - b[0]).abs() > 150.)));

    // A corner centre 1e-6 mm off its corner strays 4.1e-7 mm from the
    // exact quarter, beyond the default 1e-7 mm linear tolerance.
    let inexact = plate(1e-6);
    let scene = [&inexact, &box_behind].map(|solid| DrawingSolid { solid, hatch: None });
    let error = project_drawing(&scene, top()).unwrap_err();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert!(
        error.message.contains("cardinal quarter rounds"),
        "{}",
        error.message
    );
}

/// 100 x 60 x 10 with one R20 corner, a 15 mm rounded slot of radius 5
/// along x and an r10 hole: the shared-workflow bracket.
fn bracket() -> Solid {
    let line = |a, b| CurvedEdge2d::Line { a, b };
    let outer = vec![
        CurvedEdge2d::Arc {
            center: [80., 40.],
            a: [100., 40.],
            b: [80., 60.],
            clockwise: false,
        },
        line([80., 60.], [0., 60.]),
        line([0., 60.], [0., 0.]),
        line([0., 0.], [100., 0.]),
        line([100., 0.], [100., 40.]),
    ];
    let slot = vec![
        line([60., 35.], [75., 35.]),
        CurvedEdge2d::Arc {
            center: [75., 30.],
            a: [75., 35.],
            b: [75., 25.],
            clockwise: true,
        },
        line([75., 25.], [60., 25.]),
        CurvedEdge2d::Arc {
            center: [60., 30.],
            a: [60., 25.],
            b: [60., 35.],
            clockwise: true,
        },
    ];
    let hole = vec![CurvedEdge2d::Circle {
        center: [25., 30.],
        radius_mm: 10.,
    }];
    let origins: Vec<_> = (0..outer.len() + slot.len() + hole.len())
        .map(|_| Origin {
            source: EntityId::new(),
            role: "bracket edge".into(),
        })
        .collect();
    Solid::rust_extrude_normal_profile_named(
        &CurvedProfile {
            frame: PlaneFrame::default(),
            outer,
            holes: vec![slot, hole],
        },
        10.,
        EntityId::new(),
        &origins,
        GeometryPolicy::default(),
    )
    .unwrap()
}

#[test]
fn a_through_slot_reveals_what_lies_behind_its_whole_stadium() {
    let bracket = bracket();
    let one = Rational::new(1, 1).unwrap();
    // Under the slot's flat middle, under its semicircular end, and under
    // material just beside the slot.
    let middle = cube(1, [67., 29.5, -10.], one);
    let end = cube(2, [77.5, 29.5, -10.], one);
    let beside = cube(3, [67., 36., -10.], one);
    let vectors = project(&[&bracket, &middle, &end, &beside], top());
    assert!((visible_length(&vectors, 1) - 4.).abs() < TOL);
    assert!((visible_length(&vectors, 2) - 4.).abs() < TOL);
    assert!(vectors.iter().any(|v| v.source_index == 3));
    assert_eq!(visible_length(&vectors, 3), 0.);

    // Side-on along y the slot walls are hidden at its two ends, x = 55
    // and x = 80, and at the r10 hole's silhouettes, x = 15 and x = 35.
    let vectors: Vec<_> = project(&[&bracket, &middle], front())
        .into_iter()
        .filter(|v| v.source_index == 0)
        .collect();
    for x in [15., 35., 55., 80.] {
        let wall = |v: &&ProjectedVector| {
            matches!(v.curve, VectorCurve::Line { a, b }
                if (a[0] - x).abs() < TOL && (b[0] - x).abs() < TOL
                    && (a[1] - b[1]).abs() > 10. - TOL)
        };
        assert!(
            vectors
                .iter()
                .filter(wall)
                .any(|v| v.role == VectorRole::Hidden),
            "x = {x}"
        );
        assert!(
            !vectors
                .iter()
                .filter(wall)
                .any(|v| v.role == VectorRole::Visible)
        );
    }
    // Looking along x the slot's flat walls are edge-on at y = 25 and 35.
    let side = DrawingProjection {
        frame: PlaneFrame {
            origin: [0.; 3],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        },
        hidden_lines: true,
        ..Default::default()
    };
    let vectors = project(&[&bracket, &middle], side);
    for y in [25., 35.] {
        assert!(
            vectors.iter().any(|v| v.source_index == 0
                && v.role == VectorRole::Hidden
                && matches!(v.curve, VectorCurve::Line { a, b }
                if (a[0] - y).abs() < TOL && (b[0] - y).abs() < TOL)),
            "y = {y}"
        );
    }
}

#[test]
fn bracket_sections_cut_the_exact_sweep_shell() {
    let bracket = bracket();
    let section = |offset: f64, depth: Option<f64>| {
        project_drawing(
            &[DrawingSolid {
                solid: &bracket,
                hatch: None,
            }],
            DrawingProjection {
                section_offset_mm: Some(offset),
                section_only: depth.is_none(),
                section_depth_mm: depth,
                ..front()
            },
        )
    };
    // At y = 31 the cut crosses the r10 hole at x = 25 ± sqrt(99) and the
    // slot's semicircular ends at x = 60 - sqrt(24) and 75 + sqrt(24).
    let walls = [
        25. - 99_f64.sqrt(),
        25. + 99_f64.sqrt(),
        60. - 24_f64.sqrt(),
        75. + 24_f64.sqrt(),
    ];
    for depth in [None, Some(20.)] {
        let vectors = section(-31., depth).unwrap();
        for x in walls {
            assert!(
                vectors.iter().any(|v| v.role == VectorRole::Section
                    && matches!(v.curve, VectorCurve::Line { a, b }
                        if (a[0] - x).abs() < TOL && (b[0] - x).abs() < TOL
                            && (a[1] - b[1]).abs() > 10. - TOL)),
                "x = {x}, depth {depth:?}"
            );
        }
    }
    // At y = 30 the plane holds the hole's seam edge. The exact shell
    // section refuses a coplanar edge, so the retained view stays a typed
    // refusal rather than a guessed drawing.
    let error = section(-30., Some(20.)).unwrap_err();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
}

#[test]
fn a_finely_hatched_plate_section_stays_within_the_hatch_budget() {
    // The hatch boundary follows the exact section within the drawing chord
    // tolerance. Sampling the small holes at a quarter of the minimum
    // feature size gave each one thousands of segments and overran the
    // shared hatch budget for a 0.05 mm hatch across the whole plate.
    let board = plate(0.);
    let vectors = project_drawing(
        &[DrawingSolid {
            solid: &board,
            hatch: Some(DrawingHatch {
                angle_rad: std::f64::consts::FRAC_PI_4,
                spacing_mm: 0.05,
                line_width_mm: 0.025,
                pattern: HatchPattern::Lines,
            }),
        }],
        DrawingProjection {
            section_offset_mm: Some(0.8),
            section_only: true,
            ..Default::default()
        },
    )
    .unwrap();
    let hatch: Vec<_> = vectors
        .iter()
        .filter(|v| v.role == VectorRole::Hatch)
        .collect();
    // 45 degree lines 0.05 mm apart cross a 160 x 100 mm plate about
    // (160 + 100) / sqrt(2) / 0.05 = 3677 times; holes split a few more.
    assert!(hatch.len() > 3_600, "{}", hatch.len());
    // Holes are sampled by tangent segments, so no stroke enters a hole.
    for vector in hatch {
        let VectorCurve::Line { a, b } = vector.curve else {
            panic!("hatch {:?}", vector.curve);
        };
        for t in [0., 0.25, 0.5, 0.75, 1.] {
            let point: [f64; 2] = std::array::from_fn(|i| a[i] + (b[i] - a[i]) * t);
            for center in HOLES {
                let distance = (point[0] - center[0]).hypot(point[1] - center[1]);
                assert!(
                    distance > 1.6 - TOL,
                    "{point:?} inside the hole at {center:?}"
                );
            }
        }
    }
}
