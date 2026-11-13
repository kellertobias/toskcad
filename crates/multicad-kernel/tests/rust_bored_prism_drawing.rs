//! Rust drawings of a rectangular box with a drilled through bore.
//!
//! `box − cylinder` keeps a rational prism with an inner loop. The drawing
//! providers receive the identical point set as orthogonal stock with a
//! through bore. Expectations are analytic: the box spans [0,20]×[0,12]×[0,8]
//! and the bore of radius 2 runs along +Z through (10, 6).
use multicad_core::drawings::{VectorCurve, VectorRole};
use multicad_kernel::{
    BooleanOperation, DrawingHatch, DrawingProjection, DrawingSolid, GeometryErrorKind,
    GeometryPolicy, PlaneFrame, ProjectedVector, Solid, project_drawing,
};

#[cfg(feature = "occt")]
mod native_drawing;
#[cfg(feature = "occt")]
use native_drawing::assert_matches_native;

const TOLERANCE: f64 = 1e-9;

fn bored() -> Solid {
    use multicad_core::document::EntityId;
    use multicad_rust_kernel::{Circle, Rational, cylinder::AxialCylinder};
    let policy = GeometryPolicy::default();
    let n = |value: i64| Rational::new(value, 1).unwrap();
    let stock = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [20, 12, 8].map(n),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    // The same exact descriptor as `Solid::cylinder` in a Rust-only build,
    // spelled out so an OpenCascade build also receives Rust operands.
    let tool = EntityId::new();
    let drill = Solid::from_axial_cylinder_sweep(
        AxialCylinder::new(
            tool.uuid_bits(),
            Circle::new([n(10), n(6)], n(2)).unwrap(),
            n(-1),
            n(10),
        )
        .unwrap(),
        PlaneFrame {
            origin: [10., 6., -1.],
            ..Default::default()
        },
        [0.; 2],
        false,
        multicad_core::references::Origin {
            source: tool,
            role: "primitive.circle".into(),
        },
        policy,
    )
    .unwrap();
    let result = stock
        .boolean(&drill, BooleanOperation::Subtract, policy)
        .unwrap()
        .unwrap();
    assert!(
        result.rust_prism().is_some(),
        "box − cylinder keeps a prism"
    );
    result
}

fn hatch() -> Option<DrawingHatch> {
    Some(DrawingHatch {
        angle_rad: 0.6,
        spacing_mm: 0.5,
        line_width_mm: 0.18,
        pattern: multicad_core::materials::HatchPattern::Lines,
    })
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

fn frame(origin: [f64; 3], u: [f64; 3], v: [f64; 3]) -> PlaneFrame {
    PlaneFrame { origin, u, v }
}

fn lines(vectors: &[ProjectedVector], role: VectorRole) -> Vec<[[f64; 2]; 2]> {
    let mut lines: Vec<_> = vectors
        .iter()
        .filter(|v| v.role == role)
        .filter_map(|v| match v.curve {
            VectorCurve::Line { a, b } => Some(if a <= b { [a, b] } else { [b, a] }),
            _ => None,
        })
        .collect();
    lines.sort_by(|p, q| p.partial_cmp(q).unwrap());
    lines
}

fn circles(vectors: &[ProjectedVector], role: VectorRole) -> Vec<([f64; 2], f64, f64)> {
    vectors
        .iter()
        .filter(|v| v.role == role)
        .filter_map(|v| match v.curve {
            VectorCurve::Ellipse {
                center,
                major_mm,
                minor_mm,
                start_rad,
                end_rad,
                ..
            } => Some((center, major_mm.max(minor_mm), (end_rad - start_rad).abs())),
            _ => None,
        })
        .collect()
}

fn near(a: [f64; 2], b: [f64; 2]) -> bool {
    (a[0] - b[0]).abs() < TOLERANCE && (a[1] - b[1]).abs() < TOLERANCE
}

fn same_lines(actual: &[[[f64; 2]; 2]], expected: &[[[f64; 2]; 2]]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?}");
    for line in expected {
        let (a, b) = (line[0], line[1]);
        assert!(
            actual
                .iter()
                .any(|l| (near(l[0], a) && near(l[1], b)) || (near(l[0], b) && near(l[1], a))),
            "missing {line:?} in {actual:?}"
        );
    }
}

fn rectangle([x0, y0]: [f64; 2], [x1, y1]: [f64; 2]) -> Vec<[[f64; 2]; 2]> {
    vec![
        [[x0, y0], [x1, y0]],
        [[x1, y0], [x1, y1]],
        [[x0, y1], [x1, y1]],
        [[x0, y0], [x0, y1]],
    ]
}

#[test]
fn axial_view_publishes_the_outline_and_one_full_bore_rim() {
    let solid = bored();
    for (label, view) in [
        ("toward +Z", frame([0.; 3], [1., 0., 0.], [0., 1., 0.])),
        ("toward -Z", frame([0.; 3], [1., 0., 0.], [0., -1., 0.])),
    ] {
        let flip = if label == "toward -Z" { -1. } else { 1. };
        let vectors = draw(
            &solid,
            DrawingProjection {
                frame: view,
                ..Default::default()
            },
        );
        same_lines(
            &lines(&vectors, VectorRole::Visible),
            &rectangle([0., 0.], [20., 12.])
                .into_iter()
                .map(|l| l.map(|[x, y]| [x, flip * y]))
                .collect::<Vec<_>>(),
        );
        let rims = circles(&vectors, VectorRole::Visible);
        assert_eq!(rims.len(), 1, "{label}: {rims:?}");
        let (center, radius, sweep) = rims[0];
        assert!(near(center, [10., flip * 6.]), "{label}");
        assert!((radius - 2.).abs() < TOLERANCE && (sweep - std::f64::consts::TAU).abs() < 1e-12);
        assert!(vectors.iter().all(|v| v.role == VectorRole::Visible));

        // The far rim lies exactly behind the near rim; hidden lines add only
        // the rear outline and rim, coincident in this direction.
        let hidden = draw(
            &solid,
            DrawingProjection {
                frame: view,
                hidden_lines: true,
                ..Default::default()
            },
        );
        assert_eq!(
            lines(&hidden, VectorRole::Visible),
            lines(&vectors, VectorRole::Visible)
        );
        for (center, radius, _) in circles(&hidden, VectorRole::Hidden) {
            assert!(near(center, [10., flip * 6.]) && (radius - 2.).abs() < TOLERANCE);
        }
    }
}

#[test]
fn side_view_hides_both_bore_silhouettes_only_when_requested() {
    let solid = bored();
    let side = frame([0.; 3], [1., 0., 0.], [0., 0., 1.]);
    let visible = draw(
        &solid,
        DrawingProjection {
            frame: side,
            ..Default::default()
        },
    );
    same_lines(
        &lines(&visible, VectorRole::Visible),
        &rectangle([0., 0.], [20., 8.]),
    );
    assert!(lines(&visible, VectorRole::Hidden).is_empty());

    let hidden = draw(
        &solid,
        DrawingProjection {
            frame: side,
            hidden_lines: true,
            ..Default::default()
        },
    );
    same_lines(
        &lines(&hidden, VectorRole::Visible),
        &rectangle([0., 0.], [20., 8.]),
    );
    // Rear outline edges coincide with the visible outline, as for any box;
    // the interior hidden lines are exactly the two bore silhouettes.
    let silhouettes: Vec<_> = lines(&hidden, VectorRole::Hidden)
        .into_iter()
        .filter(|[a, b]| {
            (a[0] - b[0]).abs() < TOLERANCE && a[0] > TOLERANCE && a[0] < 20. - TOLERANCE
        })
        .collect();
    same_lines(
        &silhouettes,
        &[[[8., 0.], [8., 8.]], [[12., 0.], [12., 8.]]],
    );
}

#[test]
fn oblique_view_projects_each_bore_rim_as_an_exact_ellipse() {
    let solid = bored();
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let view = frame([0.; 3], [s, -s, 0.], [t, t, -2. * t]);
    let vectors = draw(
        &solid,
        DrawingProjection {
            frame: view,
            hidden_lines: true,
            ..Default::default()
        },
    );
    // The view normal u × v = -(1,1,1)/√3 makes |n·z| = 1/√3 with each rim.
    let mut ellipses = 0;
    for vector in &vectors {
        if let VectorCurve::Ellipse {
            major_mm, minor_mm, ..
        } = vector.curve
        {
            ellipses += 1;
            assert!((major_mm - 2.).abs() < TOLERANCE, "{vector:?}");
            assert!(
                (minor_mm - 2. / 3_f64.sqrt()).abs() < TOLERANCE,
                "{vector:?}"
            );
        }
    }
    assert!(
        ellipses >= 2,
        "both rims publish exact ellipse arcs: {vectors:?}"
    );
}

#[test]
fn axial_sections_keep_the_bore_void_and_hatch_only_material() {
    let solid = bored();
    let plane = frame([10., 6., 4.], [1., 0., 0.], [0., 1., 0.]);
    let base = DrawingProjection {
        frame: plane,
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
                section_depth_mm: Some(1.),
                ..base
            },
        ),
    ] {
        let vectors = draw(&solid, options);
        same_lines(
            &lines(&vectors, VectorRole::Section),
            &rectangle([-10., -6.], [10., 6.]),
        );
        let rims = circles(&vectors, VectorRole::Section);
        assert_eq!(rims.len(), 1, "{label}");
        assert!(near(rims[0].0, [0., 0.]) && (rims[0].1 - 2.).abs() < TOLERANCE);
        let strokes: Vec<_> = vectors
            .iter()
            .filter(|v| v.role == VectorRole::Hatch)
            .collect();
        assert!(strokes.len() > 20, "{label}: material must be hatched");
        for stroke in strokes {
            let VectorCurve::Line { a, b } = stroke.curve else {
                panic!("hatch strokes are lines");
            };
            for q in [0., 0.25, 0.5, 0.75, 1.] {
                let p = [a[0] + (b[0] - a[0]) * q, a[1] + (b[1] - a[1]) * q];
                assert!(
                    p[0].abs() <= 10. + 1e-9 && p[1].abs() <= 6. + 1e-9,
                    "{label}: {p:?}"
                );
                assert!(
                    p[0].hypot(p[1]) >= 2. - 1e-9,
                    "{label}: hatch enters the bore at {p:?}"
                );
            }
        }
        // Behind the cut, only outlines coincident with the section remain.
        for line in lines(&vectors, VectorRole::Visible) {
            assert!(
                rectangle([-10., -6.], [10., 6.])
                    .iter()
                    .any(|r| near(r[0], line[0]) && near(r[1], line[1])),
                "{label}: {line:?}"
            );
        }
        if label == "section only" {
            assert!(lines(&vectors, VectorRole::Visible).is_empty());
        }
    }
}

#[test]
fn side_sections_split_the_cut_face_at_the_bore() {
    let solid = bored();
    // The plane y = 6 contains the bore axis.
    let plane = frame([0., 6., 0.], [1., 0., 0.], [0., 0., 1.]);
    let base = DrawingProjection {
        frame: plane,
        section_offset_mm: Some(0.),
        ..Default::default()
    };
    let mut expected = rectangle([0., 0.], [8., 8.]);
    expected.extend(rectangle([12., 0.], [20., 8.]));
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
                section_depth_mm: Some(1.),
                ..base
            },
        ),
    ] {
        let vectors = draw(&solid, options);
        same_lines(&lines(&vectors, VectorRole::Section), &expected);
        for stroke in vectors.iter().filter(|v| v.role == VectorRole::Hatch) {
            let VectorCurve::Line { a, b } = stroke.curve else {
                panic!("hatch strokes are lines");
            };
            let mid = (a[0] + b[0]) / 2.;
            assert!(
                !(8. + 1e-9..12. - 1e-9).contains(&mid),
                "{label}: hatch inside the bore"
            );
        }
    }
    // A lateral window keeps only the left wall's cut face.
    let window = draw(
        &solid,
        DrawingProjection {
            section_u_mm: Some([1., 7.]),
            section_v_mm: Some([1., 7.]),
            ..base
        },
    );
    same_lines(
        &lines(&window, VectorRole::Section),
        &rectangle([1., 1.], [7., 7.]),
    );
}

#[test]
fn a_part_behind_the_bore_is_visible_only_through_the_opening() {
    let policy = GeometryPolicy::default();
    let solid = bored();
    let block = |size: [i64; 3]| {
        use multicad_rust_kernel::Rational;
        Solid::rust_box(
            multicad_core::document::EntityId::new().uuid_bits(),
            size.map(|value| Rational::new(value, 1).unwrap()),
            PlaneFrame::default(),
            policy,
        )
        .unwrap()
    };
    let behind = block([20, 12, 1])
        .transformed(multicad_kernel::RigidTransform {
            rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: [0., 0., -3.],
        })
        .unwrap();
    // A plate edge at x = 10 crosses the bore below the stock.
    let half = block([10, 12, 1])
        .transformed(multicad_kernel::RigidTransform {
            rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: [10., 0., -3.],
        })
        .unwrap();
    for (scene, crossing) in [(&behind, false), (&half, true)] {
        for options in [
            DrawingProjection::default(),
            DrawingProjection {
                frame: frame([10., 6., 4.], [1., 0., 0.], [0., 1., 0.]),
                section_offset_mm: Some(0.),
                ..Default::default()
            },
        ] {
            let vectors = project_drawing(
                &[
                    DrawingSolid {
                        solid: &solid,
                        hatch: hatch(),
                    },
                    DrawingSolid {
                        solid: scene,
                        hatch: hatch(),
                    },
                ],
                options,
            )
            .unwrap();
            let offset = options.frame.origin;
            let through: Vec<_> = vectors
                .iter()
                .filter(|v| v.source_index == 1 && v.role == VectorRole::Visible)
                .collect();
            if crossing {
                // Only the plate edge segment inside the bore is visible.
                let [segment] = through.as_slice() else {
                    panic!("{through:?}");
                };
                let VectorCurve::Line { a, b } = segment.curve else {
                    panic!("{segment:?}");
                };
                let [a, b] = [a, b].map(|p| [p[0] + offset[0], p[1] + offset[1]]);
                assert!((a[0] - 10.).abs() < TOLERANCE && (b[0] - 10.).abs() < TOLERANCE);
                let [lo, hi] = [a[1].min(b[1]), a[1].max(b[1])];
                assert!(
                    (lo - 4.).abs() < 1e-7 && (hi - 8.).abs() < 1e-7,
                    "{a:?} {b:?}"
                );
            } else {
                // The plate's outline coincides with the stock outline, which
                // the nearer stock owns. Nothing of it shows through the bore.
                assert!(through.is_empty(), "{through:?}");
            }
        }
    }
}

#[test]
fn oblique_sections_only_keep_the_exact_elliptical_void() {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let solid = bored();
    // Planes through the axis at several heights. Near the stock caps the
    // ellipse reaches a cap face and its chord there bounds no material.
    for (height, window) in [
        (4., None),
        (1.5, None),
        (0.5, None),
        (6.5, None),
        (4., Some(([-3., 6.], [-4., 2.]))),
    ] {
        let options = DrawingProjection {
            frame: frame([10., 6., height], [s, -s, 0.], [t, t, -2. * t]),
            section_offset_mm: Some(0.),
            section_only: true,
            section_u_mm: window.map(|w: ([f64; 2], [f64; 2])| w.0),
            section_v_mm: window.map(|w: ([f64; 2], [f64; 2])| w.1),
            ..Default::default()
        };
        let label = format!("height {height} window {window:?}");
        let vectors = draw(&solid, options);
        assert!(
            vectors
                .iter()
                .all(|v| matches!(v.role, VectorRole::Section | VectorRole::Hatch))
        );
        // The view normal (1, 1, 1)/√3 meets the bore axis +Z at cos 1/√3, so
        // every void arc lies on the ellipse with semiaxes 2√3 and 2 about the
        // point where the axis crosses the plane.
        for vector in &vectors {
            if let VectorCurve::Ellipse {
                center,
                major_mm,
                minor_mm,
                ..
            } = vector.curve
            {
                assert!(near(center, [0., 0.]), "{label}");
                assert!((major_mm - 2. * 3_f64.sqrt()).abs() < TOLERANCE, "{label}");
                assert!((minor_mm - 2.).abs() < TOLERANCE, "{label}");
            }
        }
        let world =
            |[x, y]: [f64; 2]| [10. + s * x + t * y, 6. - s * x + t * y, height - 2. * t * y];
        let in_window = |[x, y]: [f64; 2]| {
            window.is_none_or(|([u0, u1], [v0, v1])| {
                x >= u0 - 1e-9 && x <= u1 + 1e-9 && y >= v0 - 1e-9 && y <= v1 + 1e-9
            })
        };
        let on_boundary = |p: [f64; 2]| {
            let w = world(p);
            (0..3).any(|i| w[i].abs() < 1e-9 || (w[i] - [20., 12., 8.][i]).abs() < 1e-9)
                || window.is_some_and(|([u0, u1], [v0, v1])| {
                    [p[0] - u0, p[0] - u1, p[1] - v0, p[1] - v1]
                        .iter()
                        .any(|d| d.abs() < 1e-9)
                })
                || ((w[0] - 10.).hypot(w[1] - 6.) - 2.).abs() < 1e-9
        };
        for vector in vectors.iter().filter(|v| v.role == VectorRole::Section) {
            if let VectorCurve::Line { a, b } = vector.curve {
                for q in [0., 0.5, 1.] {
                    let p = [a[0] + (b[0] - a[0]) * q, a[1] + (b[1] - a[1]) * q];
                    assert!(on_boundary(p) && in_window(p), "{label}: {p:?}");
                }
            }
        }
        let strokes: Vec<_> = vectors
            .iter()
            .filter(|v| v.role == VectorRole::Hatch)
            .collect();
        assert!(strokes.len() > 5, "{label}");
        for stroke in strokes {
            let VectorCurve::Line { a, b } = stroke.curve else {
                panic!("hatch strokes are lines");
            };
            for q in [0., 0.25, 0.5, 0.75, 1.] {
                let p = [a[0] + (b[0] - a[0]) * q, a[1] + (b[1] - a[1]) * q];
                let w = world(p);
                assert!(in_window(p), "{label}: {p:?}");
                assert!(
                    w.iter()
                        .zip([20., 12., 8.])
                        .all(|(v, m)| *v >= -1e-9 && *v <= m + 1e-9),
                    "{label}: {w:?}"
                );
                assert!(
                    (w[0] - 10.).hypot(w[1] - 6.) >= 2. - 1e-9,
                    "{label}: hatch enters the bore at {w:?}"
                );
            }
        }
    }
}

#[test]
fn an_oblique_section_through_the_bore_is_an_explicit_finding() {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let error = project_drawing(
        &[DrawingSolid {
            solid: &bored(),
            hatch: hatch(),
        }],
        DrawingProjection {
            frame: frame([10., 6., 4.], [s, -s, 0.], [t, t, -2. * t]),
            section_offset_mm: Some(0.),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
}

#[cfg(feature = "occt")]
#[test]
fn bored_box_views_and_sections_match_native_open_cascade() {
    let solid = bored();
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let axial = frame([10., 6., 4.], [1., 0., 0.], [0., 1., 0.]);
    let side = frame([0., 6., 0.], [1., 0., 0.], [0., 0., 1.]);
    let oblique = frame([0.; 3], [s, -s, 0.], [t, t, -2. * t]);
    let section = |frame| DrawingProjection {
        frame,
        section_offset_mm: Some(0.),
        ..Default::default()
    };
    for (label, options) in [
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
        (
            "side hidden",
            DrawingProjection {
                frame: side,
                hidden_lines: true,
                ..Default::default()
            },
        ),
        (
            "oblique hidden",
            DrawingProjection {
                frame: oblique,
                hidden_lines: true,
                ..Default::default()
            },
        ),
        ("axial section", section(axial)),
        (
            "axial section only",
            DrawingProjection {
                section_only: true,
                ..section(axial)
            },
        ),
        (
            "axial finite depth",
            DrawingProjection {
                section_depth_mm: Some(1.),
                ..section(axial)
            },
        ),
        ("side section", section(side)),
        (
            "side finite depth",
            DrawingProjection {
                section_depth_mm: Some(1.),
                ..section(side)
            },
        ),
        (
            "side window",
            DrawingProjection {
                section_u_mm: Some([1., 7.]),
                section_v_mm: Some([1., 7.]),
                ..section(side)
            },
        ),
        (
            "axial detail",
            DrawingProjection {
                frame: axial,
                detail: Some(([0., 0.], 3.)),
                ..Default::default()
            },
        ),
        (
            "oblique section only",
            DrawingProjection {
                section_only: true,
                ..section(PlaneFrame {
                    origin: [10., 6., 4.],
                    ..oblique
                })
            },
        ),
        (
            "oblique section only at a cap",
            DrawingProjection {
                section_only: true,
                ..section(PlaneFrame {
                    origin: [10., 6., 1.5],
                    ..oblique
                })
            },
        ),
    ] {
        assert_matches_native(&[&solid], options, &format!("bored box {label}"));
    }
}
