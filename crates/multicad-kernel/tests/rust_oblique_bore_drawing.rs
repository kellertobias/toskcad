//! Sections of the Rust oblique-bore solid: box stock minus one finite
//! oblique cylinder (CAD-354 geometry).
//!
//! The analytic section of the stock is a polygon, and the bore removes an
//! ellipse trimmed by its cap planes. Section-only cuts take any transversal
//! plane. Retained cuts are admitted when the view runs along the bore axis
//! and all retained stock lies behind the cut face.
use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
    references::Origin,
};
use multicad_kernel::{
    BooleanOperation, DrawingHatch, DrawingProjection, DrawingSolid, GeometryPolicy, PlaneFrame,
    ProjectedVector, Solid, project_drawing,
    rust_kernel::{Circle, Rational, cylinder::AxialCylinder},
};

// The oblique bore refuses STEP exchange, so only the reference comparison
// with an independently built native solid applies.
#[cfg(feature = "occt")]
#[path = "native_drawing/reference.rs"]
mod native_drawing;

type V = [f64; 3];

fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalized(v: V) -> V {
    let length = dot(v, v).sqrt();
    v.map(|x| x / length)
}

/// A frame whose normal is the requested direction.
fn frame(origin: V, normal: V) -> PlaneFrame {
    let w = normalized(normal);
    let seed = if w[2].abs() < 0.9 {
        [0., 0., 1.]
    } else {
        [1., 0., 0.]
    };
    let u = normalized(cross(seed, w));
    PlaneFrame {
        origin,
        u,
        v: cross(w, u),
    }
}

fn stock(size: [i64; 3]) -> Solid {
    Solid::rust_box(
        EntityId::new().uuid_bits(),
        size.map(|v| Rational::new(v, 1).unwrap()),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// The exact Rust cylinder of `Solid::cylinder` in a Rust-only build.
fn cylinder(origin: V, travel: V, radius: i64, tenths: bool, height: i64) -> Solid {
    let id = EntityId::new();
    let r = if tenths {
        Rational::new(radius, 10).unwrap()
    } else {
        Rational::new(radius, 1).unwrap()
    };
    let zero = Rational::new(0, 1).unwrap();
    Solid::from_axial_cylinder_sweep(
        AxialCylinder::new(
            id.uuid_bits(),
            Circle::new([zero; 2], r).unwrap(),
            zero,
            Rational::new(height, 1).unwrap(),
        )
        .unwrap(),
        frame(origin, travel),
        [0.; 2],
        false,
        Origin {
            source: id,
            role: "primitive.circle".into(),
        },
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn cut(stock: &Solid, tool: &Solid) -> Solid {
    let result = stock
        .boolean(tool, BooleanOperation::Subtract, GeometryPolicy::default())
        .unwrap()
        .unwrap();
    assert!(
        !result.supports_step_export(),
        "the oblique-bore solid keeps its explicit exchange refusal"
    );
    result
}

fn hatch() -> Option<DrawingHatch> {
    Some(DrawingHatch {
        angle_rad: 0.4,
        spacing_mm: 0.2,
        line_width_mm: 0.18,
        pattern: multicad_core::materials::HatchPattern::Lines,
    })
}

fn draw(
    solid: &Solid,
    options: DrawingProjection,
) -> Result<Vec<ProjectedVector>, multicad_kernel::GeometryError> {
    project_drawing(
        &[DrawingSolid {
            solid,
            hatch: hatch(),
        }],
        options,
    )
}

/// Drawing coordinates of a world point in the frame of `options`.
fn screen(options: DrawingProjection, p: V) -> [f64; 2] {
    let u = normalized(options.frame.u);
    let n = normalized(cross(u, options.frame.v));
    let v = cross(n, u);
    let d = [0, 1, 2].map(|i| p[i] - options.frame.origin[i]);
    [dot(d, u), dot(d, v)]
}

fn world(options: DrawingProjection, depth: f64, [x, y]: [f64; 2]) -> V {
    let u = normalized(options.frame.u);
    let n = normalized(cross(u, options.frame.v));
    let v = cross(n, u);
    std::array::from_fn(|i| options.frame.origin[i] + x * u[i] + y * v[i] + depth * n[i])
}

/// Radial distance from the bore axis and position along it.
fn bore_coordinates(p: V, base: V, axis: V) -> (f64, f64) {
    let d = [0, 1, 2].map(|i| p[i] - base[i]);
    let along = dot(d, axis);
    let radial = [0, 1, 2].map(|i| d[i] - along * axis[i]);
    (dot(radial, radial).sqrt(), along)
}

fn hatch_points(vectors: &[ProjectedVector]) -> Vec<[f64; 2]> {
    vectors
        .iter()
        .filter(|v| v.role == VectorRole::Hatch)
        .flat_map(|v| {
            let VectorCurve::Line { a, b } = v.curve else {
                panic!("hatch strokes are lines");
            };
            [0., 0.25, 0.5, 0.75, 1.].map(|q| [a[0] + (b[0] - a[0]) * q, a[1] + (b[1] - a[1]) * q])
        })
        .collect()
}

#[test]
fn a_cut_across_an_oblique_through_bore_keeps_its_exact_ellipse() {
    let center = [5., 4., 3.];
    let axis = normalized([1., 0.5, 1.]);
    let base = [0, 1, 2].map(|i| center[i] - 10. * axis[i]);
    let solid = cut(&stock([10, 8, 6]), &cylinder(base, axis, 1, false, 20));
    let options = DrawingProjection {
        frame: PlaneFrame {
            origin: [0., 0., 3.],
            ..Default::default()
        },
        section_offset_mm: Some(0.),
        section_only: true,
        ..Default::default()
    };
    let vectors = draw(&solid, options).unwrap();
    let ellipses: Vec<_> = vectors
        .iter()
        .filter_map(|v| match v.curve {
            VectorCurve::Ellipse {
                center,
                u,
                major_mm,
                minor_mm,
                start_rad,
                end_rad,
                ..
            } => Some((center, u, major_mm, minor_mm, end_rad - start_rad)),
            _ => None,
        })
        .collect();
    let [(c, u, major, minor, sweep)] = ellipses.as_slice() else {
        panic!("{ellipses:?}");
    };
    // The axis meets z = 3 at the center; the major axis follows its
    // projection (1, 0.5) and stretches the radius by 1 / |axis_z|.
    assert!(
        (c[0] - 5.).abs() < 1e-9 && (c[1] - 4.).abs() < 1e-9,
        "{c:?}"
    );
    assert!((major - 1. / axis[2]).abs() < 1e-9 && (minor - 1.).abs() < 1e-9);
    let projected = normalized([1., 0.5, 0.]);
    assert!((u[0] * projected[0] + u[1] * projected[1]).abs() > 1. - 1e-9);
    assert!((sweep - std::f64::consts::TAU).abs() < 1e-12);
    let lines: Vec<_> = vectors
        .iter()
        .filter(|v| v.role == VectorRole::Section)
        .filter(|v| matches!(v.curve, VectorCurve::Line { .. }))
        .collect();
    assert_eq!(lines.len(), 4, "the rectangle of the stock");
    let points = hatch_points(&vectors);
    assert!(points.len() > 50);
    for p in points {
        assert!((-1e-9..=10. + 1e-9).contains(&p[0]) && (-1e-9..=8. + 1e-9).contains(&p[1]));
        let (radial, _) = bore_coordinates([p[0], p[1], 3.], base, axis);
        assert!(radial >= 1. - 1e-9, "hatch enters the bore at {p:?}");
    }
}

#[test]
fn a_cut_through_a_blind_oblique_floor_draws_its_cap_chord() {
    // An oblique blind pocket entering from the top face, 4 mm along its axis.
    let axis = normalized([1., 0., -2.]);
    let entry = [3., 4., 6.];
    let base = [0, 1, 2].map(|i| entry[i] - 2. * axis[i]);
    let solid = cut(&stock([10, 8, 6]), &cylinder(base, axis, 10, true, 6));
    let floor = [0, 1, 2].map(|i| base[i] + 6. * axis[i]);
    // A vertical plane y = 4 through the axis cuts the floor disk on a chord.
    let options = DrawingProjection {
        frame: PlaneFrame {
            origin: [0., 4., 0.],
            u: [1., 0., 0.],
            v: [0., 0., 1.],
        },
        section_offset_mm: Some(0.),
        section_only: true,
        ..Default::default()
    };
    let vectors = draw(&solid, options).unwrap();
    let chord: Vec<_> = vectors
        .iter()
        .filter(|v| v.role == VectorRole::Section)
        .filter_map(|v| match v.curve {
            VectorCurve::Line { a, b } => Some([a, b]),
            _ => None,
        })
        .filter(|[a, b]| {
            [a, b].iter().all(|p| {
                let w = world(options, 0., **p);
                let d = [0, 1, 2].map(|i| w[i] - floor[i]);
                dot(d, axis).abs() < 1e-9
            })
        })
        .collect();
    let [[a, b]] = chord.as_slice() else {
        panic!("one floor chord: {chord:?}");
    };
    // The plane contains the axis, so the chord spans the floor diameter.
    assert!(((a[0] - b[0]).hypot(a[1] - b[1]) - 2.).abs() < 1e-9);
    for p in hatch_points(&vectors) {
        let (radial, along) = bore_coordinates(world(options, 0., p), base, axis);
        assert!(
            radial >= 1. - 1e-9 || along >= 6. - 1e-9 || along <= 1e-9,
            "hatch enters the pocket at {p:?}"
        );
    }
}

/// The CAD-354 audit geometry: a radius-0.1 bore along the cube diagonal.
fn audit() -> (Solid, V, V) {
    let axis = normalized([1., 1., 1.]);
    let n = 3_f64.sqrt().recip();
    let base = [1. / 3. - n; 3];
    let id = EntityId::new();
    let tool = Solid::from_axial_cylinder_sweep(
        AxialCylinder::new(
            id.uuid_bits(),
            Circle::new(
                [Rational::new(0, 1).unwrap(); 2],
                Rational::new(1, 10).unwrap(),
            )
            .unwrap(),
            Rational::new(0, 1).unwrap(),
            Rational::new(20, 1).unwrap(),
        )
        .unwrap(),
        PlaneFrame {
            origin: base,
            u: [0.5_f64.sqrt(), -(0.5_f64.sqrt()), 0.],
            v: [
                (1_f64 / 6.).sqrt(),
                (1_f64 / 6.).sqrt(),
                -2. * (1_f64 / 6.).sqrt(),
            ],
        },
        [0.; 2],
        false,
        Origin {
            source: id,
            role: "primitive.circle".into(),
        },
        GeometryPolicy::default(),
    )
    .unwrap();
    (cut(&stock([10, 10, 10]), &tool), base, axis)
}

/// Lines as ordered end pairs and conics as center, semi-axes and total sweep,
/// rounded to a nanometre and sorted.
fn canonical(vectors: &[&ProjectedVector]) -> (Vec<[i64; 4]>, Vec<[i64; 5]>) {
    let round = |value: f64| (value * 1e9).round() as i64;
    let mut lines = Vec::new();
    let mut conics: Vec<[i64; 5]> = Vec::new();
    for vector in vectors {
        match vector.curve {
            VectorCurve::Line { a, b } => {
                let [a, b] = [a.map(round), b.map(round)];
                let [first, second] = if a <= b { [a, b] } else { [b, a] };
                lines.push([first[0], first[1], second[0], second[1]]);
            }
            VectorCurve::Ellipse {
                center,
                major_mm,
                minor_mm,
                start_rad,
                end_rad,
                ..
            } => {
                let key = [
                    round(center[0]),
                    round(center[1]),
                    round(major_mm),
                    round(minor_mm),
                ];
                // Sweeps add over the pieces of one conic; compare them to a
                // microradian.
                let sweep = ((end_rad - start_rad) * 1e12).round() as i64;
                match conics.iter_mut().find(|known| known[..4] == key) {
                    Some(known) => known[4] += sweep,
                    None => conics.push([key[0], key[1], key[2], key[3], sweep]),
                }
            }
            _ => panic!("unexpected curve"),
        }
    }
    for conic in &mut conics {
        conic[4] = (conic[4] as f64 / 1e6).round() as i64;
    }
    lines.sort_unstable();
    conics.sort_unstable();
    (lines, conics)
}

#[test]
fn corner_cuts_along_the_bore_axis_publish_retained_and_finite_sections() {
    let (solid, base, axis) = audit();
    for offset in [1., 2.] {
        // The plane x + y + z = 3 * offset retains only the corner at origin.
        let plane = PlaneFrame {
            origin: [offset; 3],
            u: [0.5_f64.sqrt(), -(0.5_f64.sqrt()), 0.],
            v: [
                (1_f64 / 6.).sqrt(),
                (1_f64 / 6.).sqrt(),
                -2. * (1_f64 / 6.).sqrt(),
            ],
        };
        let base_options = DrawingProjection {
            frame: plane,
            section_offset_mm: Some(0.),
            ..Default::default()
        };
        let only = draw(
            &solid,
            DrawingProjection {
                section_only: true,
                ..base_options
            },
        )
        .unwrap();
        for (label, options) in [
            ("retained", base_options),
            (
                "finite",
                DrawingProjection {
                    section_depth_mm: Some(0.5),
                    ..base_options
                },
            ),
        ] {
            let vectors = draw(&solid, options).unwrap();
            let sections: Vec<_> = vectors
                .iter()
                .filter(|v| v.role == VectorRole::Section)
                .collect();
            let visible: Vec<_> = vectors
                .iter()
                .filter(|v| v.role == VectorRole::Visible)
                .collect();
            // The same curves, whichever provider publishes them: equal line
            // sets up to direction and order, and equal conics by center,
            // semi-axes and total sweep.
            let only: Vec<_> = only
                .iter()
                .filter(|v| v.role == VectorRole::Section)
                .collect();
            assert_eq!(canonical(&sections), canonical(&only), "{label} {offset}");
            // Retained geometry behind the cut coincides with its boundary.
            assert_eq!(
                canonical(&visible),
                canonical(&sections),
                "{label} {offset}"
            );
        }
        let rims: Vec<_> = only
            .iter()
            .filter_map(|v| match v.curve {
                VectorCurve::Ellipse {
                    major_mm, minor_mm, ..
                } => Some((major_mm, minor_mm)),
                _ => None,
            })
            .collect();
        assert_eq!(rims.len(), 1);
        assert!((rims[0].0 - 0.1).abs() < 1e-12 && (rims[0].1 - 0.1).abs() < 1e-12);
        for p in hatch_points(&only) {
            let w = world(base_options, 0., p);
            assert!(w.iter().all(|v| *v >= -1e-9), "{w:?}");
            assert!((w.iter().sum::<f64>() - 3. * offset).abs() < 1e-9);
            let (radial, _) = bore_coordinates(w, base, axis);
            assert!(radial >= 0.1 - 1e-9, "hatch enters the bore at {w:?}");
        }
        let _ = screen(base_options, base);
    }
}

#[test]
fn views_and_general_retained_sections_draw() {
    let (solid, ..) = audit();
    let view = PlaneFrame {
        origin: [0.; 3],
        u: [0.5_f64.sqrt(), -(0.5_f64.sqrt()), 0.],
        v: [
            (1_f64 / 6.).sqrt(),
            (1_f64 / 6.).sqrt(),
            -2. * (1_f64 / 6.).sqrt(),
        ],
    };
    for options in [
        DrawingProjection {
            frame: view,
            ..Default::default()
        },
        DrawingProjection {
            frame: view,
            hidden_lines: true,
            ..Default::default()
        },
        // Through the cube center, retained corners lie outside the cut face.
        DrawingProjection {
            frame: PlaneFrame {
                origin: [5.; 3],
                ..view
            },
            section_offset_mm: Some(0.),
            ..Default::default()
        },
        // A retained cut not along the axis.
        DrawingProjection {
            frame: PlaneFrame {
                origin: [0., 0., 0.5],
                ..Default::default()
            },
            section_offset_mm: Some(0.),
            ..Default::default()
        },
    ] {
        // These draw through the scene engine (CAD-346); `rust_scene_bore.rs`
        // certifies their geometry.
        let vectors = draw(&solid, options).unwrap();
        assert!(
            vectors.iter().any(|v| v.role == VectorRole::Visible),
            "{options:?}"
        );
    }
}

#[cfg(feature = "occt")]
#[test]
fn oblique_bore_sections_match_the_native_open_cascade_cut() {
    let center = [5., 4., 3.];
    let axis = normalized([1., 0.5, 1.]);
    let base = [0, 1, 2].map(|i| center[i] - 10. * axis[i]);
    let policy = GeometryPolicy::default();
    let rust = cut(&stock([10, 8, 6]), &cylinder(base, axis, 1, false, 20));
    let native = Solid::box_solid([10., 8., 6.], policy)
        .unwrap()
        .boolean(
            &Solid::cylinder(frame(base, axis), 1., 20., policy).unwrap(),
            BooleanOperation::Subtract,
            policy,
        )
        .unwrap()
        .unwrap();
    for (label, frame) in [
        (
            "horizontal",
            PlaneFrame {
                origin: [0., 0., 3.],
                ..Default::default()
            },
        ),
        (
            "vertical",
            PlaneFrame {
                origin: [0., 4., 0.],
                u: [1., 0., 0.],
                v: [0., 0., 1.],
            },
        ),
        ("oblique", frame_from(center, normalized([1., -1., 2.]))),
    ] {
        let options = DrawingProjection {
            frame,
            section_offset_mm: Some(0.),
            section_only: true,
            ..Default::default()
        };
        let vectors = project_drawing(
            &[DrawingSolid {
                solid: &rust,
                hatch: None,
            }],
            options,
        )
        .unwrap();
        native_drawing::assert_matches_reference(
            &vectors,
            &[&native],
            options,
            &[],
            true,
            &format!("oblique bore {label} section only"),
        );
    }
}

/// Plain views and general retained sections drawn by the scene engine
/// (CAD-346), against native hidden-line removal of the native cut.
#[cfg(feature = "occt")]
#[test]
fn oblique_bore_views_and_retained_sections_match_native_open_cascade() {
    let center = [5., 4., 3.];
    let axis = normalized([1., 0.5, 1.]);
    let base = [0, 1, 2].map(|i| center[i] - 10. * axis[i]);
    let policy = GeometryPolicy::default();
    let rust = cut(&stock([10, 8, 6]), &cylinder(base, axis, 1, false, 20));
    let native = Solid::box_solid([10., 8., 6.], policy)
        .unwrap()
        .boolean(
            &Solid::cylinder(frame(base, axis), 1., 20., policy).unwrap(),
            BooleanOperation::Subtract,
            policy,
        )
        .unwrap()
        .unwrap();
    for (label, frame) in [
        (
            "top",
            PlaneFrame {
                origin: [0., 0., 3.],
                ..Default::default()
            },
        ),
        (
            "front",
            PlaneFrame {
                origin: [0., 4., 0.],
                u: [1., 0., 0.],
                v: [0., 0., 1.],
            },
        ),
        ("oblique", frame_from(center, normalized([1., -1., 2.]))),
        ("along the bore", frame_from(center, axis)),
    ] {
        let view = DrawingProjection {
            frame,
            hidden_lines: true,
            ..Default::default()
        };
        let section = DrawingProjection {
            section_offset_mm: Some(0.),
            ..view
        };
        for (case, options) in [
            ("view", view),
            (
                "detail",
                DrawingProjection {
                    detail: Some(([0.3, 0.2], 2.5)),
                    ..view
                },
            ),
            ("section", section),
            (
                "finite depth",
                DrawingProjection {
                    section_depth_mm: Some(1.5),
                    ..section
                },
            ),
            (
                "window",
                DrawingProjection {
                    section_u_mm: Some([-1.5, 2.25]),
                    section_v_mm: Some([-2., 1.75]),
                    ..section
                },
            ),
        ] {
            let vectors = project_drawing(
                &[DrawingSolid {
                    solid: &rust,
                    hatch: None,
                }],
                options,
            )
            .unwrap();
            native_drawing::assert_matches_reference(
                &vectors,
                &[&native],
                options,
                &[],
                true,
                &format!("oblique bore {label} {case}"),
            );
        }
    }
}

#[cfg(feature = "occt")]
fn frame_from(origin: V, normal: V) -> PlaneFrame {
    frame(origin, normal)
}
