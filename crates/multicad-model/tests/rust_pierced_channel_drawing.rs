use multicad_core::{
    commands::{Actor, Channel, Session},
    document::EntityId,
    drawings::{self, VectorCurve, VectorRole},
    project_file,
};
use multicad_kernel::{
    DrawingProjection, DrawingSolid, GeometryErrorKind, PlaneFrame, ProjectedVector, Solid,
    project_drawing,
};
use multicad_model::{Compiler, GeometryBackend, drawing_preview_sources};
use multicad_precision::PortableMath as _;
use std::f64::consts::{FRAC_1_SQRT_2, PI, TAU};

const OBLIQUE_V: [f64; 3] = [0., FRAC_1_SQRT_2, FRAC_1_SQRT_2];
const TOL: f64 = 1e-6;

fn fixture() -> (multicad_core::document::Document, drawings::Page) {
    let document = project_file::read(
        &mut include_bytes!("../../../examples/front-panel.tcd").as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    let page = drawings::books(&document)
        .flat_map(|id| drawings::read(&document, id).unwrap().pages)
        .find(|page| page.number == "FP-101")
        .unwrap();
    (document, page)
}

fn with_channel<T>(
    backend: GeometryBackend,
    inspect: impl FnOnce(&Solid, [DrawingProjection; 3]) -> T,
) -> T {
    let (document, page) = fixture();
    let roots = drawing_preview_sources(&document, page.id).unwrap();
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_sources_with_backend(&document, 0, &roots, backend, None)
    .unwrap();
    let leaf = model
        .snapshot()
        .leaf_occurrences
        .iter()
        .find(|(_, leaf)| document.get(leaf.producer).unwrap().name == "Front panel sheet feature")
        .unwrap();
    let options = ["Front", "Top", "Left end"].map(|name| {
        let view = page.views.iter().find(|view| view.name == name).unwrap();
        DrawingProjection {
            frame: PlaneFrame {
                origin: view.frame.origin_mm,
                u: view.frame.u,
                v: view.frame.v,
            },
            hidden_lines: view.hidden_lines,
            tangent_edges: view.tangent_edges,
            ..Default::default()
        }
    });
    inspect(model.leaf_solid(leaf.0).unwrap(), options)
}

fn project(solid: &Solid, options: DrawingProjection) -> Vec<ProjectedVector> {
    project_drawing(&[DrawingSolid { solid, hatch: None }], options).unwrap()
}

fn oblique(options: DrawingProjection) -> DrawingProjection {
    DrawingProjection {
        frame: PlaneFrame {
            v: OBLIQUE_V,
            ..options.frame
        },
        ..options
    }
}

fn merge_intervals(mut spans: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    spans.sort_unstable_by(|a, b| a[0].total_cmp(&b[0]));
    let mut merged: Vec<[f64; 2]> = Vec::new();
    for span in spans {
        if let Some(last) = merged.last_mut()
            && span[0] <= last[1] + TOL
        {
            last[1] = last[1].max(span[1]);
        } else {
            merged.push(span);
        }
    }
    merged
}

fn window_edge<'a>(
    vectors: impl IntoIterator<Item = (VectorRole, &'a VectorCurve)>,
) -> Vec<[f64; 2]> {
    let mut spans = Vec::new();
    for (role, curve) in vectors {
        if role != VectorRole::Visible {
            continue;
        }
        let line;
        let points: &[[f64; 2]] = match curve {
            VectorCurve::Line { a, b } => {
                line = [*a, *b];
                &line
            }
            VectorCurve::BSpline {
                degree: 1,
                periodic: false,
                poles,
                ..
            } => poles,
            _ => continue,
        };
        if points
            .iter()
            .all(|p| (p[1] - 26. * FRAC_1_SQRT_2).abs() < TOL)
        {
            let [a, b] = points
                .iter()
                .fold([f64::INFINITY, f64::NEG_INFINITY], |[a, b], p| {
                    [a.min(p[0]), b.max(p[0])]
                });
            let span = [a.max(34.), b.min(106.)];
            if span[1] - span[0] > TOL {
                spans.push(span);
            }
        }
    }
    merge_intervals(spans)
}

fn rim_intervals<'a>(
    vectors: impl IntoIterator<Item = (VectorRole, &'a VectorCurve)>,
    role: VectorRole,
    expected_center: [f64; 2],
    radius: f64,
) -> Vec<[f64; 2]> {
    let mut spans = Vec::new();
    for (actual, curve) in vectors {
        let VectorCurve::Ellipse {
            center,
            u,
            v,
            major_mm,
            minor_mm,
            start_rad,
            end_rad,
        } = *curve
        else {
            continue;
        };
        if actual != role
            || (center[0] - expected_center[0]).portable_hypot(center[1] - expected_center[1]) > TOL
            || (major_mm - radius).abs() > TOL
            || (minor_mm - radius * FRAC_1_SQRT_2).abs() > TOL
        {
            continue;
        }
        assert!(u[1].abs() < TOL && v[0].abs() < TOL);
        // Convert either authored ellipse basis to +X/+Y and merge trims.
        let base = if u[0] >= 0. { 0. } else { PI };
        let phase = if u[0] * v[1] > 0. {
            start_rad
        } else {
            -end_rad
        };
        let start = (base + phase).rem_euclid(TAU);
        let end = start + end_rad - start_rad;
        if TAU - start > TOL {
            spans.push([start, end.min(TAU)]);
        }
        if end - TAU > TOL {
            spans.push([0., end - TAU]);
        }
    }
    merge_intervals(spans)
}

fn same_intervals(actual: &[[f64; 2]], expected: &[[f64; 2]]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?} != {expected:?}");
    for (a, b) in actual.iter().flatten().zip(expected.iter().flatten()) {
        assert!((a - b).abs() <= TOL, "{actual:?} != {expected:?}");
    }
}

#[test]
fn saved_pierced_channel_end_view_preserves_the_hidden_cuts() {
    with_channel(GeometryBackend::Rust, |solid, [_, _, options]| {
        let vectors = project(solid, options);
        assert_eq!(vectors.len(), 88);
        let visible: Vec<_> = vectors
            .iter()
            .filter(|vector| vector.role == VectorRole::Visible)
            .map(|vector| &vector.curve)
            .collect();
        assert_eq!(visible.len(), 12);
        assert_eq!(
            visible
                .iter()
                .filter(|curve| matches!(curve, VectorCurve::Ellipse { .. }))
                .count(),
            4,
        );
        assert!(vectors.iter().any(|vector| {
            vector.role == VectorRole::Hidden
                && matches!(vector.curve, VectorCurve::Line { a, b }
                    if (a[0]-b[0]).abs() > 1. && (a[1]-31.25).abs() < 1e-8 && (b[1]-31.25).abs() < 1e-8)
        }));
        let visible_only = project(
            solid,
            DrawingProjection {
                hidden_lines: false,
                ..options
            },
        );
        assert_eq!(
            visible,
            visible_only.iter().map(|v| &v.curve).collect::<Vec<_>>()
        );

        let tiny_budget = DrawingProjection {
            maximum_curves: 12,
            ..options
        };
        let input = DrawingSolid { solid, hatch: None };
        assert_eq!(
            project_drawing(&[input], tiny_budget).unwrap_err().kind,
            GeometryErrorKind::Limit,
        );
    });
}

#[test]
fn saved_pierced_channel_projects_both_transverse_views_with_hidden_lines() {
    with_channel(GeometryBackend::Rust, |solid, [front, top, _]| {
        for (options, [visible_count, hidden_count], circular_count) in
            [(front, [18, 38], 10), (top, [6, 62], 2)]
        {
            let vectors = project(solid, options);
            let (visible, hidden): (Vec<_>, Vec<_>) = vectors
                .into_iter()
                .partition(|vector| vector.role == VectorRole::Visible);
            assert_eq!([visible.len(), hidden.len()], [visible_count, hidden_count]);
            assert!(
                hidden
                    .iter()
                    .all(|vector| vector.role == VectorRole::Hidden)
            );
            assert_eq!(
                visible
                    .iter()
                    .filter(|vector| matches!(vector.curve, VectorCurve::Ellipse { .. }))
                    .count(),
                circular_count
            );
            assert_eq!(
                visible,
                project(
                    solid,
                    DrawingProjection {
                        hidden_lines: false,
                        ..options
                    }
                )
            );
        }
    });
}

#[test]
fn saved_pierced_channel_projects_the_oblique_window_and_round_openings() {
    with_channel(GeometryBackend::Rust, |solid, [front, _, _]| {
        let options = oblique(front);
        let vectors = project(solid, options);
        let curves = || vectors.iter().map(|v| (v.role, &v.curve));
        same_intervals(&window_edge(curves()), &[[34., 106.]]);
        same_intervals(
            &rim_intervals(
                curves(),
                VectorRole::Visible,
                [150., 53.5 * FRAC_1_SQRT_2],
                6.25,
            ),
            &[[0., TAU]],
        );
        assert!(
            !rim_intervals(
                curves(),
                VectorRole::Hidden,
                [150., 55. * FRAC_1_SQRT_2],
                6.25
            )
            .is_empty()
        );
        assert_eq!(
            vectors
                .iter()
                .filter(|v| v.role == VectorRole::Visible)
                .collect::<Vec<_>>(),
            project(
                solid,
                DrawingProjection {
                    hidden_lines: false,
                    ..options
                }
            )
            .iter()
            .collect::<Vec<_>>(),
        );
    });
}

fn assembly_views(backend: GeometryBackend) -> [Vec<drawings::DrawingVector>; 3] {
    let (document, page) = fixture();
    let roots = drawing_preview_sources(&document, page.id).unwrap();
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_sources_with_backend(&document, 0, &roots, backend, None)
    .unwrap();
    let result = model
        .drawing_page_preview(&document, 0, &Default::default(), page.id)
        .unwrap();
    ["Front", "Top", "Left end"].map(|name| {
        let view = page.views.iter().find(|view| view.name == name).unwrap();
        let evaluation = &result.views[&view.id];
        let drawings::ViewEvaluation::Ready { vectors, .. } = evaluation else {
            panic!("saved {name} view: {evaluation:?}");
        };
        vectors.clone()
    })
}

#[test]
fn saved_front_panel_assembly_draws_the_near_window_and_six_buttons() {
    let views = assembly_views(GeometryBackend::Rust);
    for vectors in &views {
        assert!(vectors.iter().any(|v| v.role == VectorRole::Visible));
        assert!(vectors.iter().any(|v| v.role == VectorRole::Hidden));
    }
    let front = &views[0];
    assert!(front.iter().any(|v| v.role == VectorRole::Visible
        && matches!(
            v.curve,
            VectorCurve::Line {
                a: [34., 27.5],
                b: [106., 27.5]
            }
        )));
    assert_eq!(
        front
            .iter()
            .filter(|v| v.role == VectorRole::Visible
                && matches!(
                    v.curve,
                    VectorCurve::Ellipse {
                        major_mm: 7.5,
                        minor_mm: 7.5,
                        ..
                    }
                ))
            .count(),
        6
    );
}

#[test]
fn assembled_bent_channel_only_reveals_foreign_edges_through_material_openings() {
    use multicad_kernel::{GeometryPolicy, rust_kernel::Rational};

    with_channel(GeometryBackend::Rust, |channel, [front, _, _]| {
        let frame = front.frame;
        let n = frame.normal();
        let behind = |u| PlaneFrame {
            origin: std::array::from_fn(|i| {
                frame.origin[i] + u * frame.u[i] + 40. * frame.v[i] - 300. * n[i]
            }),
            ..frame
        };
        let make_box = |id, u| {
            Solid::rust_box(
                id,
                [2; 3].map(|value| Rational::new(value, 1).unwrap()),
                behind(u),
                GeometryPolicy::default(),
            )
            .unwrap()
        };
        let [through_window, behind_web] = [make_box(1, 50.), make_box(2, 20.)];
        let scene = [channel, &through_window, &behind_web]
            .map(|solid| DrawingSolid { solid, hatch: None });
        let vectors = project_drawing(&scene, front).unwrap();
        for (source, visible) in [(1, true), (2, false)] {
            let vectors: Vec<_> = vectors
                .iter()
                .filter(|v| v.source_index == source)
                .collect();
            assert!(!vectors.is_empty());
            assert_eq!(
                vectors.iter().any(|v| v.role == VectorRole::Visible),
                visible
            );
        }
    });
}

fn oblique_page(hidden_lines: bool) -> (multicad_core::document::Document, EntityId, EntityId) {
    let (document, saved) = fixture();
    let mut front = saved
        .views
        .iter()
        .find(|view| view.name == "Front")
        .unwrap()
        .clone();
    front.id = EntityId::new();
    front.name = "Oblique front".into();
    front.frame.v = OBLIQUE_V;
    front.hidden_lines = hidden_lines;
    front.guide = None;
    let view = front.id;
    let mut page = drawings::Page::new("FP-101-OBLIQUE");
    page.settings = saved.settings;
    page.views.push(front);
    let page_id = page.id;
    let mut session = Session::new(document).unwrap();
    let (_, prepared) = drawings::prepare_create(
        &session,
        0,
        Actor::local(Channel::Ui),
        "FP-101 oblique",
        Some(page),
    )
    .unwrap();
    session.commit(prepared).unwrap();
    (session.document().clone(), page_id, view)
}

fn oblique_assembly(
    document: &multicad_core::document::Document,
    page: EntityId,
    view: EntityId,
    backend: GeometryBackend,
) -> Vec<drawings::DrawingVector> {
    let roots = drawing_preview_sources(document, page).unwrap();
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_sources_with_backend(document, 1, &roots, backend, None)
    .unwrap();
    let result = model
        .drawing_page_preview(document, 1, &Default::default(), page)
        .unwrap();
    let drawings::ViewEvaluation::Ready { vectors, .. } = &result.views[&view] else {
        panic!("FP-101 oblique {backend:?}: {:?}", result.views[&view]);
    };
    vectors.clone()
}

#[test]
fn bundled_front_panel_assembly_projects_an_oblique_drawing_with_six_button_caps() {
    for hidden_lines in [true, false] {
        let (document, page, view) = oblique_page(hidden_lines);
        let vectors = oblique_assembly(&document, page, view, GeometryBackend::Rust);
        let curves = || vectors.iter().map(|v| (v.role, &v.curve));
        same_intervals(&window_edge(curves()), &[[34., 106.]]);
        assert_eq!(
            vectors.iter().any(|v| v.role == VectorRole::Hidden),
            hidden_lines
        );
        for x in [150., 185., 220.] {
            for z in [25., 55.] {
                let center = [x, (z - 5.5) * FRAC_1_SQRT_2];
                same_intervals(
                    &rim_intervals(curves(), VectorRole::Visible, center, 7.5),
                    &[[0., TAU]],
                );
            }
        }
    }
}

#[test]
fn transverse_projection_requires_coincident_or_separate_openings_in_two_layers() {
    use multicad_core::{document::EntityId, references::Origin};
    use multicad_kernel::{CurvedEdge2d, CurvedProfile, GeometryPolicy, SweepOpening};

    let origin = |role: &str| Origin {
        source: EntityId::new(),
        role: role.into(),
    };
    let vertices = [
        [0., 0.],
        [20., 0.],
        [20., 2.],
        [4., 2.],
        [4., 8.],
        [20., 8.],
        [20., 10.],
        [0., 10.],
    ];
    let sides = vertices.map(|_| origin("stock side"));
    let profile = CurvedProfile {
        frame: PlaneFrame::default(),
        outer: vertices
            .into_iter()
            .zip(vertices.into_iter().cycle().skip(1))
            .map(|(a, b)| CurvedEdge2d::Line { a, b })
            .collect(),
        holes: vec![],
    };
    let policy = GeometryPolicy::default();
    let stock =
        Solid::rust_extrude_normal_profile_named(&profile, 25., EntityId::new(), &sides, policy)
            .unwrap();
    let make_opening = |wall: [usize; 2], y, travel, x| SweepOpening {
        walls: wall.map(|i| sides[i].clone()),
        frame: PlaneFrame {
            origin: [0., y, 0.],
            u: [1., 0., 0.],
            v: [0., 0., 1.],
        },
        through_mm: [0., travel, 0.],
        outline: vec![CurvedEdge2d::Circle {
            center: [x, 12.],
            radius_mm: 1.,
        }],
        origins: vec![origin("round opening")],
    };
    let pierced = |x| {
        stock
            .rust_pierce_normal_profile(
                &[
                    make_opening([0, 2], 0., 2., 10.),
                    make_opening([6, 4], 10., -2., x),
                ],
                policy,
            )
            .unwrap()
    };
    let options = DrawingProjection {
        frame: PlaneFrame {
            origin: [0.; 3],
            u: [0., 0., 1.],
            v: [1., 0., 0.],
        },
        hidden_lines: true,
        ..Default::default()
    };
    for x in [10., 14.] {
        let solid = pierced(x);
        let vectors = project(&solid, options);
        let rims: Vec<_> = vectors
            .iter()
            .filter(|v| matches!(v.curve, VectorCurve::Ellipse { .. }))
            .collect();
        assert_eq!(rims.len(), 4);
        assert_eq!(
            rims.iter()
                .filter(|v| v.role == VectorRole::Visible)
                .count(),
            1
        );
    }
    let offset = pierced(10.5);
    let error = project_drawing(
        &[DrawingSolid {
            solid: &offset,
            hatch: None,
        }],
        options,
    )
    .unwrap_err();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert!(error.message.contains("overlapping projected piercings"));
}
