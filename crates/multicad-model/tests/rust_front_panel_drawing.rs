use multicad_core::{
    document::Document,
    drawings::{self, Page, VectorCurve, VectorRole},
    project_file,
};
use multicad_kernel::{
    DrawingProjection, DrawingSolid, PlaneFrame, ProjectedVector, project_drawing,
};
use multicad_model::{Compiler, GeometryBackend, drawing_preview_sources};
use std::f64::consts::TAU;

fn fixture() -> (Document, Page) {
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

fn project_selection(
    document: &Document,
    page: &Page,
    backend: GeometryBackend,
    selection: impl Fn(&str) -> bool,
    expected_count: usize,
) -> Vec<Vec<ProjectedVector>> {
    let roots = drawing_preview_sources(document, page.id).unwrap();
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_sources_with_backend(document, 0, &roots, backend, None)
    .unwrap();
    page.views
        .iter()
        .map(|view| {
            let roots = view.source.entities();
            let input: Vec<_> = model
                .snapshot()
                .leaf_occurrences
                .iter()
                .filter(|(path, leaf)| {
                    roots.contains(&path.root)
                        && selection(&document.get(leaf.producer).unwrap().name)
                })
                .map(|(path, _)| DrawingSolid {
                    solid: model.leaf_solid(path).unwrap(),
                    hatch: None,
                })
                .collect();
            assert_eq!(input.len(), expected_count);
            project_drawing(
                &input,
                DrawingProjection {
                    frame: PlaneFrame {
                        origin: view.frame.origin_mm,
                        u: view.frame.u,
                        v: view.frame.v,
                    },
                    hidden_lines: view.hidden_lines,
                    tangent_edges: view.tangent_edges,
                    ..Default::default()
                },
            )
            .unwrap()
        })
        .collect()
}

fn project_members(
    document: &Document,
    page: &Page,
    backend: GeometryBackend,
    selection: (&str, usize),
) -> Vec<Vec<ProjectedVector>> {
    project_selection(
        document,
        page,
        backend,
        |name| name.ends_with(selection.0),
        selection.1,
    )
}

fn machined_members(
    document: &Document,
    page: &Page,
    backend: GeometryBackend,
) -> Vec<Vec<ProjectedVector>> {
    project_selection(
        document,
        page,
        backend,
        |name| name.ends_with("holes") || name.contains("Hex standoff"),
        4,
    )
}

#[test]
fn saved_machined_plates_and_standoffs_project_all_three_directions() {
    let (document, page) = fixture();
    let views = machined_members(&document, &page, GeometryBackend::Rust);
    assert_eq!(
        views.iter().map(Vec::len).collect::<Vec<_>>(),
        [148, 176, 154]
    );
    assert_eq!(
        views[0]
            .iter()
            .filter(|edge| edge.role == VectorRole::Visible
                && matches!(edge.curve, VectorCurve::Ellipse { major_mm: 1.5, .. }))
            .count(),
        4
    );
    assert_eq!(
        views[2]
            .iter()
            .filter(|edge| edge.role == VectorRole::Hidden
                && matches!(edge.curve, VectorCurve::Ellipse { major_mm: 4., .. }))
            .count(),
        2
    );
}

#[test]
fn saved_front_panel_buttons_project_in_all_three_authored_directions() {
    let (document, page) = fixture();
    let views = project_members(&document, &page, GeometryBackend::Rust, ("Button union", 6));
    let counts: Vec<_> = views
        .iter()
        .map(|view| {
            let visible = view
                .iter()
                .filter(|edge| edge.role == VectorRole::Visible)
                .count();
            [visible, view.len() - visible]
        })
        .collect();
    assert_eq!(counts, [[12, 36], [48, 96], [32, 112]]);
    assert!(
        views[0]
            .iter()
            .all(|edge| matches!(edge.curve, VectorCurve::Ellipse {
        major_mm, minor_mm, start_rad, end_rad, ..
    } if major_mm == minor_mm && (end_rad - start_rad - TAU).abs() < 1e-10))
    );
    assert!(
        views[1..]
            .iter()
            .flatten()
            .all(|edge| matches!(edge.curve, VectorCurve::Line { .. }))
    );
}

#[test]
fn saved_socket_cap_envelopes_publish_authored_views_without_rear_rim_overdraw() {
    let (document, page) = fixture();
    let views = project_members(
        &document,
        &page,
        GeometryBackend::Rust,
        ("ISO 4762 M3 × 8 feature", 2),
    );
    assert_eq!(views.iter().map(Vec::len).collect::<Vec<_>>(), [48, 16, 48]);
    let top = &views[1];
    assert!(
        top.iter()
            .all(|edge| matches!(edge.curve, VectorCurve::Ellipse { .. }))
    );
    assert_eq!(
        top.iter()
            .filter(|edge| edge.role == VectorRole::Visible)
            .count(),
        2
    );
}
