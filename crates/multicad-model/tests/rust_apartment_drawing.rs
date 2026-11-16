use multicad_core::{
    document::Document,
    drawings::{self, DrawingVector, Page, VectorCurve, VectorRole, ViewEvaluation},
    project_file,
};
use multicad_model::{Compiler, GeometryBackend, drawing_preview_sources};
use std::f64::consts::{PI, TAU};

fn fixture(number: &str) -> (Document, Page) {
    let bytes = include_bytes!("../../../examples/apartment.tcd");
    let document = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    let page = drawings::books(&document)
        .flat_map(|id| drawings::read(&document, id).unwrap().pages)
        .find(|page| page.number == number)
        .unwrap();
    (document, page)
}

fn preview(document: &Document, page: &Page, backend: GeometryBackend) -> Vec<DrawingVector> {
    let [vectors] = preview_views(document, page, backend)
        .try_into()
        .expect("expected one saved drawing view");
    vectors
}

fn preview_views(
    document: &Document,
    page: &Page,
    backend: GeometryBackend,
) -> Vec<Vec<DrawingVector>> {
    let roots = drawing_preview_sources(document, page.id).unwrap();
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_sources_with_backend(document, 0, &roots, backend, None)
    .unwrap();
    let result = model
        .drawing_page_preview(document, 0, &Default::default(), page.id)
        .unwrap();
    page.views
        .iter()
        .map(|view| {
            let ViewEvaluation::Ready { vectors, .. } = &result.views[&view.id] else {
                panic!(
                    "{} / {}: {:?}",
                    page.number, view.name, result.views[&view.id]
                );
            };
            vectors.clone()
        })
        .collect()
}

#[test]
fn apartment_floor_plan_publishes_the_real_worktop_apertures() {
    let (document, page) = fixture("AP-101");
    let vectors = preview(&document, &page, GeometryBackend::Rust);
    assert_eq!(
        vectors
            .iter()
            .filter(|v| v.role == VectorRole::Visible)
            .count(),
        235
    );
    assert_eq!(
        vectors
            .iter()
            .filter(|v| v.role == VectorRole::Section)
            .count(),
        60
    );
    assert!(vectors.iter().any(|v| v.role == VectorRole::Hatch));

    let workbench = page.views[0]
        .source
        .entities()
        .into_iter()
        .find(|&id| document.get(id).unwrap().name == "Workbench")
        .unwrap();
    let worktop: Vec<_> = vectors.iter().filter(|v| v.source == workbench).collect();
    assert_eq!(worktop.len(), 104);
    assert!(worktop.iter().all(|v| v.role == VectorRole::Visible));
    let mut centers = worktop
        .iter()
        .filter_map(|v| {
            let VectorCurve::Ellipse {
                center,
                major_mm,
                minor_mm,
                start_rad,
                end_rad,
                ..
            } = v.curve
            else {
                assert!(matches!(v.curve, VectorCurve::Line { .. }));
                return None;
            };
            assert_eq!([major_mm, minor_mm], [10.; 2]);
            assert!((end_rad - start_rad - TAU).abs() < 1e-12);
            Some(center.map(|v| v.round() as i64))
        })
        .collect::<Vec<_>>();
    centers.sort_unstable();
    let expected: Vec<_> = (0..10)
        .flat_map(|x| (0..10).map(move |y| [540 + 80 * x, 2940 + 80 * y]))
        .collect();
    assert_eq!(centers, expected);
}

#[test]
fn apartment_isometric_projects_workbench_apertures_with_room_occlusion() {
    let (document, page) = fixture("AP-102");
    let source = |name| {
        page.views[0]
            .source
            .entities()
            .into_iter()
            .find(|&id| document.get(id).unwrap().name == name)
            .unwrap()
    };
    let workbench = source("Workbench");
    let cabinet = source("Cabinet");
    let vectors = preview(&document, &page, GeometryBackend::Rust);
    assert_eq!(vectors.len(), 212);
    assert!(vectors.iter().all(|v| v.role == VectorRole::Visible));
    let cabinet_vectors: Vec<_> = vectors.iter().filter(|v| v.source == cabinet).collect();
    assert_eq!(cabinet_vectors.len(), 52);
    assert_eq!(
        cabinet_vectors
            .iter()
            .filter(|v| matches!(v.curve, VectorCurve::Ellipse { major_mm: 2.5, .. }))
            .count(),
        13
    );

    let worktop: Vec<_> = vectors.iter().filter(|v| v.source == workbench).collect();
    assert_eq!(worktop.len(), 20);
    assert!(worktop.iter().all(|v| v.role == VectorRole::Visible));
    let mut centers: Vec<_> = worktop
        .iter()
        .filter_map(|v| {
            let VectorCurve::Ellipse {
                center,
                major_mm,
                minor_mm,
                start_rad,
                end_rad,
                ..
            } = v.curve
            else {
                assert!(matches!(v.curve, VectorCurve::Line { .. }));
                return None;
            };
            assert_eq!(major_mm, 10.);
            assert!((minor_mm - 10. / 3.0_f64.sqrt()).abs() < 1e-10);
            assert!((end_rad - start_rad - PI).abs() < 1e-10);
            Some(center[0])
        })
        .collect();
    centers.sort_by(f64::total_cmp);
    let expected = (0..10)
        .rev()
        .map(|index| -1187.9393923933994 - 80. * f64::from(index) / 2.0_f64.sqrt());
    assert_eq!(centers.len(), 10);
    assert!(
        centers
            .into_iter()
            .zip(expected)
            .all(|(actual, expected)| (actual - expected).abs() < 1e-9)
    );
}

#[test]
fn apartment_workshop_references_publish_both_saved_views() {
    let (document, page) = fixture("AP-103");
    let views = preview_views(&document, &page, GeometryBackend::Rust);
    assert_eq!(views.iter().map(Vec::len).collect::<Vec<_>>(), [119, 44]);
    assert!(
        views
            .iter()
            .flatten()
            .all(|vector| vector.role == VectorRole::Visible)
    );
    assert_eq!(
        views[0]
            .iter()
            .filter(|vector| matches!(
                vector.curve,
                VectorCurve::Ellipse { major_mm: 10., minor_mm: 10., start_rad, end_rad, .. }
                    if (end_rad - start_rad - TAU).abs() < 1e-10
            ))
            .count(),
        100
    );
}
