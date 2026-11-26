//! Canonical drawn-shape identity (CAD-398): recorded by the draw request,
//! stable names, shared points, pruning, save/reopen, capability and sizes.
use multicad_core::{
    constraint_sketch::{
        ConstraintKind, CurveKind, Dimension, ShapeKind, Sketch,
        binding::DimensionInput,
        shape_edit::{self, ShapeSize},
        source::{self, Source},
    },
    document::{Document, EntityData},
    project_file::{self, ProjectFile},
    sketch::Datum,
    sketch_draw::{self, DrawInput, DrawRequest, PointInput},
};
fn empty() -> Source {
    Source {
        u: [1., 0., 0.],
        sketch: Sketch::default(),
    }
}
fn p(x: f64, y: f64) -> PointInput {
    PointInput::Position([x, y])
}
fn draw(source: &Source, input: DrawInput) -> Source {
    sketch_draw::draft(
        source,
        &DrawRequest {
            input,
            construction: false,
        },
    )
    .unwrap()
    .source
}
/// Rectangle A 100 × 60 at the origin, then B from A's center to A's
/// bottom-right corner, reusing that corner by identity.
fn two_rectangles() -> Source {
    let a = draw(
        &empty(),
        DrawInput::Rectangle {
            a: p(0., 0.),
            b: p(100., 60.),
        },
    );
    let corner = a
        .sketch
        .points
        .iter()
        .find(|q| q.position_mm == [100., 0.])
        .unwrap()
        .id;
    draw(
        &a,
        DrawInput::Rectangle {
            a: p(50., 30.),
            b: PointInput::Existing(corner),
        },
    )
}
fn literal(value: f64) -> DimensionInput {
    DimensionInput::Literal {
        value,
        unit: "mm".into(),
    }
}

#[test]
fn each_draw_request_records_one_named_shape_even_when_points_are_shared() {
    let s = two_rectangles();
    let shapes = &s.sketch.shapes;
    assert_eq!(shapes.len(), 2);
    assert_eq!(
        shapes
            .iter()
            .map(|s| (s.kind, s.name.as_str()))
            .collect::<Vec<_>>(),
        [
            (ShapeKind::Rectangle, "Rectangle 1"),
            (ShapeKind::Rectangle, "Rectangle 2")
        ]
    );
    assert!(shapes.iter().all(|shape| shape.curves.len() == 4));
    // Identity is by request, not geometry: the curves are disjoint even
    // though B reuses A's corner point.
    assert!(
        shapes[0]
            .curves
            .iter()
            .all(|c| !shapes[1].curves.contains(c))
    );
    let a_points = shape_edit::points(&s.sketch, &shapes[0]);
    let b_points = shape_edit::points(&s.sketch, &shapes[1]);
    assert_eq!(a_points.iter().filter(|p| b_points.contains(p)).count(), 1);
    // Other kinds get their own counters.
    let s = draw(
        &s,
        DrawInput::Circle {
            center: p(0., 0.),
            radius_mm: 5.,
        },
    );
    let s = draw(
        &s,
        DrawInput::Polyline {
            vertices: vec![p(0., 0.), p(5., 5.), p(10., 0.)],
            closed: false,
        },
    );
    let names: Vec<_> = s.sketch.shapes.iter().map(|s| s.name.as_str()).collect();
    assert_eq!(
        names,
        ["Rectangle 1", "Rectangle 2", "Circle 1", "Polyline 1"]
    );
    assert_eq!(s.sketch.shapes[3].curves.len(), 2);
}

#[test]
fn removed_curves_are_pruned_on_write_and_numbers_are_not_reused() {
    let mut s = two_rectangles();
    let first = s.sketch.shapes[0].clone();
    s.sketch.curves.retain(|c| !first.curves.contains(&c.id));
    s.sketch.constraints.retain(|c| {
        c.kind
            .references()
            .iter()
            .all(|r| !first.curves.contains(r))
    });
    let written = Source::from_definition(&s.definition().unwrap()).unwrap();
    assert_eq!(written.sketch.shapes.len(), 1);
    assert_eq!(written.sketch.shapes[0].name, "Rectangle 2");
    let next = draw(
        &written,
        DrawInput::Rectangle {
            a: p(200., 0.),
            b: p(210., 10.),
        },
    );
    assert_eq!(next.sketch.shapes[1].name, "Rectangle 3");
}

#[test]
fn a_curve_in_two_shapes_or_a_blank_name_is_invalid() {
    let mut s = two_rectangles();
    let stolen = s.sketch.shapes[0].curves[0];
    s.sketch.shapes[1].curves.push(stolen);
    assert!(
        s.validate()
            .unwrap_err()
            .to_string()
            .contains("at most one shape")
    );
    let mut s = two_rectangles();
    s.sketch.shapes[0].name = "  ".into();
    assert!(s.validate().is_err());
}

#[test]
fn shapes_survive_save_and_reopen_and_require_their_capability() {
    let (id, edits) = Source::create(Datum::Xy, "Sketch 1").unwrap();
    let document = Document::new("Shapes").unwrap().edited(&edits).unwrap();
    assert!(!source::requires_shapes(&document));
    let drawn = two_rectangles();
    let document = document
        .edited(&drawn.replace(&document, id).unwrap())
        .unwrap();
    assert!(source::requires_shapes(&document));
    let file = ProjectFile::new(document.clone()).unwrap();
    assert!(file.capabilities().contains(source::SHAPES_CAPABILITY));
    let mut bytes = vec![];
    project_file::write(&file, &mut bytes, Default::default()).unwrap();
    let opened = project_file::read(
        &mut std::io::Cursor::new(bytes),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    assert_eq!(opened, document);
    let EntityData::Sketch { definition, .. } = &opened.get(id).unwrap().data else {
        unreachable!()
    };
    assert_eq!(
        Source::from_definition(definition).unwrap().sketch.shapes,
        drawn.sketch.shapes
    );
}

#[test]
fn sources_saved_before_shapes_still_read_with_no_shapes() {
    let mut legacy = two_rectangles();
    legacy.sketch.shapes.clear();
    let definition = legacy.definition().unwrap();
    assert!(definition.data["sketch"].get("shapes").is_none());
    let read = Source::from_definition(&definition).unwrap();
    assert!(read.sketch.shapes.is_empty());
    assert_eq!(read.sketch.curves.len(), 8);
}

#[test]
fn sizes_edit_or_add_driving_dimensions_on_the_shapes_own_curves() {
    let s = two_rectangles();
    let document = Document::new("Sizes").unwrap();
    let a = s.sketch.shapes[0].clone();
    let width =
        shape_edit::set_size(&document, &s, a.id, ShapeSize::Width, &literal(120.)).unwrap();
    let added: Vec<_> = width
        .sketch
        .constraints
        .iter()
        .filter(|c| !s.sketch.constraints.contains(c))
        .collect();
    assert_eq!(added.len(), 1);
    assert!(matches!(
        &added[0].kind,
        ConstraintKind::Length { line, dimension: Dimension::Literal { value, .. } }
            if *line == a.curves[0] && *value == 120.
    ));
    // Editing again changes that same constraint rather than adding another.
    let again =
        shape_edit::set_size(&document, &width, a.id, ShapeSize::Width, &literal(130.)).unwrap();
    assert_eq!(
        again.sketch.constraints.len(),
        width.sketch.constraints.len()
    );
    let found = shape_edit::size_constraint(&again.sketch, &a, ShapeSize::Width).unwrap();
    assert!(matches!(
        found.kind.dimension(),
        Some(Dimension::Literal { value, .. }) if *value == 130.
    ));
    let height =
        shape_edit::set_size(&document, &again, a.id, ShapeSize::Height, &literal(50.)).unwrap();
    assert!(matches!(
        shape_edit::size_constraint(&height.sketch, &a, ShapeSize::Height).unwrap().kind,
        ConstraintKind::Length { line, .. } if line == a.curves[1]
    ));
    // Shape B is untouched; a rectangle has no diameter.
    let b = &s.sketch.shapes[1];
    assert!(shape_edit::size_constraint(&height.sketch, b, ShapeSize::Width).is_none());
    assert!(shape_edit::set_size(&document, &s, a.id, ShapeSize::Diameter, &literal(5.)).is_err());
}

#[test]
fn position_moves_every_shape_point_and_fixes_the_anchor() {
    let s = two_rectangles();
    let a = s.sketch.shapes[0].clone();
    let moved = shape_edit::set_position(&s, &s.sketch, a.id, [10., 20.]).unwrap();
    let anchor = shape_edit::anchor(&s.sketch, &a).unwrap();
    let at = |id| {
        moved
            .sketch
            .points
            .iter()
            .find(|q| q.id == id)
            .unwrap()
            .position_mm
    };
    assert_eq!(at(anchor), [10., 20.]);
    for id in shape_edit::points(&s.sketch, &a) {
        let before = s
            .sketch
            .points
            .iter()
            .find(|q| q.id == id)
            .unwrap()
            .position_mm;
        assert_eq!(at(id), [before[0] + 10., before[1] + 20.]);
    }
    assert!(
        moved
            .sketch
            .constraints
            .iter()
            .any(|c| matches!(c.kind, ConstraintKind::Fix { point } if point == anchor))
    );
    // B's own corners that A does not share stay where they were.
    let b = &s.sketch.shapes[1];
    let start = match s
        .sketch
        .curves
        .iter()
        .find(|c| c.id == b.curves[0])
        .unwrap()
        .kind
    {
        CurveKind::Line { start, .. } => start,
        _ => unreachable!(),
    };
    assert_eq!(at(start), [50., 30.]);
    assert!(shape_edit::set_position(&s, &s.sketch, a.id, [f64::NAN, 0.]).is_err());
}
