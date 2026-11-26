//! Sketch element list and shape inspector (CAD-398), driven through the real
//! editor host: left list, canvas and Inspector fields with injected time.
use super::shape_list::{ElementKey, ShapeField, field_id};
use super::snap::DWELL_SECONDS;
use super::*;
use multicad_core::constraint_sketch::{
    ConstraintKind, CurveKind, Dimension, Shape,
    shape_edit::{self, ShapeSize},
};
use multicad_core::project_file::{self, ProjectFile};

struct Host {
    ctx: egui::Context,
    editor: SketchEditor,
    session: Session,
    feed: crate::viewport::SceneFeed,
    time: f64,
}
fn key(key: egui::Key, modifiers: egui::Modifiers) -> egui::Event {
    egui::Event::Key {
        key,
        physical_key: None,
        pressed: true,
        repeat: false,
        modifiers,
    }
}
impl Host {
    fn new() -> Self {
        let mut session = Session::new(Document::new("Shape list").unwrap()).unwrap();
        let mut editor = SketchEditor::default();
        editor.create(&mut session, Datum::Xy, Channel::Ui).unwrap();
        let mut host = Self {
            ctx: egui::Context::default(),
            editor,
            session,
            feed: Default::default(),
            time: 0.,
        };
        host.frame(vec![]).drop_without_applying_deltas();
        host.frame(vec![]).drop_without_applying_deltas();
        host
    }
    fn frame(&mut self, events: Vec<egui::Event>) -> egui::FullOutput {
        self.time += 0.02;
        self.frame_at(events, self.time)
    }
    fn frame_at(&mut self, events: Vec<egui::Event>, time: f64) -> egui::FullOutput {
        self.time = self.time.max(time);
        let events_empty = events.is_empty();
        let raw = egui::RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1360., 900.))),
            events,
            time: Some(self.time),
            focused: true,
            ..Default::default()
        };
        let (editor, session, feed) = (&mut self.editor, &mut self.session, &mut self.feed);
        let output = self.ctx.run_ui(raw, |ui| {
            editor.show_project(ui, session, true, feed, &Default::default())
        });
        if !events_empty {
            self.settle();
        }
        output
    }
    /// Render until the sketch worker has published the current solution.
    fn settle(&mut self) {
        let until = std::time::Instant::now() + std::time::Duration::from_secs(5);
        loop {
            self.frame_at(vec![], self.time)
                .drop_without_applying_deltas();
            if !self.editor.sketch_evaluation.pending() {
                break;
            }
            assert!(std::time::Instant::now() < until, "Sketch worker timed out");
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
        self.frame_at(vec![], self.time)
            .drop_without_applying_deltas();
    }
    fn screen(&self, p: [f64; 2]) -> Pos2 {
        let s = self.editor.camera.project(
            Datum::Xy.frame().point(p).map(|v| v as f32),
            self.editor.viewport_rect.unwrap(),
        );
        Pos2::new(s[0], s[1])
    }
    fn hover(&mut self, p: [f64; 2], offset: Vec2, time: f64) {
        let at = self.screen(p) + offset;
        self.frame_at(vec![egui::Event::PointerMoved(at)], time)
            .drop_without_applying_deltas();
    }
    fn press(&mut self, at: Pos2) {
        self.frame(vec![egui::Event::PointerMoved(at)])
            .drop_without_applying_deltas();
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(at),
                egui::Event::PointerButton {
                    pos: at,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ])
            .drop_without_applying_deltas();
        }
    }
    fn click(&mut self, p: [f64; 2], offset: Vec2) {
        let at = self.screen(p) + offset;
        self.press(at);
    }
    /// Click a painted label (list entries).
    fn click_label(&mut self, label: &str) {
        let output = self.frame(vec![]);
        let at = output
            .shapes
            .iter()
            .find_map(|shape| match &shape.shape {
                egui::Shape::Text(text) if text.galley.text() == label => {
                    Some(text.pos + text.galley.size() * 0.5)
                }
                _ => None,
            })
            .unwrap_or_else(|| panic!("Missing {label}"));
        output.drop_without_applying_deltas();
        self.press(at);
    }
    fn labels(&mut self) -> Vec<String> {
        let output = self.frame(vec![]);
        let labels = output
            .shapes
            .iter()
            .filter_map(|shape| match &shape.shape {
                egui::Shape::Text(text) => Some(text.galley.text().to_string()),
                _ => None,
            })
            .collect();
        output.drop_without_applying_deltas();
        labels
    }
    /// Replace an Inspector field's text and press Enter.
    fn enter(&mut self, field: ShapeField, text: &str) {
        self.frame(vec![]).drop_without_applying_deltas();
        let rect = self
            .ctx
            .read_response(field_id(field))
            .unwrap_or_else(|| panic!("Missing field {field:?}"))
            .rect;
        self.press(rect.center());
        for event in [
            key(egui::Key::A, egui::Modifiers::COMMAND),
            egui::Event::Text(text.into()),
            key(egui::Key::Enter, egui::Modifiers::NONE),
        ] {
            self.frame(vec![event]).drop_without_applying_deltas();
        }
        self.frame(vec![]).drop_without_applying_deltas();
    }
    fn sketch(&self) -> EntityId {
        self.editor.selected.unwrap()
    }
    fn source(&self) -> Source {
        shape_list::canonical(self.session.document(), self.sketch()).unwrap()
    }
    fn solved(&self) -> multicad_core::constraint_sketch::Sketch {
        self.editor
            .sketch_evaluation
            .current(self.sketch())
            .unwrap()
            .evaluated
            .solved()
            .unwrap()
            .clone()
    }
}
fn at(sketch: &multicad_core::constraint_sketch::Sketch, id: EntityId) -> [f64; 2] {
    sketch
        .points
        .iter()
        .find(|p| p.id == id)
        .unwrap()
        .position_mm
}
fn line(
    sketch: &multicad_core::constraint_sketch::Sketch,
    curve: EntityId,
) -> (EntityId, EntityId) {
    match sketch.curves.iter().find(|c| c.id == curve).unwrap().kind {
        CurveKind::Line { start, end } => (start, end),
        _ => panic!("not a line"),
    }
}
fn near(a: [f64; 2], b: [f64; 2]) -> bool {
    (a[0] - b[0]).abs() < 1e-6 && (a[1] - b[1]).abs() < 1e-6
}
fn length_of(sketch: &multicad_core::constraint_sketch::Sketch, curve: EntityId) -> f64 {
    let (a, b) = line(sketch, curve);
    let (a, b) = (at(sketch, a), at(sketch, b));
    (b[0] - a[0]).portable_hypot(b[1] - a[1])
}
fn reopened(document: &Document) -> Document {
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(document.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    project_file::read(
        &mut std::io::Cursor::new(bytes),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document
}
/// Steps 1–3 of the product-owner scenario through the canvas: rectangle A
/// 100 × 60, then B from the point defined by A's top-edge midpoint and A's
/// center (acquired by dwelling) to A's bottom-right corner.
fn two_rectangles() -> (Host, Shape, Shape) {
    let mut h = Host::new();
    h.editor.set_tool(DrawTool::Rectangle).unwrap();
    h.frame(vec![]).drop_without_applying_deltas();
    h.click([0., 0.], Vec2::ZERO);
    h.click([100., 60.], Vec2::ZERO);
    assert_eq!(h.source().sketch.shapes.len(), 1, "{}", h.editor.message);
    let t = h.time + 1.;
    // Acquire the top-edge midpoint and the center as references.
    h.hover([50., 60.], Vec2::new(1., 1.), t);
    h.hover([50., 60.], Vec2::new(1., 1.), t + DWELL_SECONDS + 0.01);
    h.hover([50., 30.], Vec2::new(1., -1.), t + 1.);
    h.hover(
        [50., 30.],
        Vec2::new(1., -1.),
        t + 1. + DWELL_SECONDS + 0.01,
    );
    assert_eq!(h.editor.snap_assist.references().len(), 2);
    h.click([50., 30.], Vec2::new(1., -1.));
    h.click([100., 0.], Vec2::new(-2., 2.));
    let source = h.source();
    assert_eq!(source.sketch.shapes.len(), 2, "{}", h.editor.message);
    let (a, b) = (
        source.sketch.shapes[0].clone(),
        source.sketch.shapes[1].clone(),
    );
    (h, a, b)
}

#[test]
fn the_left_list_names_drawn_shapes_and_edits_rectangle_one_without_moving_rectangle_two() {
    let (mut h, a, b) = two_rectangles();
    let solved = h.solved();
    let a_corner = line(&solved, a.curves[1]).0; // A's bottom-right corner
    assert!(near(at(&solved, a_corner), [100., 0.]));
    let (b_start, _) = line(&solved, b.curves[0]);
    assert!(
        near(at(&solved, b_start), [50., 30.]),
        "B starts where the references meet"
    );
    assert!(
        shape_edit::points(&solved, &b).contains(&a_corner),
        "B ends on A's corner by identity"
    );
    // The left panel lists two shapes, not eight lines.
    let labels = h.labels();
    assert!(labels.iter().any(|l| l == "Rectangle 1"));
    assert!(labels.iter().any(|l| l == "Rectangle 2"));
    assert!(!labels.iter().any(|l| l.starts_with("Line")), "{labels:?}");
    // List → canvas: A's four curves highlight, none of B's.
    h.click_label("Rectangle 1");
    assert_eq!(h.editor.shape, Some(ElementKey::Shape(a.id)));
    let highlight = h.editor.element_highlight(h.session.document());
    assert_eq!(highlight, a.curves.iter().copied().collect());
    // The Inspector shows A's current sizes and first corner.
    let labels = h.labels();
    for label in ["Width", "Height", "X", "Y"] {
        assert!(labels.iter().any(|l| l == label), "{label}: {labels:?}");
    }
    let constraints_before = h.source().sketch.constraints;
    let revision = h.session.revision();
    h.enter(ShapeField::Size(0), "120");
    assert_eq!(h.session.revision(), revision + 1, "{}", h.editor.message);
    h.enter(ShapeField::Size(1), "50");
    assert_eq!(h.session.revision(), revision + 2, "{}", h.editor.message);
    assert!(
        h.editor.message.contains("Rectangle 1 height"),
        "{}",
        h.editor.message
    );
    // Only A gained dimensions, on its own curves.
    let source = h.source();
    let added: Vec<_> = source
        .sketch
        .constraints
        .iter()
        .filter(|c| !constraints_before.contains(c))
        .collect();
    assert_eq!(added.len(), 2);
    for (constraint, (curve, value)) in added.iter().zip([(a.curves[0], 120.), (a.curves[1], 50.)])
    {
        assert!(
            matches!(&constraint.kind, ConstraintKind::Length { line, dimension: Dimension::Literal { value: v, .. } } if *line == curve && *v == value),
            "{constraint:?}"
        );
    }
    assert_eq!(
        source.sketch.shapes,
        vec![a.clone(), b.clone()],
        "shape records unchanged"
    );
    // Solved: A is 120 × 50, B still ends on A's moved corner and is still a
    // rectangle starting where it was placed.
    let solved = h.solved();
    assert!((length_of(&solved, a.curves[0]) - 120.).abs() < 1e-6);
    assert!((length_of(&solved, a.curves[1]) - 50.).abs() < 1e-6);
    assert!(shape_edit::points(&solved, &b).contains(&a_corner));
    let corners: Vec<_> = b
        .curves
        .iter()
        .map(|c| at(&solved, line(&solved, *c).0))
        .collect();
    for i in 0..4 {
        let (p, q) = (corners[i], corners[(i + 1) % 4]);
        assert!(
            (p[0] - q[0]).abs() < 1e-6 || (p[1] - q[1]).abs() < 1e-6,
            "B edge {i} stays horizontal or vertical: {corners:?}"
        );
    }
    assert!(
        near(at(&solved, b_start), [50., 30.]),
        "{:?}",
        at(&solved, b_start)
    );
    // The Inspector now shows the dimension values.
    let labels = h.labels();
    assert!(
        labels.iter().any(|l| l == "120") && labels.iter().any(|l| l == "50"),
        "{labels:?}"
    );
    // Each edit is one Undo step and the result survives save/reopen.
    assert_eq!(reopened(h.session.document()), *h.session.document());
    h.session
        .undo(h.session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    let undone = h.source();
    assert_eq!(
        undone.sketch.constraints.len(),
        constraints_before.len() + 1
    );
}

#[test]
fn clicking_a_curve_in_the_canvas_selects_its_shape_in_the_list() {
    let (mut h, a, b) = two_rectangles();
    h.editor.tool = None;
    h.editor.clear_drawing_progress();
    // B's top edge, away from A.
    h.click([75., 30.], Vec2::new(0., 1.));
    assert_eq!(
        h.editor.shape,
        Some(ElementKey::Shape(b.id)),
        "{}",
        h.editor.message
    );
    let highlight = h.editor.element_highlight(h.session.document());
    assert_eq!(highlight, b.curves.iter().copied().collect());
    // A's left edge selects A; empty space clears the shape selection.
    h.click([0., 40.], Vec2::new(1., 0.));
    assert_eq!(h.editor.shape, Some(ElementKey::Shape(a.id)));
    h.click([-80., -80.], Vec2::ZERO);
    assert_eq!(h.editor.shape, None);
}

#[test]
fn a_size_that_over_constrains_the_sketch_is_reported_and_changes_nothing() {
    use multicad_core::sketch_draw::{DrawInput, DrawRequest, PointInput};
    let mut h = Host::new();
    h.editor
        .add(
            &mut h.session,
            DrawTool::Rectangle,
            [0., 0.],
            [100., 60.],
            Channel::Ui,
        )
        .unwrap();
    h.frame(vec![]).drop_without_applying_deltas();
    let a = h.source().sketch.shapes[0].clone();
    let sketch = h.source().sketch;
    let bottom_left = line(&sketch, a.curves[0]).0;
    let top_right = line(&sketch, a.curves[2]).0;
    // B spans A's diagonal corners, so B's width is A's width.
    h.editor
        .draw(
            &mut h.session,
            &DrawRequest {
                input: DrawInput::Rectangle {
                    a: PointInput::Existing(bottom_left),
                    b: PointInput::Existing(top_right),
                },
                construction: false,
            },
            Channel::Ui,
        )
        .unwrap();
    h.frame(vec![]).drop_without_applying_deltas();
    h.click_label("Rectangle 1");
    h.enter(ShapeField::Size(0), "120");
    let revision = h.session.revision();
    let document = h.session.document().clone();
    h.click_label("Rectangle 2");
    h.enter(ShapeField::Size(0), "100");
    assert_eq!(h.session.revision(), revision, "nothing committed");
    assert_eq!(h.session.document(), &document);
    assert!(
        h.editor.message.contains("Over-constrained"),
        "{}",
        h.editor.message
    );
    assert!(
        h.editor.message.contains("nothing changed"),
        "{}",
        h.editor.message
    );
    let labels = h.labels();
    assert!(
        labels.iter().any(|l| l.starts_with("Over-constrained")),
        "the refusal stays visible: {labels:?}"
    );
    // The field reverts to the solved value instead of keeping the refused one.
    assert!(labels.iter().any(|l| l == "120"), "{labels:?}");
    assert!(!labels.iter().any(|l| l == "100"), "{labels:?}");
}

#[test]
fn an_under_defined_shape_says_so_and_position_moves_and_fixes_it() {
    let (mut h, a, b) = two_rectangles();
    h.click_label("Rectangle 1");
    let labels = h.labels();
    assert!(
        labels.iter().any(|l| l.starts_with("Under-defined")),
        "{labels:?}"
    );
    // Moving the undimensioned A would let attached B squeeze it: refused.
    let revision = h.session.revision();
    h.enter(ShapeField::X, "10");
    assert_eq!(h.session.revision(), revision, "{}", h.editor.message);
    assert!(
        h.editor
            .message
            .contains("Moving Rectangle 1 would change its width"),
        "{}",
        h.editor.message
    );
    // With its sizes set, A moves rigidly and B stays attached.
    h.enter(ShapeField::Size(0), "100");
    h.enter(ShapeField::Size(1), "60");
    assert_eq!(h.session.revision(), revision + 2, "{}", h.editor.message);
    h.enter(ShapeField::X, "10");
    assert_eq!(h.session.revision(), revision + 3, "{}", h.editor.message);
    let solved = h.solved();
    let anchor = shape_edit::anchor(&h.source().sketch, &a).unwrap();
    assert!(
        near(at(&solved, anchor), [10., 0.]),
        "{:?}",
        at(&solved, anchor)
    );
    let corner = line(&solved, a.curves[1]).0;
    assert!(
        near(at(&solved, corner), [110., 0.]),
        "{:?}",
        at(&solved, corner)
    );
    assert!(shape_edit::points(&solved, &b).contains(&corner));
    assert!(
        h.source()
            .sketch
            .constraints
            .iter()
            .any(|c| matches!(c.kind, ConstraintKind::Fix { point } if point == anchor))
    );
    // Invalid input is refused without a transaction.
    h.enter(ShapeField::Size(0), "wide");
    assert_eq!(h.session.revision(), revision + 3);
    assert!(
        h.editor.message.contains("Enter a number"),
        "{}",
        h.editor.message
    );
}

#[test]
fn deleting_a_listed_shape_removes_all_its_curves_and_only_those() {
    let (mut h, a, b) = two_rectangles();
    h.click_label("Rectangle 2");
    h.editor.delete(&mut h.session, Channel::Ui).unwrap();
    let source = h.source();
    assert_eq!(source.sketch.shapes, vec![a.clone()]);
    assert!(
        b.curves
            .iter()
            .all(|c| source.sketch.curves.iter().all(|k| k.id != *c))
    );
    assert_eq!(source.sketch.curves.len(), 4);
}

#[test]
fn unrecorded_curves_and_texts_are_listed_one_by_one() {
    let mut source = Source {
        u: [1., 0., 0.],
        sketch: Default::default(),
    };
    let p = source.add_point([0., 0.]);
    let q = source.add_point([10., 0.]);
    source.add_line(p, q);
    source.add_circle(p, 4.);
    let document = Document::new("Legacy list").unwrap();
    let listed = shape_list::elements(&document, EntityId::new(), &source);
    let labels: Vec<_> = listed.iter().map(|e| e.label.as_str()).collect();
    assert_eq!(labels, ["Line 1 (single curve)", "Circle 1 (single curve)"]);
    assert!(listed.iter().all(|e| matches!(e.key, ElementKey::Curve(_))));
    let _ = ShapeSize::Width;
}
