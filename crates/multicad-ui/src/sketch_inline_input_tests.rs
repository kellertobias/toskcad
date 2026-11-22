//! Canvas tests for the inline creation fields (CAD-384). Pointer and keyboard
//! input run through the real canvas; commits use the shared draw transaction.
use super::*;
use multicad_core::constraint_sketch::{ConstraintKind, Dimension, source::Source};
use multicad_core::document::Entity;
use multicad_core::project_file::{self, ProjectFile};

struct Harness {
    ctx: egui::Context,
    editor: SketchEditor,
    session: Session,
    time: f64,
}
impl Harness {
    fn new(tool: DrawTool, setup: impl FnOnce(&mut Session)) -> Self {
        let mut session = Session::new(Document::new("Inline input").unwrap()).unwrap();
        setup(&mut session);
        let mut editor = SketchEditor::default();
        editor.create(&mut session, Datum::Xy, Channel::Ui).unwrap();
        editor.snap = false;
        editor.geometry_snap = false;
        editor.set_tool(tool).unwrap();
        let mut h = Self {
            ctx: egui::Context::default(),
            editor,
            session,
            time: 0.,
        };
        h.frame(vec![]);
        h.frame(vec![]);
        h
    }
    fn settle(&mut self) {
        super::tests::settle(&mut self.editor, &self.session);
    }
    fn frame(&mut self, events: Vec<egui::Event>) {
        self.settle();
        self.time += 0.05;
        let raw = egui::RawInput {
            focused: true,
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1000., 700.))),
            events,
            time: Some(self.time),
            ..Default::default()
        };
        let (editor, session) = (&mut self.editor, &mut self.session);
        self.ctx
            .run_ui(raw, |ui| editor.canvas(ui, session, true, None, false))
            .drop_without_applying_deltas();
        self.settle();
    }
    fn screen(&self, p: [f64; 2]) -> Pos2 {
        let s = self.editor.camera.project(
            Datum::Xy.frame().point(p).map(|v| v as f32),
            self.editor.viewport_rect.unwrap(),
        );
        Pos2::new(s[0], s[1])
    }
    fn hover(&mut self, p: [f64; 2]) {
        let at = self.screen(p);
        self.frame(vec![egui::Event::PointerMoved(at)]);
    }
    fn click(&mut self, p: [f64; 2]) {
        let at = self.screen(p);
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(at),
                egui::Event::PointerButton {
                    pos: at,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ]);
        }
    }
    fn key(&mut self, key: egui::Key, modifiers: egui::Modifiers) {
        for pressed in [true, false] {
            self.frame(vec![egui::Event::Key {
                key,
                physical_key: None,
                pressed,
                repeat: false,
                modifiers,
            }]);
        }
    }
    fn type_text(&mut self, text: &str) {
        self.frame(vec![egui::Event::Text(text.into())]);
    }
    /// Tab from the canvas (the application routes it here) into a field.
    fn enter_fields(&mut self) {
        self.editor.focus_inline_field(true);
        self.frame(vec![]);
        assert!(self.editor.inline_input_active(), "a field has focus");
    }
    fn source(&self) -> Source {
        let id = self.editor.selected.unwrap();
        inline_input::authored(self.session.document(), id).unwrap()
    }
    fn dimensions(&self) -> Vec<ConstraintKind> {
        self.source()
            .sketch
            .constraints
            .into_iter()
            .map(|c| c.kind)
            .filter(|k| k.dimension().is_some())
            .collect()
    }
}
fn points(source: &Source) -> Vec<[f64; 2]> {
    source.sketch.points.iter().map(|p| p.position_mm).collect()
}
fn literal(kind: &ConstraintKind) -> f64 {
    match kind.dimension() {
        Some(Dimension::Literal { value, .. }) => *value,
        other => panic!("expected a literal dimension, got {other:?}"),
    }
}

#[test]
fn rectangle_width_and_height_fields_fix_the_size_and_commit_one_dimensioned_draw() {
    let mut h = Harness::new(DrawTool::Rectangle, |_| {});
    // The first corner is typed in sketch coordinates and placed with Enter.
    h.hover([1., 1.]);
    h.enter_fields();
    h.type_text("5");
    h.key(egui::Key::Tab, egui::Modifiers::NONE);
    h.type_text("5");
    h.key(egui::Key::Enter, egui::Modifiers::NONE);
    h.frame(vec![]);
    assert_eq!(h.editor.first, Some([5., 5.]), "{}", h.editor.message);
    assert_eq!(
        inline_input::fields(DrawTool::Rectangle, true),
        &[inline_input::Field::Width, inline_input::Field::Height]
    );
    h.enter_fields();
    h.type_text("25");
    h.key(egui::Key::Tab, egui::Modifiers::NONE);
    h.type_text("10");
    assert_eq!(h.editor.inline.texts, vec!["25".to_string(), "10".into()]);
    // The pointer still chooses the side; the entered sizes stay fixed.
    h.hover([40., 30.]);
    assert_eq!(h.editor.inline.last_point, Some([30., 15.]));
    h.hover([-20., -1.]);
    assert_eq!(h.editor.inline.last_point, Some([-20., -5.]));
    h.hover([40., 30.]);
    let before = h.session.document().clone();
    let revision = h.session.revision();
    // Enter in a field commits the shown preview exactly once.
    h.editor.focus_inline_field(true);
    h.frame(vec![]);
    h.key(egui::Key::Enter, egui::Modifiers::NONE);
    h.frame(vec![]);
    assert_eq!(h.session.revision(), revision + 1, "{}", h.editor.message);
    assert_eq!(h.editor.first, None);
    let committed = points(&h.source());
    for corner in [[5., 5.], [30., 5.], [30., 15.], [5., 15.]] {
        assert!(committed.contains(&corner), "{corner:?} in {committed:?}");
    }
    let dimensions = h.dimensions();
    assert_eq!(dimensions.len(), 2, "{dimensions:?}");
    let mut sizes: Vec<f64> = dimensions.iter().map(literal).collect();
    sizes.sort_by(f64::total_cmp);
    assert_eq!(sizes, vec![10., 25.]);
    assert!(
        dimensions
            .iter()
            .all(|k| matches!(k, ConstraintKind::Length { .. }))
    );
    // The fields are empty again for the next rectangle.
    assert!(h.editor.inline.texts.iter().all(String::is_empty));
    // One Undo removes geometry and dimensions together; Redo restores both.
    let drawn = h.session.document().clone();
    h.session
        .undo(h.session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(h.session.document(), &before);
    h.session
        .redo(h.session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(h.session.document(), &drawn);
    // The dimensions survive save and reopen.
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(drawn.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let opened = project_file::read(
        &mut std::io::Cursor::new(bytes),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    assert_eq!(opened.project.document, drawn);
}

#[test]
fn circle_center_position_and_parameter_diameter_are_entered_beside_the_preview() {
    let mut parameter = None;
    let mut h = Harness::new(DrawTool::Circle, |session| {
        let d = Entity::new(
            "d",
            EntityData::Parameter {
                expression: "6".into(),
                unit: "mm".into(),
                exposed: false,
                inputs: vec![],
            },
        );
        parameter = Some(d.id);
        session
            .apply(Request {
                expected_revision: session.revision(),
                actor: Actor::local(Channel::Ui),
                label: "Parameter".into(),
                edits: vec![Mutation::Insert(d)],
            })
            .unwrap();
    });
    // Before the first point the fields are sketch X/Y; Enter places it.
    h.hover([1., 1.]);
    assert_eq!(
        inline_input::fields(DrawTool::Circle, false),
        &[inline_input::Field::X, inline_input::Field::Y]
    );
    h.enter_fields();
    h.type_text("3");
    h.key(egui::Key::Tab, egui::Modifiers::NONE);
    h.type_text("-4");
    h.hover([8., 9.]);
    assert_eq!(h.editor.inline.last_point, Some([3., -4.]));
    h.editor.focus_inline_field(false);
    h.frame(vec![]);
    let revision = h.session.revision();
    h.key(egui::Key::Enter, egui::Modifiers::NONE);
    h.frame(vec![]);
    assert_eq!(h.editor.first, Some([3., -4.]), "{}", h.editor.message);
    assert_eq!(
        h.session.revision(),
        revision,
        "placing a point is not an edit"
    );
    // The diameter accepts a parameter expression and keeps it bound.
    h.enter_fields();
    h.type_text("=$d * 2");
    h.hover([3., 40.]);
    let shown = h.editor.inline.last_point.unwrap();
    assert!(
        // The radius is exact; the direction comes from the pointer.
        ((shown[0] - 3.).portable_hypot(shown[1] + 4.) - 6.).abs() < 1e-9
            && (shown[0] - 3.).abs() < 1e-5,
        "{shown:?} {:?} {:?}",
        h.editor.inline.texts,
        h.editor.inline_invalid()
    );
    h.editor.focus_inline_field(true);
    h.frame(vec![]);
    h.key(egui::Key::Enter, egui::Modifiers::NONE);
    h.frame(vec![]);
    assert_eq!(h.session.revision(), revision + 1, "{}", h.editor.message);
    let source = h.source();
    let circle = source
        .sketch
        .curves
        .iter()
        .find_map(|c| match c.kind {
            multicad_core::constraint_sketch::CurveKind::Circle { center, radius_mm } => {
                Some((center, radius_mm))
            }
            _ => None,
        })
        .unwrap();
    assert!((circle.1 - 6.).abs() < 1e-9, "radius {}", circle.1);
    let center = source
        .sketch
        .points
        .iter()
        .find(|p| p.id == circle.0)
        .unwrap();
    assert_eq!(center.position_mm, [3., -4.]);
    let dimensions = h.dimensions();
    assert_eq!(dimensions.len(), 1);
    let ConstraintKind::Diameter { dimension, .. } = &dimensions[0] else {
        panic!("{dimensions:?}")
    };
    let expression = dimension.expression().expect("bound expression");
    assert!(
        expression.text.contains('d'),
        "the parameter stays referenced: {expression:?}"
    );
    let _ = parameter;
}

#[test]
fn line_length_field_and_invalid_or_cancelled_entries_change_nothing() {
    let mut h = Harness::new(DrawTool::Line, |_| {});
    h.click([0., 0.]);
    h.enter_fields();
    // An invalid expression is reported and blocks the commit.
    h.type_text("2 +");
    h.hover([10., 0.]);
    let before = h.session.document().clone();
    let revision = h.session.revision();
    assert!(h.editor.inline_invalid().is_some());
    h.click([10., 0.]);
    assert_eq!(h.session.document(), &before);
    assert!(
        h.editor.message.contains("nothing was placed"),
        "{}",
        h.editor.message
    );
    // A zero or negative size is refused the same way.
    h.editor.inline.texts[0] = "0".into();
    h.hover([10., 0.]);
    assert!(h.editor.inline_invalid().unwrap().contains("positive"));
    // A valid length fixes the distance along the pointer direction.
    h.editor.inline.texts[0] = "12.5".into();
    h.hover([0., -30.]);
    assert_eq!(
        h.editor.inline.last_point,
        Some([0., -12.5]),
        "{:?} {:?} {:?}",
        h.editor.inline.texts,
        h.editor.inline_invalid(),
        h.editor.inline.last_pointer
    );
    // Escape from a field drops the pending line; the tool stays active.
    h.editor.cancel_inline_creation();
    h.frame(vec![]);
    assert_eq!(h.editor.first, None);
    assert_eq!(h.editor.tool, Some(DrawTool::Line));
    assert_eq!(h.session.document(), &before);
    assert_eq!(h.session.revision(), revision);
    // A new line with a typed length commits it as a Length dimension.
    h.click([1., 1.]);
    h.enter_fields();
    h.type_text("12.5");
    h.click([1., 30.]);
    assert_eq!(h.session.revision(), revision + 1, "{}", h.editor.message);
    let dimensions = h.dimensions();
    assert_eq!(dimensions.len(), 1);
    assert!(matches!(dimensions[0], ConstraintKind::Length { .. }));
    assert_eq!(literal(&dimensions[0]), 12.5);
    // The pointer placed the start and direction; the length is exact.
    let source = h.source();
    let length = source
        .sketch
        .curves
        .iter()
        .rev()
        .find_map(|c| match c.kind {
            multicad_core::constraint_sketch::CurveKind::Line { start, end } => {
                let at = |id| {
                    source
                        .sketch
                        .points
                        .iter()
                        .find(|p| p.id == id)
                        .unwrap()
                        .position_mm
                };
                let (a, b) = (at(start), at(end));
                Some((a[0] - b[0]).portable_hypot(a[1] - b[1]))
            }
            _ => None,
        })
        .unwrap();
    assert!((length - 12.5).abs() < 1e-9, "{length}");
}

#[test]
fn typing_a_number_enters_the_first_field_and_fields_never_block_clicks() {
    let mut h = Harness::new(DrawTool::Rectangle, |_| {});
    h.hover([2., 2.]);
    h.click([0., 0.]);
    assert!(h.editor.pending_first().is_some());
    // Without touching a field, typing goes to the width field.
    h.type_text("3");
    h.type_text("0");
    h.frame(vec![]);
    assert_eq!(h.editor.inline.texts[0], "30");
    assert!(h.editor.inline_input_active());
    // The pointer still places the rest, even across where the fields float.
    let before = h.session.revision();
    for p in [[15., 10.], [15., 20.], [5., 25.]] {
        h.hover(p);
    }
    h.click([5., 25.]);
    assert_eq!(h.session.revision(), before + 1, "{}", h.editor.message);
    let dimensions = h.dimensions();
    assert_eq!(dimensions.len(), 1);
    assert_eq!(literal(&dimensions[0]), 30.);
}
