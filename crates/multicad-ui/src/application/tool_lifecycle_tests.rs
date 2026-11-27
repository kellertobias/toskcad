//! CAD-413: sketch tool lifecycle through the full application. Single-shot
//! tools return to Select after their creation; Enter or a plain right-click
//! ends a tool, keeping what was committed and dropping the pending gesture;
//! Escape still cancels; a plain right-click in Select restarts the last
//! tool with its settings, a right-drag only navigates.
use super::sketch_select_tests::SketchHarness;
use super::*;
use crate::sketch_editor::DrawTool;
use multicad_core::constraint_sketch::{ConstraintKind, CurveKind, Dimension};

const ZERO: egui::Vec2 = egui::Vec2::ZERO;

impl SketchHarness {
    fn tool_now(&self) -> Option<DrawTool> {
        self.h.app.current().sketch.tool
    }
    fn curves(&self) -> usize {
        self.source().sketch.curves.len()
    }
    fn right_click(&mut self, p: [f64; 2]) {
        let at = self.screen(p);
        self.button(at, egui::PointerButton::Secondary, egui::Modifiers::NONE);
    }
    fn enter(&mut self) {
        self.h.key(egui::Key::Enter);
        self.settle();
    }
    /// Click the drag-value field painted beside `caption`, replace its
    /// value and press Enter, as a user types a value.
    fn type_beside(&mut self, caption: &str, text: &str) {
        let out = self.h.render(vec![]);
        let label = out
            .shapes
            .iter()
            .find_map(|s| match &s.shape {
                egui::Shape::Text(t) if t.galley.text() == caption => {
                    Some(egui::Rect::from_min_size(t.pos, t.galley.size()))
                }
                _ => None,
            })
            .unwrap_or_else(|| panic!("Missing {caption}"));
        let field = out
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .filter(|(_, n)| {
                matches!(
                    n.role(),
                    egui::accesskit::Role::SpinButton | egui::accesskit::Role::TextInput
                )
            })
            .filter_map(|(_, n)| n.bounds())
            .find(|b| {
                b.x0 as f32 >= label.max.x - 1.
                    && (b.y0 as f32) <= label.center().y
                    && (b.y1 as f32) >= label.center().y
            })
            .unwrap_or_else(|| panic!("a field beside {caption}"));
        out.drop_without_applying_deltas();
        let at = egui::pos2(
            ((field.x0 + field.x1) * 0.5) as f32,
            ((field.y0 + field.y1) * 0.5) as f32,
        );
        self.h.click_at(at);
        assert!(self.h.ctx.text_edit_focused(), "{caption} field focused");
        let key = |key, pressed, modifiers| egui::Event::Key {
            key,
            physical_key: None,
            pressed,
            repeat: false,
            modifiers,
        };
        for events in [
            vec![key(egui::Key::A, true, egui::Modifiers::COMMAND)],
            vec![key(egui::Key::A, false, egui::Modifiers::COMMAND)],
            vec![
                egui::Event::Text(text.into()),
                key(egui::Key::Enter, true, egui::Modifiers::NONE),
            ],
            vec![key(egui::Key::Enter, false, egui::Modifiers::NONE)],
        ] {
            self.h.frame(events);
        }
        self.settle();
    }
    /// Radius dimensions of every sketch arc, in mm.
    fn fillet_radii(&self) -> Vec<f64> {
        self.source()
            .sketch
            .constraints
            .iter()
            .filter_map(|c| match &c.kind {
                ConstraintKind::Radius {
                    dimension: Dimension::Literal { value, .. },
                    ..
                } => Some(*value),
                _ => None,
            })
            .collect()
    }
}

/// AC1: rectangle, circle, line, arc and spline each end after their one
/// creation; Polyline and Fillet keep going.
#[test]
fn single_shot_tools_return_to_select_and_multi_step_tools_continue() {
    let mut s = SketchHarness::new();
    for (tool, clicks) in [
        (DrawTool::Rectangle, vec![[0., 0.], [40., 30.]]),
        (DrawTool::Circle, vec![[-50., 0.], [-40., 0.]]),
        (DrawTool::Line, vec![[0., 60.], [40., 60.]]),
        (DrawTool::Arc, vec![[-50., 50.], [-35., 50.], [-50., 65.]]),
        (DrawTool::RectangleCenter, vec![[70., -30.], [80., -20.]]),
    ] {
        let before = s.curves();
        s.tool(tool);
        for p in clicks {
            s.click(p, ZERO);
        }
        assert!(s.curves() > before, "{tool:?} created: {}", s.message());
        assert_eq!(s.tool_now(), None, "{tool:?} returned to Select");
        assert!(s.h.app.current().sketch.select_active());
        assert!(
            s.message()
                .contains(&format!("right-click repeats {}", tool.label())),
            "{tool:?}: {}",
            s.message()
        );
    }
    // Spline: finishing with Enter is its one creation.
    s.tool(DrawTool::Spline);
    for p in [[0., -60.], [20., -40.], [40., -60.], [60., -40.]] {
        s.click(p, ZERO);
    }
    let before = s.curves();
    s.enter();
    assert_eq!(s.curves(), before + 1, "{}", s.message());
    assert_eq!(s.tool_now(), None);
    // Polyline continues after each committed segment.
    s.tool(DrawTool::Polyline);
    for p in [[-70., -70.], [-40., -70.], [-40., -40.]] {
        s.click(p, ZERO);
    }
    assert_eq!(s.tool_now(), Some(DrawTool::Polyline));
    s.h.key(egui::Key::Escape);
    // Fillet continues after each corner.
    s.tool(DrawTool::Fillet);
    s.click([40., 30.], ZERO);
    s.click([0., 30.], ZERO);
    assert_eq!(s.fillet_radii().len(), 2, "{}", s.message());
    assert_eq!(s.tool_now(), Some(DrawTool::Fillet));
}

/// AC2: Enter and a plain right-click end a tool; committed geometry stays,
/// the pending gesture is dropped. Escape still cancels.
#[test]
fn enter_and_right_click_end_the_tool_keeping_committed_geometry() {
    let mut s = SketchHarness::new();
    // Enter: two committed polyline segments stay, the pending vertex goes.
    s.tool(DrawTool::Polyline);
    for p in [[0., 0.], [30., 0.], [30., 30.]] {
        s.click(p, ZERO);
    }
    let revision = s.h.revision();
    assert_eq!(s.curves(), 2);
    s.enter();
    assert_eq!(s.tool_now(), None);
    assert_eq!(s.curves(), 2);
    assert_eq!(s.h.revision(), revision, "ending commits nothing");
    assert!(
        s.message().contains("right-click repeats Polyline"),
        "{}",
        s.message()
    );
    // Enter with only a first rectangle corner: nothing is created.
    s.tool(DrawTool::Rectangle);
    s.click([60., 0.], ZERO);
    s.enter();
    assert_eq!(s.tool_now(), None);
    assert_eq!(s.curves(), 2);
    assert_eq!(s.h.revision(), revision);
    // Right-click ends Fillet after one corner; the fillet stays.
    s.draw(DrawTool::Rectangle, &["-60", "-50", "-20", "-20"]);
    s.tool(DrawTool::Fillet);
    s.click([-20., -20.], ZERO);
    assert_eq!(s.fillet_radii().len(), 1);
    let revision = s.h.revision();
    s.right_click([70., 20.]);
    assert_eq!(s.tool_now(), None, "{}", s.message());
    assert_eq!(s.fillet_radii().len(), 1);
    assert_eq!(s.h.revision(), revision);
    // Right-click drops a pending line start.
    s.tool(DrawTool::Line);
    s.click([0., 60.], ZERO);
    let curves = s.curves();
    s.right_click([70., 20.]);
    assert_eq!(s.tool_now(), None);
    s.click([40., 60.], ZERO);
    assert_eq!(s.curves(), curves, "no line from the dropped start");
    // Escape still cancels the gesture and the tool.
    s.tool(DrawTool::Circle);
    s.click([0., -40.], ZERO);
    let document = s.h.app.current().session.document().clone();
    s.h.key(egui::Key::Escape);
    assert_eq!(s.tool_now(), None);
    assert_eq!(s.h.app.current().session.document(), &document);
}

/// AC3: a plain right-click in Select restarts the last tool with the fillet
/// radius and construction mode it had; a right-drag only navigates.
#[test]
fn right_click_in_select_repeats_the_last_tool_with_its_settings_and_right_drag_does_not() {
    // Repeating Fillet restores its radius of 10 mm.
    let mut s = SketchHarness::new();
    s.draw(DrawTool::Rectangle, &["0", "0", "60", "40"]);
    s.tool(DrawTool::Fillet);
    s.type_beside("Fillet radius", "10");
    s.click([60., 40.], ZERO);
    assert_eq!(s.fillet_radii(), vec![10.], "{}", s.message());
    s.enter();
    assert_eq!(s.tool_now(), None);
    s.right_click([30., 20.]);
    assert_eq!(s.tool_now(), Some(DrawTool::Fillet), "{}", s.message());
    assert!(
        s.message().contains("Sketch fillet repeated"),
        "{}",
        s.message()
    );
    s.click([0., 40.], ZERO);
    assert_eq!(s.fillet_radii(), vec![10., 10.], "{}", s.message());
    s.enter();

    // Construction mode is part of the repeated settings.
    s.h.app
        .invoke(Command::Draw(DrawTool::Line), &[], Channel::Ui)
        .unwrap();
    s.h.app
        .invoke(Command::Construction, &[], Channel::Ui)
        .unwrap();
    assert!(s.h.app.current().sketch.construction());
    s.click([0., 80.], ZERO);
    s.click([40., 80.], ZERO);
    assert_eq!(s.tool_now(), None);
    s.right_click([30., 20.]);
    assert_eq!(s.tool_now(), Some(DrawTool::Line));
    assert!(s.h.app.current().sketch.construction());
    s.click([0., -40.], ZERO);
    s.click([40., 100.], ZERO);
    let source = s.source();
    let lines: Vec<_> =
        source
            .sketch
            .curves
            .iter()
            .filter(|c| matches!(c.kind, CurveKind::Line { .. }) && c.construction)
            .filter(|c| {
                !source.sketch.shapes.iter().any(|shape| {
                    shape.name.starts_with("Rectangle") && shape.curves.contains(&c.id)
                })
            })
            .collect();
    assert!(lines.len() >= 2, "both lines are construction");

    // A right-drag navigates and does not repeat.
    assert_eq!(s.tool_now(), None);
    let camera = s.h.app.current().sketch.camera;
    let start = s.screen([30., 20.]);
    s.h.frame(vec![egui::Event::PointerMoved(start)]);
    let button = |pos, pressed| egui::Event::PointerButton {
        pos,
        button: egui::PointerButton::Secondary,
        pressed,
        modifiers: egui::Modifiers::NONE,
    };
    s.h.frame(vec![button(start, true)]);
    for i in 1..=8 {
        let p = start + egui::vec2(6. * i as f32, 3. * i as f32);
        s.h.frame(vec![egui::Event::PointerMoved(p)]);
    }
    let end = start + egui::vec2(48., 24.);
    s.h.frame(vec![button(end, false)]);
    s.settle();
    assert_eq!(s.tool_now(), None, "a right-drag does not repeat the tool");
    assert_ne!(
        s.h.app.current().sketch.camera,
        camera,
        "the drag navigated"
    );
}
