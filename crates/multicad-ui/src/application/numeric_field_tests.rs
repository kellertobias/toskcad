//! CAD-416: numeric fields in the real Model Inspector use the shared layout
//! (full-width input, label above, unit inside) and still commit edits.
use super::*;
use crate::modeling::{ModelingTool, native_fixture::NativeBackend};
use egui::accesskit::Role;
use multicad_core::{constraint_sketch::Dimension, model};
use std::{
    sync::Arc,
    time::{Duration, Instant},
};
struct Storage;
impl DocumentStorage for Storage {
    fn save(&mut self, _: &std::path::Path, _: &Document) -> Result<(), String> {
        Err("unused".into())
    }
    fn load(&mut self, _: &std::path::Path) -> Result<Document, String> {
        Err("unused".into())
    }
}
const SIZE: egui::Vec2 = egui::vec2(1360., 860.);
const NUMERIC: [&str; 6] = [
    "Origin X",
    "Origin Y",
    "Origin Z",
    "Width along U",
    "Depth along V",
    "Height along normal",
];

/// One Inspector text input: its label text, value and bounds.
struct Input {
    label: Option<String>,
    value: String,
    rect: egui::Rect,
}
struct Painted {
    inputs: Vec<Input>,
    texts: Vec<(String, egui::Rect)>,
}
impl Painted {
    fn input(&self, label: &str) -> &Input {
        self.inputs
            .iter()
            .find(|i| i.label.as_deref() == Some(label))
            .unwrap_or_else(|| panic!("Inspector input labelled {label}"))
    }
    fn text_in(&self, text: &str, within: egui::Rect) -> Option<egui::Rect> {
        self.texts
            .iter()
            .find(|(t, r)| t == text && within.expand(0.5).contains_rect(*r))
            .map(|(_, r)| *r)
    }
}

struct Harness {
    ctx: egui::Context,
    app: Application,
    time: f64,
}
impl Harness {
    fn box_tool() -> Self {
        let mut app =
            Application::with_document(Box::new(Storage), Document::new("Numeric fields").unwrap())
                .unwrap();
        app.set_model_backends(Arc::new(|| Box::new(NativeBackend::new())));
        app.mode = Mode::Model;
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        crate::theme::configure(&ctx);
        app.invoke(Command::Modeling(ModelingTool::Box), &[], Channel::Ui)
            .unwrap();
        let mut h = Self { ctx, app, time: 0. };
        h.until("native preview", |h| {
            let p = h.app.current();
            p.sketch.model_can_apply(&p.session, &p.resources)
        });
        h
    }
    fn until(&mut self, what: &str, done: impl Fn(&Self) -> bool) {
        let deadline = Instant::now() + Duration::from_secs(10);
        loop {
            self.frame(vec![]).drop_without_applying_deltas();
            if done(self) {
                return;
            }
            assert!(Instant::now() < deadline, "{what} did not arrive");
            std::thread::sleep(Duration::from_millis(1));
        }
    }
    fn frame(&mut self, events: Vec<egui::Event>) -> egui::FullOutput {
        self.time += 1. / 60.;
        let app = &mut self.app;
        self.ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, SIZE)),
                events,
                focused: true,
                time: Some(self.time),
                ..Default::default()
            },
            |ui| app.show(ui),
        )
    }
    fn paint(&mut self) -> Painted {
        let viewport = self.app.current().sketch.viewport_rect.expect("viewport");
        let out = self.frame(vec![]);
        let update = out.platform_output.accesskit_update.as_ref().unwrap();
        let rect = |b: egui::accesskit::Rect| {
            egui::Rect::from_min_max(
                egui::pos2(b.x0 as f32, b.y0 as f32),
                egui::pos2(b.x1 as f32, b.y1 as f32),
            )
        };
        let text_of = |id: &egui::accesskit::NodeId| {
            update
                .nodes
                .iter()
                .find(|(n, _)| n == id)
                .and_then(|(_, node)| node.label().or(node.value()).map(str::to_owned))
        };
        let inputs = update
            .nodes
            .iter()
            .filter(|(_, node)| node.role() == Role::TextInput)
            .filter_map(|(_, node)| {
                Some(Input {
                    label: node.labelled_by().first().and_then(text_of),
                    value: node.value()?.to_owned(),
                    rect: rect(node.bounds()?),
                })
            })
            .filter(|input| !viewport.intersects(input.rect))
            .collect();
        let texts = out
            .shapes
            .iter()
            .filter_map(|clipped| match &clipped.shape {
                egui::Shape::Text(text) => Some((
                    text.galley.text().to_owned(),
                    text.galley.rect.translate(text.pos.to_vec2()),
                )),
                _ => None,
            })
            .collect();
        out.drop_without_applying_deltas();
        Painted { inputs, texts }
    }
    fn press(&mut self, pos: egui::Pos2) {
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: Default::default(),
                },
            ])
            .drop_without_applying_deltas();
        }
        self.frame(vec![]).drop_without_applying_deltas();
    }
    fn key(&mut self, key: egui::Key, modifiers: egui::Modifiers) {
        for pressed in [true, false] {
            self.frame(vec![egui::Event::Key {
                key,
                physical_key: None,
                pressed,
                repeat: false,
                modifiers,
            }])
            .drop_without_applying_deltas();
        }
    }
    fn replace(&mut self, label: &str, text: &str) {
        let pos = self.paint().input(label).rect.center();
        self.press(pos);
        assert!(self.ctx.text_edit_focused(), "{label} focused");
        self.key(egui::Key::A, egui::Modifiers::COMMAND);
        self.frame(vec![egui::Event::Text(text.into())])
            .drop_without_applying_deltas();
    }
}

#[test]
fn box_inspector_fields_are_full_width_with_label_above_and_unit_inside() {
    let mut h = Harness::box_tool();
    let painted = h.paint();
    // The body name spans the Inspector content width (available width).
    let name = painted
        .inputs
        .iter()
        .find(|i| i.value == "Box 1")
        .expect("body name input");
    let content = name.rect;
    for label in NUMERIC {
        let input = painted.input(label);
        assert!(
            (input.rect.left() - content.left()).abs() < 0.5
                && (input.rect.width() - content.width()).abs() < 0.5,
            "{label} input {:?} spans the content width {content:?}",
            input.rect
        );
        let above = painted
            .texts
            .iter()
            .filter(|(t, r)| t == label && r.bottom() <= input.rect.top() + 0.5)
            .map(|(_, r)| *r)
            .find(|r| input.rect.top() - r.bottom() < 12.)
            .unwrap_or_else(|| panic!("{label} label directly above its input"));
        assert!(above.left() >= content.left() - 0.5);
        // The editable value is a plain number; "mm" is painted inside the
        // input at its right edge and is not part of the value.
        assert!(
            input.value.parse::<f64>().is_ok(),
            "{label}: {}",
            input.value
        );
        assert!(!input.value.contains("mm"));
        let unit = painted
            .text_in("mm", input.rect)
            .unwrap_or_else(|| panic!("{label} shows mm inside its input"));
        assert!(
            input.rect.right() - unit.right() < 12.,
            "{label} unit right-aligned"
        );
    }
    // One column: every numeric field sits below the previous one.
    let tops: Vec<f32> = NUMERIC
        .iter()
        .map(|l| painted.input(l).rect.top())
        .collect();
    assert!(tops.windows(2).all(|w| w[0] < w[1]), "{tops:?}");
}

#[test]
fn box_inspector_expression_edit_commits_on_enter_and_escape_keeps_the_tool() {
    let mut h = Harness::box_tool();
    // Escape leaves a focused field without cancelling the tool.
    h.replace("Height along normal", "35");
    h.key(egui::Key::Escape, Default::default());
    h.frame(vec![]).drop_without_applying_deltas();
    assert!(!h.ctx.text_edit_focused());
    assert!(
        h.app.current().sketch.model_active(),
        "Escape kept the tool"
    );
    assert_eq!(h.paint().input("Height along normal").value, "35");
    // An expression typed into the field (no unit text) applies with Enter.
    h.replace("Width along U", "=40 + 2");
    assert_eq!(h.paint().input("Width along U").value, "=40 + 2");
    h.until("expression preview", |h| {
        let p = h.app.current();
        p.sketch.model_can_apply(&p.session, &p.resources)
    });
    h.key(egui::Key::Enter, Default::default());
    h.until("apply", |h| h.app.current().session.revision() > 0);
    let document = h.app.current().session.document();
    let definition = document
        .entities()
        .find_map(|e| model::read_primitive(document, e.id).ok())
        .expect("applied box");
    let model::PrimitiveParameters::Box { size, .. } = definition.parameters else {
        panic!("box parameters");
    };
    let Dimension::Expression { expression } = &size[0] else {
        panic!("width is the typed expression: {:?}", size[0]);
    };
    assert_eq!(expression.text.trim(), "=40 + 2");
    assert_eq!(
        size[2],
        Dimension::Literal {
            value: 35.,
            unit: "mm".into()
        }
    );
}
