//! Full-application UI harness for modelling tool tests (CAD-408/410):
//! native backends, real key and pointer events, bounded readiness waits.
use super::*;
use crate::modeling::native_fixture::NativeBackend;
use std::{
    sync::Arc,
    time::{Duration, Instant},
};

pub(super) struct Storage;
impl DocumentStorage for Storage {
    fn save(&mut self, _: &std::path::Path, _: &Document) -> Result<(), String> {
        Err("unused".into())
    }
    fn load(&mut self, _: &std::path::Path) -> Result<Document, String> {
        Err("unused".into())
    }
}
pub(super) const SIZE: egui::Vec2 = egui::vec2(1360., 860.);
pub(super) struct Harness {
    pub(super) app: Application,
    pub(super) ctx: egui::Context,
}
impl Harness {
    pub(super) fn new(document: Document) -> Self {
        let mut app = Application::with_document(Box::new(Storage), document).unwrap();
        app.pending_files.clear();
        app.set_scene_backends(Arc::new(|| Box::new(NativeBackend::new())));
        app.set_model_backends(Arc::new(|| Box::new(NativeBackend::new())));
        app.mode = Mode::Model;
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        crate::theme::configure(&ctx);
        let mut h = Self { app, ctx };
        h.until("the native model", |app| {
            let p = app.current();
            p.scene
                .is_current(p.session.document().id(), p.session.revision())
        });
        h
    }
    pub(super) fn render(&mut self, events: Vec<egui::Event>) -> egui::FullOutput {
        self.ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, SIZE)),
                events,
                focused: true,
                ..Default::default()
            },
            |ui| self.app.show(ui),
        )
    }
    pub(super) fn frame(&mut self, events: Vec<egui::Event>) {
        self.render(events).drop_without_applying_deltas();
    }
    /// Run frames until `done` holds; fails after a bounded wait.
    pub(super) fn until(&mut self, what: &str, done: impl Fn(&Application) -> bool) {
        let deadline = Instant::now() + Duration::from_secs(20);
        while !done(&self.app) {
            assert!(Instant::now() < deadline, "timed out waiting for {what}");
            self.frame(vec![]);
            std::thread::sleep(Duration::from_millis(1));
        }
    }
    pub(super) fn key(&mut self, key: egui::Key) {
        for pressed in [true, false] {
            self.frame(vec![egui::Event::Key {
                key,
                physical_key: None,
                pressed,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            }]);
        }
    }
    pub(super) fn texts(&mut self) -> Vec<String> {
        let out = self.render(vec![]);
        let texts = out
            .shapes
            .iter()
            .filter_map(|s| match &s.shape {
                egui::Shape::Text(text) => Some(text.galley.text().to_owned()),
                _ => None,
            })
            .collect();
        out.drop_without_applying_deltas();
        texts
    }
    /// Click into the text field drawn right under `label`, replace its text
    /// and press Enter in the same frame, as typing then Enter does.
    pub(super) fn type_into_field_and_enter(&mut self, label: &str, text: &str) {
        let out = self.render(vec![]);
        let caption = out
            .shapes
            .iter()
            .find_map(|s| match &s.shape {
                egui::Shape::Text(t) if t.galley.text() == label => {
                    Some(egui::Rect::from_min_size(t.pos, t.galley.size()))
                }
                _ => None,
            })
            .unwrap_or_else(|| panic!("Missing {label}"));
        let field = out
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .filter(|(_, n)| n.role() == egui::accesskit::Role::TextInput)
            .filter_map(|(_, n)| n.bounds())
            .filter(|b| b.y0 as f32 >= caption.max.y - 1.)
            .min_by(|a, b| a.y0.total_cmp(&b.y0))
            .expect("a text field under the caption");
        out.drop_without_applying_deltas();
        let pos = egui::pos2(
            ((field.x0 + field.x1) * 0.5) as f32,
            ((field.y0 + field.y1) * 0.5) as f32,
        );
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ]);
        }
        assert!(self.ctx.text_edit_focused(), "{label} field focused");
        let key = |key, modifiers| egui::Event::Key {
            key,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers,
        };
        self.frame(vec![key(egui::Key::A, egui::Modifiers::COMMAND)]);
        self.frame(vec![egui::Event::Key {
            key: egui::Key::A,
            physical_key: None,
            pressed: false,
            repeat: false,
            modifiers: egui::Modifiers::COMMAND,
        }]);
        self.frame(vec![
            egui::Event::Text(text.into()),
            key(egui::Key::Enter, egui::Modifiers::NONE),
        ]);
        self.frame(vec![egui::Event::Key {
            key: egui::Key::Enter,
            physical_key: None,
            pressed: false,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        }]);
    }
    /// Click the control whose accessible label or painted text is `label`.
    pub(super) fn click_label(&mut self, label: &str) {
        let pos = self.label_position(label);
        self.click_at(pos);
    }
    /// Right-click the widget or text labelled `label`, e.g. to open a tree
    /// row's context menu.
    pub(super) fn right_click_label(&mut self, label: &str) {
        let pos = self.label_position(label);
        self.press_at(pos, egui::PointerButton::Secondary);
    }
    fn label_position(&mut self, label: &str) -> egui::Pos2 {
        let out = self.render(vec![]);
        let pos = out
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .find_map(|(_, n)| {
                (n.label() == Some(label))
                    .then(|| n.bounds())
                    .flatten()
                    .map(|b| egui::pos2(((b.x0 + b.x1) * 0.5) as f32, ((b.y0 + b.y1) * 0.5) as f32))
            })
            .or_else(|| {
                out.shapes.iter().find_map(|s| match &s.shape {
                    egui::Shape::Text(t) if t.galley.text() == label => {
                        Some(t.pos + t.galley.size() * 0.5)
                    }
                    _ => None,
                })
            })
            .unwrap_or_else(|| panic!("Missing {label}"));
        out.drop_without_applying_deltas();
        pos
    }
    /// Click the `nth` (top to bottom) checkbox whose label starts with `prefix`.
    pub(super) fn click_nth_checkbox(&mut self, prefix: &str, nth: usize) {
        let out = self.render(vec![]);
        let mut boxes: Vec<_> = out
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .filter(|(_, n)| {
                n.role() == egui::accesskit::Role::CheckBox
                    && n.label().is_some_and(|l| l.starts_with(prefix))
            })
            .filter_map(|(_, n)| n.bounds())
            .collect();
        out.drop_without_applying_deltas();
        boxes.sort_by(|a, b| a.y0.total_cmp(&b.y0));
        let b = boxes
            .get(nth)
            .unwrap_or_else(|| panic!("checkbox {nth} for {prefix} among {boxes:?}"));
        self.click_at(egui::pos2(
            ((b.x0 + b.x1) * 0.5) as f32,
            ((b.y0 + b.y1) * 0.5) as f32,
        ));
    }
    pub(super) fn click_at(&mut self, pos: egui::Pos2) {
        self.press_at(pos, egui::PointerButton::Primary);
    }
    fn press_at(&mut self, pos: egui::Pos2, button: egui::PointerButton) {
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ]);
        }
    }
    pub(super) fn revision(&self) -> u64 {
        self.app.current().session.revision()
    }
    pub(super) fn modeling(&self) -> bool {
        self.app.current().sketch.model_active()
    }
}
