//! CAD-371: a focused Inspector field must not block viewport drags.
use super::*;
use crate::modeling::{ModelingTool, native_fixture::NativeBackend};
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

struct Harness {
    ctx: egui::Context,
    app: Application,
    time: f64,
}
impl Harness {
    fn box_tool() -> Self {
        let mut app =
            Application::with_document(Box::new(Storage), Document::new("Viewport focus").unwrap())
                .unwrap();
        app.set_model_backends(Arc::new(|| Box::new(NativeBackend::new())));
        app.mode = Mode::Model;
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        crate::theme::configure(&ctx);
        app.invoke(Command::Modeling(ModelingTool::Box), &[], Channel::Ui)
            .unwrap();
        let mut h = Self { ctx, app, time: 0. };
        h.settle();
        h
    }
    fn settle(&mut self) {
        let until = Instant::now() + Duration::from_secs(5);
        loop {
            self.frame(vec![]).drop_without_applying_deltas();
            let p = self.app.current();
            if p.sketch.model_can_apply(&p.session, &p.resources) {
                break;
            }
            assert!(Instant::now() < until, "native preview did not arrive");
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
    fn viewport(&self) -> egui::Rect {
        self.app.current().sketch.viewport_rect.expect("viewport")
    }
    /// The centre of the first Inspector text field (outside the viewport)
    /// whose value is `value`.
    fn field(&mut self, value: &str) -> egui::Pos2 {
        self.field_where(|text| text == value)
            .unwrap_or_else(|| panic!("Inspector field {value}"))
    }
    fn field_where(&mut self, matches: impl Fn(&str) -> bool) -> Option<egui::Pos2> {
        let viewport = self.viewport();
        let out = self.frame(vec![]);
        let position = out
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .filter(|(_, node)| node.value().is_some_and(&matches))
            .filter_map(|(_, node)| node.bounds())
            .map(|b| egui::pos2(((b.x0 + b.x1) * 0.5) as f32, ((b.y0 + b.y1) * 0.5) as f32))
            .find(|pos| !viewport.contains(*pos));
        out.drop_without_applying_deltas();
        position
    }
    fn press(&mut self, pos: egui::Pos2, button: egui::PointerButton, pressed: bool) {
        self.frame(vec![
            egui::Event::PointerMoved(pos),
            egui::Event::PointerButton {
                pos,
                button,
                pressed,
                modifiers: Default::default(),
            },
        ])
        .drop_without_applying_deltas();
    }
    fn focus(&mut self, value: &str) {
        let pos = self.field(value);
        self.press(pos, egui::PointerButton::Primary, true);
        self.press(pos, egui::PointerButton::Primary, false);
        self.frame(vec![]).drop_without_applying_deltas();
        assert!(self.ctx.text_edit_focused(), "{value} field focused");
    }
    /// A drag that starts in the viewport and moves in several frames.
    fn drag(&mut self, from: egui::Pos2, to: egui::Pos2, button: egui::PointerButton) {
        self.press(from, button, true);
        for step in 1..=4 {
            let pos = from + (to - from) * (step as f32 / 4.);
            self.frame(vec![egui::Event::PointerMoved(pos)])
                .drop_without_applying_deltas();
        }
        self.press(to, button, false);
        self.frame(vec![]).drop_without_applying_deltas();
    }
    fn envelope(&self) -> [[f64; 3]; 8] {
        let p = self.app.current();
        p.sketch.primitive_envelope(&p.session).expect("box draft")
    }
}

#[test]
fn camera_drag_orbits_while_an_inspector_field_is_focused() {
    for field in ["Box 1", "100"] {
        let mut h = Harness::box_tool();
        h.focus(field);
        let before = h.app.current().sketch.camera;
        let start = h.viewport().left_top() + egui::vec2(80., 80.);
        h.drag(
            start,
            start + egui::vec2(60., 25.),
            egui::PointerButton::Secondary,
        );
        let after = h.app.current().sketch.camera;
        assert!(
            (after.yaw - before.yaw).abs() > 1e-4 || (after.pitch - before.pitch).abs() > 1e-4,
            "{field}: right-drag must orbit; {before:?} → {after:?}"
        );
        // Navigation neither applies nor cancels the tool.
        let p = h.app.current();
        assert!(p.sketch.model_active());
        assert_eq!(p.session.revision(), 0);
    }
}

#[test]
fn preview_handle_drag_resizes_while_an_inspector_field_is_focused() {
    let mut h = Harness::box_tool();
    // Stage an Inspector edit first; it must survive the viewport drag.
    h.focus("Box 1");
    h.frame(vec![egui::Event::Text("draft".into())])
        .drop_without_applying_deltas();
    let staged = |h: &mut Harness| h.field_where(|text| text.contains("draft")).is_some();
    assert!(staged(&mut h), "the typed name is staged");
    let rect = h.viewport();
    let camera = h.app.current().sketch.camera;
    let project = |p: [f32; 3]| {
        let [x, y, _] = camera.project(p, rect);
        egui::pos2(x, y)
    };
    let before = h.envelope();
    // The low-u, high-v top corner: away from the corner editor popup.
    let corner = project([0., 60., 20.]);
    let target = corner + (corner - project([100., 60., 20.])) * 0.3;
    h.drag(corner, target, egui::PointerButton::Primary);
    assert_ne!(h.envelope(), before, "the corner drag resized the draft");
    let p = h.app.current();
    assert!(p.sketch.model_active());
    assert_eq!(p.session.revision(), 0);
    assert!(staged(&mut h), "the typed name is kept through the drag");
}

#[test]
fn pan_and_wheel_work_while_focused_but_drags_inside_the_field_stay_with_it() {
    let mut h = Harness::box_tool();
    h.focus("100");
    let start = h.viewport().left_top() + egui::vec2(80., 80.);
    let before = h.app.current().sketch.camera;
    h.drag(
        start,
        start + egui::vec2(40., 30.),
        egui::PointerButton::Middle,
    );
    let panned = h.app.current().sketch.camera;
    assert_ne!(panned.target, before.target, "middle-drag pans");

    h.focus("100");
    h.frame(vec![
        egui::Event::PointerMoved(start),
        egui::Event::MouseWheel {
            unit: egui::MouseWheelUnit::Line,
            delta: egui::vec2(0., 2.),
            modifiers: Default::default(),
            phase: egui::TouchPhase::Move,
        },
    ])
    .drop_without_applying_deltas();
    h.frame(vec![]).drop_without_applying_deltas();
    let zoomed = h.app.current().sketch.camera;
    assert_ne!(
        zoomed.scale, panned.scale,
        "the wheel zooms over the viewport"
    );

    // A selection drag inside the field edits it and never moves the view,
    // even when it ends over the viewport.
    let field = h.field("100");
    let camera = h.app.current().sketch.camera;
    let envelope = h.envelope();
    h.drag(
        field - egui::vec2(8., 0.),
        field + egui::vec2(8., 0.),
        egui::PointerButton::Primary,
    );
    assert!(h.ctx.text_edit_focused(), "the field keeps focus");
    assert_eq!(h.field("100"), field, "the Inspector layout is unchanged");
    h.drag(field, start, egui::PointerButton::Secondary);
    assert_eq!(h.app.current().sketch.camera, camera);
    assert_eq!(h.envelope(), envelope);
    let p = h.app.current();
    assert!(p.sketch.model_active());
    assert_eq!(p.session.revision(), 0);
}
