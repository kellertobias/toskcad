//! CAD-436 with real geometry: New drawing on the empty sheet, Add view and
//! one click on the sheet give a correct front view of a native box.
use super::tests::{add_box, preview};
use super::*;
use multicad_core::{
    commands::{Channel, Session},
    drawings::{self, ViewEvaluation},
};
use multicad_ui::{
    drawings::{DrawingEditor, DrawingTool},
    egui,
};
use std::time::{Duration, Instant};

struct Ui {
    ctx: egui::Context,
    start: Instant,
}
impl Ui {
    fn frame(
        &self,
        editor: &mut DrawingEditor,
        session: &mut Session,
        events: Vec<egui::Event>,
    ) -> egui::FullOutput {
        self.ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1200., 850.),
                )),
                time: Some(self.start.elapsed().as_secs_f64()),
                focused: true,
                events,
                ..Default::default()
            },
            |ui| editor.show(ui, session),
        )
    }
    fn click(&self, editor: &mut DrawingEditor, session: &mut Session, at: egui::Pos2) {
        for pressed in [true, false] {
            self.frame(
                editor,
                session,
                vec![
                    egui::Event::PointerMoved(at),
                    egui::Event::PointerButton {
                        pos: at,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: Default::default(),
                    },
                ],
            )
            .drop_without_applying_deltas();
        }
    }
    fn until(
        &self,
        editor: &mut DrawingEditor,
        session: &mut Session,
        what: &str,
        done: impl Fn(&DrawingEditor) -> bool,
    ) {
        let deadline = Instant::now() + Duration::from_secs(120);
        while !done(editor) {
            assert!(Instant::now() < deadline, "Timed out waiting for {what}");
            self.frame(editor, session, vec![])
                .drop_without_applying_deltas();
            std::thread::sleep(Duration::from_millis(2));
        }
    }
}

#[test]
fn new_drawing_add_view_and_one_click_give_a_correct_native_front_view() {
    let mut session = Session::new(Document::new("Guided front view").unwrap()).unwrap();
    // A 20 × 10 × 5 mm box at the origin: its front view is 20 wide, 5 high.
    add_box(&mut session, "Block", 0.);
    let mut editor = DrawingEditor::default();
    editor.set_backend_factory(Arc::new(|| {
        Box::new(Backend::with_cache_and_backend(None, GeometryBackend::Rust))
    }));
    let ui = Ui {
        ctx: egui::Context::default(),
        start: Instant::now(),
    };
    ui.ctx.enable_accesskit();
    let start = session.revision();
    // 1. New drawing on the empty sheet area.
    let output = ui.frame(&mut editor, &mut session, vec![]);
    let button = output
        .platform_output
        .accesskit_update
        .as_ref()
        .unwrap()
        .nodes
        .iter()
        .find_map(|(_, n)| {
            (n.label() == Some("New drawing"))
                .then(|| n.bounds())
                .flatten()
        })
        .expect("New drawing on the empty sheet");
    output.drop_without_applying_deltas();
    ui.click(
        &mut editor,
        &mut session,
        egui::pos2(
            ((button.x0 + button.x1) / 2.) as f32,
            ((button.y0 + button.y1) / 2.) as f32,
        ),
    );
    assert_eq!(session.revision(), start + 1);
    // 2. Add view (the ribbon's command).
    editor
        .invoke(&mut session, DrawingTool::AddView, Channel::Ui)
        .unwrap();
    ui.until(
        &mut editor,
        &mut session,
        "the fitted native preview",
        |e| {
            e.next_step_hint()
                .is_some_and(|h| h.starts_with("Click the sheet to place the front view"))
        },
    );
    // 3. One click on the sheet.
    let at = editor.sheet_point([130., 90.]).unwrap();
    ui.click(&mut editor, &mut session, at);
    ui.until(&mut editor, &mut session, "the applied view", |e| {
        !e.view_placement_active()
    });
    assert_eq!(session.revision(), start + 2);
    assert!(!editor.dirty());
    let book = drawings::books(session.document()).next().unwrap();
    let page = drawings::read(session.document(), book).unwrap().pages[0].clone();
    let view = page.views[0].clone();
    assert!((view.position_mm[0] - 130.).abs() < 0.5 && (view.position_mm[1] - 90.).abs() < 0.5);
    // The native front view is the box's X × Z silhouette, centered on it.
    let mut backend = Backend::with_cache_and_backend(None, GeometryBackend::Rust);
    let evaluation = preview(&mut backend, &session, page.id);
    let Some(ViewEvaluation::Ready { vectors, .. }) = evaluation.views.get(&view.id) else {
        panic!("front view resolves")
    };
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for p in vectors.iter().flat_map(|v| &v.preview_mm) {
        for i in 0..2 {
            lo[i] = lo[i].min(p[i]);
            hi[i] = hi[i].max(p[i]);
        }
    }
    let near = |a: f64, b: f64| (a - b).abs() < 1e-6;
    assert!(
        near(hi[0] - lo[0], 20.) && near(hi[1] - lo[1], 5.),
        "{lo:?} {hi:?}"
    );
    assert!(
        near(lo[0] + hi[0], 0.) && near(lo[1] + hi[1], 0.),
        "centered"
    );
    // At the largest standard scale that fits the sheet.
    let layout = drawings::layout_page(&page, 1, 1).unwrap();
    let fit = drawings::fit_view(&view, &layout, lo, hi, drawings::FitMode::Standard).unwrap();
    assert_eq!(view.scale, fit.scale);
    assert_eq!(drawings::scale_text(view.scale), "10:1");
}
