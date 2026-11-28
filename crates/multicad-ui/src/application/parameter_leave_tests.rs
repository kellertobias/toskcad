//! CAD-417 UI tests: real clicks on the mode tabs, Project tabs and the
//! prompt's buttons. Window close is exempt (user decision 2026-11-27).
use super::*;

struct Storage;
impl DocumentStorage for Storage {
    fn save(&mut self, _: &std::path::Path, _: &Document) -> Result<(), String> {
        Err("No file writes in parameter-leave tests".into())
    }
    fn load(&mut self, _: &std::path::Path) -> Result<Document, String> {
        Err("unused".into())
    }
}

const SIZE: egui::Vec2 = egui::vec2(1360., 860.);

struct Ui {
    app: Application,
    ctx: egui::Context,
}

impl Ui {
    /// A saved project in Project → Parameters.
    fn new() -> Self {
        let mut app =
            Application::with_document(Box::new(Storage), Document::new("Leave").unwrap()).unwrap();
        app.pending_files.clear();
        app.current_mut().session.mark_saved(0).unwrap();
        app.mode = Mode::Parameters;
        app.tab = project::PARAMETERS_TAB;
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        crate::theme::configure(&ctx);
        let mut ui = Self { app, ctx };
        ui.frame(vec![]);
        ui
    }
    fn run(&mut self, events: Vec<egui::Event>, close: bool) -> egui::FullOutput {
        let mut input = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, SIZE)),
            events,
            focused: true,
            ..Default::default()
        };
        if close {
            input
                .viewports
                .entry(egui::ViewportId::ROOT)
                .or_default()
                .events
                .push(egui::ViewportEvent::Close);
        }
        self.ctx.run_ui(input, |ui| self.app.show(ui))
    }
    fn frame(&mut self, events: Vec<egui::Event>) {
        self.run(events, false).drop_without_applying_deltas();
    }
    /// Type a pending parameter row, as the table does.
    fn edit_parameter(&mut self) {
        let editor = &mut self.app.current_mut().editor;
        editor.ensure_blank_row();
        let row = editor.rows.last_mut().unwrap();
        row.name = "leave_width".into();
        row.unit = "mm".into();
        row.expression = "42".into();
        editor.ensure_blank_row();
        // The table refreshes its preview over the next frames.
        for _ in 0..3 {
            self.frame(vec![]);
        }
        assert!(self.app.current().editor.rows_dirty());
    }
    fn node(&mut self, label: &str) -> Option<(egui::Rect, bool)> {
        let out = self.run(vec![], false);
        let found = out
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .find_map(|(_, n)| {
                let b = n.bounds()?;
                (n.label() == Some(label)).then(|| {
                    (
                        egui::Rect::from_min_max(
                            egui::pos2(b.x0 as f32, b.y0 as f32),
                            egui::pos2(b.x1 as f32, b.y1 as f32),
                        ),
                        !n.is_disabled(),
                    )
                })
            });
        out.drop_without_applying_deltas();
        found
    }
    fn click(&mut self, label: &str) {
        let (rect, enabled) = self
            .node(label)
            .unwrap_or_else(|| panic!("Missing {label}"));
        assert!(enabled, "{label} disabled");
        let pos = rect.center();
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
        self.frame(vec![]);
    }
    fn prompt_painted(&mut self) -> bool {
        let out = self.run(vec![], false);
        let visible = out
            .shapes
            .iter()
            .any(|s| matches!(&s.shape, egui::Shape::Text(t) if t.galley.text() == PROMPT_TITLE));
        out.drop_without_applying_deltas();
        visible
    }
    /// The prompt is open and painted; its area needs one sizing pass
    /// before it draws, so wait (bounded) for the painted title.
    fn prompt_visible(&mut self) -> bool {
        if self.app.parameter_leave.is_none() {
            return self.prompt_painted();
        }
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
        while !self.prompt_painted() {
            if std::time::Instant::now() > deadline || self.app.parameter_leave.is_none() {
                return false;
            }
        }
        true
    }
    /// The pending dot painted inside the Project mode tab.
    fn marker_visible(&mut self) -> bool {
        let (tab, _) = self.node("Project").expect("Project mode tab");
        let out = self.run(vec![], false);
        let found = out.shapes.iter().any(|s| {
            matches!(&s.shape, egui::Shape::Circle(c)
                if c.fill == PENDING_MARKER && tab.contains(c.center))
        });
        out.drop_without_applying_deltas();
        found
    }
    fn revision(&self) -> u64 {
        self.app.current().session.revision()
    }
}

fn parameter_count(app: &Application) -> usize {
    app.current()
        .session
        .document()
        .entities()
        .filter(|e| matches!(e.data, EntityData::Parameter { .. }))
        .count()
}

#[test]
fn switching_modes_without_pending_edits_never_prompts() {
    let mut ui = Ui::new();
    assert!(!ui.marker_visible());
    ui.click("Model");
    assert_eq!(ui.app.mode, Mode::Model);
    assert!(ui.app.parameter_leave.is_none());
    assert!(!ui.prompt_visible());
    ui.click("Project");
    assert_eq!(ui.app.mode, Mode::Parameters);
    ui.click("Drawings");
    assert_eq!(ui.app.mode, Mode::Drawings);
    assert!(!ui.prompt_visible());
}

#[test]
fn pending_marker_shows_on_the_project_tab_until_the_edit_is_resolved() {
    let mut ui = Ui::new();
    assert!(!ui.marker_visible());
    ui.edit_parameter();
    assert!(ui.marker_visible());
    ui.click("Discard Edit");
    assert!(!ui.app.current().editor.rows_dirty());
    assert!(!ui.marker_visible());
}

#[test]
fn apply_and_continue_applies_one_undo_step_then_switches_mode() {
    let mut ui = Ui::new();
    let before = ui.app.current().session.document().clone();
    ui.edit_parameter();
    ui.click("Model");
    assert_eq!(ui.app.mode, Mode::Parameters, "the prompt holds the switch");
    assert_eq!(
        ui.app.parameter_leave,
        Some(ParameterLeave::Command(Command::Mode(Mode::Model)))
    );
    assert!(ui.prompt_visible());
    for label in [APPLY_LABEL, DISCARD_LABEL, STAY_LABEL] {
        assert!(
            ui.node(label).is_some_and(|(_, enabled)| enabled),
            "{label}"
        );
    }
    ui.click(APPLY_LABEL);
    assert_eq!(ui.app.mode, Mode::Model);
    assert!(!ui.prompt_visible());
    assert!(!ui.app.current().editor.rows_dirty());
    assert_eq!(ui.revision(), 1);
    assert_eq!(parameter_count(&ui.app), 1);
    ui.app.invoke(Command::Undo, &[], Channel::Ui).unwrap();
    assert_eq!(ui.app.current().session.document(), &before);
}

#[test]
fn discard_changes_restores_the_table_then_switches_mode() {
    let mut ui = Ui::new();
    ui.edit_parameter();
    ui.click("Motion");
    assert!(ui.prompt_visible());
    ui.click(DISCARD_LABEL);
    assert_eq!(ui.app.mode, Mode::Motion);
    assert_eq!(ui.revision(), 0);
    assert_eq!(parameter_count(&ui.app), 0);
    assert!(!ui.app.current().editor.rows_dirty());
    assert!(!ui.prompt_visible());
}

#[test]
fn stay_keeps_the_draft_in_the_parameter_table() {
    let mut ui = Ui::new();
    ui.edit_parameter();
    ui.click("Model");
    assert!(ui.prompt_visible());
    ui.click(STAY_LABEL);
    assert_eq!(ui.app.mode, Mode::Parameters);
    assert_eq!(ui.app.tab, project::PARAMETERS_TAB);
    assert!(ui.app.current().editor.rows_dirty());
    assert_eq!(ui.revision(), 0);
    assert!(!ui.prompt_visible());
    assert!(ui.marker_visible());
    // Escape stays too.
    ui.click("Drawings");
    assert!(ui.prompt_visible());
    for pressed in [true, false] {
        ui.frame(vec![egui::Event::Key {
            key: egui::Key::Escape,
            physical_key: None,
            pressed,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        }]);
    }
    assert!(!ui.prompt_visible());
    assert_eq!(ui.app.mode, Mode::Parameters);
    assert!(ui.app.current().editor.rows_dirty());
}

#[test]
fn another_project_tab_prompts_and_continues_there() {
    let mut ui = Ui::new();
    ui.edit_parameter();
    ui.click("General");
    assert!(ui.prompt_visible());
    assert_eq!(ui.app.tab, project::PARAMETERS_TAB);
    ui.click(DISCARD_LABEL);
    assert_eq!(ui.app.tab, 0);
    ui.click("Parameters");
    ui.edit_parameter();
    ui.click("Materials");
    assert!(
        ui.prompt_visible(),
        "Materials is reachable through the prompt"
    );
    ui.click(APPLY_LABEL);
    assert_eq!(ui.app.tab, project::MATERIALS_TAB);
    assert_eq!(ui.revision(), 1);
}

#[test]
fn closing_the_project_with_pending_edits_prompts_first() {
    let mut ui = Ui::new();
    ui.edit_parameter();
    ui.app
        .invoke(Command::CloseProject, &[], Channel::Ui)
        .unwrap();
    assert_eq!(
        ui.app.parameter_leave,
        Some(ParameterLeave::Command(Command::CloseProject))
    );
    assert!(ui.prompt_visible());
    assert!(
        !ui.app.welcome,
        "the project stays open while the prompt asks"
    );
    ui.click(STAY_LABEL);
    assert!(!ui.app.welcome);
    assert!(ui.app.current().editor.rows_dirty());
    ui.app
        .invoke(Command::CloseProject, &[], Channel::Ui)
        .unwrap();
    ui.click(DISCARD_LABEL);
    // The saved, discarded project closes to the start screen.
    assert!(ui.app.welcome);
}

#[test]
fn closing_the_window_skips_the_parameter_prompt_for_the_quit_dialog() {
    // User decision 2026-11-27: window close is exempt from the
    // parameter-leave prompt and goes straight to the Quit dialog.
    let mut ui = Ui::new();
    ui.edit_parameter();
    ui.run(vec![], true).drop_without_applying_deltas();
    ui.frame(vec![]);
    assert!(ui.app.parameter_leave.is_none());
    assert!(!ui.prompt_painted());
    assert!(ui.app.close_confirm, "the Quit dialog asks about the draft");
    assert!(ui.app.current().editor.rows_dirty());
}
