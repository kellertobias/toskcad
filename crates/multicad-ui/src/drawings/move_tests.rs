//! CAD-437: placed views move by dragging them on the sheet, and Apply page
//! is a live control only while page edits are pending. Driven through the
//! real editor frame with headless pointer input.
use super::tests::{publish_editor, view_editor};
use super::*;
use multicad_core::project_file::{self, ProjectFile};

pub(super) fn run(
    ctx: &egui::Context,
    editor: &mut DrawingEditor,
    session: &mut Session,
    events: Vec<egui::Event>,
) -> egui::FullOutput {
    ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1200., 850.),
            )),
            focused: true,
            events,
            ..Default::default()
        },
        |ui| editor.show(ui, session),
    )
}
pub(super) fn texts(output: &egui::FullOutput) -> Vec<String> {
    output
        .shapes
        .iter()
        .filter_map(|shape| match &shape.shape {
            egui::Shape::Text(t) => Some(t.galley.text().to_owned()),
            _ => None,
        })
        .collect()
}
pub(super) fn frame(
    ctx: &egui::Context,
    editor: &mut DrawingEditor,
    session: &mut Session,
    events: Vec<egui::Event>,
) -> Vec<String> {
    let output = run(ctx, editor, session, events);
    let shown = texts(&output);
    output.drop_without_applying_deltas();
    shown
}
/// Accessible controls with this label, as (enabled, center).
pub(super) fn controls(output: &egui::FullOutput, label: &str) -> Vec<(bool, egui::Pos2)> {
    output
        .platform_output
        .accesskit_update
        .as_ref()
        .expect("accesskit enabled")
        .nodes
        .iter()
        .filter(|(_, n)| n.label() == Some(label))
        .filter_map(|(_, n)| {
            let b = n.bounds()?;
            Some((
                !n.is_disabled(),
                egui::pos2(((b.x0 + b.x1) / 2.) as f32, ((b.y0 + b.y1) / 2.) as f32),
            ))
        })
        .collect()
}
pub(super) fn context() -> egui::Context {
    let ctx = egui::Context::default();
    ctx.enable_accesskit();
    crate::theme::configure(&ctx);
    ctx
}
pub(super) fn press(pressed: bool, at: egui::Pos2) -> egui::Event {
    egui::Event::PointerButton {
        pos: at,
        button: egui::PointerButton::Primary,
        pressed,
        modifiers: Default::default(),
    }
}
pub(super) fn click_at(
    ctx: &egui::Context,
    editor: &mut DrawingEditor,
    session: &mut Session,
    at: egui::Pos2,
) {
    for pressed in [true, false] {
        frame(
            ctx,
            editor,
            session,
            vec![egui::Event::PointerMoved(at), press(pressed, at)],
        );
    }
}
/// Run frames until the page preview is current or the budget is spent.
pub(super) fn settle(ctx: &egui::Context, editor: &mut DrawingEditor, session: &mut Session) {
    for _ in 0..40 {
        frame(ctx, editor, session, vec![]);
        if editor.preview_outcome(session).is_some() {
            break;
        }
        publish_editor(editor, session);
    }
    frame(ctx, editor, session, vec![]);
}
/// Press at `from`, move in small steps by `by` screen points, release.
pub(super) fn drag(
    ctx: &egui::Context,
    editor: &mut DrawingEditor,
    session: &mut Session,
    from: egui::Pos2,
    by: egui::Vec2,
) {
    frame(
        ctx,
        editor,
        session,
        vec![egui::Event::PointerMoved(from), press(true, from)],
    );
    for step in 1..=8 {
        let at = from + by * (step as f32 / 8.);
        frame(ctx, editor, session, vec![egui::Event::PointerMoved(at)]);
    }
    frame(ctx, editor, session, vec![press(false, from + by)]);
}
fn screen(editor: &DrawingEditor, mm: [f64; 2]) -> egui::Pos2 {
    let (paper, scale) = editor.last_paper.expect("paper drawn");
    paper.min + egui::vec2(mm[0] as f32, mm[1] as f32) * scale
}
fn saved(session: &Session, editor: &DrawingEditor) -> drawings::Page {
    drawings::read(session.document(), editor.book.unwrap())
        .unwrap()
        .pages
        .into_iter()
        .find(|p| Some(p.id) == editor.selected)
        .unwrap()
}
fn position(page: &drawings::Page, id: EntityId) -> [f64; 2] {
    page.views.iter().find(|v| v.id == id).unwrap().position_mm
}
/// A saved projection at `position` drawing the test line (0,0)–(600,400)
/// model mm at 1:10, so it covers 60 × 40 paper mm up and right of it.
fn placed(
    editor: &mut DrawingEditor,
    session: &mut Session,
    source: EntityId,
    position: [f64; 2],
) -> EntityId {
    editor.add_view(drawings::ViewKind::Projection).unwrap();
    let view = editor.draft.as_mut().unwrap().views.last_mut().unwrap();
    view.source = drawings::ViewSource::Model {
        entities: vec![source],
    };
    view.position_mm = position;
    view.scale = drawings::Scale {
        numerator: 1,
        denominator: 10,
    };
    let id = view.id;
    publish_editor(editor, session);
    editor
        .invoke(session, DrawingTool::Apply, Channel::Ui)
        .unwrap();
    id
}

#[test]
fn dragging_a_placed_view_moves_it_as_one_undo_step_that_survives_save_and_reopen() {
    let (mut editor, mut session, source) = view_editor(false);
    let id = placed(&mut editor, &mut session, source, [80., 120.]);
    let ctx = context();
    settle(&ctx, &mut editor, &mut session);
    let before = session.document().clone();
    let revision = session.revision();
    // Grab the middle of the drawn line and drag it 40 × 20 paper mm.
    let (_, scale) = editor.last_paper.unwrap();
    let grab = screen(&editor, [110., 100.]);
    drag(
        &ctx,
        &mut editor,
        &mut session,
        grab,
        egui::vec2(40., 20.) * scale,
    );
    settle(&ctx, &mut editor, &mut session);
    let moved = position(&saved(&session, &editor), id);
    assert!(
        (moved[0] - 120.).abs() < 0.6 && (moved[1] - 140.).abs() < 0.6,
        "Dragging the view must move its saved anchor, got {moved:?}"
    );
    assert_eq!(session.revision(), revision + 1, "One canonical edit");
    assert!(!editor.dirty(), "The move is applied, not left as a draft");
    // The new position survives save and reopen.
    let project = ProjectFile::new(session.document().clone()).unwrap();
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let reopened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let book = drawings::read(&reopened.project.document, editor.book.unwrap()).unwrap();
    assert_eq!(position(&book.pages[0], id), moved);
    // One Undo restores the original placement.
    session
        .undo(session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(session.document(), &before);
    settle(&ctx, &mut editor, &mut session);
    assert_eq!(position(editor.draft.as_ref().unwrap(), id), [80., 120.]);
    // A plain click still selects the view without moving it.
    let at = screen(&editor, [110., 100.]);
    click_at(&ctx, &mut editor, &mut session, at);
    assert_eq!(editor.selected_view, Some(id));
    assert_eq!(session.document(), &before);
}

#[test]
fn projected_views_stay_aligned_with_their_base_while_dragging() {
    let (mut editor, mut session, source) = view_editor(false);
    let front = placed(&mut editor, &mut session, source, [40., 170.]);
    let top = placed(&mut editor, &mut session, source, [40., 90.]);
    // Top shares the front's right axis: it is projected from the front.
    let page = editor.draft.as_mut().unwrap();
    page.views[0].frame.v = [0., 0., 1.];
    page.projections.push(drawings::ProjectedView {
        view: top,
        base: front,
    });
    publish_editor(&mut editor, &session);
    editor
        .invoke(&mut session, DrawingTool::Apply, Channel::Ui)
        .unwrap();
    let ctx = context();
    settle(&ctx, &mut editor, &mut session);
    let (_, scale) = editor.last_paper.unwrap();
    // Dragging the projected view sideways keeps it in the front's column.
    let grab = screen(&editor, [70., 70.]);
    drag(
        &ctx,
        &mut editor,
        &mut session,
        grab,
        egui::vec2(50., 10.) * scale,
    );
    settle(&ctx, &mut editor, &mut session);
    let page = saved(&session, &editor);
    let t = position(&page, top);
    assert!((t[0] - 40.).abs() < 1e-9, "Top left its column: {t:?}");
    assert!(
        (t[1] - 100.).abs() < 0.6,
        "Top slides along its free axis: {t:?}"
    );
    // Dragging the base carries the projected view along.
    let grab = screen(&editor, [70., 150.]);
    drag(
        &ctx,
        &mut editor,
        &mut session,
        grab,
        egui::vec2(30., -5.) * scale,
    );
    settle(&ctx, &mut editor, &mut session);
    let page = saved(&session, &editor);
    let (f, t) = (position(&page, front), position(&page, top));
    assert!(
        (f[0] - 70.).abs() < 0.6 && (f[1] - 165.).abs() < 0.6,
        "{f:?}"
    );
    assert_eq!(t[0], f[0], "Projection keeps the base column");
    assert!((t[1] - 95.).abs() < 0.6, "{t:?}");
}

#[test]
fn apply_page_is_shown_only_for_pending_edits_and_confirms_what_it_applied() {
    let (mut editor, mut session, _) = view_editor(false);
    let ctx = context();
    settle(&ctx, &mut editor, &mut session);
    // Nothing pending: no Apply page control at all, live or dead.
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    assert!(
        controls(&output, "Apply page").is_empty(),
        "Apply page must not be offered without pending page settings"
    );
    output.drop_without_applying_deltas();
    // A pending paper change shows a live Apply page and says what is pending.
    editor.draft.as_mut().unwrap().settings.paper = drawings::PaperSize::A3;
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    let apply = controls(&output, "Apply page");
    assert_eq!(apply.len(), 1, "One Apply page control");
    assert!(apply[0].0, "Apply page is enabled for a valid pending edit");
    assert!(
        texts(&output)
            .iter()
            .any(|t| t.contains("Unapplied page edits"))
    );
    output.drop_without_applying_deltas();
    let revision = session.revision();
    click_at(&ctx, &mut editor, &mut session, apply[0].1);
    let shown = frame(&ctx, &mut editor, &mut session, vec![]);
    assert_eq!(session.revision(), revision + 1);
    assert_eq!(
        saved(&session, &editor).settings.paper,
        drawings::PaperSize::A3
    );
    assert!(
        shown.iter().any(|t| t.contains("Page applied")),
        "Apply page gives visible feedback: {shown:?}"
    );
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    assert!(controls(&output, "Apply page").is_empty());
    output.drop_without_applying_deltas();
}
