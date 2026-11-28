//! CAD-438: the Drawings tree (drawings › sheets › views › shown parts),
//! selection between tree and sheet, and per-view part filters driven
//! through the real editor frame.
use super::guided_tests::{box_factory, saved, screen, tagged_document, until};
use super::move_tests::{click_at, context, run};
use super::*;
use multicad_core::{
    drawings::ViewEvaluation,
    project_file::{self, ProjectFile},
};
use std::collections::BTreeSet;

struct Node {
    rect: egui::Rect,
    selected: bool,
    checked: Option<bool>,
}
fn nodes(output: &egui::FullOutput, label: &str) -> Vec<Node> {
    output
        .platform_output
        .accesskit_update
        .as_ref()
        .unwrap()
        .nodes
        .iter()
        .filter(|(_, n)| n.label() == Some(label))
        .filter_map(|(_, n)| {
            let b = n.bounds()?;
            Some(Node {
                rect: egui::Rect::from_min_max(
                    egui::pos2(b.x0 as f32, b.y0 as f32),
                    egui::pos2(b.x1 as f32, b.y1 as f32),
                ),
                // egui reports both a selected row and a ticked box as toggled.
                selected: n.toggled() == Some(egui::accesskit::Toggled::True),
                checked: n.toggled().map(|t| t == egui::accesskit::Toggled::True),
            })
        })
        .collect()
}
/// The one tree row with this label in the left sidebar.
fn tree_row(output: &egui::FullOutput, label: &str) -> Node {
    let mut rows: Vec<_> = nodes(output, label)
        .into_iter()
        .filter(|n| n.rect.left() < 340.)
        .collect();
    assert_eq!(rows.len(), 1, "one tree row {label}");
    rows.remove(0)
}
/// Two tagged blocks, a drawing and the standard views of both.
fn fixture() -> (DrawingEditor, Session, Vec<EntityId>, egui::Context) {
    let (doc, parts, boxes) = tagged_document(&[
        ([200., 50., 80.], &["frame"]),
        ([40., 30., 20.], &["hardware"]),
    ]);
    let mut editor = DrawingEditor::default();
    editor.set_backend_factory(box_factory(boxes));
    let mut session = Session::new(doc).unwrap();
    editor
        .invoke(&mut session, DrawingTool::NewDocument, Channel::Ui)
        .unwrap();
    let ctx = context();
    editor
        .invoke(&mut session, DrawingTool::StandardViews, Channel::Ui)
        .unwrap();
    until(
        &ctx,
        &mut editor,
        &mut session,
        "the standard views",
        |e, _| !e.view_placement_active(),
    );
    settle(&ctx, &mut editor, &mut session);
    (editor, session, parts, ctx)
}
fn settle(ctx: &egui::Context, editor: &mut DrawingEditor, session: &mut Session) {
    until(ctx, editor, session, "a current preview", |e, s| {
        e.preview_outcome(s).is_some()
    });
}
/// The parts drawn in a view of the current preview.
fn drawn(editor: &DrawingEditor, view: EntityId) -> BTreeSet<EntityId> {
    match editor.feed.result.as_ref().unwrap().views.get(&view) {
        Some(ViewEvaluation::Ready { vectors, .. }) => vectors.iter().map(|v| v.source).collect(),
        _ => BTreeSet::new(),
    }
}
fn view_id(editor: &DrawingEditor, name: &str) -> EntityId {
    editor
        .draft
        .as_ref()
        .unwrap()
        .views
        .iter()
        .find(|v| v.name == name)
        .unwrap()
        .id
}

#[test]
fn the_sidebar_is_a_tree_of_drawings_sheets_views_and_shown_parts() {
    let (mut editor, mut session, _, ctx) = fixture();
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    let book = tree_row(&output, "Drawing 1");
    let sheet = tree_row(&output, "1 · Untitled sheet");
    let views: Vec<_> = ["Front view", "Top view", "Right view", "Iso view"]
        .map(|name| tree_row(&output, name))
        .into();
    // Each view lists both parts as toggles: four rows per part name.
    let mut parts: Vec<_> = nodes(&output, "Block 1")
        .into_iter()
        .filter(|n| n.rect.left() < 340.)
        .collect();
    parts.sort_by(|a, b| a.rect.top().total_cmp(&b.rect.top()));
    assert_eq!(parts.len(), 4);
    assert_eq!(
        nodes(&output, "Block 2")
            .iter()
            .filter(|n| n.rect.left() < 340.)
            .count(),
        4
    );
    output.drop_without_applying_deltas();
    // Nested by indent and in order: drawing, sheet, its views, their parts.
    assert!(book.rect.left() < sheet.rect.left());
    assert!(sheet.rect.left() < views[0].rect.left());
    assert!(views[0].rect.left() < parts[0].rect.left());
    assert!(book.rect.top() < sheet.rect.top() && sheet.rect.top() < views[0].rect.top());
    for pair in views.windows(2) {
        assert!(pair[0].rect.top() < pair[1].rect.top());
    }
    assert!(parts[0].rect.top() > views[0].rect.top() && parts[0].rect.top() < views[1].rect.top());
    // The drawing renders view rows with icons: an icon shape precedes each label.
    assert!(
        views
            .iter()
            .all(|v| v.rect.left() > sheet.rect.left() + 10.)
    );
}

#[test]
fn selecting_in_the_tree_selects_on_the_sheet_and_back() {
    let (mut editor, mut session, _, ctx) = fixture();
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    let top = tree_row(&output, "Top view");
    output.drop_without_applying_deltas();
    click_at(&ctx, &mut editor, &mut session, top.rect.center());
    let top_id = view_id(&editor, "Top view");
    assert_eq!(editor.selected_view, Some(top_id));
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    assert!(tree_row(&output, "Top view").selected);
    assert!(!tree_row(&output, "Front view").selected);
    // The inspector opened for that view.
    assert!(
        super::move_tests::texts(&output)
            .iter()
            .any(|t| t == "Parts shown")
    );
    output.drop_without_applying_deltas();
    // Clicking the front view on the sheet selects its tree row.
    let front = editor
        .draft
        .as_ref()
        .unwrap()
        .views
        .iter()
        .find(|v| v.name == "Front view")
        .unwrap()
        .clone();
    let corner = front.paper_point([-100., -40.]);
    let at = screen(&editor, corner);
    click_at(&ctx, &mut editor, &mut session, at);
    assert_eq!(editor.selected_view, Some(front.id));
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    assert!(tree_row(&output, "Front view").selected);
    assert!(!tree_row(&output, "Top view").selected);
    output.drop_without_applying_deltas();
}

#[test]
fn a_view_shows_everything_by_default() {
    let (mut editor, mut session, parts, ctx) = fixture();
    let page = saved(&session, &editor);
    assert!(
        page.part_filters.is_empty(),
        "no filter is stored by default"
    );
    for view in &page.views {
        let listed: BTreeSet<_> = drawings::view_parts(session.document(), view)
            .into_iter()
            .collect();
        assert_eq!(listed, parts.iter().copied().collect());
        assert_eq!(drawn(&editor, view.id), parts.iter().copied().collect());
    }
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    assert!(
        nodes(&output, "Block 1")
            .iter()
            .chain(&nodes(&output, "Block 2"))
            .all(|n| n.checked == Some(true)),
        "every part is ticked"
    );
    output.drop_without_applying_deltas();
}

#[test]
fn tag_include_and_exclude_filter_what_a_view_draws_as_undoable_edits() {
    let (mut editor, mut session, parts, ctx) = fixture();
    let front = view_id(&editor, "Front view");
    let top = view_id(&editor, "Top view");
    // Select the front view; its inspector offers each model tag twice.
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    let row = tree_row(&output, "Front view");
    output.drop_without_applying_deltas();
    click_at(&ctx, &mut editor, &mut session, row.rect.center());
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    let mut frame_chips = nodes(&output, "frame");
    assert_eq!(frame_chips.len(), 2, "include and exclude chips");
    frame_chips.sort_by(|a, b| a.rect.top().total_cmp(&b.rect.top()));
    output.drop_without_applying_deltas();
    // Include "frame": one click, one undoable edit, hardware is not drawn.
    let revision = session.revision();
    click_at(
        &ctx,
        &mut editor,
        &mut session,
        frame_chips[0].rect.center(),
    );
    assert_eq!(session.revision(), revision + 1);
    let filter = saved(&session, &editor)
        .part_filter(front)
        .cloned()
        .unwrap();
    assert_eq!(filter.include_tags, BTreeSet::from(["frame".to_owned()]));
    settle(&ctx, &mut editor, &mut session);
    assert_eq!(drawn(&editor, front), BTreeSet::from([parts[0]]));
    assert_eq!(
        drawn(&editor, top),
        parts.iter().copied().collect(),
        "only this view"
    );
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    let block2: Vec<_> = nodes(&output, "Block 2")
        .into_iter()
        .filter(|n| n.rect.left() < 340.)
        .collect();
    assert_eq!(
        block2.iter().filter(|n| n.checked == Some(false)).count(),
        1
    );
    output.drop_without_applying_deltas();
    session
        .undo(session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    settle(&ctx, &mut editor, &mut session);
    assert!(saved(&session, &editor).part_filters.is_empty());
    assert_eq!(drawn(&editor, front), parts.iter().copied().collect());
    // Exclude "hardware" through the shared command.
    editor
        .set_view_tag_filter(
            &mut session,
            front,
            BTreeSet::new(),
            BTreeSet::from(["hardware".to_owned()]),
        )
        .unwrap();
    settle(&ctx, &mut editor, &mut session);
    assert_eq!(drawn(&editor, front), BTreeSet::from([parts[0]]));
    // Persisted through save and reopen.
    let project = ProjectFile::new(session.document().clone()).unwrap();
    assert!(
        project
            .capabilities()
            .contains(drawings::PART_FILTER_CAPABILITY)
    );
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let reopened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let page = drawings::read(&reopened.project.document, editor.book.unwrap())
        .unwrap()
        .pages
        .remove(0);
    assert_eq!(
        page.part_filter(front).unwrap().exclude_tags,
        BTreeSet::from(["hardware".to_owned()])
    );
    // Editing a part's canonical tags in Model changes what the view shows.
    let prepared = multicad_core::materials::prepare_tags(
        &session,
        session.revision(),
        Actor::local(Channel::Ui),
        &[parts[0]],
        &multicad_core::materials::TagEdit::Add(BTreeSet::from(["hardware".to_owned()])),
    )
    .unwrap();
    session.commit(prepared).unwrap();
    settle(&ctx, &mut editor, &mut session);
    assert!(matches!(
        editor.feed.result.as_ref().unwrap().views.get(&front),
        Some(ViewEvaluation::Unresolved { .. })
    ));
}

#[test]
fn toggling_a_part_in_the_tree_hides_it_from_that_view_as_one_undo_step() {
    let (mut editor, mut session, parts, ctx) = fixture();
    let front = view_id(&editor, "Front view");
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    let mut toggles: Vec<_> = nodes(&output, "Block 2")
        .into_iter()
        .filter(|n| n.rect.left() < 340.)
        .collect();
    toggles.sort_by(|a, b| a.rect.top().total_cmp(&b.rect.top()));
    output.drop_without_applying_deltas();
    // The first Block 2 row belongs to the first view, the front view.
    let revision = session.revision();
    let before = session.document().clone();
    click_at(&ctx, &mut editor, &mut session, toggles[0].rect.center());
    assert_eq!(session.revision(), revision + 1, "one undoable edit");
    let filter = saved(&session, &editor)
        .part_filter(front)
        .cloned()
        .unwrap();
    assert_eq!(filter.hidden, BTreeSet::from([parts[1]]));
    settle(&ctx, &mut editor, &mut session);
    assert_eq!(
        drawn(&editor, front),
        BTreeSet::from([parts[0]]),
        "hidden part not drawn"
    );
    assert_eq!(editor.selected_view, Some(front));
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    let mut toggles: Vec<_> = nodes(&output, "Block 2")
        .into_iter()
        .filter(|n| n.rect.left() < 340.)
        .collect();
    toggles.sort_by(|a, b| a.rect.top().total_cmp(&b.rect.top()));
    assert_eq!(toggles[0].checked, Some(false));
    assert!(toggles[1..].iter().all(|n| n.checked == Some(true)));
    output.drop_without_applying_deltas();
    // Ticking it again shows it again and stores nothing.
    click_at(&ctx, &mut editor, &mut session, toggles[0].rect.center());
    assert!(saved(&session, &editor).part_filters.is_empty());
    // Undo twice restores the original page.
    for _ in 0..2 {
        session
            .undo(session.revision(), Actor::local(Channel::Ui))
            .unwrap();
    }
    assert_eq!(session.document(), &before);
    // A part hidden by its tag can be ticked on against the filter.
    editor
        .set_view_tag_filter(
            &mut session,
            front,
            BTreeSet::new(),
            BTreeSet::from(["hardware".to_owned()]),
        )
        .unwrap();
    editor
        .set_part_shown(&mut session, front, parts[1], true)
        .unwrap();
    settle(&ctx, &mut editor, &mut session);
    assert_eq!(drawn(&editor, front), parts.iter().copied().collect());
}
