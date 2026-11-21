use super::*;
use multicad_core::{document::Document, project_file::ProjectResources};
use std::time::{Duration, Instant};
#[allow(dead_code)]
#[path = "../../../../multicad-app/examples/shared_materials_review.rs"]
mod shared_fixture;

#[test]
fn hex_color_editor_keeps_valid_rgb_when_switching_modes_or_typing_a_partial_value() {
    let ctx = egui::Context::default();
    let mut rgb = [12, 34, 56];
    for (hex, input, expected) in [
        (true, Some("#336699"), [51, 102, 153]),
        (true, Some("#33"), [51, 102, 153]),
        (false, None, [51, 102, 153]),
        (true, None, [51, 102, 153]),
    ] {
        ctx.run_ui(egui::RawInput::default(), |ui| {
            let mode = ui.id().with("rgb-mode");
            let field = mode.with("hex-value");
            ctx.data_mut(|data| {
                data.insert_temp(mode, hex);
                if let Some(value) = input {
                    data.insert_temp(field, value.to_owned());
                } else {
                    data.remove::<String>(field);
                }
            });
            if input.is_some() {
                ctx.memory_mut(|memory| memory.request_focus(field));
            }
            rgb_fields(ui, &mut rgb);
        })
        .drop_without_applying_deltas();
        assert_eq!(rgb, expected);
    }
}

fn editor(s: &Session) -> LibraryEditor {
    let mut editor = LibraryEditor::default();
    editor.set_backend(Box::new(
        crate::modeling::native_fixture::NativeBackend::new(),
    ));
    editor.open(s);
    editor
}
fn preview(editor: &mut LibraryEditor, s: &Session) -> Result<(), String> {
    let resources = ProjectResources::default();
    editor.refresh_preview(s, &resources, Channel::Ui)?;
    let deadline = Instant::now() + Duration::from_secs(15);
    while editor.pending() {
        assert!(
            Instant::now() < deadline,
            "native definition preview timed out"
        );
        editor.poll(s, &resources);
        std::thread::sleep(Duration::from_millis(1));
    }
    if let Some(error) = &editor.error {
        return Err(error.clone());
    }
    assert!(editor.can_apply(s, &resources));
    Ok(())
}
fn apply(editor: &mut LibraryEditor, s: &mut Session) -> Result<(), String> {
    preview(editor, s)?;
    editor.apply(s, &Default::default(), Channel::Ui)
}
fn set_stock(editor: &mut LibraryEditor, f: impl FnOnce(&mut Stock)) {
    let d = editor.draft.as_mut().unwrap();
    let EntityData::Extension {
        definition: wire,
        inputs,
    } = &mut d.entity.data
    else {
        panic!("stock draft");
    };
    let LibraryDefinition::Stock(mut stock) = serde_json::from_value(wire.data.clone()).unwrap()
    else {
        panic!("stock draft");
    };
    f(&mut stock);
    let definition = LibraryDefinition::Stock(stock);
    wire.data = serde_json::to_value(&definition).unwrap();
    *inputs = definition.links();
}
#[test]
fn inline_stock_rejects_invalid_dimensions_and_keeps_the_draft_editable() {
    let session = Session::new(Document::new("Stock corrections").unwrap()).unwrap();
    let mut editor = LibraryEditor::for_model_stock(&session, ToolStock::Sheet);
    set_stock(&mut editor, |stock| {
        if let StockKind::Sheet { thickness_mm, .. } = &mut stock.specification {
            *thickness_mm = 0.;
        }
    });
    assert!(matches!(
        editor.prepare_model_stock(&session),
        Err(error) if error.contains("positive sheet dimensions")
    ));
    assert_eq!(session.document().entities().count(), 0);
    // Keep the same draft and correct the field.
    set_stock(&mut editor, |stock| {
        if let StockKind::Sheet { thickness_mm, .. } = &mut stock.specification {
            *thickness_mm = 1.5;
        }
    });
    let (_, prepared) = editor.prepare_model_stock(&session).unwrap();
    assert_eq!(prepared.document().entities().count(), 2);
}
#[test]
fn native_material_preview_is_neutral_and_apply_and_history_keep_one_identity() {
    let mut s = Session::new(Document::new("Project materials").unwrap()).unwrap();
    let before = s.document().clone();
    let mut editor = editor(&s);
    editor.create(&s, 0);
    let id = editor.draft.as_ref().unwrap().entity.id;
    editor.draft.as_mut().unwrap().entity.name = "Reusable oak".into();
    assert!(
        editor
            .apply(&mut s, &Default::default(), Channel::Ui)
            .is_err()
    );
    preview(&mut editor, &s).unwrap();
    assert_eq!(s.document(), &before);
    assert_eq!(s.revision(), 0);
    editor
        .apply(&mut s, &Default::default(), Channel::Ui)
        .unwrap();
    assert_eq!(s.revision(), 1);
    assert_eq!(s.document().get(id).unwrap().name, "Reusable oak");
    assert!(!editor.dirty());
    assert_eq!(s.activity().last().unwrap().actor.channel, Channel::Ui);
    s.undo(1, Actor::local(Channel::Ui)).unwrap();
    assert_eq!(s.document(), &before);
    editor.open(&s);
    assert!(
        editor.draft.is_none(),
        "Undo must not retain an authorable stale definition"
    );
    s.redo(s.revision(), Actor::local(Channel::Ui)).unwrap();
    assert!(s.document().get(id).is_some());
}
#[test]
fn new_stock_in_empty_project_creates_companion_material_in_one_native_checked_edit() {
    let mut s = Session::new(Document::new("Empty stock").unwrap()).unwrap();
    let before = s.document().clone();
    let mut editor = editor(&s);
    editor.create_stock(&s, 1);
    let d = editor.draft.as_ref().unwrap();
    let stock = d.entity.id;
    let material = d.companion.as_ref().unwrap().id;
    assert_ne!(stock, material);
    assert_eq!(stock_material(&d.entity), Some(material));
    apply(&mut editor, &mut s).unwrap();
    assert_eq!(s.revision(), 1);
    assert_eq!(s.document().entities().count(), 2);
    assert!(diagnostics(s.document()).unwrap().is_empty());
    s.undo(1, Actor::local(Channel::Ui)).unwrap();
    assert_eq!(s.document(), &before);
}
#[test]
fn material_first_sequence_creates_panel_sheet_and_profile_in_the_same_catalog() {
    let mut s = Session::new(Document::new("Shared catalog").unwrap()).unwrap();
    let mut editor = editor(&s);
    editor.create(&s, 0);
    let material = editor.draft.as_ref().unwrap().entity.id;
    apply(&mut editor, &mut s).unwrap();
    for kind in 0..3 {
        editor.create_stock(&s, kind);
        assert!(editor.draft.as_ref().unwrap().companion.is_none());
        let id = editor.draft.as_ref().unwrap().entity.id;
        apply(&mut editor, &mut s).unwrap();
        let Some(LibraryDefinition::Stock(stock)) = read(s.document().get(id).unwrap()).unwrap()
        else {
            panic!("committed stock")
        };
        assert_eq!(stock.material, material);
        assert_eq!(stock_index(&stock.specification), kind);
    }
    assert_eq!(s.revision(), 4);
    assert_eq!(s.document().entities().count(), 4);
    assert!(diagnostics(s.document()).unwrap().is_empty());
}
#[test]
fn invalid_missing_and_stale_definitions_refuse_without_source_edits() {
    let mut s = Session::new(Document::new("Invalid definitions").unwrap()).unwrap();
    let mut editor = editor(&s);
    editor.create(&s, 0);
    if let EntityData::Material { density_kg_m3, .. } =
        &mut editor.draft.as_mut().unwrap().entity.data
    {
        *density_kg_m3 = Some(-10.);
    }
    let before = s.document().clone();
    assert!(preview(&mut editor, &s).is_err());
    assert_eq!(s.document(), &before);
    editor.create_stock(&s, 0);
    let stock = editor.draft.as_ref().unwrap().entity.id;
    editor.draft.as_mut().unwrap().entity.name = "Birch ply".into();
    editor.draft.as_mut().unwrap().companion = None;
    set_stock(&mut editor, |s| s.material = stock);
    let problem = preview(&mut editor, &s).unwrap_err();
    assert!(
        problem.contains("Birch ply") && problem.contains("manufacturing material"),
        "{problem}"
    );
    assert_eq!(s.document(), &before);
    editor.create(&s, 0);
    s.apply(multicad_core::commands::Request {
        expected_revision: 0,
        actor: Actor::local(Channel::Ui),
        label: "Rename".into(),
        edits: vec![multicad_core::document::Mutation::RenameDocument(
            "New source".into(),
        )],
    })
    .unwrap();
    assert!(
        preview(&mut editor, &s)
            .unwrap_err()
            .to_lowercase()
            .contains("stale")
    );
    assert!(s.document().entities().next().is_none());
    editor.discard(&s);
    assert!(!editor.dirty());
}
#[test]
fn checked_preview_cannot_authorize_changed_draft_resources_session_or_channel() {
    let mut s = Session::new(Document::new("Preview guards").unwrap()).unwrap();
    let mut editor = editor(&s);
    editor.create(&s, 0);
    preview(&mut editor, &s).unwrap();
    let original = s.document().clone();
    editor.draft.as_mut().unwrap().entity.name = "Changed since preview".into();
    assert!(!editor.can_apply(&s, &Default::default()));
    assert!(editor.preview_scene(&s, &Default::default()).is_none());
    editor.poll(&s, &Default::default());
    assert!(!editor.pending());
    preview(&mut editor, &s).unwrap();
    let mut resources = ProjectResources::default();
    resources
        .metadata
        .insert("review".into(), serde_json::json!("changed"));
    assert!(!editor.can_apply(&s, &resources));
    let foreign = Session::new(s.document().clone()).unwrap();
    assert!(!editor.can_apply(&foreign, &Default::default()));
    assert!(
        editor
            .apply(&mut s, &Default::default(), Channel::Mcp)
            .is_err()
    );
    assert_eq!(s.document(), &original);
    editor.discard(&s);
    assert_eq!(s.document(), &original);
}
#[test]
fn applying_used_panel_stock_validates_geometry_and_commits_without_a_second_click() {
    let (mut s, part) = crate::modeling::panel_fixture::flat_panel();
    let stock = s
        .document()
        .entities()
        .find(|e| matches!(read(e), Ok(Some(LibraryDefinition::Stock(_)))))
        .unwrap()
        .clone();
    let stock_id = stock.id;
    let before = s.document().clone();
    let revision = s.revision();
    let mut editor = editor(&s);
    editor.select(&s, stock);
    set_stock(&mut editor, |s| {
        if let StockKind::Panel { thickness_mm, .. } = &mut s.specification {
            *thickness_mm = 24.;
        }
    });
    let resources = ProjectResources::default();
    editor.confirm(&mut s, &resources, Channel::Ui).unwrap();
    assert_eq!(editor.confirm_on_ready, Some(Channel::Ui));
    assert!(editor.pending());
    let deadline = Instant::now() + Duration::from_secs(15);
    while editor.pending() {
        assert!(
            Instant::now() < deadline,
            "native definition check timed out"
        );
        editor.poll(&s, &resources);
        std::thread::sleep(Duration::from_millis(1));
    }
    let bounds = editor
        .preview_scene(&s, &Default::default())
        .unwrap()
        .native_instance_bounds(part)
        .unwrap();
    assert!((bounds[1][2] - bounds[0][2] - 24.).abs() < 1e-6);
    assert_eq!(s.document(), &before);
    egui::Context::default()
        .run_ui(Default::default(), |ui| {
            editor.show_inline(ui, &mut s, &resources, true);
        })
        .drop_without_applying_deltas();
    assert_eq!(
        s.revision(),
        revision + 1,
        "one source edit after confirmation"
    );
    assert!(!editor.dirty());
    assert!(s.document().get(part).is_some());
    assert!(s.document().get(stock_id).is_some());
    s.undo(s.revision(), Actor::local(Channel::Ui)).unwrap();
    assert_eq!(s.document(), &before);
    editor.open(&s);
    assert_eq!(
        &editor.draft.as_ref().unwrap().entity,
        s.document().get(stock_id).unwrap()
    );
    assert!(!editor.dirty());
}
#[test]
fn a_shared_stock_edit_that_breaks_dependent_geometry_cannot_apply() {
    let (s, _, _) = crate::modeling::panel_fixture::carcass();
    let stock = s
        .document()
        .entities()
        .find(|e| matches!(read(e), Ok(Some(LibraryDefinition::Stock(_)))))
        .unwrap()
        .clone();
    let before = s.document().clone();
    let mut editor = editor(&s);
    editor.select(&s, stock);
    set_stock(&mut editor, |s| {
        if let StockKind::Panel { thickness_mm, .. } = &mut s.specification {
            *thickness_mm = 400.;
        }
    });
    let error = preview(&mut editor, &s).unwrap_err();
    assert!(!error.is_empty());
    assert!(!editor.can_apply(&s, &Default::default()));
    assert_eq!(s.document(), &before);
}
#[test]
fn a_display_scene_without_native_source_metadata_cannot_authorize_apply() {
    use crate::viewport::{Scene, SceneBackend, SceneReply, SceneStamp};
    struct NoReceipt(Option<SceneReply>);
    impl SceneBackend for NoReceipt {
        fn request(&mut self, _: &Document, stamp: SceneStamp) -> Result<(), String> {
            self.0 = Some(SceneReply {
                stamp,
                result: Ok(std::sync::Arc::new(Scene::new(stamp, vec![], vec![])?)),
            });
            Ok(())
        }
        fn poll(&mut self) -> Option<SceneReply> {
            self.0.take()
        }
    }
    let s = Session::new(Document::new("No receipt").unwrap()).unwrap();
    let mut editor = editor(&s);
    editor.set_backend(Box::new(NoReceipt(None)));
    editor.create(&s, 0);
    let error = preview(&mut editor, &s).unwrap_err();
    assert!(error.contains("source-bound"), "{error}");
    assert!(!editor.can_apply(&s, &Default::default()));
    assert_eq!(s.revision(), 0);
}

#[test]
fn project_definition_edits_drive_wood_sheet_and_profile_geometry_and_shared_drawing_hatches() {
    let mut fixture = shared_fixture::fixture().unwrap();
    let original = fixture.session.document().clone();
    let count = original.entities().count();
    shared_fixture::check(&fixture.session, &fixture, false).unwrap();
    let mut editor = editor(&fixture.session);
    for id in fixture.stocks.into_iter().chain([fixture.hatch]) {
        let entity = fixture.session.document().get(id).unwrap().clone();
        editor.select(&fixture.session, entity);
        editor.draft.as_mut().unwrap().entity =
            shared_fixture::changed_definition(&fixture.session, id).unwrap();
        let source = fixture.session.document().clone();
        preview(&mut editor, &fixture.session).unwrap();
        assert_eq!(fixture.session.document(), &source);
        editor
            .apply(&mut fixture.session, &Default::default(), Channel::Ui)
            .unwrap();
        assert_eq!(fixture.session.document().entities().count(), count);
    }
    shared_fixture::check(&fixture.session, &fixture, true).unwrap();
    for _ in 0..4 {
        fixture
            .session
            .undo(fixture.session.revision(), Actor::local(Channel::Ui))
            .unwrap();
    }
    assert_eq!(fixture.session.document(), &original);
    shared_fixture::check(&fixture.session, &fixture, false).unwrap();
    for _ in 0..4 {
        fixture
            .session
            .redo(fixture.session.revision(), Actor::local(Channel::Ui))
            .unwrap();
    }
    shared_fixture::check(&fixture.session, &fixture, true).unwrap();
}

#[test]
fn changing_a_used_stocks_material_reports_the_dependent_assignment_instead_of_breaking_it() {
    let fixture = shared_fixture::fixture().unwrap();
    let original = fixture.session.document().clone();
    let mut editor = editor(&fixture.session);
    let stock = fixture
        .session
        .document()
        .get(fixture.stocks[0])
        .unwrap()
        .clone();
    editor.select(&fixture.session, stock);
    let Some(LibraryDefinition::Stock(other)) =
        read(fixture.session.document().get(fixture.stocks[1]).unwrap()).unwrap()
    else {
        panic!("sheet stock")
    };
    set_stock(&mut editor, |stock| stock.material = other.material);
    let error = preview(&mut editor, &fixture.session).unwrap_err();
    assert!(error.contains("material"), "{error}");
    assert!(!editor.can_apply(&fixture.session, &Default::default()));
    assert_eq!(fixture.session.document(), &original);
}

#[test]
fn hollow_section_widget_can_author_a_real_solid_bar_without_replacing_stock_or_part() {
    let mut fixture = shared_fixture::fixture().unwrap();
    let mut editor = editor(&fixture.session);
    editor.select(
        &fixture.session,
        fixture
            .session
            .document()
            .get(fixture.stocks[2])
            .unwrap()
            .clone(),
    );
    let ctx = egui::Context::default();
    ctx.enable_accesskit();
    crate::theme::configure(&ctx);
    let mut frame = |editor: &mut LibraryEditor, events| {
        ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1000., 800.),
                )),
                events,
                focused: true,
                ..Default::default()
            },
            |ui| editor.show_inline(ui, &mut fixture.session, &Default::default(), true),
        )
    };
    for _ in 0..3 {
        frame(&mut editor, vec![]).drop_without_applying_deltas();
    }
    let output = frame(&mut editor, vec![]);
    let bounds = output
        .platform_output
        .accesskit_update
        .as_ref()
        .unwrap()
        .nodes
        .iter()
        .filter_map(|(_, node)| {
            (node.label() == Some("Hollow section"))
                .then(|| node.bounds())
                .flatten()
        })
        .max_by(|a, b| (a.y1 - a.y0).total_cmp(&(b.y1 - b.y0)))
        .expect("hollow checkbox");
    let pos = egui::pos2(
        ((bounds.x0 + bounds.x1) / 2.) as f32,
        ((bounds.y0 + bounds.y1) / 2.) as f32,
    );
    output.drop_without_applying_deltas();
    for pressed in [true, false] {
        frame(
            &mut editor,
            vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: Default::default(),
                },
            ],
        )
        .drop_without_applying_deltas();
    }
    let Some(LibraryDefinition::Stock(stock)) =
        read(&editor.draft.as_ref().unwrap().entity).unwrap()
    else {
        panic!("profile")
    };
    assert!(matches!(
        stock.specification,
        StockKind::Profile { wall_mm: None, .. }
    ));
    apply(&mut editor, &mut fixture.session).unwrap();
    let native = multicad_model::Compiler::default()
        .compile(fixture.session.document(), fixture.session.revision())
        .unwrap();
    let volume = native
        .solid(fixture.parts[2])
        .unwrap()
        .volume_mm3()
        .unwrap();
    assert!((volume - 120000.).abs() < 1e-5);
    assert!(fixture.session.document().get(fixture.stocks[2]).is_some());
}

/// Render the inline editor once and return (accessible labels, unclipped
/// painted texts with their rects, focused accessible label).
fn render(
    ctx: &egui::Context,
    editor: &mut LibraryEditor,
    s: &mut Session,
    size: egui::Vec2,
    events: Vec<egui::Event>,
) -> (Vec<String>, Vec<(String, egui::Rect)>, Option<String>) {
    let output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size)),
            events,
            focused: true,
            ..Default::default()
        },
        |ui| editor.show_inline(ui, s, &Default::default(), true),
    );
    let update = output.platform_output.accesskit_update.as_ref().unwrap();
    let labels = update
        .nodes
        .iter()
        .filter_map(|(_, node)| node.label().map(str::to_owned))
        .collect();
    let focused = update
        .nodes
        .iter()
        .find(|(id, _)| *id == update.focus)
        .and_then(|(_, node)| node.label().map(str::to_owned));
    let texts = output
        .shapes
        .iter()
        .filter_map(|shape| match &shape.shape {
            egui::Shape::Text(text) => {
                let rect = egui::Rect::from_min_size(text.pos, text.galley.size());
                (shape.clip_rect.contains_rect(rect)
                    && egui::Rect::from_min_size(egui::Pos2::ZERO, size).contains_rect(rect))
                .then(|| (text.galley.text().to_owned(), rect))
            }
            _ => None,
        })
        .collect();
    output.drop_without_applying_deltas();
    (labels, texts, focused)
}

#[test]
fn materials_tree_groups_stock_with_its_material_and_says_metal_sheet_stock() {
    let mut fixture = shared_fixture::fixture().unwrap();
    let before = fixture.session.document().clone();
    let mut editor = editor(&fixture.session);
    let ctx = egui::Context::default();
    ctx.enable_accesskit();
    crate::theme::configure(&ctx);
    // A narrow window: the tree keeps its minimum width and wraps rows.
    let size = egui::vec2(640., 560.);
    let sheet = fixture
        .session
        .document()
        .get(fixture.stocks[1])
        .unwrap()
        .clone();
    editor.select(&fixture.session, sheet);
    let mut out = Default::default();
    for _ in 0..3 {
        out = render(&ctx, &mut editor, &mut fixture.session, size, vec![]);
    }
    let (labels, texts, _) = out;
    for heading in [
        "Materials",
        "Panel stock",
        "Metal sheet stock",
        "Profile stock",
        "Print stock",
        "Appearances",
        "Hatches",
    ] {
        assert!(labels.iter().any(|l| l == heading), "missing {heading}");
    }
    assert!(!labels.iter().any(|l| l == "Sheet stock"));
    let painted = |needle: &str| texts.iter().find(|(t, _)| t == needle).map(|(_, r)| *r);
    let tree = size.x * 0.24 + 1.;
    for text in [
        "Stock definitions",
        "Appearance and drawing",
        "Metal sheet stock",
        "Steel sheet 2",
        "Material: Shared steel",
        "Used by 2 stock definitions",
        "Used by 1 stock definition",
    ] {
        let rect = painted(text).unwrap_or_else(|| panic!("{text} not painted unclipped"));
        assert!(rect.left() < tree, "{text} belongs in the tree");
    }
    // The stock editor uses the same wording for its kind.
    assert!(
        texts
            .iter()
            .any(|(t, r)| t == "Metal sheet" && r.left() > tree)
    );
    assert_eq!(fixture.session.document(), &before, "navigation only");
}

#[test]
fn keyboard_tab_and_enter_select_a_tree_definition() {
    let mut fixture = shared_fixture::fixture().unwrap();
    let mut editor = editor(&fixture.session);
    let ctx = egui::Context::default();
    ctx.enable_accesskit();
    crate::theme::configure(&ctx);
    let size = egui::vec2(900., 640.);
    for _ in 0..3 {
        render(&ctx, &mut editor, &mut fixture.session, size, vec![]);
    }
    let key = |key| egui::Event::Key {
        key,
        physical_key: None,
        pressed: true,
        repeat: false,
        modifiers: Default::default(),
    };
    let mut reached = false;
    for _ in 0..12 {
        let (_, _, focused) = render(
            &ctx,
            &mut editor,
            &mut fixture.session,
            size,
            vec![key(egui::Key::Tab)],
        );
        if focused.as_deref() == Some("Birch") {
            reached = true;
            break;
        }
    }
    assert!(reached, "Tab reaches the Birch material row");
    render(
        &ctx,
        &mut editor,
        &mut fixture.session,
        size,
        vec![key(egui::Key::Enter)],
    );
    assert_eq!(
        editor.draft.as_ref().map(|d| d.entity.name.as_str()),
        Some("Birch")
    );
}

#[test]
fn draft_stock_rows_name_a_new_companion_and_broken_links_as_missing() {
    let s = Session::new(Document::new("Rows").unwrap()).unwrap();
    let mut editor = editor(&s);
    editor.create_stock(&s, 1);
    let draft = &editor.draft.as_ref().unwrap().entity;
    assert_eq!(draft.name, "New metal sheet stock");
    assert_eq!(
        row_detail(s.document(), &BTreeMap::new(), draft).as_deref(),
        Some("Material: new (applied with this stock)")
    );
    let committed = BTreeMap::from([(draft.id, stock_material(draft).unwrap())]);
    assert_eq!(
        row_detail(s.document(), &committed, draft).as_deref(),
        Some("Material: missing")
    );
    assert_eq!(ToolStock::Sheet.label(), "Create metal sheet stock…");
}

fn click_label(
    ctx: &egui::Context,
    editor: &mut LibraryEditor,
    s: &mut Session,
    size: egui::Vec2,
    label: &str,
) {
    let output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size)),
            focused: true,
            ..Default::default()
        },
        |ui| editor.show_inline(ui, s, &Default::default(), true),
    );
    let bounds = output
        .platform_output
        .accesskit_update
        .as_ref()
        .unwrap()
        .nodes
        .iter()
        .find_map(|(_, node)| {
            (node.label() == Some(label))
                .then(|| node.bounds())
                .flatten()
        })
        .unwrap_or_else(|| panic!("no control {label}"));
    output.drop_without_applying_deltas();
    let pos = egui::pos2(
        ((bounds.x0 + bounds.x1) / 2.) as f32,
        ((bounds.y0 + bounds.y1) / 2.) as f32,
    );
    for pressed in [true, false] {
        render(
            ctx,
            editor,
            s,
            size,
            vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: Default::default(),
                },
            ],
        );
    }
}

#[test]
fn stock_routes_to_its_material_color_which_drives_preview_swatch_and_model_and_round_trips() {
    use multicad_core::project_file::{self, ProjectFile};
    let mut fixture = shared_fixture::fixture().unwrap();
    let s = &mut fixture.session;
    let ctx = egui::Context::default();
    ctx.enable_accesskit();
    crate::theme::configure(&ctx);
    let size = egui::vec2(1100., 800.);
    let mut editor = editor(s);
    let panel = s.document().get(fixture.stocks[0]).unwrap().clone();
    let birch = stock_material(&panel).unwrap();
    editor.select(s, panel);
    for _ in 0..3 {
        render(&ctx, &mut editor, s, size, vec![]);
    }
    let (labels, texts, _) = render(&ctx, &mut editor, s, size, vec![]);
    assert!(texts.iter().any(|(t, _)| t == "Color (from material)"));
    assert!(
        texts.iter().any(|(t, _)| t == "#D2AD78"),
        "Birch color shown on stock"
    );
    assert!(labels.iter().any(|l| l == "Edit material color: Birch"));
    // A dirty stock draft cannot jump away and lose its edits.
    editor.draft.as_mut().unwrap().entity.name = "Renamed panel".into();
    click_label(&ctx, &mut editor, s, size, "Edit material color: Birch");
    assert_eq!(editor.draft.as_ref().unwrap().entity.name, "Renamed panel");
    editor.discard(s);
    click_label(&ctx, &mut editor, s, size, "Edit material color: Birch");
    assert_eq!(editor.draft.as_ref().unwrap().entity.id, birch);
    let (labels, texts, _) = render(&ctx, &mut editor, s, size, vec![]);
    assert!(texts.iter().any(|(t, _)| t == "Appearance color"));
    assert!(labels.iter().any(|l| l == "Hex") && labels.iter().any(|l| l == "RGB"));
    assert!(
        texts
            .iter()
            .any(|(t, _)| t.starts_with("Shared: changing this color updates 1 stock definition")),
        "shared effect is explained"
    );
    // Edit the color and preview it natively before Apply.
    let new_color = [40, 120, 60, 255];
    if let EntityData::Material { color_rgba, .. } = &mut editor.draft.as_mut().unwrap().entity.data
    {
        *color_rgba = new_color;
    }
    let before = s.document().clone();
    preview(&mut editor, s).unwrap();
    assert!(editor.preview_scene(s, &Default::default()).is_some());
    assert_eq!(s.document(), &before);
    editor.apply(s, &Default::default(), Channel::Ui).unwrap();
    // The model scene colors each Part with `Resolver::resolve(part)`
    // (multicad-app scene publication); the Birch panel part follows.
    assert_eq!(
        resolve(s.document(), fixture.parts[0])
            .unwrap()
            .appearance
            .color_rgba,
        new_color
    );
    // Swatch and stock preview resolve the same color; stock has no own color.
    editor.open(s);
    render(&ctx, &mut editor, s, size, vec![]);
    assert_eq!(editor.swatches.get(&birch), Some(&new_color));
    let Some(LibraryDefinition::Stock(stock)) =
        read(s.document().get(fixture.stocks[0]).unwrap()).unwrap()
    else {
        panic!("stock")
    };
    assert_eq!(stock.material, birch);
    let applied = s.document().clone();
    s.undo(s.revision(), Actor::local(Channel::Ui)).unwrap();
    assert_eq!(s.document(), &before);
    s.redo(s.revision(), Actor::local(Channel::Ui)).unwrap();
    assert_eq!(s.document(), &applied);
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(applied.clone()).unwrap(),
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
    assert!(matches!(
        opened.project.document.get(birch).unwrap().data,
        EntityData::Material { color_rgba, .. } if color_rgba == new_color
    ));
}

#[test]
fn linked_appearance_overrides_are_named_and_shown_in_stock_color() {
    let mut s = Session::new(Document::new("Linked appearance").unwrap()).unwrap();
    let mut editor = editor(&s);
    editor.create(&s, 2);
    if let EntityData::Extension { definition, .. } =
        &mut editor.draft.as_mut().unwrap().entity.data
    {
        definition.data = serde_json::to_value(LibraryDefinition::Appearance(Appearance {
            color_rgba: [10, 20, 30, 255],
            ..Appearance::default()
        }))
        .unwrap();
    }
    editor.draft.as_mut().unwrap().entity.name = "Lacquer".into();
    let appearance = editor.draft.as_ref().unwrap().entity.id;
    apply(&mut editor, &mut s).unwrap();
    editor.create(&s, 0);
    let material = editor.draft.as_ref().unwrap().entity.id;
    editor.draft.as_mut().unwrap().style = Some((
        EntityId::new(),
        MaterialStyle {
            material,
            appearance: Some(LibraryRef::capture(s.document(), appearance).unwrap()),
            hatch: None,
        },
    ));
    apply(&mut editor, &mut s).unwrap();
    let ctx = egui::Context::default();
    ctx.enable_accesskit();
    crate::theme::configure(&ctx);
    let size = egui::vec2(1100., 800.);
    editor.select(&s, s.document().get(material).unwrap().clone());
    for _ in 0..2 {
        render(&ctx, &mut editor, &mut s, size, vec![]);
    }
    let (_, texts, _) = render(&ctx, &mut editor, &mut s, size, vec![]);
    assert!(
        texts
            .iter()
            .any(|(t, _)| t.contains("linked appearance 'Lacquer'")),
        "the override is not silent"
    );
    assert_eq!(editor.swatches.get(&material), Some(&[10, 20, 30, 255]));
    editor.create_stock(&s, 0);
    apply(&mut editor, &mut s).unwrap();
    let (_, texts, _) = render(&ctx, &mut editor, &mut s, size, vec![]);
    assert!(
        texts.iter().any(|(t, _)| t == "#0A141E"),
        "stock shows the canonical color"
    );
}
