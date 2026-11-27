//! CAD-433: documented single-key tool shortcuts start their tools through the
//! full application key handling, and never fire while a text field is focused.
use super::ui_harness::Harness;
use super::*;
use crate::modeling::ModelingTool;
use multicad_core::model::{self, PrimitiveParameters, PrimitiveRequest};
use std::collections::BTreeSet;

fn box_part() -> (Document, EntityId) {
    let mut session = Session::new(Document::new("Shortcuts").unwrap()).unwrap();
    let request = PrimitiveRequest::new(
        "Panel",
        PrimitiveParameters::Box {
            size: [model::length(30.), model::length(20.), model::length(10.)],
            frame: Default::default(),
        },
    );
    session
        .commit(
            model::prepare_primitive(&session, 0, Actor::local(Channel::Ui), &request)
                .unwrap()
                .prepared,
        )
        .unwrap();
    (session.document().clone(), request.body)
}
fn selected(h: &mut Harness, part: EntityId) {
    h.app
        .current_mut()
        .sketch
        .set_tree_selection_for_tests(BTreeSet::from([part]));
    h.frame(vec![]);
}

#[test]
fn m_starts_move_on_the_selected_part() {
    let (document, part) = box_part();
    let mut h = Harness::new(document);
    selected(&mut h, part);
    h.key(egui::Key::M);
    let model = h.app.current().sketch.model_for_tests();
    assert!(model.placement_active(), "M must start Move / Copy");
    assert_eq!(model.placement_targets(), Some(&BTreeSet::from([part])));
    assert_eq!(h.revision(), 0, "starting Move commits nothing");
}

/// Existing Model bindings after the ways a user actually reaches the canvas.
#[test]
fn model_keys_dispatch_after_clicking_the_viewport_or_a_tree_row() {
    let (document, _) = box_part();
    for (how, key) in [
        ("viewport", egui::Key::B),
        ("viewport", egui::Key::C),
        ("tree", egui::Key::B),
    ] {
        let mut h = Harness::new(document.clone());
        match how {
            "viewport" => h.click_at(egui::pos2(700., 700.)),
            _ => h.click_label("Panel"),
        }
        h.frame(vec![]);
        assert!(!h.modeling(), "{how}: no tool before the key");
        h.key(key);
        assert!(h.modeling(), "{how}: {key:?} must start its tool");
    }
}

fn chord(h: &mut Harness, key: egui::Key, modifiers: egui::Modifiers) {
    for pressed in [true, false] {
        h.frame(vec![egui::Event::Key {
            key,
            physical_key: None,
            pressed,
            repeat: false,
            modifiers,
        }]);
    }
}

#[test]
fn m_without_a_selection_asks_for_targets_instead_of_doing_nothing() {
    let (document, _) = box_part();
    let mut h = Harness::new(document);
    h.key(egui::Key::M);
    assert!(!h.app.current().sketch.model_for_tests().placement_active());
    let (message, _) = h.app.status();
    assert!(message.to_lowercase().contains("select"), "{message}");
    assert_eq!(h.revision(), 0);
}

#[test]
fn model_tool_keys_start_their_tools() {
    let (document, part) = box_part();
    for (key, tool) in [
        (egui::Key::B, ModelingTool::Box),
        (egui::Key::C, ModelingTool::Cylinder),
        (egui::Key::M, ModelingTool::Move),
    ] {
        let mut h = Harness::new(document.clone());
        selected(&mut h, part);
        h.key(key);
        assert!(h.modeling(), "{key:?} must start {tool:?}");
        assert!(
            h.app.tool_active(Command::Modeling(tool)),
            "{key:?} must start {tool:?}"
        );
        assert_eq!(h.revision(), 0, "{key:?} commits nothing");
    }
}

#[test]
fn sketch_tool_keys_start_their_drawing_tools() {
    let (document, _) = box_part();
    let mut h = Harness::new(document);
    h.app
        .invoke(Command::NewSketch, &["xy"], Channel::Ui)
        .unwrap();
    h.frame(vec![]);
    for (key, tool) in [
        (egui::Key::E, crate::sketch_editor::DrawTool::Line),
        (egui::Key::R, crate::sketch_editor::DrawTool::Rectangle),
        (egui::Key::F, crate::sketch_editor::DrawTool::Circle),
    ] {
        h.key(key);
        assert_eq!(h.app.current().sketch.tool, Some(tool), "{key:?}");
    }
}

#[test]
fn tool_keys_never_fire_while_a_text_field_has_focus() {
    let (document, part) = box_part();
    let mut h = Harness::new(document);
    selected(&mut h, part);
    // ⌘F focuses the tool search field.
    chord(&mut h, egui::Key::F, egui::Modifiers::COMMAND);
    h.frame(vec![]);
    assert!(h.ctx.text_edit_focused(), "search field focused");
    for key in [egui::Key::M, egui::Key::B, egui::Key::C] {
        h.frame(vec![
            egui::Event::Key {
                key,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            },
            egui::Event::Text(key.name().to_lowercase()),
        ]);
        h.frame(vec![egui::Event::Key {
            key,
            physical_key: None,
            pressed: false,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        }]);
        assert!(
            !h.modeling(),
            "{key:?} typed into the field must not start a tool"
        );
    }
    assert!(h.ctx.text_edit_focused());
    assert_eq!(h.revision(), 0);
}

#[test]
fn tooltips_and_the_shortcut_list_show_the_keys() {
    let (document, _) = box_part();
    let mut h = Harness::new(document);
    h.frame(vec![]);
    let help = h
        .app
        .command_help(Command::Modeling(ModelingTool::Move), "Move / Copy", &h.ctx);
    assert_eq!(help, "Move / Copy · M");
    assert!(
        h.app
            .command_help(Command::Modeling(ModelingTool::Box), "Box", &h.ctx)
            .ends_with("· B")
    );
    assert!(h.app.shortcuts.bindings().iter().any(|b| {
        b.chord == crate::shortcuts::Chord::plain(egui::Key::M)
            && b.action == crate::shortcuts::Action::Command(Command::Modeling(ModelingTool::Move))
    }));
}

#[test]
fn an_unavailable_tool_key_reports_why_nothing_started() {
    let (document, part) = box_part();
    let mut h = Harness::new(document);
    selected(&mut h, part);
    // An unapplied parameter draft blocks modelling tools.
    h.app.current_mut().editor.rows[0].name = "width".into();
    h.app.current_mut().editor.rows[0].expression = "600".into();
    h.app.current_mut().editor.ensure_blank_row();
    assert!(h.app.current().editor.dirty());
    h.key(egui::Key::M);
    assert!(!h.modeling());
    let (message, _) = h.app.status();
    assert!(
        message.contains("not available") && message.contains("(M)"),
        "{message}"
    );
}
