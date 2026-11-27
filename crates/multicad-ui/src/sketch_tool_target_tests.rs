//! Opening cabinet machining, joint and print-aid tools keeps the selected
//! targets; an unsupported or ambiguous selection is explained and leaves
//! the choice open instead of falling back to an unrelated Part.
use super::*;
use crate::{
    modeling::{ModelingTool, native_fixture::NativeBackend, panel_fixture},
    viewport::SceneFeed,
};
use multicad_core::{
    commands::{Actor, Channel, Session},
    model::{self, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, length},
};
use std::collections::BTreeSet;

fn editor() -> SketchEditor {
    let mut editor = SketchEditor::default();
    editor.set_model_backend(Box::new(NativeBackend::new()));
    editor
}
fn open(editor: &mut SketchEditor, session: &Session, tool: ModelingTool) {
    editor
        .activate_model_tool(
            session,
            &Default::default(),
            &SceneFeed::default(),
            tool,
            Channel::Ui,
        )
        .unwrap();
}
fn add_box(session: &mut Session, name: &str, x: f64) -> EntityId {
    let request = PrimitiveRequest::new(
        name,
        PrimitiveParameters::Box {
            size: [length(30.), length(30.), length(10.)],
            frame: PrimitiveFrame {
                origin: [length(x), length(0.), length(0.)],
                ..Default::default()
            },
        },
    );
    let p = model::prepare_primitive(
        session,
        session.revision(),
        Actor::local(Channel::Ui),
        &request,
    )
    .unwrap();
    session.commit(p.prepared).unwrap();
    request.body
}

#[test]
fn wood_presets_seed_the_selected_feature_panel_and_cancel_source_neutrally() {
    use crate::wood_navigation::Preset;
    let (session, left, _) = panel_fixture::carcass();
    let before = session.document().clone();
    let feature = model::panel_feature(&before, left).unwrap();
    let mut editor = editor();
    for preset in [Preset::Drill, Preset::Pocket] {
        editor.tree_selection = BTreeSet::from([feature]);
        editor
            .begin_wood_preset(&session, &Default::default(), Channel::Ui, preset)
            .unwrap();
        assert_eq!(editor.model.form_targets(), Some(vec![left]));
        assert_eq!(editor.wood_preset(), Some(preset));
        assert!(editor.message.contains("Cabinet left side"));
        editor.cancel_model();
        assert_eq!(editor.wood_preset(), None);
        assert_eq!(session.document(), &before);
    }
}
#[test]
fn a_selected_panel_seeds_drill_pocket_and_joint_targets() {
    let (session, left, right) = panel_fixture::carcass();
    let mut editor = editor();
    // The side's Part or its panel feature both stand for the panel.
    let feature = model::panel_feature(session.document(), left).unwrap();
    editor.tree_selection = BTreeSet::from([feature]);
    open(&mut editor, &session, ModelingTool::Machining);
    assert_eq!(editor.model.form_targets(), Some(vec![left]));
    assert!(
        editor.message.contains("Cabinet left side"),
        "{}",
        editor.message
    );
    // Explicit multi-panel changes stay possible: the viewport/Inspector
    // add more targets to the same form.
    editor.cancel_model();

    // Two selected panels seed the joint's first pair.
    editor.tree_selection = BTreeSet::from([left, right]);
    open(&mut editor, &session, ModelingTool::Joint);
    let mut seeded = editor.model.form_targets().unwrap();
    seeded.sort();
    let mut expected = vec![left, right];
    expected.sort();
    assert_eq!(seeded, expected);
    editor.cancel_model();

    // System 32 takes both sides from one selected side.
    editor.tree_selection = BTreeSet::from([right]);
    open(&mut editor, &session, ModelingTool::ShelfPins);
    assert_eq!(editor.model.form_targets(), Some(vec![right, left]));
    editor.cancel_model();
}

#[test]
fn a_selection_that_is_not_a_panel_is_explained_and_seeds_nothing() {
    let (mut session, left, _) = panel_fixture::carcass();
    let knob = add_box(&mut session, "Printed control knob", 900.);
    let mut editor = editor();
    editor.tree_selection = BTreeSet::from([knob]);
    open(&mut editor, &session, ModelingTool::Machining);
    assert_eq!(editor.model.form_targets(), Some(vec![]));
    assert!(
        editor
            .message
            .contains("Printed control knob is not a cabinet panel"),
        "{}",
        editor.message
    );
    editor.cancel_model();
    // A panel and a non-panel: the panel is kept, the other is named.
    editor.tree_selection = BTreeSet::from([knob, left]);
    open(&mut editor, &session, ModelingTool::Machining);
    assert_eq!(editor.model.form_targets(), Some(vec![left]));
    assert!(editor.message.contains("not used"), "{}", editor.message);
    editor.cancel_model();
}

#[test]
fn a_selected_printed_part_seeds_the_print_aid_not_the_first_part() {
    let mut session = Session::new(Document::new("Front panel").unwrap()).unwrap();
    let sheet = add_box(&mut session, "Sheet part 1", 0.);
    let knob = add_box(&mut session, "Printed control knob", 100.);
    let mut editor = editor();
    for tool in [
        ModelingTool::FdmAid,
        ModelingTool::FdmBridge,
        ModelingTool::FdmSupport,
    ] {
        editor.tree_selection = BTreeSet::from([knob]);
        open(&mut editor, &session, tool);
        assert_eq!(editor.model.form_targets(), Some(vec![knob]));
        assert_eq!(editor.fdm_aid_kind(), tool.fdm_kind());
        assert!(
            editor.message.starts_with(tool.label()),
            "{}",
            editor.message
        );
        assert!(
            editor.message.contains("Printed control knob"),
            "{}",
            editor.message
        );
        editor.cancel_model();
        // Nothing selected with two Parts: no silent fallback to the first.
        editor.tree_selection.clear();
        open(&mut editor, &session, tool);
        assert_eq!(editor.model.form_targets(), Some(vec![]));
        assert!(
            editor.message.contains("Select the finished Part"),
            "{}",
            editor.message
        );
        editor.cancel_model();
        // Both selected: ambiguous, explained.
        editor.tree_selection = BTreeSet::from([sheet, knob]);
        open(&mut editor, &session, tool);
        assert_eq!(editor.model.form_targets(), Some(vec![]));
        assert!(
            editor.message.contains("2 Parts are selected"),
            "{}",
            editor.message
        );
        editor.cancel_model();
    }
}

/// CAD-392 in the real editor frame: with the joint tool active, clicking
/// panel A then panel B fills both inputs, the corner auto-previews, Escape
/// cancels without geometry and Enter applies one undoable joint.
#[test]
fn joint_tool_clicks_fill_a_then_b_preview_and_enter_or_escape_finish_atomically() {
    use std::time::{Duration, Instant};
    let (mut session, side, bottom) = panel_fixture::corner(false);
    let mut feed = SceneFeed::default();
    feed.set_backend(Box::new(NativeBackend::new()));
    let until = Instant::now() + Duration::from_secs(10);
    while !feed.is_current(session.document().id(), session.revision()) && Instant::now() < until {
        feed.tick(session.document(), session.revision());
        std::thread::sleep(Duration::from_millis(2));
    }
    assert!(feed.is_current(session.document().id(), session.revision()));
    let before = session.document().clone();
    let revision = session.revision();
    let ctx = egui::Context::default();
    let mut t = 0.;
    let mut frame =
        |editor: &mut SketchEditor, session: &mut Session, feed: &mut SceneFeed, events| {
            t += 0.1;
            ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(1360., 860.),
                    )),
                    time: Some(t),
                    events,
                    ..Default::default()
                },
                |ui| editor.show_project(ui, session, false, feed, &Default::default()),
            )
            .drop_without_applying_deltas();
            editor.poll_model(session, &Default::default());
        };
    let top = crate::camera::Camera {
        roll: 0.0,
        yaw: 0.,
        pitch: std::f32::consts::FRAC_PI_2,
        scale: 0.5,
        target: [300., 200., 0.],
        perspective: false,
    };
    let mut editor = editor();
    for confirm in [false, true] {
        editor
            .activate_model_tool(
                &session,
                &Default::default(),
                &feed,
                ModelingTool::Joint,
                Channel::Ui,
            )
            .unwrap();
        frame(&mut editor, &mut session, &mut feed, vec![]);
        editor.camera = top;
        frame(&mut editor, &mut session, &mut feed, vec![]);
        let rect = editor.viewport_rect.unwrap();
        for (point, expected) in [
            ([9f32, 200., 300.], vec![side]),
            ([300., 200., 18.], vec![side, bottom]),
        ] {
            let p = editor.camera.project(point, rect);
            let pos = egui::pos2(p[0], p[1]);
            for pressed in [true, false] {
                frame(
                    &mut editor,
                    &mut session,
                    &mut feed,
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
            assert_eq!(
                editor.model.form_targets(),
                Some(expected),
                "{}",
                editor.message
            );
        }
        // The staged preview becomes ready without touching the document.
        let until = Instant::now() + Duration::from_secs(20);
        while !editor.model_can_apply(&session, &Default::default()) && Instant::now() < until {
            frame(&mut editor, &mut session, &mut feed, vec![]);
            std::thread::sleep(Duration::from_millis(10));
        }
        assert!(
            editor.model_can_apply(&session, &Default::default()),
            "{}",
            editor.message
        );
        assert_eq!(session.document(), &before);
        if !confirm {
            // Escape runs the shared cancel: no partial geometry remains.
            editor.cancel();
            assert!(!editor.model.active());
            assert_eq!(session.document(), &before);
            assert_eq!(session.revision(), revision);
            continue;
        }
        assert!(editor.arm_model_confirmation(&session, &Default::default()));
        assert!(
            editor
                .finish_model_confirmation(&mut session, &Default::default())
                .unwrap()
        );
        assert_eq!(session.revision(), revision + 1);
        let joints: Vec<_> = model::panel_joints(session.document()).collect();
        assert_eq!(joints.len(), 1);
        let joint = model::read_panel_joint(session.document(), joints[0].id).unwrap();
        assert_eq!(
            joint.pairs,
            vec![model::PanelJointPair { a: side, b: bottom }]
        );
    }
}

/// CAD-393 in the real editor frame: a butted corner shows the corner
/// block and the bottom's missing 18 mm in the viewport and the exact fix
/// in the Inspector, and stays refused without changing the document.
#[test]
fn a_butted_joint_corner_shows_the_missing_material_and_the_fix() {
    use std::time::{Duration, Instant};
    let (mut session, side, bottom) = panel_fixture::corner(true);
    let mut feed = SceneFeed::default();
    feed.set_backend(Box::new(NativeBackend::new()));
    let until = Instant::now() + Duration::from_secs(10);
    while !feed.is_current(session.document().id(), session.revision()) && Instant::now() < until {
        feed.tick(session.document(), session.revision());
        std::thread::sleep(Duration::from_millis(2));
    }
    let before = session.document().clone();
    let mut editor = editor();
    editor.tree_selection = BTreeSet::from([side, bottom]);
    open(&mut editor, &session, ModelingTool::Joint);
    // A two-panel selection is unordered; read which role each panel got.
    let seeded = editor.model.form_targets().unwrap();
    assert_eq!(seeded.len(), 2);
    let (side_role, bottom_role) = if seeded[0] == side {
        ("A", "B")
    } else {
        ("B", "A")
    };
    let ctx = egui::Context::default();
    let wanted = [
        format!("Pair 1 · {side_role} · Side"),
        format!("Pair 1 · {bottom_role} · Bottom"),
        "Corner block both panels must fill".into(),
        format!("+18 mm length of {bottom_role}"),
        format!("lengthen Bottom ({bottom_role}) by 18 mm at its start (length 582 → 600 mm)"),
    ];
    let mut texts = String::new();
    let until = Instant::now() + Duration::from_secs(20);
    for i in 0.. {
        let out = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1360., 1400.),
                )),
                time: Some(i as f64 * 0.1),
                ..Default::default()
            },
            |ui| editor.show_project(ui, &mut session, false, &mut feed, &Default::default()),
        );
        editor.poll_model(&session, &Default::default());
        if i == 0 {
            editor.camera = crate::camera::Camera {
                roll: 0.0,
                yaw: 0.6,
                pitch: 0.5,
                scale: 0.5,
                target: [300., 200., 100.],
                perspective: false,
            };
        }
        texts.clear();
        for shape in &out.shapes {
            if let egui::Shape::Text(text) = &shape.shape {
                texts.push_str(text.galley.text());
                texts.push('\n');
            }
        }
        out.drop_without_applying_deltas();
        if i > 2 && wanted.iter().all(|w| texts.contains(w.as_str())) || Instant::now() > until {
            break;
        }
        std::thread::sleep(Duration::from_millis(5));
    }
    for wanted in wanted {
        assert!(texts.contains(&wanted), "missing {wanted:?} in\n{texts}");
    }
    assert!(!editor.model_can_apply(&session, &Default::default()));
    assert_eq!(session.document(), &before);
    editor.cancel();
}
