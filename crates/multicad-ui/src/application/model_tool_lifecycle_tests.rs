//! CAD-413: Model tool lifecycle through the full application. A plain
//! viewport right-click ends the active tool (a complete form applies as
//! Apply does, an incomplete one ends unchanged), and in Select restarts the
//! last tool with its settings; a right-drag only navigates.
use super::ui_harness::Harness;
use super::*;
use crate::modeling::ModelingTool;
use multicad_core::{
    constraint_sketch::Dimension,
    model::{self, PrimitiveParameters, PrimitiveRequest},
    references::{GeometryIntent, TopologyKind},
};
use std::collections::BTreeSet;

fn box_document() -> (Document, model::BodyInput) {
    let mut session = Session::new(Document::new("Model tool lifecycle").unwrap()).unwrap();
    let request = PrimitiveRequest::new(
        "Stock",
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
    let body = model::BodyInput {
        body: request.body,
        producer: request.feature,
    };
    (session.document().clone(), body)
}
fn select_part(h: &mut Harness, body: EntityId) {
    h.app
        .current_mut()
        .sketch
        .set_tree_selection_for_tests(BTreeSet::from([body]));
}
fn open_fillet(h: &mut Harness, body: EntityId) {
    select_part(h, body);
    h.app
        .invoke(
            Command::Modeling(ModelingTool::EdgeTreatment),
            &[],
            Channel::Ui,
        )
        .unwrap();
    h.frame(vec![]);
    assert!(fillet_active(h));
}
fn fillet_active(h: &Harness) -> bool {
    h.app
        .current()
        .sketch
        .model_for_tests()
        .edge_treatment_active()
}
/// Pick the first edge of the body's current producer and wait for the
/// native preview.
fn pick_edge(h: &mut Harness, producer: EntityId) {
    let p = h.app.current_mut();
    let scene = p.scene.visible().unwrap().clone();
    let topology = scene.geometry_snapshot(producer).unwrap();
    let edge = topology
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Edge)
        .unwrap();
    let capture = GeometryIntent::capture(&topology, edge.key).unwrap();
    p.sketch
        .model_for_tests_mut()
        .set_edge_treatment_edges(&p.session, &p.resources, &scene, vec![capture])
        .unwrap();
    h.until("the fillet preview", |app| {
        let p = app.current();
        p.sketch.model_can_apply(&p.session, &p.resources)
    });
}
/// An empty corner of the 3D viewport.
fn viewport_corner(h: &Harness) -> egui::Pos2 {
    let rect = h.app.current().sketch.viewport_rect.unwrap();
    rect.left_top() + egui::vec2(60., 60.)
}
fn right_click(h: &mut Harness) {
    let at = viewport_corner(h);
    h.frame(vec![egui::Event::PointerMoved(at)]);
    for pressed in [true, false] {
        h.frame(vec![
            egui::Event::PointerMoved(at),
            egui::Event::PointerButton {
                pos: at,
                button: egui::PointerButton::Secondary,
                pressed,
                modifiers: egui::Modifiers::NONE,
            },
        ]);
    }
}
fn wait_current(h: &mut Harness) {
    h.until("the current model", |app| {
        let p = app.current();
        p.scene
            .is_current(p.session.document().id(), p.session.revision())
    });
}
fn last_radius_mm(document: &Document, body: EntityId) -> f64 {
    let EntityData::Part { features, .. } = &document.get(body).unwrap().data else {
        panic!("Part expected")
    };
    let (saved, _) = model::read_edge_treatment(document, *features.last().unwrap()).unwrap();
    match saved.parameters.size {
        Dimension::Literal { value, .. } => value,
        other => panic!("literal size expected, got {other:?}"),
    }
}

fn size_text(h: &Harness) -> Option<String> {
    h.app
        .current()
        .sketch
        .model_for_tests()
        .edge_treatment_size_text()
}

/// AC2 for Model tools: a right-click on a complete Fillet applies it exactly
/// as Apply does and ends the tool.
#[test]
fn right_click_applies_a_complete_model_fillet_and_ends_the_tool() {
    let (document, body) = box_document();
    let mut h = Harness::new(document);
    open_fillet(&mut h, body.body);
    pick_edge(&mut h, body.producer);
    assert_eq!(h.revision(), 0);
    right_click(&mut h);
    h.until("the right-click to apply the fillet", |app| {
        app.current().session.revision() == 1
    });
    assert!(!h.modeling(), "the tool ended");
    assert!(
        h.app.message.contains("right-click repeats"),
        "{}",
        h.app.message
    );
    let document = h.app.current().session.document().clone();
    let EntityData::Part { features, .. } = &document.get(body.body).unwrap().data else {
        panic!("Part expected")
    };
    model::read_edge_treatment(&document, *features.last().unwrap()).unwrap();
}

/// AC3 for Model tools: in Select, a plain right-click restarts the last
/// tool with the 3 mm radius it was applied with.
#[test]
fn right_click_in_select_repeats_the_model_fillet_with_its_radius() {
    let (document, body) = box_document();
    let mut h = Harness::new(document);
    open_fillet(&mut h, body.body);
    assert_ne!(size_text(&h).as_deref(), Some("3"), "the default differs");
    pick_edge(&mut h, body.producer);
    // Enter in the Radius field applies, as Apply does (CAD-410).
    h.type_into_field_and_enter("Radius", "3");
    h.until("the 3 mm fillet", |app| {
        app.current().session.revision() == 1
    });
    assert!(!h.modeling());
    let document = h.app.current().session.document().clone();
    assert!((last_radius_mm(&document, body.body) - 3.).abs() < 1e-9);
    wait_current(&mut h);
    select_part(&mut h, body.body);
    right_click(&mut h);
    h.frame(vec![]);
    assert!(fillet_active(&h), "{}", h.app.message);
    assert!(h.app.message.contains("repeated"), "{}", h.app.message);
    assert_eq!(size_text(&h).as_deref(), Some("3"), "the radius is kept");
    assert_eq!(h.revision(), 1, "repeating commits nothing");
    // A fresh Fillet from the ribbon starts from the default again.
    h.key(egui::Key::Escape);
    open_fillet(&mut h, body.body);
    assert_ne!(size_text(&h).as_deref(), Some("3"));
}

/// AC2: an incomplete Model tool ends on right-click or Enter without
/// changing the document; Escape still cancels.
#[test]
fn right_click_or_enter_ends_an_incomplete_model_tool_unchanged() {
    let (document, body) = box_document();
    let mut h = Harness::new(document.clone());
    open_fillet(&mut h, body.body);
    right_click(&mut h);
    assert!(!h.modeling());
    assert_eq!(h.app.current().session.document(), &document);
    assert!(h.app.message.contains("ended"), "{}", h.app.message);
    open_fillet(&mut h, body.body);
    h.key(egui::Key::Enter);
    assert!(!h.modeling());
    assert_eq!(h.app.current().session.document(), &document);
    open_fillet(&mut h, body.body);
    h.key(egui::Key::Escape);
    assert!(!h.modeling());
    assert_eq!(h.revision(), 0);
}

/// AC3: a right-drag in the Model viewport navigates and does not repeat.
#[test]
fn a_right_drag_in_the_model_viewport_does_not_repeat_the_tool() {
    let (document, body) = box_document();
    let mut h = Harness::new(document);
    open_fillet(&mut h, body.body);
    h.key(egui::Key::Escape);
    assert!(!h.modeling());
    select_part(&mut h, body.body);
    let camera = h.app.current().sketch.camera;
    let start = viewport_corner(&h);
    let button = |pos, pressed| egui::Event::PointerButton {
        pos,
        button: egui::PointerButton::Secondary,
        pressed,
        modifiers: egui::Modifiers::NONE,
    };
    h.frame(vec![egui::Event::PointerMoved(start)]);
    h.frame(vec![button(start, true)]);
    for i in 1..=8 {
        h.frame(vec![egui::Event::PointerMoved(
            start + egui::vec2(6. * i as f32, 4. * i as f32),
        )]);
    }
    h.frame(vec![button(start + egui::vec2(48., 32.), false)]);
    h.frame(vec![]);
    assert!(!h.modeling(), "a right-drag does not repeat the tool");
    assert_ne!(h.app.current().sketch.camera, camera, "the drag navigated");
}
