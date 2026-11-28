//! CAD-418: Model's Common group has an Add Sketch button that runs the
//! global Sketch command; clicking it and then a datum plane in the real
//! viewport opens the new sketch in one Undo step (CAD-397 behaviour).
use super::ui_harness::Harness;
use super::*;

fn sketches(document: &Document) -> usize {
    document
        .entities()
        .filter(|e| matches!(e.data, EntityData::Sketch { .. }))
        .count()
}

#[test]
fn model_common_add_sketch_picks_a_plane_and_opens_the_sketch() {
    let mut h = Harness::new(Document::new("Add sketch").unwrap());
    let before = h.app.current().session.document().clone();
    assert_eq!(h.app.mode, Mode::Model);
    assert_eq!(h.app.effective_ribbon_tab(), 0, "Common");
    h.click_label("Add Sketch");
    h.until("plane picking", |app| app.current().sketch.picking);
    assert_eq!(h.app.mode.top_level(), Mode::Model);
    assert!(!h.app.current().sketch.editing);
    assert_eq!(h.app.current().session.document(), &before, "nothing yet");
    // Look straight down and click the XY datum plane in the viewport.
    let camera = &mut h.app.current_mut().sketch.camera;
    camera.target = [0.; 3];
    camera.yaw = 0.;
    camera.pitch = std::f32::consts::FRAC_PI_2;
    camera.roll = 0.;
    h.frame(vec![]);
    let rect = h
        .app
        .current()
        .sketch
        .viewport_rect
        .expect("viewport drawn");
    let p = h.app.current().sketch.camera.project([10., 10., 0.], rect);
    h.click_at(egui::pos2(p[0], p[1]));
    h.until("the new sketch", |app| app.current().sketch.editing);
    assert!(!h.app.current().sketch.picking);
    assert_eq!(h.app.mode, Mode::Sketch);
    assert_eq!(sketches(h.app.current().session.document()), 1);
    assert_eq!(h.revision(), 1);
    // One transaction: a single Undo restores the document before Add Sketch.
    let session = &mut h.app.current_mut().session;
    session
        .undo(session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(session.document(), &before, "one Undo step");
    assert!(!session.can_undo());
}

#[test]
fn model_common_add_sketch_and_the_global_action_share_one_command() {
    let registry = crate::registry::Registry::default();
    let common = &registry.workspace(Mode::Model).tabs[0];
    assert_eq!(common.title, "Common");
    let tool = registry.tool("new-sketch").unwrap();
    assert!(common.tools.contains(&tool.id));
    // The global right-side action issues the same command.
    let h = Harness::new(Document::new("Shared").unwrap());
    assert_eq!(tool.action, Command::NewSketch);
    assert!(h.app.enabled(Command::NewSketch));
}
