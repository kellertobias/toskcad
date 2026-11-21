//! New-sketch entry, cancellation and face targeting (CAD-380, CAD-381).
use super::tests::input;
use super::*;
use crate::viewport::{Geometry, Instance, Scene, SceneStamp, Transform};
use multicad_core::document::{Definition, Entity, Placement};
use multicad_core::project_file::{self, ProjectFile};
use multicad_core::sketch::Datum;
use multicad_kernel::{MeshOptions, Solid};
use std::sync::Arc;

/// A 100 × 60 × 18 mm box placed rotated and translated, with its native
/// topology, so face hits resolve exactly like the application's scene.
pub(super) struct Panel {
    pub document: Document,
    pub producer: EntityId,
    pub instance: EntityId,
    pub pose: Transform,
    geometry: Arc<Geometry>,
}
impl Panel {
    pub fn new() -> Self {
        Self::placed(Placement {
            translation_mm: [140., 25., 40.],
            rotation_xyzw: [
                0.,
                0.,
                std::f64::consts::FRAC_1_SQRT_2,
                std::f64::consts::FRAC_1_SQRT_2,
            ],
        })
    }
    pub fn placed(placement: Placement) -> Self {
        let producer = Entity::new(
            "Panel source",
            EntityData::Feature {
                definition: Definition {
                    type_id: "test.box".into(),
                    schema_version: 1,
                    data: serde_json::json!({}),
                },
                inputs: vec![],
                suppressed: false,
            },
        );
        let part = Entity::new(
            "Panel",
            EntityData::Part {
                features: vec![producer.id],
                material: None,
                tags: Default::default(),
            },
        );
        let instance = Entity::new(
            "Placed panel",
            EntityData::Instance {
                component: part.id,
                placement,
                placement_source: None,
            },
        );
        let mut document = Document::new("Sketch entry").unwrap();
        for entity in [producer.clone(), part, instance.clone()] {
            document.insert(entity).unwrap();
        }
        let native = Solid::box_named([100., 60., 18.], producer.id, Default::default()).unwrap();
        let mesh = native.mesh(MeshOptions::default()).unwrap();
        let geometry = Arc::new(
            Geometry::new(
                producer.id,
                mesh.vertices.as_chunks::<3>().0.to_vec(),
                mesh.indices,
                mesh.face_indices,
                native.topology_snapshot(producer.id, 0).unwrap().elements,
            )
            .unwrap(),
        );
        Self {
            document,
            producer: producer.id,
            instance: instance.id,
            pose: Transform::from_placement(placement).unwrap(),
            geometry,
        }
    }
    pub fn scene(&self, revision: u64) -> Arc<Scene> {
        Arc::new(
            Scene::new(
                SceneStamp {
                    document: self.document.id(),
                    revision,
                    generation: 1,
                },
                vec![self.geometry.clone()],
                vec![Instance {
                    id: self.instance,
                    geometry: self.producer,
                    transform: self.pose,
                    color: [140, 165, 190, 255],
                }],
            )
            .unwrap(),
        )
    }
    /// An editor looking straight at the panel along `-normal` (a world
    /// direction), so a click at the panel-local `point` hits that face.
    pub fn editor(&self, revision: u64, normal: [f64; 3]) -> SketchEditor {
        let mut editor = SketchEditor {
            display_scene: Some(self.scene(revision)),
            ..Default::default()
        };
        editor.camera.target = self.pose.point([50., 30., 9.]).map(|v| v as f32);
        editor.camera.yaw = normal[0].atan2(normal[1]) as f32;
        editor.camera.pitch = normal[2].clamp(-1., 1.).asin() as f32;
        editor
    }
    pub fn screen(&self, editor: &SketchEditor, local: [f64; 3]) -> Pos2 {
        let rect = Rect::from_min_size(Pos2::ZERO, Vec2::new(1000., 700.));
        let p = editor
            .camera
            .project(self.pose.point(local).map(|v| v as f32), rect);
        Pos2::new(p[0], p[1])
    }
}
/// Hover then click once at `point`, rendering the canvas for each frame.
pub(super) fn click(
    editor: &mut SketchEditor,
    session: &mut Session,
    scene: Arc<Scene>,
    sketch_mode: bool,
    point: Pos2,
) {
    let ctx = egui::Context::default();
    for events in [
        vec![],
        vec![egui::Event::PointerMoved(point)],
        vec![egui::Event::PointerButton {
            pos: point,
            button: egui::PointerButton::Primary,
            pressed: true,
            modifiers: egui::Modifiers::NONE,
        }],
        vec![egui::Event::PointerButton {
            pos: point,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: egui::Modifiers::NONE,
        }],
    ] {
        ctx.run_ui(input(events), |ui| {
            editor.canvas(ui, session, sketch_mode, Some(scene.clone()), true)
        })
        .drop_without_applying_deltas();
    }
}
fn sketches(document: &Document) -> usize {
    document
        .entities()
        .filter(|e| matches!(e.data, EntityData::Sketch { .. }))
        .count()
}
fn reopened(document: &Document) -> Document {
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(document.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    project_file::read(
        &mut std::io::Cursor::new(bytes),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document
}

#[test]
fn cancelling_an_untouched_new_datum_sketch_withdraws_it_without_history() {
    let panel = Panel::new();
    let mut session = Session::new(panel.document.clone()).unwrap();
    // A genuine pre-existing unsaved edit must survive the cancellation.
    session
        .apply(Request {
            expected_revision: 0,
            actor: Actor::local(Channel::Ui),
            label: "Rename".into(),
            edits: vec![Mutation::RenameDocument("Edited".into())],
        })
        .unwrap();
    let before = session.document().clone();
    let mut editor = SketchEditor::default();
    editor.create(&mut session, Datum::Xy, Channel::Ui).unwrap();
    assert!(editor.editing);
    assert_eq!(sketches(session.document()), 1);
    editor.abort(&mut session, Channel::Ui).unwrap();
    assert!(!editor.editing && editor.selected.is_none());
    assert_eq!(session.document(), &before);
    assert_eq!(sketches(session.document()), 0);
    assert!(session.dirty(), "the earlier rename is still unsaved");
    assert!(!session.can_redo(), "no redoable phantom creation");
    session
        .undo(session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(session.document(), &panel.document);
    assert!(!session.can_undo(), "only the rename was ever undoable");
    session
        .redo(session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(reopened(session.document()), before);
    assert!(editor.message.contains("cancelled"), "{}", editor.message);
}

#[test]
fn cancelling_an_untouched_face_sketch_keeps_the_part_and_its_references() {
    let panel = Panel::new();
    let mut session = Session::new(panel.document.clone()).unwrap();
    session.mark_saved(0).unwrap();
    let top = panel.pose.vector([0., 0., 1.]);
    let mut editor = panel.editor(0, top);
    let point = panel.screen(&editor, [30., 20., 18.]);
    click(&mut editor, &mut session, panel.scene(0), false, point);
    assert!(editor.solid.selected_hit().is_some(), "face selected first");
    editor
        .begin_model_sketch(&mut session, Channel::Ui)
        .unwrap();
    assert!(editor.editing && session.revision() == 1);
    assert_eq!(sketches(session.document()), 1);
    editor.abort(&mut session, Channel::Ui).unwrap();
    assert_eq!(session.document(), &panel.document);
    assert!(!session.can_undo() && !session.can_redo());
    assert!(!session.dirty());
    let opened = reopened(session.document());
    assert_eq!(opened, panel.document);
    assert!(opened.get(panel.producer).is_some() && opened.get(panel.instance).is_some());
}

#[test]
fn cancelling_an_existing_or_edited_sketch_keeps_it_and_only_drops_pending_work() {
    let panel = Panel::new();
    let mut session = Session::new(panel.document.clone()).unwrap();
    let mut editor = SketchEditor::default();
    // An empty sketch that was finished once is an existing sketch.
    editor.create(&mut session, Datum::Xz, Channel::Ui).unwrap();
    let empty = editor.selected.unwrap();
    editor.finish().unwrap();
    let kept = session.document().clone();
    editor.open(session.document(), empty).unwrap();
    editor.abort(&mut session, Channel::Ui).unwrap();
    assert_eq!(session.document(), &kept);
    assert!(session.document().get(empty).is_some());
    // A new sketch with committed geometry keeps that geometry on cancel.
    editor.create(&mut session, Datum::Xy, Channel::Ui).unwrap();
    editor
        .add(
            &mut session,
            DrawTool::Rectangle,
            [0., 0.],
            [10., 5.],
            Channel::Ui,
        )
        .unwrap();
    editor.set_tool(DrawTool::Line).unwrap();
    let drawn = session.document().clone();
    let revision = session.revision();
    editor.abort(&mut session, Channel::Ui).unwrap();
    assert_eq!(session.document(), &drawn);
    assert_eq!(session.revision(), revision);
    assert_eq!(sketches(session.document()), 2);
    // Cancelling plane picking before any plane is chosen changes nothing.
    editor.begin();
    editor.abort(&mut session, Channel::Ui).unwrap();
    assert_eq!(session.document(), &drawn);
    assert_eq!(session.revision(), revision);
}
