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
        editor.camera.yaw = normal[0].portable_atan2(normal[1]) as f32;
        editor.camera.pitch = normal[2].clamp(-1., 1.).portable_asin() as f32;
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

#[test]
fn hovering_a_face_while_picking_shows_the_ordinary_face_highlight() {
    let panel = Panel::new();
    let mut session = Session::new(panel.document.clone()).unwrap();
    let top = panel.pose.vector([0., 0., 1.]);
    let mut editor = panel.editor(0, top);
    editor.begin();
    let point = panel.screen(&editor, [30., 20., 18.]);
    let ctx = egui::Context::default();
    for events in [vec![], vec![egui::Event::PointerMoved(point)], vec![]] {
        ctx.run_ui(input(events), |ui| {
            editor.canvas(ui, &mut session, true, Some(panel.scene(0)), true)
        })
        .drop_without_applying_deltas();
    }
    assert!(
        matches!(
            editor.solid.reference_highlight(),
            Some(crate::viewport::Selection::Face { instance, .. }) if *instance == panel.instance
        ),
        "hovered face is highlighted"
    );
    assert!(editor.picking && !editor.editing);
    assert_eq!(session.document(), &panel.document);
    // Away from the part nothing is highlighted.
    let empty = Pos2::new(5., 5.);
    for events in [vec![egui::Event::PointerMoved(empty)], vec![]] {
        ctx.run_ui(input(events), |ui| {
            editor.canvas(ui, &mut session, true, Some(panel.scene(0)), true)
        })
        .drop_without_applying_deltas();
    }
    assert!(editor.solid.reference_highlight().is_none());
}

#[test]
fn datum_planes_still_need_use_plane_and_invalid_faces_change_nothing() {
    let panel = Panel::new();
    let mut session = Session::new(panel.document.clone()).unwrap();
    let mut editor = SketchEditor {
        display_scene: Some(panel.scene(0)),
        ..Default::default()
    };
    editor.camera.target = [0., 0., 0.];
    editor.camera.pitch = std::f32::consts::FRAC_PI_2;
    editor.camera.yaw = 0.;
    editor.begin();
    let rect = Rect::from_min_size(Pos2::ZERO, Vec2::new(1000., 700.));
    let p = editor.camera.project([10., 10., 0.], rect);
    click(
        &mut editor,
        &mut session,
        panel.scene(0),
        true,
        Pos2::new(p[0], p[1]),
    );
    assert!(editor.picking && !editor.editing);
    assert!(editor.message.contains("Use plane"), "{}", editor.message);
    assert_eq!(session.document(), &panel.document);
    assert_eq!(session.revision(), 0);
    // A face that cannot host a sketch reports why and leaves the source alone.
    let EntityData::Feature {
        definition, inputs, ..
    } = panel.document.get(panel.producer).unwrap().data.clone()
    else {
        unreachable!()
    };
    let mut suppressed = panel.document.get(panel.producer).unwrap().clone();
    suppressed.data = EntityData::Feature {
        definition,
        inputs,
        suppressed: true,
    };
    let document = panel
        .document
        .edited(&[Mutation::Replace(suppressed)])
        .unwrap();
    let mut session = Session::new(document.clone()).unwrap();
    let top = panel.pose.vector([0., 0., 1.]);
    let mut editor = panel.editor(0, top);
    editor.begin();
    let point = panel.screen(&editor, [30., 20., 18.]);
    click(&mut editor, &mut session, panel.scene(0), true, point);
    assert!(editor.picking && !editor.editing);
    assert!(
        editor.message.contains("active modeling feature"),
        "{}",
        editor.message
    );
    assert_eq!(session.document(), &document);
    assert_eq!(session.revision(), 0);
}

fn toward_eye(camera: &crate::camera::Camera) -> [f64; 3] {
    let (s, c) = (camera.yaw as f64).portable_sin_cos();
    let (sp, cp) = (camera.pitch as f64).portable_sin_cos();
    [cp * s, cp * c, sp]
}
fn near(a: [f64; 3], b: [f64; 3], tolerance: f64) -> bool {
    (0..3).all(|i| (a[i] - b[i]).abs() < tolerance)
}

#[test]
fn face_sketches_look_at_their_face_from_outside_on_a_moved_and_rotated_part() {
    // Top, bottom and two sides of a translated, quarter-turned box. The
    // bottom is the prism start cap, whose semantic frame normal points into
    // the part; new face sketches orient it outward at attach time.
    for (local_point, local_normal) in [
        ([30., 20., 18.], [0., 0., 1.]),
        ([30., 20., 0.], [0., 0., -1.]),
        ([100., 30., 9.], [1., 0., 0.]),
        ([50., 0., 9.], [0., -1., 0.]),
    ] {
        let panel = Panel::new();
        let mut session = Session::new(panel.document.clone()).unwrap();
        let normal = panel.pose.vector(local_normal);
        let mut editor = panel.editor(0, normal);
        let point = panel.screen(&editor, local_point);
        click(&mut editor, &mut session, panel.scene(0), false, point);
        let hit = editor.solid.selected_hit().cloned().expect("face hit");
        assert_eq!(hit.instance, panel.instance);
        editor
            .begin_model_sketch(&mut session, Channel::Ui)
            .unwrap();
        assert!(editor.editing, "{local_normal:?}: {}", editor.message);
        // The view looks along the face normal from the outward side.
        assert!(
            near(toward_eye(&editor.camera), normal, 1e-6),
            "{local_normal:?}: eye {:?} vs {normal:?}",
            toward_eye(&editor.camera)
        );
        let id = editor.selected.unwrap();
        let revision = session.revision();
        let (_, world, _) = editor
            .read_source(session.document(), id, revision)
            .map(|(s, frame, plane)| (s, editor.display_frame(id, frame), plane))
            .unwrap();
        // The authored plane is the selected face in world space.
        assert!(
            near(world.n, normal, 1e-9),
            "{local_normal:?}: {:?}",
            world.n
        );
        let face_point = panel.pose.point(local_point);
        let offset: f64 = (0..3)
            .map(|i| (world.origin[i] - face_point[i]) * normal[i])
            .sum();
        assert!(
            offset.abs() < 1e-9,
            "{local_normal:?}: plane offset {offset}"
        );
        // Reopening restores the same view of the same face.
        let entered = editor.camera;
        editor.finish().unwrap();
        editor.camera = Default::default();
        editor
            .open_at(session.document(), id, session.revision())
            .unwrap();
        assert!(
            near(toward_eye(&editor.camera), toward_eye(&entered), 1e-6),
            "{local_normal:?}: reopened view differs"
        );
        // Sketch coordinates land on that face, inside its outline. The
        // rectangle starts at the clicked point: where the face lies relative
        // to the frame origin depends on which side the frame is seen from.
        let at = world.coordinates(face_point);
        let far = [at[0] + 2., at[1] + 1.];
        editor
            .add(&mut session, DrawTool::Rectangle, at, far, Channel::Ui)
            .unwrap();
        // The quarter-turned box stays axis aligned, so its world bounds are
        // exact; the drawn corner sits on the face, inside the part outline.
        let corners: Vec<[f64; 3]> = (0..8)
            .map(|k| {
                panel.pose.point(std::array::from_fn(|i| {
                    if k & (1 << i) == 0 {
                        0.
                    } else {
                        [100., 60., 18.][i]
                    }
                }))
            })
            .collect();
        let corner = world.point(far);
        let on_face: f64 = (0..3)
            .map(|i| (corner[i] - face_point[i]) * normal[i])
            .sum();
        assert!(on_face.abs() < 1e-9, "{local_normal:?}: {on_face}");
        for axis in 0..3 {
            let low = corners
                .iter()
                .map(|c| c[axis])
                .fold(f64::INFINITY, f64::min);
            let high = corners
                .iter()
                .map(|c| c[axis])
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(
                corner[axis] > low - 1e-9 && corner[axis] < high + 1e-9,
                "{local_normal:?}: {corner:?}"
            );
        }
        // Entering, reopening and drawing never move the part.
        assert_eq!(
            session.document().get(panel.instance),
            panel.document.get(panel.instance)
        );
    }
}

#[test]
fn clicking_a_face_while_picking_enters_its_sketch_without_use_plane() {
    let panel = Panel::new();
    let mut session = Session::new(panel.document.clone()).unwrap();
    let top = panel.pose.vector([0., 0., 1.]);
    let mut editor = panel.editor(0, top);
    editor.begin();
    assert!(editor.picking && !editor.editing);
    let point = panel.screen(&editor, [30., 20., 18.]);
    // Hovering shows the ordinary face highlight before anything is created.
    let ctx = egui::Context::default();
    for events in [vec![], vec![egui::Event::PointerMoved(point)], vec![]] {
        ctx.run_ui(input(events), |ui| {
            editor.canvas(ui, &mut session, true, Some(panel.scene(0)), true)
        })
        .drop_without_applying_deltas();
    }
    assert!(
        matches!(
            editor.solid.reference_highlight(),
            Some(crate::viewport::Selection::Face { instance, .. }) if *instance == panel.instance
        ),
        "hovered face is highlighted"
    );
    assert_eq!(session.document(), &panel.document);
    click(&mut editor, &mut session, panel.scene(0), true, point);
    assert!(editor.editing && !editor.picking, "{}", editor.message);
    assert_eq!(session.revision(), 1);
    assert_eq!(sketches(session.document()), 1);
    let id = editor.selected.unwrap();
    // A repeated click on the same face never adds a second sketch.
    click(&mut editor, &mut session, panel.scene(0), true, point);
    click(&mut editor, &mut session, panel.scene(1), true, point);
    assert_eq!(sketches(session.document()), 1);
    assert_eq!(editor.selected, Some(id));
    // The face attachment is canonical and persists through save/reopen.
    let opened = reopened(session.document());
    let plane = opened
        .entities()
        .find_map(|e| match &e.data {
            EntityData::Plane {
                anchor: Some(anchor),
                ..
            } => Some(anchor.target),
            _ => None,
        })
        .expect("face attachment saved");
    assert_eq!(plane, panel.producer);
    // The cancellation contract still holds for a face-click sketch.
    editor.abort(&mut session, Channel::Ui).unwrap();
    assert_eq!(session.document(), &panel.document);
    assert!(!session.can_undo());
}
