//! Add Sketch creates and opens the sketch on the first valid viewer pick
//! (CAD-397): datum planes and planar faces, one Undo step, clear feedback for
//! invalid picks, and correct frames on a placed, rotated occurrence.
use super::sketch_entry_tests::{Panel, click};
use super::*;
use crate::viewport::{Geometry, Instance, Scene, SceneStamp, Transform};
use multicad_core::document::{Definition, Entity, Placement};
use multicad_kernel::{MeshOptions, PlaneFrame, Solid};
use std::sync::Arc;

fn sketches(document: &Document) -> usize {
    document
        .entities()
        .filter(|e| matches!(e.data, EntityData::Sketch { .. }))
        .count()
}
fn undo(session: &mut Session) {
    session
        .undo(session.revision(), Actor::local(Channel::Ui))
        .unwrap();
}
fn top_down(editor: &mut SketchEditor) -> Pos2 {
    editor.camera.target = [0., 0., 0.];
    editor.camera.pitch = std::f32::consts::FRAC_PI_2;
    editor.camera.yaw = 0.;
    let rect = Rect::from_min_size(Pos2::ZERO, Vec2::new(1000., 700.));
    let p = editor.camera.project([10., 10., 0.], rect);
    Pos2::new(p[0], p[1])
}

#[test]
fn picking_a_datum_plane_in_the_viewer_creates_and_opens_the_sketch_in_one_undo_step() {
    let panel = Panel::new();
    let mut session = Session::new(panel.document.clone()).unwrap();
    let mut editor = SketchEditor {
        display_scene: Some(panel.scene(0)),
        ..Default::default()
    };
    editor.begin();
    assert!(editor.picking && !editor.editing);
    let point = top_down(&mut editor);
    click(&mut editor, &mut session, panel.scene(0), true, point);
    // No Use plane / Apply step: the pick alone creates and enters the sketch.
    assert!(editor.editing && !editor.picking, "{}", editor.message);
    assert_eq!(session.revision(), 1);
    assert_eq!(sketches(session.document()), 1);
    let id = editor.selected.expect("new sketch is open");
    let (_, frame, _) = editor
        .read_source(session.document(), id, session.revision())
        .unwrap();
    assert_eq!(frame.n, [0., 0., 1.], "XY datum frame");
    assert_eq!(frame.origin, [0., 0., 0.]);
    // One transaction: a single Undo restores the source before Add Sketch.
    undo(&mut session);
    assert_eq!(session.document(), &panel.document);
    assert!(!session.can_undo());
    assert!(session.can_redo());
}

#[test]
fn picking_a_face_creates_one_transaction_that_one_undo_removes() {
    let panel = Panel::new();
    let mut session = Session::new(panel.document.clone()).unwrap();
    let top = panel.pose.vector([0., 0., 1.]);
    let mut editor = panel.editor(0, top);
    editor.begin();
    let point = panel.screen(&editor, [30., 20., 18.]);
    click(&mut editor, &mut session, panel.scene(0), true, point);
    assert!(editor.editing && !editor.picking, "{}", editor.message);
    assert_eq!(session.revision(), 1);
    assert_eq!(session.activity().count(), 1, "one recorded transaction");
    undo(&mut session);
    assert_eq!(session.document(), &panel.document);
    assert!(!session.can_undo());
}

#[test]
fn empty_space_and_escape_create_nothing_and_say_why() {
    let panel = Panel::new();
    let mut session = Session::new(panel.document.clone()).unwrap();
    let top = panel.pose.vector([0., 0., 1.]);
    let mut editor = panel.editor(0, top);
    editor.begin();
    // Far away from the part and every datum plane's pick area.
    click(
        &mut editor,
        &mut session,
        panel.scene(0),
        true,
        Pos2::new(3., 3.),
    );
    assert!(editor.picking && !editor.editing);
    assert!(
        editor.message.contains("nothing was created"),
        "{}",
        editor.message
    );
    assert!(editor.message.contains("Escape"), "{}", editor.message);
    assert_eq!(session.document(), &panel.document);
    assert_eq!(session.revision(), 0);
    // Escape (AbortSketch) while picking leaves the document untouched.
    editor.abort(&mut session, Channel::Ui).unwrap();
    assert!(!editor.picking && !editor.editing);
    assert_eq!(session.document(), &panel.document);
    assert!(!session.can_undo());
}

#[test]
fn a_curved_face_is_refused_with_feedback_and_creates_nothing() {
    let producer = Entity::new(
        "Round source",
        EntityData::Feature {
            definition: Definition {
                type_id: "test.cylinder".into(),
                schema_version: 1,
                data: serde_json::json!({}),
            },
            inputs: vec![],
            suppressed: false,
        },
    );
    let part = Entity::new(
        "Round",
        EntityData::Part {
            features: vec![producer.id],
            material: None,
            tags: Default::default(),
        },
    );
    let placement = Placement::default();
    let instance = Entity::new(
        "Placed round",
        EntityData::Instance {
            component: part.id,
            placement,
            placement_source: None,
        },
    );
    let mut document = Document::new("Curved pick").unwrap();
    for entity in [producer.clone(), part, instance.clone()] {
        document.insert(entity).unwrap();
    }
    let native = Solid::cylinder(PlaneFrame::default(), 20., 30., Default::default()).unwrap();
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
    let scene = Arc::new(
        Scene::new(
            SceneStamp {
                document: document.id(),
                revision: 0,
                generation: 1,
            },
            vec![geometry],
            vec![Instance {
                id: instance.id,
                geometry: producer.id,
                transform: Transform::from_placement(placement).unwrap(),
                color: [140, 165, 190, 255],
            }],
        )
        .unwrap(),
    );
    let mut session = Session::new(document.clone()).unwrap();
    let mut editor = SketchEditor {
        display_scene: Some(scene.clone()),
        ..Default::default()
    };
    // Look horizontally at the cylinder's side wall (+Y side, mid height).
    editor.camera.target = [0., 0., 15.];
    editor.camera.pitch = 0.;
    editor.camera.yaw = 0.;
    editor.begin();
    let rect = Rect::from_min_size(Pos2::ZERO, Vec2::new(1000., 700.));
    let p = editor.camera.project([0., 20., 15.], rect);
    click(
        &mut editor,
        &mut session,
        scene.clone(),
        true,
        Pos2::new(p[0], p[1]),
    );
    assert!(editor.picking && !editor.editing, "{}", editor.message);
    assert!(
        editor.message.contains("planar") || editor.message.contains("choose another face"),
        "{}",
        editor.message
    );
    assert_eq!(session.document(), &document);
    assert_eq!(session.revision(), 0);
}

#[test]
fn picking_a_face_on_a_tilted_placed_occurrence_authors_that_face_frame() {
    // A part instance turned 30° about a skew axis and moved away from the
    // origin: the new sketch plane must be the picked face in world space.
    let axis: [f64; 3] = [1., 2., 2.];
    let len = (axis.iter().map(|v| v * v).sum::<f64>()).sqrt();
    let (s, c) = (std::f64::consts::PI / 12.).portable_sin_cos();
    let placement = Placement {
        translation_mm: [-60., 35., 80.],
        rotation_xyzw: [axis[0] / len * s, axis[1] / len * s, axis[2] / len * s, c],
    };
    let panel = Panel::placed(placement);
    let mut session = Session::new(panel.document.clone()).unwrap();
    let normal = panel.pose.vector([0., 0., 1.]);
    let mut editor = panel.editor(0, normal);
    editor.begin();
    let local = [30., 20., 18.];
    let point = panel.screen(&editor, local);
    click(&mut editor, &mut session, panel.scene(0), true, point);
    assert!(editor.editing && !editor.picking, "{}", editor.message);
    assert_eq!(session.revision(), 1);
    let id = editor.selected.unwrap();
    let (_, frame, _) = editor
        .read_source(session.document(), id, session.revision())
        .unwrap();
    let world = editor.display_frame(id, frame);
    for (n, expected) in world.n.iter().zip(normal) {
        assert!((n - expected).abs() < 1e-9, "{:?}", world.n);
    }
    let face_point = panel.pose.point(local);
    let offset: f64 = (0..3)
        .map(|i| (world.origin[i] - face_point[i]) * normal[i])
        .sum();
    assert!(offset.abs() < 1e-9, "plane offset {offset}");
    // The part placement itself is untouched by the pick.
    assert_eq!(
        session.document().get(panel.instance),
        panel.document.get(panel.instance)
    );
    undo(&mut session);
    assert_eq!(session.document(), &panel.document);
}
