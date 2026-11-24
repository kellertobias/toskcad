//! CAD-383 snap-time projection, against the real native compiler: the
//! supporting face's corners are snap targets, a placement on them captures
//! canonical projected references together with the drawing in one checked
//! transaction, and the result survives Undo/Redo, save/reopen and an
//! upstream edit.
use super::*;
use crate::{
    modeling::native_fixture,
    viewport::{Scene, SceneBackend, SceneReply, SceneStamp},
};
use multicad_core::{
    model::{self, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, length},
    references::{AnalyticGeometry, ReferenceReport, TopologyKey, TopologyKind},
    sketch::{CreationPlane, PlaneOptions},
    sketch_draw::PointInput,
    sketch_external,
};
use multicad_model::{Compiler, ModelError, ModelSnapshot, Runtime};
use std::{
    sync::Arc,
    time::{Duration, Instant},
};

struct NativeBackend {
    runtime: Runtime<Arc<Scene>>,
    stamp: Option<(u64, SceneStamp)>,
    fail: bool,
}
fn publish(
    document: &Document,
    snapshot: &ModelSnapshot,
    generation: u64,
) -> Result<Arc<Scene>, ModelError> {
    native_fixture::scene(document, snapshot, generation)
}
impl NativeBackend {
    fn new(fail: bool) -> Self {
        Self {
            runtime: Runtime::with_publisher(Compiler::default(), publish).unwrap(),
            stamp: None,
            fail,
        }
    }
}
impl SceneBackend for NativeBackend {
    fn request(&mut self, document: &Document, stamp: SceneStamp) -> Result<(), String> {
        let generation = self.runtime.request_evaluation(
            document,
            stamp.revision,
            stamp.generation,
            &Default::default(),
        )?;
        self.stamp = Some((generation, stamp));
        Ok(())
    }
    fn poll(&mut self) -> Option<SceneReply> {
        let (generation, stamp) = self.stamp?;
        let reply = self.runtime.poll()?;
        if reply.generation != generation {
            return None;
        }
        self.stamp = None;
        Some(SceneReply {
            stamp,
            result: if self.fail {
                Err("Native evaluation failed".into())
            } else {
                reply.result.map_err(|e| e.to_string())
            },
        })
    }
    fn cancel(&mut self) {
        self.stamp = None;
    }
}

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
/// A 30 × 20 × 10 mm box, quarter-turned about Z and moved to (40, 10, 5).
fn turned() -> PrimitiveFrame {
    PrimitiveFrame {
        origin: [length(40.), length(10.), length(5.)],
        u: [0., 1., 0.],
        v: [-1., 0., 0.],
    }
}
fn box_request(size: [f64; 3]) -> PrimitiveParameters {
    PrimitiveParameters::Box {
        size: size.map(length),
        frame: turned(),
    }
}
fn native(session: &Session) -> (Arc<Scene>, ReferenceReport) {
    let compiled = Compiler::default()
        .compile(session.document(), session.revision())
        .unwrap();
    let scene = native_fixture::scene(session.document(), compiled.snapshot(), 1).unwrap();
    (scene, compiled.snapshot().references.clone())
}
/// The box's bottom corners in producer coordinates: z = 5 and the box
/// spans 30 mm along +Y and 20 mm along -X from (40, 10).
fn bottom_corners(length_mm: f64) -> Vec<[f64; 3]> {
    vec![
        [40., 10., 5.],
        [40., 10. + length_mm, 5.],
        [20., 10., 5.],
        [20., 10. + length_mm, 5.],
    ]
}
fn near(a: [f64; 3], b: [f64; 3]) -> bool {
    (0..3).all(|i| (a[i] - b[i]).abs() < 1e-9)
}

struct Fixture {
    session: Session,
    stock: PrimitiveRequest,
    sketch: EntityId,
    editor: SketchEditor,
}
fn fixture(fail: bool) -> Fixture {
    let mut session = Session::new(Document::new("Snap-time projection").unwrap()).unwrap();
    let stock = PrimitiveRequest::new("Box", box_request([30., 20., 10.]));
    let created = model::prepare_primitive(&session, 0, actor(), &stock).unwrap();
    session.commit(created.prepared).unwrap();
    let (scene, _) = native(&session);
    let snapshot = scene.geometry_snapshot(stock.feature).unwrap();
    let bottom = snapshot
        .elements
        .iter()
        .find(|e| matches!(e.geometry, Some(AnalyticGeometry::Plane { frame }) if frame.normal()[2] < -0.99))
        .unwrap()
        .key;
    let created = sketch::prepare_new(
        &session,
        session.revision(),
        actor(),
        "Bottom",
        CreationPlane::Face {
            snapshot: &snapshot,
            key: bottom,
            anchor: None,
        },
        PlaneOptions::default(),
    )
    .unwrap();
    let sketch = created.sketch;
    session.commit(created.prepared).unwrap();
    let editor = editor(&session, sketch, fail);
    Fixture {
        session,
        stock,
        sketch,
        editor,
    }
}
fn editor(session: &Session, sketch: EntityId, fail: bool) -> SketchEditor {
    let (scene, report) = native(session);
    let mut editor = SketchEditor {
        display_scene: Some(scene),
        evaluated: Some(report.clone()),
        ..Default::default()
    };
    editor
        .open_at(session.document(), sketch, session.revision())
        .unwrap();
    editor.evaluated = Some(report);
    editor.set_external_backend(Box::new(NativeBackend::new(fail)));
    editor.tool = Some(DrawTool::Rectangle);
    editor
}
/// Polls the native validation until the projection finishes.
fn settle(editor: &mut SketchEditor, session: &mut Session) {
    let deadline = Instant::now() + Duration::from_secs(20);
    while !editor.finish_projection(session) {
        assert!(Instant::now() < deadline, "native check never finished");
        editor.poll_model(session, &Default::default());
        std::thread::sleep(Duration::from_millis(2));
    }
}
/// Solved positions of each external reference's output point, in producer
/// coordinates, from a fresh native evaluation.
fn projected(session: &Session, sketch: EntityId) -> Vec<[f64; 3]> {
    let (_, report) = native(session);
    assert!(report.for_model(session.document()).unwrap().all_resolved());
    let evaluated = multicad_solver::evaluate_document_sketch_with_parameters(
        session.document(),
        session.revision(),
        sketch,
        Some(&report),
        session.parameters(),
    )
    .unwrap();
    let solved = evaluated.solved().expect("solved sketch");
    let (_, frame, _) = Source::read_for_evaluation(
        session.document(),
        sketch,
        session.revision(),
        Some(&report),
    )
    .unwrap();
    sketch_external::bindings(session.document(), sketch)
        .unwrap()
        .iter()
        .flat_map(|binding| binding.definition.output.points())
        .map(|id| {
            let p = solved.points.iter().find(|p| p.id == id).unwrap();
            frame.point(p.position_mm)
        })
        .collect()
}
fn corners(f: &Fixture) -> Vec<(TopologyKey, [f64; 2])> {
    f.editor.support_corners(&f.session, f.sketch)
}

#[test]
fn rectangle_between_face_corners_projects_both_and_commits_once() {
    let mut f = fixture(false);
    let (_, frame, _) = f
        .editor
        .read_source(f.session.document(), f.sketch, f.session.revision())
        .unwrap();
    // The four bottom corners of the turned, moved box, exactly.
    let offered = corners(&f);
    assert_eq!(offered.len(), 4);
    for corner in bottom_corners(30.) {
        assert!(
            offered.iter().any(|(_, at)| near(frame.point(*at), corner)),
            "{corner:?} not offered: {offered:?}"
        );
    }
    // Hovering a corner shows a face-corner target with the vertex marker.
    let geometry = snap::SnapGeometry::default().with_support(offered.clone());
    let to_screen = |p: [f64; 2]| Some([p[0] as f32 * 10., p[1] as f32 * 10.]);
    let at = offered[0].1;
    let outcome = snap::SnapAssist::default().update(
        &geometry,
        &snap::Query {
            pointer_mm: [at[0] + 0.2, at[1] - 0.1],
            pointer_px: [(at[0] + 0.2) as f32 * 10., (at[1] - 0.1) as f32 * 10.],
            to_screen: &to_screen,
            from: None,
        },
        0.,
    );
    let target = outcome.target.unwrap();
    assert_eq!(target.key, snap::TargetKey::Support(offered[0].0));
    assert_eq!(target.kind, snap::TargetKind::FaceCorner);
    assert_eq!(target.kind.marker(), snap::Marker::Vertex);
    assert_eq!(target.position, at);
    assert_eq!(target.anchor, None);

    // Draw a rectangle between two diagonal corners.
    let start = offered[0];
    let end = *offered
        .iter()
        .max_by(|a, b| {
            let d = |c: &(TopologyKey, [f64; 2])| {
                (c.1[0] - start.1[0]).powi(2) + (c.1[1] - start.1[1]).powi(2)
            };
            d(a).total_cmp(&d(b))
        })
        .unwrap();
    let before = f.session.revision();
    let original = f.session.document().clone();
    f.editor
        .draw_on_corners(
            &f.session,
            DrawTool::Rectangle,
            (PointInput::Position(start.1), Some(start.0)),
            (PointInput::Position(end.1), Some(end.0)),
            0.,
            Channel::Ui,
        )
        .unwrap();
    // Nothing is committed before the native check passes.
    assert!(f.editor.projection_pending() && f.editor.has_pending_edits());
    assert_eq!(f.session.revision(), before);
    assert_eq!(*f.session.document(), original);
    settle(&mut f.editor, &mut f.session);
    assert_eq!(f.session.revision(), before + 1, "{}", f.editor.message);
    let bindings = sketch_external::bindings(f.session.document(), f.sketch).unwrap();
    assert_eq!(bindings.len(), 2);
    let (_, report) = native(&f.session);
    let (source, _, _) = Source::read_for_evaluation(
        f.session.document(),
        f.sketch,
        f.session.revision(),
        Some(&report),
    )
    .unwrap();
    assert_eq!(source.sketch.curves.len(), 4);
    // The rectangle reuses the projected corner points.
    for binding in &bindings {
        let point = binding.definition.output.points()[0];
        assert!(source.sketch.curves.iter().any(|curve| matches!(
            curve.kind,
            multicad_core::constraint_sketch::CurveKind::Line { start, end }
                if start == point || end == point
        )));
    }
    let projected_now = projected(&f.session, f.sketch);
    for (corner, expected) in [(start, true), (end, true)] {
        assert_eq!(
            projected_now
                .iter()
                .any(|p| near(*p, frame.point(corner.1))),
            expected
        );
    }

    // One Undo step removes the references and the rectangle; Redo restores.
    f.session.undo(f.session.revision(), actor()).unwrap();
    assert_eq!(*f.session.document(), original);
    f.session.redo(f.session.revision(), actor()).unwrap();
    assert_eq!(
        sketch_external::bindings(f.session.document(), f.sketch)
            .unwrap()
            .len(),
        2
    );

    // Save and reopen keep the bindings and their exact resolution.
    let saved = serde_json::to_string(f.session.document()).unwrap();
    let reopened: Document = serde_json::from_str(&saved).unwrap();
    assert_eq!(reopened, *f.session.document());
    let reopened = Session::new(reopened).unwrap();
    assert_eq!(projected(&reopened, f.sketch).len(), 2);

    // The references follow an upstream edit: the box grows to 36 mm.
    let edited = model::prepare_primitive_edit(
        &f.session,
        f.session.revision(),
        actor(),
        f.stock.feature,
        &box_request([36., 20., 10.]),
    )
    .unwrap();
    f.session.commit(edited.prepared).unwrap();
    let moved = projected(&f.session, f.sketch);
    let grown = bottom_corners(36.);
    for p in &moved {
        assert!(grown.iter().any(|c| near(*p, *c)), "{p:?} left the face");
    }
    assert!(moved.iter().any(|p| (p[1] - 46.).abs() < 1e-9));
}

#[test]
fn a_placed_occurrence_projects_the_same_local_corners() {
    let mut f = fixture(false);
    let local = corners(&f);
    let mut placed = model::ComponentPlacementRequest::new(f.stock.body, Default::default());
    placed.parameters.translation = [length(200.), length(-50.), length(15.)];
    f.session
        .commit(
            model::prepare_place_component(&f.session, f.session.revision(), actor(), &placed)
                .unwrap()
                .prepared,
        )
        .unwrap();
    f.editor = editor(&f.session, f.sketch, false);
    let (scene, _) = native(&f.session);
    assert!(
        scene
            .instances()
            .iter()
            .any(|i| i.geometry == f.stock.feature && i.transform != Default::default()),
        "the occurrence is placed"
    );
    // The occurrence pose moves the display only; sketch and producer share
    // local coordinates, so the offered corners are unchanged.
    let offered = corners(&f);
    assert_eq!(offered.len(), 4);
    for (key, at) in &local {
        assert!(offered.iter().any(|(k, a)| k == key && a == at));
    }
    let before = f.session.revision();
    f.editor
        .draw_on_corners(
            &f.session,
            DrawTool::Line,
            (PointInput::Position(offered[0].1), Some(offered[0].0)),
            (PointInput::Position(offered[1].1), None),
            0.,
            Channel::Ui,
        )
        .unwrap();
    settle(&mut f.editor, &mut f.session);
    assert_eq!(f.session.revision(), before + 1, "{}", f.editor.message);
    assert_eq!(projected(&f.session, f.sketch).len(), 1);
}

#[test]
fn failed_native_check_unsupported_targets_and_stale_scenes_change_nothing() {
    // A failing native check reports and leaves the document unchanged.
    let mut f = fixture(true);
    let offered = corners(&f);
    let before = f.session.document().clone();
    f.editor
        .draw_on_corners(
            &f.session,
            DrawTool::Line,
            (PointInput::Position(offered[0].1), Some(offered[0].0)),
            (PointInput::Position(offered[1].1), Some(offered[1].0)),
            0.,
            Channel::Ui,
        )
        .unwrap();
    settle(&mut f.editor, &mut f.session);
    assert!(
        f.editor.message.contains("nothing was placed"),
        "{}",
        f.editor.message
    );
    assert_eq!(*f.session.document(), before);
    assert!(!f.editor.projection_pending());

    // A face is not a corner: nothing is guessed or prepared.
    let mut f = fixture(false);
    let (scene, _) = native(&f.session);
    let face = scene
        .geometry_snapshot(f.stock.feature)
        .unwrap()
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face)
        .unwrap()
        .key;
    let before = f.session.document().clone();
    assert!(
        f.editor
            .draw_on_corners(
                &f.session,
                DrawTool::Line,
                (PointInput::Position([0., 0.]), Some(face)),
                (PointInput::Position([1., 1.]), None),
                0.,
                Channel::Ui,
            )
            .is_err()
    );
    assert!(!f.editor.projection_pending());
    assert_eq!(*f.session.document(), before);

    // Tools without corner placement and stale scenes offer no corners.
    f.editor.tool = Some(DrawTool::Fillet);
    assert!(corners(&f).is_empty());
    f.editor.tool = Some(DrawTool::Rectangle);
    f.session
        .commit(
            f.session
                .prepare(Request {
                    expected_revision: f.session.revision(),
                    actor: actor(),
                    label: "Rename".into(),
                    edits: vec![Mutation::RenameDocument("Renamed".into())],
                })
                .unwrap(),
        )
        .unwrap();
    assert!(corners(&f).is_empty());
}

/// Runs one canvas frame after the sketch solve is current; returns the
/// pointer shape the canvas requested.
fn frame(
    ctx: &egui::Context,
    editor: &mut SketchEditor,
    session: &mut Session,
    events: Vec<egui::Event>,
    time: f64,
) -> egui::CursorIcon {
    if let Some(report) = editor
        .evaluated
        .clone()
        .filter(|r| r.revision == session.revision())
    {
        editor.sketch_evaluation.sync_creation_report(
            session.document(),
            session.revision(),
            report,
        );
    }
    let deadline = Instant::now() + Duration::from_secs(20);
    while editor.sketch_evaluation.pending() {
        editor.sketch_evaluation.poll();
        assert!(Instant::now() < deadline, "sketch solve timed out");
        std::thread::sleep(Duration::from_millis(1));
    }
    let raw = egui::RawInput {
        screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1000., 700.))),
        events,
        time: Some(time),
        ..Default::default()
    };
    let output = ctx.run_ui(raw, |ui| editor.canvas(ui, session, true, None, false));
    let cursor = output.platform_output.cursor_icon;
    output.drop_without_applying_deltas();
    cursor
}

#[test]
fn hovering_and_clicking_two_face_corners_on_the_canvas_draws_onto_projected_corners() {
    let mut f = fixture(false);
    let ctx = egui::Context::default();
    f.editor.tool = None;
    f.editor.set_tool(DrawTool::Rectangle).unwrap();
    let offered = corners(&f);
    let (_, local, _) = f
        .editor
        .read_source(f.session.document(), f.sketch, f.session.revision())
        .unwrap();
    let display = f.editor.display_frame(f.sketch, local);
    f.editor.camera.target = display.point([0., 0.]).map(|v| v as f32);
    f.editor.camera.scale = 8.;
    frame(&ctx, &mut f.editor, &mut f.session, vec![], 0.);
    let screen = |editor: &SketchEditor, p: [f64; 2]| {
        let s = editor.camera.project(
            display.point(p).map(|v| v as f32),
            editor.viewport_rect.unwrap(),
        );
        Pos2::new(s[0], s[1])
    };
    let start = offered[0];
    let end = *offered
        .iter()
        .find(|c| (c.1[0] - start.1[0]).abs() > 1. && (c.1[1] - start.1[1]).abs() > 1.)
        .unwrap();
    assert!(screen(&f.editor, start.1).distance(screen(&f.editor, end.1)) > 40.);
    // Hovering near a corner shows the face-corner target and the point cursor.
    let hover_at = screen(&f.editor, start.1) + Vec2::new(3., -2.);
    let cursor = frame(
        &ctx,
        &mut f.editor,
        &mut f.session,
        vec![egui::Event::PointerMoved(hover_at)],
        0.02,
    );
    let shown = f.editor.snap_assist.shown().cloned().unwrap();
    assert_eq!(shown.kind, snap::TargetKind::FaceCorner);
    assert_eq!(shown.position, start.1);
    assert_eq!(cursor, egui::CursorIcon::Crosshair);
    let before = f.session.revision();
    let mut time = 0.1;
    for at in [hover_at, screen(&f.editor, end.1) + Vec2::new(-2., 2.)] {
        for pressed in [true, false] {
            frame(
                &ctx,
                &mut f.editor,
                &mut f.session,
                vec![
                    egui::Event::PointerMoved(at),
                    egui::Event::PointerButton {
                        pos: at,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::NONE,
                    },
                ],
                time,
            );
            time += 0.01;
        }
    }
    assert!(f.editor.projection_pending(), "{}", f.editor.message);
    assert_eq!(f.session.revision(), before);
    // The canvas commits once the native check passes.
    let deadline = Instant::now() + Duration::from_secs(20);
    while f.session.revision() == before {
        assert!(Instant::now() < deadline, "{}", f.editor.message);
        f.editor.poll_model(&f.session, &Default::default());
        frame(&ctx, &mut f.editor, &mut f.session, vec![], time);
        std::thread::sleep(Duration::from_millis(2));
    }
    assert_eq!(f.session.revision(), before + 1, "{}", f.editor.message);
    let projected = projected(&f.session, f.sketch);
    assert_eq!(projected.len(), 2);
    for corner in [start, end] {
        assert!(projected.iter().any(|p| near(*p, local.point(corner.1))));
    }
}
