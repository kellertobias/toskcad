//! Snap-time projection (CAD-383) through the real editor and the
//! production native scene backend: face corners and edges, every draw tool
//! that places points, entered sizes, merged captures, explicit refusals and
//! cancellation. Every wait is on a readiness condition (scene current,
//! target shown, native validation finished), never on a frame count.
//! Ported from the parallel CAD-383 implementation onto the merged design.
use super::*;
use crate::{modeling::native_fixture::NativeBackend, viewport::SceneFeed};
use multicad_core::{
    constraint_sketch::{ConstraintKind, CurveKind, Dimension, source::Source},
    model::{
        self, ComponentPlacementRequest, NumericTransform, PrimitiveFrame, PrimitiveParameters,
        PrimitiveRequest, length,
    },
    project_file::{self, ProjectFile},
    references::{AnalyticGeometry, Resolution, TopologyKey, TopologyKind},
    sketch::{CreationPlane, PlaneOptions},
    sketch_draw::PointInput,
    sketch_external,
};
use multicad_model::Compiler;
use snap::TargetKind;
use std::collections::{BTreeMap, BTreeSet};
use std::time::{Duration, Instant};

const SIZE: [f64; 3] = [30., 20., 10.];

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    a.map(|v| v * s)
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// A box whose own frame is rotated and translated, in the Part's coordinates.
struct Stock {
    request: PrimitiveRequest,
    origin: [f64; 3],
    u: [f64; 3],
    v: [f64; 3],
}
impl Stock {
    fn new(session: &mut Session) -> Self {
        let (s, c) = (0.5f64, 0.75f64.sqrt());
        let origin = [12., -7., 4.];
        let (u, v) = ([c, s, 0.], [-s, c, 0.]);
        let request = PrimitiveRequest::new(
            "Rotated stock",
            PrimitiveParameters::Box {
                size: SIZE.map(length),
                frame: PrimitiveFrame {
                    origin: origin.map(length),
                    u,
                    v,
                },
            },
        );
        session
            .commit(
                model::prepare_primitive(session, session.revision(), actor(), &request)
                    .unwrap()
                    .prepared,
            )
            .unwrap();
        Self {
            request,
            origin,
            u,
            v,
        }
    }
    /// A top-face corner in Part coordinates.
    fn corner(&self, x: f64, y: f64) -> [f64; 3] {
        let n = cross(self.u, self.v);
        add(
            add(self.origin, scale(self.u, x)),
            add(scale(self.v, y), scale(n, SIZE[2])),
        )
    }
}
/// A sketch attached to the stock's top face.
fn top_face_sketch(session: &mut Session, stock: &Stock) -> (EntityId, Frame) {
    let native = Compiler::default()
        .compile(session.document(), session.revision())
        .unwrap();
    let snapshot = &native.snapshot().shapes[&stock.request.feature].topology;
    let n = cross(stock.u, stock.v);
    let element = snapshot
        .elements
        .iter()
        .find(|e| {
            e.key.kind == TopologyKind::Face
                && matches!(e.geometry, Some(AnalyticGeometry::Plane { frame }) if {
                    let normal = cross(frame.u, frame.v);
                    let along: f64 = (0..3).map(|i| (frame.origin_mm[i] - stock.origin[i]) * n[i]).sum();
                    (0..3).map(|i| normal[i] * n[i]).sum::<f64>().abs() > 0.999 && (along - SIZE[2]).abs() < 1e-7
                })
        })
        .expect("top face");
    let created = sketch::prepare_new(
        session,
        session.revision(),
        actor(),
        "Face sketch",
        CreationPlane::Face {
            snapshot,
            key: element.key,
            anchor: Some(&element.frames[0].origin),
        },
        PlaneOptions::default(),
    )
    .unwrap();
    let (id, frame) = (created.sketch, created.frame);
    session.commit(created.prepared).unwrap();
    (id, frame)
}
fn authored(document: &Document, sketch: EntityId) -> Source {
    inline_input::authored(document, sketch).unwrap()
}
/// Solved point positions of the sketch under a fresh native compile.
fn solved(document: &Document, revision: u64, sketch: EntityId) -> BTreeMap<EntityId, [f64; 2]> {
    let native = Compiler::default().compile(document, revision).unwrap();
    let evaluated = multicad_solver::evaluate_document_sketch(
        document,
        revision,
        sketch,
        Some(&native.snapshot().references),
    )
    .unwrap();
    evaluated
        .solved()
        .expect("the sketch solves")
        .points
        .iter()
        .map(|p| (p.id, p.position_mm))
        .collect()
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
fn close(a: [f64; 2], b: [f64; 2]) -> bool {
    (a[0] - b[0]).abs() < 1e-7 && (a[1] - b[1]).abs() < 1e-7
}

struct Host {
    session: Session,
    editor: SketchEditor,
    feed: SceneFeed,
    ctx: egui::Context,
    sketch: EntityId,
    time: f64,
}
impl Host {
    fn new(session: Session, sketch: EntityId) -> Self {
        Self::in_context(session, sketch, None)
    }
    /// Opens the sketch, optionally in the context of one placed occurrence.
    fn in_context(session: Session, sketch: EntityId, occurrence: Option<EntityId>) -> Self {
        let mut feed = SceneFeed::default();
        feed.set_backend(Box::new(NativeBackend::new()));
        let mut editor = SketchEditor::default();
        editor.set_external_backend(Box::new(NativeBackend::new()));
        editor.grid = false;
        editor.snap = false;
        let mut host = Self {
            session,
            editor,
            feed,
            ctx: Default::default(),
            sketch,
            time: 0.,
        };
        host.wait("the native scene", |h| {
            h.feed
                .is_current(h.session.document().id(), h.session.revision())
        });
        if let Some(instance) = occurrence {
            host.editor.tree_selection = std::collections::BTreeSet::from([instance]);
        }
        host.wait("the face sketch to open", |h| {
            h.editor.editing
                || h.editor
                    .open_at(h.session.document(), h.sketch, h.session.revision())
                    .is_ok()
        });
        host.editor.set_tool(DrawTool::Rectangle).unwrap();
        host
    }
    fn frame(&mut self, events: Vec<egui::Event>) {
        self.time += 0.02;
        self.feed
            .tick(self.session.document(), self.session.revision());
        let raw = egui::RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1200., 800.))),
            events,
            time: Some(self.time),
            ..Default::default()
        };
        let (editor, session, feed) = (&mut self.editor, &mut self.session, &mut self.feed);
        self.ctx
            .run_ui(raw, |ui| {
                editor.show_project(ui, session, true, feed, &Default::default())
            })
            .drop_without_applying_deltas();
        self.editor.poll_model(&self.session, &Default::default());
    }
    fn wait(&mut self, what: &str, mut ready: impl FnMut(&mut Self) -> bool) {
        let deadline = Instant::now() + Duration::from_secs(30);
        while !ready(self) {
            assert!(
                Instant::now() < deadline,
                "Timed out waiting for {what}: {} · {}",
                self.editor.message,
                self.diagnostics()
            );
            self.frame(vec![]);
            std::thread::sleep(Duration::from_millis(1));
        }
    }
    fn screen(&self, world: [f64; 3]) -> Pos2 {
        let p = self
            .editor
            .camera
            .project(world.map(|v| v as f32), self.editor.viewport_rect.unwrap());
        Pos2::new(p[0], p[1])
    }
    /// Readiness of the drawing canvas, reported when a wait times out.
    fn diagnostics(&self) -> String {
        let (document, revision) = (self.session.document(), self.session.revision());
        format!(
            "revision {revision}; scene {:?}; feed current {}; sketch {:?}; pose wait {}; \
             {} corners, {} edges; shown {:?}; pending edits {}; tool {:?}; camera {:?}; \
             viewport {:?}; last pointer {:?}",
            self.editor
                .display_scene
                .as_ref()
                .map(|scene| scene.stamp()),
            self.feed.is_current(document.id(), revision),
            self.editor
                .read_sketch(document, self.sketch, revision)
                .err()
                .map(|e| e.to_string()),
            self.editor.waiting_for_sketch_pose(self.sketch),
            self.editor
                .support_corners(&self.session, self.sketch)
                .len(),
            self.editor.support_edges(&self.session, self.sketch).len(),
            self.editor
                .snap_assist
                .shown()
                .map(|target| (target.kind, target.position)),
            self.editor.has_pending_edits(),
            self.editor.tool,
            self.editor.camera,
            self.editor.viewport_rect,
            self.editor.inline.last_pointer,
        )
    }
    /// Hovers `world` until the canvas shows a support snap there. It waits
    /// for the canvas to draw (current solve, current scene, offered face
    /// targets), then projects `world` with the camera of every frame: the
    /// editor may still reframe the model once the scene and viewport exist.
    fn hover_support(&mut self, world: [f64; 3], kind: TargetKind) {
        self.wait("a drawable face sketch with face targets", |h| {
            let (document, revision) = (h.session.document(), h.session.revision());
            h.editor.read_sketch(document, h.sketch, revision).is_ok()
                && !(h.editor.support_corners(&h.session, h.sketch).is_empty()
                    && h.editor.support_edges(&h.session, h.sketch).is_empty())
        });
        self.wait("the face snap target", |h| {
            let at = h.screen(world);
            h.frame(vec![egui::Event::PointerMoved(at)]);
            h.editor
                .snap_assist
                .shown()
                .is_some_and(|target| target.kind == kind)
        });
    }
    fn click(&mut self, world: [f64; 3]) {
        let at = self.screen(world);
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(at),
                egui::Event::PointerButton {
                    pos: at,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ]);
        }
    }
    fn settle_support(&mut self) {
        self.wait("the projected reference and its drawing", |h| {
            !h.editor.projection_pending() && !h.editor.external_active()
        });
    }
}

/// The rectangle's corner points and the reference outputs they reuse.
fn rectangle_corners(document: &Document, sketch: EntityId) -> Vec<EntityId> {
    let source = authored(document, sketch);
    let mut corners: Vec<EntityId> = source
        .sketch
        .curves
        .iter()
        .filter(|c| !c.construction)
        .flat_map(|c| match c.kind {
            CurveKind::Line { start, end } => [start, end],
            _ => panic!("rectangle sides are lines"),
        })
        .collect();
    corners.sort();
    corners.dedup();
    corners
}

#[test]
fn rectangle_on_rotated_part_face_corners_projects_both_references_in_one_undoable_edit() {
    let mut session = Session::new(Document::new("Snap-time projection").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, frame) = top_face_sketch(&mut session, &stock);
    let mut host = Host::new(session, sketch);
    let (a, c) = (stock.corner(0., 0.), stock.corner(SIZE[0], SIZE[1]));
    let before = host.session.document().clone();
    let revision = host.session.revision();

    host.hover_support(a, TargetKind::FaceCorner);
    assert_eq!(host.session.document(), &before, "hovering changes nothing");
    host.click(a);
    assert_eq!(
        host.session.revision(),
        revision,
        "the first corner only starts the gesture"
    );
    host.hover_support(c, TargetKind::FaceCorner);
    host.click(c);
    host.settle_support();
    assert!(
        host.editor.message.starts_with("Face reference projected"),
        "{}",
        host.editor.message
    );
    assert_eq!(
        host.session.revision(),
        revision + 1,
        "both references and the rectangle are one transaction"
    );
    let document = host.session.document().clone();
    let bindings = sketch_external::bindings(&document, sketch).unwrap();
    assert_eq!(bindings.len(), 2);
    assert!(
        bindings
            .iter()
            .all(|b| b.intent.kind == TopologyKind::Vertex
                && b.definition.behavior == sketch_external::Behavior::ReferenceOnly)
    );
    let outputs: Vec<EntityId> = bindings
        .iter()
        .flat_map(|b| b.definition.output.points())
        .collect();
    let corners = rectangle_corners(&document, sketch);
    assert_eq!(corners.len(), 4);
    assert!(
        outputs.iter().all(|id| corners.contains(id)),
        "the rectangle reuses both projected corner points by identity"
    );
    let positions = solved(&document, host.session.revision(), sketch);
    let expected: Vec<[f64; 2]> = [a, c].map(|p| frame.coordinates(p)).to_vec();
    for output in &outputs {
        assert!(
            expected.iter().any(|e| close(positions[output], *e)),
            "projected corner {:?} is exactly a face corner {expected:?}",
            positions[output]
        );
    }

    // Save and reopen keeps both bindings and the exact corners.
    let reloaded = reopened(&document);
    assert_eq!(
        sketch_external::bindings(&reloaded, sketch).unwrap(),
        bindings
    );
    let again = solved(&reloaded, 0, sketch);
    for output in &outputs {
        assert!(close(again[output], positions[output]));
    }

    // One Undo removes both references and the rectangle; Redo restores them.
    let after = document.clone();
    host.session.undo(host.session.revision(), actor()).unwrap();
    assert_eq!(host.session.document(), &before);
    host.session.redo(host.session.revision(), actor()).unwrap();
    assert_eq!(host.session.document(), &after);

    // A supported upstream edit moves the projected corner and the rectangle with it.
    let grown = model::prepare_primitive_edit(
        &host.session,
        host.session.revision(),
        actor(),
        stock.request.feature,
        &PrimitiveParameters::Box {
            size: [length(40.), length(25.), length(SIZE[2])],
            frame: PrimitiveFrame {
                origin: stock.origin.map(length),
                u: stock.u,
                v: stock.v,
            },
        },
    )
    .unwrap();
    host.session.commit(grown.prepared).unwrap();
    let moved = solved(host.session.document(), host.session.revision(), sketch);
    let far = frame.coordinates(stock.corner(40., 25.));
    assert!(
        outputs.iter().any(|id| close(moved[id], far)),
        "the far projected corner follows the grown face: {:?}",
        outputs.iter().map(|id| moved[id]).collect::<Vec<_>>()
    );

    // Replacing the box removes its corners: the references report it.
    let replaced = model::prepare_primitive_edit(
        &host.session,
        host.session.revision(),
        actor(),
        stock.request.feature,
        &PrimitiveParameters::Cylinder {
            radius: length(5.),
            height: length(SIZE[2]),
            frame: Default::default(),
        },
    )
    .unwrap();
    host.session.commit(replaced.prepared).unwrap();
    let native = Compiler::default()
        .compile(host.session.document(), host.session.revision())
        .unwrap();
    let report = &native.snapshot().references;
    let references: BTreeSet<EntityId> = bindings.iter().map(|b| b.reference).collect();
    let lost: Vec<_> = report
        .entries
        .iter()
        .filter(|d| references.contains(&d.address.owner))
        .map(|d| d.resolution.clone())
        .collect();
    assert!(!lost.is_empty(), "both references are inspected");
    assert!(
        lost.iter()
            .all(|r| matches!(r, Resolution::Missing | Resolution::Ambiguous { .. })),
        "lost corners are reported explicitly: {lost:?}"
    );
}

#[test]
fn a_project_change_while_projecting_creates_nothing_and_datum_sketches_offer_no_face_targets() {
    let mut session = Session::new(Document::new("Snap-time refusal").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, _) = top_face_sketch(&mut session, &stock);
    let mut host = Host::new(session, sketch);
    let (a, c) = (stock.corner(0., 0.), stock.corner(SIZE[0], SIZE[1]));
    host.hover_support(a, TargetKind::FaceCorner);
    host.click(a);
    host.hover_support(c, TargetKind::FaceCorner);
    // The release frame defers the drawing; the references are captured
    // from the next frame on. An unrelated commit lands in between.
    host.click(c);
    assert!(host.editor.projection_pending(), "{}", host.editor.message);
    let mut renamed = host
        .session
        .document()
        .get(stock.request.feature)
        .unwrap()
        .clone();
    renamed.name = "Renamed stock".into();
    host.session
        .apply(Request {
            expected_revision: host.session.revision(),
            actor: actor(),
            label: "Rename".into(),
            edits: vec![multicad_core::document::Mutation::Replace(renamed)],
        })
        .unwrap();
    let changed = host.session.document().clone();
    host.settle_support();
    assert!(
        host.editor.message.ends_with("nothing was placed"),
        "{}",
        host.editor.message
    );
    assert_eq!(
        host.session.document(),
        &changed,
        "no reference and no rectangle"
    );
    assert!(
        sketch_external::bindings(&changed, sketch)
            .unwrap()
            .is_empty()
    );

    // A datum sketch in the same plane is not face-attached: no face targets.
    let mut session = Session::new(Document::new("Datum sketch").unwrap()).unwrap();
    let flat = PrimitiveRequest::new(
        "Flat stock",
        PrimitiveParameters::Box {
            size: SIZE.map(length),
            frame: Default::default(),
        },
    );
    session
        .commit(
            model::prepare_primitive(&session, 0, actor(), &flat)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let (datum, edits) = Source::create(multicad_core::sketch::Datum::Xy, "Datum").unwrap();
    session
        .apply(Request {
            expected_revision: session.revision(),
            actor: actor(),
            label: "Datum sketch".into(),
            edits,
        })
        .unwrap();
    let mut host = Host::new(session, datum);
    let at = host.screen([0., 0., 0.]);
    host.frame(vec![egui::Event::PointerMoved(at)]);
    host.frame(vec![egui::Event::PointerMoved(at)]);
    assert!(host.editor.support_corners(&host.session, datum).is_empty());
    assert!(host.editor.support_edges(&host.session, datum).is_empty());
}

#[test]
fn placed_occurrence_face_corners_project_in_producer_coordinates_in_one_undoable_edit() {
    let mut session = Session::new(Document::new("Placed snap projection").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, _) = top_face_sketch(&mut session, &stock);
    let placement = ComponentPlacementRequest::new(
        stock.request.body,
        NumericTransform {
            translation: [length(150.), length(40.), length(25.)],
            axis: [0., 1., 1.],
            angle: model::angle(35.),
            ..Default::default()
        },
    );
    session
        .commit(
            model::prepare_place_component(&session, session.revision(), actor(), &placement)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let mut host = Host::in_context(session, sketch, Some(placement.instance));
    // The occurrence pose from the evaluated assembly, as the editor shows it.
    host.wait("the placed sketch context", |h| {
        h.editor.sketch_pose(sketch).is_some()
    });
    let pose = host.editor.sketch_pose(sketch).unwrap();
    assert_ne!(pose, crate::viewport::Transform::default());
    let (a, c) = (
        pose.point(stock.corner(0., 0.)),
        pose.point(stock.corner(SIZE[0], SIZE[1])),
    );
    let before = host.session.document().clone();
    let revision = host.session.revision();
    host.hover_support(a, TargetKind::FaceCorner);
    host.click(a);
    host.hover_support(c, TargetKind::FaceCorner);
    host.click(c);
    host.settle_support();
    assert!(
        host.editor.message.starts_with("Face reference projected"),
        "{}",
        host.editor.message
    );
    assert_eq!(host.session.revision(), revision + 1);
    let document = host.session.document().clone();
    // Corner and sketch share producer-local coordinates (merged design), so
    // the placed occurrence projects ordinary references.
    let bindings = sketch_external::bindings(&document, sketch).unwrap();
    assert_eq!(bindings.len(), 2, "two projected references");
    let outputs: Vec<EntityId> = bindings
        .iter()
        .flat_map(|b| b.definition.output.points())
        .collect();
    let corners = rectangle_corners(&document, sketch);
    assert!(outputs.iter().all(|id| corners.contains(id)));
    let reloaded = reopened(&document);
    assert_eq!(
        sketch_external::bindings(&reloaded, sketch).unwrap(),
        bindings
    );
    host.session.undo(host.session.revision(), actor()).unwrap();
    assert_eq!(host.session.document(), &before);
    host.session.redo(host.session.revision(), actor()).unwrap();
    assert_eq!(host.session.document(), &document);
}

#[test]
fn a_line_ending_on_a_face_edge_projects_that_edge_and_ends_exactly_on_it() {
    let mut session = Session::new(Document::new("Edge snap projection").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, frame) = top_face_sketch(&mut session, &stock);
    let mut host = Host::new(session, sketch);
    host.editor.set_tool(DrawTool::Line).unwrap();
    let revision = host.session.revision();
    let a = stock.corner(0., 0.);
    // A point on the far long edge (y = SIZE[1]), away from its corners.
    let on_edge = stock.corner(17., SIZE[1]);
    host.hover_support(a, TargetKind::FaceCorner);
    host.click(a);
    host.hover_support(on_edge, TargetKind::FaceEdge);
    let shown = host.editor.snap_assist.shown().unwrap().position;
    host.click(on_edge);
    host.settle_support();
    assert!(
        host.editor.message.starts_with("Face reference projected"),
        "{}",
        host.editor.message
    );
    assert_eq!(host.session.revision(), revision + 1);
    let document = host.session.document().clone();
    let bindings = sketch_external::bindings(&document, sketch).unwrap();
    let kinds: BTreeSet<_> = bindings.iter().map(|b| b.intent.kind).collect();
    assert_eq!(
        kinds,
        BTreeSet::from([TopologyKind::Vertex, TopologyKind::Edge]),
        "the corner and the edge are both projected"
    );
    let positions = solved(&document, host.session.revision(), sketch);
    let source = authored(&document, sketch);
    let line = source
        .sketch
        .curves
        .iter()
        .find(|c| !c.construction)
        .unwrap();
    let CurveKind::Line { start, end } = line.kind else {
        panic!("a line")
    };
    let corner = bindings
        .iter()
        .find(|b| b.intent.kind == TopologyKind::Vertex)
        .unwrap()
        .definition
        .output
        .points()[0];
    assert_eq!(start, corner, "the line starts on the projected corner");
    let edge_far = frame.coordinates(stock.corner(SIZE[0], SIZE[1]));
    let edge_near = frame.coordinates(stock.corner(0., SIZE[1]));
    let p = positions[&end];
    let cross = (edge_far[0] - edge_near[0]) * (p[1] - edge_near[1])
        - (edge_far[1] - edge_near[1]) * (p[0] - edge_near[0]);
    assert!(cross.abs() < 1e-7, "the end lies exactly on the face edge");
    assert!(close(p, shown), "the committed end is the shown snap");
}

#[test]
fn escape_while_the_face_reference_is_projected_creates_nothing() {
    let mut session = Session::new(Document::new("Snap-time cancel").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, _) = top_face_sketch(&mut session, &stock);
    let mut host = Host::new(session, sketch);
    let before = host.session.document().clone();
    let (a, c) = (stock.corner(0., 0.), stock.corner(SIZE[0], SIZE[1]));
    host.hover_support(a, TargetKind::FaceCorner);
    host.click(a);
    host.hover_support(c, TargetKind::FaceCorner);
    host.click(c);
    assert!(host.editor.projection_pending(), "{}", host.editor.message);
    assert!(
        !host.editor.message.contains("committed"),
        "nothing is reported committed before the commit: {}",
        host.editor.message
    );
    host.editor.cancel();
    host.settle_support();
    assert_eq!(
        host.editor.message,
        "Face reference projection cancelled; nothing was placed"
    );
    for _ in 0..3 {
        host.frame(vec![]);
    }
    assert_eq!(host.session.document(), &before);
}

/// Output points of the sketch's projected corner references.
fn corner_outputs(document: &Document, sketch: EntityId) -> Vec<EntityId> {
    sketch_external::bindings(document, sketch)
        .unwrap()
        .iter()
        .filter(|b| b.intent.kind == TopologyKind::Vertex)
        .flat_map(|b| b.definition.output.points())
        .collect()
}
fn drawn_lines(document: &Document, sketch: EntityId) -> Vec<(EntityId, EntityId)> {
    authored(document, sketch)
        .sketch
        .curves
        .iter()
        .filter(|c| !c.construction)
        .filter_map(|c| match c.kind {
            CurveKind::Line { start, end } => Some((start, end)),
            _ => None,
        })
        .collect()
}

#[test]
fn a_polyline_segment_between_face_corners_commits_once_and_the_chain_continues() {
    let mut session = Session::new(Document::new("Polyline projection").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, frame) = top_face_sketch(&mut session, &stock);
    let mut host = Host::new(session, sketch);
    host.editor.set_tool(DrawTool::Polyline).unwrap();
    let revision = host.session.revision();
    let (a, c) = (stock.corner(0., 0.), stock.corner(SIZE[0], SIZE[1]));
    host.hover_support(a, TargetKind::FaceCorner);
    host.click(a);
    host.hover_support(c, TargetKind::FaceCorner);
    host.click(c);
    assert!(host.editor.projection_pending(), "{}", host.editor.message);
    assert!(
        !host.editor.message.contains("committed"),
        "{}",
        host.editor.message
    );
    host.settle_support();
    assert_eq!(
        host.session.revision(),
        revision + 1,
        "{}",
        host.editor.message
    );
    let document = host.session.document().clone();
    let outputs = corner_outputs(&document, sketch);
    assert_eq!(outputs.len(), 2);
    let lines = drawn_lines(&document, sketch);
    assert_eq!(lines.len(), 1);
    assert!(outputs.contains(&lines[0].0) && outputs.contains(&lines[0].1));
    // The chain continues from the projected far corner.
    let next = host.editor.pending_first().expect("the polyline continues");
    assert!(close(next, frame.coordinates(c)), "{next:?}");
    // The canvas places the next vertex under the pointer before the click.
    let free = stock.corner(10., 15.);
    let at = host.screen(free);
    let expected = frame.coordinates(free);
    host.wait("the next polyline vertex under the pointer", |h| {
        h.frame(vec![egui::Event::PointerMoved(at)]);
        h.editor
            .inline
            .last_point
            .is_some_and(|p| (p[0] - expected[0]).portable_hypot(p[1] - expected[1]) < 1e-3)
    });
    host.click(free);
    assert_eq!(
        host.session.revision(),
        revision + 2,
        "{}",
        host.editor.message
    );
    let lines = drawn_lines(host.session.document(), sketch);
    assert_eq!(lines.len(), 2);
    assert_eq!(
        corner_outputs(host.session.document(), sketch).len(),
        2,
        "the free vertex projects nothing"
    );
    assert!(
        lines
            .iter()
            .any(|(start, _)| *start == lines[0].1 || *start == lines[0].0),
        "the second segment starts on the projected corner: {lines:?}"
    );
}

#[test]
fn an_arc_centred_on_a_face_corner_reuses_both_projected_corners() {
    let mut session = Session::new(Document::new("Arc projection").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, _) = top_face_sketch(&mut session, &stock);
    let mut host = Host::new(session, sketch);
    host.editor.set_tool(DrawTool::Arc).unwrap();
    let revision = host.session.revision();
    let (a, c) = (stock.corner(0., 0.), stock.corner(SIZE[0], SIZE[1]));
    host.hover_support(a, TargetKind::FaceCorner);
    host.click(a);
    host.hover_support(c, TargetKind::FaceCorner);
    host.click(c);
    let radius = SIZE[0].portable_hypot(SIZE[1]);
    let angle = 80f64.to_radians();
    host.click(stock.corner(radius * angle.portable_cos(), radius * angle.portable_sin()));
    host.settle_support();
    assert_eq!(
        host.session.revision(),
        revision + 1,
        "{}",
        host.editor.message
    );
    let document = host.session.document().clone();
    let outputs = corner_outputs(&document, sketch);
    assert_eq!(outputs.len(), 2);
    let arc = authored(&document, sketch)
        .sketch
        .curves
        .iter()
        .find_map(|c| match c.kind {
            CurveKind::Arc { center, start, .. } if !c.construction => Some((center, start)),
            _ => None,
        })
        .expect("an arc");
    assert!(outputs.contains(&arc.0) && outputs.contains(&arc.1));
    host.session.undo(host.session.revision(), actor()).unwrap();
    assert!(
        sketch_external::bindings(host.session.document(), sketch)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn a_three_point_rectangle_on_two_face_corners_is_one_transaction() {
    let mut session = Session::new(Document::new("Multi-click projection").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, _) = top_face_sketch(&mut session, &stock);
    let mut host = Host::new(session, sketch);
    host.editor.set_tool(DrawTool::RectangleThree).unwrap();
    let revision = host.session.revision();
    let (a, b) = (stock.corner(0., 0.), stock.corner(SIZE[0], 0.));
    host.hover_support(a, TargetKind::FaceCorner);
    host.click(a);
    host.hover_support(b, TargetKind::FaceCorner);
    host.click(b);
    host.click(stock.corner(15., 8.));
    assert!(
        !host.editor.message.contains("committed"),
        "{}",
        host.editor.message
    );
    host.settle_support();
    assert_eq!(
        host.session.revision(),
        revision + 1,
        "{}",
        host.editor.message
    );
    let document = host.session.document().clone();
    let outputs = corner_outputs(&document, sketch);
    assert_eq!(outputs.len(), 2);
    let corners = rectangle_corners(&document, sketch);
    assert_eq!(corners.len(), 4);
    assert!(outputs.iter().all(|id| corners.contains(id)));
}

#[test]
fn an_entered_width_joins_the_projection_transaction() {
    let mut session = Session::new(Document::new("Sized projection").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, _) = top_face_sketch(&mut session, &stock);
    let mut host = Host::new(session, sketch);
    let revision = host.session.revision();
    let a = stock.corner(0., 0.);
    host.hover_support(a, TargetKind::FaceCorner);
    host.click(a);
    host.wait("the size fields", |h| {
        h.editor.inline_stage() == Some((DrawTool::Rectangle, true))
            && h.editor.inline.texts.len() == 2
    });
    host.editor.inline.texts[0] = "12".into();
    host.click(stock.corner(20., 15.));
    host.settle_support();
    assert_eq!(
        host.session.revision(),
        revision + 1,
        "references, rectangle and size are one commit: {}",
        host.editor.message
    );
    let document = host.session.document().clone();
    let outputs = corner_outputs(&document, sketch);
    assert_eq!(outputs.len(), 1);
    assert!(rectangle_corners(&document, sketch).contains(&outputs[0]));
    let source = authored(&document, sketch);
    let sized = source.sketch.constraints.iter().any(|c| {
        matches!(&c.kind, ConstraintKind::Length {
            dimension: Dimension::Literal { value, .. }, ..
        } if (value - 12.).abs() < 1e-12)
    });
    assert!(sized, "the entered width is a driving dimension");
    host.session.undo(host.session.revision(), actor()).unwrap();
    assert!(
        authored(host.session.document(), sketch)
            .sketch
            .curves
            .is_empty(),
        "one Undo removes the reference, the rectangle and its size"
    );
}

#[test]
fn missing_and_ambiguous_face_references_are_refused_explicitly() {
    let mut session = Session::new(Document::new("Refused projection").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, _) = top_face_sketch(&mut session, &stock);
    let mut host = Host::new(session, sketch);
    host.wait("the face corners", |h| {
        !h.editor.support_corners(&h.session, h.sketch).is_empty()
    });
    let before = host.session.document().clone();
    let gone = TopologyKey {
        kind: TopologyKind::Vertex,
        ordinal: 9_999,
    };
    let error = host
        .editor
        .draw_on_corners(
            &host.session,
            DrawTool::Line,
            (PointInput::Position([0., 0.]), Some(gone)),
            (PointInput::Position([5., 5.]), None),
            0.,
            Channel::Ui,
        )
        .unwrap_err();
    assert!(error.starts_with("Missing face reference"), "{error}");
    assert!(!host.editor.projection_pending());
    assert_eq!(host.session.document(), &before);
    let key = |ordinal| TopologyKey {
        kind: TopologyKind::Edge,
        ordinal,
    };
    let ambiguous = support_projection::resolution_error(&Resolution::Ambiguous {
        candidates: vec![key(1), key(2)],
    })
    .unwrap();
    assert!(
        ambiguous.starts_with("Ambiguous face reference (2 projections)"),
        "{ambiguous}"
    );
    assert!(
        support_projection::resolution_error(&Resolution::Missing)
            .unwrap()
            .starts_with("Missing face reference")
    );
}

#[test]
fn spline_poles_on_face_corners_reuse_the_projected_points_in_one_transaction() {
    let mut session = Session::new(Document::new("Spline projection").unwrap()).unwrap();
    let stock = Stock::new(&mut session);
    let (sketch, _) = top_face_sketch(&mut session, &stock);
    let mut host = Host::new(session, sketch);
    host.editor.set_tool(DrawTool::Spline).unwrap();
    let revision = host.session.revision();
    let (a, c) = (stock.corner(0., 0.), stock.corner(SIZE[0], SIZE[1]));
    host.hover_support(a, TargetKind::FaceCorner);
    host.click(a);
    for free in [stock.corner(10., 12.), stock.corner(20., 6.)] {
        let poles = host.editor.spline_poles.len();
        host.click(free);
        assert_eq!(host.editor.spline_poles.len(), poles + 1);
    }
    host.hover_support(c, TargetKind::FaceCorner);
    host.click(c);
    assert_eq!(host.editor.spline_poles.len(), 4);
    assert_eq!(host.session.revision(), revision, "poles change nothing");
    host.editor
        .finish_spline(&mut host.session, false, Channel::Ui)
        .unwrap();
    assert!(host.editor.projection_pending(), "{}", host.editor.message);
    assert!(
        !host.editor.message.contains("committed"),
        "{}",
        host.editor.message
    );
    host.settle_support();
    assert_eq!(
        host.session.revision(),
        revision + 1,
        "{}",
        host.editor.message
    );
    let document = host.session.document().clone();
    let outputs = corner_outputs(&document, sketch);
    assert_eq!(outputs.len(), 2);
    let poles = authored(&document, sketch)
        .sketch
        .curves
        .iter()
        .find_map(|c| match &c.kind {
            CurveKind::BSpline { poles, .. } if !c.construction => Some(poles.clone()),
            _ => None,
        })
        .expect("a spline");
    assert_eq!(poles.len(), 4);
    assert!(
        outputs.contains(&poles[0]) && outputs.contains(&poles[3]),
        "the end poles are the projected corners: {poles:?} {outputs:?}"
    );
    host.session.undo(host.session.revision(), actor()).unwrap();
    assert!(
        sketch_external::bindings(host.session.document(), sketch)
            .unwrap()
            .is_empty(),
        "one Undo removes the references and the spline"
    );
}
