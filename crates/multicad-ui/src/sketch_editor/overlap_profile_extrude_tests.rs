//! CAD-400: the product owner's sketch. Rectangle A is 100 × 60 mm.
//! Rectangle B starts at the point given by A's top-edge midpoint (x) and
//! A's center (y), and ends on A's bottom-right corner. B's top and left
//! edges end on A's edges, and its right and bottom edges lie on A's edges.
//! Free corners are filleted. After Finish Sketch, both closed areas must be
//! pickable in the 3D view and extrudable in either order.
use super::tests::{click_sketch_modifiers, input, session, settle};
use super::*;
use crate::modeling::native_fixture::{NativeBackend, scene};
use multicad_core::{
    constraint_sketch::CurveKind,
    sketch_draw::{DrawInput, DrawRequest, PointInput},
};
use multicad_kernel::NativeRegionQueries;
use multicad_model::Compiler;
use std::{
    f64::consts::PI,
    sync::Arc,
    time::{Duration, Instant},
};

const RADIUS: f64 = 5.;
/// Material a 90° fillet of `RADIUS` removes from a corner.
const CORNER: f64 = RADIUS * RADIUS * (1. - PI / 4.);
/// L-shaped part of A outside B: three filleted corners of A, and it gains
/// what the fillet on B's inner corner removes from B.
const L_AREA: f64 = 100. * 60. - 50. * 30. - 3. * CORNER + CORNER;
/// B's quarter, with its free inner corner filleted.
const B_AREA: f64 = 50. * 30. - CORNER;
const L_POINT: [f64; 2] = [25., 30.];
const B_POINT: [f64; 2] = [75., 15.];

fn corner_id(session: &Session, sketch: EntityId, at: [f64; 2]) -> EntityId {
    let (source, _, _) = Source::read(session.document(), sketch).unwrap();
    source
        .sketch
        .points
        .iter()
        .find(|p| {
            (p.position_mm[0] - at[0]).abs() < 1e-6 && (p.position_mm[1] - at[1]).abs() < 1e-6
        })
        .unwrap()
        .id
}

struct Reported {
    session: Session,
    editor: SketchEditor,
    sketch: EntityId,
    ctx: egui::Context,
}

impl Reported {
    fn new() -> Self {
        let mut session = session();
        let mut editor = SketchEditor::default();
        editor.set_query_factory(Arc::new(|| Ok(Box::new(NativeRegionQueries::default()))));
        editor.set_model_backend(Box::new(NativeBackend::new()));
        editor.create(&mut session, Datum::Xy, Channel::Ui).unwrap();
        editor.snap = false;
        editor
            .add(
                &mut session,
                DrawTool::Rectangle,
                [0., 0.],
                [100., 60.],
                Channel::Ui,
            )
            .unwrap();
        let sketch = editor.selected.unwrap();
        // The second corner snaps onto A's existing corner point.
        let shared = corner_id(&session, sketch, [100., 0.]);
        settle(&mut editor, &session);
        editor
            .draw(
                &mut session,
                &DrawRequest {
                    input: DrawInput::Rectangle {
                        a: PointInput::Position([50., 30.]),
                        b: PointInput::Existing(shared),
                    },
                    construction: false,
                },
                Channel::Ui,
            )
            .unwrap();
        let ctx = egui::Context::default();
        ctx.all_styles_mut(|style| style.animation_time = 0.);
        editor.set_tool(DrawTool::Fillet).unwrap();
        editor.corner_size_mm = RADIUS;
        for _ in 0..2 {
            settle(&mut editor, &session);
            ctx.run_ui(input(vec![]), |ui| {
                editor.canvas(ui, &mut session, true, None, false)
            })
            .drop_without_applying_deltas();
        }
        for corner in [[0., 0.], [0., 60.], [100., 60.], [50., 30.]] {
            let before = session.revision();
            click_sketch_modifiers(
                &ctx,
                &mut editor,
                &mut session,
                corner,
                egui::Modifiers::NONE,
            );
            assert_eq!(
                session.revision(),
                before + 1,
                "fillet at {corner:?}: {}",
                editor.message
            );
        }
        let (source, _, _) = Source::read(session.document(), sketch).unwrap();
        let arcs = source
            .sketch
            .curves
            .iter()
            .filter(|c| !c.construction && matches!(c.kind, CurveKind::Arc { .. }))
            .count();
        assert_eq!(arcs, 4);
        editor.tool = None;
        settle(&mut editor, &session);
        editor.finish().unwrap();
        Self {
            session,
            editor,
            sketch,
            ctx,
        }
    }
    fn scene(&self) -> Arc<crate::viewport::Scene> {
        let native = Compiler::default()
            .compile(self.session.document(), self.session.revision())
            .unwrap();
        scene(self.session.document(), native.snapshot(), 1).unwrap()
    }
    fn frame(&mut self, events: Vec<egui::Event>) {
        let scene = self.scene();
        self.editor.display_scene = Some(scene.clone());
        settle(&mut self.editor, &self.session);
        self.ctx
            .run_ui(input(events), |ui| {
                self.editor
                    .canvas(ui, &mut self.session, false, Some(scene.clone()), true)
            })
            .drop_without_applying_deltas();
    }
    fn screen(&mut self, point: [f64; 2]) -> Pos2 {
        self.frame(vec![]);
        let frame = self
            .editor
            .read_sketch(
                self.session.document(),
                self.sketch,
                self.session.revision(),
            )
            .unwrap()
            .1;
        let p = self.editor.camera.project(
            frame.point(point).map(|v| v as f32),
            self.editor.viewport_rect.unwrap(),
        );
        Pos2::new(p[0], p[1])
    }
    fn hover(&mut self, point: [f64; 2]) {
        let pos = self.screen(point);
        self.frame(vec![egui::Event::PointerMoved(pos)]);
        self.drain();
    }
    fn click(&mut self, point: [f64; 2], shift: bool) {
        let pos = self.screen(point);
        let modifiers = egui::Modifiers {
            shift,
            ..Default::default()
        };
        for pressed in [true, false] {
            let mut events = vec![
                egui::Event::ModifiersChanged(modifiers),
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers,
                },
            ];
            if !pressed {
                events.push(egui::Event::ModifiersChanged(egui::Modifiers::NONE));
            }
            self.frame(events);
        }
        self.drain();
    }
    fn drain(&mut self) {
        let deadline = Instant::now() + Duration::from_secs(10);
        loop {
            self.editor.sketch_evaluation.poll();
            self.editor
                .poll_profile_pick(self.session.document(), self.session.revision());
            if self.editor.pending_profile_pick.is_none()
                && self.editor.queued_profile_clicks.is_empty()
                && !self.editor.sketch_evaluation.pending()
            {
                break;
            }
            assert!(Instant::now() < deadline, "Profile click did not drain");
            std::thread::sleep(Duration::from_millis(1));
        }
    }
    fn areas(&self) -> Vec<f64> {
        let mut areas: Vec<_> = self
            .editor
            .sketch_evaluation
            .current(self.sketch)
            .unwrap()
            .regions
            .regions
            .iter()
            .map(|r| r.area_mm2)
            .collect();
        areas.sort_by(f64::total_cmp);
        areas
    }
    fn selected_areas(&self) -> Vec<f64> {
        let presentation = self.editor.sketch_evaluation.current(self.sketch).unwrap();
        let mut areas: Vec<_> = self
            .editor
            .profiles
            .selected()
            .iter()
            .map(|intent| presentation.regions.resolve(intent).unwrap().area_mm2)
            .collect();
        areas.sort_by(f64::total_cmp);
        areas
    }
    /// Preview, then Apply; returns the volume of the new solids.
    fn preview_and_apply(&mut self) -> f64 {
        let resources = Default::default();
        self.editor
            .preview_model(&self.session, &resources, Channel::Ui)
            .unwrap();
        let wait = Instant::now() + Duration::from_secs(30);
        while !self.editor.model.can_apply(&self.session, &resources) {
            self.editor.poll_model(&self.session, &resources);
            assert!(
                self.editor.model.error().is_none(),
                "{:?}",
                self.editor.model.error()
            );
            assert!(Instant::now() < wait, "Extrude preview timed out");
            std::thread::sleep(Duration::from_millis(1));
        }
        let before = self.session.revision();
        let previous: BTreeSet<_> = self.session.document().entities().map(|e| e.id).collect();
        self.editor
            .apply_model(&mut self.session, &resources)
            .unwrap();
        assert_eq!(
            self.session.revision(),
            before + 1,
            "one Extrude transaction"
        );
        let native = Compiler::default()
            .compile(self.session.document(), self.session.revision())
            .unwrap();
        native
            .snapshot()
            .shapes
            .iter()
            .filter(|(id, _)| !previous.contains(id))
            .map(|(_, shape)| shape.volume_mm3)
            .sum()
    }
}

use std::collections::BTreeSet;

fn near(actual: f64, expected: f64, what: &str) {
    assert!(
        (actual - expected).abs() < 1e-3,
        "{what}: {actual} != {expected}"
    );
}

#[test]
fn reported_overlapping_filleted_rectangles_offer_both_closed_areas_after_finish() {
    let mut h = Reported::new();
    h.frame(vec![]);
    let presentation = h.editor.sketch_evaluation.current(h.sketch).unwrap();
    assert_eq!(
        h.areas().len(),
        2,
        "closed areas {:?}, diagnostics {:?}",
        h.areas(),
        presentation.regions.diagnostics
    );
    near(h.areas()[0], B_AREA, "B area");
    near(h.areas()[1], L_AREA, "L area");
    // Hover gives visible feedback in the 3D view before any click.
    h.hover(B_POINT);
    assert!(
        h.editor
            .hovered_profile_area()
            .is_some_and(|a| (a - B_AREA).abs() < 1e-3),
        "{}",
        h.editor.message
    );
    h.click(L_POINT, false);
    assert_eq!(h.selected_areas().len(), 1, "{}", h.editor.message);
    h.click(B_POINT, true);
    let selected = h.selected_areas();
    assert_eq!(selected.len(), 2, "{}", h.editor.message);
    near(selected[0], B_AREA, "selected B");
    near(selected[1], L_AREA, "selected L");
    assert!(h.editor.solid.selected_id().is_none());
}

#[test]
fn reported_sketch_extrudes_selection_first_and_tool_first_with_one_undo_each() {
    let mut h = Reported::new();
    let resources = Default::default();
    // Select, then Extrude.
    h.click(L_POINT, false);
    h.click(B_POINT, true);
    assert_eq!(h.selected_areas().len(), 2, "{}", h.editor.message);
    h.editor
        .begin_extrude(&h.session, &resources, Channel::Ui)
        .unwrap();
    let distance = h
        .editor
        .model
        .extrude_distance_for_tests(h.session.document())
        .unwrap();
    let volume = h.preview_and_apply();
    near(
        volume,
        (L_AREA + B_AREA) * distance.abs(),
        "selection-first volume",
    );
    let applied = h.session.document().clone();
    h.session
        .undo(h.session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert!(
        h.session.document().entities().count() < applied.entities().count(),
        "one Undo removes the Extrude"
    );

    // Extrude first, then pick one area in the 3D view.
    h.editor.profiles.clear();
    h.editor.tree_selection.clear();
    h.frame(vec![]);
    h.editor
        .activate_model_tool(
            &h.session,
            &resources,
            &Default::default(),
            crate::modeling::ModelingTool::Extrude,
            Channel::Ui,
        )
        .unwrap();
    h.click(B_POINT, false);
    let wait = Instant::now() + Duration::from_secs(10);
    while !h.editor.model.active() {
        h.frame(vec![]);
        assert!(
            Instant::now() < wait,
            "Extrude did not start: {}",
            h.editor.message
        );
        std::thread::sleep(Duration::from_millis(1));
    }
    let distance = h
        .editor
        .model
        .extrude_distance_for_tests(h.session.document())
        .unwrap();
    let volume = h.preview_and_apply();
    near(volume, B_AREA * distance.abs(), "tool-first volume");

    // The applied document survives save and reopen with the same solid.
    let project = multicad_core::project_file::ProjectFile {
        document: h.session.document().clone(),
        resources: Default::default(),
        recovery: None,
    };
    let mut bytes = vec![];
    multicad_core::project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let reopened = multicad_core::project_file::read(
        &mut &bytes[..],
        Default::default(),
        &multicad_core::project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    let before = Compiler::default()
        .compile(h.session.document(), h.session.revision())
        .unwrap();
    let after = Compiler::default().compile(&reopened.document, 0).unwrap();
    let total = |native: &multicad_model::Model| -> f64 {
        native
            .snapshot()
            .shapes
            .values()
            .map(|s| s.volume_mm3)
            .sum()
    };
    near(total(&after), total(&before), "reopened volume");
}
