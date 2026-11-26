//! Primitive corner drags driven through real egui input and the production
//! native model worker (CAD-401). Every drag state the pointer can reach must
//! leave a draft that the native validation accepts.
use super::*;
use crate::{camera::Camera, modeling::native_fixture::NativeBackend};
use egui::{Event, Pos2, Rect, Vec2};
use multicad_core::{document::Document, model::ResolvedPrimitive};
use std::time::{Duration, Instant};

const SIZE: Vec2 = Vec2::new(800., 600.);

pub(super) struct PrimitiveHarness {
    pub session: Session,
    pub controller: ModelToolController,
    ctx: egui::Context,
    pub camera: Camera,
    pub rect: Rect,
    pub modifiers: egui::Modifiers,
}

impl PrimitiveHarness {
    pub fn new(camera: Camera) -> Self {
        let session = Session::new(Document::new("Primitive drag limits").unwrap()).unwrap();
        let mut controller = ModelToolController::default();
        controller.set_backend(Box::new(NativeBackend::new()));
        controller
            .begin_primitive(
                &session,
                &Default::default(),
                ModelingTool::Box,
                Channel::Ui,
            )
            .unwrap();
        let mut harness = Self {
            session,
            controller,
            ctx: egui::Context::default(),
            camera,
            rect: Rect::from_min_size(Pos2::ZERO, SIZE),
            modifiers: egui::Modifiers::NONE,
        };
        harness.wait();
        harness
    }
    pub fn screen(&self, p: [f64; 3]) -> Pos2 {
        let s = self.camera.project(p.map(|v| v as f32), self.rect);
        Pos2::new(s[0], s[1])
    }
    pub fn frame(&mut self, mut events: Vec<Event>) -> bool {
        let (camera, rect) = (self.camera, self.rect);
        events.insert(0, Event::ModifiersChanged(self.modifiers));
        let mut owned = false;
        self.ctx
            .run_ui(
                egui::RawInput {
                    focused: true,
                    screen_rect: Some(rect),
                    events,
                    ..Default::default()
                },
                |ui| {
                    ui.interact(
                        rect,
                        ui.id().with("viewport"),
                        egui::Sense::click_and_drag(),
                    );
                    owned = self.controller.primitive_viewport(
                        ui,
                        &self.session,
                        &Default::default(),
                        camera,
                        rect,
                        true,
                    );
                },
            )
            .drop_without_applying_deltas();
        owned
    }
    pub fn press(&mut self, pos: Pos2) {
        let modifiers = self.modifiers;
        assert!(
            self.frame(vec![
                Event::PointerMoved(pos),
                Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed: true,
                    modifiers,
                },
            ]),
            "the corner handle must own the press"
        );
    }
    pub fn drag(&mut self, pos: Pos2) {
        self.frame(vec![Event::PointerMoved(pos)]);
    }
    pub fn release(&mut self, pos: Pos2) {
        let modifiers = self.modifiers;
        self.frame(vec![Event::PointerButton {
            pos,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers,
        }]);
    }
    /// Wait for the native receipt of the current draft.
    pub fn wait(&mut self) {
        let until = Instant::now() + Duration::from_secs(30);
        while self.controller.pending() && Instant::now() < until {
            self.controller.poll(&self.session, &Default::default());
            std::thread::sleep(Duration::from_millis(2));
        }
        assert!(!self.controller.pending(), "native preview did not settle");
    }
    /// Preview the current draft natively and require that it validates.
    pub fn assert_valid_preview(&mut self, state: &str) {
        self.controller
            .preview(&self.session, &Default::default(), Channel::Ui)
            .unwrap_or_else(|error| panic!("{state}: draft request refused: {error}"));
        self.wait();
        assert_eq!(self.controller.error, None, "{state}");
        assert!(
            self.controller
                .preview_visible(&self.session, &Default::default()),
            "{state}: native preview did not validate"
        );
        if !self.controller.transform_dragging() {
            assert!(
                self.controller
                    .can_apply(&self.session, &Default::default()),
                "{state}: validated draft cannot be applied"
            );
        }
    }
    pub fn shape(&self) -> ResolvedPrimitive {
        let Some(ModelRequest::Primitive(request)) =
            self.controller.request_for_tests(&self.session)
        else {
            panic!("primitive request");
        };
        request
            .parameters
            .resolve(self.session.parameters())
            .unwrap()
    }
    pub fn box_shape(&self) -> ([f64; 3], [f64; 3]) {
        let ResolvedPrimitive::Box { size_mm, frame } = self.shape() else {
            panic!("box remains a box");
        };
        (size_mm, frame.origin)
    }
}

/// Top view on the XY plane with an irregular zoom so pointer positions map
/// to long decimal millimetre values.
pub(super) fn top_camera(scale: f32) -> Camera {
    Camera {
        yaw: 0.,
        pitch: std::f32::consts::FRAC_PI_2,
        scale,
        target: [53.17, 29.41, 0.],
        perspective: false,
    }
}

#[test]
fn every_free_primitive_drag_state_previews_natively_without_arithmetic_errors() {
    let mut h = PrimitiveHarness::new(top_camera(3.713));
    h.assert_valid_preview("opened");
    let corner = h.screen([100., 60., 20.]);
    let low = h.screen([0., 0., 20.]);
    h.press(corner);
    // Start of the drag: below the movement threshold, then a first move.
    h.drag(corner + egui::vec2(1.3, 0.7));
    h.assert_valid_preview("drag start");
    h.drag(corner + egui::vec2(37.31, -11.93));
    h.assert_valid_preview("unrounded drag");
    // Through zero size onto the opposite corner, then across it.
    h.drag(low);
    h.assert_valid_preview("zero size");
    let (size, _) = h.box_shape();
    assert!(size.iter().all(|v| *v > 0.), "{size:?}");
    h.drag(low - egui::vec2(41.77, -23.19));
    h.assert_valid_preview("across the start point");
    let (size, _) = h.box_shape();
    assert!(size.iter().all(|v| *v > 0.), "{size:?}");
    // Far outside the viewport and outside model bounds.
    h.drag(Pos2::new(1.0e7, -3.3e6));
    h.assert_valid_preview("far away");
    h.drag(Pos2::new(f32::MAX, f32::MIN));
    h.assert_valid_preview("extreme pointer");
    // Back into a regular size and release.
    let back = h.screen([123.456_789, 71.234_567, 20.]);
    h.drag(back);
    h.release(back);
    h.assert_valid_preview("released");
}

#[test]
fn off_plane_perspective_drag_keeps_the_last_valid_draft() {
    let mut camera = Camera {
        yaw: 0.6,
        pitch: 0.9,
        scale: 3.1,
        target: [50., 30., 10.],
        perspective: true,
    };
    camera.fit(
        &[[0., 0., 0.], [100., 60., 20.]],
        Rect::from_min_size(Pos2::ZERO, SIZE),
    );
    let mut h = PrimitiveHarness::new(camera);
    let corner = h.screen([100., 60., 20.]);
    h.press(corner);
    h.drag(corner + egui::vec2(17.3, 9.1));
    h.assert_valid_preview("perspective drag");
    let before = h.box_shape();
    // Far above the horizon the pointer ray misses the drag plane.
    h.drag(Pos2::new(corner.x, -5.0e5));
    let feedback = h.controller.primitive_drag_feedback().unwrap().to_owned();
    assert!(
        feedback.contains("keeping the last valid size"),
        "{feedback}"
    );
    assert_eq!(h.box_shape(), before);
    h.assert_valid_preview("off plane");
    h.drag(Pos2::new(corner.x + 3.0e5, 4.0e5));
    h.assert_valid_preview("far below");
    h.drag(corner + egui::vec2(17.3, 9.1));
    assert_eq!(h.box_shape(), before);
    h.release(corner + egui::vec2(17.3, 9.1));
    h.assert_valid_preview("released");
}

#[test]
fn a_dragged_primitive_commits_as_one_undo_step_with_the_shown_dimensions() {
    let mut h = PrimitiveHarness::new(top_camera(3.713));
    let corner = h.screen([100., 60., 20.]);
    h.press(corner);
    let target = corner + egui::vec2(51.7, -33.3);
    h.drag(target);
    h.release(target);
    h.assert_valid_preview("released");
    let (size, origin) = h.box_shape();
    let revision = h.session.revision();
    let parts = h
        .controller
        .apply(&mut h.session, &Default::default())
        .unwrap();
    assert_eq!(parts.len(), 1);
    assert_eq!(h.session.revision(), revision + 1);
    let stored = model::read_primitive(h.session.document(), parts[0].feature).unwrap();
    let ResolvedPrimitive::Box {
        size_mm,
        frame: stored_frame,
    } = stored.parameters.resolve(h.session.parameters()).unwrap()
    else {
        panic!("box")
    };
    assert_eq!((size_mm, stored_frame.origin), (size, origin));
    h.session
        .undo(h.session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert!(h.session.document().entities().next().is_none());
}
