//! Form routing checks; a pending scripted lane is never treated as native proof.
use super::*;
use multicad_core::{
    document::{Document, EntityData},
    model,
};
use std::sync::Mutex;

#[derive(Default)]
struct Capture(Arc<Mutex<Option<Document>>>);
impl SceneBackend for Capture {
    fn request(&mut self, document: &Document, _: SceneStamp) -> Result<(), String> {
        *self.0.lock().unwrap() = Some(document.clone());
        Ok(())
    }
    fn poll(&mut self) -> Option<crate::viewport::SceneReply> {
        None
    }
}
/// The Wall tool draws (CAD-427): click `corners` on the plan at `z`, then
/// Enter ends the chain, keeping the path for Preview and Apply.
fn draw_wall(
    controller: &mut ModelToolController,
    session: &multicad_core::commands::Session,
    z: f32,
    corners: &[[f32; 2]],
) {
    let ctx = crate::egui::Context::default();
    let rect =
        crate::egui::Rect::from_min_size(crate::egui::Pos2::ZERO, crate::egui::vec2(800., 600.));
    let camera = crate::camera::Camera {
        roll: 0.0,
        yaw: 0.,
        pitch: std::f32::consts::FRAC_PI_2,
        target: [2000., 0., z],
        scale: 0.1,
        perspective: false,
    };
    for corner in corners {
        let mut output = ctx.run_ui(Default::default(), |ui| {
            let p = camera.project([corner[0], corner[1], z], rect);
            assert!(controller.outline_viewport(
                session,
                super::super::TerminationViewport {
                    scene: None,
                    report: None,
                    camera,
                    rect,
                    pointer: Some(crate::egui::pos2(p[0], p[1])),
                    clicked: true,
                    painter: ui.painter(),
                    section: None,
                }
            ));
        });
        output.textures_delta.clear();
    }
    assert_eq!(
        controller.outline_key(super::super::OutlineKey::Finish),
        Some(true)
    );
}
fn controller() -> (ModelToolController, Arc<Mutex<Option<Document>>>) {
    let capture = Capture::default();
    let sent = capture.0.clone();
    let mut controller = ModelToolController::default();
    controller.set_backend(Box::new(capture));
    (controller, sent)
}

#[test]
fn repeated_equipment_keeps_its_chosen_floor_when_the_active_level_changes() {
    use crate::egui::{self, Event, Pos2, Rect, Vec2};
    use multicad_core::equipment;
    let (mut session, levels, _, _) = crate::architecture_view::tests::building();
    let resources = ProjectResources::default();
    let mut controller = ModelToolController::default();
    controller.set_backend(Box::new(super::super::native_fixture::NativeBackend::new()));
    controller.set_architecture_level(Some(levels[0]));
    controller
        .begin_equipment_preset(
            &session,
            &resources,
            Channel::Ui,
            equipment::Kind::Refrigerator,
        )
        .unwrap();
    let ctx = egui::Context::default();
    let mut paint = |controller: &mut ModelToolController, events: Vec<Event>| {
        ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(600., 2400.))),
                events,
                ..Default::default()
            },
            |ui| {
                controller.show(ui, &mut session, &resources).unwrap();
            },
        )
    };
    paint(&mut controller, vec![]).drop_without_applying_deltas();
    let output = paint(&mut controller, vec![]);
    let at = output
        .shapes
        .iter()
        .find_map(|shape| match &shape.shape {
            egui::Shape::Text(text) if text.galley.text() == "Keep placing this type" => {
                Some(text.pos + text.galley.size() * 0.5)
            }
            _ => None,
        })
        .expect("repeat checkbox is painted");
    output.drop_without_applying_deltas();
    for pressed in [true, false] {
        paint(
            &mut controller,
            vec![
                Event::PointerMoved(at),
                Event::PointerButton {
                    pos: at,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ],
        )
        .drop_without_applying_deltas();
    }
    for (index, x) in [400., 800.].into_iter().enumerate() {
        assert!(
            controller
                .preview(&session, &resources, Channel::Ui)
                .is_err(),
            "every repeated occurrence requires a fresh viewport pick"
        );
        ctx.run_ui(Default::default(), |ui| {
            let camera = crate::camera::Camera {
                roll: 0.0,
                yaw: 0.,
                pitch: std::f32::consts::FRAC_PI_2,
                target: [x, 400., 0.],
                scale: 0.1,
                perspective: false,
            };
            let rect = Rect::from_min_size(Pos2::ZERO, Vec2::new(800., 600.));
            assert!(controller.outline_viewport(
                &session,
                super::super::TerminationViewport {
                    scene: None,
                    report: None,
                    camera,
                    rect,
                    pointer: Some(rect.center()),
                    clicked: true,
                    painter: ui.painter(),
                    section: None,
                }
            ));
        })
        .drop_without_applying_deltas();
        controller
            .preview(&session, &resources, Channel::Ui)
            .unwrap();
        let until = std::time::Instant::now() + std::time::Duration::from_secs(10);
        while !controller.can_apply(&session, &resources) && std::time::Instant::now() < until {
            controller.poll(&session, &resources);
            std::thread::sleep(std::time::Duration::from_millis(2));
        }
        assert!(controller.can_apply(&session, &resources));
        controller.set_architecture_level(Some(levels[1]));
        controller.apply(&mut session, &resources).unwrap();
        let attachments: Vec<_> = session
            .document()
            .entities()
            .filter_map(|entity| equipment::read_attachment(session.document(), entity.id).ok())
            .collect();
        assert_eq!(attachments.len(), index + 1);
        assert!(
            attachments
                .iter()
                .all(|attachment| attachment.floor == levels[0]),
            "repetition preserves the chosen floor instead of reseeding from the active level"
        );
    }
    session
        .undo(session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(
        session
            .document()
            .entities()
            .filter(|entity| equipment::read_attachment(session.document(), entity.id).is_ok())
            .count(),
        1,
        "the repeated occurrence has a separate undo step"
    );
}
#[test]
fn active_level_seeds_new_wall_stair_slab_and_equipment_but_never_reassigns_saved_features() {
    let (session, levels, _, features) = crate::architecture_view::tests::building();
    let resources = ProjectResources::default();
    let before = session.document().clone();
    let (mut controller, sent) = controller();
    controller.set_architecture_level(Some(levels[1]));
    controller
        .begin_wall(&session, &resources, Channel::Ui)
        .unwrap();
    draw_wall(&mut controller, &session, 3000., &[[0., 0.], [4000., 0.]]);
    controller
        .preview(&session, &resources, Channel::Ui)
        .unwrap();
    let wall_document = sent.lock().unwrap().clone().unwrap();
    let wall = model::walls(&wall_document).next().unwrap();
    assert_eq!(
        model::read_wall(&wall_document, wall.id).unwrap().level,
        Some(levels[1])
    );
    assert!(
        !controller.can_apply(&session, &resources),
        "a captured request is not a native receipt"
    );
    controller.cancel();
    controller
        .begin_stairs(&session, &resources, Channel::Ui)
        .unwrap();
    let ctx = crate::egui::Context::default();
    for point in [[1000., 500., 3000.], [2000., 500., 3000.]] {
        let mut output = ctx.run_ui(Default::default(), |ui| {
            let camera = crate::camera::Camera {
                roll: 0.0,
                yaw: 0.,
                pitch: std::f32::consts::FRAC_PI_2,
                target: [0., 0., 3000.],
                scale: 0.1,
                perspective: false,
            };
            let rect = crate::egui::Rect::from_min_size(
                crate::egui::Pos2::ZERO,
                crate::egui::vec2(800., 600.),
            );
            let p = camera.project(point, rect);
            assert!(controller.outline_viewport(
                &session,
                super::super::TerminationViewport {
                    scene: None,
                    report: None,
                    camera,
                    rect,
                    pointer: Some(crate::egui::pos2(p[0], p[1])),
                    clicked: true,
                    painter: ui.painter(),
                    section: None,
                }
            ));
        });
        output.textures_delta.clear();
    }
    let Some(ModelForm::Architecture(panel)) = &controller.form else {
        panic!("stair form")
    };
    let super::super::architecture::ArchitectureCommand::Stair { stair, .. } =
        panel.request(session.document(), true).unwrap()
    else {
        panic!("stair request")
    };
    assert_eq!(stair.level, levels[1]);
    controller.cancel();
    controller
        .begin_architecture(&session, &resources, Channel::Ui)
        .unwrap();
    controller
        .preview(&session, &resources, Channel::Ui)
        .unwrap();
    let slab_document = sent.lock().unwrap().clone().unwrap();
    let slab = model::slabs(&slab_document)
        .find(|slab| session.document().get(slab.id).is_none())
        .unwrap();
    assert_eq!(
        model::read_slab(&slab_document, slab.id).unwrap().level,
        levels[1]
    );
    controller.cancel();
    assert!(ModelToolController::editable(
        session.document(),
        features[0]
    ));
    controller
        .begin_existing(&session, &resources, features[0], Channel::Ui)
        .unwrap();
    controller
        .preview(&session, &resources, Channel::Ui)
        .unwrap();
    assert_eq!(
        model::read_slab(sent.lock().unwrap().as_ref().unwrap(), features[0])
            .unwrap()
            .level,
        levels[0]
    );
    controller.cancel();
    controller
        .begin_equipment_preset(
            &session,
            &resources,
            Channel::Ui,
            multicad_core::equipment::Kind::Refrigerator,
        )
        .unwrap();
    // The existing Inspector must still require a graphical placement before Preview.
    assert!(
        controller
            .preview(&session, &resources, Channel::Ui)
            .is_err()
    );
    let ctx = crate::egui::Context::default();
    let mut output = ctx.run_ui(Default::default(), |ui| {
        let rect = crate::egui::Rect::from_min_size(
            crate::egui::Pos2::ZERO,
            crate::egui::vec2(800., 600.),
        );
        let camera = crate::camera::Camera {
            roll: 0.0,
            yaw: 0.,
            pitch: std::f32::consts::FRAC_PI_2,
            target: [0., 0., 3000.],
            scale: 0.1,
            perspective: false,
        };
        assert!(controller.outline_viewport(
            &session,
            super::super::TerminationViewport {
                scene: None,
                report: None,
                camera,
                rect,
                pointer: Some(rect.center()),
                clicked: true,
                painter: ui.painter(),
                section: None,
            }
        ));
    });
    output.textures_delta.clear();
    controller
        .preview(&session, &resources, Channel::Ui)
        .unwrap();
    let equipment_document = sent.lock().unwrap().clone().unwrap();
    let attachment = equipment_document.entities().find_map(|entity|
        matches!(&entity.data,EntityData::Extension{definition,..} if definition.type_id == multicad_core::equipment::ATTACHMENT_TYPE)
            .then(||multicad_core::equipment::read_attachment(&equipment_document,entity.id).unwrap())).unwrap();
    assert_eq!(attachment.floor, levels[1]);
    multicad_model::Compiler::default()
        .compile(&equipment_document, session.revision())
        .unwrap();
    assert_eq!(session.document(), &before);
}

#[test]
fn native_candidate_level_filter_uses_prepared_membership_and_keeps_canonical_context() {
    let (mut session, levels, bodies, _) = crate::architecture_view::tests::building();
    let source = session.document().clone();
    let resources = ProjectResources::default();
    let mut controller = ModelToolController::default();
    controller.set_backend(Box::new(super::super::native_fixture::NativeBackend::new()));
    controller.set_architecture_level(Some(levels[1]));
    controller
        .begin_wall(&session, &resources, Channel::Ui)
        .unwrap();
    draw_wall(&mut controller, &session, 3000., &[[0., 0.], [4000., 0.]]);
    controller
        .preview(&session, &resources, Channel::Ui)
        .unwrap();
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
    while !controller.can_apply(&session, &resources) && std::time::Instant::now() < deadline {
        controller.poll(&session, &resources);
        std::thread::sleep(std::time::Duration::from_millis(2));
    }
    assert!(controller.can_apply(&session, &resources));
    let native = controller.display_scene(&session, &resources).unwrap();
    let wall = controller.candidate.as_ref().unwrap().outputs[0].body;
    let prepared = controller
        .architecture_preview_levels(&session, &native)
        .unwrap();
    assert_eq!(prepared.level_for(wall), Some(Some(levels[1])));
    let mut context = crate::architecture_view::Context::default();
    context.sync(session.document(), session.revision());
    let mut camera = crate::camera::Camera::default();
    context.set_enabled(true, &mut camera);
    context.set_plan(true, &mut camera);
    let ground = context.preview_scene(&native, prepared).unwrap();
    assert_eq!(
        ground
            .instances()
            .iter()
            .map(|instance| instance.id)
            .collect::<Vec<_>>(),
        vec![bodies[0]]
    );
    assert!(Arc::ptr_eq(
        &ground,
        &context.preview_scene(&native, prepared).unwrap()
    ));
    context.recall.level = Some(levels[1]);
    let upper = context.preview_scene(&native, prepared).unwrap();
    let visible: std::collections::BTreeSet<_> = upper
        .instances()
        .iter()
        .map(|instance| instance.id)
        .collect();
    assert_eq!(visible, std::collections::BTreeSet::from([bodies[1], wall]));
    assert!(
        context.levels.as_ref().unwrap().level_for(wall).is_none(),
        "preview facts must not replace canonical level ownership"
    );
    for (id, geometry) in native.geometries() {
        assert!(Arc::ptr_eq(geometry, &upper.geometries()[id]));
    }
    assert_eq!(session.document(), &source);
    controller.apply(&mut session, &resources).unwrap();
    context.sync(session.document(), session.revision());
    let handoff = controller.display_scene(&session, &resources).unwrap();
    let levels = controller
        .architecture_preview_levels(&session, &handoff)
        .unwrap();
    let view = context.preview_scene(&handoff, levels).unwrap();
    assert_eq!(
        view.instances()
            .iter()
            .map(|instance| instance.id)
            .collect::<std::collections::BTreeSet<_>>(),
        visible
    );
    session
        .undo(session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(session.document(), &source);
}

#[test]
fn path_keys_belong_to_the_wall_drawing_only_while_it_captures_the_plan() {
    use crate::egui;
    let session =
        multicad_core::commands::Session::new(Document::new("Wall keys").unwrap()).unwrap();
    let resources = ProjectResources::default();
    let (mut controller, _) = controller();
    controller
        .begin_wall(&session, &resources, Channel::Ui)
        .unwrap();
    // The Wall tool draws at once (CAD-427), so the path keys are the
    // drawing's from the start; nothing is placed yet to remove.
    assert_eq!(
        controller.outline_key(super::super::OutlineKey::RemoveLast),
        Some(false)
    );
    let ctx = egui::Context::default();
    let rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(800., 600.));
    let camera = crate::camera::Camera {
        roll: 0.0,
        yaw: 0.,
        pitch: 1.5,
        scale: 0.1,
        target: [0., 2000., 0.],
        perspective: false,
    };
    let click = |controller: &mut ModelToolController, at: [f32; 2]| {
        let mut output = ctx.run_ui(Default::default(), |ui| {
            let s = camera.project([at[0], at[1], 0.], rect);
            controller.outline_viewport(
                &session,
                super::super::TerminationViewport {
                    scene: None,
                    report: None,
                    camera,
                    rect,
                    pointer: Some(egui::pos2(s[0], s[1])),
                    clicked: true,
                    painter: ui.painter(),
                    section: None,
                },
            );
        });
        output.textures_delta.clear();
    };
    // The first click sets where the wall starts.
    click(&mut controller, [-1000., 3000.]);
    click(&mut controller, [1000., 3000.]);
    click(&mut controller, [1000., 4000.]);
    assert_eq!(
        controller.outline_key(super::super::OutlineKey::RemoveLast),
        Some(true)
    );
    assert_eq!(
        controller.outline_key(super::super::OutlineKey::Finish),
        Some(true)
    );
    assert!(controller.active(), "Enter finished the path, not the form");
    let Some(ModelRequest::Wall(command)) = controller.request_for_tests(&session) else {
        panic!("wall request")
    };
    assert_eq!(
        command.wall.path,
        vec![
            [model::length(-1000.), model::length(3000.)],
            [model::length(1000.), model::length(3000.)]
        ]
    );
    // After the first drawing an empty click selects nothing and Escape is
    // the form's again.
    click(&mut controller, [5000., 5000.]);
    let Some(ModelRequest::Wall(drawn)) = controller.request_for_tests(&session) else {
        panic!("wall request")
    };
    assert_eq!(
        drawn.wall.path, command.wall.path,
        "no click lands while idle"
    );
    assert_eq!(
        controller.outline_key(super::super::OutlineKey::Cancel),
        None
    );
    assert!(controller.active());
}

#[test]
fn wall_tool_starts_empty_and_every_end_gesture_keeps_the_chain_as_one_undoable_wall() {
    use crate::egui;
    let mut session =
        multicad_core::commands::Session::new(Document::new("Wall chain").unwrap()).unwrap();
    let source = session.document().clone();
    let resources = ProjectResources::default();
    let mut controller = ModelToolController::default();
    controller.set_backend(Box::new(super::super::native_fixture::NativeBackend::new()));
    let ctx = egui::Context::default();
    let rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(800., 600.));
    let camera = crate::camera::Camera {
        roll: 0.0,
        yaw: 0.,
        pitch: std::f32::consts::FRAC_PI_2,
        scale: 0.1,
        target: [2000., 1500., 0.],
        perspective: false,
    };
    let frame = |controller: &mut ModelToolController,
                 session: &multicad_core::commands::Session,
                 at: [f32; 2],
                 clicked: bool,
                 events: Vec<egui::Event>| {
        let s = camera.project([at[0], at[1], 0.], rect);
        let mut output = ctx.run_ui(
            egui::RawInput {
                events,
                ..Default::default()
            },
            |ui| {
                controller.outline_viewport(
                    session,
                    super::super::TerminationViewport {
                        scene: None,
                        report: None,
                        camera,
                        rect,
                        pointer: Some(egui::pos2(s[0], s[1])),
                        clicked,
                        painter: ui.painter(),
                        section: None,
                    },
                );
            },
        );
        output.textures_delta.clear();
    };
    let corners = [[0., 0.], [4000., 0.], [4000., 3000.]];
    let expected: Vec<_> = corners
        .iter()
        .map(|[x, y]| [model::length(*x as f64), model::length(*y as f64)])
        .collect();
    for gesture in ["Enter", "Escape", "right-click"] {
        controller
            .begin_wall(&session, &resources, Channel::Ui)
            .unwrap();
        // Choosing Wall places nothing: the path is empty until clicked.
        match controller.request_for_tests(&session) {
            Some(ModelRequest::Wall(command)) => {
                assert!(command.wall.path.is_empty(), "{gesture}")
            }
            None => {}
            Some(_) => panic!("wall request ({gesture})"),
        }
        assert_eq!(
            controller.outline_key(super::super::OutlineKey::RemoveLast),
            Some(false),
            "drawing from the start ({gesture})"
        );
        for at in corners {
            frame(&mut controller, &session, at, true, vec![]);
        }
        match gesture {
            "Enter" => assert_eq!(
                controller.outline_key(super::super::OutlineKey::Finish),
                Some(true)
            ),
            "Escape" => assert_eq!(
                controller.outline_key(super::super::OutlineKey::Cancel),
                Some(true)
            ),
            _ => {
                let s = camera.project([2000., 2000., 0.], rect);
                let pos = egui::pos2(s[0], s[1]);
                for pressed in [true, false] {
                    frame(
                        &mut controller,
                        &session,
                        [2000., 2000.],
                        false,
                        vec![
                            egui::Event::PointerMoved(pos),
                            egui::Event::PointerButton {
                                pos,
                                button: egui::PointerButton::Secondary,
                                pressed,
                                modifiers: egui::Modifiers::NONE,
                            },
                        ],
                    );
                }
            }
        }
        // The chain ended: the drawing stopped, the form stays open, every
        // placed segment is kept and the document is untouched.
        assert_eq!(
            controller.outline_key(super::super::OutlineKey::RemoveLast),
            None,
            "{gesture}"
        );
        assert!(controller.active(), "{gesture}");
        let Some(ModelRequest::Wall(command)) = controller.request_for_tests(&session) else {
            panic!("wall request ({gesture})")
        };
        assert_eq!(command.wall.path, expected, "{gesture}");
        assert!(!command.wall.closed, "{gesture}");
        assert_eq!(session.document(), &source, "{gesture}");
        if gesture != "right-click" {
            controller.cancel();
            assert_eq!(session.document(), &source);
        }
    }
    // Apply builds the whole chain as one wall feature; one Undo removes
    // it and one Redo restores it.
    controller
        .preview(&session, &resources, Channel::Ui)
        .unwrap();
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
    while !controller.can_apply(&session, &resources) && std::time::Instant::now() < deadline {
        controller.poll(&session, &resources);
        std::thread::sleep(std::time::Duration::from_millis(2));
    }
    let applied = controller.apply(&mut session, &resources).unwrap();
    assert_eq!(applied.len(), 1);
    let after = session.document().clone();
    assert_eq!(model::walls(&after).count(), 1);
    assert_eq!(
        model::read_wall(&after, applied[0].feature).unwrap().path,
        expected
    );
    session
        .undo(session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(session.document(), &source, "one Undo removes the chain");
    session
        .redo(session.revision(), Actor::local(Channel::Ui))
        .unwrap();
    assert_eq!(session.document(), &after, "one Redo restores it");
}
