//! CAD-375 point-to-point placement of a new panel through the production
//! viewport pick and native preview lane.
use super::*;
use crate::{
    camera::Camera,
    egui::{self, Context, Pos2, Rect},
    modeling::{TerminationViewport, cabinet::PanelPointPick, native_fixture},
};
use multicad_core::{
    model::{PanelSource, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest},
    project_file::{self, ProjectFile},
};
use multicad_model::Compiler;
use multicad_precision::PortableMath as _;
use std::time::{Duration, Instant};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
/// Panel stock plus a 100 × 60 × 40 box rotated 30° about Z at (1000, 0, 0).
fn fixture() -> (Session, [f64; 3]) {
    let mut session = super::super::cabinet::tests::session();
    let (c, s) = (
        30_f64.to_radians().portable_cos(),
        30_f64.to_radians().portable_sin(),
    );
    let other = PrimitiveRequest::new(
        "Rotated block",
        PrimitiveParameters::Box {
            size: [length(100.), length(60.), length(40.)],
            frame: PrimitiveFrame {
                origin: [length(1000.), length(0.), length(0.)],
                u: [c, s, 0.],
                v: [-s, c, 0.],
            },
        },
    );
    let prepared = model::prepare_primitive(&session, session.revision(), actor(), &other).unwrap();
    session.commit(prepared.prepared).unwrap();
    (session, [1000. + 100. * c, 100. * s, 40.])
}
fn canonical(session: &Session) -> Arc<Scene> {
    let model = Compiler::default()
        .compile(session.document(), session.revision())
        .unwrap();
    native_fixture::scene(session.document(), model.snapshot(), 0).unwrap()
}
fn wait(controller: &mut ModelToolController, session: &Session) {
    let until = Instant::now() + Duration::from_secs(5);
    while controller.pending() && Instant::now() < until {
        controller.poll(session, &Default::default());
        std::thread::sleep(Duration::from_millis(1));
    }
    assert!(!controller.pending(), "native worker did not reply");
}
fn camera() -> Camera {
    Camera {
        yaw: -0.6,
        pitch: 0.6,
        scale: 0.6,
        target: [600., 150., 0.],
        perspective: false,
    }
}
fn rect() -> Rect {
    Rect::from_min_max(Pos2::new(0., 0.), Pos2::new(1200., 800.))
}
fn screen(p: [f64; 3]) -> Pos2 {
    let [x, y, _] = camera().project(p.map(|v| v as f32), rect());
    Pos2::new(x, y)
}
fn click(
    c: &mut ModelToolController,
    s: &Session,
    scene: Option<&Scene>,
    pointer: Pos2,
) -> Result<Option<String>, String> {
    let ctx = Context::default();
    let mut result = None;
    ctx.run_ui(egui::RawInput::default(), |ui| {
        result = Some(c.panel_point_viewport(
            s,
            &Default::default(),
            TerminationViewport {
                scene,
                report: None,
                camera: camera(),
                rect: rect(),
                pointer: Some(pointer),
                clicked: true,
                painter: ui.painter(),
                section: None,
            },
        ));
    })
    .drop_without_applying_deltas();
    result.unwrap()
}
fn form(c: &mut ModelToolController) -> &mut crate::modeling::cabinet::CabinetPanel {
    let Some(ModelForm::Cabinet(form)) = &mut c.form else {
        panic!("panel form")
    };
    form
}
fn origin(
    document: &Document,
    feature: EntityId,
    parameters: &multicad_core::parameters::EvaluatedParameters,
) -> [f64; 3] {
    let panel = model::read_panel(document, feature).unwrap();
    let PanelSource::Rectangle { origin, .. } = panel.source else {
        panic!("rectangle panel")
    };
    origin.map(|d| {
        d.value_in(
            multicad_core::parameters::Unit::Millimeter,
            Some(parameters),
        )
        .unwrap()
    })
}
fn close(actual: [f64; 3], expected: [f64; 3]) {
    for (a, b) in actual.into_iter().zip(expected) {
        assert!((a - b).abs() < 1e-6, "{actual:?} != {expected:?}");
    }
}

#[test]
fn panel_corner_snaps_to_a_rotated_part_and_survives_apply_undo_redo_and_reopen() {
    let (mut session, target) = fixture();
    let original = session.document().clone();
    let revision = session.revision();
    let scene = canonical(&session);
    let mut c = ModelToolController::default();
    c.set_backend(Box::new(native_fixture::NativeBackend::new()));
    c.begin_panel(&session, &Default::default(), Channel::Ui)
        .unwrap();
    form(&mut c).toggle_point_pick();
    assert!(c.panel_point_picking());
    assert!(!c.arm_confirmation(&session, &Default::default()));

    // The draft's top, low-v, high-u corner is the source.
    let source = [800., 0., 18.];
    let message = click(&mut c, &session, Some(&scene), screen(source)).unwrap();
    assert!(message.unwrap().contains("panel corner"));
    assert!(matches!(
        form(&mut c).point_pick(),
        Some(PanelPointPick::Target { .. })
    ));

    // Invalid targets keep the last valid draft: empty space and a stale or
    // missing model. Apply is refused while a pick is pending.
    let staged: Vec<_> = form(&mut c).origin_texts();
    assert!(click(&mut c, &session, Some(&scene), Pos2::new(3., 3.)).is_err());
    assert!(click(&mut c, &session, None, screen(target)).is_err());
    assert_eq!(form(&mut c).origin_texts(), staged);
    assert!(c.apply(&mut session, &Default::default()).is_err());
    assert_eq!(session.document(), &original);

    let message = click(&mut c, &session, Some(&scene), screen(target)).unwrap();
    assert!(message.unwrap().contains("vertex"));
    assert!(!c.panel_point_picking());
    let corners = c.primitive_envelope(&session).unwrap();
    // The chosen source coincides with the rotated part's corner; size,
    // plane and stock thickness are unchanged.
    close(corners[5], target);
    close(
        std::array::from_fn(|i| corners[7][i] - corners[0][i]),
        [800., 400., 18.],
    );
    assert_eq!(session.document(), &original);

    wait(&mut c, &session);
    assert!(
        c.can_apply(&session, &Default::default()),
        "{:?}",
        c.error()
    );
    let created = c.apply(&mut session, &Default::default()).unwrap();
    assert_eq!(session.revision(), revision + 1);
    let feature = created[0].feature;
    let expected: [f64; 3] = std::array::from_fn(|i| target[i] - source[i]);
    close(
        origin(session.document(), feature, session.parameters()),
        expected,
    );
    let native = Compiler::default()
        .compile(session.document(), session.revision())
        .unwrap();
    let bounds = native.snapshot().shapes[&feature].bounds;
    close(bounds.min, expected);
    close(
        bounds.max,
        std::array::from_fn(|i| expected[i] + [800., 400., 18.][i]),
    );

    let applied = session.document().clone();
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &original);
    session.redo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &applied);
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(applied).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let reopened = Session::new(opened.project.document).unwrap();
    close(
        origin(reopened.document(), feature, reopened.parameters()),
        expected,
    );
}

#[test]
fn other_panel_source_points_and_cancel_leave_the_document_untouched() {
    let (session, target) = fixture();
    let original = session.document().clone();
    let scene = canonical(&session);
    let mut c = ModelToolController::default();
    c.set_backend(Box::new(native_fixture::NativeBackend::new()));
    c.begin_panel(&session, &Default::default(), Channel::Ui)
        .unwrap();
    // Clicking a draft corner before Choose Points does nothing here.
    assert_eq!(
        click(&mut c, &session, Some(&scene), screen([0.; 3])),
        Ok(None)
    );
    form(&mut c).toggle_point_pick();
    // An edge midpoint is a valid source too.
    let midpoint = [400., 400., 18.];
    let message = click(&mut c, &session, Some(&scene), screen(midpoint)).unwrap();
    assert!(message.unwrap().contains("edge midpoint"));
    // Cancel Points drops the pending source without moving the draft.
    let before = form(&mut c).origin_texts();
    form(&mut c).toggle_point_pick();
    assert!(!c.panel_point_picking());
    assert_eq!(form(&mut c).origin_texts(), before);
    // The world origin is always a target.
    form(&mut c).toggle_point_pick();
    click(&mut c, &session, Some(&scene), screen(midpoint)).unwrap();
    let message = click(&mut c, &session, Some(&scene), screen([0.; 3])).unwrap();
    assert!(message.unwrap().contains("world origin"));
    let corners = c.primitive_envelope(&session).unwrap();
    close(corners[0], [-400., -400., -18.]);
    close(corners[7], [400., 0., 0.]);
    // A later pick from another corner follows the draft where it is now.
    form(&mut c).toggle_point_pick();
    click(&mut c, &session, Some(&scene), screen(corners[7])).unwrap();
    click(&mut c, &session, Some(&scene), screen(target)).unwrap();
    close(c.primitive_envelope(&session).unwrap()[7], target);
    c.cancel();
    assert_eq!(session.document(), &original);
}
