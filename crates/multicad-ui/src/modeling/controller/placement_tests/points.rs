//! CAD-391 point-to-point Move through the production viewport and native lane.
use super::*;
use crate::{
    camera::Camera,
    egui::{self, Context, Pos2, Rect},
    modeling::{
        TerminationViewport,
        placement::points::{PointPick, moving_instances, snap_points},
    },
};
use multicad_core::{
    model::PrimitiveFrame,
    project_file::{self, ProjectFile},
};
use multicad_precision::PortableMath as _;

struct Fixture {
    session: Session,
    moving: PrimitiveRequest,
    other: PrimitiveRequest,
}
/// A 20 × 10 × 5 box at the origin and a 30 × 18 × 4 panel rotated 30°
/// about Z with its frame origin at (50, 0, 0).
fn fixture() -> Fixture {
    let mut session = Session::new(Document::new("Point to point").unwrap()).unwrap();
    let moving = PrimitiveRequest::new(
        "Moving",
        PrimitiveParameters::Box {
            size: [length(20.), length(10.), length(5.)],
            frame: Default::default(),
        },
    );
    let (c, s) = (
        30_f64.to_radians().portable_cos(),
        30_f64.to_radians().portable_sin(),
    );
    let other = PrimitiveRequest::new(
        "Rotated panel",
        PrimitiveParameters::Box {
            size: [length(30.), length(18.), length(4.)],
            frame: PrimitiveFrame {
                origin: [length(50.), length(0.), length(0.)],
                u: [c, s, 0.],
                v: [-s, c, 0.],
            },
        },
    );
    for request in [&moving, &other] {
        let prepared =
            model::prepare_primitive(&session, session.revision(), actor(), request).unwrap();
        session.commit(prepared.prepared).unwrap();
    }
    Fixture {
        session,
        moving,
        other,
    }
}
fn camera() -> Camera {
    Camera {
        yaw: -0.6,
        pitch: 0.6,
        scale: 6.,
        target: [25., 5., 0.],
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
    s: &mut Session,
    scene: Option<&Scene>,
    pointer: Pos2,
) -> Result<Option<String>, String> {
    let ctx = Context::default();
    let mut result = None;
    ctx.run_ui(egui::RawInput::default(), |ui| {
        result = Some(c.placement_point_viewport(
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
fn panel(c: &mut ModelToolController) -> &mut crate::modeling::placement::PlacementPanel {
    let Some(ModelForm::Placement(panel)) = &mut c.form else {
        panic!()
    };
    panel
}
fn close(actual: [f64; 3], expected: [f64; 3]) {
    for (a, b) in actual.into_iter().zip(expected) {
        assert!((a - b).abs() < 1e-6, "{actual:?} != {expected:?}");
    }
}

#[test]
fn moving_part_snap_points_follow_the_staged_pose_and_targets_exclude_moving_geometry() {
    let f = fixture();
    let scene = canonical(&f.session);
    let selected = BTreeSet::from([f.moving.body]);
    let moving = moving_instances(f.session.document(), &scene, &selected);
    assert_eq!(moving, selected);
    let sources = snap_points(
        &scene,
        |id| moving.contains(&id),
        |p| [p[0] + 1., p[1], p[2]],
    );
    assert!(sources.iter().all(|p| p.owner == Some(f.moving.body)));
    // Every box corner is offered, shown at the staged position.
    for corner in [[0., 0., 0.], [20., 10., 5.], [20., 0., 5.]] {
        assert!(
            sources.iter().any(|p| p.canonical == corner
                && p.shown == [corner[0] + 1., corner[1], corner[2]]
                && p.label == "Vertex"),
            "{corner:?} missing from {sources:?}"
        );
    }
    assert!(sources.iter().any(|p| p.label == "Edge midpoint"));
    let targets = snap_points(&scene, |id| !moving.contains(&id), |p| p);
    assert!(!targets.is_empty());
    assert!(targets.iter().all(|p| p.owner == Some(f.other.body)));
}

#[test]
fn point_to_point_move_aligns_a_rotated_stage_to_a_rotated_panel_and_commits_once() {
    let mut f = fixture();
    let original = f.session.document().clone();
    let revision = f.session.revision();
    let scene = canonical(&f.session);
    let selected = BTreeSet::from([f.moving.body]);
    let mut c = ModelToolController::default();
    c.set_backend(Box::new(native_fixture::NativeBackend::new()));
    c.begin_placement(
        &f.session,
        &Default::default(),
        &scene,
        &selected,
        PlacementMode::Move,
        Channel::Ui,
    )
    .unwrap();
    // A staged 90° turn about Z through the world origin before choosing points.
    let staged = NumericTransform {
        angle: model::angle(90.),
        ..Default::default()
    };
    panel(&mut c).set_drag_parameters(&staged).unwrap();
    panel(&mut c).toggle_point_pick();
    assert!(c.placement_point_picking());
    // Enter cannot confirm while points are being chosen.
    assert!(!c.arm_confirmation(&f.session, &Default::default()));

    // The source corner (20, 10, 5) is shown at its staged position (-10, 20, 5).
    let source = [20., 10., 5.];
    let shown = [-10., 20., 5.];
    let message = click(&mut c, &mut f.session, Some(&scene), screen(shown)).unwrap();
    assert!(message.unwrap().contains("vertex"));
    assert_eq!(
        panel(&mut c).point_pick(),
        Some(PointPick::Target {
            source,
            label: "Vertex"
        })
    );

    // Invalid targets keep the last valid state: empty space, the moving
    // part itself, and a stale or missing model.
    let staged_before = panel(&mut c).parameters(f.session.document()).unwrap();
    assert!(click(&mut c, &mut f.session, Some(&scene), Pos2::new(3., 3.)).is_err());
    assert!(click(&mut c, &mut f.session, Some(&scene), screen([-10., 0., 5.])).is_err());
    assert!(click(&mut c, &mut f.session, None, screen(shown)).is_err());
    assert_eq!(
        panel(&mut c).parameters(f.session.document()).unwrap(),
        staged_before
    );
    assert!(matches!(
        panel(&mut c).point_pick(),
        Some(PointPick::Target { .. })
    ));

    // The far bottom corner of the rotated panel along its u direction.
    let (cos, sin) = (
        30_f64.to_radians().portable_cos(),
        30_f64.to_radians().portable_sin(),
    );
    let target = [50. + 30. * cos, 30. * sin, 0.];
    let message = click(&mut c, &mut f.session, Some(&scene), screen(target)).unwrap();
    assert!(message.unwrap().contains("staged"));
    assert!(!c.placement_point_picking());
    let aligned = panel(&mut c).parameters(f.session.document()).unwrap();
    let pose = aligned.resolve(f.session.parameters()).unwrap();
    close(pose.point(source), target);
    assert_eq!(aligned.angle, staged.angle);
    assert_eq!(aligned.pivot, staged.pivot);
    let alignment = panel(&mut c).alignment().unwrap();
    assert_eq!((alignment.source, alignment.target), (source, target));
    // Choosing points never replaced the moving-part selection.
    assert_eq!(c.placement_targets(), Some(&selected));
    assert_eq!(f.session.document(), &original);
    assert_eq!(f.session.revision(), revision);

    wait(&mut c, &f.session);
    assert!(
        c.can_apply(&f.session, &Default::default()),
        "{:?}",
        c.error()
    );
    let parts = c.apply(&mut f.session, &Default::default()).unwrap();
    assert_eq!(f.session.revision(), revision + 1);
    let delta: [f64; 3] = std::array::from_fn(|i| target[i] - shown[i]);
    let check = |document: &Document, revision| {
        let native = Compiler::default().compile(document, revision).unwrap();
        let bounds = native.snapshot().shapes[&parts[0].feature].bounds;
        // Rotated 90°: x ∈ [-10, 0], y ∈ [0, 20], z ∈ [0, 5]; then translated.
        close(
            bounds.min,
            std::array::from_fn(|i| [-10., 0., 0.][i] + delta[i]),
        );
        close(
            bounds.max,
            std::array::from_fn(|i| [0., 20., 5.][i] + delta[i]),
        );
    };
    check(f.session.document(), f.session.revision());
    let applied = f.session.document().clone();
    f.session.undo(f.session.revision(), actor()).unwrap();
    assert_eq!(f.session.document(), &original);
    f.session.redo(f.session.revision(), actor()).unwrap();
    assert_eq!(f.session.document(), &applied);
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
    check(&opened.project.document, 0);
}

#[test]
fn cancelling_points_and_escape_leave_the_source_document_untouched() {
    let mut f = fixture();
    let original = f.session.document().clone();
    let scene = canonical(&f.session);
    let mut c = ModelToolController::default();
    c.set_backend(Box::new(native_fixture::NativeBackend::new()));
    c.begin_placement(
        &f.session,
        &Default::default(),
        &scene,
        &BTreeSet::from([f.moving.body]),
        PlacementMode::Move,
        Channel::Ui,
    )
    .unwrap();
    let before = panel(&mut c).parameters(f.session.document()).unwrap();
    panel(&mut c).toggle_point_pick();
    click(&mut c, &mut f.session, Some(&scene), screen([20., 10., 5.])).unwrap();
    // Cancel Points drops the pending source without staging anything.
    panel(&mut c).toggle_point_pick();
    assert!(!c.placement_point_picking());
    assert_eq!(
        panel(&mut c).parameters(f.session.document()).unwrap(),
        before
    );
    // A complete pair, then Escape: nothing is committed.
    panel(&mut c).toggle_point_pick();
    click(&mut c, &mut f.session, Some(&scene), screen([20., 10., 5.])).unwrap();
    click(&mut c, &mut f.session, Some(&scene), screen([50., 0., 0.])).unwrap();
    assert_ne!(
        panel(&mut c).parameters(f.session.document()).unwrap(),
        before
    );
    c.cancel();
    assert!(!c.placement_active());
    assert_eq!(f.session.document(), &original);
}
