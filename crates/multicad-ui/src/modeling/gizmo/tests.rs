use super::*;
use crate::egui::{Context, Event, Modifiers, PointerButton, RawInput};
use multicad_core::{document::Document, model};

fn camera() -> Camera {
    Camera {
        roll: 0.0,
        yaw: -0.6,
        pitch: 0.6,
        scale: 4.,
        target: [0.; 3],
        perspective: false,
    }
}
fn rect() -> Rect {
    Rect::from_min_max(Pos2::new(10., 10.), Pos2::new(810., 610.))
}
fn parameters() -> EvaluatedParameters {
    multicad_core::parameters::evaluate(&Document::new("Pointer math").unwrap()).unwrap()
}
fn frame(
    ctx: &Context,
    state: &mut Interaction,
    source: &NumericTransform,
    camera: Camera,
    basis: GizmoBasis,
    events: Vec<Event>,
    modifiers: Modifiers,
) -> Gesture {
    // Pointer-math tests drag freely; snapping tests opt into increments.
    snapped_frame(
        ctx,
        state,
        source,
        camera,
        basis,
        events,
        modifiers,
        TransformSnap::OFF,
    )
}
#[allow(clippy::too_many_arguments)]
fn snapped_frame(
    ctx: &Context,
    state: &mut Interaction,
    source: &NumericTransform,
    camera: Camera,
    basis: GizmoBasis,
    mut events: Vec<Event>,
    modifiers: Modifiers,
    snap: TransformSnap,
) -> Gesture {
    let mut result = None;
    events.insert(0, Event::ModifiersChanged(modifiers));
    let output = ctx.run_ui(
        RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, egui::vec2(900., 700.))),
            events,
            ..Default::default()
        },
        |ui| {
            result = Some(
                state
                    .show(
                        ui,
                        &PlacementViewport {
                            scene: None,
                            camera,
                            rect: rect(),
                            basis,
                            start_allowed: true,
                            snap,
                        },
                        source,
                        &parameters(),
                    )
                    .unwrap(),
            );
        },
    );
    output.drop_without_applying_deltas();
    result.unwrap()
}
fn press(p: Pos2) -> Event {
    Event::PointerButton {
        pos: p,
        button: PointerButton::Primary,
        pressed: true,
        modifiers: Modifiers::NONE,
    }
}
fn release(p: Pos2) -> Event {
    Event::PointerButton {
        pos: p,
        button: PointerButton::Primary,
        pressed: false,
        modifiers: Modifiers::NONE,
    }
}
fn point(camera: Camera, p: Point) -> Pos2 {
    projected(camera, rect(), p).unwrap()
}

#[test]
fn actual_axis_pointer_drag_uses_basis_and_original_controls_without_shift_copy() {
    let ctx = Context::default();
    let mut interaction = Interaction::default();
    let source = NumericTransform::default();
    let camera = camera();
    let basis = GizmoBasis {
        axes: [[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]],
    };
    let start = point(camera, [0., 18., 0.]);
    let end = point(camera, [0., 25., 0.]);
    frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![Event::PointerMoved(start)],
        Modifiers::NONE,
    );
    let first = frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![press(start)],
        Modifiers::NONE,
    );
    assert!(first.started && first.input_owned);
    let moved = frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![Event::PointerMoved(end)],
        Modifiers {
            shift: true,
            ..Modifiers::NONE
        },
    );
    assert!(moved.input_owned && !moved.copy_modifier);
    let output = moved.changed.unwrap().resolve(&parameters()).unwrap();
    assert!(output.translation_mm[0].abs() < 1e-6);
    assert!((output.translation_mm[1] - 7.).abs() < 1e-5);
    assert_eq!(source, NumericTransform::default());
    let last = frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![release(end)],
        Modifiers::NONE,
    );
    assert!(last.released && last.moved && last.input_owned);
    assert!(!interaction.active());
}

#[test]
fn actual_plane_drag_and_rotation_ring_have_geometric_pointer_solutions() {
    let parameters = parameters();
    let source = NumericTransform::default();
    let camera = camera();
    let basis = GizmoBasis::default();
    let ctx = Context::default();
    let mut interaction = Interaction::default();
    // Interior of the XY square, away from axes and rotation arcs.
    let start = point(camera, [5.1, 5.1, 0.]);
    let end = point(camera, [11.1, 1.1, 0.]);
    frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![Event::PointerMoved(start)],
        Modifiers::NONE,
    );
    frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![press(start)],
        Modifiers::NONE,
    );
    assert!(matches!(
        interaction.drag.as_ref().unwrap().handle,
        Handle::Plane(2)
    ));
    let output = frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![Event::PointerMoved(end)],
        Modifiers::NONE,
    )
    .changed
    .unwrap()
    .resolve(&parameters)
    .unwrap();
    assert!((output.translation_mm[0] - 6.).abs() < 1e-5);
    assert!((output.translation_mm[1] + 4.).abs() < 1e-5);
    assert!(output.translation_mm[2].abs() < 1e-5);
    frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![release(end)],
        Modifiers::NONE,
    );
    let radius = 82. / camera.scale as f64 * 0.74;
    // The Z ring at 45 degrees is distinct from all projected cardinal axes.
    let t = std::f64::consts::FRAC_PI_4;
    let start = point(
        camera,
        [radius * t.portable_cos(), radius * t.portable_sin(), 0.],
    );
    let end = point(
        camera,
        [-radius * t.portable_sin(), radius * t.portable_cos(), 0.],
    );
    frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![Event::PointerMoved(start)],
        Modifiers::NONE,
    );
    frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![press(start)],
        Modifiers::NONE,
    );
    assert!(matches!(
        interaction.drag.as_ref().unwrap().handle,
        Handle::Rotate(2)
    ));
    let output = frame(
        &ctx,
        &mut interaction,
        &source,
        camera,
        basis,
        vec![Event::PointerMoved(end)],
        Modifiers::NONE,
    )
    .changed
    .unwrap();
    let pose = output.resolve(&parameters).unwrap();
    let expected = NumericTransform {
        angle: model::angle(90.),
        ..Default::default()
    }
    .resolve(&parameters)
    .unwrap();
    for (a, b) in pose
        .rotation
        .into_iter()
        .flatten()
        .zip(expected.rotation.into_iter().flatten())
    {
        assert!((a - b).abs() < 1e-5);
    }
}

#[test]
fn empty_space_does_not_own_selection_and_invalid_coordinate_frames_fail() {
    let ctx = Context::default();
    let mut interaction = Interaction::default();
    let camera = camera();
    let outcome = frame(
        &ctx,
        &mut interaction,
        &NumericTransform::default(),
        camera,
        GizmoBasis::default(),
        vec![press(Pos2::new(50., 50.))],
        Modifiers::NONE,
    );
    assert!(!outcome.input_owned && !interaction.active());
    for axes in [
        [[1., 0., 0.], [1., 0., 0.], [0., 0., 1.]],
        [[1., 0., 0.], [0., 1., 0.], [0., 0., -1.]],
        [[f64::NAN, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    ] {
        assert!(GizmoBasis { axes }.validate().is_err());
    }
}

fn translation(gesture: &Gesture) -> Point {
    gesture
        .changed
        .as_ref()
        .unwrap()
        .resolve(&parameters())
        .unwrap()
        .translation_mm
}

#[test]
fn move_snaps_to_the_grid_and_shift_frees_it_without_ending_the_gesture() {
    let ctx = Context::default();
    let mut interaction = Interaction::default();
    let source = NumericTransform::default();
    let camera = camera();
    let basis = GizmoBasis::default();
    let snap = TransformSnap::grid(10.);
    let mut run = |source: &NumericTransform, events, shift| {
        snapped_frame(
            &ctx,
            &mut interaction,
            source,
            camera,
            basis,
            events,
            Modifiers {
                shift,
                ..Modifiers::NONE
            },
            snap,
        )
    };
    let start = point(camera, [0., 18., 0.]);
    run(&source, vec![Event::PointerMoved(start)], false);
    assert!(run(&source, vec![press(start)], false).started);
    // 7 mm along Y stages exactly one 10 mm grid step.
    let snapped = run(
        &source,
        vec![Event::PointerMoved(point(camera, [0., 25., 0.]))],
        false,
    );
    assert_eq!(translation(&snapped), [0., 10., 0.]);
    let staged = snapped.changed.clone().unwrap();
    // Holding Shift re-solves the same pointer without increments.
    let free = run(&staged, vec![], true);
    assert!(free.input_owned && !free.released);
    let t = translation(&free);
    assert!(t[0].abs() < 1e-9 && (t[1] - 7.).abs() < 1e-5 && t[2].abs() < 1e-9);
    // Releasing Shift restores the grid; the gesture is still captured.
    let restored = run(&free.changed.clone().unwrap(), vec![], false);
    assert!(restored.input_owned && !restored.released);
    assert_eq!(translation(&restored), [0., 10., 0.]);
    // 14 mm rounds down to one step; the pivot is never changed by snapping.
    let end = point(camera, [0., 32., 0.]);
    let last = run(&staged, vec![Event::PointerMoved(end)], false);
    assert!(last.changed.is_none(), "still one grid step");
    let done = run(&staged, vec![release(end)], false);
    assert!(done.released && done.moved);
    assert!(!interaction.active());
    assert_eq!(staged.pivot, source.pivot);
    assert_eq!(source, NumericTransform::default());
}

#[test]
fn plane_snapping_uses_the_active_rotated_coordinate_frame() {
    let ctx = Context::default();
    let mut interaction = Interaction::default();
    let source = NumericTransform::default();
    let camera = camera();
    let (c, s) = (
        30_f64.to_radians().portable_cos(),
        30_f64.to_radians().portable_sin(),
    );
    let (u, v) = ([c, s, 0.], [-s, c, 0.]);
    let basis = GizmoBasis {
        axes: [u, v, [0., 0., 1.]],
    };
    let at = |a: f64, b: f64| point(camera, add(mul(u, a), mul(v, b)));
    let snap = TransformSnap::grid(10.);
    let run = |interaction: &mut Interaction, events| {
        snapped_frame(
            &ctx,
            interaction,
            &source,
            camera,
            basis,
            events,
            Modifiers::NONE,
            snap,
        )
    };
    run(&mut interaction, vec![Event::PointerMoved(at(5.1, 5.1))]);
    run(&mut interaction, vec![press(at(5.1, 5.1))]);
    assert!(matches!(
        interaction.drag.as_ref().unwrap().handle,
        Handle::Plane(2)
    ));
    let moved = run(&mut interaction, vec![Event::PointerMoved(at(22.1, 9.1))]);
    let t = translation(&moved);
    let expected = mul(u, 20.);
    for i in 0..3 {
        assert!((t[i] - expected[i]).abs() < 1e-9, "{t:?} != {expected:?}");
    }
    let readout = interaction.drag.as_ref().unwrap().readout.clone().unwrap();
    assert!(readout.contains("grid 10"), "{readout}");
}

#[test]
fn rotate_snaps_to_45_degree_steps_shows_the_angle_and_shift_frees_it() {
    let parameters = parameters();
    let ctx = Context::default();
    let mut interaction = Interaction::default();
    let source = NumericTransform::default();
    let camera = camera();
    let snap = TransformSnap::grid(10.);
    let radius = 82. / camera.scale as f64 * 0.74;
    let ring = |degrees: f64| {
        let t = degrees.to_radians();
        point(
            camera,
            [radius * t.portable_cos(), radius * t.portable_sin(), 0.],
        )
    };
    let run = |interaction: &mut Interaction, source: &NumericTransform, events, shift| {
        snapped_frame(
            &ctx,
            interaction,
            source,
            camera,
            GizmoBasis::default(),
            events,
            Modifiers {
                shift,
                ..Modifiers::NONE
            },
            snap,
        )
    };
    run(
        &mut interaction,
        &source,
        vec![Event::PointerMoved(ring(45.))],
        false,
    );
    run(&mut interaction, &source, vec![press(ring(45.))], false);
    // Sweep 80 degrees in small samples; the staged angle is 90.
    let samples = (1..=8).map(|k| Event::PointerMoved(ring(45. + 10. * k as f64)));
    let snapped = run(&mut interaction, &source, samples.collect(), false);
    let staged = snapped.changed.clone().unwrap();
    let angle =
        |t: &NumericTransform| Values::read(t, &parameters).unwrap().angle_deg * t.axis[2].signum();
    assert!((angle(&staged) - 90.).abs() < 1e-9, "{staged:?}");
    let readout = interaction.drag.as_ref().unwrap().readout.clone().unwrap();
    assert!(
        readout.contains("90.0°") && readout.contains("45° steps"),
        "{readout}"
    );
    let free = run(&mut interaction, &staged, vec![], true)
        .changed
        .unwrap();
    assert!((angle(&free) - 80.).abs() < 1e-3, "{free:?}");
    let readout = interaction.drag.as_ref().unwrap().readout.clone().unwrap();
    assert!(readout.contains("free (Shift)"), "{readout}");
    let restored = run(&mut interaction, &free, vec![], false);
    assert!(restored.input_owned && !restored.released);
    assert!((angle(&restored.changed.unwrap()) - 90.).abs() < 1e-9);
    assert_eq!(staged.pivot, source.pivot);
}
