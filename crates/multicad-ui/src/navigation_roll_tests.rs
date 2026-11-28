//! CAD-419: roll arrows on the navigation widget, driven by real clicks.
use super::*;
use std::time::{Duration, Instant};

struct Widget {
    ctx: egui::Context,
    camera: Camera,
    axis: Option<[f32; 3]>,
    id: egui::Id,
}
impl Widget {
    fn new(camera: Camera, axis: Option<[f32; 3]>) -> Self {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let mut widget = Self {
            ctx,
            camera,
            axis,
            id: egui::Id::NULL,
        };
        widget.frame(vec![]).drop_without_applying_deltas();
        widget
    }
    fn frame(&mut self, events: Vec<egui::Event>) -> egui::FullOutput {
        let Self {
            ctx,
            camera,
            axis,
            id,
        } = self;
        ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(600.0, 400.0))),
                events,
                ..Default::default()
            },
            |ui| {
                egui::Area::new(ui.id().with("cross"))
                    .pivot(Align2::RIGHT_TOP)
                    .fixed_pos(Pos2::new(500.0, 10.0))
                    .movable(false)
                    .show(ui.ctx(), |ui| {
                        *id = orientation_about(ui, camera, *axis).id;
                    });
            },
        )
    }
    fn button(&mut self, roll: Roll) -> (Rect, String) {
        let out = self.frame(vec![]);
        let found = out
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .find_map(|(_, n)| {
                (n.label() == Some(roll.label())).then(|| {
                    let b = n.bounds().unwrap();
                    Rect::from_min_max(
                        Pos2::new(b.x0 as f32, b.y0 as f32),
                        Pos2::new(b.x1 as f32, b.y1 as f32),
                    )
                })
            })
            .unwrap_or_else(|| panic!("Missing {}", roll.label()));
        out.drop_without_applying_deltas();
        // Hover long enough for the tooltip, which must explain the action.
        let pos = found.center();
        let deadline = Instant::now() + Duration::from_secs(5);
        self.frame(vec![egui::Event::PointerMoved(pos)])
            .drop_without_applying_deltas();
        let tooltip = loop {
            let out = self.frame(vec![]);
            let text = out.shapes.iter().find_map(|s| match &s.shape {
                egui::Shape::Text(t) if t.galley.text().starts_with("Roll the view 90°") => {
                    Some(t.galley.text().to_owned())
                }
                _ => None,
            });
            out.drop_without_applying_deltas();
            if let Some(text) = text {
                break text;
            }
            assert!(Instant::now() < deadline, "no tooltip for {}", roll.label());
            std::thread::sleep(Duration::from_millis(5));
        };
        (found, tooltip)
    }
    /// Click a roll arrow and wait until its animation has landed.
    fn click(&mut self, roll: Roll) {
        let (rect, _) = self.button(roll);
        let pos = rect.center();
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ])
            .drop_without_applying_deltas();
        }
        assert!(roll_animating(&self.ctx, self.id), "the roll animates");
        self.settle();
    }
    fn settle(&mut self) {
        let deadline = Instant::now() + Duration::from_secs(5);
        while roll_animating(&self.ctx, self.id) {
            assert!(Instant::now() < deadline, "roll animation never finished");
            self.frame(vec![]).drop_without_applying_deltas();
            std::thread::sleep(Duration::from_millis(2));
        }
    }
}

fn neg(v: [f32; 3]) -> [f32; 3] {
    v.map(|c| -c)
}

#[test]
fn roll_arrows_are_labelled_and_explained() {
    let mut widget = Widget::new(Camera::default(), None);
    for roll in Roll::ALL {
        let (rect, tooltip) = widget.button(roll);
        assert!(rect.width() >= 14.0, "{} is a usable target", roll.label());
        assert!(tooltip.contains(match roll {
            Roll::Clockwise => "clockwise",
            Roll::CounterClockwise => "counter-clockwise",
        }));
    }
    let (ccw, _) = widget.button(Roll::CounterClockwise);
    let (cw, _) = widget.button(Roll::Clockwise);
    assert!(
        ccw.right() <= cw.left() + 0.5,
        "counter-clockwise sits left"
    );
    // Both sit in the top-right corner of the 104 px widget (500 - 104 .. 500).
    for rect in [ccw, cw] {
        assert!(rect.top() <= 12.0 && rect.right() <= 500.5 && rect.left() >= 500.0 - 52.0);
    }
}

#[test]
fn one_click_rolls_a_quarter_turn_keeping_target_and_zoom() {
    for (roll, expected_up) in [
        (Roll::Clockwise, neg as fn([f32; 3]) -> [f32; 3]),
        (Roll::CounterClockwise, |v| v),
    ] {
        let start = Camera::default();
        let [right, up, toward_eye] = start.basis();
        let mut widget = Widget::new(start, None);
        widget.click(roll);
        let camera = widget.camera;
        let [new_right, new_up, new_toward_eye] = camera.basis();
        // Clockwise shows the former up direction at the right.
        assert_eq!(new_up, expected_up(right), "{roll:?}");
        assert_eq!(
            new_right,
            match roll {
                Roll::Clockwise => up,
                Roll::CounterClockwise => neg(up),
            }
        );
        assert_eq!(new_toward_eye, toward_eye, "rolls about the view direction");
        assert_eq!(camera.target, start.target);
        assert_eq!(camera.scale, start.scale);
        assert_eq!((camera.yaw, camera.pitch), (start.yaw, start.pitch));
    }
}

#[test]
fn four_clicks_return_the_exact_up_vector_in_both_directions() {
    for roll in Roll::ALL {
        let start = Camera::default();
        let mut widget = Widget::new(start, None);
        let mut ups = vec![];
        for _ in 0..4 {
            widget.click(roll);
            ups.push(widget.camera.basis()[1]);
        }
        assert_eq!(widget.camera.basis(), start.basis(), "{roll:?}");
        assert_eq!(widget.camera.roll, 0.0);
        assert_eq!(widget.camera, start);
        // Each step was a distinct quarter turn.
        for pair in ups.windows(2) {
            assert_ne!(pair[0], pair[1]);
        }
    }
    // Opposite clicks cancel exactly too.
    let start = Camera::default();
    let mut widget = Widget::new(start, None);
    widget.click(Roll::Clockwise);
    widget.click(Roll::CounterClockwise);
    assert_eq!(widget.camera, start);
}

#[test]
fn the_roll_animates_through_intermediate_angles() {
    let mut widget = Widget::new(Camera::default(), None);
    let (rect, _) = widget.button(Roll::Clockwise);
    let pos = rect.center();
    for pressed in [true, false] {
        widget
            .frame(vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ])
            .drop_without_applying_deltas();
    }
    let mut seen = vec![];
    let deadline = Instant::now() + Duration::from_secs(5);
    while roll_animating(&widget.ctx, widget.id) {
        assert!(Instant::now() < deadline);
        widget.frame(vec![]).drop_without_applying_deltas();
        seen.push(widget.camera.roll);
        std::thread::sleep(Duration::from_millis(2));
    }
    assert!(
        seen.iter()
            .any(|r| *r > 0.0 && *r < std::f32::consts::FRAC_PI_2),
        "{seen:?}"
    );
    assert_eq!(widget.camera.roll, std::f32::consts::FRAC_PI_2);
}

#[test]
fn a_roll_axis_turns_about_the_sketch_normal() {
    // An XZ sketch seen from +Y side normal (0, -1, 0), with the view orbited
    // away from it.
    let normal = [0.0, -1.0, 0.0];
    let start = Camera {
        yaw: 0.4,
        pitch: 0.3,
        ..Default::default()
    };
    let mut widget = Widget::new(start, Some(normal));
    widget.click(Roll::Clockwise);
    let [_, up, toward_eye] = widget.camera.basis();
    for i in 0..3 {
        assert!((toward_eye[i] - normal[i]).abs() < 1e-6, "{toward_eye:?}");
    }
    assert_eq!(widget.camera.target, start.target);
    assert_eq!(widget.camera.scale, start.scale);
    // A clockwise roll puts the former right direction at the bottom.
    let mut level = widget.camera;
    level.roll = 0.0;
    let [level_right, ..] = level.basis();
    assert_eq!(up, neg(level_right));
    for _ in 0..3 {
        widget.click(Roll::Clockwise);
    }
    assert_eq!(widget.camera.roll, 0.0);
    assert_eq!(widget.camera.basis(), level.basis());
}
