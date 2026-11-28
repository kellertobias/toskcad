//! CAD-419 in the real application viewport: the roll arrows turn the Model
//! view about its viewing direction and a sketch view about the sketch normal.
use super::ui_harness::Harness;
use super::*;
use std::f32::consts::FRAC_PI_2;

fn close(a: [f32; 3], b: [f32; 3]) -> bool {
    (0..3).all(|i| (a[i] - b[i]).abs() < 1e-5)
}

#[test]
fn model_roll_arrows_turn_the_view_about_its_direction() {
    let mut h = Harness::new(Document::new("Roll model").unwrap());
    let start = h.app.current().sketch.camera;
    let [right, up, toward_eye] = start.basis();
    h.click_label("Roll view counter-clockwise");
    h.until("the roll to land", |app| {
        app.current().sketch.camera.roll == 3. * FRAC_PI_2
    });
    let camera = h.app.current().sketch.camera;
    let [_, new_up, new_toward_eye] = camera.basis();
    assert_eq!(
        new_up, right,
        "counter-clockwise shows the former right at the top"
    );
    assert_eq!(new_toward_eye, toward_eye);
    assert_eq!((camera.target, camera.scale), (start.target, start.scale));
    for _ in 0..3 {
        let quarters = h.app.current().sketch.camera.roll_quarters();
        h.click_label("Roll view counter-clockwise");
        h.until("the next quarter turn", |app| {
            let c = app.current().sketch.camera;
            c.roll == crate::camera::Camera::quarter_roll(quarters - 1)
        });
    }
    let camera = h.app.current().sketch.camera;
    assert_eq!(
        camera.basis()[1],
        up,
        "four clicks restore the exact up vector"
    );
    assert_eq!(camera, start);
}

#[test]
fn sketch_roll_arrows_turn_about_the_sketch_normal() {
    let mut h = Harness::new(Document::new("Roll sketch").unwrap());
    h.app
        .invoke(Command::NewSketch, &["xz"], Channel::Ui)
        .unwrap();
    h.until("the sketch editor", |app| {
        app.current().sketch.editing && app.mode == Mode::Sketch
    });
    let level = h.app.current().sketch.camera;
    let normal = level.basis()[2];
    let frame = multicad_core::sketch::Datum::Xz.frame();
    assert!(
        close(normal, frame.n.map(|v| v as f32)),
        "the sketch opens facing its plane"
    );
    // Orbit away from the plane; the roll still turns about its normal.
    let camera = &mut h.app.current_mut().sketch.camera;
    camera.yaw += 0.5;
    camera.pitch = 0.3;
    let orbited = *camera;
    h.click_label("Roll view clockwise");
    h.until("the roll to land", |app| {
        app.current().sketch.camera.roll == FRAC_PI_2
    });
    let camera = h.app.current().sketch.camera;
    let [_, up, toward_eye] = camera.basis();
    assert!(
        close(toward_eye, normal),
        "{toward_eye:?} looks along the sketch normal"
    );
    assert!(
        close(up, level.basis()[0].map(|v| -v)),
        "clockwise puts the plane's right at the bottom"
    );
    assert_eq!(
        (camera.target, camera.scale),
        (orbited.target, orbited.scale)
    );
    for quarters in [2, 3, 0] {
        h.click_label("Roll view clockwise");
        h.until("the next quarter turn", |app| {
            app.current().sketch.camera.roll == crate::camera::Camera::quarter_roll(quarters)
        });
    }
    let camera = h.app.current().sketch.camera;
    assert_eq!(camera.roll, 0.);
    assert!(close(camera.basis()[1], level.basis()[1]));
}
