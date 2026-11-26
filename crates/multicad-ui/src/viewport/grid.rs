use crate::camera::Camera;
use eframe::egui::Rect;
use multicad_precision::PortableMath as _;
/// Spacing of the floor grid at this zoom: the configured base spacing
/// coarsened by powers of ten so lines stay at least 8 px apart and at most
/// about 100 fit across the view. Dimension drags snap to this same spacing
/// (CAD-402), also in edge-on views where the grid lines themselves are hidden.
pub fn grid_step(camera: Camera, rect: Rect, base_mm: f64) -> Option<f64> {
    let facing = (camera.pitch.portable_sin().abs() as f64).max(0.05);
    if !base_mm.is_finite()
        || base_mm <= 0.
        || !camera.scale.is_finite()
        || camera.scale <= 0.
        || !facing.is_finite()
        || !rect.is_finite()
    {
        return None;
    }
    let scale = camera.scale as f64;
    let factor = (8. / (base_mm * scale)).max(1.).portable_log10().ceil();
    let mut step = base_mm * 10_f64.portable_powf(factor);
    let extent = rect.width().max(rect.height()) as f64 / scale / facing;
    while extent / step > 100. {
        step *= 10.;
    }
    step.is_finite().then_some(step)
}

/// Coarsen display lines at distant zoom without changing the configured sketch
/// snap spacing. Center the finite display patch on the current view.
pub fn floor_grid(camera: Camera, rect: Rect, base_mm: f64) -> Vec<[[f64; 3]; 2]> {
    let facing = camera.pitch.portable_sin().abs() as f64;
    if facing < 0.05 || !facing.is_finite() || camera.target.iter().any(|v| !v.is_finite()) {
        return vec![];
    }
    let Some(step) = grid_step(camera, rect, base_mm) else {
        return vec![];
    };
    let scale = camera.scale as f64;
    let extent = rect.width().max(rect.height()) as f64 / scale / facing;
    let count = (extent / step).ceil().clamp(4., 100.) as i32;
    let center = [
        (camera.target[0] as f64 / step).round() * step,
        (camera.target[1] as f64 / step).round() * step,
    ];
    let half = count as f64 * step;
    let mut lines = vec![];
    for i in -count..=count {
        let offset = i as f64 * step;
        lines.extend([
            [
                [center[0] + offset, center[1] - half, 0.],
                [center[0] + offset, center[1] + half, 0.],
            ],
            [
                [center[0] - half, center[1] + offset, 0.],
                [center[0] + half, center[1] + offset, 0.],
            ],
        ]);
    }
    lines
}

#[cfg(test)]
mod tests {
    use super::*;
    use eframe::egui::{Pos2, vec2};

    fn camera(scale: f32, pitch: f32) -> Camera {
        Camera {
            yaw: 0.,
            pitch,
            scale,
            target: [0.; 3],
            perspective: false,
        }
    }

    #[test]
    fn drag_step_is_the_displayed_floor_grid_spacing() {
        let rect = Rect::from_min_size(Pos2::ZERO, vec2(800., 600.));
        let top = std::f32::consts::FRAC_PI_2;
        for (scale, expected) in [(3.713, 10.), (0.5, 100.), (0.05, 1000.)] {
            let camera = camera(scale, top);
            assert_eq!(grid_step(camera, rect, 10.), Some(expected), "{scale}");
            let lines = floor_grid(camera, rect, 10.);
            let spacing = lines[2][0][0] - lines[0][0][0];
            assert_eq!(spacing, expected, "{scale}");
        }
        // Edge-on views hide the lines but keep the zoom's step for drags.
        assert!(floor_grid(camera(3.713, 0.), rect, 10.).is_empty());
        assert_eq!(grid_step(camera(3.713, 0.), rect, 10.), Some(100.));
        assert_eq!(grid_step(camera(3.713, top), rect, 0.), None);
        assert_eq!(grid_step(camera(f32::NAN, top), rect, 10.), None);
    }
}
