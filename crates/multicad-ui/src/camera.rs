//! Camera math is independent of native windows and model ownership.
use eframe::egui::{Rect, Vec2};
use multicad_precision::PortableMath as _;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Camera {
    pub yaw: f32,
    pub pitch: f32,
    pub scale: f32,
    pub target: [f32; 3],
    pub perspective: bool,
    /// Turn of the image about the view direction, in radians; positive
    /// turns the image clockwise on screen (CAD-419). Roll buttons keep it
    /// on exact quarter turns in `[0, 2π)`.
    pub roll: f32,
}
impl Default for Camera {
    fn default() -> Self {
        Self {
            yaw: -0.6,
            pitch: 0.6,
            scale: 2.4,
            target: [90.0, 40.0, 12.0],
            perspective: false,
            roll: 0.0,
        }
    }
}
impl Camera {
    /// Yaw at which a straight top-down view shows +X to the right and +Y up
    /// (and a bottom-up view +X right, +Y down), as in a standard plan view.
    pub const PLAN_YAW: f32 = std::f32::consts::PI;
    /// Screen right, up and toward-the-eye directions in world space. The
    /// basis is right-handed (right × up = toward the eye), so the image is
    /// never a mirror of the model (CAD-382).
    pub fn basis(&self) -> [[f32; 3]; 3] {
        // Axis headings (0, ±π/2, π) are exact, so plan views map model
        // coordinates without f32 rounding of sin(π) and friends.
        let exact = |(s, c): (f32, f32)| {
            if s.abs() < 1e-6 {
                (0., c.signum())
            } else if c.abs() < 1e-6 {
                (s.signum(), 0.)
            } else {
                (s, c)
            }
        };
        let (s, c) = exact(self.yaw.portable_sin_cos());
        let (sp, cp) = exact(self.pitch.portable_sin_cos());
        let right = [-c, s, 0.];
        let up = [-sp * s, -sp * c, cp];
        let toward_eye = [cp * s, cp * c, sp];
        if self.roll == 0. {
            return [right, up, toward_eye];
        }
        // Rolling the image clockwise shows the former up direction at the
        // right; quarter turns stay exact like the axis headings above.
        let (sr, cr) = exact(self.roll.portable_sin_cos());
        [
            std::array::from_fn(|i| cr * right[i] + sr * up[i]),
            std::array::from_fn(|i| cr * up[i] - sr * right[i]),
            toward_eye,
        ]
    }
    /// The roll as whole quarter turns, `0..4`.
    pub fn roll_quarters(&self) -> i32 {
        ((self.roll / std::f32::consts::FRAC_PI_2).round() as i32).rem_euclid(4)
    }
    /// The exact roll angle of `quarters` quarter turns, in `[0, 2π)`.
    pub fn quarter_roll(quarters: i32) -> f32 {
        quarters.rem_euclid(4) as f32 * std::f32::consts::FRAC_PI_2
    }
    /// Look along `toward_eye` (a unit vector from the target to the eye)
    /// without changing target, zoom or roll. A vertical direction keeps the
    /// current heading.
    pub fn look_from(&mut self, toward_eye: [f32; 3]) {
        let horizontal = toward_eye[0].portable_hypot(toward_eye[1]);
        if horizontal > 1e-6 {
            self.yaw = toward_eye[0].portable_atan2(toward_eye[1]);
        }
        self.pitch = toward_eye[2].clamp(-1., 1.).portable_asin();
    }
    pub fn orbit(&mut self, delta: Vec2) {
        self.yaw += delta.x * 0.01;
        self.pitch = (self.pitch + delta.y * 0.01).clamp(-1.57, 1.57);
    }
    pub fn pan(&mut self, delta: Vec2) {
        let [right, up_axis, _] = self.basis();
        let x = delta.x / self.scale;
        let up = -delta.y / self.scale;
        for i in 0..3 {
            self.target[i] -= right[i] * x + up_axis[i] * up;
        }
    }
    /// Keyboard panning uses screen units per second, independent of frame rate.
    pub fn pan_keyboard(&mut self, direction: Vec2, seconds: f32) -> bool {
        if !direction.is_finite()
            || !seconds.is_finite()
            || seconds <= 0.
            || direction == Vec2::ZERO
        {
            return false;
        }
        let direction = if direction.length_sq() > 1. {
            direction.normalized()
        } else {
            direction
        };
        self.pan(direction * 360.0 * seconds.min(0.05));
        true
    }
    pub fn zoom(&mut self, scroll: f32) {
        self.scale = (self.scale * (scroll * 0.002).portable_exp()).clamp(0.02, 100.0);
    }
    pub fn project(&self, point: [f32; 3], rect: Rect) -> [f32; 3] {
        let relative: [f32; 3] = std::array::from_fn(|i| point[i] - self.target[i]);
        let [right, up, depth] = self
            .basis()
            .map(|axis| (0..3).map(|i| axis[i] * relative[i]).sum::<f32>());
        let distance = rect.height().max(1.0) / self.scale * 2.0;
        let z = distance - depth;
        let factor = if self.perspective {
            distance / z.max(0.01)
        } else {
            1.0
        };
        let p = rect.center() + Vec2::new(right, -up) * self.scale * factor;
        [p.x, p.y, z]
    }
    /// A normalized world ray shared by surface picking and datum intersections.
    pub fn ray(&self, point: eframe::egui::Pos2, rect: Rect) -> Option<crate::viewport::Ray> {
        if !self.scale.is_finite() || self.scale <= 0. || !point.is_finite() || !rect.is_finite() {
            return None;
        }
        let [right, up, forward] = self.basis().map(|axis| axis.map(f64::from));
        let x = ((point.x - rect.center().x) / self.scale) as f64;
        let y = ((rect.center().y - point.y) / self.scale) as f64;
        let d = (rect.height().max(1.) / self.scale * 2.) as f64;
        let offset: [f64; 3] = std::array::from_fn(|i| right[i] * x + up[i] * y);
        let origin: [f64; 3] = std::array::from_fn(|i| {
            self.target[i] as f64 + forward[i] * d + if self.perspective { 0. } else { offset[i] }
        });
        let direction = std::array::from_fn(|i| {
            if self.perspective {
                offset[i] - forward[i] * d
            } else {
                -forward[i]
            }
        });
        crate::viewport::Ray::new(origin, direction).ok()
    }
    /// Inverse projection against a datum; edge-on planes have no unique hit.
    pub fn on_plane(
        &self,
        point: eframe::egui::Pos2,
        rect: Rect,
        frame: multicad_core::sketch::Frame,
    ) -> Option<[f64; 2]> {
        let ray = self.ray(point, rect)?;
        let origin = ray.origin();
        let direction = ray.direction();
        let dot = multicad_core::sketch::dot;
        let denominator = dot(direction, frame.n);
        if denominator.abs() < 1e-8 {
            return None;
        }
        let t = dot(
            std::array::from_fn(|i| frame.origin[i] - origin[i]),
            frame.n,
        ) / denominator;
        if t < 0. || !t.is_finite() {
            return None;
        }
        let hit = std::array::from_fn(|i| origin[i] + direction[i] * t);
        Some(frame.coordinates(hit))
    }
    /// Clip source segments before perspective division so geometry behind the
    /// eye never folds across the viewport.
    pub fn segment(&self, a: [f64; 3], b: [f64; 3], rect: Rect) -> Option<[eframe::egui::Pos2; 2]> {
        self.clip_segment(a, b, rect).map(|points| {
            points.map(|p| {
                let p = self.project(p.map(|v| v as f32), rect);
                eframe::egui::Pos2::new(p[0], p[1])
            })
        })
    }
    pub fn clip_segment(
        &self,
        mut a: [f64; 3],
        mut b: [f64; 3],
        rect: Rect,
    ) -> Option<[[f64; 3]; 2]> {
        let depth = |p: [f64; 3]| self.project(p.map(|v| v as f32), rect)[2];
        let za = depth(a);
        let zb = depth(b);
        let near = 0.1;
        if self.perspective {
            if za < near && zb < near {
                return None;
            }
            if za < near || zb < near {
                let t = ((near - za) / (zb - za)) as f64;
                let clipped = std::array::from_fn(|i| a[i] + (b[i] - a[i]) * t);
                if za < near {
                    a = clipped;
                } else {
                    b = clipped;
                }
            }
        }
        Some([a, b])
    }
    pub fn fit(&mut self, points: &[[f32; 3]], rect: Rect) {
        if points.is_empty() {
            return;
        }
        let mut low = [f32::INFINITY; 3];
        let mut high = [f32::NEG_INFINITY; 3];
        for p in points {
            for i in 0..3 {
                low[i] = low[i].min(p[i]);
                high[i] = high[i].max(p[i]);
            }
        }
        self.target = std::array::from_fn(|i| (low[i] + high[i]) * 0.5);
        // Bound a sphere in either projection, leaving space for camera controls.
        let radius = (0..3)
            .map(|i| (high[i] - low[i]).powi(2))
            .sum::<f32>()
            .sqrt()
            * 0.5;
        self.scale = (rect.width().min(rect.height()) * 0.32 / radius.max(1.0)).clamp(0.02, 100.0);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use eframe::egui::Pos2;
    fn rect() -> Rect {
        Rect::from_min_size(Pos2::ZERO, Vec2::new(800.0, 600.0))
    }
    #[test]
    fn pan_tracks_pointer_and_zoom_is_bounded() {
        for perspective in [false, true] {
            let mut c = Camera {
                perspective,
                ..Default::default()
            };
            let p = c.target;
            let a = c.project(p, rect());
            c.pan(Vec2::new(30.0, -20.0));
            let b = c.project(p, rect());
            assert!((b[0] - a[0] - 30.0).abs() < 0.001);
            assert!((b[1] - a[1] + 20.0).abs() < 0.001);
            c.zoom(100000.0);
            assert_eq!(c.scale, 100.0);
            c.zoom(-100000.0);
            assert_eq!(c.scale, 0.02);
        }
    }
    #[test]
    fn near_surface_orbit_tracks_drag_and_fit_contains_scene() {
        let mut c = Camera {
            yaw: 0.0,
            pitch: 0.0,
            target: [0.0; 3],
            ..Default::default()
        };
        let p = [0.0, 20.0, 0.0];
        let a = c.project(p, rect());
        c.orbit(Vec2::new(10.0, 10.0));
        let b = c.project(p, rect());
        assert!(b[0] > a[0] && b[1] > a[1]);
        let points = [[-500.0, -200.0, -50.0], [1000.0, 200.0, 150.0]];
        for perspective in [false, true] {
            c.perspective = perspective;
            c.fit(&points, rect());
            for p in points {
                let q = c.project(p, rect());
                assert!(rect().contains(Pos2::new(q[0], q[1])));
            }
        }
    }
}
#[cfg(test)]
mod handedness_tests {
    use super::*;
    #[test]
    fn every_view_is_right_handed_so_the_model_is_never_mirrored() {
        let cross = |a: [f32; 3], b: [f32; 3]| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };
        for yaw in [-3.0, -0.6, 0.0, 0.4, Camera::PLAN_YAW] {
            for pitch in [-1.5707964, -0.6, 0.0, 0.6, 1.5707964] {
                let camera = Camera {
                    yaw,
                    pitch,
                    ..Default::default()
                };
                let [right, up, toward_eye] = camera.basis();
                let normal = cross(right, up);
                for i in 0..3 {
                    assert!((normal[i] - toward_eye[i]).abs() < 1e-5, "{yaw} {pitch}");
                }
            }
        }
        // A plan view shows +X right and +Y up.
        let rect = Rect::from_min_size(eframe::egui::Pos2::ZERO, Vec2::new(800., 600.));
        let mut top = Camera {
            target: [0.; 3],
            ..Default::default()
        };
        top.orient(crate::navigation::AxisView::PositiveZ);
        let origin = top.project([0.; 3], rect);
        let x = top.project([10., 0., 0.], rect);
        let y = top.project([0., 10., 0.], rect);
        assert!(x[0] > origin[0] && (x[1] - origin[1]).abs() < 1e-3);
        assert!(y[1] < origin[1] && (y[0] - origin[0]).abs() < 1e-3);
    }
}
#[cfg(test)]
mod plane_tests {
    use super::*;
    use multicad_core::sketch::Datum;
    #[test]
    fn ray_inverse_matches_projection_in_both_modes_and_rejects_edge_on_planes() {
        let rect = Rect::from_min_size(eframe::egui::Pos2::ZERO, Vec2::new(1000., 700.));
        for perspective in [false, true] {
            for datum in Datum::ALL {
                let camera = Camera {
                    target: [0.; 3],
                    perspective,
                    ..Default::default()
                };
                let frame = datum.frame();
                let expected = [20., 30.];
                let p = camera.project(frame.point(expected).map(|v| v as f32), rect);
                let actual = camera
                    .on_plane(eframe::egui::Pos2::new(p[0], p[1]), rect, frame)
                    .unwrap();
                assert!((actual[0] - expected[0]).abs() < 0.001);
                assert!((actual[1] - expected[1]).abs() < 0.001);
            }
        }
        let camera = Camera {
            pitch: 0.,
            yaw: 0.,
            target: [0.; 3],
            ..Default::default()
        };
        assert!(
            camera
                .on_plane(rect.center(), rect, Datum::Xy.frame())
                .is_none()
        );
    }
}
#[cfg(test)]
mod keyboard_tests {
    use super::*;
    #[test]
    fn keyboard_pan_speed_is_independent_of_frame_rate_and_bounds_bad_input() {
        let mut a = Camera::default();
        let mut b = a;
        for _ in 0..60 {
            a.pan_keyboard(Vec2::new(1., -1.), 1. / 60.);
        }
        for _ in 0..120 {
            b.pan_keyboard(Vec2::new(1., -1.), 1. / 120.);
        }
        for (a, b) in a.target.into_iter().zip(b.target) {
            assert!((a - b).abs() < 0.002);
        }
        let before = a;
        assert!(!a.pan_keyboard(Vec2::X, f32::NAN));
        assert!(!a.pan_keyboard(Vec2::ZERO, 1.));
        assert_eq!(a, before);
        let mut normal = before;
        normal.pan_keyboard(Vec2::X, 0.05);
        a.pan_keyboard(Vec2::X, 10.);
        assert_eq!(a, normal);
    }
}
