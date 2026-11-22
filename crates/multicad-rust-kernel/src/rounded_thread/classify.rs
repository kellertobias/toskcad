//! Whole-box material classification for a rounded external thread.
//!
//! The material is `rho <= radius(phase)` for `0 <= z <= length`, with the
//! phase `z - rate * theta`. A box is proved only when its whole radial range
//! lies beyond every radius its phase range can reach; anything near the
//! thread surface or an end cap stays unresolved.
use super::{Thread, ThroughBore};
use crate::helical_thread::BoxClassification;

impl Thread {
    /// Classify the closed local box `[[x0, x1], [y0, y1], [z0, z1]]` given in
    /// the thread chart (axis along z from the start cap).
    pub fn classify_box(&self, bounds: [[f64; 2]; 3]) -> BoxClassification {
        use BoxClassification::*;
        if bounds.iter().flatten().any(|value| !value.is_finite())
            || bounds.iter().any(|[lo, hi]| lo > hi)
        {
            return UnresolvedBoundary;
        }
        let [[x0, x1], [y0, y1], [z0, z1]] = bounds;
        let spec = self.spec;
        let magnitude = bounds
            .iter()
            .flatten()
            .map(|value| value.abs())
            .fold(spec.major_radius_mm.max(spec.length_mm).max(1.), f64::max);
        let rounding = 64. * f64::EPSILON * magnitude;
        if z1 < -rounding || z0 > spec.length_mm + rounding {
            return Outside;
        }
        // Radial range of the box about the axis.
        let near = |lo: f64, hi: f64| {
            if lo > 0. {
                lo
            } else if hi < 0. {
                -hi
            } else {
                0.
            }
        };
        let far = |lo: f64, hi: f64| lo.abs().max(hi.abs());
        let rho_lo = near(x0, x1).hypot(near(y0, y1)) * (1. - 4. * f64::EPSILON);
        let rho_hi = far(x0, x1).hypot(far(y0, y1)) * (1. + 4. * f64::EPSILON);
        if rho_lo > spec.major_radius_mm + rounding {
            return Outside;
        }
        // The phase range: an angular enclosure of the box about its centre.
        let center = [0.5 * (x0 + x1), 0.5 * (y0 + y1)];
        let reach = 0.5 * (x1 - x0).hypot(y1 - y0);
        let rho_center = center[0].hypot(center[1]);
        let [low, high] = if rho_center > 2. * reach && reach.is_finite() {
            let theta = center[1].atan2(center[0]);
            // A point within `reach` of the centre deviates in angle by at
            // most asin(reach / rho) <= reach / (rho - reach).
            let angular = reach / (rho_center - reach) + 16. * f64::EPSILON;
            let rate = self.helical_rate();
            let phase = 0.5 * (z0 + z1) - rate * theta;
            let half = 0.5 * (z1 - z0) + rate.abs() * angular;
            if 2. * half >= spec.pitch_mm {
                self.radius_range()
            } else {
                let [radius, _] = self.radius_and_slope_at_phase(phase);
                let slope = self.profile_limits().0;
                [radius - slope * half, radius + slope * half]
            }
        } else {
            self.radius_range()
        };
        let [low, high] = [low - rounding, high + rounding];
        if rho_lo > high {
            Outside
        } else if rho_hi < low && z0 > rounding && z1 < spec.length_mm - rounding {
            Inside
        } else {
            UnresolvedBoundary
        }
    }

    /// Root and crest radii: the meridian's whole range.
    fn radius_range(&self) -> [f64; 2] {
        [self.root_radius_mm(), self.spec.major_radius_mm]
    }
}

impl ThroughBore {
    /// Classify a closed local box against the receiver: stock cylinder
    /// material outside the thread cavity, within the receiver's length.
    pub fn classify_box(&self, bounds: [[f64; 2]; 3]) -> BoxClassification {
        use BoxClassification::*;
        let cavity = self.thread().classify_box(bounds);
        if cavity == Inside {
            return Outside;
        }
        if bounds.iter().flatten().any(|value| !value.is_finite()) {
            return UnresolvedBoundary;
        }
        let [[x0, x1], [y0, y1], [z0, z1]] = bounds;
        let length = self.length_mm();
        let rounding = 64.
            * f64::EPSILON
            * bounds
                .iter()
                .flatten()
                .map(|value| value.abs())
                .fold(self.radius_mm().max(length).max(1.), f64::max);
        let near = |lo: f64, hi: f64| {
            if lo > 0. {
                lo
            } else if hi < 0. {
                -hi
            } else {
                0.
            }
        };
        let far = |lo: f64, hi: f64| lo.abs().max(hi.abs());
        let rho_lo = near(x0, x1).hypot(near(y0, y1)) * (1. - 4. * f64::EPSILON);
        let rho_hi = far(x0, x1).hypot(far(y0, y1)) * (1. + 4. * f64::EPSILON);
        let radius = self.radius_mm();
        if rho_lo > radius + rounding || z1 < -rounding || z0 > length + rounding {
            Outside
        } else if cavity == Outside
            && rho_hi < radius - rounding
            && z0 > rounding
            && z1 < length - rounding
        {
            Inside
        } else {
            UnresolvedBoundary
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::Spec;
    use super::*;
    use BoxClassification::*;

    fn thread(left_hand: bool) -> Thread {
        let angle = 27.5_f64.to_radians();
        Thread::new(Spec {
            major_radius_mm: 5.,
            pitch_mm: 1.,
            length_mm: 4.,
            depth_mm: 1. / (2. * angle.tan()) - 0.274 * (1. / angle.sin() - 1.),
            crest_radius_mm: 0.137,
            root_radius_mm: 0.137,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand,
        })
        .unwrap()
    }

    fn cube(center: [f64; 3], half: f64) -> [[f64; 2]; 3] {
        center.map(|c| [c - half, c + half])
    }

    #[test]
    fn core_crest_and_exterior_boxes_are_proved_for_both_hands() {
        for left in [false, true] {
            let t = thread(left);
            let root = t.root_radius_mm();
            // Well inside the core, beyond the crest, beyond each cap.
            assert_eq!(t.classify_box(cube([root - 0.5, 0., 2.], 0.05)), Inside);
            assert_eq!(t.classify_box(cube([0., 0., 2.], 0.05)), Inside);
            assert_eq!(t.classify_box(cube([5.5, 0., 2.], 0.05)), Outside);
            assert_eq!(t.classify_box(cube([1., 1., -1.], 0.05)), Outside);
            assert_eq!(t.classify_box(cube([1., 1., 5.], 0.05)), Outside);
            // At a crest phase the tooth is material just under the crest; at
            // a root phase the same radius is empty.
            let rate = t.helical_rate();
            let theta = 0.3_f64;
            let crest = rate * theta + 2.; // phase 0 at the crest centre
            let root_phase = crest + 0.5;
            let r = 4.95;
            let point = |z: f64| [r * theta.cos(), r * theta.sin(), z];
            assert_eq!(t.classify_box(cube(point(crest), 1e-4)), Inside);
            assert_eq!(t.classify_box(cube(point(root_phase), 1e-4)), Outside);
        }
    }

    #[test]
    fn receiver_material_lies_between_the_cavity_and_the_stock_wall() {
        let bore = ThroughBore::new(thread(false), 8.).unwrap();
        let rate = bore.thread().helical_rate();
        let theta = 0.3_f64;
        let point = |r: f64, phase: f64| {
            let z = rate * theta + 2. + phase;
            [r * theta.cos(), r * theta.sin(), z]
        };
        // Stock wall material, the empty core, beyond the wall and the caps.
        assert_eq!(bore.classify_box(cube(point(6.5, 0.), 0.05)), Inside);
        assert_eq!(bore.classify_box(cube([0., 0., 2.], 0.05)), Outside);
        assert_eq!(bore.classify_box(cube(point(9., 0.), 0.05)), Outside);
        assert_eq!(bore.classify_box(cube([6.5, 0., -1.], 0.05)), Outside);
        // Just outside the crest at a root phase is receiver tooth material;
        // the same radius at the crest phase is inside the cavity.
        assert_eq!(bore.classify_box(cube(point(4.95, 0.5), 1e-4)), Inside);
        assert_eq!(bore.classify_box(cube(point(4.95, 0.), 1e-4)), Outside);
        // On the stock wall, a cap, or invalid input.
        assert_eq!(
            bore.classify_box(cube(point(8., 0.), 0.01)),
            UnresolvedBoundary
        );
        assert_eq!(
            bore.classify_box(cube([6.5, 0., 4.], 0.01)),
            UnresolvedBoundary
        );
        assert_eq!(
            bore.classify_box([[f64::NAN, 1.], [0., 1.], [0., 1.]]),
            UnresolvedBoundary
        );
    }

    #[test]
    fn surface_cap_and_invalid_boxes_stay_unresolved() {
        let t = thread(false);
        // Straddling the crest cylinder, an end cap, and boxes too wide to
        // resolve the phase.
        assert_eq!(t.classify_box(cube([5., 0., 2.], 0.01)), UnresolvedBoundary);
        assert_eq!(t.classify_box(cube([1., 0., 0.], 0.01)), UnresolvedBoundary);
        assert_eq!(t.classify_box(cube([4.6, 0., 2.], 0.6)), UnresolvedBoundary);
        assert_eq!(
            t.classify_box([[f64::NAN, 1.], [0., 1.], [0., 1.]]),
            UnresolvedBoundary
        );
        assert_eq!(
            t.classify_box([[2., 1.], [0., 1.], [0., 1.]]),
            UnresolvedBoundary
        );
    }
}
