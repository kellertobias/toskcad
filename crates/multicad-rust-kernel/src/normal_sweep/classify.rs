//! Whole-box material classification for an unpierced normal sweep.
//!
//! This is a proof, not a distance query: a box is `Inside` or `Outside` only
//! when no boundary of the actual profile or caps can reach it. Anything near
//! a wall, cap or profile edge stays explicitly unresolved.
use super::{Sweep, inside};
use crate::helical_thread::BoxClassification;

impl Sweep {
    /// Classify the closed local box `[[x0, x1], [y0, y1], [z0, z1]]` given in
    /// the sweep chart (profile in x/y, material between z = 0 and the height).
    pub fn classify_box(&self, bounds: [[f64; 2]; 3]) -> BoxClassification {
        use BoxClassification::*;
        if bounds.iter().flatten().any(|value| !value.is_finite())
            || bounds.iter().any(|[lo, hi]| lo > hi)
        {
            return UnresolvedBoundary;
        }
        let [[x0, x1], [y0, y1], [z0, z1]] = bounds;
        let height = self.height_mm;
        let magnitude = bounds
            .iter()
            .flatten()
            .map(|value| value.abs())
            .fold(height.max(1.), f64::max);
        // Rounding of the analytic edge distances and of the chart itself.
        let margin = 64. * f64::EPSILON * magnitude;
        if z1 < -margin || z0 > height + margin {
            return Outside;
        }
        let center = [0.5 * (x0 + x1), 0.5 * (y0 + y1)];
        // A disc around the centre that encloses the whole rectangle.
        let reach = 0.5 * (x1 - x0).hypot(y1 - y0) * (1. + 4. * f64::EPSILON) + margin;
        if self.rings.iter().flatten().any(|edge| {
            edge.distance(center).partial_cmp(&reach) != Some(std::cmp::Ordering::Greater)
        }) {
            return UnresolvedBoundary;
        }
        // No profile edge reaches the rectangle, so it lies wholly on one side.
        let material = inside(&self.rings[0], center)
            && self.rings[1..].iter().all(|ring| !inside(ring, center));
        if !material {
            Outside
        } else if z0 > margin && z1 < height - margin {
            Inside
        } else {
            UnresolvedBoundary
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normal_sweep::Edge;
    use BoxClassification::*;

    /// A 10 x 10 square with a radius-2 hole at its centre, 5 high.
    fn sweep() -> Sweep {
        let corners = [[0., 0.], [10., 0.], [10., 10.], [0., 10.]];
        let outer = (0..4)
            .map(|i| Edge::Line {
                a: corners[i],
                b: corners[(i + 1) % 4],
            })
            .collect();
        let hole = vec![Edge::Circle {
            center: [5., 5.],
            radius_mm: 2.,
            clockwise: true,
        }];
        Sweep::new(vec![outer, hole], 5., 1e-7, 1e-5).unwrap()
    }

    fn cube(center: [f64; 3], half: f64) -> [[f64; 2]; 3] {
        center.map(|c| [c - half, c + half])
    }

    #[test]
    fn boxes_away_from_every_boundary_are_proved() {
        let s = sweep();
        assert_eq!(s.classify_box(cube([1., 1., 2.5], 0.1)), Inside);
        assert_eq!(s.classify_box(cube([8.5, 5., 1.], 0.4)), Inside);
        // In the hole, beside the profile, below and above the caps.
        assert_eq!(s.classify_box(cube([5., 5., 2.5], 0.5)), Outside);
        assert_eq!(s.classify_box(cube([12., 5., 2.5], 0.5)), Outside);
        assert_eq!(s.classify_box(cube([5., 1., -1.], 0.5)), Outside);
        assert_eq!(s.classify_box(cube([5., 1., 6.], 0.5)), Outside);
        // Outside the profile but straddling a cap plane is still outside.
        assert_eq!(s.classify_box(cube([12., 5., 5.], 0.5)), Outside);
    }

    #[test]
    fn boxes_touching_walls_caps_or_hole_stay_unresolved() {
        let s = sweep();
        // Crossing an outer wall, the hole wall, and each cap.
        assert_eq!(
            s.classify_box(cube([10., 5., 2.5], 0.1)),
            UnresolvedBoundary
        );
        assert_eq!(s.classify_box(cube([7., 5., 2.5], 0.1)), UnresolvedBoundary);
        assert_eq!(s.classify_box(cube([1., 1., 0.], 0.1)), UnresolvedBoundary);
        assert_eq!(s.classify_box(cube([1., 1., 5.], 0.1)), UnresolvedBoundary);
        // The rectangle corner reaches the hole even though the centre is clear.
        assert_eq!(s.classify_box(cube([7.2, 7.2, 2.5], 0.2)), Inside);
        assert_eq!(
            s.classify_box(cube([6.3, 6.3, 2.5], 0.6)),
            UnresolvedBoundary
        );
        // Invalid or inverted input is never classified.
        assert_eq!(
            s.classify_box([[f64::NAN, 1.], [0., 1.], [0., 1.]]),
            UnresolvedBoundary
        );
        assert_eq!(
            s.classify_box([[2., 1.], [1., 2.], [1., 2.]]),
            UnresolvedBoundary
        );
    }
}
