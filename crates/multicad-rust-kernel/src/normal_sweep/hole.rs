//! Add a straight circular through hole to an unpierced normal sweep.
use super::{Edge, Error, Sweep, inside};

impl Sweep {
    /// The same sweep with one more clockwise circular hole ring. The circle
    /// must lie strictly inside the material: its centre is material and no
    /// existing profile edge comes within `radius + rounding` of it, so the
    /// new wall touches no other wall or hole.
    pub fn with_circular_hole(
        &self,
        center: [f64; 2],
        radius_mm: f64,
        tolerance_mm: f64,
        minimum_mm: f64,
    ) -> Result<Self, Error> {
        if !center.iter().all(|value| value.is_finite())
            || !radius_mm.is_finite()
            || radius_mm <= 0.
        {
            return Err(Error::InvalidProfile);
        }
        let magnitude = center
            .into_iter()
            .map(f64::abs)
            .fold(radius_mm.max(1.), f64::max);
        let clearance = radius_mm + tolerance_mm.max(64. * f64::EPSILON * magnitude);
        let material = inside(&self.rings[0], center)
            && self.rings[1..].iter().all(|ring| !inside(ring, center));
        if !material
            || self.rings.iter().flatten().any(|edge| {
                edge.distance(center).partial_cmp(&clearance) != Some(std::cmp::Ordering::Greater)
            })
        {
            return Err(Error::InvalidProfile);
        }
        let mut rings = self.rings.clone();
        rings.push(vec![Edge::Circle {
            center,
            radius_mm,
            clockwise: true,
        }]);
        Self::new(rings, self.height_mm, tolerance_mm, minimum_mm)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plate() -> Sweep {
        let corners = [[0., 0.], [10., 0.], [10., 6.], [0., 6.]];
        let outer = (0..4)
            .map(|i| Edge::Line {
                a: corners[i],
                b: corners[(i + 1) % 4],
            })
            .collect();
        let hole = vec![Edge::Circle {
            center: [3., 3.],
            radius_mm: 1.,
            clockwise: true,
        }];
        Sweep::new(vec![outer, hole], 2., 1e-7, 1e-5).unwrap()
    }

    #[test]
    fn a_clear_hole_removes_exactly_its_disc() {
        let source = plate();
        let drilled = source
            .with_circular_hole([7., 3.], 1.5, 1e-7, 1e-5)
            .unwrap();
        assert_eq!(drilled.rings().len(), 3);
        let disc = std::f64::consts::PI * 2.25;
        assert!((source.area_mm2() - drilled.area_mm2() - disc).abs() < 1e-9);
        assert!((source.volume_mm3() - drilled.volume_mm3() - 2. * disc).abs() < 1e-9);
        assert!(!drilled.contains([7., 3., 1.], 0.));
        assert!(drilled.contains([9., 5.5, 1.], 0.));
    }

    #[test]
    fn holes_reaching_walls_or_other_holes_or_outside_are_refused() {
        let source = plate();
        for (center, radius) in [
            ([9., 3.], 1.),  // tangent to the outer wall
            ([5., 3.], 1.),  // tangent to the existing hole
            ([3., 3.], 0.5), // inside the existing hole
            ([12., 3.], 1.), // outside the profile
            ([7., 3.], 0.),  // degenerate
            ([f64::NAN, 3.], 1.),
        ] {
            assert_eq!(
                source.with_circular_hole(center, radius, 1e-7, 1e-5),
                Err(Error::InvalidProfile),
                "{center:?} {radius}"
            );
        }
    }
}
