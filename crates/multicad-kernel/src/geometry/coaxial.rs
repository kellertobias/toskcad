//! Relative placement of two charts that share one axis direction.
//!
//! Helical sources are periodic in the phase `z - rate * theta`, so a second
//! chart rotated about or translated along the first chart's axis only shifts
//! that phase. Rounding of the placements is carried as explicit enclosures;
//! nothing is snapped to a nearby exact placement.
use super::*;

#[derive(Clone, Copy, Debug)]
pub(in crate::geometry) struct Coaxial {
    /// Signed distance of the second origin along the first axis.
    pub(in crate::geometry) axial_mm: f64,
    /// Rotation of the second chart's u axis about the first axis.
    pub(in crate::geometry) angle_rad: f64,
    /// Largest radial displacement between the two axes over `span_mm`.
    pub(in crate::geometry) radial_mm: f64,
    /// Rounding allowance for `axial_mm` and `angle_rad` (as a length).
    pub(in crate::geometry) error_mm: f64,
}

impl Coaxial {
    /// Charts with the same axis direction up to rounding, whose axes stay
    /// within a rounding-level distance of each other over `span_mm`.
    pub(in crate::geometry) fn of(a: PlaneFrame, b: PlaneFrame, span_mm: f64) -> Option<Self> {
        let [na, nb] = [a.normal(), b.normal()];
        let tilt = dot(cross(na, nb), cross(na, nb)).sqrt();
        let magnitude = a
            .origin
            .into_iter()
            .chain(b.origin)
            .map(f64::abs)
            .fold(span_mm.abs().max(1.), f64::max);
        let rounding = 64. * f64::EPSILON * magnitude;
        let parallel = dot(na, nb);
        if parallel.is_nan()
            || parallel <= 0.
            || tilt.is_nan()
            || tilt > 64. * f64::EPSILON
            || !span_mm.is_finite()
        {
            return None;
        }
        let delta: [f64; 3] = std::array::from_fn(|i| b.origin[i] - a.origin[i]);
        let axial_mm = dot(delta, na);
        let offset: [f64; 3] = std::array::from_fn(|i| delta[i] - na[i] * axial_mm);
        let radial_mm = dot(offset, offset).sqrt() + tilt * (span_mm.abs() + axial_mm.abs());
        if radial_mm.is_nan() || radial_mm > rounding {
            return None;
        }
        let angle_rad = dot(b.u, a.v).atan2(dot(b.u, a.u));
        Some(Self {
            axial_mm,
            angle_rad,
            radial_mm,
            error_mm: rounding,
        })
    }

    /// Enclosure of the phase shift `axial - rate * angle` for a helix of
    /// signed helical rate `rate` (lead per radian) whose cavity lies at a
    /// radius of at least `radius_mm` from its axis.
    pub(in crate::geometry) fn phase_shift(self, rate: f64, radius_mm: f64) -> Option<[f64; 2]> {
        let shift = self.axial_mm - rate * self.angle_rad;
        // Axis displacement changes the angle seen from the other axis by at
        // most radial / radius; the angle itself carries a few ulps.
        let angular = self.radial_mm / radius_mm + 16. * f64::EPSILON * std::f64::consts::PI;
        let error = self.error_mm + rate.abs() * angular + 16. * f64::EPSILON * shift.abs();
        (shift.is_finite() && error.is_finite() && radius_mm > 0.)
            .then_some([shift - error, shift + error])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn rotation_about_and_travel_along_the_axis_become_a_phase_shift() {
        let a = PlaneFrame {
            origin: [10., 10., 9.7],
            u: [1., 0., 0.],
            v: [0., -1., 0.],
        };
        let b = PlaneFrame {
            origin: [10., 10., 0.],
            u: [0., -1., 0.],
            v: [-1., 0., 0.],
        };
        let placement = Coaxial::of(a, b, 30.).unwrap();
        // The first axis is -Z, so the second origin lies 9.7 along it.
        assert!((placement.axial_mm - 9.7).abs() < 1e-12);
        assert!((placement.angle_rad - FRAC_PI_2).abs() < 1e-12);
        assert_eq!(placement.radial_mm, 0.);
        let [lo, hi] = placement.phase_shift(0.5, 1.).unwrap();
        let expected = 9.7 - 0.5 * FRAC_PI_2;
        assert!(lo <= expected && expected <= hi && hi - lo < 1e-9);
        assert!(placement.phase_shift(0.5, 0.).is_none());
    }

    #[test]
    fn opposed_tilted_or_offset_axes_are_not_coaxial() {
        let a = PlaneFrame::default();
        let flipped = PlaneFrame {
            v: [0., -1., 0.],
            ..a
        };
        assert!(Coaxial::of(a, flipped, 1.).is_none());
        let tilted = PlaneFrame {
            u: [1., 0., 0.],
            v: [0., 1., 1e-9].map(|x| x / (1f64 + 1e-18).sqrt()),
            ..a
        };
        assert!(Coaxial::of(a, tilted, 1.).is_none());
        let offset = PlaneFrame {
            origin: [1e-6, 0., 0.],
            ..a
        };
        assert!(Coaxial::of(a, offset, 1.).is_none());
        assert!(Coaxial::of(a, a, f64::NAN).is_none());
    }
}
