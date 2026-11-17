//! Certified radial clearance between an external thread and a coaxial
//! threaded cavity of the same lead and hand.
//!
//! Both meridians are periodic functions of the helical phase alone, so on a
//! common axis the material of the external thread (`rho <= outer(phase)`)
//! and the material around the cavity (`rho >= cavity(phase - shift)`) share
//! no interior exactly when `cavity(phase - shift) - outer(phase)` is never
//! negative. That one-dimensional bound is proved over a whole period by
//! subdivision with the meridians' certified slope limits; no sample is
//! trusted without its Lipschitz margin.
use super::{Error, Thread};

/// Cells admitted before an unresolved (near-touching) clearance is reported.
const MAXIMUM_CELLS: usize = 1 << 16;

impl Thread {
    /// Upper bound of the meridian's slope magnitude over a whole period.
    pub fn meridian_slope_limit(&self) -> f64 {
        self.profile_limits().0
    }

    /// Prove that this external thread lies within `cavity`, with the cavity
    /// phase shifted by `shift_mm` along the axis. Returns a certified lower
    /// bound on the radial clearance, which may be zero within rounding.
    ///
    /// `Invalid` reports incompatible leads, hands or starts; `Unresolved`
    /// reports measurable overlap or contact that cannot be separated from
    /// overlap within the cell budget.
    pub fn prove_clearance(&self, cavity: &Thread, shift_mm: f64) -> Result<f64, Error> {
        let (a, b) = (self.spec, cavity.spec);
        if !shift_mm.is_finite()
            || a.pitch_mm != b.pitch_mm
            || a.starts != b.starts
            || a.left_hand != b.left_hand
        {
            return Err(Error::Invalid);
        }
        let pitch = a.pitch_mm;
        let lipschitz = self.profile_limits().0 + cavity.profile_limits().0;
        if !lipschitz.is_finite() {
            return Err(Error::Unresolved);
        }
        let shift = shift_mm.rem_euclid(pitch);
        // Rounding of each meridian evaluation and of the phase arithmetic.
        let rounding = 64.
            * f64::EPSILON
            * (a.major_radius_mm + b.major_radius_mm + lipschitz * (pitch + shift_mm.abs()));
        let gap = |phase: f64| cavity.radius_at_phase(phase - shift) - self.radius_at_phase(phase);
        let mut lower = f64::INFINITY;
        let mut cells = vec![(0., pitch)];
        let mut visited = 0usize;
        while let Some((start, end)) = cells.pop() {
            visited += 1;
            if visited > MAXIMUM_CELLS {
                return Err(Error::Unresolved);
            }
            let middle = 0.5 * (start + end);
            let value = gap(middle);
            if !value.is_finite() || value < -rounding {
                return Err(Error::Unresolved);
            }
            let bound = value - 0.5 * (end - start) * lipschitz - rounding;
            if bound >= -rounding {
                lower = lower.min(bound.max(0.));
            } else if end - start <= 64. * f64::EPSILON * pitch {
                return Err(Error::Unresolved);
            } else {
                cells.push((start, middle));
                cells.push((middle, end));
            }
        }
        Ok(lower)
    }
}

#[cfg(test)]
mod tests {
    use super::super::Spec;
    use super::*;

    /// A Whitworth-like 55 degree form scaled by its major radius.
    fn thread(major: f64, depth: f64, crest: f64, root: f64) -> Thread {
        Thread::new(spec(major, depth, crest, root)).unwrap()
    }

    fn spec(major: f64, depth: f64, crest: f64, root: f64) -> Spec {
        Spec {
            major_radius_mm: major,
            pitch_mm: 1.,
            length_mm: 4.,
            depth_mm: depth,
            crest_radius_mm: crest,
            root_radius_mm: root,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand: false,
        }
    }

    /// Depth consistent with the pitch, flank angle and roundings.
    fn depth(crest: f64, root: f64) -> f64 {
        let angle = 27.5_f64.to_radians();
        1. / (2. * angle.tan()) - (crest + root) * (1. / angle.sin() - 1.)
    }

    #[test]
    fn an_enlarged_cavity_clears_unless_its_phase_is_shifted() {
        let d = depth(0.137, 0.137);
        let external = thread(5., d, 0.137, 0.137);
        let cavity = thread(5.05, d, 0.137, 0.137);
        let clearance = external.prove_clearance(&cavity, 0.).unwrap();
        // A certified lower bound of the uniform 0.05 mm radial clearance.
        assert!(clearance > 0. && clearance <= 0.05 + 1e-12, "{clearance}");
        // A whole period shift is the same helix.
        let shifted = external.prove_clearance(&cavity, 3.).unwrap();
        assert!(shifted > 0. && shifted <= 0.05 + 1e-12);
        // A half-period shift puts crest into crest.
        assert_eq!(
            external.prove_clearance(&cavity, 0.5),
            Err(Error::Unresolved)
        );
        // The cavity must be outside, not inside.
        assert_eq!(
            cavity.prove_clearance(&external, 0.),
            Err(Error::Unresolved)
        );
    }

    #[test]
    fn incompatible_helices_are_not_compared() {
        let d = depth(0.137, 0.137);
        let external = thread(5., d, 0.137, 0.137);
        let mut left = spec(5.05, d, 0.137, 0.137);
        left.left_hand = true;
        let left = Thread::new(left).unwrap();
        assert_eq!(external.prove_clearance(&left, 0.), Err(Error::Invalid));
        let mut coarse = spec(5.05, 0., 0.137, 0.137);
        coarse.pitch_mm = 2.;
        let angle = 27.5_f64.to_radians();
        coarse.depth_mm = 2. / (2. * angle.tan()) - 0.274 * (1. / angle.sin() - 1.);
        let coarse = Thread::new(coarse).unwrap();
        assert_eq!(external.prove_clearance(&coarse, 0.), Err(Error::Invalid));
        let cavity = thread(5.05, d, 0.137, 0.137);
        assert_eq!(
            external.prove_clearance(&cavity, f64::NAN),
            Err(Error::Invalid)
        );
    }
}
