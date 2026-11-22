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
/// Relative pitch difference admitted as decode roundoff by
/// `prove_clearance_over`; larger differences are distinct threads.
const MAXIMUM_PITCH_DRIFT: f64 = 1e-9;

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
        self.gap_lower_bound(cavity, shift_mm)
    }

    /// `prove_clearance` for decoded cavities whose pitch differs from this
    /// thread's by at most `MAXIMUM_PITCH_DRIFT` relative roundoff, over this
    /// thread's finite `external_length_mm`.
    ///
    /// The proof runs in this thread's period. Writing an axial station as
    /// `z = k * pitch + r` with `0 <= k <= ceil(length / pitch)`, the cavity
    /// argument differs from the period-reduced one by `(k - m) * delta`,
    /// where `m` folds the shift, plus at most `starts * delta` from the
    /// differing helical rates over one revolution. That drift is charged at
    /// the cavity meridian's certified slope limit; a negative remainder is
    /// `Unresolved`. Equal pitches reduce to `prove_clearance` exactly.
    pub fn prove_clearance_over(
        &self,
        cavity: &Thread,
        shift_mm: f64,
        external_length_mm: f64,
    ) -> Result<f64, Error> {
        let (a, b) = (self.spec, cavity.spec);
        if a.pitch_mm == b.pitch_mm {
            return self.prove_clearance(cavity, shift_mm);
        }
        let delta = (b.pitch_mm - a.pitch_mm).abs();
        if !shift_mm.is_finite()
            || !external_length_mm.is_finite()
            || external_length_mm < 0.
            || delta.is_nan()
            || delta > MAXIMUM_PITCH_DRIFT * a.pitch_mm
            || a.starts != b.starts
            || a.left_hand != b.left_hand
        {
            return Err(Error::Invalid);
        }
        let turns = (external_length_mm / a.pitch_mm).ceil()
            + (shift_mm / a.pitch_mm).floor().abs()
            + f64::from(a.starts)
            + 1.;
        let drift = (turns * delta * (1. + 16. * f64::EPSILON)).next_up();
        let slope = cavity.profile_limits().0;
        if !slope.is_finite() || !drift.is_finite() {
            return Err(Error::Unresolved);
        }
        let lower = self.gap_lower_bound(cavity, shift_mm)?;
        let remaining = lower - (drift * slope * (1. + 16. * f64::EPSILON)).next_up();
        if remaining >= 0. {
            Ok(remaining)
        } else {
            Err(Error::Unresolved)
        }
    }

    fn gap_lower_bound(&self, cavity: &Thread, shift_mm: f64) -> Result<f64, Error> {
        let (a, b) = (self.spec, cavity.spec);
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
    use multicad_precision::PortableMath as _;

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
        1. / (2. * angle.portable_tan()) - (crest + root) * (1. / angle.portable_sin() - 1.)
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
    fn roundoff_pitch_drift_is_charged_over_the_external_length() {
        let d = depth(0.137, 0.137);
        let external = thread(5., d, 0.137, 0.137);
        let mut decoded = spec(5.05, d, 0.137, 0.137);
        decoded.pitch_mm = 1_f64.next_down();
        decoded.depth_mm = d.next_down();
        let decoded = Thread::new(decoded).unwrap();
        // Unequal leads are never compared by the equal-pitch proof.
        assert_eq!(external.prove_clearance(&decoded, 0.), Err(Error::Invalid));
        let exact = external
            .prove_clearance(&thread(5.05, d, 0.137, 0.137), 0.)
            .unwrap();
        let drifted = external.prove_clearance_over(&decoded, 0., 4.).unwrap();
        assert!(drifted > 0. && drifted < exact, "{drifted} vs {exact}");
        // Equal pitches take the exact path unchanged.
        let cavity = thread(5.05, d, 0.137, 0.137);
        assert_eq!(
            external.prove_clearance_over(&cavity, 0., 4.),
            external.prove_clearance(&cavity, 0.)
        );
        // A half-period shift still overlaps.
        assert_eq!(
            external.prove_clearance_over(&decoded, 0.5, 4.),
            Err(Error::Unresolved)
        );
        // A real pitch difference is a different thread, not roundoff.
        let closed = |pitch: f64| {
            let angle = 27.5_f64.to_radians();
            pitch / (2. * angle.portable_tan()) - 0.274 * (1. / angle.portable_sin() - 1.)
        };
        let mut other = spec(5.05, closed(1.0001), 0.137, 0.137);
        other.pitch_mm = 1.0001;
        let other = Thread::new(other).unwrap();
        assert_eq!(
            external.prove_clearance_over(&other, 0., 4.),
            Err(Error::Invalid)
        );
        // Over a huge length the charged drift consumes the clearance.
        let mut coarse = spec(5.05, closed(1. + 5e-10), 0.137, 0.137);
        coarse.pitch_mm = 1. + 5e-10;
        let coarse = Thread::new(coarse).unwrap();
        assert_eq!(
            external.prove_clearance_over(&coarse, 0., 1e9),
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
        coarse.depth_mm =
            2. / (2. * angle.portable_tan()) - 0.274 * (1. / angle.portable_sin() - 1.);
        let coarse = Thread::new(coarse).unwrap();
        assert_eq!(external.prove_clearance(&coarse, 0.), Err(Error::Invalid));
        let cavity = thread(5.05, d, 0.137, 0.137);
        assert_eq!(
            external.prove_clearance(&cavity, f64::NAN),
            Err(Error::Invalid)
        );
    }
}
