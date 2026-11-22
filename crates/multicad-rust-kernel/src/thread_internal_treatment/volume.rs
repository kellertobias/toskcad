//! Literal affine source provider and its independent angular-sliver correction.
use super::*;
use crate::tapered_thread::interval::I;
use crate::thread_cavity_volume::{self, Source};
impl Cavity {
    pub fn removed_volume_mm3(&self, options: QueryOptions) -> Result<Enclosure, Error> {
        self.removed_volume_with_interrupt(options, None)
    }
    pub fn removed_volume_with_interrupt(
        &self,
        options: QueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Enclosure, Error> {
        thread_cavity_volume::integrate(
            self,
            options,
            self.phase_advance_volume_reserve(),
            interrupt,
        )
    }
    pub(super) fn phase_advance_volume_reserve(&self) -> f64 {
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let p = I::point(self.thread.spec().pitch_mm);
        let lead = p * I::point(f64::from(self.thread.spec().starts));
        let true_rate = lead / (I::point(2.) * pi);
        let true_rate = if self.thread.spec().left_hand {
            -true_rate
        } else {
            true_rate
        };
        let evaluated = self.thread.helical_rate();
        let delta = (I::point(evaluated) - true_rate).magnitude();
        let knots = self.thread.profile_breaks();
        let mut radius = self.entry_radius_mm().unwrap_or(0.);
        for band in 0..5 {
            radius = radius.max(
                self.profile(band, I::new(knots[band], knots[band + 1]))
                    .0
                    .hi,
            );
        }
        // At each z the two phase integrals differ by their normalization
        // and a phase-domain sliver. Bounding both by R^2 gives this outward
        // (conservative factor-two) reserve, independent of source continuity.
        (I::point(2.)
            * pi
            * I::point(self.terminal_station)
            * I::point(radius).square()
            * I::point(delta)
            / I::point(evaluated.abs()))
        .hi
    }
    pub(super) fn profile(&self, band: usize, q: I) -> (I, I) {
        let spec = self.thread.spec();
        let knots = self.thread.profile_breaks();
        if band == 0 || band == 4 {
            return (I::point(spec.major_radius_mm), I::point(0.));
        }
        if band == 2 {
            return (I::point(self.thread.root_radius_mm()), I::point(0.));
        }
        let anchor = if band == 1 { knots[1] } else { knots[4] };
        let slope = I::point(self.thread.slope_at_phase(anchor, band));
        (
            I::point(spec.major_radius_mm) + slope * (q - I::point(anchor)),
            slope,
        )
    }
}
impl Source for Cavity {
    fn phase_knots(&self) -> [f64; 6] {
        self.thread.profile_breaks()
    }
    fn profile(&self, band: usize, q: I) -> (I, I) {
        Cavity::profile(self, band, q)
    }
    fn pitch_mm(&self) -> f64 {
        self.thread.spec().pitch_mm
    }
    fn full_length_mm(&self) -> f64 {
        Cavity::full_length_mm(self)
    }
    fn prefix_length_mm(&self) -> f64 {
        self.terminal_station
    }
    fn terminal_radius_mm(&self) -> Option<f64> {
        self.intent.terminal.map(|_| self.terminal_radius)
    }
    fn entry(&self) -> Option<thread_cavity_volume::Entry> {
        self.intent
            .entry_depth_mm
            .map(|depth| thread_cavity_volume::Entry {
                radius: self.entry_radius,
                depth,
                slope: self.entry_slope.expect("validated enabled cone"),
            })
    }
}
