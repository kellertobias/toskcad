//! Direct normalized-source integration; no legacy cavity proxy or rate reserve.
use super::*;
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
        thread_cavity_volume::integrate(self, options, 0., interrupt)
    }
}
impl Source for Cavity {
    fn phase_knots(&self) -> [f64; 6] {
        let bands = self.thread.bands();
        [
            bands[0].phase_mm[0],
            bands[0].phase_mm[1],
            bands[1].phase_mm[1],
            bands[2].phase_mm[1],
            bands[3].phase_mm[1],
            bands[4].phase_mm[1],
        ]
    }
    fn profile(&self, band: usize, q: I) -> (I, I) {
        let b = self.thread.bands()[band];
        if b.radius_mm[0] == b.radius_mm[1] {
            return (I::point(b.radius_mm[0]), I::point(0.));
        }
        let slope = (I::point(b.radius_mm[1]) - I::point(b.radius_mm[0]))
            / (I::point(b.phase_mm[1]) - I::point(b.phase_mm[0]));
        let r = I::point(b.radius_mm[0]) + slope * (q - I::point(b.phase_mm[0]));
        // The actual source is convex interpolation of its shared endpoints.
        // Intersecting with that independent exact range loses no material.
        (
            I::new(
                r.lo.max(b.radius_mm[0].min(b.radius_mm[1])),
                r.hi.min(b.radius_mm[0].max(b.radius_mm[1])),
            ),
            slope,
        )
    }
    fn pitch_mm(&self) -> f64 {
        self.thread.spec().pitch_mm
    }
    fn full_length_mm(&self) -> f64 {
        Cavity::full_length_mm(self)
    }
    fn full_length_interval(&self) -> I {
        self.end_interval()
    }
    fn prefix_length_mm(&self) -> f64 {
        self.tail.map_or(self.full_length_mm(), |t| t.station)
    }
    fn prefix_length_interval(&self) -> I {
        self.tail
            .map_or(self.end_interval(), |t| I::point(t.station))
    }
    fn terminal_radius_mm(&self) -> Option<f64> {
        Cavity::terminal_radius_mm(self)
    }
    fn entry(&self) -> Option<thread_cavity_volume::Entry> {
        self.entry.map(|e| thread_cavity_volume::Entry {
            radius: e.radius,
            depth: e.depth,
            slope: e.slope,
        })
    }
}
