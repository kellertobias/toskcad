//! Exact shared-crest radial domain, independent of original-file ownership.
//! r(u,t)=(1-t)*thread_r(u)+t*head_r on bands0,1,2. The periodic
//! source has ordered real phases and 0<root<major; hence these three ruled
//! patches fill the head disk outside the thread star section. At nodes0/3
//! the radial interval collapses; constant-major crest band3 is omitted.
use super::*;
pub(super) struct Domain {
    pub source: PeriodicThreadSource,
    pub bottom: [f64; 2],
}
impl Domain {
    pub(super) fn new(
        source: PeriodicThreadSource,
        bottom: [f64; 2],
        stations: [f64; 2],
        work: &mut Work<'_>,
    ) -> FResult<Self> {
        work.step(64)?;
        // Equality is equality of defining real constants, not tolerance-based
        // projection. Raw topology and literal support transfer are separate.
        if source.radius_nodes_mm()[1] >= source.radius_nodes_mm()[0]
            || bottom[0] != source.radius_nodes_mm()[0]
            || bottom[1] != stations[0]
            || !stations[0].is_finite()
            || !stations[1].is_finite()
            || stations[0] < 0.
            || stations[1] <= stations[0]
            || stations[1] > source.length_mm()
        {
            return Err(PhillipsReadError::Unsupported);
        }
        Ok(Self { source, bottom })
    }
    pub(super) fn radial(&self, band: usize, u: I, t: I, work: &mut Work<'_>) -> FResult<I> {
        work.step(128)?;
        if band >= 3 || !u.valid() || u.lo < 0. || u.hi > 1. || !t.valid() || t.lo < 0. || t.hi > 1.
        {
            return Err(PhillipsReadError::Invalid);
        }
        Ok((I::point(1.) - t) * self.source.radii(band, u) + t * I::point(self.bottom[0]))
    }
}
