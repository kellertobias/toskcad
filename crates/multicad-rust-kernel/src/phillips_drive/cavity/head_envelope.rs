//! Whole-source containment in a conservative retained-profile support envelope.
//! This is not a closed outer BRep, raw material or mass certificate.
use super::*;
use crate::head_profile::Profile;
#[cfg(test)]
mod tests;
const MAX_CELLS: usize = 65536;
const MAX_WORK: usize = 16_000_000;
const MAX_BYTES: usize = 64 * 1024 * 1024;
#[derive(Clone, Copy)]
pub struct HeadEnvelopeOptions {
    pub minimum_clearance_mm: f64,
    pub maximum_cells: usize,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for HeadEnvelopeOptions {
    fn default() -> Self {
        Self {
            minimum_clearance_mm: 1e-7,
            maximum_cells: MAX_CELLS,
            maximum_work: MAX_WORK,
            maximum_bytes: MAX_BYTES,
        }
    }
}
/// Source-chart radial/axial margins to the supported head envelope. They are
/// not Euclidean normal-offset wall thickness under an arbitrary affine frame.
pub struct HeadEnvelopeClearance<'s, 'p> {
    source: &'s Cutter,
    profile: &'p Profile,
    cells: usize,
    radial: f64,
    floor: f64,
    work: usize,
}
impl<'s, 'p> HeadEnvelopeClearance<'s, 'p> {
    pub fn source(&self) -> &'s Cutter {
        self.source
    }
    pub fn profile_proposal(&self) -> &'p Profile {
        self.profile
    }
    pub fn cells(&self) -> usize {
        self.cells
    }
    pub fn minimum_radial_clearance_mm(&self) -> f64 {
        self.radial
    }
    pub fn axial_floor_clearance_mm(&self) -> f64 {
        self.floor
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + std::mem::size_of::<Cutter>() + std::mem::size_of::<Profile>()
    }
    pub(crate) fn scratch_budget_bytes(cells: usize) -> usize {
        4 * (cells + 1) * std::mem::size_of::<[f64; 2]>() + 8192
    }
}
impl<'s> ClosedBoundaryCertificate<'s> {
    /// Evaluate every depth using the actual line/circle supports and the
    /// profile's separately declared reconstruction reserve. No volume cache or
    /// guessed head catalog is constructed. The supplied profile stays explicitly
    /// a proposal; a raw-file caller must earn its own whole finite matching.
    pub fn certify_head_envelope<'p>(
        &self,
        profile: &'p Profile,
        options: HeadEnvelopeOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HeadEnvelopeClearance<'s, 'p>, Error> {
        let checkpoint = || {
            if interrupt.is_some_and(|f| f()) {
                Err(Error::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        if !options.minimum_clearance_mm.is_finite()
            || options.minimum_clearance_mm <= 0.
            || options.maximum_cells == 0
            || options.maximum_cells > MAX_CELLS
            || options.maximum_work > MAX_WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > MAX_BYTES
        {
            return Err(Error::Invalid);
        }
        let mut remaining = options
            .maximum_work
            .checked_sub(self.work_used())
            .ok_or(Error::Limit)?;
        let peak = std::mem::size_of::<HeadEnvelopeClearance<'_, '_>>()
            + std::mem::size_of::<Cutter>()
            + std::mem::size_of::<Profile>()
            + HeadEnvelopeClearance::scratch_budget_bytes(options.maximum_cells);
        if peak > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let source = self.source();
        let reserve = I::point(2.) * I::point(profile.rounding_error_mm());
        let top = I::point(profile.top()[1]);
        let tip = source.envelope.tip;
        let floor = (I::point(profile.end()[1]) - top - tip - reserve).lo;
        if floor <= options.minimum_clearance_mm
            || (I::point(profile.top()[0]) - source.envelope.entry - reserve).lo
                <= options.minimum_clearance_mm
        {
            return Err(Error::Unresolved);
        }
        let mut pending = vec![[0., tip.hi]];
        let mut cells = 0;
        let mut radial = f64::INFINITY;
        while let Some([a, b]) = pending.pop() {
            checkpoint()?;
            if cells == options.maximum_cells {
                return Err(Error::Limit);
            }
            // At most five retained meridian segments (<=10 interval operations
            // each) and three analytic suffix pieces (<20), plus fixed range,
            // gap and split work. Charge the full bounded evaluation, not samples
            // or measured CPU time; no hidden profile/volume traversal occurs.
            remaining = remaining.checked_sub(256).ok_or(Error::Limit)?;
            cells += 1;
            let lower = host::head_radius_lower(profile, top + I::new(a, b));
            let radius = source.radial_suffix_bound_mm(a)?;
            let gap = (I::point(lower) - I::point(radius) - reserve).lo;
            if gap > options.minimum_clearance_mm {
                radial = radial.min(gap)
            } else {
                let mid = a + (b - a) / 2.;
                if mid <= a || mid >= b {
                    return Err(Error::Unresolved);
                }
                pending.push([mid, b]);
                pending.push([a, mid]);
            }
        }
        checkpoint()?;
        remaining = remaining.checked_sub(1).ok_or(Error::Limit)?;
        Ok(HeadEnvelopeClearance {
            source,
            profile,
            cells,
            radial,
            floor,
            work: options.maximum_work - remaining,
        })
    }
}
