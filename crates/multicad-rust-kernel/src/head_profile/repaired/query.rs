//! Whole-box queries of the retained monotone meridian. Candidate parameter
//! cells are excluded only by outward source enclosures, never a sampled fit.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RepairedHeadLocation {
    Inside,
    Outside,
    Unresolved,
}
#[derive(Clone, Copy)]
pub struct RepairedHeadQueryOptions {
    pub maximum_work: usize,
    pub maximum_cells: usize,
}
impl Default for RepairedHeadQueryOptions {
    fn default() -> Self {
        Self {
            maximum_work: WORK,
            maximum_cells: 4096,
        }
    }
}
/// A predicate result for this exact source and input coordinate box. This
/// does not prove a mesh, manufacturing boundary, or Boolean fragment.
pub struct RepairedHeadQuery<'q, 's> {
    source: &'q RepairedHeadMeridian<'s>,
    bounds: [[f64; 2]; 3],
    location: RepairedHeadLocation,
    used: usize,
    remaining: usize,
}
impl<'s> RepairedHeadQuery<'_, 's> {
    pub fn source(&self) -> &RepairedHeadMeridian<'s> {
        self.source
    }
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        self.bounds
    }
    pub fn location(&self) -> RepairedHeadLocation {
        self.location
    }
    pub fn work_used(&self) -> usize {
        self.used
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
}
struct QueryWork<'a> {
    remaining: usize,
    cells: usize,
    maximum_cells: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl QueryWork<'_> {
    fn charge(&mut self, amount: usize) -> Result<(), RepairedHeadError> {
        poll(self.interrupt)?;
        self.remaining = self
            .remaining
            .checked_sub(amount)
            .ok_or(RepairedHeadError::Limit)?;
        Ok(())
    }
    fn cell(&mut self) -> Result<(), RepairedHeadError> {
        self.charge(512)?;
        self.cells += 1;
        if self.cells > self.maximum_cells {
            return Err(RepairedHeadError::Limit);
        }
        Ok(())
    }
}
fn interval(v: [f64; 2]) -> I {
    I::new(v[0], v[1])
}
impl<'s> RepairedHeadMeridian<'s> {
    pub fn classify_box<'q>(
        &'q self,
        bounds: [[f64; 2]; 3],
        options: RepairedHeadQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RepairedHeadQuery<'q, 's>, RepairedHeadError> {
        let mut w = QueryWork {
            remaining: options.maximum_work,
            cells: 0,
            maximum_cells: options.maximum_cells,
            interrupt,
        };
        w.charge(1)?;
        if options.maximum_work > WORK
            || options.maximum_cells == 0
            || options.maximum_cells > 1_000_000
            || bounds.into_iter().map(interval).any(|v| !v.valid())
        {
            return Err(RepairedHeadError::Invalid);
        }
        let location = self.locate(bounds.map(interval), &mut w)?;
        w.charge(1)?;
        Ok(RepairedHeadQuery {
            source: self,
            bounds,
            location,
            used: options.maximum_work - w.remaining,
            remaining: w.remaining,
        })
    }
    fn locate(
        &self,
        p: [I; 3],
        w: &mut QueryWork<'_>,
    ) -> Result<RepairedHeadLocation, RepairedHeadError> {
        use RepairedHeadLocation::*;
        let [x, y, z] = p;
        let start = self.profile.top()[1];
        let end = self.profile.end()[1];
        if z.hi < start || z.lo > end {
            return Ok(Outside);
        }
        let mut radial_upper = 0_f64;
        for segment in 0..self.count {
            w.cell()?;
            let b = self.certified_interval_bounds(segment, [0., 1.]);
            radial_upper = radial_upper.max(b.position_mm[0][1]);
        }
        if [x, y]
            .into_iter()
            .any(|v| v.lo > radial_upper || v.hi < -radial_upper)
        {
            return Ok(Outside);
        }
        let rho = (x.square() + y.square()).sqrt();
        if !rho.valid() {
            return Ok(Unresolved);
        }
        if rho.lo > radial_upper {
            return Ok(Outside);
        }
        if z.lo <= start || z.hi >= end {
            return Ok(Unresolved);
        }
        if self.prove_radial(true, rho, z, w)? {
            Ok(Inside)
        } else if self.prove_radial(false, rho, z, w)? {
            Ok(Outside)
        } else {
            Ok(Unresolved)
        }
    }
    fn prove_radial(
        &self,
        inside: bool,
        rho: I,
        z: I,
        w: &mut QueryWork<'_>,
    ) -> Result<bool, RepairedHeadError> {
        for segment in 0..self.count {
            // Depth-first bisection needs at most65 simultaneously live cells.
            let mut pending = [([0., 0.], 0_u8); 65];
            pending[0] = ([0., 1.], 0);
            let mut count = 1;
            while count > 0 {
                count -= 1;
                let (range, depth) = pending[count];
                w.cell()?;
                let b = self.certified_interval_bounds(segment, range);
                let [r, axial] = b.position_mm.map(interval);
                if axial.hi < z.lo || axial.lo > z.hi {
                    continue;
                }
                let radially_inside = rho.hi < r.lo;
                let radially_outside = rho.lo > r.hi;
                if if inside {
                    radially_inside
                } else {
                    radially_outside
                } {
                    continue;
                }
                if if inside {
                    radially_outside
                } else {
                    radially_inside
                } {
                    // A broad interval cell may only APPEAR to meet the query
                    // height. Opposite signs from that false candidate cannot
                    // disprove a classification until actual endpoint order
                    // proves its continuous axial image intersects the box.
                    w.charge(1024)?;
                    let a = self
                        .certified_interval_bounds(segment, [range[0]; 2])
                        .position_mm[1];
                    let b = self
                        .certified_interval_bounds(segment, [range[1]; 2])
                        .position_mm[1];
                    if a[1] <= z.hi && b[0] >= z.lo {
                        return Ok(false);
                    }
                }
                let mid = range[0] + (range[1] - range[0]) / 2.;
                if depth == 64 || mid <= range[0] || mid >= range[1] {
                    return Ok(false);
                }
                pending[count] = ([mid, range[1]], depth + 1);
                pending[count + 1] = ([range[0], mid], depth + 1);
                count += 2;
            }
        }
        // Every source cell has been excluded axially or proved radially.
        // The constructor's connected ordered chain covers the entire strict
        // finite axial interval, so every actual crossing has been checked.
        Ok(true)
    }
}
#[cfg(test)]
mod tests;
