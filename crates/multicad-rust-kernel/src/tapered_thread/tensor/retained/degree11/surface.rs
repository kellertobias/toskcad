//! Source surface integrals; no occupied-solid or manufacturing admission.
use super::*;
mod integrate;
mod quadrature;
#[cfg(test)]
mod tests;
type Values = [I; 4];
const MAX_WORK: usize = 512_000_000;
const MAX_BYTES: usize = 64 * 1024 * 1024;
const MAX_CELLS: usize = 8192;
const MAX_DEPTH: u8 = 48;
/// Fresh calculation limits and absolute FULL-width output requirements.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceOptions11 {
    pub area_width_mm2: f64,
    pub first_moment_width_mm3: f64,
    pub centroid_width_mm: f64,
    pub max_work: usize,
    pub max_bytes: usize,
    pub max_cells: usize,
    pub max_depth: u8,
}
impl Default for SurfaceOptions11 {
    fn default() -> Self {
        Self {
            area_width_mm2: 1e-5,
            first_moment_width_mm3: 1e-5,
            centroid_width_mm: 1e-5,
            max_work: MAX_WORK,
            max_bytes: MAX_BYTES,
            max_cells: MAX_CELLS,
            max_depth: MAX_DEPTH,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceProperties11 {
    pub area_mm2: [f64; 2],
    pub first_moment_mm3: [[f64; 2]; 3],
    pub centroid_mm: [[f64; 2]; 3],
    pub centered_first_moment_mm3: [[f64; 2]; 3],
    pub center_mm: [f64; 3],
    pub work_used: usize,
    pub remaining_work: usize,
    /// Conservative live source/capacity/workspace reservation, not RSS.
    pub accounted_peak_bytes: usize,
    pub evaluated_cells: usize,
}
struct Budget<'a> {
    remaining: usize,
    peak: usize,
    evaluated: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Budget<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        Ok(())
    }
    fn bytes(&mut self, n: usize, limit: usize) -> Result<(), Error> {
        if n > limit {
            return Err(Error::Limit);
        }
        self.peak = self.peak.max(n);
        Ok(())
    }
}
fn pair(x: I) -> [f64; 2] {
    [x.lo, x.hi]
}
impl RetainedTensor11 {
    /// Integrate the unchanged finite tensor cells. Normals come from their
    /// actual degree-(21,21) polynomial, never a tessellation or source recipe.
    pub fn surface_properties(
        &self,
        options: SurfaceOptions11,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<SurfaceProperties11, Error> {
        let positive = [
            options.area_width_mm2,
            options.first_moment_width_mm3,
            options.centroid_width_mm,
        ];
        if positive.iter().any(|x| !x.is_finite() || *x <= 0.)
            || options.max_work == 0
            || options.max_work > MAX_WORK
            || options.max_bytes == 0
            || options.max_bytes > MAX_BYTES
            || options.max_cells == 0
            || options.max_cells > MAX_CELLS
            || options.max_depth > MAX_DEPTH
        {
            return Err(Error::Invalid);
        }
        let mut budget = Budget {
            remaining: options.max_work,
            peak: 0,
            evaluated: 0,
            interrupt,
        };
        budget.charge(1)?;
        integrate::calculate(self, options, &mut budget)
    }
}
