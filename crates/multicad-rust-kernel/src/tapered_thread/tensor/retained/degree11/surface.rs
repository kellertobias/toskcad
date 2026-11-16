//! Source surface integrals; no occupied-solid or manufacturing admission.
use super::normal::{SQUARED11_WORK, SquaredNormal11};
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
/// Fixed ledger indices for actual successful operation reservations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(usize)]
pub enum SurfaceStage11 {
    Setup,
    NormalConstruction,
    Fallback,
    Taylor,
    Gauss2,
    Gauss4,
    Subdivision,
    Reduction,
    Publication,
    Storage,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceFailure11 {
    pub error: Error,
    /// Invalid options have false, an empty ledger and no attempted stage.
    pub execution_started: bool,
    pub attempted_stage: Option<SurfaceStage11>,
    pub attempted_work: usize,
    pub attempted_bytes: usize,
    /// Successful debits only. The failed attempted charge is separate.
    pub work_by_stage: [usize; 10],
    pub work_used: usize,
    pub remaining_work: usize,
    pub accounted_peak_bytes: usize,
    pub evaluated_cells: usize,
    pub axis_samples: [usize; 2],
    pub axis_unresolved: [usize; 2],
    /// Max observed Gauss4 remainder radii: area, centered moments XYZ.
    /// An axis with no finite samples has no observed remainder bound.
    pub axis_remainder_max: [[f64; 4]; 2],
    /// Last COMPLETE source reduction: area, raw moments XYZ, centroid XYZ.
    /// None means no complete source reduction was reached. Infinite centroid
    /// widths mean the aggregate area was not proved strictly positive.
    pub last_aggregate_widths: Option<[f64; 7]>,
}
impl SurfaceFailure11 {
    fn invalid() -> Self {
        Self {
            error: Error::Invalid,
            execution_started: false,
            attempted_stage: None,
            attempted_work: 0,
            attempted_bytes: 0,
            work_by_stage: [0; 10],
            work_used: 0,
            remaining_work: 0,
            accounted_peak_bytes: 0,
            evaluated_cells: 0,
            axis_samples: [0; 2],
            axis_unresolved: [0; 2],
            axis_remainder_max: [[0.; 4]; 2],
            last_aggregate_widths: None,
        }
    }
}
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
    pub work_by_stage: [usize; 10],
    pub axis_samples: [usize; 2],
    pub axis_unresolved: [usize; 2],
    pub axis_remainder_max: [[f64; 4]; 2],
    pub last_aggregate_widths: Option<[f64; 7]>,
}
struct Budget<'a> {
    initial: usize,
    remaining: usize,
    peak: usize,
    evaluated: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
    work_by_stage: [usize; 10],
    attempted_stage: Option<SurfaceStage11>,
    attempted_work: usize,
    attempted_bytes: usize,
    axis_samples: [usize; 2],
    axis_unresolved: [usize; 2],
    axis_remainder_max: [[f64; 4]; 2],
    last_aggregate_widths: Option<[f64; 7]>,
}
impl Budget<'_> {
    fn charge(&mut self, stage: SurfaceStage11, n: usize) -> Result<(), Error> {
        self.attempted_stage = Some(stage);
        self.attempted_work = n;
        self.attempted_bytes = 0;
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        self.work_by_stage[stage as usize] += n;
        Ok(())
    }
    fn bytes(&mut self, n: usize, limit: usize) -> Result<(), Error> {
        self.attempted_stage = Some(SurfaceStage11::Storage);
        self.attempted_work = 0;
        self.attempted_bytes = n;
        if n > limit {
            return Err(Error::Limit);
        }
        self.peak = self.peak.max(n);
        Ok(())
    }
    fn failure(&self, error: Error) -> SurfaceFailure11 {
        SurfaceFailure11 {
            error,
            execution_started: true,
            attempted_stage: self.attempted_stage,
            attempted_work: self.attempted_work,
            attempted_bytes: self.attempted_bytes,
            work_by_stage: self.work_by_stage,
            work_used: self.initial - self.remaining,
            remaining_work: self.remaining,
            accounted_peak_bytes: self.peak,
            evaluated_cells: self.evaluated,
            axis_samples: self.axis_samples,
            axis_unresolved: self.axis_unresolved,
            axis_remainder_max: self.axis_remainder_max,
            last_aggregate_widths: self.last_aggregate_widths,
        }
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
        self.surface_properties_with_diagnostics(options, interrupt)
            .map_err(|f| f.error)
    }
    /// Same calculation and limits, retaining the real resource ledger on
    /// refusal. No geometry, parser owner or mutable budget escapes.
    #[allow(
        clippy::result_large_err,
        reason = "Allocation-free bounded diagnostic reporting after resource failure"
    )]
    pub fn surface_properties_with_diagnostics(
        &self,
        options: SurfaceOptions11,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<SurfaceProperties11, SurfaceFailure11> {
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
            return Err(SurfaceFailure11::invalid());
        }
        let mut budget = Budget {
            initial: options.max_work,
            remaining: options.max_work,
            peak: 0,
            evaluated: 0,
            interrupt,
            work_by_stage: [0; 10],
            attempted_stage: None,
            attempted_work: 0,
            attempted_bytes: 0,
            axis_samples: [0; 2],
            axis_unresolved: [0; 2],
            axis_remainder_max: [[0.; 4]; 2],
            last_aggregate_widths: None,
        };
        let result = (|| {
            budget.charge(SurfaceStage11::Setup, 1)?;
            integrate::calculate(self, options, &mut budget)
        })();
        result.map_err(|error| budget.failure(error))
    }
}
