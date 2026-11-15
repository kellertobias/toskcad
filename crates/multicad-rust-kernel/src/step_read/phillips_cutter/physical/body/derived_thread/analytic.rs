//! A DISTINCT analytic thread repair, qualified against every decoded tensor
//! coefficient cell. Header data is an untrusted proposal only. This receipt
//! does not admit a complete head/cavity material, raw-file mass or STEP export.
use super::*;
mod common_chart;
mod curve;
mod finite;
pub use common_chart::CommonChartPeriodicThreadRepair;
pub use finite::{
    PeriodicFiniteCurve, PeriodicFiniteCurves, PeriodicFiniteKind, PeriodicFiniteOptions,
    PeriodicShoulderOptions, PeriodicShoulderRepair,
};
mod periodic;
mod rotation;
pub use periodic::{PeriodicThreadBounds, PeriodicThreadRepair, PeriodicThreadSource};
const TOTAL: usize = 32_000_000;
const SCRATCH_ANALYTIC: usize = 2 * 1024 * 1024;
#[derive(Clone, Copy)]
pub struct AnalyticThreadOptions {
    pub maximum_distance_mm: f64,
    /// Derivative discrepancy in each ORIGINAL normalized Bernstein cell.
    pub maximum_cell_derivative_mm: f64,
    /// New work only, after the separately reported closed-head prerequisite.
    pub maximum_work: usize,
    pub maximum_total_work: usize,
    pub maximum_bytes: usize,
}
impl Default for AnalyticThreadOptions {
    fn default() -> Self {
        Self {
            maximum_distance_mm: 1e-7,
            maximum_cell_derivative_mm: 1e-7,
            maximum_work: WORK,
            maximum_total_work: TOTAL,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct AnalyticThreadTransfer {
    pub face_entity: u32,
    pub surface_entity: u32,
    pub original_rail_entities: [u32; 2],
    pub original_v_start_cell: usize,
    pub original_cells: [usize; 2],
    pub distance_mm: f64,
    pub derivative_mm: [f64; 2],
}
pub struct AnalyticThreadRepair<'a, 'input> {
    body: &'a Body<'input>,
    thread: crate::helical_thread::Thread,
    periodic: Option<PeriodicThreadSource>,
    chart: Option<CommonMouthChart>,
    frame: [Point; 3],
    transfers: Vec<AnalyticThreadTransfer>,
    new_work: usize,
    inherited_work: usize,
    retained: usize,
}
impl<'a, 'input> AnalyticThreadRepair<'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.body
    }
    pub fn thread(&self) -> crate::helical_thread::Thread {
        self.thread
    }
    pub fn frame(&self) -> [Point; 3] {
        self.frame
    }
    pub fn transfers(&self) -> &[AnalyticThreadTransfer] {
        &self.transfers
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn inherited_head_work_used(&self) -> usize {
        self.inherited_work
    }
    pub fn total_work_used(&self) -> usize {
        self.new_work + self.inherited_work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
fn valid(o: AnalyticThreadOptions) -> bool {
    o.maximum_distance_mm.is_finite()
        && o.maximum_distance_mm > 0.
        && o.maximum_distance_mm <= 1e-7
        && o.maximum_cell_derivative_mm.is_finite()
        && o.maximum_cell_derivative_mm > 0.
        && o.maximum_cell_derivative_mm <= 1e-7
        && o.maximum_work <= WORK
        && o.maximum_total_work <= TOTAL
        && o.maximum_bytes <= BYTES
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn analytic_thread_repair(
        &self,
        options: AnalyticThreadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<AnalyticThreadRepair<'a, 'input>, PhillipsReadError> {
        self.build_thread_repair(options, interrupt, false, None)
    }
    fn build_thread_repair(
        &self,
        options: AnalyticThreadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
        periodic: bool,
        chart: Option<CommonMouthChart>,
    ) -> std::result::Result<AnalyticThreadRepair<'a, 'input>, PhillipsReadError> {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !valid(options) {
                return Err(ReadError::Invalid);
            }
            let allowance = options.maximum_work.min(
                options
                    .maximum_total_work
                    .checked_sub(self.work_used())
                    .ok_or(ReadError::Limit)?,
            );
            work.remaining = allowance.checked_sub(1).ok_or(ReadError::Limit)?;
            let body = self.body();
            let thread = body.thread;
            // Rounded roots require their own non-affine radial Taylor law;
            // no flat-root surrogate or template substitution is admitted.
            if thread.spec().root_radius_mm.is_some() || thread.spec().starts != 1 {
                return Err(ReadError::Unsupported);
            }
            let base =
                self.retained_budget_bytes() + std::mem::size_of::<AnalyticThreadRepair<'_, '_>>();
            if base.saturating_add(SCRATCH_ANALYTIC) > options.maximum_bytes {
                return Err(ReadError::Limit);
            }
            let mut shell = Shell::new(&body.entities, body.tolerance, true);
            let raw = raw(
                self,
                &mut shell,
                &mut work,
                options.maximum_bytes - base - SCRATCH_ANALYTIC,
            )?;
            let count: usize = raw.iter().map(|r| r.ids.len() * r.ids[0].len()).sum();
            let factors: usize =
                raw.iter().map(|r| r.ids.len()).sum::<usize>() + raw[0].ids[0].len();
            if count > MAX_POINTS
                || base
                    .saturating_add(SCRATCH_ANALYTIC)
                    .saturating_add(count.saturating_mul(160))
                    .saturating_add(factors.saturating_mul(64))
                    > options.maximum_bytes
            {
                return Err(ReadError::Limit);
            }
            // Fixed untrusted factor-generation work; every actual field still
            // receives the complete measured transfer below. Counts are earned
            // by the SAME original grids, so no guessed cheap source wins quota.
            work.step(factors.saturating_mul(1024))?;
            let reference = crate::step::helical::prepare_retained_reference(
                &thread,
                body.visible[0],
                body.tolerance,
            )
            .map_err(helical::step_error)?;
            // Body::outer_carrier_bounds already prepared this SAME immutable
            // Body thread/visible/tolerance triple (outer.rs116), checked the
            // band count against tensor_face_count and validated every grid via
            // helical::tensor_surface. One-start non-rounded profiles have four
            // phase bands. A discrepancy here is an internal implementation bug.
            assert_eq!(raw.len(), reference.bands.len());
            assert_eq!(reference.bands.len(), 4);
            let rotation = rotation::Rotation::new(&mut work)?;
            let source = if periodic {
                Some(PeriodicThreadSource::new(thread, &mut work)?)
            } else {
                None
            };
            let geometry =
                curve::factors(&thread, &reference, &rotation, source.as_ref(), &mut work)?;
            let frame = body.frame;
            let target = if let Some(chart) = chart {
                Target {
                    columns: chart.columns,
                    z_shift: chart.outer_top_z,
                }
            } else {
                Target::raw(frame)
            };
            let mut transfers = Vec::with_capacity(raw.len());
            for (index, r) in raw.iter().enumerate() {
                let distance =
                    curve::literal_transfer(r, &reference, index, &target, &shell, &mut work)?;
                let bounds = curve::world_transfer(index, &target, distance, &geometry, &mut work)?;
                if bounds[0] > options.maximum_distance_mm
                    || bounds[1] > options.maximum_cell_derivative_mm
                    || bounds[2] > options.maximum_cell_derivative_mm
                {
                    return Err(ReadError::Unsupported);
                }
                transfers.push(AnalyticThreadTransfer {
                    face_entity: r.face,
                    surface_entity: r.surface,
                    original_rail_entities: r.rails,
                    original_v_start_cell: r.entry / 11,
                    original_cells: [(r.ids.len() - 1) / 11, (r.ids[0].len() - 1) / 11],
                    distance_mm: bounds[0],
                    derivative_mm: [bounds[1], bounds[2]],
                });
            }
            work.step(1)?;
            let retained =
                base + transfers.capacity() * std::mem::size_of::<AnalyticThreadTransfer>();
            Ok(AnalyticThreadRepair {
                body,
                thread,
                periodic: source,
                chart,
                frame,
                transfers,
                new_work: allowance - work.remaining,
                inherited_work: self.work_used(),
                retained,
            })
        })();
        result.map_err(|e| bound_error(e, &work))
    }
}
#[cfg(test)]
mod tests;

#[derive(Clone, Copy)]
struct Target {
    columns: [Point; 4],
    z_shift: f64,
}
impl Target {
    fn raw([o, u, v]: [Point; 3]) -> Self {
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        Self {
            columns: [o, u, v, n],
            z_shift: 0.,
        }
    }
}
