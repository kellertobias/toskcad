//! Actual DISTINCT repair-source mass: repaired meridian integral + exact
//! real periodic radial section integral - complete Cutter section integral.
//! No old Profile volume, original-file mass or mesh enters this result.
use super::*;
use crate::{
    head_profile::{RepairedHeadVolume, RepairedHeadVolumeOptions},
    phillips_drive::Cutter,
};
#[cfg(test)]
mod tests;
const BASE: usize = 16 * 1024;
const STREAM: usize = 128 * 1024;
const META: usize = 16384;
#[derive(Clone, Copy)]
pub struct PeriodicVolumeOptions {
    pub maximum_width_mm3: f64,
    pub maximum_head_cells: usize,
    pub maximum_cavity_evaluations: usize,
    pub maximum_work: usize,
    pub maximum_total_work: usize,
    pub maximum_bytes: usize,
}
impl Default for PeriodicVolumeOptions {
    fn default() -> Self {
        Self {
            maximum_width_mm3: 1e-6,
            maximum_head_cells: 65536,
            maximum_cavity_evaluations: 1_000_000,
            maximum_work: WORK,
            maximum_total_work: TOTAL,
            maximum_bytes: BYTES,
        }
    }
}
pub struct PeriodicMaterialVolume<'v, 's, 'input> {
    source: &'v PeriodicClosedMaterial<'s, 'input>,
    local: [f64; 2],
    world: [f64; 2],
    components: [[f64; 2]; 3],
    head_cells: usize,
    cavity_evaluations: usize,
    head_work: usize,
    thread_work: usize,
    new_work: usize,
    retained: usize,
}
impl<'v, 's, 'input> PeriodicMaterialVolume<'v, 's, 'input> {
    pub fn source(&self) -> &'v PeriodicClosedMaterial<'s, 'input> {
        self.source
    }
    pub fn local_volume_mm3(&self) -> [f64; 2] {
        self.local
    }
    pub fn world_volume_mm3(&self) -> [f64; 2] {
        self.world
    }
    /// Actual source head, thread and full cavity in the common local chart.
    pub fn component_volumes_mm3(&self) -> [[f64; 2]; 3] {
        self.components
    }
    pub fn head_cells(&self) -> usize {
        self.head_cells
    }
    pub fn cavity_evaluations(&self) -> usize {
        self.cavity_evaluations
    }
    pub fn head_work_used(&self) -> usize {
        self.head_work
    }
    pub fn thread_work_used(&self) -> usize {
        self.thread_work
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn inherited_work_used(&self) -> usize {
        self.source.total_work_used()
    }
    pub fn total_work_used(&self) -> usize {
        self.new_work + self.inherited_work_used()
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
fn interval([a, b]: [f64; 2]) -> I {
    I::new(a, b)
}
fn certified(v: I, width: f64) -> FResult<[f64; 2]> {
    if v.valid()
        && v.lo.is_finite()
        && v.hi.is_finite()
        && v.lo > 0.
        && (I::point(v.hi) - I::point(v.lo)).hi <= width
    {
        Ok([v.lo, v.hi])
    } else {
        Err(PhillipsReadError::Unsupported)
    }
}
/// Caller has already validated maximum>=1 and available>=STREAM; the two
/// bounded monotone component scratch functions at1 fit that reserve.
fn fitting_limit(maximum: usize, available: usize, bytes: fn(usize) -> usize) -> usize {
    let (mut lo, mut hi) = (1, maximum);
    while lo < hi {
        let mid = lo + (hi - lo).div_ceil(2);
        if bytes(mid) <= available {
            lo = mid;
        } else {
            hi = mid - 1;
        }
    }
    lo
}
/// The source owns ordered REAL phase nodes, positive radii and a full wrap by
/// identity. Each band's r is affine in phase, so integrating r² is exact with
/// the symmetric quadratic Bernstein average. The trim interval is supplied
/// only by the validated same-source finite/shoulder receipt, never a header.
fn periodic_volume(
    source: &PeriodicThreadSource,
    stations: [f64; 2],
    work: &mut Work<'_>,
) -> FResult<I> {
    work.step(1024)?;
    let rotation = rotation::Rotation::new(work)?;
    let pi = rotation.pi(work)?;
    let p = source.phases();
    let radii = source.radius_nodes_mm();
    let mut radial = I::point(0.);
    for band in 0..4 {
        work.step(256)?;
        let a = I::point(radii[band]);
        let b = I::point(radii[band + 1]);
        radial =
            radial + (p[band + 1] - p[band]) * (a.square() + a * b + b.square()) / I::point(3.);
    }
    Ok(pi * radial / I::point(source.pitch_mm()) * (I::point(stations[1]) - I::point(stations[0])))
}
impl<'s, 'input> PeriodicClosedMaterial<'s, 'input> {
    pub fn volume_bounds<'v>(
        &'v self,
        o: PeriodicVolumeOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<PeriodicMaterialVolume<'v, 's, 'input>> {
        let mut work = Work {
            remaining: o.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !o.maximum_width_mm3.is_finite()
                || o.maximum_width_mm3 <= 0.
                || o.maximum_head_cells == 0
                || o.maximum_head_cells > 65536
                || o.maximum_cavity_evaluations == 0
                || o.maximum_cavity_evaluations > 1_000_000
                || o.maximum_work > WORK
                || o.maximum_total_work > TOTAL
                || o.maximum_bytes > BYTES
            {
                return Err(PhillipsReadError::Invalid);
            }
            let allowance = o.maximum_work.min(
                o.maximum_total_work
                    .checked_sub(self.total_work_used())
                    .ok_or(PhillipsReadError::Limit)?,
            );
            work.remaining = allowance.checked_sub(1).ok_or(PhillipsReadError::Limit)?;
            work.step(META)?;
            let retained = self
                .retained_budget_bytes()
                .saturating_add(std::mem::size_of::<PeriodicMaterialVolume<'_, '_, '_>>());
            let prefix = retained.saturating_add(BASE);
            let available = o
                .maximum_bytes
                .checked_sub(prefix)
                .filter(|n| *n >= STREAM)
                .ok_or(PhillipsReadError::Limit)?;
            let before = work.remaining;
            let thread = periodic_volume(
                self.join.periodic_source(),
                self.join.finite_curves().cap_stations_mm(),
                &mut work,
            )?;
            let thread_work = before - work.remaining;
            let determinant = interval(self.determinant_bounds());
            let preliminary =
                (I::point(o.maximum_width_mm3) / (I::point(8.) * I::point(determinant.hi))).lo;
            certified(thread, preliminary)?;
            let head = self.join.head_source();
            let mut radius = 0_f64;
            for i in 0..head.segment_count() {
                radius = radius.max(head.certified_interval_bounds(i, [0., 1.]).position_mm[0][1]);
            }
            let top = head.station(0).unwrap()[1];
            let end = head.station(head.segment_count()).unwrap()[1];
            // This finite cylinder contains the entire ACTUAL repaired head, not its
            // legacy Profile; it is used only to allocate the affine numerical reserve.
            let pi = I::new(
                std::f64::consts::PI.next_down(),
                std::f64::consts::PI.next_up(),
            );
            let outer = pi * I::point(radius).square() * (I::point(end) - I::point(top))
                + I::point(thread.hi);
            let reserve = ((I::point(determinant.hi) - I::point(determinant.lo)) * outer).hi;
            let remaining_width = I::point(o.maximum_width_mm3) - I::point(reserve);
            if !remaining_width.lo.is_finite() || remaining_width.lo <= 0. {
                return Err(PhillipsReadError::Unsupported);
            }
            let component_width =
                (I::point(remaining_width.lo) / (I::point(8.) * I::point(determinant.hi))).lo;
            certified(thread, component_width)?;
            let head_limit = fitting_limit(
                o.maximum_head_cells,
                available,
                RepairedHeadVolume::scratch_budget_bytes,
            );
            let h = head
                .volume_bounds(
                    RepairedHeadVolumeOptions {
                        maximum_width_mm3: component_width,
                        maximum_cells: head_limit,
                        maximum_work: work.remaining.saturating_add(head.work_used()).min(WORK),
                        maximum_bytes: o.maximum_bytes,
                    },
                    interrupt,
                )
                .map_err(
                    crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error,
                )?;
            let head_volume = interval(h.volume_mm3());
            let head_cells = h.cells();
            let head_work = h
                .work_used()
                .checked_sub(head.work_used())
                .expect("head integral receipt includes its immutable source prerequisite");
            work.step(head_work)?;
            let cavity_max = o.maximum_cavity_evaluations.min(work.remaining);
            if cavity_max == 0 {
                return Err(PhillipsReadError::Limit);
            }
            let cavity_limit =
                fitting_limit(cavity_max, available, Cutter::volume_scratch_budget_bytes);
            let c = self
                .containment
                .source()
                .volume_bounds_with_cancel(component_width, cavity_limit, &|| {
                    interrupt.is_some_and(|f| f())
                })
                .map_err(cavity_error)?;
            let cavity = interval(c.volume_mm3);
            let cavity_evaluations = c.section_evaluations;
            work.step(cavity_evaluations)?;
            let mut local = head_volume + thread - cavity;
            // Complete repaired-head containment and its strict floor leave the entire
            // thread in material, while the source cylinder supplies an upper bound.
            local.lo = local.lo.max(thread.lo);
            local.hi = local.hi.min(outer.hi);
            let world = certified(local * determinant, o.maximum_width_mm3)?;
            work.step(1)?;
            Ok(PeriodicMaterialVolume {
                source: self,
                local: [local.lo, local.hi],
                world,
                components: [
                    [head_volume.lo, head_volume.hi],
                    [thread.lo, thread.hi],
                    [cavity.lo, cavity.hi],
                ],
                head_cells,
                cavity_evaluations,
                head_work,
                thread_work,
                new_work: allowance - work.remaining,
                retained,
            })
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result
        }
    }
}
