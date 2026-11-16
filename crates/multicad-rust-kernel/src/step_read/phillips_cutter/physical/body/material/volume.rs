//! Volume of the ACTUAL distinct closed material, not the unsewn raw STEP.
//! New head integral + literal new tensor flux - actual Cutter section integral.
//! The common affine map contributes its independently certified determinant.
use super::*;
use crate::head_profile::{RepairedHeadVolume, RepairedHeadVolumeOptions};
mod tensor;
#[cfg(test)]
mod tests;
const WORK_LIMIT: usize = 512_000_000;
const BASE_SCRATCH: usize = 16 * 1024;
const STREAM_SCRATCH: usize = 128 * 1024;
const META_WORK: usize = 16384;
const CELL_WORK: usize = 32768;
#[derive(Clone, Copy)]
pub struct MaterialVolumeOptions {
    pub maximum_width_mm3: f64,
    pub maximum_head_cells: usize,
    pub maximum_cavity_evaluations: usize,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for MaterialVolumeOptions {
    fn default() -> Self {
        Self {
            maximum_width_mm3: 1e-6,
            maximum_head_cells: 65536,
            maximum_cavity_evaluations: 1_000_000,
            maximum_work: WORK_LIMIT,
            maximum_bytes: BYTES,
        }
    }
}
pub struct MaterialVolume<'v, 's, 'input> {
    source: &'v DerivedClosedMaterial<'s, 'input>,
    local: [f64; 2],
    world: [f64; 2],
    components: [[f64; 2]; 3],
    thread_cells: usize,
    head_cells: usize,
    cavity_evaluations: usize,
    thread_work: usize,
    head_work: usize,
    work: usize,
}
impl<'v, 's, 'input> MaterialVolume<'v, 's, 'input> {
    pub fn source(&self) -> &'v DerivedClosedMaterial<'s, 'input> {
        self.source
    }
    pub fn local_volume_mm3(&self) -> [f64; 2] {
        self.local
    }
    pub fn world_volume_mm3(&self) -> [f64; 2] {
        self.world
    }
    /// Head, thread and cavity, all in the common local source chart.
    pub fn component_volumes_mm3(&self) -> [[f64; 2]; 3] {
        self.components
    }
    pub fn thread_cells(&self) -> usize {
        self.thread_cells
    }
    pub fn head_cells(&self) -> usize {
        self.head_cells
    }
    pub fn cavity_evaluations(&self) -> usize {
        self.cavity_evaluations
    }
    pub fn thread_work(&self) -> usize {
        self.thread_work
    }
    pub fn head_work(&self) -> usize {
        self.head_work
    }
    /// Known material prerequisites +all fresh numerical query charges. Existing
    /// occupied-provider admission was a separate operation, not rerun/renewed.
    pub fn work_used(&self) -> usize {
        self.work
    }
}
fn interval([lo, hi]: [f64; 2]) -> I {
    I::new(lo, hi)
}
fn certified(v: I, width: f64) -> std::result::Result<[f64; 2], PhillipsReadError> {
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
/// Largest requested cap fitting actual logical workspace. Caller validates
/// maximum>=1 and available>=128KiB; both supplied scratch functions at1 fit
/// that minimum, and are monotone for their bounded integer input domains.
fn fitting_limit(maximum: usize, available: usize, bytes: fn(usize) -> usize) -> usize {
    let (mut lo, mut hi) = (1, maximum);
    while lo < hi {
        let mid = lo + (hi - lo).div_ceil(2);
        if bytes(mid) <= available {
            lo = mid
        } else {
            hi = mid - 1
        }
    }
    lo
}
impl<'s, 'input> DerivedClosedMaterial<'s, 'input> {
    fn thread_flux(
        &self,
        work: &mut Work<'_>,
    ) -> std::result::Result<(I, usize), PhillipsReadError> {
        let input = self.thread.input();
        let weights = tensor::weights();
        let mut total = I::point(0.);
        let mut cells = 0;
        // This private material source can only contain the occupied candidate
        // built by derived_thread::graph:2 EXACT horizontal cap planes and all
        // remaining faces degree9 tensors with identical literal shared rails.
        // F=(x/2,y/2,0) has identically zero cap flux. Source coefficients below
        // are the explicitly chosen NEW binary64 geometry, not cached original
        // decimal leaves. All physical tensor cells are traversed exactly once.
        for face in &input.faces {
            if let crate::radial_tensor::Surface::Tensor(t) = &input.surfaces[face.surface] {
                let [nu, nv] = t.spans();
                let h = 9 * nv + 1;
                for u in 0..nu {
                    for v in 0..nv {
                        // <=8100 fixed3-vector normal terms +10728 contracted
                        // component terms +100 pole copies/remaining fixed setup.
                        // No sampled normal, adaptive integration or mesh occurs.
                        work.step(CELL_WORK).map_err(|e| read_error(e, work))?;
                        let p = std::array::from_fn(|i| {
                            t.poles()[(9 * u + i / 10) * h + 9 * v + i % 10]
                        });
                        let value = tensor::cell(p, &weights);
                        work.step(0).map_err(|e| read_error(e, work))?;
                        total = total + if face.forward { value } else { -value };
                        cells += 1;
                    }
                }
            }
        }
        Ok((total, cells))
    }
    pub fn volume_bounds<'v>(
        &'v self,
        options: MaterialVolumeOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<MaterialVolume<'v, 's, 'input>, PhillipsReadError> {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        work.step(1).map_err(|e| read_error(e, &work))?;
        if !options.maximum_width_mm3.is_finite()
            || options.maximum_width_mm3 <= 0.
            || options.maximum_head_cells == 0
            || options.maximum_head_cells > 65536
            || options.maximum_cavity_evaluations == 0
            || options.maximum_cavity_evaluations > 1_000_000
            || options.maximum_work > WORK_LIMIT
            || options.maximum_bytes > BYTES
        {
            return Err(PhillipsReadError::Invalid);
        }
        work.step(self.work_used())
            .and_then(|_| work.step(META_WORK))
            .map_err(|e| read_error(e, &work))?;
        let prefix = self
            .retained_budget_bytes()
            .saturating_add(std::mem::size_of::<MaterialVolume<'_, '_, '_>>())
            .saturating_add(BASE_SCRATCH);
        let available = options
            .maximum_bytes
            .checked_sub(prefix)
            .filter(|n| *n >= STREAM_SCRATCH)
            .ok_or(PhillipsReadError::Limit)?;
        // ONE cell at a time: calculation grid~5KiB, normals/derivatives<32KiB,
        // two contracted fields<20KiB, weights/stack/setup<16KiB. The128KiB
        // reserve also covers transient old/new Vec buffers. Nothing is cached.
        let before = work.remaining;
        let (thread, thread_cells) = self.thread_flux(&mut work)?;
        let thread_work = before - work.remaining;
        let determinant = interval(self.determinant_bounds());
        let preliminary =
            (I::point(options.maximum_width_mm3) / (I::point(8.) * I::point(determinant.hi))).lo;
        certified(thread, preliminary)?;
        // A direct whole-source cylinder bound, solely for the affine numerical
        // reserve. This does NOT use old Profile volume or infer raw STEP mass.
        let mut radius = 0f64;
        for i in 0..self.head.segment_count() {
            radius = radius.max(self.head.certified_interval_bounds(i, [0., 1.]).position_mm[0][1]);
        }
        let top = self.head.station(0).unwrap()[1];
        let end = self.head.station(self.head.segment_count()).unwrap()[1];
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let outer =
            pi * I::point(radius).square() * (I::point(end) - I::point(top)) + I::point(thread.hi);
        let reserve = ((I::point(determinant.hi) - I::point(determinant.lo)) * outer).hi;
        let remaining_width = I::point(options.maximum_width_mm3) - I::point(reserve);
        if !remaining_width.lo.is_finite() || remaining_width.lo <= 0. {
            return Err(PhillipsReadError::Unsupported);
        }
        let component_width =
            (I::point(remaining_width.lo) / (I::point(8.) * I::point(determinant.hi))).lo;
        certified(thread, component_width)?;
        // Operations are sequential. The complete shared live source is in
        // prefix; only the largest active component workspace is added to it.
        let head_limit = fitting_limit(
            options.maximum_head_cells,
            available,
            RepairedHeadVolume::scratch_budget_bytes,
        );
        let (head, head_cells, head_work) = {
            let result = self
                .head
                .volume_bounds(
                    RepairedHeadVolumeOptions {
                        maximum_width_mm3: component_width,
                        maximum_cells: head_limit,
                        maximum_work: work.remaining.min(16_000_000),
                        maximum_bytes: options.maximum_bytes,
                    },
                    interrupt,
                )
                .map_err(super::super::closed_head::repaired_error)?;
            (
                interval(result.volume_mm3()),
                result.cells(),
                result.work_used(),
            )
        };
        work.step(head_work).map_err(|e| read_error(e, &work))?;
        let cavity_max = options.maximum_cavity_evaluations.min(work.remaining);
        if cavity_max == 0 {
            return Err(PhillipsReadError::Limit);
        }
        let cavity_limit =
            fitting_limit(cavity_max, available, Cutter::volume_scratch_budget_bytes);
        let result = self
            .cavity
            .volume_bounds_with_cancel(component_width, cavity_limit, &|| {
                interrupt.is_some_and(|f| f())
            })
            .map_err(cavity_error)?;
        let cavity = interval(result.volume_mm3);
        let cavity_evaluations = result.section_evaluations;
        work.step(cavity_evaluations)
            .map_err(|e| read_error(e, &work))?;
        let mut local = head + thread - cavity;
        // Actual containment and strict floor prove the ENTIRE occupied thread
        // belongs to this material. Intersect with that earned positive volume
        // lower bound; no tolerance zeroing or change of geometry is performed.
        local.lo = local.lo.max(thread.lo);
        // The source cylinder also bounds the true material from above. Keep
        // the determinant reserve tied to this earned upper bound even when a
        // deliberately loose component query returned a wider enclosure.
        local.hi = local.hi.min(outer.hi);
        let world = local * determinant;
        let result = certified(world, options.maximum_width_mm3);
        work.step(1).map_err(|e| read_error(e, &work))?;
        Ok(MaterialVolume {
            source: self,
            local: [local.lo, local.hi],
            world: result?,
            components: [
                [head.lo, head.hi],
                [thread.lo, thread.hi],
                [cavity.lo, cavity.hi],
            ],
            thread_cells,
            head_cells,
            cavity_evaluations,
            thread_work,
            head_work,
            work: options.maximum_work - work.remaining,
        })
    }
}
