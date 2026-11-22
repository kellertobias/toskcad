//! Independent volume of the earned closed material: actual head integral plus
//! original D11 signed flux minus the actual contained closed Cutter volume.
use super::*;
use crate::head_profile::{RepairedHeadVolume, RepairedHeadVolumeOptions};
use crate::tapered_thread::{NORMAL11_WORK, NormalPolynomial11};
use std::sync::Arc;
mod supplied;
mod tensor;
pub use supplied::{SuppliedDegree11MaterialMatch, SuppliedDegree11MaterialOptions};
#[cfg(test)]
mod tests;
const WORK_LIMIT: usize = 512_000_000;
const BASE_SCRATCH: usize = 16 * 1024;
const STREAM_SCRATCH: usize = 128 * 1024;
const META_WORK: usize = 16384;
// Actual normal construction includes its derivative caches. The additional
//49152 credits cover37841 directed contraction operations, original cell/pole
// setup, Arc publication and all loop/cancellation bookkeeping.
const CELL_WORK: usize = NORMAL11_WORK + 49_152;
const NORMAL_ARC_BYTES: usize = 64;
fn normal_error(e: crate::tapered_thread::Error, work: &Work<'_>) -> ReadError {
    use crate::tapered_thread::Error as E;
    match e {
        E::Invalid => ReadError::Invalid,
        E::Limit => ReadError::Limit,
        E::Cancelled => {
            work.cancelled.set(true);
            ReadError::Limit
        }
        E::Unresolved | E::Unsupported => ReadError::Unsupported,
    }
}
fn cached_payload(normals: &[Arc<NormalPolynomial11>]) -> usize {
    normals
        .iter()
        .map(|n| n.retained_budget_bytes() + NORMAL_ARC_BYTES)
        .sum()
}
fn cached_bytes(normals: &[Arc<NormalPolynomial11>], capacity: usize) -> usize {
    capacity * std::mem::size_of::<Arc<NormalPolynomial11>>() + cached_payload(normals)
}
#[derive(Clone, Copy)]
pub struct Degree11MaterialVolumeOptions {
    pub maximum_width_mm3: f64,
    pub maximum_head_cells: usize,
    pub maximum_cavity_evaluations: usize,
    /// Whole material aggregate plus fresh numerical work, inherited once.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for Degree11MaterialVolumeOptions {
    fn default() -> Self {
        Self {
            maximum_width_mm3: 1e-6,
            maximum_head_cells: 65536,
            maximum_cavity_evaluations: 1_000_000,
            maximum_work: WORK_LIMIT,
            maximum_bytes: MAX_BYTES,
        }
    }
}
pub struct Degree11MaterialVolume<'v, 's, 'a, 'input> {
    source: &'v Degree11ClosedMaterial<'s, 'a, 'input>,
    local: [f64; 2],
    world: [f64; 2],
    components: [[f64; 2]; 3],
    thread_cells: usize,
    head_cells: usize,
    cavity_evaluations: usize,
    thread_work: usize,
    head_work: usize,
    new_work: usize,
    peak: usize,
    normals: Vec<Arc<NormalPolynomial11>>,
}
impl<'v, 's, 'a, 'input> Degree11MaterialVolume<'v, 's, 'a, 'input> {
    pub fn source(&self) -> &'v Degree11ClosedMaterial<'s, 'a, 'input> {
        self.source
    }
    pub fn local_volume_mm3(&self) -> [f64; 2] {
        self.local
    }
    pub fn world_volume_mm3(&self) -> [f64; 2] {
        self.world
    }
    /// Actual head, full original D11 thread and closed cavity local integrals.
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
    /// Fresh head integration only; its retained construction is inherited once.
    pub fn head_work(&self) -> usize {
        self.head_work
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn work_used(&self) -> usize {
        self.source.work_used() + self.new_work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.source.retained_budget_bytes()
            + std::mem::size_of::<Self>()
            + cached_bytes(&self.normals, self.normals.capacity())
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
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
impl<'s, 'a, 'input> Degree11ClosedMaterial<'s, 'a, 'input> {
    fn thread_flux(
        &self,
        work: &mut Work<'_>,
        count: usize,
    ) -> Result<(I, usize, Vec<Arc<NormalPolynomial11>>)> {
        // Bounded integer binomial loops plus396 directed weight divisions.
        work.step(32_768)?;
        let weights = tensor::Quadrature::new();
        let mut normals = Vec::new();
        normals
            .try_reserve_exact(count)
            .map_err(|_| ReadError::Limit)?;
        if normals.capacity() > count {
            return Err(ReadError::Limit);
        }

        let mut total = I::point(0.);
        let mut cells = 0;
        // Every retained face is degree11 and every actual cell is visited once.
        // Horizontal tip/join/top faces have identically zero Fxy flux.
        for face in self.thread.carrier().faces() {
            let source = face.source();
            let [nu, nv] = source.spans();
            let height = source.dimensions()[1];
            for u in 0..nu {
                for v in 0..nv {
                    work.step(CELL_WORK)?;
                    let p = std::array::from_fn(|i| {
                        source.poles()[(11 * u + i / 12) * height + 11 * v + i % 12].map(I::point)
                    });
                    let cell = source.cell([u, v]).map_err(|e| normal_error(e, work))?;
                    let normal = Arc::new(
                        NormalPolynomial11::new(&cell).map_err(|e| normal_error(e, work))?,
                    );
                    if normal.retained_budget_bytes()
                        > NormalPolynomial11::retained_capacity_bound_bytes()
                    {
                        return Err(ReadError::Limit);
                    }
                    let value = tensor::cached_cell(p, &weights, &normal);
                    normals.push(normal);
                    total = total + if face.forward() { value } else { -value };
                    work.step(0)?;
                    cells += 1;
                }
            }
        }
        Ok((total, cells, normals))
    }
    pub fn volume_bounds<'v>(
        &'v self,
        options: Degree11MaterialVolumeOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<Degree11MaterialVolume<'v, 's, 'a, 'input>, PhillipsReadError> {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !options.maximum_width_mm3.is_finite()
                || options.maximum_width_mm3 <= 0.
                || options.maximum_head_cells == 0
                || options.maximum_head_cells > 65536
                || options.maximum_cavity_evaluations == 0
                || options.maximum_cavity_evaluations > 1_000_000
                || options.maximum_work > WORK_LIMIT
                || options.maximum_bytes > MAX_BYTES
            {
                return Err(ReadError::Invalid);
            }
            work.step(self.work_used())?;
            work.step(META_WORK)?;
            let count = self
                .thread
                .carrier()
                .faces()
                .iter()
                .try_fold(0usize, |sum, f| {
                    let [u, v] = f.source().spans();
                    sum.checked_add(u.checked_mul(v).ok_or(ReadError::Limit)?)
                        .ok_or(ReadError::Limit)
                })?;
            let cache_bound = count
                .checked_mul(
                    NormalPolynomial11::retained_capacity_bound_bytes()
                        + NORMAL_ARC_BYTES
                        + std::mem::size_of::<Arc<NormalPolynomial11>>(),
                )
                .ok_or(ReadError::Limit)?;
            let prefix = self
                .retained_budget_bytes()
                .checked_add(std::mem::size_of::<Degree11MaterialVolume<'_, '_, '_, '_>>())
                .and_then(|n| n.checked_add(BASE_SCRATCH))
                .and_then(|n| n.checked_add(cache_bound))
                .ok_or(ReadError::Limit)?;
            let available = options
                .maximum_bytes
                .checked_sub(prefix)
                .filter(|n| *n >= STREAM_SCRATCH)
                .ok_or(ReadError::Limit)?;
            // Fixed stack arrays: source6912B, derivatives12672B, fullnormal
            //23232B, weights6336B and contraction4224B.128KiB covers all live
            // copies, return temporaries, helper frames and source-cell setup.
            let before = work.remaining;
            let (thread, thread_cells, normals) = self.thread_flux(&mut work, count)?;
            let thread_work = before - work.remaining;
            Ok((
                prefix,
                available,
                thread,
                thread_cells,
                thread_work,
                normals,
            ))
        })();
        let (prefix, available, thread, thread_cells, thread_work, normals) =
            result.map_err(|e| bound_error(e, &work))?;
        let determinant = interval(self.determinant_bounds());
        let preliminary =
            (I::point(options.maximum_width_mm3) / (I::point(8.) * I::point(determinant.hi))).lo;
        certified(thread, preliminary)?;
        let mut radius = 0f64;
        for i in 0..self.head.segment_count() {
            radius = radius.max(self.head.certified_interval_bounds(i, [0., 1.]).position_mm[0][1]);
        }
        // Material construction earned at least one original head segment.
        let top = self.head.station(0).expect("closed head top station")[1];
        let end = self
            .head
            .station(self.head.segment_count())
            .expect("closed head bottom station")[1];
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let outer =
            pi * I::point(radius).square() * (I::point(end) - I::point(top)) + I::point(thread.hi);
        let component_width =
            allocate_component_width(options.maximum_width_mm3, determinant, outer)?;
        certified(thread, component_width)?;
        let head_limit = fitting_limit(
            options.maximum_head_cells,
            available,
            RepairedHeadVolume::scratch_budget_bytes,
        );
        let head_prefix = self.head.work_used();
        let head_maximum = head_prefix
            .checked_add(work.remaining)
            .ok_or(PhillipsReadError::Limit)?
            .min(16_000_000);
        let (head, head_cells, head_work) = {
            let result = self
                .head
                .volume_bounds(
                    RepairedHeadVolumeOptions {
                        maximum_width_mm3: component_width,
                        maximum_cells: head_limit,
                        maximum_work: head_maximum,
                        maximum_bytes: options.maximum_bytes,
                    },
                    interrupt,
                )
                .map_err(super::super::super::super::closed_head::repaired_error)?;
            (
                interval(result.volume_mm3()),
                result.cells(),
                result
                    .work_used()
                    .checked_sub(head_prefix)
                    .ok_or(PhillipsReadError::Invalid)?,
            )
        };
        work.step(head_work).map_err(|e| bound_error(e, &work))?;
        let cavity_max = options.maximum_cavity_evaluations.min(work.remaining);
        if cavity_max == 0 {
            return Err(PhillipsReadError::Limit);
        }
        let cavity_limit =
            fitting_limit(cavity_max, available, Cutter::volume_scratch_budget_bytes);
        let result = self
            .cavity
            .source()
            .volume_bounds_with_cancel(component_width, cavity_limit, &|| {
                interrupt.is_some_and(|f| f())
            })
            .map_err(cavity_error)?;
        let cavity = interval(result.volume_mm3);
        let cavity_evaluations = result.section_evaluations;
        work.step(cavity_evaluations)
            .map_err(|e| bound_error(e, &work))?;
        let mut local = head + thread - cavity;
        // Earned axial cavity exclusion leaves the entire thread occupied.
        local.lo = local.lo.max(thread.lo);
        local.hi = local.hi.min(outer.hi);
        let world = certified(local * determinant, options.maximum_width_mm3)?;
        work.step(1).map_err(|e| bound_error(e, &work))?;
        let peak = prefix
            + STREAM_SCRATCH
                .max(RepairedHeadVolume::scratch_budget_bytes(head_limit))
                .max(Cutter::volume_scratch_budget_bytes(cavity_limit));
        Ok(Degree11MaterialVolume {
            source: self,
            local: [local.lo, local.hi],
            world,
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
            new_work: options.maximum_work - work.remaining - self.work_used(),
            peak,
            normals,
        })
    }
}

// Allocate the requested world interval width after the unavoidable affine
// determinant uncertainty. Pure numerical contract, independent of receipts.
fn allocate_component_width(
    width: f64,
    determinant: I,
    outer: I,
) -> std::result::Result<f64, PhillipsReadError> {
    let reserve = ((I::point(determinant.hi) - I::point(determinant.lo)) * outer).hi;
    let remaining = I::point(width) - I::point(reserve);
    if !remaining.lo.is_finite() || remaining.lo <= 0. {
        return Err(PhillipsReadError::Unsupported);
    }
    Ok((I::point(remaining.lo) / (I::point(8.) * I::point(determinant.hi))).lo)
}

mod owned;
pub use owned::{
    Degree11EndDiskOptions, Degree11EndDiskProperties, Degree11EndDisks, Degree11FaceMeasurement,
    Degree11FacePropertyOptions, Degree11MaterialQuery, Degree11MaterialQueryError,
    Degree11MaterialQueryOptions, Degree11OwnedMaterialOptions, Degree11ReplayMatch,
    Degree11ReplayMatchOptions, OwnedDegree11Material,
};
