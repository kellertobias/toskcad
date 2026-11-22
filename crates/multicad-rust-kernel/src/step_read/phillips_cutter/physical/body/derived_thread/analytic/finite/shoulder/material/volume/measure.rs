//! Component integration shared by imported and directly authored sources.
//! This routine proves numeric enclosures only. Its caller separately earns
//! whole-source closure/containment and credits the already-paid head prefix.
use super::*;
pub(super) struct MeasureInput<'q, 'p> {
    pub head: &'q crate::head_profile::RepairedHeadMeridian<'p>,
    pub cavity: &'q Cutter,
    pub thread: &'q PeriodicThreadSource,
    pub stations: [[f64; 2]; 2],
    pub determinant: [f64; 2],
}
pub(super) struct Measured {
    pub local: [f64; 2],
    pub world: [f64; 2],
    pub components: [[f64; 2]; 3],
    pub head_cells: usize,
    pub cavity_evaluations: usize,
    pub head_work: usize,
    pub thread_work: usize,
}
pub(super) fn measure(
    input: MeasureInput<'_, '_>,
    o: PeriodicVolumeOptions,
    available: usize,
    work: &mut Work<'_>,
) -> FResult<Measured> {
    let before = work.remaining;
    let thread = periodic_volume(input.thread, input.stations, work)?;
    let thread_work = before - work.remaining;
    let determinant = interval(input.determinant);
    let preliminary =
        (I::point(o.maximum_width_mm3) / (I::point(8.) * I::point(determinant.hi))).lo;
    certified(thread, preliminary)?;
    let head = input.head;
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
    let outer =
        pi * I::point(radius).square() * (I::point(end) - I::point(top)) + I::point(thread.hi);
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
            work.interrupt,
        )
        .map_err(crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error)?;
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
    let cavity_limit = fitting_limit(cavity_max, available, Cutter::volume_scratch_budget_bytes);
    let c = input
        .cavity
        .volume_bounds_with_cancel(component_width, cavity_limit, &|| {
            work.interrupt.is_some_and(|f| f())
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
    Ok(Measured {
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
    })
}
