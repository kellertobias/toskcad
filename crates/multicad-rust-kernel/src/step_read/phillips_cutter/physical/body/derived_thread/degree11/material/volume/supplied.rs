//! Complete supplied-author/material correspondence, distinct from raw identity.
//! Both original rounded-root laws and original D11 coefficients remain intact.
use super::super::super::super::analytic::{AnalyticThreadOptions, SuppliedThreadCorrespondence};
use super::*;
use crate::{head_profile::RepairedHeadOptions, phillips_drive::PhillipsScrew};
mod chart;
type PResult<T> = std::result::Result<T, PhillipsReadError>;
const EXTRA_SCRATCH: usize = 4 * 1024 * 1024;
const ATLAS_BYTES: usize = 1024 * 1024;
const HEAD_BYTES: usize = 1024 * 1024;
// Static primitive bounds for the ACTUAL order8 paired chart implementation.
// Series9: add72/sub90/mul1536/recip1536/div3072/sqrt1536;
// constant multiply1024/constant reciprocal768. Entry-roundover sine inputs
// have only c0,c1: their full sin/cos recurrence (including all zero checks)
// costs<=2304. reduced entry/envelope therefore costs<39,000 per evaluation.
// Four evaluations, three paired degree7+remainder differences, both complete
// derivative lifts, validity/zero checks, array/DFS/callback bookkeeping fit
// under163840. Rounded inner walls add two square roots/divisions per reduced
// cone expression;196608 covers BOTH operands and those extra coefficient
// loops. No chart constructor allocation occurs inside a Taylor attempt.
fn cavity_cell_work(a: &Cutter, b: &Cutter) -> usize {
    if [a, b]
        .into_iter()
        .any(|c| matches!(c.spec().inner, crate::phillips_drive::Inner::Rounded { .. }))
    {
        196_608
    } else {
        163_840
    }
}
#[derive(Clone, Copy)]
pub struct SuppliedDegree11MaterialOptions {
    pub maximum_distance_mm: f64,
    /// Each corresponding source chart, including original thread unit cells.
    pub maximum_derivative_mm: f64,
    pub maximum_volume_discrepancy_mm3: f64,
    /// Complete already-earned decoded volume plus fresh comparison work.
    /// Prior author construction is excluded, its reachable storage is not.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub maximum_cavity_cells: usize,
}
impl Default for SuppliedDegree11MaterialOptions {
    fn default() -> Self {
        Self {
            maximum_distance_mm: 1e-7,
            maximum_derivative_mm: 1e-7,
            maximum_volume_discrepancy_mm3: 1e-5,
            maximum_work: WORK_LIMIT,
            maximum_bytes: MAX_BYTES,
            maximum_cavity_cells: 262_144,
        }
    }
}
pub struct SuppliedDegree11MaterialMatch<'m, 'v, 's, 'author, 'a, 'input> {
    decoded: &'m Degree11MaterialVolume<'v, 's, 'a, 'input>,
    thread: SuppliedThreadCorrespondence<'author, 'a, 'input>,
    distance: f64,
    derivative: f64,
    volume_difference: [f64; 2],
    proposal: [f64; 2],
    fresh_work: usize,
    peak: usize,
    retained: usize,
    cavity_cells: usize,
    cavity_accuracy: [f64; 2],
}
impl SuppliedDegree11MaterialMatch<'_, '_, '_, '_, '_, '_> {
    pub fn author(&self) -> &PhillipsScrew {
        self.thread.source()
    }
    pub fn decoded_volume(&self) -> &Degree11MaterialVolume<'_, '_, '_, '_> {
        self.decoded
    }
    pub fn supplied_thread(&self) -> &SuppliedThreadCorrespondence<'_, '_, '_> {
        &self.thread
    }
    pub fn maximum_world_distance_mm(&self) -> f64 {
        self.distance
    }
    pub fn maximum_world_chart_derivative_mm(&self) -> f64 {
        self.derivative
    }
    /// Author minus independently decoded material. Overlap is not required.
    pub fn independent_volume_difference_mm3(&self) -> [f64; 2] {
        self.volume_difference
    }
    /// Earned original author-head proposal to its repaired-head C0/C1.
    pub fn original_head_proposal_bounds_mm(&self) -> [f64; 2] {
        self.proposal
    }
    pub fn new_work_used(&self) -> usize {
        self.fresh_work
    }
    pub fn total_work_used(&self) -> usize {
        self.decoded.work_used() + self.fresh_work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
    pub fn cavity_cells_used(&self) -> usize {
        self.cavity_cells
    }
    /// Actual paired C0/C1 allocation before the complete downstream ray lift.
    pub fn cavity_accuracy_mm(&self) -> [f64; 2] {
        self.cavity_accuracy
    }
}
fn compatible(ok: bool) -> PResult<()> {
    if ok {
        Ok(())
    } else {
        Err(PhillipsReadError::Unsupported)
    }
}
fn step(work: &mut Work<'_>, n: usize) -> PResult<()> {
    work.step(n).map_err(|e| bound_error(e, work))
}
fn maximum(values: impl IntoIterator<Item = I>) -> I {
    let mut upper = 0_f64;
    for value in values {
        if !value.valid() {
            return I::new(f64::NEG_INFINITY, f64::INFINITY);
        }
        upper = upper.max(value.hi);
    }
    I::point(upper)
}
fn abs(x: I) -> I {
    I::point(x.magnitude())
}
fn valid(o: SuppliedDegree11MaterialOptions) -> bool {
    [o.maximum_distance_mm, o.maximum_derivative_mm]
        .into_iter()
        .all(|x| x.is_finite() && x > 0. && x <= 1e-7)
        && o.maximum_volume_discrepancy_mm3.is_finite()
        && o.maximum_volume_discrepancy_mm3 > 0.
        && o.maximum_volume_discrepancy_mm3 <= 1e-5
        && o.maximum_work <= WORK_LIMIT
        && o.maximum_bytes <= MAX_BYTES
        && o.maximum_cavity_cells > 0
        && o.maximum_cavity_cells <= 262_144
}
struct ThreadBounds {
    delta: I,
    derivative: [I; 2],
    speed: [I; 2],
    rims: [RimBounds; 4],
}
#[derive(Clone, Copy)]
struct RimBounds {
    delta: I,
    derivative: I,
    speed: I,
}
fn thread_bounds(
    material: &Degree11ClosedMaterial<'_, '_, '_>,
    supplied: &SuppliedThreadCorrespondence<'_, '_, '_>,
    world: &chart::ChartBounds,
    extent: I,
    work: &mut Work<'_>,
) -> PResult<ThreadBounds> {
    let faces = material.thread().carrier().faces();
    compatible(
        faces.len() == 4
            && faces.len() == supplied.transfers().len()
            && faces.len() == supplied.entry_rim_transfers().len(),
    )?;
    let mut delta = I::point(0.);
    let mut derivative = [I::point(0.); 2];
    let mut speed = [I::point(0.); 2];
    let mut rims = [RimBounds {
        delta: I::point(0.),
        derivative: I::point(0.),
        speed: I::point(0.),
    }; 4];
    for (index, (face, author)) in faces.iter().zip(supplied.transfers()).enumerate() {
        let transfer = face.transfer();
        compatible(
            transfer.face_entity == author.face_entity
                && transfer.surface_entity == author.surface_entity
                && transfer.original_rail_entities == author.original_rail_entities
                && transfer.original_v_start_cell == author.original_v_start_cell
                && transfer.original_cells[0] == author.original_cells[0]
                && transfer.original_cells[1] + transfer.original_v_start_cell
                    == author.original_cells[1],
        )?;
        let source = face.source();
        // Each retained pole has at most two outgoing degree11 differences.
        // 512 steps cover interval vector differences/scaling/norm and loops.
        step(work, source.poles().len().saturating_mul(512))?;
        let [width, height] = source.dimensions();
        let mut face_speed = [I::point(0.); 2];
        let mut rim_speed = I::point(0.);
        for u in 0..width {
            for v in 0..height {
                let p = source.poles()[u * height + v];
                for (axis, next) in [
                    (0, (u + 1 < width).then_some((u + 1) * height + v)),
                    (1, (v + 1 < height).then_some(u * height + v + 1)),
                ] {
                    if let Some(next) = next {
                        let q = source.poles()[next];
                        let d = chart::norm(
                            (0..3).map(|i| I::point(11.) * (I::point(q[i]) - I::point(p[i]))),
                        );
                        face_speed[axis] = maximum([face_speed[axis], d]);
                        if axis == 0 && v == 0 {
                            rim_speed = maximum([rim_speed, d]);
                        }
                    }
                }
            }
        }
        // World source/raw error plus actual D11/raw transfer, then the
        // actual chart difference: no world-mm bound is treated as local-mm.
        delta = maximum([
            delta,
            world.inverse_author_norm
                * (I::point(author.distance_mm)
                    + I::point(transfer.distance_mm)
                    + world.translation
                    + world.linear * extent),
        ]);
        for axis in 0..2 {
            derivative[axis] = maximum([
                derivative[axis],
                world.inverse_author_norm
                    * (I::point(author.derivative_mm[axis])
                        + I::point(transfer.derivative_mm[axis])
                        + world.linear * face_speed[axis]),
            ]);
            speed[axis] = maximum([speed[axis], face_speed[axis]]);
        }
        let rim = supplied.entry_rim_transfers()[index];
        compatible(
            rim.face_entity == transfer.face_entity
                && rim.surface_entity == transfer.surface_entity
                && rim.original_rail_entities == transfer.original_rail_entities
                && rim.original_v_start_cell == transfer.original_v_start_cell
                && rim.original_u_cells == transfer.original_cells[0],
        )?;
        rims[index] = RimBounds {
            delta: world.inverse_author_norm
                * (I::point(rim.distance_mm)
                    + I::point(transfer.distance_mm)
                    + world.translation
                    + world.linear * extent),
            derivative: world.inverse_author_norm
                * (I::point(rim.derivative_mm)
                    + I::point(transfer.derivative_mm[0])
                    + world.linear * rim_speed),
            speed: rim_speed,
        };
    }
    Ok(ThreadBounds {
        delta,
        derivative,
        speed,
        rims,
    })
}
struct CavityBudget<'w> {
    world: &'w chart::ChartBounds,
    extent: I,
    other_target_speed: I,
    options: SuppliedDegree11MaterialOptions,
}
struct CavityBounds {
    position: I,
    derivative: I,
    top_position: I,
    top_derivative: I,
    cells: usize,
    speed: I,
    accuracy: [f64; 2],
}
fn whole_patch_speeds(
    atlas: &crate::phillips_drive::PatchAtlas<'_>,
    work: &mut Work<'_>,
) -> PResult<[I; 2]> {
    let mut maximum_speed = [I::point(0.); 2];
    for i in 0..atlas.patches().len() {
        // Fixed second-order J state, source law and normal; no adaptive
        // storage.65536 covers all branch arithmetic/checks and fixed reads.
        step(work, 65536)?;
        let bound = atlas
            .patch_bounds(i, [[0., 1.], [0., 1.]])
            .map_err(cavity_error)?;
        for (axis, derivative) in [bound.du_mm, bound.dv_mm].into_iter().enumerate() {
            let speed = derivative
                .into_iter()
                .fold(I::point(0.), |s, x| s + abs(interval(x)));
            compatible(speed.valid())?;
            maximum_speed[axis] = maximum([maximum_speed[axis], speed]);
        }
    }
    Ok(maximum_speed)
}
fn cavity_allocation(
    budget: &CavityBudget<'_>,
    rho: I,
    target_radius: I,
    radius_delta: I,
    mouth_shift: I,
    speeds: [I; 2],
) -> PResult<[f64; 2]> {
    compatible(rho.valid() && rho.lo > 0.)?;
    let target_speed = maximum([
        budget.other_target_speed,
        speeds[0],
        speeds[1],
        target_radius * speeds[0] / rho,
        I::point(2.) * target_radius,
    ]);
    let norm = I::point(budget.world.author_norm.hi);
    let position_budget = (I::point(budget.options.maximum_distance_mm)
        - budget.world.translation
        - budget.world.linear * budget.extent)
        / norm;
    let derivative_budget = (I::point(budget.options.maximum_derivative_mm)
        - budget.world.linear * target_speed)
        / norm;
    let position_free = position_budget - radius_delta - mouth_shift;
    let derivative_free =
        derivative_budget - radius_delta * maximum([I::point(1.), speeds[0] / rho]);
    compatible(
        position_free.valid()
            && position_free.lo > 0.
            && derivative_free.valid()
            && derivative_free.lo > 0.,
    )?;
    // Top C0 <= dr+dz+max(1,2R/rho)*ep.
    // Top dv <= dr+(1+2R/rho)*ep.
    // Top du <= dr*M/rho +(R+dr)*ed/rho +5R*M*ep/rho²:
    // the paired receipt also earns sourceSpeed<=targetSpeed+ed.
    // Half the residual goes to ed, a quarter to both ep contributions;
    // the final complete actual world C0/D1 gates are never bypassed.
    let position_amplification = maximum([I::point(1.), I::point(2.) * target_radius / rho]);
    let derivative_position_amplification = I::point(1.)
        + I::point(2.) * target_radius / rho
        + I::point(5.) * target_radius * speeds[0] / rho.square();
    let derivative_amplification = maximum([I::point(1.), (target_radius + radius_delta) / rho]);
    let ep = (I::point(position_free.lo) / (I::point(2.) * position_amplification))
        .lo
        .min((I::point(derivative_free.lo) / (I::point(4.) * derivative_position_amplification)).lo)
        .min(1e-7);
    let ed = (I::point(derivative_free.lo) / (I::point(2.) * derivative_amplification))
        .lo
        .min((I::point(derivative_budget.lo) / I::point(2.)).lo)
        .min(1e-7);
    compatible([ep, ed].into_iter().all(|x| x.is_finite() && x > 0.))?;
    Ok([ep, ed])
}
fn cavity_maps(
    author: &PhillipsScrew,
    material: &Degree11ClosedMaterial<'_, '_, '_>,
    source_head: &RepairedHeadMeridian<'_>,
    budget: CavityBudget<'_>,
    work: &mut Work<'_>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> PResult<CavityBounds> {
    // Each constructor has <=128nodes/256edges/88patches. Both full payloads
    // and fixed source construction are reserved before either allocation.
    step(work, 1_000_000)?;
    let a = author
        .cavity()
        .source()
        .surface_patches_with_cancel(&|| interrupt.is_some_and(|f| f()))
        .map_err(cavity_error)?;
    step(work, 1_000_000)?;
    let b = material
        .cavity()
        .surface_patches_with_cancel(&|| interrupt.is_some_and(|f| f()))
        .map_err(cavity_error)?;
    if a.retained_budget_bytes() > ATLAS_BYTES || b.retained_budget_bytes() > ATLAS_BYTES {
        return Err(PhillipsReadError::Limit);
    }
    step(work, 32768)?;
    let source_top = source_head.station(0).unwrap();
    let target_top = material.head().station(0).unwrap();
    let mouth_shift = abs(I::point(source_top[1])
        - I::point(author.spec().shank_length_mm)
        - I::point(target_top[1]));
    let rho = I::point(a.radial_margin_mm_bounds()[0].min(b.radial_margin_mm_bounds()[0]));
    let radius_delta = abs(I::point(source_top[0]) - I::point(target_top[0]));
    let speeds = whole_patch_speeds(&b, work)?;
    let accuracy = cavity_allocation(
        &budget,
        rho,
        I::point(target_top[0]),
        radius_delta,
        mouth_shift,
        speeds,
    )?;
    // Reserve the actual-family operation upper bound before any attempt;
    // independently capped Taylor attempts are not treated as free work.
    let cell_work = cavity_cell_work(author.cavity().source(), material.cavity());
    let maximum_cells = budget
        .options
        .maximum_cavity_cells
        .min(work.remaining / cell_work);
    if maximum_cells == 0 {
        return Err(PhillipsReadError::Limit);
    }
    let c = a
        .source_correspondence(
            &b,
            crate::phillips_drive::PatchCorrespondenceOptions {
                maximum_position_mm: accuracy[0],
                maximum_derivative_mm: accuracy[1],
                maximum_cells,
            },
            interrupt,
        )
        .map_err(cavity_error)?;
    step(work, c.cells_used * cell_work)?;
    let mouth_delta = I::point(c.position_mm);
    let mouth_du = I::point(c.du_mm);
    // Both independent bounds cover the SAME complete source charts. Taking
    // their minimum is safe; the paired C1 also bounds source speed by M+ed.
    let target_du = I::point(speeds[0].hi.min(c.maximum_target_du_mm));
    let target_dv = I::point(speeds[1].hi.min(c.maximum_target_dv_mm));
    let source_du = I::point(c.maximum_source_du_mm.min((target_du + mouth_du).hi));
    let direction = chart::ray(mouth_delta, mouth_du, target_du, rho)?;
    let outer_delta = radius_delta + I::point(target_top[0]) * direction[0];
    let outer_du = radius_delta * source_du / rho + I::point(target_top[0]) * direction[1];
    let top_position = maximum([mouth_delta, outer_delta]) + mouth_shift;
    let top_derivative = maximum([mouth_du, outer_du, mouth_delta + outer_delta]);
    let target_speed = maximum([
        target_du,
        target_dv,
        I::point(target_top[0]) * target_du / rho,
        I::point(2.) * I::point(target_top[0]),
    ]);
    Ok(CavityBounds {
        position: mouth_delta + mouth_shift,
        derivative: maximum([I::point(c.du_mm), I::point(c.dv_mm)]),
        top_position,
        top_derivative,
        cells: c.cells_used,
        speed: target_speed,
        accuracy,
    })
}
impl<'v, 's, 'a, 'input> Degree11MaterialVolume<'v, 's, 'a, 'input> {
    /// Match every complete physical face to this actual supplied author.
    /// No manufacturing capability or exact raw-file identity is granted.
    pub fn match_supplied_author<'m, 'author>(
        &'m self,
        author: &'author PhillipsScrew,
        frame: [Point; 3],
        options: SuppliedDegree11MaterialOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> PResult<SuppliedDegree11MaterialMatch<'m, 'v, 's, 'author, 'a, 'input>> {
        if !valid(options) {
            return Err(PhillipsReadError::Invalid);
        }
        let allowance = options
            .maximum_work
            .checked_sub(self.work_used())
            .ok_or(PhillipsReadError::Limit)?;
        let mut work = Work {
            remaining: allowance,
            interrupt,
            cancelled: Cell::new(false),
        };
        step(&mut work, 1)?;
        let material = self.source();
        let origin = material.thread().carrier().originating_head();
        // The COMPLETE originating receipt is borrowed by the material. Only
        // this exact shared allocation/proof prefix is subtracted, never a
        // separately reconstructed equivalent head or its buffers/options.
        let outer_resident = material
            .retained_budget_bytes()
            .checked_add(std::mem::size_of_val(self))
            .and_then(|n| {
                n.checked_add(std::mem::size_of::<
                    SuppliedDegree11MaterialMatch<'_, '_, '_, '_, '_, '_>,
                >())
            })
            .ok_or(PhillipsReadError::Limit)?;
        let other_resident = outer_resident
            .checked_sub(origin.retained_budget_bytes())
            .ok_or(PhillipsReadError::Limit)?;
        let child_bytes = options
            .maximum_bytes
            .checked_sub(other_resident)
            .and_then(|n| n.checked_sub(EXTRA_SCRATCH))
            .ok_or(PhillipsReadError::Limit)?;
        let thread = origin.supplied_thread_correspondence(
            author,
            frame,
            AnalyticThreadOptions {
                maximum_bytes: child_bytes,
                maximum_work: work.remaining.min(WORK),
                maximum_total_work: origin
                    .work_used()
                    .saturating_add(work.remaining)
                    .min(32_000_000),
                ..Default::default()
            },
            interrupt,
        )?;
        step(&mut work, thread.new_work_used())?;
        let retained = other_resident.saturating_add(thread.retained_budget_bytes());
        if retained > options.maximum_bytes {
            return Err(PhillipsReadError::Limit);
        }
        // EXTRA_SCRATCH precharges one repaired author head, two bounded
        // canonical atlases, the fixed64-cell Taylor stack and cell temporaries.
        let peak = other_resident
            .saturating_add(child_bytes)
            .saturating_add(EXTRA_SCRATCH);
        step(&mut work, 32768)?;
        let world = chart::compare(frame, material.chart())?;
        let source_head = author
            .profile()
            .closed_repaired_meridian(
                RepairedHeadOptions {
                    maximum_work: work.remaining.min(WORK),
                    maximum_bytes: HEAD_BYTES,
                    ..Default::default()
                },
                interrupt,
            )
            .map_err(
                crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error,
            )?;
        step(&mut work, source_head.work_used())?;
        compatible(source_head.retained_budget_bytes() <= HEAD_BYTES)?;
        let proposal = source_head
            .proposal_transfer(work.remaining.min(WORK), interrupt)
            .map_err(
                crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error,
            )?;
        step(&mut work, proposal.work_used)?;
        let target_head = material.head();
        let h = source_head
            .translated_correspondence(
                target_head,
                author.spec().shank_length_mm,
                0.,
                work.remaining.min(WORK),
                interrupt,
            )
            .map_err(
                crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error,
            )?;
        step(&mut work, h.work_used)?;
        let head_position = I::point(proposal.position_mm) + I::point(h.position_mm);
        let head_derivative = I::point(proposal.derivative_mm) + I::point(h.derivative_mm);
        let mut radius = I::point(0.);
        let mut z = I::point(0.);
        let mut target_speed = I::point(0.);
        for i in 0..target_head.segment_count() {
            step(&mut work, 4096)?;
            let bound = target_head.certified_interval_bounds(i, [0., 1.]);
            radius = maximum([radius, abs(interval(bound.position_mm[0]))]);
            z = maximum([z, abs(interval(bound.position_mm[1]))]);
            target_speed = maximum([
                target_speed,
                abs(interval(bound.derivative_mm[0])) + abs(interval(bound.derivative_mm[1])),
            ]);
        }
        let thread_box = material.thread().local_enclosure_bounds_mm();
        let extent = maximum([
            chart::norm([radius, radius, z]),
            chart::norm(thread_box.map(|x| I::point(x[0].abs().max(x[1].abs())))),
        ]);
        let t = thread_bounds(material, &thread, &world, extent, &mut work)?;
        let source_radius = if let Some(root) = author.thread().root_circle() {
            I::point(root.center_radius) - I::point(root.radius)
        } else {
            I::point(author.thread().spec().major_radius_mm)
                - I::point(author.thread().spec().depth_mm)
        };
        let source_bottom = source_head.station(source_head.segment_count()).unwrap();
        let target_bottom = target_head.station(target_head.segment_count()).unwrap();
        let radius_delta = abs(I::point(source_bottom[0]) - I::point(target_bottom[0]));
        let cap_delta = abs(I::point(source_bottom[1])
            - I::point(author.spec().shank_length_mm)
            - I::point(material.thread().cap_stations_mm()[0]));
        let mut join_position = I::point(0.);
        let mut join_du = I::point(0.);
        let mut join_dv = I::point(0.);
        let mut join_speed = I::point(0.);
        for rim in t.rims {
            // The whole source/rim correspondence proves target norm >=
            // author minimum radius minus C0, sharper than a separately
            // sufficient occupied-core cylinder. Both norms are still earned.
            let rho = I::point((source_radius - rim.delta).lo);
            let ray = chart::ray(rim.delta, rim.derivative, rim.speed, rho)?;
            let outer_delta = radius_delta + I::point(target_bottom[0]) * ray[0] + cap_delta;
            let outer_du = radius_delta * (rim.speed + rim.derivative) / rho
                + I::point(target_bottom[0]) * ray[1];
            join_position = maximum([join_position, rim.delta, outer_delta]);
            join_du = maximum([join_du, rim.derivative, outer_du]);
            join_dv = maximum([join_dv, rim.delta + outer_delta]);
            join_speed = maximum([join_speed, I::point(target_bottom[0]) * rim.speed / rho]);
        }
        // Complete disk((1-v)*rimXY(u),capZ), in each own local chart. Both
        // actual D11 cap rows were pinned BEFORE its measured raw transfer.
        // capZs may differ: their constant difference cancels from dv, leaving
        // only -deltaXY, whose norm<=full rim delta (not twice that delta).
        // du is the contracted rim derivative; all cap corrections are already
        // included in the full supplied/raw plus D11/raw C0/Du transfers.
        let tip_position = t.delta;
        let tip_derivative = maximum([t.derivative[0], t.delta]);
        let tau = I::point(2.)
            * I::new(
                std::f64::consts::PI.next_down(),
                std::f64::consts::PI.next_up(),
            );
        let other_speed = maximum([
            target_speed,
            tau * radius,
            radius,
            t.speed[0],
            t.speed[1],
            join_speed,
            I::point(target_bottom[0]) + radius,
        ]);
        let cavity = cavity_maps(
            author,
            material,
            &source_head,
            CavityBudget {
                world: &world,
                extent,
                other_target_speed: other_speed,
                options,
            },
            &mut work,
            interrupt,
        )?;
        let local_position = maximum([
            head_position,
            t.delta,
            join_position,
            tip_position,
            cavity.position,
            cavity.top_position,
        ]);
        let local_derivative = maximum([
            head_derivative,
            tau * head_position,
            head_position,
            t.derivative[0],
            t.derivative[1],
            join_du,
            join_dv,
            tip_derivative,
            cavity.derivative,
            cavity.top_derivative,
        ]);
        let target_speed = maximum([other_speed, cavity.speed]);
        let distance =
            world.translation + world.linear * extent + world.author_norm * local_position;
        let derivative = world.linear * target_speed + world.author_norm * local_derivative;
        // Fixed analytic head/root integrals, cached cavity bounds and bounded
        // retained-period precision loops; no mesh volume or nominal scalar.
        step(&mut work, 1_000_000)?;
        let author_mass = author
            .source_volume_bounds()
            .map_err(cavity_error)?
            .volume_mm3();
        let difference =
            interval(author_mass) * world.author_determinant - interval(self.world_volume_mm3());
        compatible(
            [distance, derivative, difference]
                .into_iter()
                .all(|x| x.valid() && x.lo.is_finite() && x.hi.is_finite())
                && distance.hi <= options.maximum_distance_mm
                && derivative.hi <= options.maximum_derivative_mm
                && difference.magnitude() <= options.maximum_volume_discrepancy_mm3,
        )?;
        step(&mut work, 1)?;
        Ok(SuppliedDegree11MaterialMatch {
            decoded: self,
            thread,
            distance: distance.hi,
            derivative: derivative.hi,
            volume_difference: [difference.lo, difference.hi],
            proposal: [proposal.position_mm, proposal.derivative_mm],
            fresh_work: allowance - work.remaining,
            peak,
            retained,
            cavity_cells: cavity.cells,
            cavity_accuracy: cavity.accuracy,
        })
    }
}
#[cfg(test)]
mod tests;
