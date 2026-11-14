//! Actual regular screw with a general measured clock, never a fitted source.
//! The similarity is an auxiliary enclosure of the literal affine map only.
use super::*;
use crate::headed_screw::Screw;

#[derive(Clone, Copy)]
pub struct ClockedMatingOptions {
    pub requirements: MatingRequirements,
    /// Logical whole-phase cell evaluations, including rejected parents.
    pub maximum_cells: usize,
    pub linear_tolerance_mm: f64,
}
pub struct ClockedInsertMating {
    pub mating: InsertMating,
    pub radial_map_error_mm: f64,
    pub axial_map_error_mm: f64,
    pub angle_error_rad: f64,
    pub cells_used: usize,
}
struct Alignment {
    projection: Projection,
    scale: I,
    angle: f64,
    ray_error: f64,
    radial_error: f64,
    axial_error: f64,
}
fn exact(v: &Q) -> R<I> {
    let f = v.to_f64().ok_or(Error::Unsupported)?;
    if q(f)? == *v { Ok(ip(f)) } else { enclosure(v) }
}
fn tau() -> I {
    ip(2.)
        * I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        )
}
fn abs(v: &Q) -> R<I> {
    exact(&v.abs())
}
// Only these constructors supply boxes: all dimensions come from admitted
// immutable sources. No generic numeric box can enter Projection::enclose.
struct SourceBox([[f64; 2]; 3]);
impl SourceBox {
    fn prefix(s: &Screw) -> R<Self> {
        Ok(Self(prefix_bounds(s)?))
    }
    fn thread(s: &crate::helical_thread::Thread) -> Self {
        Self(s.source_coordinate_bounds_mm())
    }
    fn tip(s: &crate::helical_thread::Thread) -> Self {
        let b = s.source_coordinate_bounds_mm();
        Self([b[0], b[1], [s.spec().length_mm, s.spec().length_mm]])
    }
    fn cavity(s: &Insert) -> Self {
        let r = s.inner().maximum_radius_mm();
        Self([[-r, r], [-r, r], [0., s.full_length_mm()]])
    }
}
struct Projection {
    origin: [Q; 3],
    columns: [[Q; 3]; 3],
    origin_i: [I; 3],
    columns_i: [[I; 3]; 3],
    metric: Metric,
}
impl Projection {
    fn new(from: Frame, to: Frame, interrupt: &dyn Fn() -> bool) -> R<Self> {
        // Validation is owned HERE, not a caller convention. These exact
        // original frames bound cofactor ratios/entries below2, origin below
        // 4e7, and retain a positive near-unit target metric.
        scale(from)?;
        let metric = scale(to)?;
        let (origin, columns) = affine(to, from, interrupt)?;
        let mut origin_i = [ip(0.); 3];
        let mut columns_i = [[ip(0.); 3]; 3];
        for i in 0..3 {
            origin_i[i] = exact(&origin[i])?;
            for j in 0..3 {
                columns_i[j][i] = exact(&columns[j][i])?;
            }
        }
        Ok(Self {
            origin,
            columns,
            origin_i,
            columns_i,
            metric,
        })
    }
    fn enclose(&self, bounds: SourceBox) -> [I; 3] {
        // Head dimensions <=1e7; the whole retained root-circle hull <1e9;
        // insert dimensions <=1e7. Together with constructor bounds every
        // product and accumulated sum is <1e11 and finite. Each valid input
        // interval is multiplied/added outward; no repeated Q conversion.
        std::array::from_fn(|i| {
            (0..3).fold(self.origin_i[i], |sum, j| {
                sum + self.columns_i[j][i] * I::new(bounds.0[j][0], bounds.0[j][1])
            })
        })
    }
}
fn prefix_bounds(screw: &Screw) -> R<[[f64; 2]; 3]> {
    let radius = head_radii(screw)?.1.max(screw.spec().shank_radius_mm);
    let bottom = finite(-ip(screw.spec().shank_length_mm) - ip(screw.spec().head.height_mm()))?.lo;
    Ok([[-radius, radius], [-radius, radius], [bottom, 0.]])
}
/// Certify an actual ray enclosed by c/s and its positive norm enclosure.
/// The proposed atan2 value is checked by directed trig; uncertain rays refuse.
fn ray_angle(c: I, s: I, scale: I) -> R<(f64, f64)> {
    let c = finite(c)?;
    let s = finite(s)?;
    let scale = positive(scale)?;
    let angle = s.midpoint().atan2(c.midpoint());
    let chord = finite(
        ip((c / scale - ip(angle).cos()).magnitude())
            + ip((s / scale - ip(angle).sin()).magnitude()),
    )?;
    if chord.hi > 1. {
        return Err(Error::Unsupported);
    }
    Ok((angle, (ip(2.) * chord).hi))
}
fn align(screw: &Screw, from: Frame, to: Frame, interrupt: &dyn Fn() -> bool) -> R<Alignment> {
    let projection = Projection::new(from, to, interrupt)?;
    let origin = &projection.origin;
    let columns = &projection.columns;
    let c = projection.columns_i[0][0];
    let s = projection.columns_i[0][1];
    let scale = positive(c.square() + s.square())?.sqrt();
    let scale = positive(scale)?;
    let (angle, ray_error) = ray_angle(c, s, scale)?;
    let zero = Q::zero();
    let one = q(1.)?;
    let ideal = [
        [columns[0][0].clone(), columns[0][1].clone(), zero.clone()],
        [-&columns[0][1], columns[0][0].clone(), zero.clone()],
        [zero.clone(), zero, one],
    ];
    let radius = ip(screw.thread().source_coordinate_bounds_mm()[0][1]);
    let length = ip(screw.thread().spec().length_mm);
    let mut error = [ip(0.); 3];
    for i in 0..3 {
        error[i] = (abs(&(&columns[0][i] - &ideal[0][i]))?
            + abs(&(&columns[1][i] - &ideal[1][i]))?)
            * radius
            + abs(&(&columns[2][i] - &ideal[2][i]))? * length;
        if i < 2 {
            error[i] = error[i] + abs(&origin[i])?;
        }
    }
    Ok(Alignment {
        projection,
        scale,
        angle,
        ray_error,
        radial_error: finite(error[0] + error[1])?.hi,
        axial_error: finite(error[2])?.hi,
    })
}
/// Complete finite-cavity mating of an actual regular screw. Partial axial
/// overlaps retain the existing quarter-clock provider; this operation refuses
/// unless the WHOLE actual insert cavity has an original male thread domain.
/// At most 4096 interval cells and a 4097-entry frontier are retained. Each
/// processed cell is charged once; finite source phase accounting is separate.
pub fn certify_clocked_insert_mating(
    insert: &Insert,
    insert_frame: Frame,
    screw: &Screw,
    screw_frame: Frame,
    options: ClockedMatingOptions,
    interrupt: &dyn Fn() -> bool,
) -> R<ClockedInsertMating> {
    checkpoint(interrupt)?;
    let req = options.requirements;
    if [req.minimum_engagement_mm, req.minimum_clearance_mm]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0. || x > 1e7)
        || options.maximum_cells == 0
        || options.maximum_cells > 4096
        || !(1e-7..=1e-3).contains(&options.linear_tolerance_mm)
    {
        return Err(Error::Invalid);
    }
    let internal = insert.inner();
    let external = screw.thread();
    let i = internal.spec();
    let e = external.spec();
    if i.pitch_mm != e.pitch_mm || i.starts != e.starts || i.left_hand != e.left_hand {
        return Err(Error::Unsupported);
    }
    let map = align(screw, screw_frame, insert_frame, interrupt)?;
    let metric = &map.projection.metric;
    if ((ip(map.radial_error) + ip(map.axial_error)) * metric.maximum_scale).hi
        > options.linear_tolerance_mm
    {
        return Err(Error::Unsupported);
    }
    // Full source head and plain interval, not a centerline prefix assertion.
    let prefix = map.projection.enclose(SourceBox::prefix(screw)?);
    if prefix[2].hi > 0. {
        return Err(Error::Unsupported);
    }
    // Every insert cavity point is inside the finite original male chart.
    // This inverse test keeps true cap domains even for a slightly tilted map.
    let inverse = Projection::new(insert_frame, screw_frame, interrupt)?;
    let l = insert.full_length_mm();
    let cavity = inverse.enclose(SourceBox::cavity(insert));
    if cavity[2].lo < 0. || cavity[2].hi > e.length_mm {
        return Err(Error::Unsupported);
    }
    let norm = vector(insert_frame.axes[2])?;
    if bounded(q(l)?.pow(2) * dot(&norm, &norm)?)? < q(req.minimum_engagement_mm)?.pow(2) {
        return Err(Error::Unsupported);
    }
    // A source core cylinder is handled without polar coordinates. The rest
    // has radius >=core, so measured point error yields a finite angle error.
    let source_min = external
        .radius_at_phase_interval(I::new(0., e.pitch_mm))
        .ok_or(Error::Unsupported)?;
    let core = positive(ip(source_min.lo) * ip(0.5))?.lo;
    let min_bore = ip(internal.minimum_radius_mm());
    if (map.scale * ip(core) + ip(map.radial_error)).hi >= min_bore.lo {
        return Err(Error::Unsupported);
    }
    let radial = positive(map.scale * ip(core))?;
    if map.radial_error > (radial * ip(0.5)).lo {
        return Err(Error::Unsupported);
    }
    let angle_error = finite(ip(map.ray_error) + ip(2.) * ip(map.radial_error) / radial)?.hi;
    let rate = ip(1.) / internal.angular_gain();
    let allowance = mating::phase_error::measure(external, interrupt)?;
    // The polar guard gives turn_error<=1 (+rounding), the ray gives <=2,
    // hence angle_error<4. Admitted source rates <1.4e7, allowance<.001,
    // and the measured map gate <.002 bound this valid sum below6e7.
    let pad =
        (ip(map.axial_error) + ip(rate.magnitude()) * ip(angle_error) + ip(allowance.total_mm)).hi;
    // Exact authored turns remove integer angular wraps. Literal external rate
    // drift and original indexed/evaluator error remain explicitly enclosed.
    // Projection owns |origin|<4e7; finite atan2 has |angle|<=pi. Together
    // with the source rate/tau and pad bounds every term/sum is <1e9.
    let offset = exact(&map.projection.origin[2])? - rate * ip(map.angle)
        + (ip(external.helical_rate()) - rate) * I::new(0., tau().hi)
        + I::new(-pad, pad);
    let mut pending = Vec::with_capacity(options.maximum_cells + 1);
    pending.push([0., e.pitch_mm]);
    let mut used = 0;
    let mut clearance = f64::INFINITY;
    let mut inner_envelope = false;
    let mut outer_envelope = false;
    while let Some([lo, hi]) = pending.pop() {
        checkpoint(interrupt)?;
        if used == options.maximum_cells {
            return Err(Error::Limit);
        }
        used += 1;
        let phase = I::new(lo, hi);
        let (bore, ib) = mating::authored_periodic_bound(internal, phase + offset);
        let (shaft, eb) = mating::external_periodic_bound(external, phase);
        let gap = finite((bore - map.scale * shaft - ip(map.radial_error)) * metric.minimum_scale)?;
        if gap.lo >= req.minimum_clearance_mm {
            clearance = clearance.min(gap.lo);
            inner_envelope |= ib;
            outer_envelope |= eb;
            continue;
        }
        let mid = lo.midpoint(hi);
        if mid <= lo || mid >= hi {
            return Err(Error::Unsupported);
        }
        pending.push([mid, hi]);
        pending.push([lo, mid]);
    }
    let length = insert_frame.axes[2]
        .map(ip)
        .into_iter()
        .fold(ip(0.), |a, b| a + b.square())
        .sqrt();
    let engagement = finite(ip(l) * length)?;
    let prefix_gap = finite(-prefix[2] * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi))?;
    let exit = map.projection.enclose(SourceBox::tip(external));
    let whole = map.projection.enclose(SourceBox::thread(external));
    let exit_gap =
        finite((ip(l) - exit[2]) * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi))?;
    checkpoint(interrupt)?;
    Ok(ClockedInsertMating {
        mating: InsertMating {
            engagement_mm: [engagement.lo.max(req.minimum_engagement_mm), engagement.hi],
            minimum_radial_clearance_lower_mm: clearance,
            source_phase_offset_mm: [offset.lo, offset.hi],
            prefix_entry_clearance_mm: [prefix_gap.lo, prefix_gap.hi],
            thread_exit_clearance_mm: [exit_gap.lo, exit_gap.hi],
            entire_thread_within_insert: whole[2].lo >= 0. && whole[2].hi <= l,
            external_indexed_phase_error_mm: allowance.indexed_mm,
            external_evaluator_phase_error_mm: allowance.evaluator_mm,
            inner_whole_period_envelope_used: inner_envelope,
            outer_whole_source_envelope_used: outer_envelope,
        },
        radial_map_error_mm: map.radial_error,
        axial_map_error_mm: map.axial_error,
        angle_error_rad: angle_error,
        cells_used: used,
    })
}

/// Full actual regular head/plain/thread versus the SAME actual blind pilot.
/// This uses a conservative transformed source box, not a point-axis fit or
/// the insert's shorter cavity. No wood-stock clearance follows from mating.
pub fn certify_clocked_insert_screw_clearance(
    stock: Stock<'_>,
    screw: &Screw,
    screw_frame: Frame,
    req: ScrewClearanceRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<ScrewClearance> {
    checkpoint(interrupt)?;
    if [
        req.minimum_tip_clearance_mm,
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    ]
    .into_iter()
    .any(|v| !v.is_finite() || !(0. ..=1e7).contains(&v))
        || !(1e-7..=1e-3).contains(&req.linear_tolerance_mm)
        || req.minimum_feature_mm == 0.
    {
        return Err(Error::Invalid);
    }
    if stock.material.cells().len() > 4096 {
        return Err(Error::Limit);
    }
    if stock.material.solid_count() != 1 {
        return Err(Error::Unsupported);
    }
    let boundary = orthogonal_boundary::build(stock.material.quads()).map_err(boundary_error)?;
    orthogonal_pockets::validate(
        stock.material,
        &boundary,
        &[stock.pilot],
        req.linear_tolerance_mm,
        req.minimum_feature_mm,
    )
    .map_err(pocket_error)?;
    let face = boundary.faces.get(stock.pilot.face).ok_or(Error::Invalid)?;
    let (axis, sign) = inward_axis(face.normal)?;
    receiving_halfspace(&stock, axis, sign, interrupt)?;
    let projection = Projection::new(screw_frame, stock.frame, interrupt)?;
    let metric = &projection.metric;
    let prefix = projection.enclose(SourceBox::prefix(screw)?);
    let prefix_depth = (prefix[axis] - ip(stock.pilot.center_mm[axis])) * ip(sign);
    if prefix_depth.hi > 0. {
        return Err(Error::Unsupported);
    }
    let thread = projection.enclose(SourceBox::thread(screw.thread()));
    let depth = (thread[axis] - ip(stock.pilot.center_mm[axis])) * ip(sign);
    let radial = [(axis + 1) % 3, (axis + 2) % 3].map(|i| thread[i] - ip(stock.pilot.center_mm[i]));
    let reach = finite(radial[0].square() + radial[1].square())?.sqrt();
    let radial_clearance = finite((ip(stock.pilot.radius_mm) - reach) * metric.minimum_scale)?;
    // Validated actual pilot dimensions and the constructor-owned source-box
    // projection keep all endpoints/products <1e11 and ordered, even when a
    // real negative floor gap must subsequently refuse.
    let floor = (ip(stock.pilot.depth_mm) - depth)
        * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi);
    if radial_clearance.lo < req.minimum_feature_mm || floor.lo < req.minimum_tip_clearance_mm {
        return Err(Error::Unsupported);
    }
    let prefix_gap =
        finite(-prefix_depth * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi))?;
    checkpoint(interrupt)?;
    Ok(ScrewClearance {
        tip_clearance_mm: [floor.lo, floor.hi],
        prefix_entry_clearance_mm: [prefix_gap.lo, prefix_gap.hi],
        minimum_radial_clearance_mm: radial_clearance.lo,
    })
}

#[cfg(test)]
mod tests;
