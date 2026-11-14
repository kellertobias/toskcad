//! Independent clearance between actual authored cavity bands and original
//! external literal flank/root supports, over the whole physical phase period.
use super::*;
pub(in crate::wood_thread_contact) mod phase_error;
use crate::{authored_flat_thread::AuthoredThread, helical_thread::Thread};
#[derive(Clone, Copy)]
pub enum MatingScrew<'a> {
    Regular(&'a crate::headed_screw::Screw),
    Profiled(&'a crate::profiled_screw::Screw),
    Phillips(&'a crate::phillips_drive::PhillipsScrew),
}
#[derive(Clone, Copy)]
pub struct MatingExternal<'a> {
    pub source: MatingScrew<'a>,
    pub frame: Frame,
}
impl<'a> MatingScrew<'a> {
    pub(super) fn parts(self) -> R<(&'a Thread, [f64; 2], f64)> {
        match self {
            Self::Regular(s) => Ok((s.thread(), [0., s.thread().spec().length_mm], 0.)),
            Self::Profiled(s) => Ok((
                s.thread(),
                s.visible_thread_interval_mm(),
                profile_prefix(s.profile(), s.spec().shank_length_mm)?,
            )),
            Self::Phillips(s) => Ok((
                s.thread(),
                s.visible_thread_interval_mm(),
                profile_prefix(s.profile(), s.spec().shank_length_mm)?,
            )),
        }
    }
}
fn profile_prefix(profile: &crate::head_profile::Profile, plain: f64) -> R<f64> {
    // A whole-circle axial support is intentionally conservative. Literal
    // support + source reserve qualifies the actual head/neck, not an inferred
    // head recipe or a nominal planar socket-floor approximation.
    let mut upper = f64::NEG_INFINITY;
    for segment in profile.segments() {
        let z = match segment.curve {
            crate::head_profile::Curve::Line => ip(segment.start[1].max(segment.end[1])),
            crate::head_profile::Curve::Arc {
                center, radius_mm, ..
            } => ip(center[1]) + ip(radius_mm),
        };
        upper = upper.max(finite(z)?.hi);
    }
    Ok(finite(ip(upper) - ip(plain) + ip(profile.rounding_error_mm()))?.hi)
}
#[derive(Clone, Copy)]
pub struct MatingRequirements {
    pub minimum_engagement_mm: f64,
    pub minimum_clearance_mm: f64,
}
#[derive(Clone, Debug)]
pub struct InsertMating {
    pub engagement_mm: [f64; 2],
    /// Uniform physical clearance lower bound, not a sampled minimum.
    pub minimum_radial_clearance_lower_mm: f64,
    pub source_phase_offset_mm: [f64; 2],
    pub prefix_entry_clearance_mm: [f64; 2],
    pub thread_exit_clearance_mm: [f64; 2],
    pub entire_thread_within_insert: bool,
    pub external_indexed_phase_error_mm: f64,
    pub external_evaluator_phase_error_mm: f64,
    /// True only when an inner bound used the earned whole-period envelope
    /// instead of resolving its narrower indexed phase interval.
    pub inner_whole_period_envelope_used: bool,
    pub outer_whole_source_envelope_used: bool,
}
fn exact_interval(v: &Q) -> R<I> {
    let x = v.to_f64().ok_or(Error::Unsupported)?;
    if q(x)? == *v { Ok(ip(x)) } else { enclosure(v) }
}
pub fn certify_insert_mating(
    insert: &Insert,
    insert_frame: Frame,
    external: MatingExternal<'_>,
    req: MatingRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<InsertMating> {
    checkpoint(interrupt)?;
    let (outer, [start, end], prefix_end_upper_mm) = external.source.parts()?;
    if [
        start,
        end,
        prefix_end_upper_mm,
        req.minimum_engagement_mm,
        req.minimum_clearance_mm,
    ]
    .into_iter()
    .any(|x| !x.is_finite() || x.abs() > 1e7)
        || start < 0.
        || start >= end
        || end > outer.spec().length_mm
        || req.minimum_engagement_mm <= 0.
        || req.minimum_clearance_mm <= 0.
    {
        return Err(Error::Invalid);
    }
    let internal = insert.inner();
    let i = internal.spec();
    let e = outer.spec();
    if i.pitch_mm != e.pitch_mm || i.starts != e.starts || i.left_hand != e.left_hand {
        return Err(Error::Unsupported);
    }
    let metric = scale(insert_frame)?;
    scale(external.frame)?;
    let (origin, columns) = relative(insert_frame, external.frame, interrupt)?;
    if columns[2] != [0, 0, 1]
        || !origin[0].is_zero()
        || !origin[1].is_zero()
        || &origin[2] + q(prefix_end_upper_mm)? > Q::zero()
    {
        return Err(Error::Unsupported);
    }
    // relative() earned a positive signed-permutation matrix and +Z axis.
    let quarter: i32 = if columns[0][0] == -1 {
        2
    } else {
        i32::from(columns[0][1])
    };
    let prefix_gap = bounded(-&origin[2] - q(prefix_end_upper_mm)?)?;
    let exit_gap = bounded(q(insert.full_length_mm())? - &origin[2] - q(end)?)?;
    let begin = (&origin[2] + q(start)?).max(Q::zero());
    let finish = (&origin[2] + q(end)?).min(q(insert.full_length_mm())?);
    let overlap = bounded(finish - begin)?;
    let norm = vector(insert_frame.axes[2])?;
    if overlap <= Q::zero()
        || bounded(overlap.pow(2) * dot(&norm, &norm)?)? < q(req.minimum_engagement_mm)?.pow(2)
    {
        return Err(Error::Unsupported);
    }
    let length_norm = insert_frame.axes[2]
        .map(ip)
        .into_iter()
        .fold(ip(0.), |s, v| s + v.square())
        .sqrt();
    let mut engagement = exact_interval(&overlap)? * length_norm;
    engagement.lo = engagement.lo.max(req.minimum_engagement_mm);
    let signed = if i.left_hand { -1. } else { 1. };
    // A quarter-turn's exact true-angle phase is signed lead/4. It does not
    // inherit independent PI/TAU rounding or silently rotate either source.
    // Exact frame leaves have common denominator D=2^1074. The inverse
    // translation's integer adjugate/determinant gives <6500 combined bits.
    // Admitted pitch>=.001 has dyadic denominator<=2^62; starts<=8 and the
    // quarter/4 term leave this expression <7000 bits, below the 16384 guard
    // already earned by relative(). No second rational-budget gate is needed.
    let base = &origin[2]
        - q(signed)?
            * q(i.pitch_mm)?
            * Q::from_integer(i.starts.into())
            * Q::from_integer(quarter.into())
            / Q::from_integer(4.into());
    let true_rate = ip(1.) / internal.angular_gain();
    let tau = ip(2.)
        * I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
    // The original external evaluator measures its actual polar angle in
    // [0,2pi]. Retain its literal helical rate and cover true-rate drift.
    let correspondence = phase_error::measure(outer, interrupt)?;
    // Positive finite overlap bounds |originZ|<=1e7; lead<=8e7 and quarter<=2
    // bound |base|<=5e7. Rates are <1.4e7, tau<7 and measure()'s admitted
    // finite-domain allowance is <.001, so this offset is finite (<2.5e8).
    // Every periodic-radius quotient is still independently checked before
    // iteration and before this interval may be published in a receipt.
    let offset = exact_interval(&base)?
        + (ip(outer.helical_rate()) - true_rate) * I::new(0., tau.hi)
        + I::new(-correspondence.total_mm, correspondence.total_mm);
    let mut clearance = f64::INFINITY;
    let mut inner_envelope_used = false;
    let mut outer_envelope_used = false;
    // Between either source's joins the cavity radius is affine and the
    // original external radius is affine/constant or a convex circular root.
    // Their difference is concave. BOTH one-sided laws and periodic endpoints
    // are evaluated by the interval radius routines at every candidate join.
    for phase in outer.profile_breaks().map(ip).into_iter().chain(
        internal
            .bands()
            .into_iter()
            .flat_map(|b| b.phase_mm)
            .map(|p| ip(p) - offset),
    ) {
        checkpoint(interrupt)?;
        let (bore, used_envelope) = authored_periodic_bound(internal, phase + offset);
        inner_envelope_used |= used_envelope;
        let (shaft, used_envelope) = external_periodic_bound(outer, phase);
        outer_envelope_used |= used_envelope;
        clearance = clearance.min(finite(bore - shaft)?.lo);
    }
    let physical = ip(clearance) * metric.minimum_scale;
    if !physical.valid() || physical.lo < req.minimum_clearance_mm {
        return Err(Error::Unsupported);
    }
    checkpoint(interrupt)?;
    Ok(InsertMating {
        engagement_mm: [engagement.lo, engagement.hi],
        minimum_radial_clearance_lower_mm: physical.lo,
        source_phase_offset_mm: [offset.lo, offset.hi],
        prefix_entry_clearance_mm: physical_gap(&prefix_gap, &metric)?,
        thread_exit_clearance_mm: physical_gap(&exit_gap, &metric)?,
        entire_thread_within_insert: exit_gap >= Q::zero(),
        external_indexed_phase_error_mm: correspondence.indexed_mm,
        external_evaluator_phase_error_mm: correspondence.evaluator_mm,
        inner_whole_period_envelope_used: inner_envelope_used,
        outer_whole_source_envelope_used: outer_envelope_used,
    })
}
pub(super) fn physical_gap(gap: &Q, metric: &Metric) -> R<[f64; 2]> {
    let v = exact_interval(gap)? * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi);
    let lower = if gap >= &Q::zero() {
        v.lo.max(0.)
    } else {
        v.lo
    };
    let upper = if gap <= &Q::zero() {
        v.hi.min(0.)
    } else {
        v.hi
    };
    Ok([lower, upper])
}

fn external_periodic_bound(source: &Thread, phase: I) -> (I, bool) {
    match periodic_radius(phase, source.spec().pitch_mm, &mut |p| {
        source.radius_at_phase_interval(p).ok_or(Error::Unsupported)
    }) {
        Ok(value) => (value, false),
        Err(_) => {
            // Actual source_support bounds retain the WHOLE literal root circle.
            // Closed flank trims are monotone from the literal major, hence
            // never exceed it. The evaluator's signed FMA flank updates cannot
            // exceed that same representable major; root evaluation subtracts
            // a nonnegative height from its literal center. Thus this upper
            // encloses both all actual supports and relevant evaluator paths.
            (I::new(0., source.source_coordinate_bounds_mm()[0][1]), true)
        }
    }
}
fn authored_periodic_bound(source: &AuthoredThread, phase: I) -> (I, bool) {
    match periodic_radius(phase, source.spec().pitch_mm, &mut |p| {
        authored_radius(source, p)
    }) {
        Ok(value) => (value, false),
        Err(_) => {
            // These private evaluators have only numeric Unsupported failures:
            // cancellation/resource publication is owned by the outer query.
            // Every exact normalized affine source band lies between its sealed
            // endpoint radii for ALL phases. This is an earned cavity bound,
            // not a nominal diameter or a claim that an envelope is material.
            (
                I::new(source.minimum_radius_mm(), source.maximum_radius_mm()),
                true,
            )
        }
    }
}
fn authored_radius(source: &AuthoredThread, phase: I) -> R<I> {
    let mut result = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for band in source.bands() {
        let p = I::new(
            phase.lo.max(band.phase_mm[0]),
            phase.hi.min(band.phase_mm[1]),
        );
        if p.lo > p.hi {
            continue;
        }
        let [a, b] = band.radius_mm;
        let u = (p - ip(band.phase_mm[0])) / (ip(band.phase_mm[1]) - ip(band.phase_mm[0]));
        let r = ip(a) + (ip(b) - ip(a)) * u;
        let r = finite(r)?;
        result.lo = result.lo.min(r.lo);
        result.hi = result.hi.max(r.hi);
    }
    finite(result)
}
fn periodic_radius(phase: I, pitch: f64, evaluate: &mut dyn FnMut(I) -> R<I>) -> R<I> {
    let quotient = finite(phase / ip(pitch))?;
    if quotient.width() >= 1. {
        return evaluate(I::new(0., pitch));
    }
    if quotient.lo.abs().max(quotient.hi.abs()) > 4096. {
        return Err(Error::Unsupported);
    }
    let mut result = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for k in quotient.lo.floor() as i64 - 1..=quotient.hi.floor() as i64 + 1 {
        let p = phase - ip(k as f64) * ip(pitch);
        let p = I::new(p.lo.max(0.), p.hi.min(pitch));
        if p.lo > p.hi {
            continue;
        }
        let r = evaluate(p)?;
        result.lo = result.lo.min(r.lo);
        result.hi = result.hi.max(r.hi);
    }
    finite(result)
}
#[cfg(test)]
mod tests;
