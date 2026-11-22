//! Whole-period clearance of the NEW continuous authored metric source inside
//! the ORIGINAL retained nut. The auxiliary similarity only bounds the actual
//! affine map; it never changes either frame, profile, phase or finite domain.
use super::super::insert::mating::phase_error;
use super::*;
use crate::{
    authored_metric_thread::{AuthoredMetricThread, Knot},
    helical_thread::Thread,
    hex_nut::Nut,
};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HangerMatingRequirements {
    pub minimum_engagement_mm: f64,
    pub minimum_radial_clearance_mm: f64,
    pub linear_tolerance_mm: f64,
    pub maximum_cells: usize,
}
impl Default for HangerMatingRequirements {
    fn default() -> Self {
        Self {
            minimum_engagement_mm: 1.,
            minimum_radial_clearance_mm: 1e-7,
            linear_tolerance_mm: 1e-7,
            maximum_cells: 4096,
        }
    }
}
#[derive(Clone, Debug)]
pub struct HangerMating {
    /// Actual normal separation of the nut cap planes. Full nut containment
    /// makes this a conservative geometric thread-engagement measure.
    pub engagement_mm: [f64; 2],
    pub minimum_radial_clearance_lower_mm: f64,
    pub metric_tip_clearance_mm: [f64; 2],
    pub plain_shoulder_clearance_mm: [f64; 2],
    pub relative_radial_scale: [f64; 2],
    pub finite_radial_map_error_mm: f64,
    pub finite_axial_map_error_mm: f64,
    pub relative_angle_error_rad: f64,
    pub source_phase_offset_mm: [f64; 2],
    pub nut_indexed_phase_error_mm: f64,
    pub nut_evaluator_phase_error_mm: f64,
    pub phase_cells_used: usize,
    pub proved_phase_cells: usize,
}
#[derive(Clone, Copy)]
struct Alignment {
    scale: I,
    angle: f64,
    angle_error: f64,
    radial: f64,
    axial: f64,
    origin_z: I,
}
fn tau() -> I {
    ip(2.)
        * I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        )
}
fn absq(x: &Q) -> R<I> {
    finite(exact_interval(&x.abs())?)
}
/// atan2 proposes a principal angle for a finite nonzero enclosed ray. Its
/// bounded output is not an accuracy assertion: true interval trig earns that
/// separately. A broad ray enclosure may honestly fail this certificate.
fn ray_angle(ci: I, si: I, scale: I) -> R<(f64, f64)> {
    let ci = finite(ci)?;
    let si = finite(si)?;
    let scale = positive(scale)?;
    let angle = si.midpoint().portable_atan2(ci.midpoint());
    let ray = [ci / scale, si / scale];
    let unit = [ip(angle).cos(), ip(angle).sin()];
    let chord = finite(ip((ray[0] - unit[0]).magnitude()) + ip((ray[1] - unit[1]).magnitude()))?;
    if chord.hi > 1. {
        return Err(Error::Unsupported);
    }
    Ok((angle, (ip(2.) * chord).hi))
}
fn phase_midpoint(a: f64, b: f64) -> R<f64> {
    let mid = a + (b - a) * 0.5;
    if mid <= a || mid >= b {
        Err(Error::Unsupported)
    } else {
        Ok(mid)
    }
}
fn align(source: &Hanger, from: Frame, to: Frame, interrupt: &dyn Fn() -> bool) -> R<Alignment> {
    let (mut o, a) = affine(to, from, interrupt)?;
    let shift = bounded(-q(source.spec().plain_length_mm)? - q(source.metric().full_length_mm())?)?;
    for i in 0..3 {
        o[i] = bounded(&o[i] + &a[2][i] * &shift)?;
    }
    let (c, s) = (&a[0][0], &a[0][1]);
    let ci = exact_interval(c)?;
    let si = exact_interval(s)?;
    let scale = positive(ci.square() + si.square())?.sqrt();
    let scale = positive(scale)?;
    let (angle, angle_error) = ray_angle(ci, si, scale)?;
    let zero = Q::zero();
    let one = q(1.)?;
    let ideal = [
        [c.clone(), s.clone(), zero.clone()],
        [-s.clone(), c.clone(), zero.clone()],
        [zero.clone(), zero, one],
    ];
    let residual: [[Q; 3]; 3] =
        std::array::from_fn(|j| std::array::from_fn(|i| &a[j][i] - &ideal[j][i]));
    let r = ip(source.metric().maximum_radius_mm());
    let l = ip(source.metric().full_length_mm());
    let mut error = [ip(0.); 3];
    for i in 0..3 {
        error[i] =
            (absq(&residual[0][i])? + absq(&residual[1][i])?) * r + absq(&residual[2][i])? * l;
        if i < 2 {
            error[i] = error[i] + absq(&o[i])?;
        }
    }
    let radial = finite(error[0] + error[1])?.hi;
    let axial = finite(error[2])?.hi;
    Ok(Alignment {
        scale,
        angle,
        angle_error,
        radial,
        axial,
        origin_z: exact_interval(&o[2])?,
    })
}
fn metric_error(e: crate::authored_metric_thread::Error) -> Error {
    match e {
        crate::authored_metric_thread::Error::Invalid => Error::Invalid,
        crate::authored_metric_thread::Error::Cancelled => Error::Cancelled,
        crate::authored_metric_thread::Error::Limit => Error::Limit,
        crate::authored_metric_thread::Error::Unresolved => Error::Unsupported,
    }
}
fn continuous_radius(source: &AuthoredMetricThread, phase: I) -> R<I> {
    let knots = [
        Knot::Begin,
        Knot::CrestEnd,
        Knot::UpperRoot,
        Knot::LowerRoot,
        Knot::CrestStart,
        Knot::End,
    ]
    .map(|k| source.knot(k));
    let mut value = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for band in 0..5 {
        let a = I::new(knots[band].phase_mm[0], knots[band].phase_mm[1]);
        let b = I::new(knots[band + 1].phase_mm[0], knots[band + 1].phase_mm[1]);
        let p = I::new(phase.lo.max(a.lo), phase.hi.min(b.hi));
        if p.lo > p.hi {
            continue;
        }
        let u = (p - a) / (b - a);
        let bound = source
            .profile_bounds(band, [u.lo.max(0.), u.hi.min(1.)])
            .map_err(metric_error)?;
        value.lo = value.lo.min(bound.radius_mm[0]);
        value.hi = value.hi.max(bound.radius_mm[1]);
    }
    finite(value)
}
/// Actual one-sided coefficient/support evaluation on ALL potentially selected
/// finite trims. Enlarging the trims also retains tiny global-station extensions
/// instead of silently projecting an indexed flank onto a canonical endpoint.
fn retained_radius(source: &Thread, phase: I, pad: f64) -> R<I> {
    let breaks = source.profile_breaks();
    let spec = source.spec();
    let mut value = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for band in 0..5 {
        let lo = (ip(breaks[band]) - ip(pad)).lo;
        let hi = (ip(breaks[band + 1]) + ip(pad)).hi;
        let p = I::new(phase.lo.max(lo), phase.hi.min(hi));
        if p.lo > p.hi {
            continue;
        }
        let r = match band {
            1 | 3 => {
                ip(spec.major_radius_mm)
                    + ip(source.slope_at_phase(0., band))
                        * (p - ip(breaks[if band == 1 { 1 } else { 4 }]))
            }
            2 => {
                if let Some(root) = source.root_circle() {
                    let x = p - ip(root.phase);
                    let a = ip(root.radius) - x;
                    let b = ip(root.radius) + x;
                    if a.lo < 0. || b.lo < 0. {
                        return Err(Error::Unsupported);
                    }
                    ip(root.center_radius) - (a * b).sqrt()
                } else {
                    ip(spec.major_radius_mm) - ip(spec.depth_mm)
                }
            }
            _ => ip(spec.major_radius_mm),
        };
        let r = finite(r)?;
        value.lo = value.lo.min(r.lo);
        value.hi = value.hi.max(r.hi);
    }
    finite(value)
}
fn retained_periodic(source: &Thread, phase: I, pad: f64) -> R<I> {
    let pitch = source.spec().pitch_mm;
    let quotient = finite(phase / ip(pitch))?;
    if quotient.width() >= 1. {
        return retained_radius(source, I::new(-pad, (ip(pitch) + ip(pad)).hi), pad);
    }
    if quotient.magnitude() > 4096. {
        return Err(Error::Unsupported);
    }
    let mut result = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for k in quotient.lo.floor() as i64 - 1..=quotient.hi.floor() as i64 + 1 {
        let local = phase - ip(k as f64) * ip(pitch);
        let local = I::new(local.lo.max(-pad), local.hi.min((ip(pitch) + ip(pad)).hi));
        if local.lo > local.hi {
            continue;
        }
        let r = retained_radius(source, local, pad)?;
        result.lo = result.lo.min(r.lo);
        result.hi = result.hi.max(r.hi);
    }
    finite(result)
}
/// Exact support of a typed finite Nut corner after `scale` and `affine`.
/// All binary64 fields share D=2^1074. Axes<=2 give integer determinant
/// <=3231 bits; relative-column numerators<=3231 and origin numerators<=3254.
/// Finite Nut corners and the symbolic shoulder shift are <=2e7. On common
/// denominator det*D, the complete corner has <8704 combined numerator and
/// denominator bits, below bounded()'s 16384 cap. Rational normalization cannot
/// increase that common representation. The affine constructor's checks and
/// the caller's finite coordinate conversions remain authoritative; repeating
/// a generic overflow check on this earned domain would be unreachable.
fn axial_corner(o: &Q, a: [&Q; 3], xyz: [Q; 3], shift: &Q) -> Q {
    o + a[0] * &xyz[0] + a[1] * &xyz[1] + a[2] * &xyz[2] - shift
}
/// Complete directed phase uncertainty, including original finite evaluation
/// and the possible true-angle period wrap. No angle normalization changes law.
fn phase_offset(
    map: &Alignment,
    he: I,
    hn: I,
    lead: I,
    left: bool,
    allowance: phase_error::Allowance,
) -> R<I> {
    let wrap = (hn * tau() - ip(if left { -1. } else { 1. }) * lead).magnitude();
    let error = finite(
        ip(map.axial)
            + ip(hn.magnitude()) * ip(map.angle_error)
            + ip(wrap)
            + ip(allowance.total_mm),
    )?
    .hi;
    finite(
        map.origin_z - hn * ip(map.angle)
            + (he - hn) * I::new(0., tau().hi)
            + I::new(-error, error),
    )
}
/// Radial correspondence on checked actual source intervals, not minimum
/// Euclidean flank distance. Thread admission gives slopes<58, root-center
/// radius<1e9 and dimensions<=1e7. The finite original phase reserve is<.001;
/// retained_radius checks the whole extended circle domain and all finite
/// bounds. Authored radii are positive. Checked near-isometric frames give
/// scale<2; the independently admitted residual is<.002. Thus all products
/// below stay below 1e12, including directed rounding, and are finite. Positive
/// radius/scale ordering gives low<=high. No generic unknown interval reaches
/// this source-owned handoff; source/frame/phase failures still propagate.
fn physical_gap(inner: I, external: I, map: &Alignment, metric: &Metric) -> I {
    let low = ip(inner.lo) - ip(map.scale.hi) * ip(external.hi) - ip(map.radial);
    let high = ip(inner.hi) - ip(map.scale.lo) * ip(external.lo) + ip(map.radial);
    I::new(low.lo, high.hi) * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi)
}
fn finite_nut_span(
    source: &Hanger,
    frame: Frame,
    nut: &Nut,
    nut_frame: Frame,
    interrupt: &dyn Fn() -> bool,
) -> R<(Q, Q)> {
    let (o, a) = affine(frame, nut_frame, interrupt)?;
    let shift = bounded(-q(source.spec().plain_length_mm)? - q(source.metric().full_length_mm())?)?;
    let mut bounds: Option<(Q, Q)> = None;
    for xy in crate::headed_screw::hex_vertices(nut.across_flats_mm()) {
        for z in [0., nut.thread().spec().length_mm] {
            checkpoint(interrupt)?;
            let xyz = [q(xy[0])?, q(xy[1])?, q(z)?];
            let p = axial_corner(&o[2], [&a[0][2], &a[1][2], &a[2][2]], xyz, &shift);
            bounds = Some(bounds.map_or_else(
                || (p.clone(), p.clone()),
                |(lo, hi)| (lo.min(p.clone()), hi.max(p.clone())),
            ));
        }
    }
    // Six actual source corners times two noncoincident finite axial stations.
    let (lo, hi) = bounds.unwrap();
    let length = q(source.metric().full_length_mm())?;
    if lo < Q::zero() || hi >= length {
        return Err(Error::Unsupported);
    }
    Ok((lo, bounded(length - hi)?))
}
fn engagement(nut: &Nut, frame: Frame, minimum: f64) -> R<I> {
    let a = [
        vector(frame.axes[0])?,
        vector(frame.axes[1])?,
        vector(frame.axes[2])?,
    ];
    let normal = cross(&a[0], &a[1])?;
    let denominator = dot(&normal, &normal)?;
    let numerator = bounded(q(nut.thread().spec().length_mm)? * dot(&a[2], &normal)?)?;
    if denominator <= Q::zero()
        || numerator <= Q::zero()
        || bounded(numerator.pow(2))? < bounded(q(minimum)?.pow(2) * &denominator)?
    {
        return Err(Error::Unsupported);
    }
    let mut value =
        finite(exact_interval(&numerator)? / positive(exact_interval(&denominator)?)?.sqrt())?;
    value.lo = value.lo.max(minimum);
    Ok(value)
}
/// All accepted cells cover the complete physical phase period. Corner samples
/// never certify clearance, and no legacy transfer certificate substitutes for
/// either actual profile. Other components/receivers require their own checks.
pub fn certify_hanger_nut(
    source: &Hanger,
    frame: Frame,
    nut: &Nut,
    nut_frame: Frame,
    req: HangerMatingRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<HangerMating> {
    checkpoint(interrupt)?;
    if [req.minimum_engagement_mm, req.minimum_radial_clearance_mm]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0. || x > 1e7)
        || !(1e-7..=1e-3).contains(&req.linear_tolerance_mm)
        || req.maximum_cells == 0
        || req.maximum_cells > 4096
    {
        return Err(Error::Invalid);
    }
    let e = source.metric();
    let i = nut.thread();
    let es = e.spec();
    let is = i.spec();
    if es.pitch_mm != is.pitch_mm || es.starts != is.starts || es.left_hand != is.left_hand {
        return Err(Error::Unsupported);
    }
    let source_metric = scale(frame)?;
    let nut_metric = scale(nut_frame)?;
    let (tip, shoulder) = finite_nut_span(source, frame, nut, nut_frame, interrupt)?;
    let engaged = engagement(nut, nut_frame, req.minimum_engagement_mm)?;
    let map = align(source, frame, nut_frame, interrupt)?;
    if ((ip(map.radial) + ip(map.axial)) * nut_metric.maximum_scale).hi > req.linear_tolerance_mm {
        return Err(Error::Unsupported);
    }
    let allowance = phase_error::measure(i, interrupt)?;
    let bore = positive(retained_periodic(
        i,
        I::new(0., is.pitch_mm),
        allowance.total_mm,
    )?)?;
    let core = positive(ip(e.minimum_radius_mm()[0]) * ip(0.5))?.lo;
    if (map.scale * ip(core) + ip(map.radial)).hi >= bore.lo {
        return Err(Error::Unsupported);
    }
    let radial = positive(map.scale * ip(core))?;
    if ip(map.radial).hi > (radial * ip(0.5)).lo {
        return Err(Error::Unsupported);
    }
    let turn_error = finite(ip(2.) * ip(map.radial) / radial)?.hi;
    let angle_error = finite(ip(map.angle_error) + ip(turn_error))?.hi;
    // Principal |angle|<=pi, certified ray error<=2 and turn error<=1
    // give total<pi+3<true2pi: at most one actual angular seam wraps.

    let lead = ip(es.pitch_mm) * ip(f64::from(es.starts));
    let he = lead / tau() * ip(if es.left_hand { -1. } else { 1. });
    let hn = ip(i.helical_rate());
    let phase_map = Alignment { angle_error, ..map };
    let offset = phase_offset(&phase_map, he, hn, lead, is.left_hand, allowance)?;
    let mut pending = Vec::with_capacity(req.maximum_cells + 1);
    pending.push((0., es.pitch_mm));
    let mut used = 0;
    let mut proved = 0;
    let mut clearance = f64::INFINITY;
    while let Some((a, b)) = pending.pop() {
        checkpoint(interrupt)?;
        if used == req.maximum_cells {
            return Err(Error::Limit);
        }
        used += 1;
        let phase = I::new(a, b);
        let external = continuous_radius(e, phase)?;
        let inner = retained_periodic(i, phase + offset, allowance.total_mm)?;
        let physical = physical_gap(inner, external, &map, &nut_metric);
        if physical.lo >= req.minimum_radial_clearance_mm {
            clearance = clearance.min(physical.lo);
            proved += 1;
            continue;
        }
        if physical.hi < req.minimum_radial_clearance_mm {
            return Err(Error::Unsupported);
        }
        let mid = phase_midpoint(a, b)?;
        if used == req.maximum_cells {
            return Err(Error::Limit);
        }
        pending.push((a, mid));
        pending.push((mid, b));
    }
    let tip = exact_interval(&tip)?
        * I::new(
            source_metric.minimum_scale.lo,
            source_metric.maximum_scale.hi,
        );
    let shoulder = exact_interval(&shoulder)?
        * I::new(
            source_metric.minimum_scale.lo,
            source_metric.maximum_scale.hi,
        );
    checkpoint(interrupt)?;
    Ok(HangerMating {
        engagement_mm: [engaged.lo, engaged.hi],
        minimum_radial_clearance_lower_mm: clearance,
        metric_tip_clearance_mm: [tip.lo.max(0.), tip.hi],
        plain_shoulder_clearance_mm: [shoulder.lo.max(0.), shoulder.hi],
        relative_radial_scale: [map.scale.lo, map.scale.hi],
        finite_radial_map_error_mm: map.radial,
        finite_axial_map_error_mm: map.axial,
        relative_angle_error_rad: angle_error,
        source_phase_offset_mm: [offset.lo, offset.hi],
        nut_indexed_phase_error_mm: allowance.indexed_mm,
        nut_evaluator_phase_error_mm: allowance.evaluator_mm,
        phase_cells_used: used,
        proved_phase_cells: proved,
    })
}
