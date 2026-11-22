use super::*;
const MAX_MM: f64 = 1e7;
fn cancel(cancelled: &dyn Fn() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
struct Angle {
    radians: I,
    sine: I,
    cosine: I,
}
fn angle(degrees: f64, divisor: f64) -> Angle {
    let pi = I::new(PI.next_down(), PI.next_up());
    let radians = (I::point(degrees) / I::point(divisor)) * pi / I::point(180.);
    // Callers prove a strictly acute authored angle. Intersecting its trig
    // enclosure with [0,1] uses that analytic domain, not a geometry tolerance.
    let positive = |v: I| I::new(v.lo.max(0.), v.hi.min(1.));
    Angle {
        radians,
        sine: positive(radians.sin()),
        cosine: positive(radians.cos()),
    }
}
fn valid(value: I) -> Result<I, Error> {
    if !value.valid() {
        Err(Error::Unresolved)
    } else {
        Ok(value)
    }
}
fn bounded(value: I) -> Result<I, Error> {
    let value = valid(value)?;
    if value.magnitude() > MAX_MM {
        Err(Error::Limit)
    } else {
        Ok(value)
    }
}
fn strictly_less(a: I, b: I) -> Result<(), Error> {
    if a.lo >= b.hi {
        Err(Error::Invalid)
    } else if a.hi >= b.lo {
        Err(Error::Unresolved)
    } else {
        Ok(())
    }
}
fn arc(center: [I; 2], radius: f64, sweep: I) -> Curve {
    Curve::CircularArc {
        center_mm_bounds: point2(center),
        radius_mm: radius,
        sweep_rad_bounds: bounds(sweep),
    }
}
fn section(role: SectionRole, a: [I; 2], b: [I; 2], curve: Curve) -> Option<SectionSegment> {
    Some(SectionSegment {
        role,
        start_mm_bounds: point2(a),
        end_mm_bounds: point2(b),
        curve,
    })
}
pub(super) fn construct(
    spec: Spec,
    domain: Domain,
    cancelled: &dyn Fn() -> bool,
) -> Result<Cutter, Error> {
    cancel(cancelled)?;
    let acute = |v: f64| v > 0. && v < 90.;
    let lengths = [
        spec.mouth_diameter_mm,
        spec.mouth_radius_mm,
        spec.b_width_mm,
        spec.g_diameter_mm,
    ];
    let angles = [
        spec.tilt_angle_deg,
        spec.outer_flank_angle_deg,
        spec.floor_angle_deg,
        spec.outer_corner_angle_deg,
    ];
    let domains = [
        domain.corner_start_depth_mm,
        domain.corner_end_depth_mm,
        domain.corner_reach_mm,
    ];
    if !lengths
        .into_iter()
        .chain(angles)
        .chain(domains)
        .all(f64::is_finite)
        || lengths.into_iter().any(|v| v <= 0.)
        || spec.b_width_mm >= spec.g_diameter_mm
        || spec.g_diameter_mm >= spec.mouth_diameter_mm
        || 2. * spec.mouth_radius_mm >= spec.mouth_diameter_mm
        || !acute(spec.tilt_angle_deg)
        || !acute(spec.outer_flank_angle_deg)
        || !acute(spec.floor_angle_deg)
        || spec.outer_corner_angle_deg <= 0.
        || spec.outer_corner_angle_deg >= 180.
        || domain.corner_start_depth_mm >= domain.corner_end_depth_mm
        || domain.corner_reach_mm <= 0.
    {
        return Err(Error::Invalid);
    }
    let inner_length = match spec.inner {
        Inner::Rounded { radius_mm } => radius_mm,
        Inner::ShallowV {
            angle_deg,
            e_width_mm,
        } => {
            if !angle_deg.is_finite()
                || angle_deg <= 90.
                || angle_deg >= 180.
                || e_width_mm >= spec.b_width_mm
            {
                return Err(Error::Invalid);
            }
            e_width_mm
        }
    };
    if !inner_length.is_finite() || inner_length <= 0. {
        return Err(Error::Invalid);
    }
    if lengths
        .into_iter()
        .chain(domains)
        .chain([inner_length])
        .any(|v| v.abs() > MAX_MM)
    {
        return Err(Error::Limit);
    }
    cancel(cancelled)?;
    let [a, c, beta] = [
        spec.outer_flank_angle_deg,
        spec.floor_angle_deg,
        spec.tilt_angle_deg,
    ]
    .map(|degrees| angle(degrees, 1.));
    let theta = angle(spec.outer_corner_angle_deg, 2.);
    let m = I::point(spec.mouth_diameter_mm);
    let g = I::point(spec.g_diameter_mm);
    let b = I::point(spec.b_width_mm);
    let r = I::point(spec.mouth_radius_mm);
    let two = I::point(2.);
    // cos/(1+sin) and cos²/(1+sin) retain tiny positive joins without
    // subtracting two nearly equal evaluations of the tangent construction.
    let entry = bounded(m / two + r * a.cosine / (I::point(1.) + a.sine))?;
    let join = bounded(r * a.cosine.square() / (I::point(1.) + a.sine))?;
    let flank = bounded((m - g) * a.cosine / (two * a.sine))?;
    strictly_less(join, flank)?;
    let tan_c = valid(c.sine / c.cosine)?;
    let tip = bounded(flank + g / two * tan_c)?;
    let b_depth = bounded(flank + (g - b) / two * tan_c)?;
    let zero = I::point(0.);
    let mouth = [entry, zero];
    let tangent = [entry - r * a.cosine, r - r * a.sine];
    let floor_join = [g / two, flank];
    let pi = I::new(PI.next_down(), PI.next_up());
    let segments = [
        MeridianSegment {
            role: EnvelopeRole::EntryRoundover,
            start_mm_bounds: point2(mouth),
            end_mm_bounds: point2(tangent),
            curve: arc([entry, r], spec.mouth_radius_mm, a.radians - pi / two),
        },
        MeridianSegment {
            role: EnvelopeRole::OuterFlank,
            start_mm_bounds: point2(tangent),
            end_mm_bounds: point2(floor_join),
            curve: Curve::Line,
        },
        MeridianSegment {
            role: EnvelopeRole::Floor,
            start_mm_bounds: point2(floor_join),
            end_mm_bounds: point2([zero, tip]),
            curve: Curve::Line,
        },
    ];
    let envelope = Envelope {
        segments,
        entry,
        join,
        flank,
        tip,
        sin_a: a.sine,
        cos_a: a.cosine,
        sin_c: c.sine,
        cos_c: c.cosine,
    };
    cancel(cancelled)?;
    let (inner, join_s, join_y) = match spec.inner {
        Inner::Rounded { radius_mm } => {
            let radius = I::point(radius_mm);
            let y = radius * theta.cosine;
            strictly_less(two * y, b)?;
            (
                SourceInner::Rounded { radius_mm },
                radius * theta.cosine.square() / (I::point(1.) + theta.sine),
                y,
            )
        }
        Inner::ShallowV {
            angle_deg,
            e_width_mm,
        } => {
            let half = angle(angle_deg, 2.);
            let y = I::point(e_width_mm) / two;
            (
                SourceInner::ShallowV {
                    sine: half.sine,
                    cosine: half.cosine,
                },
                y * half.cosine / half.sine,
                y,
            )
        }
    };
    let join_s = valid(join_s)?;
    let join_y = valid(join_y)?;
    let far_s = valid(join_s + I::point(domain.corner_reach_mm) * theta.cosine)?;
    let far_y = valid(join_y + I::point(domain.corner_reach_mm) * theta.sine)?;
    let upper = [join_s, join_y];
    let lower = [join_s, -join_y];
    let far_upper = [far_s, far_y];
    let far_lower = [far_s, -far_y];
    let mut section_segments = [None; 5];
    match spec.inner {
        Inner::Rounded { radius_mm } => {
            section_segments[..4].copy_from_slice(&[
                section(
                    SectionRole::InnerRound,
                    lower,
                    upper,
                    arc(
                        [I::point(radius_mm), zero],
                        radius_mm,
                        two * theta.radians - pi,
                    ),
                ),
                section(SectionRole::OuterUpper, upper, far_upper, Curve::Line),
                section(SectionRole::Closure, far_upper, far_lower, Curve::Line),
                section(SectionRole::OuterLower, far_lower, lower, Curve::Line),
            ]);
        }
        Inner::ShallowV { .. } => {
            section_segments = [
                section(SectionRole::InnerUpper, [zero; 2], upper, Curve::Line),
                section(SectionRole::OuterUpper, upper, far_upper, Curve::Line),
                section(SectionRole::Closure, far_upper, far_lower, Curve::Line),
                section(SectionRole::OuterLower, far_lower, lower, Curve::Line),
                section(SectionRole::InnerLower, lower, [zero; 2], Curve::Line),
            ];
        }
    }
    cancel(cancelled)?;
    let corner = Corner {
        quarter_turn: 0,
        domain,
        section: section_segments,
        sin_beta: beta.sine,
        cos_beta: beta.cosine,
        sin_theta: theta.sine,
        cos_theta: theta.cosine,
        b_depth,
        join_s,
        join_y,
        far_s,
        b_half: b / two,
        inner,
    };
    let mut corners = std::array::from_fn(|_| corner.clone());
    for (i, corner) in corners.iter_mut().enumerate() {
        cancel(cancelled)?;
        corner.quarter_turn = i as u8;
    }
    Ok(Cutter {
        spec,
        domain,
        envelope,
        corners,
    })
}
