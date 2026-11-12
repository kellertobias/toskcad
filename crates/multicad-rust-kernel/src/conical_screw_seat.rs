//! Source-proven proposed seating of an actual retained screw-head meridian.
//!
//! This constructs a conservative local placement; it does not validate an
//! arbitrary existing pose. The positive bearing interval matches within the
//! reported radial gap. Clearance never spends that gap as intrusion tolerance.
use crate::head_profile::{Curve, Profile, Role, Segment};
use crate::tapered_thread::interval::I;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Seat {
    pub bore_radius_mm: f64,
    pub mouth_radius_mm: f64,
    pub depth_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Options {
    pub tolerance_mm: f64,
    /// Minimum guaranteed radial bearing width, not axial contact length.
    pub minimum_land_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Certificate {
    pub ideal_neck_z_mm: [f64; 2],
    pub safe_neck_z_mm: f64,
    pub certified_neck_z_mm: [f64; 2],
    pub placement_deviation_mm: f64,
    pub bearing_axial_mm: [f64; 2],
    pub bearing_radial_mm: [f64; 2],
    pub maximum_radial_gap_mm: f64,
    pub source_rounding_error_mm: f64,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Error {
    Invalid,
    NoCone,
    IncompatibleCone,
    InsufficientLand,
    Interference,
    Unresolved,
}
fn checked(v: I) -> Result<I, Error> {
    v.checked().map_err(|_| Error::Unresolved)
}
fn magnitude(v: I) -> Result<f64, Error> {
    Ok(checked(v)?.magnitude())
}
fn max(a: I, b: I) -> I {
    I::new(a.lo.max(b.lo), a.hi.max(b.hi))
}
fn trig(t: I) -> Result<(I, I), Error> {
    if !t.valid() || t.magnitude() > 2. {
        return Err(Error::Unresolved);
    }
    let square = -t.square();
    let (mut s, mut c) = (t, I::point(1.));
    let (mut st, mut ct) = (s, c);
    for n in 1..=16 {
        st = st * square / I::point(f64::from(2 * n * (2 * n + 1)));
        ct = ct * square / I::point(f64::from((2 * n - 1) * 2 * n));
        s = s + st;
        c = c + ct;
    }
    // Alternating Taylor tails decrease for |t|<=2. The next term bounds
    // every point remainder; interval evaluation already encloses arithmetic.
    let se = magnitude(st * square / I::point(34. * 35.))?;
    let ce = magnitude(ct * square / I::point(33. * 34.))?;
    Ok((checked(s + I::new(-se, se))?, checked(c + I::new(-ce, ce))?))
}
fn arc_bounds(segment: Segment) -> Result<[I; 2], Error> {
    let Curve::Arc { center, sweep } = segment.curve else {
        return Err(Error::Unresolved);
    };
    let [x, y] = std::array::from_fn(|i| I::point(segment.start[i]) - I::point(center[i]));
    let mut bounds = segment.start.map(I::point);
    for i in 0..256 {
        let a = I::point(sweep) * I::point(f64::from(i)) / I::point(256.);
        let b = I::point(sweep) * I::point(f64::from(i + 1)) / I::point(256.);
        let (s, c) = trig(I::new(a.lo.min(b.lo), a.hi.max(b.hi)))?;
        let p = [
            I::point(center[0]) + x * c - y * s,
            I::point(center[1]) + x * s + y * c,
        ];
        for (bound, p) in bounds.iter_mut().zip(p) {
            let p = checked(p)?;
            *bound = I::new(bound.lo.min(p.lo), bound.hi.max(p.hi));
        }
    }
    for (bound, end) in bounds.iter_mut().zip(segment.end) {
        *bound = I::new(bound.lo.min(end), bound.hi.max(end));
    }
    Ok(bounds)
}
fn allowed(seat: Seat, slope: I, z: I) -> I {
    max(
        I::point(seat.bore_radius_mm),
        I::point(seat.mouth_radius_mm) - slope * z,
    )
}
fn clear(radial: I, cavity: I) -> Result<(), Error> {
    let residual = checked(radial - cavity)?;
    if residual.hi > 0. {
        Err(if residual.lo > 0. {
            Error::Interference
        } else {
            Error::Unresolved
        })
    } else {
        Ok(())
    }
}
fn clear_segment(
    segment: Segment,
    shift: f64,
    reserve: f64,
    seat: Seat,
    slope: I,
) -> Result<(), Error> {
    match segment.curve {
        Curve::Arc { .. } => {
            let [radial, axial] = arc_bounds(segment)?;
            if (axial + I::point(shift) + I::point(reserve)).hi < 0. {
                return Ok(());
            }
            // The finite actual arc, not its full circle or a recipe radius,
            // must lie inside the bore. Otherwise this bounded proof refuses.
            if (radial + I::point(reserve)).hi > seat.bore_radius_mm {
                return Err(Error::Unresolved);
            }
        }
        Curve::Line => {
            let [r, z] = segment.start.map(I::point);
            if segment.start[1] == segment.end[1] {
                let far = z + I::point(shift) + I::point(reserve);
                if checked(far)?.hi < 0. {
                    return Ok(());
                }
                let radial = I::point(segment.start[0].max(segment.end[0])) + I::point(reserve);
                return clear(
                    radial,
                    allowed(seat, slope, z + I::point(shift) + I::point(reserve)),
                );
            }
            let dz = I::point(segment.end[1]) - z;
            if dz.lo <= 0. {
                return Err(Error::Unresolved);
            }
            let dr = I::point(segment.end[0]) - r;
            let low = segment.start[1].max((-I::point(shift) - I::point(reserve)).lo);
            let high = segment.end[1];
            if low > high {
                return Ok(());
            }
            let transition =
                checked(I::point(seat.depth_mm) - I::point(shift) - I::point(reserve))?;
            for at in [
                I::point(low),
                I::point(high),
                I::new(
                    transition.lo.clamp(low, high),
                    transition.hi.clamp(low, high),
                ),
            ] {
                let radial = r + dr * (at - z) / dz + I::point(reserve);
                let cavity = allowed(seat, slope, at + I::point(shift) + I::point(reserve));
                clear(radial, cavity)?;
            }
        }
    }
    Ok(())
}
fn supporting_cone(segments: &[Segment]) -> Result<(Segment, I, I), Error> {
    let mut cones = segments.iter().copied().filter(|s| s.role == Role::Cone);
    let cone = cones.next().ok_or(Error::NoCone)?;
    if cones.next().is_some() || cone.curve != Curve::Line {
        return Err(Error::IncompatibleCone);
    }
    let [r, z] = cone.start.map(I::point);
    let dz = checked(I::point(cone.end[1]) - z)?;
    if dz.lo <= 0. {
        return Err(Error::IncompatibleCone);
    }
    let k = checked((r - I::point(cone.end[0])) / dz)?;
    if k.lo <= 0. {
        return Err(Error::IncompatibleCone);
    }
    let intercept = checked(r + k * z)?;
    Ok((cone, k, intercept))
}
struct Proposal {
    segments: Vec<Segment>,
    cone: Segment,
    k: I,
    intercept: I,
    slope: I,
    ideal: I,
    safe: f64,
    reserve: f64,
}
fn propose(profile: &Profile, seat: Seat, options: Options) -> Result<Proposal, Error> {
    let reserve = profile.rounding_error_mm();
    if ![
        seat.bore_radius_mm,
        seat.mouth_radius_mm,
        seat.depth_mm,
        options.tolerance_mm,
        options.minimum_land_mm,
        reserve,
    ]
    .into_iter()
    .all(f64::is_finite)
        || seat.bore_radius_mm <= 0.
        || seat.mouth_radius_mm <= seat.bore_radius_mm
        || seat.depth_mm <= 0.
        || options.tolerance_mm <= 0.
        || options.minimum_land_mm <= 0.
        || reserve < 0.
    {
        return Err(Error::Invalid);
    }
    if reserve >= options.tolerance_mm {
        return Err(Error::Unresolved);
    }
    let segments: Vec<_> = profile.segments().collect();
    let (cone, k, intercept) = supporting_cone(&segments)?;
    let slope =
        (I::point(seat.mouth_radius_mm) - I::point(seat.bore_radius_mm)) / I::point(seat.depth_mm);
    let slope = checked(slope)?;
    let ideal = checked((I::point(seat.mouth_radius_mm) - intercept) / k)?;
    let variation = I::point(magnitude(k - slope)?) * I::point(seat.depth_mm);
    if variation.hi > options.tolerance_mm * 0.25 {
        return Err(Error::IncompatibleCone);
    }
    let source_margin = I::point(4. * reserve) * (I::point(1.) + slope);
    let placement_margin = I::point(8. * ideal.width()) * k;
    let margin = variation + source_margin + placement_margin;
    let margin = checked(margin + I::point(options.tolerance_mm) * I::point(0.25))?;
    let shift = checked(ideal - margin / k)?.lo;
    Ok(Proposal {
        segments,
        cone,
        k,
        intercept,
        slope,
        ideal,
        safe: shift,
        reserve,
    })
}

fn at(
    proposal: Proposal,
    seat: Seat,
    options: Options,
    neck: I,
    extra: f64,
) -> Result<Certificate, Error> {
    if !neck.valid() || !extra.is_finite() || extra < 0. {
        return Err(Error::Invalid);
    }
    let Proposal {
        segments,
        cone,
        k,
        intercept,
        slope,
        ideal,
        safe,
        reserve: source_reserve,
    } = proposal;
    let reserve = checked(I::point(source_reserve) + I::point(extra))?.hi;
    if reserve >= options.tolerance_mm {
        return Err(Error::Unresolved);
    }
    let z = I::point(cone.start[1]);
    let mouth_gap = I::point(seat.mouth_radius_mm) - intercept - k * neck;
    let departure = (k - slope) * I::new(0., seat.depth_mm);
    let source_gap = I::point(reserve) * (I::point(1.) + slope);
    let gap = checked(mouth_gap + departure + source_gap)?;
    if gap.hi > options.tolerance_mm {
        return Err(Error::Unresolved);
    }
    let radial_error = I::point(reserve) * (I::point(1.) + k);
    let low = checked(z + neck + I::point(reserve))?.hi.max(0.);
    let end = checked(I::point(cone.end[1]) + neck - I::point(reserve))?;
    let bore_end = (intercept - I::point(seat.bore_radius_mm) - radial_error) / k + neck;
    let high = end.lo.min(seat.depth_mm).min(checked(bore_end)?.lo);
    if high <= low {
        return Err(Error::InsufficientLand);
    }
    let inner = checked(intercept - k * (I::point(high) - neck) + radial_error)?.hi;
    let outer = checked(intercept - k * (I::point(low) - neck) - radial_error)?.lo;
    if checked(I::point(outer) - I::point(inner))?.lo < options.minimum_land_mm {
        return Err(Error::InsufficientLand);
    }
    for segment in segments {
        clear_segment(segment, neck.hi, reserve, seat, slope)?;
    }
    Ok(Certificate {
        ideal_neck_z_mm: [ideal.lo, ideal.hi],
        safe_neck_z_mm: safe,
        certified_neck_z_mm: [neck.lo, neck.hi],
        placement_deviation_mm: extra,
        bearing_axial_mm: [low, high],
        bearing_radial_mm: [inner, outer],
        maximum_radial_gap_mm: gap.hi,
        source_rounding_error_mm: source_reserve,
    })
}

/// Propose a seated source-origin station with a reserved radial gap. The
/// caller must physically apply this placement and prove its world frame.
pub fn certify(profile: &Profile, seat: Seat, options: Options) -> Result<Certificate, Error> {
    let proposal = propose(profile, seat, options)?;
    let neck = I::point(proposal.safe);
    at(proposal, seat, options, neck, 0.)
}
/// Certify an actual placed head, using its source-origin axial enclosure.
/// `placement_deviation_mm` must bound every head point's additional spatial
/// displacement from coaxial placement, including frame/axis/rounding errors.
/// All poses in the enclosure must clear; proximity to a proposal is not proof.
pub fn certify_at_neck(
    profile: &Profile,
    seat: Seat,
    options: Options,
    neck_z_mm: [f64; 2],
    placement_deviation_mm: f64,
) -> Result<Certificate, Error> {
    let proposal = propose(profile, seat, options)?;
    at(
        proposal,
        seat,
        options,
        I::new(neck_z_mm[0], neck_z_mm[1]),
        placement_deviation_mm,
    )
}

#[cfg(test)]
mod tests;
