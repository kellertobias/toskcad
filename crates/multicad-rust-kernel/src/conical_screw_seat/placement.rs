//! Actual finite placement admission for a retained axisymmetric head.
use super::*;
mod between;
pub use between::{BetweenFrameCertificate, Frame, certify_between_frames};
type Vector = [I; 3];
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Placement {
    pub origin_mm: [f64; 3],
    /// Actual columns, including the literal emitted normal, without repair.
    pub axes: [[f64; 3]; 3],
    /// P = origin + axes0*x + axes1*y + axes2*(profile_z + offset).
    pub axial_offset_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Axis {
    /// Receiving stock entry and cone mouth center.
    pub origin_mm: [f64; 3],
    /// Positive direction points into the receiving stock.
    pub direction: [f64; 3],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlacementCertificate {
    pub bearing: Certificate,
    pub neck_z_mm: [f64; 2],
    pub radial_offset_mm: f64,
    pub maximum_point_deviation_mm: f64,
    pub frame_gram_error: f64,
}
fn dot(a: Vector, b: Vector) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |sum, (a, b)| sum + a * b)
}
fn cross(a: Vector, b: Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn length(a: Vector) -> Result<I, Error> {
    checked(
        a.into_iter()
            .fold(I::point(0.), |sum, a| sum + a.square())
            .sqrt(),
    )
}
fn unit(a: Vector) -> Result<Vector, Error> {
    let norm = length(a)?;
    if norm.lo <= 0. {
        return Err(Error::Placement);
    }
    Ok(a.map(|a| a / norm))
}
fn frame(axes: [Vector; 3]) -> Result<f64, Error> {
    let mut gram = I::point(0.);
    for (i, &a) in axes.iter().enumerate() {
        for (j, &b) in axes.iter().enumerate() {
            let delta = dot(a, b) - I::point(if i == j { 1. } else { 0. });
            gram = gram + I::point(magnitude(delta)?);
        }
    }
    let determinant = checked(dot(cross(axes[0], axes[1]), axes[2]))?;
    if checked(gram)?.hi > 1e-6 || determinant.lo <= 0. {
        return Err(Error::Placement);
    }
    Ok(gram.hi)
}
fn extents(segments: &[Segment], reserve: f64) -> Result<(I, I), Error> {
    let mut radial = 0_f64;
    let mut axial = 0_f64;
    for &segment in segments {
        let bounds = match segment.curve {
            Curve::Line => std::array::from_fn(|i| {
                I::new(
                    segment.start[i].min(segment.end[i]),
                    segment.start[i].max(segment.end[i]),
                )
            }),
            Curve::Arc { .. } => arc_bounds(segment)?,
        };
        radial = radial.max(magnitude(bounds[0])?);
        axial = axial.max(magnitude(bounds[1])?);
    }
    Ok((
        I::point(radial) + I::point(reserve),
        I::point(axial) + I::point(reserve),
    ))
}
#[derive(Clone, Copy)]
struct IntervalPlacement {
    origin: Vector,
    axes: [Vector; 3],
    offset: I,
}
fn interval(placement: Placement) -> Result<IntervalPlacement, Error> {
    if !placement
        .origin_mm
        .into_iter()
        .chain(placement.axes.into_iter().flatten())
        .chain([placement.axial_offset_mm])
        .all(f64::is_finite)
    {
        return Err(Error::Invalid);
    }
    Ok(IntervalPlacement {
        origin: placement.origin_mm.map(I::point),
        axes: placement.axes.map(|a| a.map(I::point)),
        offset: I::point(placement.axial_offset_mm),
    })
}

fn floating_bound(placement: IntervalPlacement, radius: I, height: I) -> Result<I, Error> {
    // Eight-operation gamma includes source-coordinate offset plus ordinary
    // affine evaluation (including FMA/reassociation), with underflow reserve.
    let unit_roundoff = I::point(f64::EPSILON * 0.5);
    let gamma = I::point(8.) * unit_roundoff / (I::point(1.) - I::point(8.) * unit_roundoff);
    let z = height + I::point(magnitude(placement.offset)?);
    let mut floating = I::point(0.);
    for i in 0..3 {
        let bound = I::point(magnitude(placement.origin[i])?)
            + I::point(magnitude(placement.axes[0][i])?) * radius
            + I::point(magnitude(placement.axes[1][i])?) * radius
            + I::point(magnitude(placement.axes[2][i])?) * z;
        floating = floating + gamma * bound + I::point(8. * f64::from_bits(1));
    }
    checked(floating)
}
fn admit(
    proposal: Proposal,
    seat: Seat,
    options: Options,
    placement: IntervalPlacement,
    axis_origin: Vector,
    direction: Vector,
    extra: I,
) -> Result<PlacementCertificate, Error> {
    let axes = placement.axes;
    let gram = frame(axes)?;
    if checked(dot(axes[2], direction))?.lo <= 0. {
        return Err(Error::Placement);
    }
    let radial_axial = dot(axes[0], direction);
    let u = unit(std::array::from_fn(|i| {
        axes[0][i] - direction[i] * radial_axial
    }))?;
    let reference = [u, cross(direction, u), direction];
    let (radius, height) = extents(&proposal.segments, proposal.reserve)?;
    let mut deformation = I::point(0.);
    for (i, (&actual, &ideal)) in axes.iter().zip(&reference).enumerate() {
        let error = length(std::array::from_fn(|j| actual[j] - ideal[j]))?;
        deformation = deformation + error * if i == 2 { height } else { radius };
    }
    // The outward world evaluation contains the exact affine neck and its
    // rounded multiply/add variants. Subtract stock origin before projection.
    let relative = std::array::from_fn(|i| {
        placement.origin[i] + axes[2][i] * placement.offset - axis_origin[i]
    });
    let neck = checked(dot(relative, direction))?;
    let radial_offset = length(std::array::from_fn(|i| relative[i] - direction[i] * neck))?;
    let floating = floating_bound(placement, radius, height)?;
    let deviation = checked(radial_offset + deformation + floating + extra)?.hi;
    if deviation > options.tolerance_mm {
        return Err(Error::Placement);
    }
    let bearing = at(proposal, seat, options, neck, deviation)?;
    Ok(PlacementCertificate {
        bearing,
        neck_z_mm: [neck.lo, neck.hi],
        radial_offset_mm: radial_offset.hi,
        maximum_point_deviation_mm: deviation,
        frame_gram_error: gram,
    })
}

/// Prove the actual supplied head frame against an actual axisymmetric cavity.
/// Axis must describe the cavity itself; this does not infer receiving shear.
/// No retained input geometry is normalized or translated by this admission.
pub fn certify_placement(
    profile: &Profile,
    seat: Seat,
    options: Options,
    placement: Placement,
    axis: Axis,
) -> Result<PlacementCertificate, Error> {
    let proposal = propose(profile, seat, options)?;
    let placement = interval(placement)?;
    if !axis
        .origin_mm
        .into_iter()
        .chain(axis.direction)
        .all(f64::is_finite)
    {
        return Err(Error::Invalid);
    }
    let direction_scale = axis.direction.into_iter().map(f64::abs).fold(0., f64::max);
    if direction_scale == 0. {
        return Err(Error::Placement);
    }
    let direction = unit(
        axis.direction
            .map(|v| I::point(v) / I::point(direction_scale)),
    )?;
    admit(
        proposal,
        seat,
        options,
        placement,
        axis.origin_mm.map(I::point),
        direction,
        I::point(0.),
    )
}

#[cfg(test)]
mod tests;
