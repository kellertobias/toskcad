//! Analytic whole-curve second-derivative and actual endpoint-transfer bounds.
use super::*;
fn curvature(body: &Chamfer, source: SourceCurve, t: I) -> Result<I, Error> {
    let zero = I::point(0.);
    let rate = body.rate_interval()?;
    // r, dr/dt, d²r/dt², dz/dt, d²z/dt² and dq/dt.
    let (r, r1, r2, z1, z2, q1) = match source {
        SourceCurve::Intersection { band, phase, .. }
        | SourceCurve::EndProfile { band, phase, .. } => {
            let dq = phase[1].range() - phase[0].range();
            let q = curve::interpolate(phase.map(PhaseEndpoint::range), t);
            let p = body.profile(band, q)?;
            let first = p.first * dq;
            let second = p.second * dq.square();
            let (z1, z2) = if active_kind(source) {
                (first / body.cone_slope()?, second / body.cone_slope()?)
            } else {
                (zero, zero)
            };
            (p.radius, first, second, z1, z2, dq)
        }
        SourceCurve::EntryArc { phase, .. } => (
            I::point(body.entry_radius_mm),
            zero,
            zero,
            zero,
            zero,
            phase[1].range() - phase[0].range(),
        ),
        SourceCurve::WallRail {
            band,
            phase,
            active,
            ..
        } => {
            let (r, g) = curve::endpoint_join(body, band, phase, active)?;
            (
                r,
                zero,
                zero,
                I::point(body.material_length_mm()) - g,
                zero,
                zero,
            )
        }
        SourceCurve::ConeRail { band, phase, .. } => {
            let (_, g) = curve::endpoint_join(body, band, phase, true)?;
            let k = body.cone_slope()?;
            (
                I::point(body.entry_radius_mm) + k * g * t,
                k * g,
                zero,
                g,
                zero,
                zero,
            )
        }
        SourceCurve::Circle { radius_mm, .. } => {
            return (I::point(2.) * I::point(radius_mm).abs() * (I::point(2.) * I::pi()).square())
                .checked();
        }
        SourceCurve::NeckSeam { .. } => return Ok(zero),
    };
    let theta1 = (z1 - q1) / rate;
    let theta2 = z2 / rate;
    // Rotation preserves the two radial components' Euclidean norm. Twice
    // their triangle-inequality bound also encloses the full x/y L1 norm.
    (I::point(2.)
        * (r2.abs()
            + r.abs() * theta1.square()
            + I::point(2.) * (r1 * theta1).abs()
            + (r * theta2).abs())
        + z2.abs())
    .checked()
}
fn active_kind(source: SourceCurve) -> bool {
    matches!(source, SourceCurve::Intersection { active: true, .. })
}
pub(super) fn error(
    body: &Chamfer,
    source: SourceCurve,
    fraction: [f64; 2],
    points: [Point; 2],
) -> Result<f64, Error> {
    if points.into_iter().flatten().any(|p| !p.is_finite()) {
        return Err(Error::Invalid);
    }
    // Reuse the range and source-domain validation before differentiating.
    curve::bounds(body, source, fraction)?;
    let mut endpoint = 0_f64;
    for (t, point) in fraction.into_iter().zip(points) {
        let value = curve::bounds(body, source, [t; 2])?;
        let error = (0..3)
            .fold(I::point(0.), |s, k| {
                s + I::point((value[k] - I::point(point[k])).magnitude())
            })
            .checked()?;
        endpoint = endpoint.max(error.hi);
    }
    let width = I::point(fraction[1]) - I::point(fraction[0]);
    let remainder =
        curvature(body, source, I::new(fraction[0], fraction[1]))? * width.square() / I::point(8.);
    Ok((remainder + I::point(endpoint)).checked()?.hi)
}
