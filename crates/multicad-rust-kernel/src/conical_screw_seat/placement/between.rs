//! Actual head and receiving affine charts, without midpoint or axis repair.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    pub origin_mm: [f64; 3],
    pub axes: [[f64; 3]; 3],
    /// World P = origin + axes * (offset + receiving_chart_P).
    pub offset_mm: [f64; 3],
    /// Total local per-point cavity uncertainty from retained coordinate/radius
    /// conversion, including all independent effects. It is spent as strict
    /// clearance erosion and included in the maximum gap, never waived.
    pub source_error_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BetweenFrameCertificate {
    /// Coordinates and error budgets in the actual receiving source chart.
    pub chart: PlacementCertificate,
    pub maximum_physical_gap_mm: f64,
    pub minimum_physical_land_mm: f64,
    pub receiver_operator_norm_upper: f64,
    pub receiver_inverse_norm_upper: f64,
}
// Exact finite add residual; retain its sign so large cancellations do not
// turn several known rounding corrections into a fictitious positional tube.
fn add(sum: &mut f64, correction: &mut I, term: f64) -> Result<(), Error> {
    let next = *sum + term;
    let bridge = next - *sum;
    let residual = (*sum - (next - bridge)) + (term - bridge);
    checked(I::point(next))?;
    *correction = checked(*correction + I::point(residual))?;
    *sum = next;
    Ok(())
}
fn affine_relative(head: Placement, receiver: Frame) -> Result<Vector, Error> {
    let mut result = [I::point(0.); 3];
    for (i, out) in result.iter_mut().enumerate() {
        let mut sum = 0.;
        let mut correction = I::point(0.);
        add(&mut sum, &mut correction, head.origin_mm[i])?;
        add(&mut sum, &mut correction, -receiver.origin_mm[i])?;
        for j in 0..3 {
            let a = receiver.axes[j][i];
            let b = receiver.offset_mm[j];
            let product = a * b;
            checked(I::point(product))?;
            let residual = a.mul_add(b, -product);
            // FMA TwoProduct residual is exact away from underflow; one least
            // subnormal on each side also encloses rounded tiny residuals.
            let tiny = f64::from_bits(1);
            let error = checked(I::point(residual) + I::new(-tiny, tiny))?;
            add(&mut sum, &mut correction, -product)?;
            correction = checked(correction - error)?;
        }
        *out = checked(I::point(sum) + correction)?;
    }
    Ok(result)
}
/// Certify against the complete actual receiving chart, including its local
/// offset and source uncertainty. Both retained affine maps remain unchanged.
pub fn certify_between_frames(
    profile: &Profile,
    seat: Seat,
    options: Options,
    head: Placement,
    receiver: Frame,
) -> Result<BetweenFrameCertificate, Error> {
    propose(profile, seat, options)?;
    let literal_head = head;
    let head = interval(head)?;
    if !receiver
        .origin_mm
        .into_iter()
        .chain(receiver.axes.into_iter().flatten())
        .chain(receiver.offset_mm)
        .chain([receiver.source_error_mm])
        .all(f64::is_finite)
        || receiver.source_error_mm < 0.
    {
        return Err(Error::Invalid);
    }
    let uncertainty = I::point(receiver.source_error_mm) * I::point(2.);
    let radial_span = I::point(seat.mouth_radius_mm) - I::point(seat.bore_radius_mm);
    if uncertainty.hi >= seat.bore_radius_mm
        || uncertainty.hi >= seat.depth_mm
        || uncertainty.hi >= radial_span.lo
    {
        return Err(Error::Unresolved);
    }
    let receiver_axes = receiver.axes.map(|a| a.map(I::point));
    let gram = frame(receiver_axes)?;
    let operator = checked((I::point(1.) + I::point(gram)).sqrt())?.hi;
    let minimum = checked((I::point(1.) - I::point(gram)).sqrt())?.lo;
    let inverse_norm = checked(I::point(1.) / I::point(minimum))?.hi;
    let cofactors = [
        cross(receiver_axes[1], receiver_axes[2]),
        cross(receiver_axes[2], receiver_axes[0]),
        cross(receiver_axes[0], receiver_axes[1]),
    ];
    let determinant = checked(dot(receiver_axes[0], cofactors[0]))?;
    let inverse = cofactors.map(|row| row.map(|v| v / determinant));
    // Exact cancellation is qualified before applying the uncertain inverse.
    // No midpoint or newly rounded world receiving origin is introduced.
    let relative = affine_relative(literal_head, receiver)?;
    let origin = std::array::from_fn(|i| dot(inverse[i], relative));
    let axes = head
        .axes
        .map(|column| std::array::from_fn(|i| dot(inverse[i], column)));
    let chart_head = IntervalPlacement {
        origin,
        axes,
        offset: head.offset,
    };
    let chart_options = Options {
        tolerance_mm: checked(I::point(options.tolerance_mm) / I::point(operator))?.lo,
        minimum_land_mm: checked(I::point(options.minimum_land_mm) / I::point(minimum))?.hi,
    };
    let proposal = propose(profile, seat, chart_options)?;
    let (radius, height) = extents(&proposal.segments, proposal.reserve)?;
    let world_rounding = floating_bound(head, radius, height)? * I::point(inverse_norm);
    let extra = checked(world_rounding + I::point(receiver.source_error_mm))?;
    let chart = admit(
        proposal,
        seat,
        chart_options,
        chart_head,
        [I::point(0.); 3],
        [I::point(0.), I::point(0.), I::point(1.)],
        extra,
    )?;
    let maximum_physical_gap_mm =
        checked(I::point(chart.bearing.maximum_radial_gap_mm) * I::point(operator))?.hi;
    let land =
        I::point(chart.bearing.bearing_radial_mm[1]) - I::point(chart.bearing.bearing_radial_mm[0]);
    let minimum_physical_land_mm = checked(land * I::point(minimum))?.lo;
    Ok(BetweenFrameCertificate {
        chart,
        maximum_physical_gap_mm,
        minimum_physical_land_mm,
        receiver_operator_norm_upper: operator,
        receiver_inverse_norm_upper: inverse_norm,
    })
}

#[cfg(test)]
mod tests;
