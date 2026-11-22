//! ORIGINAL-REAL finite two-pole rational isocurves on analytic surfaces.
//! Positive weights earn the complete finite UV segment, not an infinite
//! support. A finite master window is mapped into that segment before PASS.
use super::*;
use crate::tapered_thread::interval::I;
mod endpoints;
/// Complete original finite-endpoint transfer, separate from curve coefficients.
/// For circles/ellipses the result also bounds `(a+b)*delta`, where delta is
/// each original endpoint's angular error from the retained parameter window.
/// The cylinder/ellipse certificate relies on this stronger contract when
/// transferring tangent bounds; a tighter physical-only reserve must expose
/// the angular bound separately before replacing this formula.
pub(super) fn endpoint_reserve(
    source: &Source<'_>,
    curve_id: u32,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<f64> {
    endpoints::reserve(source, curve_id, edge, budget)
}

fn finish(proof: Proof, target: f64, tol: f64) -> Proof {
    if !target.is_finite() || target < 0. {
        return Proof::Unsupported;
    }
    match proof {
        Proof::Certified(upper) => {
            let value = I::point(upper) + I::point(target);
            if value.valid() && value.hi <= tol {
                Proof::Certified(value.hi)
            } else {
                Proof::Unsupported
            }
        }
        Proof::Disagreement(lower) => {
            let value = I::point(lower) - I::point(target);
            if value.valid() && value.lo > tol {
                Proof::Disagreement(value.lo)
            } else {
                Proof::Unsupported
            }
        }
        Proof::Unsupported => Proof::Unsupported,
    }
}
fn periodic(
    image: analytic::Image,
    segment: &precision::AxisSegment,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    budget.charge(2048)?;
    let delta = segment.endpoints[1][segment.moving] - segment.endpoints[0][segment.moving];
    let sign = if delta.valid() && delta.lo > 0. {
        1.
    } else if delta.valid() && delta.hi < 0. {
        -1.
    } else {
        return Err(ReadError::Unsupported);
    };
    // actual_image orients SINE for positive UV curve parameter. Thus its
    // angular argument is sign*the ORIGINAL moving surface coordinate.
    let first = segment.endpoints[0][segment.moving] * I::point(sign);
    let last = segment.endpoints[1][segment.moving] * I::point(sign);
    if !first.valid() || !last.valid() {
        return Err(ReadError::Unsupported);
    }
    // Choosing the inward endpoints earns membership for EVERY actual REAL
    // endpoint enclosed here, rather than extending a rounded finite domain.
    let (mut low, mut high) = (first.hi, last.lo);
    if matches!(surface, Surface::Sphere { .. }) && segment.moving == 1 {
        let h = std::f64::consts::FRAC_PI_2.next_down();
        low = low.max(-h);
        high = high.min(h);
    }
    let mut best = f64::INFINITY;
    if low.is_finite()
        && high.is_finite()
        && low < high
        && let Some(base) = analytic::periodic_map(image, edge, 0., budget)?
    {
        // A rounded integer shift is only a proposal. Every shifted operator
        // is independently checked by the WHOLE periodic coefficient proof.
        let center = low / 2. + high / 2.;
        let target_center = edge.range[0] / 2. + edge.range[1] / 2.;
        let turns = ((center - target_center - base.phase) / std::f64::consts::TAU).round();
        if turns.is_finite() {
            for adjacent in [-1., 0., 1.] {
                budget.charge(512)?;
                let shift = (turns + adjacent) * std::f64::consts::TAU;
                let Some(map) = analytic::periodic_map(image, edge, shift, budget)? else {
                    continue;
                };
                let desired = edge.range.map(|t| I::point(t) + I::point(map.phase));
                if desired.into_iter().any(|p| !p.valid()) {
                    continue;
                }
                let chosen = edge.range.map(|t| (t + map.phase).clamp(low, high));
                if chosen.into_iter().any(|p| !p.is_finite()) || chosen[0] >= chosen[1] {
                    continue;
                }
                // Both endpoint maps are affine in target angle. Their
                // difference attains its maximum absolute value at an end;
                // the complete source derivative is bounded by map.speed.
                let differences = [
                    I::point(chosen[0]) - desired[0],
                    I::point(chosen[1]) - desired[1],
                ];
                if differences.into_iter().any(|d| !d.valid()) {
                    continue;
                }
                let angular = differences[0].magnitude().max(differences[1].magnitude());
                let error = I::point(map.error) + I::point(angular) * I::point(map.speed);
                if error.valid() {
                    best = best.min(error.hi)
                }
            }
        }
    }
    if best.is_finite() && best <= tol {
        return Ok(Proof::Certified(best));
    }
    // A point-to-ENTIRE continuation lower remains a lower to its finite
    // subset. A continuation upper, or a failed domain estimate, never PASS.
    Ok(match analytic::certify_image(image, edge, tol, budget)? {
        Proof::Disagreement(lower) => Proof::Disagreement(lower),
        _ => Proof::Unsupported,
    })
}
/// FINAL raw-source certificate: the coordinator must NOT apply another
/// decoded-image reserve. Its exact SI-radian ownership gate runs first.
#[allow(
    clippy::too_many_arguments,
    reason = "Actual source owner, identities, finite domain and original limits are independent proof inputs."
)]
pub(super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    definition: &uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    if !matches!(
        surface,
        Surface::Cylinder { .. }
            | Surface::Cone { .. }
            | Surface::Sphere { .. }
            | Surface::Torus { .. }
    ) {
        return Ok(None);
    }
    let Curve::Nurbs {
        curve,
        closed: false,
    } = &definition.geometry
    else {
        return Ok(None);
    };
    if (curve.degree() == 1 && (curve.points().len() != 2 || curve.knots().1 != [2, 2]))
        || curve.degree() > 11
        || curve.points().len() > 1024
    {
        return Ok(None);
    }
    budget.charge(512)?;
    if !tol.is_finite()
        || tol <= 0.
        || edge.range.into_iter().any(|t| !t.is_finite())
        || edge.range[0] >= edge.range[1]
    {
        return Err(ReadError::Invalid);
    }
    let Some(segment) = precision::axis_segment(source, definition, budget)? else {
        return Ok(None);
    };
    let image = precision::analytic_image(source, definition, surface, budget)?;
    let coefficients = precision::analytic_curve_reserve(source, curve_id, edge, budget)?;
    let target =
        I::point(coefficients) + I::point(endpoint_reserve(source, curve_id, edge, budget)?);
    if !target.valid() {
        return Err(ReadError::Unsupported);
    }
    let target = target.hi;
    let proof = match image {
        analytic::Image::Line { origin, direction } => {
            budget.charge(512)?;
            let end = std::array::from_fn(|i| origin[i] + direction[i]);
            affine::nurbs(
                1,
                (&[0., 1.], &[2, 2]),
                None,
                &[origin, end],
                true,
                edge,
                tol,
                budget,
            )?
        }
        analytic::Image::Periodic { .. } => periodic(image, &segment, surface, edge, tol, budget)?,
    };
    let proof = finish(proof, target, tol);
    if matches!(proof, Proof::Certified(_)) {
        precision::sense(source, curve_id, definition, surface, edge, budget)?;
    }
    budget.charge(64)?;
    Ok(Some(proof))
}
#[cfg(test)]
mod tests;
