//! A whole-set separating slab for constant ORIGINAL cone V traces.
//! Positive rational basis controls bound EVERY source/target point; a finite
//! used subset cannot leave those slabs. No sampled point or failed upper
//! creates Invalid, and overlap grants no positive geometry permission.
use super::super::Proof;
use super::*;
use crate::brep::NurbsCurve;
use crate::step_read::brep::check;
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn basis<'s, 'i>(
    source: &'s Source<'i>,
    id: u32,
    curve: &NurbsCurve,
    dim: usize,
    budget: &mut Budget<'_>,
) -> Result<Fields<'s, 'i>> {
    let field = fields(source, id, 1, budget)?;
    // This existing proof earns original degree/multiplicity/order, all
    // original positive weights and finite control bounds. Natural periodic
    // wrapping, if declared, remains inside the SAME convex control hull.
    nurbs_error(
        source,
        id,
        &field,
        &[curve.degree()],
        &[curve.knots()],
        curve.points(),
        curve.weights(),
        dim,
        None,
        budget,
    )?;
    Ok(field)
}
fn trig(value: I, cosine: bool) -> Result<I> {
    check(value.valid(), ReadError::Unsupported)?;
    let middle = value.midpoint();
    let radius = value - I::point(middle);
    check(middle.is_finite() && radius.valid(), ReadError::Unsupported)?;
    let at = if cosine {
        I::point(middle).cos()
    } else {
        I::point(middle).sin()
    };
    let result = I::new(
        (at.lo - radius.magnitude()).next_down().max(-1.),
        (at.hi + radius.magnitude()).next_up().min(1.),
    );
    check(result.valid(), ReadError::Unsupported)?;
    Ok(result)
}
/// Lower-only: None means no separating certificate, never successful C0.
#[allow(
    clippy::too_many_arguments,
    reason = "Actual source owner, original identities, both curves and original limits remain explicit."
)]
pub(super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    let (
        Surface::Cone {
            frame: cached_frame,
            ..
        },
        Curve::Nurbs { curve: uv, .. },
        Curve::Nurbs { curve: target, .. },
    ) = (surface, &definition.geometry, &edge.curve)
    else {
        return Ok(None);
    };
    budget.charge(8192)?;
    check(tol.is_finite() && tol > 0., ReadError::Invalid)?;
    check(
        source.angle == 1. && super::super::si_radians(source, budget)?,
        ReadError::Unsupported,
    )?;
    let uv_fields = basis(source, definition.curve_id, uv, 2, budget)?;
    let uv_points = list(uv_fields.controls)?;
    let first = reference(uv_points.first().ok_or(ReadError::Unsupported)?)?;
    let [_, coords] = source.args(first, "CARTESIAN_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let coords = list(coords)?;
    check(coords.len() == 2, ReadError::Unsupported)?;
    let value = &coords[1];
    let v = scalar(source, first, value, budget)?;
    for reference_value in uv_points {
        budget.charge(128)?;
        let id = reference(reference_value)?;
        let [_, coords] = source.args(id, "CARTESIAN_POINT")? else {
            return Err(ReadError::Invalid);
        };
        let coords = list(coords)?;
        check(coords.len() == 2, ReadError::Unsupported)?;
        if !exact_same(source, (first, value), (id, &coords[1]), budget)? {
            return Ok(None);
        }
    }
    let [_, placement, radius, angle] = source.args(definition.surface, "CONICAL_SURFACE")? else {
        return Err(ReadError::Invalid);
    };
    let frame = frame(source, reference(placement)?, 3, budget)?;
    let radius = scalar(source, definition.surface, radius, budget)?;
    let angle = scalar(source, definition.surface, angle, budget)?;
    check(
        radius.valid()
            && radius.lo >= 0.
            && angle.valid()
            && angle.lo > 0.
            && angle.hi < std::f64::consts::FRAC_PI_2,
        ReadError::Unsupported,
    )?;
    let cosine = trig(angle, true)?;
    check(cosine.lo > 0., ReadError::Unsupported)?;
    let rho = radius + v * (trig(angle, false)? / cosine);
    check(rho.valid(), ReadError::Unsupported)?;
    // q is only a proposed separating direction. All projections use actual
    // REAL frame/control enclosures; no orthogonality of cached q is assumed.
    let q = cached(cached_frame.z);
    let length = norm(q);
    check(length.valid() && length.lo > 0., ReadError::Unsupported)?;
    let source_slab = dot(q, frame.origin)
        + v * dot(q, frame.axes[2])
        + rho * (dot(q, frame.axes[0]) * I::new(-1., 1.) + dot(q, frame.axes[1]) * I::new(-1., 1.));
    check(source_slab.valid(), ReadError::Unsupported)?;
    let target_id = unwrapped(source, curve_id, budget)?;
    let target_fields = basis(source, target_id, target, 3, budget)?;
    let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in list(target_fields.controls)? {
        budget.charge(256)?;
        let projection = dot(q, point(source, reference(p)?, 3, budget)?);
        check(projection.valid(), ReadError::Unsupported)?;
        low = low.min(projection.lo);
        high = high.max(projection.hi);
    }
    check(
        low.is_finite() && high.is_finite() && low <= high,
        ReadError::Unsupported,
    )?;
    let above = I::point(source_slab.lo) - I::point(high);
    let below = I::point(low) - I::point(source_slab.hi);
    check(above.valid() && below.valid(), ReadError::Unsupported)?;
    let gap = above.lo.max(below.lo).max(0.);
    let lower = I::point(gap) / I::point(length.hi);
    check(lower.valid(), ReadError::Unsupported)?;
    budget.charge(64)?;
    Ok((lower.lo > tol).then_some(Proof::Disagreement(lower.lo)))
}
#[cfg(test)]
mod tests;
