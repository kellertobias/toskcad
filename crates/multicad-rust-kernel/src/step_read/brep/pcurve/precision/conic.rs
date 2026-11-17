//! Whole original rational quadratic plane curves against an actual conic.
//! At most four positive, open-half-plane Bezier arcs earn exactly one
//! co-oriented winding. Neither sampled points nor a failed upper admit it.
use super::*;
use crate::brep::Curve;
const MAX_CONTROLS: usize = 9;
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn multiply<const N: usize>(
    a: &[I],
    b: &[I],
    bin_a: &[f64],
    bin_b: &[f64],
    bin_out: &[f64],
) -> [I; N] {
    let mut result = [I::point(0.); N];
    for (i, &a) in a.iter().enumerate() {
        for (j, &b) in b.iter().enumerate() {
            let weight = I::point(bin_a[i]) * I::point(bin_b[j]) / I::point(bin_out[i + j]);
            result[i + j] = result[i + j] + weight * a * b;
        }
    }
    result
}
fn closure(source: &Source<'_>, controls: &[Value<'_>], budget: &mut Budget<'_>) -> Result<bool> {
    let (first, last) = (
        reference(&controls[0])?,
        reference(&controls[controls.len() - 1])?,
    );
    if first == last {
        return Ok(true);
    }
    let [_, a] = source.args(first, "CARTESIAN_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let [_, b] = source.args(last, "CARTESIAN_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let (a, b) = (list(a)?, list(b)?);
    if a.len() != 2 || b.len() != 2 {
        return Err(ReadError::Unsupported);
    }
    for i in 0..2 {
        if rational(source, first, &a[i], budget)? != rational(source, last, &b[i], budget)? {
            return Ok(false);
        }
    }
    Ok(true)
}
fn segments(
    source: &Source<'_>,
    id: u32,
    field: &Fields<'_, '_>,
    curve: &crate::brep::NurbsCurve,
    budget: &mut Budget<'_>,
) -> Result<Vec<usize>> {
    integer(source, id, field.degrees[0], 2, budget)?;
    let Some((multiplicities, values)) = field.knots[0] else {
        return if curve.points().len() == 3 {
            Ok(vec![0])
        } else {
            Err(ReadError::Unsupported)
        };
    };
    let (multiplicities, values) = (list(multiplicities)?, list(values)?);
    let (_, counts) = curve.knots();
    if values.len() != counts.len() || multiplicities.len() != counts.len() {
        return Err(ReadError::Unsupported);
    }
    let mut expanded = Vec::new();
    expanded
        .try_reserve_exact(curve.points().len() + 3)
        .map_err(|_| ReadError::Limit)?;
    let mut previous: Option<I> = None;
    for i in 0..values.len() {
        budget.charge(256)?;
        integer(source, id, &multiplicities[i], counts[i], budget)?;
        let value = scalar(source, id, &values[i], budget)?;
        if let Some(previous) = previous {
            let gap = value - previous;
            if !gap.valid() || gap.lo <= 0. {
                return Err(ReadError::Unsupported);
            }
        }
        previous = Some(value);
        if counts[i] > curve.points().len() + 3 - expanded.len() {
            return Err(ReadError::Unsupported);
        }
        expanded.extend(std::iter::repeat_n(i, counts[i]));
    }
    if expanded.len() != curve.points().len() + 3 {
        return Err(ReadError::Unsupported);
    }
    let mut result = Vec::new();
    result.try_reserve_exact(4).map_err(|_| ReadError::Limit)?;
    for k in 2..curve.points().len() {
        if expanded[k] == expanded[k + 1] {
            continue;
        }
        if counts[expanded[k]] < 2 || counts[expanded[k + 1]] < 2 {
            return Err(ReadError::Unsupported);
        }
        if k > 2 && counts[expanded[k]] != 2 {
            return Err(ReadError::Unsupported);
        }
        let start = k - 2;
        if result.last().is_some_and(|&p| p + 2 != start) {
            return Err(ReadError::Unsupported);
        }
        result.push(start);
    }
    if result.is_empty()
        || result.len() > 4
        || result[0] != 0
        || result[result.len() - 1] + 2 != curve.points().len() - 1
    {
        return Err(ReadError::Unsupported);
    }
    Ok(result)
}
/// An actual coefficient certificate, already in ORIGINAL REAL geometry.
/// The target's raw-to-retained transfer is included exactly once.
pub(super) fn certify(
    source: &Source<'_>,
    curve3d_id: u32,
    definition: &super::super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::super::Proof>> {
    let (
        Surface::Plane(_),
        Curve::Nurbs {
            curve,
            closed: false,
        },
    ) = (surface, &definition.geometry)
    else {
        return Ok(None);
    };
    if !matches!(edge.curve, Curve::Circle { .. } | Curve::Ellipse { .. }) {
        return Ok(None);
    }
    budget.charge(2048)?;
    if edge.range.into_iter().any(|t| !t.is_finite())
        || edge.range[0] >= edge.range[1]
        || edge.range[1] - edge.range[0] > std::f64::consts::TAU.next_up()
    {
        return Err(ReadError::Unsupported);
    }
    if !tol.is_finite()
        || tol <= 0.
        || curve.degree() != 2
        || !(3..=MAX_CONTROLS).contains(&curve.points().len())
    {
        return Ok(Some(super::super::Proof::Unsupported));
    }
    let id = unwrapped(source, curve3d_id, budget)?;
    let entity = source.simple(id)?;
    let (axis, a, b) = match (entity.name, entity.args.as_slice()) {
        ("CIRCLE", [_, axis, r]) => {
            let r = scalar(source, id, r, budget)?;
            (reference(axis)?, r, r)
        }
        ("ELLIPSE", [_, axis, a, b]) => (
            reference(axis)?,
            scalar(source, id, a, budget)?,
            scalar(source, id, b, budget)?,
        ),
        _ => return Err(ReadError::Unsupported),
    };
    if !a.valid() || !b.valid() || a.lo <= 0. || b.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let target = frame(source, axis, 3, budget)?;
    let field = fields(source, definition.curve_id, 1, budget)?;
    let controls = list(field.controls)?;
    if controls.len() != curve.points().len() || !closure(source, controls, budget)? {
        return Ok(Some(super::super::Proof::Unsupported));
    }
    let spans = segments(source, definition.curve_id, &field, curve, budget)?;
    let mut weights = [I::point(1.); MAX_CONTROLS];
    if let Some(values) = field.weights {
        let values = list(values)?;
        if values.len() != controls.len() {
            return Err(ReadError::Unsupported);
        }
        for (i, value) in values.iter().enumerate() {
            weights[i] = scalar(source, definition.curve_id, value, budget)?;
            if !weights[i].valid() || weights[i].lo <= 0. {
                return Err(ReadError::Unsupported);
            }
        }
    }
    let maximum = weights[..controls.len()]
        .iter()
        .map(|x| x.hi)
        .fold(0_f64, f64::max);
    if !maximum.is_finite() || maximum <= 0. {
        return Err(ReadError::Unsupported);
    }
    let mut homogeneous = [[I::point(0.); 4]; MAX_CONTROLS];
    for (i, control) in controls.iter().enumerate() {
        // Fixed-size frame projection, normalized homogeneous coefficients,
        // positive-weight checks/copies: <160 interval/control actions.
        budget.charge(512)?;
        let world =
            super::analytic::plane_point(source, definition.surface, reference(control)?, budget)?;
        let delta = difference(world, target.origin);
        let w = weights[i] / I::point(maximum);
        homogeneous[i] = [
            dot(delta, target.axes[0]) * w / a,
            dot(delta, target.axes[1]) * w / b,
            dot(delta, target.axes[2]) * w,
            w,
        ];
        if !homogeneous[i].into_iter().all(I::valid) || w.lo <= 0. {
            return Err(ReadError::Unsupported);
        }
    }
    let mut upper = 0_f64;
    for start in spans {
        // <=36 quadratic-product pairs, <=12 derivative-product pairs,
        // coefficients/normalizations, positive-ray and interval checks.
        budget.charge(16_384)?;
        let h = &homogeneous[start..start + 3];
        let x = [h[0][0], h[1][0], h[2][0]];
        let y = [h[0][1], h[1][1], h[2][1]];
        let w = [h[0][3], h[1][3], h[2][3]];
        let min_w = w.into_iter().map(|w| w.lo).fold(f64::INFINITY, f64::min);
        if !min_w.is_finite() || min_w <= 0. {
            return Err(ReadError::Unsupported);
        }
        let xx = multiply::<5>(&x, &x, &[1., 2., 1.], &[1., 2., 1.], &[1., 4., 6., 4., 1.]);
        let yy = multiply::<5>(&y, &y, &[1., 2., 1.], &[1., 2., 1.], &[1., 4., 6., 4., 1.]);
        let ww = multiply::<5>(&w, &w, &[1., 2., 1.], &[1., 2., 1.], &[1., 4., 6., 4., 1.]);
        let mut residual = 0_f64;
        for i in 0..5 {
            let f = xx[i] + yy[i] - ww[i];
            if !f.valid() {
                return Err(ReadError::Unsupported);
            }
            residual = residual.max(f.magnitude());
        }
        let radial_squared = I::point(residual) / I::point(min_w).square();
        if !radial_squared.valid() || radial_squared.hi >= 1. {
            return Ok(Some(super::super::Proof::Unsupported));
        }
        let rho = (I::point(1.) - radial_squared).sqrt();
        if !rho.valid() || rho.lo <= 0. {
            return Ok(Some(super::super::Proof::Unsupported));
        }
        let radial = radial_squared / (I::point(1.) + I::point(rho.lo));
        let z = h.iter().map(|p| p[2].magnitude()).fold(0_f64, f64::max);
        let error = I::point(z) / I::point(min_w) + I::point(a.hi.max(b.hi)) * radial;
        if !error.valid() {
            return Err(ReadError::Unsupported);
        }
        upper = upper.max(error.hi);
        let dx = [I::point(2.) * (x[1] - x[0]), I::point(2.) * (x[2] - x[1])];
        let dy = [I::point(2.) * (y[1] - y[0]), I::point(2.) * (y[2] - y[1])];
        let xy = multiply::<4>(&x, &dy, &[1., 2., 1.], &[1., 1.], &[1., 3., 3., 1.]);
        let yx = multiply::<4>(&y, &dx, &[1., 2., 1.], &[1., 1.], &[1., 3., 3., 1.]);
        if (0..4).any(|i| !(xy[i] - yx[i]).valid() || (xy[i] - yx[i]).lo <= 0.) {
            return Ok(Some(super::super::Proof::Unsupported));
        }
        // An open half-plane containing all three actual homogeneous XY
        // controls contains every ray of this arc. Its positive angle gain
        // is strictly <pi. At most four such arcs, exact closure and omega>0
        // imply 0<total<4pi: exactly ONE winding, not merely some cover.
        let proposal = [
            (x[0] * I::point(0.25) + x[1] * I::point(0.5) + x[2] * I::point(0.25)).midpoint(),
            (y[0] * I::point(0.25) + y[1] * I::point(0.5) + y[2] * I::point(0.25)).midpoint(),
        ];
        if proposal.iter().any(|v| !v.is_finite()) {
            return Err(ReadError::Unsupported);
        }
        for i in 0..3 {
            let dot = x[i] * I::point(proposal[0]) + y[i] * I::point(proposal[1]);
            if !dot.valid() || dot.lo <= 0. {
                return Ok(Some(super::super::Proof::Unsupported));
            }
        }
    }
    // Every finite positive target edge window is covered by this once-wound
    // closed geometric trace. A seam crossing is two in-domain UV windows;
    // no extrapolation, sampled closure or equal-parameter claim is made.
    let target_error = super::analytic::curve_reserve(source, id, edge, budget)?;
    let total = I::point(upper) + I::point(target_error);
    if !total.valid() {
        return Err(ReadError::Unsupported);
    }
    budget.charge(64)?;
    Ok(Some(if total.hi <= tol {
        super::super::Proof::Certified(total.hi)
    } else {
        super::super::Proof::Unsupported
    }))
}
#[cfg(test)]
mod tests;
