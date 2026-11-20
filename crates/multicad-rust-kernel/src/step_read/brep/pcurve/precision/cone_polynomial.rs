//! Full ORIGINAL polynomial UV traces on an actual axial cone. Both the
//! position Taylor remainder and an INDEPENDENT derivative remainder are
//! checked over every complete cell. A recipe or sampled point grants nothing.
use super::super::Proof;
use super::polynomial_trace::Polynomial as P;
use super::*;
use crate::step_read::brep::check;
const SCRATCH: usize = 256 * 1024;
const MAX_BYTES: usize = 64 * 1024 * 1024;
type VecP = [P; 3];
fn exact_zero(v: I) -> bool {
    v.lo == 0. && v.hi == 0.
}
fn scaled(p: &P, v: I, b: &mut Budget<'_>) -> Result<P> {
    if v.lo == 1. && v.hi == 1. {
        b.charge(128)?;
        Ok(*p)
    } else {
        p.scale(v, b)
    }
}
fn coordinate_polynomial(controls: &[[I; 4]], axis: usize, b: &mut Budget<'_>) -> Result<P> {
    b.charge(256)?;
    check(
        !controls.is_empty() && controls.len() <= 12,
        ReadError::Unsupported,
    )?;
    let mut out = [I::point(0.); 12];
    for (i, p) in controls.iter().enumerate() {
        out[i] = p[axis];
    }
    P::new(&out[..controls.len()], b)
}
fn unit(controls: &[[I; 4]]) -> bool {
    controls.iter().all(|p| p[3].lo == 1. && p[3].hi == 1.)
}
fn power(x: I, n: usize) -> I {
    (0..n).fold(I::point(1.), |v, _| v * x)
}
fn factorial(n: usize) -> f64 {
    (1..=n).fold(1., |a, i| a * i as f64)
}
/// Taylor polynomials of sin and cos of `theta` about its hull midpoint to
/// order `n`, with the whole-field tail bound. Each series is evaluated by
/// Horner's rule in the offset, so no power is raised to the final degree.
fn taylor(theta: &P, n: usize, b: &mut Budget<'_>) -> Result<(P, P, f64)> {
    b.charge(4096)?;
    let center = theta.hull().midpoint();
    check(
        center.is_finite() && center.abs() <= 1e6,
        ReadError::Unsupported,
    )?;
    let delta = theta.add(&P::constant(I::point(-center), b)?, b)?;
    let range = delta.hull();
    check(range.valid(), ReadError::Unsupported)?;
    let (s, c) = (I::point(center).sin(), I::point(center).cos());
    check(s.valid() && c.valid(), ReadError::Unsupported)?;
    check(
        delta.degree().checked_mul(n).is_some_and(|d| d <= 16),
        ReadError::Unsupported,
    )?;
    // The j-th derivatives of sin and cos at the center.
    let coefficient = |j: usize| match j % 4 {
        0 => (s, c),
        1 => (c, -s),
        2 => (-s, -c),
        _ => (-c, s),
    };
    let series = |sine: bool, b: &mut Budget<'_>| -> Result<P> {
        let term = |j: usize| {
            let (a, d) = coefficient(j);
            (if sine { a } else { d }) / I::point(factorial(j))
        };
        let mut value = P::constant(term(n), b)?;
        for j in (0..n).rev() {
            b.charge(128)?;
            value = value
                .multiply(&delta, b)?
                .add(&P::constant(term(j), b)?, b)?;
        }
        Ok(value)
    };
    let sine = series(true, b)?;
    let cosine = series(false, b)?;
    // Scalar sine/cosine derivatives have absolute value<=1 on ALL real
    // angles; this is an actual whole-field tail, not a declared tolerance.
    let tail = power(I::point(range.magnitude()), n + 1) / I::point(factorial(n + 1));
    check(tail.valid() && tail.lo >= 0., ReadError::Unsupported)?;
    Ok((sine, cosine, tail.hi))
}
/// The cone-frame coordinates R^T (p - o) of a world polynomial curve. The
/// interval axes enclose the exact orthonormal frame, so distances measured
/// here are the world distances.
fn local(world: &VecP, frame: RawFrame, position: bool, b: &mut Budget<'_>) -> Result<VecP> {
    let zero = P::constant(I::point(0.), b)?;
    let mut out = [zero; 3];
    for (k, local) in out.iter_mut().enumerate() {
        b.charge(256)?;
        let mut value: Option<P> = None;
        let mut offset = I::point(0.);
        for (i, p) in world.iter().enumerate() {
            let coefficient = frame.axes[k][i];
            if exact_zero(coefficient) {
                continue;
            }
            offset = offset + coefficient * frame.origin[i];
            let term = scaled(p, coefficient, b)?;
            value = Some(match value {
                Some(v) => v.add(&term, b)?,
                None => term,
            });
        }
        let mut value = value.ok_or(ReadError::Unsupported)?;
        if position && !exact_zero(offset) {
            check(offset.valid(), ReadError::Unsupported)?;
            value = value.add(&P::constant(-offset, b)?, b)?;
        }
        *local = value;
    }
    Ok(out)
}
/// An upper bound on |a - c| over the cell: each component's Bernstein
/// residual bounds that component, and the components combine in quadrature.
fn component_error(a: &VecP, c: &VecP, b: &mut Budget<'_>) -> Result<f64> {
    let mut squares = I::point(0.);
    for k in 0..3 {
        let residual = a[k].add(&c[k].scale(I::point(-1.), b)?, b)?;
        b.charge((residual.degree() + 1) * 64 + 64)?;
        let bound = residual.hull().magnitude();
        check(bound.is_finite(), ReadError::Unsupported)?;
        squares = squares + I::point(bound).square();
    }
    let error = squares.sqrt();
    check(error.valid(), ReadError::Unsupported)?;
    Ok(error.hi)
}
fn derivative(p: VecP, b: &mut Budget<'_>) -> Result<VecP> {
    Ok([
        p[0].derivative(b)?,
        p[1].derivative(b)?,
        p[2].derivative(b)?,
    ])
}
struct Cone {
    frame: RawFrame,
    radius: I,
    slope: I,
}
/// The position bound at Taylor order `n`, in the cone frame.
fn position_error(
    rho: &P,
    theta: &P,
    z: &P,
    actual: &VecP,
    n: usize,
    b: &mut Budget<'_>,
) -> Result<f64> {
    let (s, c, tail) = taylor(theta, n, b)?;
    let expected = [rho.multiply(&c, b)?, rho.multiply(&s, b)?, *z];
    let coefficient = component_error(&expected, actual, b)?;
    // The two radial scalar errors give this Euclidean bound; no one common
    // mean-value point is assumed.
    let radial_tail = I::point(rho.hull().magnitude()) * I::point(tail);
    let position = I::point(coefficient) + norm([radial_tail, radial_tail, I::point(0.)]);
    check(position.valid(), ReadError::Unsupported)?;
    Ok(position.hi)
}
/// Whether the analytic first derivative at Taylor order `n` agrees in
/// direction with the target derivative everywhere on the cell.
fn tangent_agrees(
    rho: &P,
    theta: &P,
    z: &P,
    actual: &VecP,
    n: usize,
    b: &mut Budget<'_>,
) -> Result<bool> {
    let (s, c, dtail) = taylor(theta, n, b)?;
    let rho_derivative = rho.derivative(b)?;
    let theta_derivative = theta.derivative(b)?;
    let radial_angular = rho.multiply(&theta_derivative, b)?;
    let dx = rho_derivative
        .multiply(&c, b)?
        .add(&radial_angular.multiply(&s, b)?.scale(I::point(-1.), b)?, b)?;
    let dy = rho_derivative
        .multiply(&s, b)?
        .add(&radial_angular.multiply(&c, b)?, b)?;
    let dz = z.derivative(b)?;
    let actual_derivative = derivative(*actual, b)?;
    let derivative_coefficient = component_error(&[dx, dy, dz], &actual_derivative, b)?;
    let amplitude =
        I::point(rho_derivative.hull().magnitude()) + I::point(radial_angular.hull().magnitude());
    let radial_tail = amplitude * I::point(dtail);
    let error = I::point(derivative_coefficient) + norm([radial_tail, radial_tail, I::point(0.)]);
    // This tail approximates the ANALYTIC FIRST derivative independently.
    // It is never obtained by differentiating an unknown position remainder.
    let source_speed =
        (rho_derivative.hull().square() + radial_angular.hull().square() + dz.hull().square())
            .sqrt();
    let target_speed = norm(actual_derivative.map(|p| p.hull()));
    check(
        error.valid() && source_speed.valid() && target_speed.valid(),
        ReadError::Unsupported,
    )?;
    // If ||s-t|| is strictly below EITHER proved speed lower, then
    // s·t>0 by Cauchy-Schwarz. Taking the larger independent lower bound
    // chooses that sufficient proof; no sampled direction is trusted.
    let lower = source_speed.lo.max(target_speed.lo);
    Ok(lower.is_finite() && lower > 0. && error.hi < lower)
}
/// The value of a Bernstein polynomial at `t` by de Casteljau.
fn at(p: &P, t: I, b: &mut Budget<'_>) -> Result<I> {
    let controls = p.controls();
    b.charge(controls.len() * controls.len() * 16 + 64)?;
    let mut values = controls.to_vec();
    for n in (1..values.len()).rev() {
        for i in 0..n {
            values[i] = (I::point(1.) - t) * values[i] + t * values[i + 1];
        }
    }
    check(values[0].valid(), ReadError::Unsupported)?;
    Ok(values[0])
}
/// Direction agreement from the cell midpoint and second-derivative hulls:
/// with f = E' - A', |f(t)| <= |f(1/2)| + sup|f'| / 2 on [0, 1], and each
/// speed is at least its midpoint value less half its second-derivative
/// hull. A cheap sufficient proof; `false` leaves the full expansion to
/// decide.
fn tangent_by_hulls(rho: &P, theta: &P, z: &P, actual: &VecP, b: &mut Budget<'_>) -> Result<bool> {
    b.charge(2048)?;
    let m = I::point(0.5);
    let d1 = |p: &P, b: &mut Budget<'_>| p.derivative(b);
    let (rho1, theta1, z1) = (d1(rho, b)?, d1(theta, b)?, d1(z, b)?);
    let (rho2, theta2, z2) = (d1(&rho1, b)?, d1(&theta1, b)?, d1(&z1, b)?);
    let (r, t, r1, t1, dz) = (
        at(rho, m, b)?,
        at(theta, m, b)?,
        at(&rho1, m, b)?,
        at(&theta1, m, b)?,
        at(&z1, m, b)?,
    );
    let (sin, cos) = (t.sin(), t.cos());
    let angular = r * t1;
    let expected = [r1 * cos - angular * sin, r1 * sin + angular * cos, dz];
    let actual1 = derivative(*actual, b)?;
    let actual2 = derivative(actual1, b)?;
    let mut target = [I::point(0.); 3];
    for (k, value) in target.iter_mut().enumerate() {
        *value = at(&actual1[k], m, b)?;
    }
    let difference = norm(std::array::from_fn(|k| expected[k] - target[k]));
    // In the rotating frame the radial and angular second-derivative parts
    // are (rho'' - rho theta'^2, 2 rho' theta' + rho theta'').
    let h = |p: &P| I::point(p.hull().magnitude());
    let (hr, ht1) = (h(rho), h(&theta1));
    let radial = h(&rho2) + hr * ht1 * ht1;
    let tangential = I::point(2.) * h(&rho1) * ht1 + hr * h(&theta2);
    let source_second = (radial.square() + tangential.square() + h(&z2).square()).sqrt();
    let target_second = norm(actual2.map(|p| h(&p)));
    let half = I::point(0.5);
    let error = difference + half * (source_second + target_second);
    let source_speed = norm(expected) - half * source_second;
    let target_speed = norm(target) - half * target_second;
    check(
        error.valid() && source_speed.valid() && target_speed.valid(),
        ReadError::Unsupported,
    )?;
    // As in the full proof, ||s - t|| below either proved speed gives s.t > 0.
    let lower = source_speed.lo.max(target_speed.lo);
    Ok(lower.is_finite() && lower > 0. && error.hi < lower)
}
/// The lowest Taylor order, up to `order`, whose whole-field tail alone
/// leaves at least half the tolerance for the coefficient residual.
fn first_order(theta: &P, rho: &P, order: usize, tol: f64) -> usize {
    let range = theta.hull();
    let range = I::point(0.5 * (range.hi - range.lo));
    let radius = I::point(rho.hull().magnitude());
    (1..order)
        .find(|&n| {
            let tail = radius * power(range, n + 1) / I::point(factorial(n + 1));
            tail.valid() && tail.hi <= 0.5 * tol
        })
        .unwrap_or(order)
}
fn cell(
    cone: &Cone,
    uv: &[[I; 4]],
    target: &[[I; 4]],
    tol: f64,
    b: &mut Budget<'_>,
) -> Result<Proof> {
    b.charge(8192)?;
    check(unit(uv) && unit(target), ReadError::Unsupported)?;
    let (order, derivative_order) = match (uv.len(), target.len()) {
        (4, 4) => (4, 3),
        (2, 12) => (10, 9),
        _ => return Err(ReadError::Unsupported),
    };
    let theta = coordinate_polynomial(uv, 0, b)?;
    let z = coordinate_polynomial(uv, 1, b)?;
    let rho = z
        .scale(cone.slope, b)?
        .add(&P::constant(cone.radius, b)?, b)?;
    let rh = rho.hull();
    check(rh.valid() && rh.lo > 0., ReadError::Unsupported)?;
    let world = [
        coordinate_polynomial(target, 0, b)?,
        coordinate_polynomial(target, 1, b)?,
        coordinate_polynomial(target, 2, b)?,
    ];
    let actual = local(&world, cone.frame, true, b)?;
    // Lower orders are tried first; the original order remains the last
    // attempt, so a cell certified before is still certified.
    let mut position = None;
    for n in first_order(&theta, &rho, order, tol)..=order {
        let bound = position_error(&rho, &theta, &z, &actual, n, b)?;
        if bound <= tol {
            position = Some(bound);
            break;
        }
    }
    let Some(position) = position else {
        return Ok(Proof::Unsupported);
    };
    // Direction agreement needs only an error below the proved speed: a
    // first-order expansion usually suffices.
    let mut agrees = tangent_by_hulls(&rho, &theta, &z, &actual, b)?;
    for n in [1, derivative_order] {
        if agrees {
            break;
        }
        if tangent_agrees(&rho, &theta, &z, &actual, n, b)? {
            agrees = true;
            break;
        }
    }
    if !agrees {
        return Ok(Proof::Unsupported);
    }
    b.charge(128)?;
    Ok(Proof::Certified(position))
}
fn full_target(source: &Source<'_>, id: u32, edge: &Edge, b: &mut Budget<'_>) -> Result<()> {
    b.charge(256)?;
    let f = fields(source, id, 1, b)?;
    let controls = list(f.controls)?;
    let ends = [
        controls.first().ok_or(ReadError::Unsupported)?,
        controls.last().ok_or(ReadError::Unsupported)?,
    ];
    for (i, end) in ends.into_iter().enumerate() {
        let vertex = u32::try_from(edge.vertices[i]).map_err(|_| ReadError::Unsupported)?;
        let [_, point] = source.args(vertex, "VERTEX_POINT")? else {
            return Err(ReadError::Invalid);
        };
        check(reference(point)? == reference(end)?, ReadError::Unsupported)?;
    }
    // Clamping and exact original end multiplicities are earned by the
    // source-bound extractor. SAME original point IDs bind the entire actual
    // topological edge to that natural curve, with no cached-window oracle.
    Ok(())
}
fn knot_values<'a, 'i>(
    source: &'a Source<'i>,
    id: u32,
    b: &mut Budget<'_>,
) -> Result<&'a [Value<'i>]> {
    let f = fields(source, id, 1, b)?;
    let Some((_, k)) = f.knots[0] else {
        return Err(ReadError::Unsupported);
    };
    list(k)
}
fn equal_charts(
    source: &Source<'_>,
    uv: u32,
    target: u32,
    linear: bool,
    b: &mut Budget<'_>,
) -> Result<()> {
    let a = knot_values(source, uv, b)?;
    let c = knot_values(source, target, b)?;
    if linear {
        check(a.len() == 2 && c.len() >= 2, ReadError::Unsupported)?;
        for (u, t) in [(0, 0), (1, c.len() - 1)] {
            b.charge(128)?;
            check(
                exact_same(source, (uv, &a[u]), (target, &c[t]), b)?,
                ReadError::Unsupported,
            )?;
        }
    } else {
        check(a.len() == c.len(), ReadError::Unsupported)?;
        for (u, t) in a.iter().zip(c) {
            b.charge(128)?;
            check(
                exact_same(source, (uv, u), (target, t), b)?,
                ReadError::Unsupported,
            )?;
        }
    }
    Ok(())
}
fn restrict_linear(
    uv: &[[I; 4]],
    range: [I; 2],
    full: [I; 2],
    b: &mut Budget<'_>,
) -> Result<[[I; 4]; 2]> {
    b.charge(1024)?;
    check(uv.len() == 2 && unit(uv), ReadError::Unsupported)?;
    let width = full[1] - full[0];
    check(width.valid() && width.lo > 0., ReadError::Unsupported)?;
    let mut q = [I::point(0.); 2];
    for (i, t) in range.into_iter().enumerate() {
        let value = (t - full[0]) / width;
        check(value.valid(), ReadError::Unsupported)?;
        // Actual original knot order and exact matched global endpoints
        // prove membership. The interval is intersected, never the source.
        q[i] = I::new(value.lo.max(0.), value.hi.min(1.));
        check(q[i].valid(), ReadError::Unsupported)?;
    }
    check((q[1] - q[0]).lo > 0., ReadError::Unsupported)?;
    Ok(q.map(|q| {
        std::array::from_fn(|i| {
            if i == 3 {
                I::point(1.)
            } else {
                uv[0][i] + (uv[1][i] - uv[0][i]) * q
            }
        })
    }))
}
#[allow(
    clippy::too_many_arguments,
    reason = "Original owner, identities, finite domains and original limits are independent proof inputs."
)]
pub(super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    b: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    let (
        Surface::Cone { .. },
        Curve::Nurbs {
            curve: uv,
            closed: false,
        },
        Curve::Nurbs {
            curve: target,
            closed: false,
        },
    ) = (surface, &definition.geometry, &edge.curve)
    else {
        return Ok(None);
    };
    let linear = match (uv.degree(), target.degree()) {
        (3, 3) => false,
        (1, 11) => true,
        _ => return Ok(None),
    };
    b.charge(4096)?;
    check(
        uv.points().len() <= 1024 && target.points().len() <= 1024,
        ReadError::Unsupported,
    )?;
    // Two bounded8MiB extraction phases (first result retained), fixed
    // numeric stack, decoded argument vectors and caller cache preflight.
    const ENVELOPE: usize = 2 * 8 * 1024 * 1024 + 4 * 1024 * 1024 + 256 * 1024 + SCRATCH;
    check(ENVELOPE <= MAX_BYTES, ReadError::Limit)?;
    check(tol.is_finite() && tol > 0., ReadError::Invalid)?;
    check(super::super::si_radians(source, b)?, ReadError::Unsupported)?;
    let target_id = unwrapped(source, curve_id, b)?;
    full_target(source, target_id, edge, b)?;
    equal_charts(source, definition.curve_id, target_id, linear, b)?;
    let [_, placement, radius, angle] = source.args(definition.surface, "CONICAL_SURFACE")? else {
        return Err(ReadError::Invalid);
    };
    let frame = frame(source, reference(placement)?, 3, b)?;
    let radius = scalar(source, definition.surface, radius, b)?;
    let angle = scalar(source, definition.surface, angle, b)?;
    check(
        radius.valid()
            && radius.lo > 0.
            && angle.valid()
            && angle.lo > 0.
            && angle.hi < std::f64::consts::FRAC_PI_2,
        ReadError::Unsupported,
    )?;
    let cosine = angle.cos();
    check(cosine.valid() && cosine.lo > 0., ReadError::Unsupported)?;
    let slope = angle.sin() / cosine;
    check(slope.valid(), ReadError::Unsupported)?;
    let cone = Cone {
        frame,
        radius,
        slope,
    };
    let uv_data = super::bezier::original_curve(
        source,
        definition.curve_id,
        uv,
        super::bezier::Coordinates::Raw(2),
        b,
    )?;
    let target_data = super::bezier::original_curve(
        source,
        target_id,
        target,
        super::bezier::Coordinates::Raw(3),
        b,
    )?;
    check(
        uv_data.matches(source, definition.curve_id)
            && matches!(uv_data.coordinates(), super::bezier::Coordinates::Raw(2))
            && target_data.matches(source, target_id)
            && matches!(
                target_data.coordinates(),
                super::bezier::Coordinates::Raw(3)
            ),
        ReadError::Invalid,
    )?;
    let resident = uv_data
        .retained_bytes()
        .checked_add(target_data.retained_bytes())
        .and_then(|v| v.checked_add(SCRATCH))
        .ok_or(ReadError::Limit)?;
    let peak = uv_data
        .peak_bytes()
        .max(
            uv_data
                .retained_bytes()
                .checked_add(target_data.peak_bytes())
                .ok_or(ReadError::Limit)?,
        )
        .max(resident);
    check(peak <= MAX_BYTES, ReadError::Limit)?;
    if linear {
        check(uv_data.cells().len() == 1, ReadError::Unsupported)?;
    } else {
        check(
            uv_data.cells().len() == target_data.cells().len(),
            ReadError::Unsupported,
        )?;
    }
    let first = target_data
        .cells()
        .first()
        .ok_or(ReadError::Unsupported)?
        .domain()[0];
    let last = target_data
        .cells()
        .last()
        .ok_or(ReadError::Unsupported)?
        .domain()[1];
    let mut maximum = 0_f64;
    for (i, target) in target_data.cells().iter().enumerate() {
        b.charge(256)?;
        let result = if linear {
            let uv = restrict_linear(
                uv_data.cells()[0].controls(),
                target.domain(),
                [first, last],
                b,
            )?;
            cell(&cone, &uv, target.controls(), tol, b)?
        } else {
            let u = &uv_data.cells()[i];
            check(
                u.knot_indices() == target.knot_indices(),
                ReadError::Invalid,
            )?;
            cell(&cone, u.controls(), target.controls(), tol, b)?
        };
        match result {
            Proof::Certified(value) => maximum = maximum.max(value),
            _ => return Ok(Some(Proof::Unsupported)),
        }
    }
    b.charge(128)?;
    Ok(Some(Proof::Certified(maximum)))
}
#[cfg(test)]
mod tests;
