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
    let mut sine = P::constant(I::point(0.), b)?;
    let mut cosine = sine;
    let mut term = P::constant(I::point(1.), b)?;
    for j in 0..=n {
        b.charge(128)?;
        let (a, d) = match j % 4 {
            0 => (s, c),
            1 => (c, -s),
            2 => (-s, -c),
            _ => (-c, s),
        };
        let divisor = I::point(factorial(j));
        sine = sine.add(&term.scale(a / divisor, b)?, b)?;
        cosine = cosine.add(&term.scale(d / divisor, b)?, b)?;
        if j < n {
            term = term.multiply(&delta, b)?;
        }
    }
    // Scalar sine/cosine derivatives have absolute value<=1 on ALL real
    // angles; this is an actual whole-field tail, not a declared tolerance.
    let tail = power(I::point(range.magnitude()), n + 1) / I::point(factorial(n + 1));
    check(tail.valid() && tail.lo >= 0., ReadError::Unsupported)?;
    Ok((sine, cosine, tail.hi))
}
fn world(local: VecP, frame: RawFrame, position: bool, b: &mut Budget<'_>) -> Result<VecP> {
    let zero = P::constant(I::point(0.), b)?;
    let mut out = [zero; 3];
    for (i, world) in out.iter_mut().enumerate() {
        b.charge(256)?;
        let mut value = None;
        for (axis, p) in local.iter().enumerate() {
            let coefficient = frame.axes[axis][i];
            if exact_zero(coefficient) {
                continue;
            }
            let term = scaled(p, coefficient, b)?;
            value = Some(if let Some(v) = value {
                P::add(&v, &term, b)?
            } else {
                term
            });
        }
        let mut value = value.ok_or(ReadError::Unsupported)?;
        if position && !exact_zero(frame.origin[i]) {
            value = value.add(&P::constant(frame.origin[i], b)?, b)?;
        }
        *world = value;
    }
    Ok(out)
}
fn vector_error(a: VecP, c: VecP, b: &mut Budget<'_>) -> Result<f64> {
    let degree = a
        .iter()
        .chain(c.iter())
        .map(P::degree)
        .max()
        .ok_or(ReadError::Invalid)?;
    let zero = P::constant(I::point(0.), b)?;
    let mut residual = [zero; 3];
    for i in 0..3 {
        residual[i] = a[i]
            .add(&c[i].scale(I::point(-1.), b)?, b)?
            .elevate(degree, b)?;
    }
    b.charge((degree + 1) * 128 + 128)?;
    let mut maximum = 0_f64;
    for i in 0..=degree {
        let error = norm(std::array::from_fn(|j| residual[j].controls()[i]));
        check(error.valid(), ReadError::Unsupported)?;
        maximum = maximum.max(error.hi);
    }
    Ok(maximum)
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
    let actual = [
        coordinate_polynomial(target, 0, b)?,
        coordinate_polynomial(target, 1, b)?,
        coordinate_polynomial(target, 2, b)?,
    ];
    let (s, c, tail) = taylor(&theta, order, b)?;
    let expected = world(
        [rho.multiply(&c, b)?, rho.multiply(&s, b)?, z],
        cone.frame,
        true,
        b,
    )?;
    let coefficient = vector_error(expected, actual, b)?;
    // Exact raw orthonormal frame: the two radial scalar errors give this
    // Euclidean bound; no one common mean-value point is assumed.
    let radial_tail = I::point(rh.magnitude()) * I::point(tail);
    let position = I::point(coefficient) + norm([radial_tail, radial_tail, I::point(0.)]);
    check(position.valid(), ReadError::Unsupported)?;
    if position.hi > tol {
        return Ok(Proof::Unsupported);
    }
    let (s, c, dtail) = taylor(&theta, derivative_order, b)?;
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
    let expected_derivative = world([dx, dy, dz], cone.frame, false, b)?;
    let actual_derivative = derivative(actual, b)?;
    let derivative_coefficient = vector_error(expected_derivative, actual_derivative, b)?;
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
    if !lower.is_finite() || lower <= 0. || error.hi >= lower {
        return Ok(Proof::Unsupported);
    }
    b.charge(128)?;
    Ok(Proof::Certified(position.hi))
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
