//! Whole analytic cylinder composition on one earned common parameter cell.
//! Binding original records and finite windows is a separate obligation.
use super::polynomial_trace::Polynomial as P;
use super::*;
use crate::step_read::brep::check;
pub(super) struct CellProof {
    pub(super) position: f64,
    pub(super) derivative: f64,
    pub(super) speed: f64,
}
fn power(x: I, n: usize) -> I {
    (0..n).fold(I::point(1.), |v, _| v * x)
}
fn factorial(n: usize) -> f64 {
    (1..=n).product::<usize>() as f64
}
fn zero(value: I) -> bool {
    value.lo == 0. && value.hi == 0.
}
fn trig_polynomial(angle: &P, order: usize, b: &mut Budget<'_>) -> Result<(P, P, I)> {
    b.charge(4096)?;
    check(order <= 4, ReadError::Unsupported)?;
    let hull = angle.hull();
    check(
        hull.valid() && hull.magnitude() <= 1e6,
        ReadError::Unsupported,
    )?;
    let center = I::point(hull.midpoint());
    let delta = angle.add(&P::constant(-center, b)?, b)?;
    let (s, c) = (center.sin(), center.cos());
    let mut sine = P::constant(I::point(0.), b)?;
    let mut cosine = sine;
    let mut term = P::constant(I::point(1.), b)?;
    for n in 0..=order {
        let (cs, cc) = match n % 4 {
            0 => (s, c),
            1 => (c, -s),
            2 => (-s, -c),
            _ => (-c, s),
        };
        let denominator = I::point(factorial(n));
        sine = sine.add(&term.scale(cs / denominator, b)?, b)?;
        cosine = cosine.add(&term.scale(cc / denominator, b)?, b)?;
        if n < order {
            term = term.multiply(&delta, b)?;
        }
    }
    let remainder =
        power(I::point(delta.hull().magnitude()), order + 1) / I::point(factorial(order + 1));
    check(remainder.valid(), ReadError::Unsupported)?;
    Ok((sine, cosine, remainder))
}
fn vector_error(a: &[P; 3], c: &[P; 3], b: &mut Budget<'_>) -> Result<f64> {
    let n = a
        .iter()
        .chain(c.iter())
        .map(P::degree)
        .max()
        .ok_or(ReadError::Unsupported)?;
    let mut difference = [P::constant(I::point(0.), b)?; 3];
    for i in 0..3 {
        difference[i] = a[i].add(&c[i].scale(I::point(-1.), b)?, b)?.elevate(n, b)?;
    }
    let mut error = 0f64;
    for j in 0..=n {
        b.charge(128)?;
        let e = norm(std::array::from_fn(|i| difference[i].controls()[j]));
        check(e.valid(), ReadError::Unsupported)?;
        error = error.max(e.hi);
    }
    Ok(error)
}
pub(super) fn cell(
    frame: &RawFrame,
    radius: I,
    uv: &[P; 2],
    target: &[P; 3],
    b: &mut Budget<'_>,
) -> Result<CellProof> {
    b.charge(4096)?;
    check(radius.valid() && radius.lo > 0., ReadError::Unsupported)?;
    check(
        uv.iter().all(|p| p.degree() <= 3) && target.iter().all(|p| p.degree() <= 4),
        ReadError::Unsupported,
    )?;
    let (sine, cosine, tail) = trig_polynomial(&uv[0], 2, b)?;
    let mut image = [P::constant(I::point(0.), b)?; 3];
    for (i, value) in image.iter_mut().enumerate() {
        *value = P::constant(frame.origin[i], b)?;
        if !zero(frame.axes[0][i]) {
            *value = value.add(&cosine.scale(radius * frame.axes[0][i], b)?, b)?;
        }
        if !zero(frame.axes[1][i]) {
            *value = value.add(&sine.scale(radius * frame.axes[1][i], b)?, b)?;
        }
        if !zero(frame.axes[2][i]) {
            *value = value.add(&uv[1].scale(frame.axes[2][i], b)?, b)?;
        }
    }
    let axes = norm(frame.axes[0]) + norm(frame.axes[1]);
    let position = I::point(vector_error(&image, target, b)?) + radius * axes * tail;
    // Independent analytic FIRST derivative approximation. No derivative of
    // the position remainder is taken.
    let (sine, cosine, tail) = trig_polynomial(&uv[0], 1, b)?;
    let du = uv[0].derivative(b)?;
    let dv = uv[1].derivative(b)?;
    let sine = du.multiply(&sine, b)?;
    let cosine = du.multiply(&cosine, b)?;
    let mut derivative = [P::constant(I::point(0.), b)?; 3];
    let mut actual = [P::constant(I::point(0.), b)?; 3];
    for i in 0..3 {
        if !zero(frame.axes[0][i]) {
            derivative[i] = derivative[i].add(&sine.scale(-radius * frame.axes[0][i], b)?, b)?;
        }
        if !zero(frame.axes[1][i]) {
            derivative[i] = derivative[i].add(&cosine.scale(radius * frame.axes[1][i], b)?, b)?;
        }
        if !zero(frame.axes[2][i]) {
            derivative[i] = derivative[i].add(&dv.scale(frame.axes[2][i], b)?, b)?;
        }
        actual[i] = target[i].derivative(b)?;
    }
    let first = I::point(vector_error(&derivative, &actual, b)?)
        + radius * axes * tail * I::point(du.hull().magnitude());
    let speed = norm(actual.map(|p| p.hull()));
    check(
        position.valid() && first.valid() && speed.valid(),
        ReadError::Unsupported,
    )?;
    b.charge(64)?;
    Ok(CellProof {
        position: position.hi,
        derivative: first.hi,
        speed: speed.lo,
    })
}

mod binding;
pub(super) use binding::certify;

#[cfg(test)]
mod tests;
