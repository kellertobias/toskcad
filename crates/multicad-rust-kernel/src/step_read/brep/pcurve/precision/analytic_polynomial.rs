//! ORIGINAL polynomial UV traces of degree<=7 on an actual circular cylinder
//! or axial cone, proved in the directed Bernstein algebra through degree80.
//! Each complete common cell earns a whole-cell POSITION Taylor remainder and
//! an INDEPENDENT analytic first-derivative remainder; neither is obtained by
//! differentiating the other. A cell that cannot earn both is bisected
//! exactly (dyadic, in both curves' shared chart) up to a fixed depth. No
//! sample or proposal grants anything, and nothing is ever repaired.
use super::bernstein80::{self, Poly};
use super::cylinder_polynomial::CellProof;
use super::polynomial_trace::Polynomial as P;
use super::*;
use crate::step_read::brep::check;
/// Highest UV or 3D Bezier degree admitted here; the legacy Polynomial16
/// routes keep their own lower guards and unchanged work.
pub(super) const MAX_DEGREE: usize = 7;
const MAX_ORDER: usize = 12;
const MAX_DEPTH: usize = 6;
/// Cylinder: `slope` is exactly zero. Cone: `radius + slope*v` with
/// `slope=tan(semi-angle)`, the STEP/ISO 10303-42 parameterisation.
pub(super) struct Support {
    pub(super) frame: RawFrame,
    pub(super) radius: I,
    pub(super) slope: I,
}
fn exact_zero(v: I) -> bool {
    v.lo == 0. && v.hi == 0.
}
fn hull(p: Poly) -> I {
    let c = &p.c[..=p.degree];
    I::new(
        c.iter().map(|x| x.lo).fold(f64::INFINITY, f64::min),
        c.iter().map(|x| x.hi).fold(f64::NEG_INFINITY, f64::max),
    )
}
fn lift(p: &P) -> Result<Poly> {
    let mut q = Poly::zero(p.degree())?;
    q.c[..=p.degree()].copy_from_slice(p.controls());
    check(q.valid(), ReadError::Unsupported)?;
    Ok(q)
}
fn factorial(n: usize) -> I {
    // n<=13 so n! < 2^53 and every partial product is an exact f64.
    I::point((1..=n).fold(1., |a, i| a * i as f64))
}
/// Taylor polynomials of sin/cos(theta) about the hull midpoint, evaluated
/// by Horner in Bernstein form, and the whole-field Lagrange tail
/// `r^(n+1)/(n+1)!` (all sin/cos derivatives have magnitude<=1).
fn trig(theta: Poly, order: usize, b: &mut Budget<'_>) -> Result<(Poly, Poly, I)> {
    b.charge(4096)?;
    check(order <= MAX_ORDER, ReadError::Unsupported)?;
    let h = hull(theta);
    check(h.valid() && h.magnitude() <= 1e6, ReadError::Unsupported)?;
    let center = I::point(h.midpoint());
    let delta = theta.add(Poly::constant(theta.degree, center, b)?, true, b)?;
    let (s, c) = (center.sin(), center.cos());
    check(s.valid() && c.valid(), ReadError::Unsupported)?;
    let coefficient = |k: usize, sine: bool| -> I {
        let value = match (k % 4, sine) {
            (0, true) | (3, false) => s,
            (1, true) | (0, false) => c,
            (2, true) | (1, false) => -s,
            _ => -c,
        };
        value / factorial(k)
    };
    let mut sine = Poly::constant(0, coefficient(order, true), b)?;
    let mut cosine = Poly::constant(0, coefficient(order, false), b)?;
    for k in (0..order).rev() {
        sine = bernstein80::product(sine, delta, b)?;
        sine = sine.add(
            Poly::constant(sine.degree, coefficient(k, true), b)?,
            false,
            b,
        )?;
        cosine = bernstein80::product(cosine, delta, b)?;
        cosine = cosine.add(
            Poly::constant(cosine.degree, coefficient(k, false), b)?,
            false,
            b,
        )?;
    }
    let r = I::point(hull(delta).magnitude());
    let tail = (0..=order).fold(I::point(1.), |v, _| v * r) / factorial(order + 1);
    check(tail.valid() && tail.lo >= 0., ReadError::Unsupported)?;
    Ok((sine, cosine, tail))
}
fn same_degree(p: [Poly; 3], b: &mut Budget<'_>) -> Result<[Poly; 3]> {
    let n = p.iter().map(|q| q.degree).max().unwrap_or(0);
    Ok([
        p[0].elevate(n, b)?,
        p[1].elevate(n, b)?,
        p[2].elevate(n, b)?,
    ])
}
/// `origin*(position) + sum_k axes[k] * local[k]` in the RAW file frame.
fn world(
    local: [Poly; 3],
    frame: &RawFrame,
    position: bool,
    b: &mut Budget<'_>,
) -> Result<[Poly; 3]> {
    let local = same_degree(local, b)?;
    let n = local[0].degree;
    let mut out = [Poly::zero(n)?; 3];
    for (i, value) in out.iter_mut().enumerate() {
        b.charge(256)?;
        if position && !exact_zero(frame.origin[i]) {
            *value = Poly::constant(n, frame.origin[i], b)?;
        }
        for (axis, p) in local.iter().enumerate() {
            let coefficient = frame.axes[axis][i];
            if !exact_zero(coefficient) {
                *value = value.add(p.scale(coefficient, b)?, false, b)?;
            }
        }
    }
    Ok(out)
}
/// Convex-hull bound of `sup ||a-c||` over the whole cell, and a LOWER
/// bound of `||a-c||` at the two cell ends. Bernstein end coefficients are
/// the exact end values, so that lower bound needs no hull argument.
fn residual_bounds(a: [Poly; 3], c: [Poly; 3], b: &mut Budget<'_>) -> Result<(f64, f64)> {
    let n = a.iter().chain(&c).map(|p| p.degree).max().unwrap_or(0);
    let mut residual = [Poly::zero(n)?; 3];
    for (i, r) in residual.iter_mut().enumerate() {
        *r = a[i].elevate(n, b)?.add(c[i].elevate(n, b)?, true, b)?;
    }
    b.charge((n + 1) * 64 + 128)?;
    let mut maximum = 0_f64;
    for j in 0..=n {
        let e = norm(std::array::from_fn(|i| residual[i].c[j]));
        check(e.valid(), ReadError::Unsupported)?;
        maximum = maximum.max(e.hi);
    }
    let ends = [0, n].map(|j| norm(std::array::from_fn(|i| residual[i].c[j])));
    check(ends.iter().all(|e| e.valid()), ReadError::Unsupported)?;
    Ok((maximum, ends[0].lo.max(ends[1].lo)))
}
fn vector_error(a: [Poly; 3], c: [Poly; 3], b: &mut Budget<'_>) -> Result<f64> {
    Ok(residual_bounds(a, c, b)?.0)
}
fn polynomial_norm_hull(p: [Poly; 3]) -> I {
    norm(p.map(hull))
}
enum Cell {
    /// This cell's Taylor order would exceed the degree-80 algebra or
    /// `MAX_ORDER`; a smaller piece may still be proved.
    Undecided,
    /// At a cell END the actual deviation already exceeds `tol`, so no
    /// subdivision can certify the pair. This only refuses; it never
    /// reports a disagreement or grants anything.
    Exceeds,
    Proved(CellProof),
}
fn proof(
    support: &Support,
    uv: [Poly; 2],
    target: [Poly; 3],
    tol: f64,
    b: &mut Budget<'_>,
) -> Result<Cell> {
    b.charge(8192)?;
    let [theta, v] = uv;
    let rho =
        Poly::constant(v.degree, support.radius, b)?.add(v.scale(support.slope, b)?, false, b)?;
    let rho_hull = hull(rho);
    check(rho_hull.valid() && rho_hull.lo > 0., ReadError::Unsupported)?;
    let axes = norm(support.frame.axes[0]) + norm(support.frame.axes[1]);
    check(axes.valid(), ReadError::Unsupported)?;
    let r = I::point(
        hull(theta.add(
            Poly::constant(theta.degree, I::point(hull(theta).midpoint()), b)?,
            true,
            b,
        )?)
        .magnitude(),
    );
    // Smallest order whose WHOLE-FIELD position tail is at most tol/16 and
    // whose image degree fits the algebra. Exact and data independent.
    let mut order = None;
    let mut power = r;
    for n in 0..=MAX_ORDER {
        b.charge(64)?;
        let tail = power / factorial(n + 1) * I::point(rho_hull.magnitude()) * axes;
        let degree = rho.degree + theta.degree * n.max(1) + theta.degree;
        if degree <= 80 && tail.valid() && tail.hi <= tol / 16. {
            order = Some(n.max(1));
            break;
        }
        power = power * r;
    }
    let Some(order) = order else {
        return Ok(Cell::Undecided);
    };
    let (sine, cosine, tail) = trig(theta, order, b)?;
    let image = world(
        [
            bernstein80::product(rho, cosine, b)?,
            bernstein80::product(rho, sine, b)?,
            v,
        ],
        &support.frame,
        true,
        b,
    )?;
    let (upper, end) = residual_bounds(image, target, b)?;
    let taylor = I::point(rho_hull.magnitude()) * axes * tail;
    let position = I::point(upper) + taylor;
    let actual_end = I::point(end) - taylor;
    check(actual_end.valid(), ReadError::Unsupported)?;
    if actual_end.lo > tol {
        return Ok(Cell::Exceeds);
    }
    // Independent analytic FIRST derivative of the ORIGINAL composition:
    // d/dt(rho*cos) = rho'*cos - rho*theta'*sin, d/dt(rho*sin) =
    // rho'*sin + rho*theta'*cos, d/dt(v) = v'. Its own trig remainder
    // multiplies the actual amplitudes, never a differentiated tail.
    let rho_derivative = rho.derivative(b)?;
    let angular = bernstein80::product(rho, theta.derivative(b)?, b)?;
    let dx = bernstein80::product(rho_derivative, cosine, b)?;
    let dx = same_degree([dx, bernstein80::product(angular, sine, b)?, dx], b)?;
    let dx = dx[0].add(dx[1], true, b)?;
    let dy = bernstein80::product(rho_derivative, sine, b)?;
    let dy = same_degree([dy, bernstein80::product(angular, cosine, b)?, dy], b)?;
    let dy = dy[0].add(dy[1], false, b)?;
    let expected = world([dx, dy, v.derivative(b)?], &support.frame, false, b)?;
    let actual = [
        target[0].derivative(b)?,
        target[1].derivative(b)?,
        target[2].derivative(b)?,
    ];
    let amplitude =
        I::point(hull(rho_derivative).magnitude()) + I::point(hull(angular).magnitude());
    let derivative = I::point(vector_error(expected, actual, b)?) + amplitude * axes * tail;
    let speed = polynomial_norm_hull(actual);
    check(
        position.valid() && derivative.valid() && speed.valid(),
        ReadError::Unsupported,
    )?;
    b.charge(64)?;
    Ok(Cell::Proved(CellProof {
        position: position.hi,
        derivative: derivative.hi,
        speed: speed.lo,
    }))
}
/// Directed de Casteljau value of a Bernstein polynomial at `t`.
fn value(p: &P, t: I, b: &mut Budget<'_>) -> Result<I> {
    let n = p.degree();
    b.charge((n + 1) * (n + 1) * 8 + 64)?;
    let mut row = [I::point(0.); MAX_DEGREE + 1];
    row[..=n].copy_from_slice(p.controls());
    for level in 1..=n {
        for i in 0..=n - level {
            row[i] = (I::point(1.) - t) * row[i] + t * row[i + 1];
        }
    }
    check(row[0].valid(), ReadError::Unsupported)?;
    Ok(row[0])
}
/// Points scanned per common cell before any whole-cell proof is built.
const SCAN: usize = 32;
/// Cheap REFUSAL scan: a lower bound of the actual deviation at `SCAN`
/// interior points of the cell, using directed interval sin/cos of the
/// actual angle. `false` means some point is already farther than `tol`,
/// so no proof can certify this pair. `true` grants nothing.
pub(super) fn scan(
    support: &Support,
    uv: &[P; 2],
    target: &[P; 3],
    tol: f64,
    b: &mut Budget<'_>,
) -> Result<bool> {
    check(
        uv.iter().chain(target).all(|p| p.degree() <= MAX_DEGREE),
        ReadError::Unsupported,
    )?;
    for j in 0..SCAN {
        b.charge(512)?;
        let t = I::point((2 * j + 1) as f64) / I::point((2 * SCAN) as f64);
        let theta = value(&uv[0], t, b)?;
        let v = value(&uv[1], t, b)?;
        let rho = support.radius + support.slope * v;
        let local = [rho * theta.cos(), rho * theta.sin(), v];
        let image: [I; 3] = std::array::from_fn(|i| {
            support.frame.origin[i]
                + (0..3).fold(I::point(0.), |a, k| a + support.frame.axes[k][i] * local[k])
        });
        let actual = [
            value(&target[0], t, b)?,
            value(&target[1], t, b)?,
            value(&target[2], t, b)?,
        ];
        let gap = norm(std::array::from_fn(|i| image[i] - actual[i]));
        check(gap.valid(), ReadError::Unsupported)?;
        if gap.lo > tol {
            return Ok(false);
        }
    }
    Ok(true)
}
fn halves(p: Poly, b: &mut Budget<'_>) -> Result<[Poly; 2]> {
    p.split(b)
}
/// Prove one common cell, bisecting exactly at local 1/2 in BOTH curves
/// when `accept` rejects a sub-cell. Every reported proof is expressed in
/// the caller's cell-local parameter: a dyadic piece of width `2^-k` has
/// derivative and speed bounds multiplied by `2^k` (exact).
pub(super) fn certify(
    support: &Support,
    uv: &[P; 2],
    target: &[P; 3],
    tol: f64,
    accept: &mut dyn FnMut(&CellProof, &mut Budget<'_>) -> Result<bool>,
    b: &mut Budget<'_>,
) -> Result<bool> {
    check(
        uv.iter().chain(target).all(|p| p.degree() <= MAX_DEGREE),
        ReadError::Unsupported,
    )?;
    let uv = [lift(&uv[0])?, lift(&uv[1])?];
    let target = [lift(&target[0])?, lift(&target[1])?, lift(&target[2])?];
    piece(support, uv, target, tol, 0, accept, b)
}
fn piece(
    support: &Support,
    uv: [Poly; 2],
    target: [Poly; 3],
    tol: f64,
    depth: usize,
    accept: &mut dyn FnMut(&CellProof, &mut Budget<'_>) -> Result<bool>,
    b: &mut Budget<'_>,
) -> Result<bool> {
    b.charge(1024)?;
    let local = match proof(support, uv, target, tol, b)? {
        Cell::Exceeds => return Ok(false),
        Cell::Undecided => None,
        Cell::Proved(local) => Some(local),
    };
    if let Some(local) = local {
        let scale = I::point((1_u64 << depth) as f64);
        let scaled = CellProof {
            position: local.position,
            derivative: (I::point(local.derivative) * scale).hi,
            speed: (I::point(local.speed) * scale).lo,
        };
        if accept(&scaled, b)? {
            return Ok(true);
        }
    }
    if depth == MAX_DEPTH {
        return Ok(false);
    }
    let [u0, u1] = halves(uv[0], b)?;
    let [v0, v1] = halves(uv[1], b)?;
    let [x0, x1] = halves(target[0], b)?;
    let [y0, y1] = halves(target[1], b)?;
    let [z0, z1] = halves(target[2], b)?;
    Ok(
        piece(support, [u0, v0], [x0, y0, z0], tol, depth + 1, accept, b)?
            && piece(support, [u1, v1], [x1, y1, z1], tol, depth + 1, accept, b)?,
    )
}

#[cfg(test)]
mod tests;
