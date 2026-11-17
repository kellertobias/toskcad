//! Direct ORIGINAL-REAL bilinear surface/finite rational trace composition.
//! A chart is a proposal until every homogeneous coefficient and orientation
//! polynomial earns its bound; no sampled point publishes agreement.
use super::super::Proof;
use super::*;
use bezier::{Coordinates, CurveCells};
const N: usize = 65;
const SCRATCH: usize = 2 * 1024 * 1024;
#[derive(Clone, Copy)]
struct Poly {
    c: [I; N],
    degree: usize,
}
impl Poly {
    fn zero(degree: usize) -> Result<Self> {
        if degree >= N {
            return Err(ReadError::Unsupported);
        }
        Ok(Self {
            c: [I::point(0.); N],
            degree,
        })
    }
    fn from_controls(values: &[[I; 4]], axis: usize) -> Result<Self> {
        let mut p = Self::zero(values.len().checked_sub(1).ok_or(ReadError::Invalid)?)?;
        for (i, v) in values.iter().enumerate() {
            p.c[i] = v[axis];
        }
        Ok(p)
    }
    fn valid(self) -> bool {
        self.c[..=self.degree].iter().all(|x| x.valid())
    }
    fn scale(self, value: I, budget: &mut Budget<'_>) -> Result<Self> {
        budget.charge((self.degree + 1) * 16 + 512)?;
        let mut p = self;
        for c in &mut p.c[..=p.degree] {
            *c = *c * value;
        }
        if p.valid() {
            Ok(p)
        } else {
            Err(ReadError::Unsupported)
        }
    }
    fn add(self, other: Self, subtract: bool, budget: &mut Budget<'_>) -> Result<Self> {
        if self.degree != other.degree {
            return Err(ReadError::Unsupported);
        }
        budget.charge((self.degree + 1) * 16 + 512)?;
        let mut p = self;
        for (i, c) in p.c[..=p.degree].iter_mut().enumerate() {
            *c = if subtract {
                *c - other.c[i]
            } else {
                *c + other.c[i]
            };
        }
        if p.valid() {
            Ok(p)
        } else {
            Err(ReadError::Unsupported)
        }
    }
    fn derivative(self, budget: &mut Budget<'_>) -> Result<Self> {
        budget.charge(512)?;
        if self.degree == 0 {
            return Self::zero(0);
        }
        budget.charge(self.degree * 24 + 64)?;
        let mut p = Self::zero(self.degree - 1)?;
        for (i, c) in p.c[..=p.degree].iter_mut().enumerate() {
            *c = I::point(self.degree as f64) * (self.c[i + 1] - self.c[i]);
        }
        if p.valid() {
            Ok(p)
        } else {
            Err(ReadError::Unsupported)
        }
    }
    fn split(self, budget: &mut Budget<'_>) -> Result<[Self; 2]> {
        budget.charge(self.degree * (self.degree + 1) / 2 * 16 + (self.degree + 1) * 32 + 512)?;
        let (mut a, mut b) = (Self::zero(self.degree)?, Self::zero(self.degree)?);
        let mut row = self.c;
        a.c[0] = row[0];
        b.c[self.degree] = row[self.degree];
        for level in 1..=self.degree {
            for i in 0..=self.degree - level {
                row[i] = (row[i] + row[i + 1]) * I::point(0.5);
            }
            a.c[level] = row[0];
            b.c[self.degree - level] = row[self.degree - level];
        }
        if a.valid() && b.valid() {
            Ok([a, b])
        } else {
            Err(ReadError::Unsupported)
        }
    }
}
fn binomial_row(n: usize, budget: &mut Budget<'_>) -> Result<[I; N]> {
    if n >= N {
        return Err(ReadError::Unsupported);
    }
    budget.charge((n + 1) * 32 + 512)?;
    let mut row = [I::point(0.); N];
    let mut value = 1_u128;
    for (i, slot) in row[..=n].iter_mut().enumerate() {
        if i > 0 {
            value = value
                .checked_mul((n - i + 1) as u128)
                .ok_or(ReadError::Unsupported)?
                / (i as u128);
        }
        let f = value as f64;
        *slot = if f as u128 == value {
            I::point(f)
        } else {
            I::new(f.next_down(), f.next_up())
        };
    }
    Ok(row)
}
fn product(a: Poly, b: Poly, budget: &mut Budget<'_>) -> Result<Poly> {
    let degree = a
        .degree
        .checked_add(b.degree)
        .ok_or(ReadError::Unsupported)?;
    budget.charge((a.degree + 1) * (b.degree + 1) * 16 + (degree + 1) * 8 + 512)?;
    let mut p = Poly::zero(degree)?;
    let ar = binomial_row(a.degree, budget)?;
    let br = binomial_row(b.degree, budget)?;
    let cr = binomial_row(degree, budget)?;
    for i in 0..=a.degree {
        for j in 0..=b.degree {
            let weight = ar[i] * br[j] / cr[i + j];
            p.c[i + j] = p.c[i + j] + weight * a.c[i] * b.c[j];
        }
    }
    if p.valid() {
        Ok(p)
    } else {
        Err(ReadError::Unsupported)
    }
}
fn dot_polynomial(a: [Poly; 3], b: [Poly; 3], budget: &mut Budget<'_>) -> Result<Poly> {
    product(a[0], b[0], budget)?
        .add(product(a[1], b[1], budget)?, false, budget)?
        .add(product(a[2], b[2], budget)?, false, budget)
}
fn derivative_numerator(h: [Poly; 4], budget: &mut Budget<'_>) -> Result<[Poly; 3]> {
    let dw = h[3].derivative(budget)?;
    let mut out = [Poly::zero(0)?; 3];
    for i in 0..3 {
        out[i] = product(h[i].derivative(budget)?, h[3], budget)?.add(
            product(h[i], dw, budget)?,
            true,
            budget,
        )?;
    }
    Ok(out)
}
fn positive(p: Poly, budget: &mut Budget<'_>) -> Result<bool> {
    budget.charge((p.degree + 1) * 8 + 64)?;
    Ok(p.c[..=p.degree].iter().all(|x| x.valid() && x.lo > 0.))
}
#[derive(Clone, Copy)]
struct Node {
    p: [Poly; 3],
    depth: usize,
}
fn vector_upper(p: [Poly; 3], budget: &mut Budget<'_>) -> Result<f64> {
    budget.charge((p[0].degree + 1) * 64 + 64)?;
    let mut value = 0_f64;
    for i in 0..=p[0].degree {
        let n = norm([p[0].c[i], p[1].c[i], p[2].c[i]]);
        if !n.valid() {
            return Err(ReadError::Unsupported);
        }
        value = value.max(n.hi);
    }
    Ok(value)
}
fn bound(mut p: [Poly; 3], denominator: f64, allowed: f64, budget: &mut Budget<'_>) -> Result<f64> {
    if !denominator.is_finite() || denominator <= 0. || !allowed.is_finite() || allowed <= 0. {
        return Err(ReadError::Unsupported);
    }
    // Fixed depth<=24 has at most25 pending siblings. The caller preflights
    // the entire 2MiB numerical stack before entering this implementation.
    let mut stack: [Option<Node>; 26] = std::array::from_fn(|_| None);
    let (mut len, mut depth, mut maximum) = (0_usize, 0_usize, 0_f64);
    loop {
        budget.charge(64)?;
        let upper = (I::point(vector_upper(p, budget)?) / I::point(denominator)).hi;
        if !upper.is_finite() {
            return Err(ReadError::Unsupported);
        }
        if upper <= allowed {
            maximum = maximum.max(upper);
            if len == 0 {
                return Ok(maximum);
            }
            len -= 1;
            let node = stack[len].take().ok_or(ReadError::Invalid)?;
            p = node.p;
            depth = node.depth;
            continue;
        }
        if depth == 24 {
            return Err(ReadError::Unsupported);
        }
        let x = p[0].split(budget)?;
        let y = p[1].split(budget)?;
        let z = p[2].split(budget)?;
        if len >= stack.len() {
            return Err(ReadError::Limit);
        }
        stack[len] = Some(Node {
            p: [x[1], y[1], z[1]],
            depth: depth + 1,
        });
        len += 1;
        p = [x[0], y[0], z[0]];
        depth += 1;
    }
}
struct Bilinear {
    h: [[I; 4]; 4],
    minimum_weight: f64,
    min: [BigRational; 2],
    max: [BigRational; 2],
}
fn source_surface(
    source: &Source<'_>,
    id: u32,
    surface: &crate::brep::NurbsSurface,
    budget: &mut Budget<'_>,
) -> Result<Bilinear> {
    budget.charge(4096)?;
    let f = fields(source, id, 2, budget)?;
    let mut min = std::array::from_fn(|_| BigRational::from_integer(0.into()));
    let mut max = std::array::from_fn(|_| BigRational::from_integer(1.into()));
    for axis in 0..2 {
        integer(source, id, f.degrees[axis], 1, budget)?;
        if let Some((m, k)) = f.knots[axis] {
            let (m, k) = (list(m)?, list(k)?);
            if m.len() != 2 || k.len() != 2 {
                return Err(ReadError::Unsupported);
            }
            for v in m {
                integer(source, id, v, 2, budget)?;
            }
            min[axis] = rational(source, id, &k[0], budget)?;
            max[axis] = rational(source, id, &k[1], budget)?;
        }
        let bits = min[axis]
            .numer()
            .bits()
            .max(min[axis].denom().bits())
            .max(max[axis].numer().bits())
            .max(max[axis].denom().bits());
        if bits > 4096 {
            return Err(ReadError::Unsupported);
        }
        budget.charge(bigint_work(
            (bits as usize).saturating_mul(2).saturating_add(2),
        ))?;
        if min[axis] >= max[axis] {
            return Err(ReadError::Unsupported);
        }
    }
    if surface.degrees() != [1, 1] || surface.counts() != [2, 2] {
        return Err(ReadError::Unsupported);
    }
    let rows = list(f.controls)?;
    if rows.len() != 2 {
        return Err(ReadError::Unsupported);
    }
    let mut ids = [0_u32; 4];
    for i in 0..2 {
        let row = list(&rows[i])?;
        if row.len() != 2 {
            return Err(ReadError::Unsupported);
        }
        for j in 0..2 {
            ids[2 * i + j] = reference(&row[j])?;
        }
    }
    let mut weights = [I::point(1.); 4];
    if let Some(w) = f.weights {
        let rows = list(w)?;
        if rows.len() != 2 {
            return Err(ReadError::Unsupported);
        }
        for i in 0..2 {
            let row = list(&rows[i])?;
            if row.len() != 2 {
                return Err(ReadError::Unsupported);
            }
            for j in 0..2 {
                weights[2 * i + j] = scalar(source, id, &row[j], budget)?;
            }
        }
    }
    if weights.iter().any(|x| !x.valid() || x.lo <= 0.) {
        return Err(ReadError::Unsupported);
    }
    let maximum = weights.iter().map(|x| x.hi).fold(0_f64, f64::max);
    let mut h = [[I::point(0.); 4]; 4];
    let mut minimum = f64::INFINITY;
    for i in 0..4 {
        budget.charge(128)?;
        let p = point(source, ids[i], 3, budget)?;
        let w = weights[i] / I::point(maximum);
        minimum = minimum.min(w.lo);
        h[i] = [p[0] * w, p[1] * w, p[2] * w, w];
        if h[i].into_iter().any(|x| !x.valid()) || w.lo <= 0. {
            return Err(ReadError::Unsupported);
        }
    }
    Ok(Bilinear {
        h,
        minimum_weight: minimum,
        min,
        max,
    })
}
fn composition(
    uv: &[[I; 4]],
    surface: &Bilinear,
    budget: &mut Budget<'_>,
) -> Result<([Poly; 4], f64)> {
    budget.charge(4096)?;
    let u = Poly::from_controls(uv, 0)?;
    let v = Poly::from_controls(uv, 1)?;
    let w = Poly::from_controls(uv, 3)?;
    let uc = w.add(u, true, budget)?;
    let vc = w.add(v, true, budget)?;
    let basis = [
        product(uc, vc, budget)?,
        product(uc, v, budget)?,
        product(u, vc, budget)?,
        product(u, v, budget)?,
    ];
    let mut h = [Poly::zero(2 * u.degree)?; 4];
    for (axis, row) in h.iter_mut().enumerate() {
        for (i, b) in basis.iter().enumerate() {
            *row = row.add(b.scale(surface.h[i][axis], budget)?, false, budget)?;
        }
    }
    // UnitBox2 exact membership and positive UV/surface weights prove this
    // lower independently of interval cancellation in W-U or W-V.
    let minimum = uv.iter().map(|p| p[3].lo).fold(f64::INFINITY, f64::min);
    let lower = (I::point(minimum).square() * I::point(surface.minimum_weight)).lo;
    if !lower.is_finite() || lower <= 0. {
        return Err(ReadError::Unsupported);
    }
    Ok((h, lower))
}
fn centered(mut h: [Poly; 4], origin: V, budget: &mut Budget<'_>) -> Result<[Poly; 4]> {
    for i in 0..3 {
        h[i] = h[i].add(h[3].scale(origin[i], budget)?, true, budget)?;
    }
    Ok(h)
}
fn endpoint(h: [I; 4], origin: V, direction: V, den: I) -> Result<I> {
    let p: V = std::array::from_fn(|i| h[i] / h[3] - origin[i]);
    let n = (0..3).fold(I::point(0.), |s, i| s + p[i] * direction[i]) / den;
    if n.valid() {
        Ok(n)
    } else {
        Err(ReadError::Unsupported)
    }
}
fn line_certificate(
    source: &Source<'_>,
    id: u32,
    edge: &Edge,
    curve: &CurveCells<'_, '_, '_>,
    surface: &Bilinear,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    let Curve::Line { origin, direction } = &edge.curve else {
        return Err(ReadError::Unsupported);
    };
    let origin = cached(*origin);
    let direction = cached(*direction);
    let den = (0..3).fold(I::point(0.), |s, i| s + direction[i].square());
    let length = den.sqrt();
    if !den.valid() || den.lo <= 0. || !length.valid() || length.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let raw_id = unwrapped(source, id, budget)?;
    let [_, _, raw_vector] = source.args(raw_id, "LINE")? else {
        return Err(ReadError::Unsupported);
    };
    let raw_direction = vector(source, reference(raw_vector)?, 3, budget)?;
    let target = (I::point(analytic::curve_reserve(source, id, edge, budget)?)
        + I::point(super::super::finite_axis::endpoint_reserve(
            source, id, edge, budget,
        )?))
    .hi;
    if !target.is_finite() || target >= tol {
        return Ok(Proof::Unsupported);
    }
    let cells = curve.cells();
    let first = cells.first().ok_or(ReadError::Unsupported)?;
    let last = cells.last().ok_or(ReadError::Unsupported)?;
    let (first_h, _) = composition(first.controls(), surface, budget)?;
    let (last_h, _) = composition(last.controls(), surface, budget)?;
    let q0 = endpoint(
        std::array::from_fn(|i| first_h[i].c[0]),
        origin,
        direction,
        den,
    )?;
    let q1 = endpoint(
        std::array::from_fn(|i| last_h[i].c[last_h[i].degree]),
        origin,
        direction,
        den,
    )?;
    let (low, high) = (q0.hi, q1.lo);
    if !low.is_finite() || !high.is_finite() || low >= high {
        return Ok(Proof::Unsupported);
    }
    let chosen = edge.range.map(|t| t.clamp(low, high));
    if chosen[0] >= chosen[1] {
        return Ok(Proof::Unsupported);
    }
    let gap = (I::point(chosen[0]) - I::point(edge.range[0]))
        .magnitude()
        .max((I::point(chosen[1]) - I::point(edge.range[1])).magnitude());
    let fixed = (I::point(target) + I::point(gap) * length).hi;
    let allowance = (I::point(tol) - I::point(fixed)).lo;
    if !fixed.is_finite() || allowance <= 0. {
        return Ok(Proof::Unsupported);
    }
    let mut error = 0_f64;
    for cell in cells {
        budget.charge(16_384)?;
        let (h, lower) = composition(cell.controls(), surface, budget)?;
        let centered = centered(h, origin, budget)?;
        // Both reference and ORIGINAL 3D directions must have positive entire
        // projected derivative. This is orientation, never a sampled arrow.
        let q = derivative_numerator(centered, budget)?;
        for d in [direction, raw_direction] {
            let p = q[0]
                .scale(d[0], budget)?
                .add(q[1].scale(d[1], budget)?, false, budget)?
                .add(q[2].scale(d[2], budget)?, false, budget)?;
            if !positive(p, budget)? {
                return Ok(Proof::Unsupported);
            }
        }
        let cross = [
            centered[1].scale(direction[2], budget)?.add(
                centered[2].scale(direction[1], budget)?,
                true,
                budget,
            )?,
            centered[2].scale(direction[0], budget)?.add(
                centered[0].scale(direction[2], budget)?,
                true,
                budget,
            )?,
            centered[0].scale(direction[1], budget)?.add(
                centered[1].scale(direction[0], budget)?,
                true,
                budget,
            )?,
        ];
        let divisor = (I::point(lower) * I::point(length.lo)).lo;
        let upper = bound(cross, divisor, allowance, budget)?;
        error = error.max(upper);
    }
    let total = (I::point(error) + I::point(fixed)).hi;
    Ok(if total.is_finite() && total <= tol {
        Proof::Certified(total)
    } else {
        Proof::Unsupported
    })
}
fn spline_certificate(
    source: &Source<'_>,
    id: u32,
    edge: &Edge,
    curve: &CurveCells<'_, '_, '_>,
    surface: &Bilinear,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    let Curve::Nurbs {
        curve: target,
        closed: false,
    } = &edge.curve
    else {
        return Err(ReadError::Unsupported);
    };
    let target = bezier::original_curve(
        source,
        unwrapped(source, id, budget)?,
        target,
        Coordinates::Raw(3),
        budget,
    )?;
    if !target.matches(source, unwrapped(source, id, budget)?)
        || !matches!(target.coordinates(), Coordinates::Raw(3))
    {
        return Err(ReadError::Invalid);
    }
    if curve.cells().len() != target.cells().len() {
        return Ok(Proof::Unsupported);
    }
    // A separately proposed increasing map joins each ACTUAL source and
    // target Bezier knot cell, regardless of unequal knot spacing. Complete
    // coefficient/error/orientation checks below earn that map's authority.
    let mut maximum = 0_f64;
    for (s, t) in curve.cells().iter().zip(target.cells()) {
        budget.charge(16_384)?;
        let (sh, slower) = composition(s.controls(), surface, budget)?;
        let th = std::array::from_fn(|i| Poly::from_controls(t.controls(), i));
        let [tx, ty, tz, tw] = th;
        let th = [tx?, ty?, tz?, tw?];
        let tlower = t
            .controls()
            .iter()
            .map(|p| p[3].lo)
            .fold(f64::INFINITY, f64::min);
        if !tlower.is_finite() || tlower <= 0. {
            return Err(ReadError::Unsupported);
        }
        let origin =
            std::array::from_fn(|i| I::point((t.controls()[0][i] / t.controls()[0][3]).midpoint()));
        let shc = centered(sh, origin, budget)?;
        let thc = centered(th, origin, budget)?;
        let mut residual = [Poly::zero(0)?; 3];
        for i in 0..3 {
            residual[i] = product(shc[i], th[3], budget)?.add(
                product(thc[i], sh[3], budget)?,
                true,
                budget,
            )?;
        }
        let divisor = (I::point(slower) * I::point(tlower)).lo;
        let upper = bound(residual, divisor, tol, budget)?;
        let dot = dot_polynomial(
            derivative_numerator(shc, budget)?,
            derivative_numerator(thc, budget)?,
            budget,
        )?;
        if !positive(dot, budget)? {
            return Ok(Proof::Unsupported);
        }
        maximum = maximum.max(upper);
    }
    // Both ORIGINAL natural curve domains were proved in full. Any actual
    // finite topological restriction is a subset; no cached endpoint oracle
    // or equal-natural-parameter assertion is used here.
    Ok(Proof::Certified(maximum))
}
#[allow(
    clippy::too_many_arguments,
    reason = "Original owner, record identities, actual finite curves and resource policy are independent inputs."
)]
pub(super) fn certify(
    source: &Source<'_>,
    id: u32,
    definition: &super::super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    let (
        Surface::Nurbs { surface, closed },
        Curve::Nurbs {
            curve,
            closed: false,
        },
    ) = (surface, &definition.geometry)
    else {
        return Ok(None);
    };
    if surface.degrees() != [1, 1] || surface.counts() != [2, 2] || closed.iter().any(|x| *x) {
        return Ok(None);
    }
    if curve.degree() > 11 || curve.points().len() > 1024 {
        return Ok(None);
    }
    match &edge.curve {
        Curve::Line { .. } => {}
        Curve::Nurbs {
            curve,
            closed: false,
        } if curve.degree() <= 11 && curve.points().len() <= 1024 => {}
        _ => return Ok(None),
    }
    // Fresh decoded arguments are bounded above by4MiB in this class. Two
    // 8MiB extraction envelopes (the earlier result remains resident), the
    // complete2MiB numerical stack, domain rationals and the caller's256KiB
    // cache all fit inside the unchanged64MiB fresh checking workspace.
    const ENVELOPE: usize = 2 * 8 * 1024 * 1024 + SCRATCH + 4 * 1024 * 1024 + 256 * 1024;
    if ENVELOPE > 64 * 1024 * 1024 {
        return Err(ReadError::Limit);
    }
    budget.charge(4096)?;
    let surface = source_surface(source, definition.surface, surface, budget)?;
    let trace = bezier::original_curve(
        source,
        definition.curve_id,
        curve,
        Coordinates::UnitBox2 {
            min: &surface.min,
            max: &surface.max,
        },
        budget,
    )?;
    if !trace.matches(source, definition.curve_id) {
        return Err(ReadError::Invalid);
    }
    match trace.coordinates() {
        Coordinates::UnitBox2 { min, max }
            if std::ptr::eq(min, &surface.min) && std::ptr::eq(max, &surface.max) => {}
        _ => return Err(ReadError::Invalid),
    }
    if trace
        .peak_bytes()
        .saturating_add(SCRATCH)
        .saturating_add(4 * 1024 * 1024)
        > 64 * 1024 * 1024
        || trace.retained_bytes() > 8 * 1024 * 1024
    {
        return Err(ReadError::Limit);
    }
    let result = match edge.curve {
        Curve::Line { .. } => line_certificate(source, id, edge, &trace, &surface, tol, budget),
        _ => spline_certificate(source, id, edge, &trace, &surface, tol, budget),
    };
    budget.charge(64)?;
    match result {
        Ok(p) => Ok(Some(p)),
        Err(ReadError::Unsupported) => Ok(Some(Proof::Unsupported)),
        Err(e) => Err(e),
    }
}
#[cfg(test)]
mod oracle_tests;
#[cfg(test)]
mod tests;
