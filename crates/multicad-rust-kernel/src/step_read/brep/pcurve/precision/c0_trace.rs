//! Straight ISO UV traces on an ORIGINAL polynomial tensor surface whose
//! interior knots may be only C0 (multiplicity equal to the degree), which
//! the C1 tensor-trace proof refuses by contract.
//!
//! An iso line at a clamped end or at a knot of multiplicity equal to the
//! degree is EXACTLY one control row of the original net: a 1D polynomial
//! B-spline in the other parameter. Each cell between exact rational breaks
//! (target knots, the trace's knot-line crossings) lies in one span of that
//! row; blossoming gives its exact Bezier controls, and the Bernstein
//! residual against the restricted target bounds position and tangent
//! directly. The target is an original polynomial B-spline on the SAME
//! parameter, or a CIRCLE whose angle is the trace parameter (proved with
//! the analytic support proof). Nothing is sampled or repaired.
use super::super::Proof;
use super::analytic_polynomial::{self, Support};
use super::bezier::{self, Coordinates, CurveCells};
use super::cylinder_polynomial::{CellProof, compare, operation};
use super::polynomial_trace::Polynomial as P;
use super::*;
use crate::step_read::brep::check;
const MAX_DEGREE: usize = 11;
const MAX_KNOTS: usize = 1025;
const MAX_DEPTH: usize = 40;
const MAX_ROW: usize = 1024;
type V3 = [I; 3];
type Q2 = [BigRational; 2];

/// Distinct, strictly increasing original knot rationals of one axis.
fn distinct_knots(
    source: &Source<'_>,
    id: u32,
    dimension: usize,
    axis: usize,
    b: &mut Budget<'_>,
) -> Result<Vec<BigRational>> {
    let field = fields(source, id, dimension, b)?;
    let Some((_, values)) = field.knots[axis] else {
        return Err(ReadError::Unsupported);
    };
    let values = list(values)?;
    check(
        (2..=MAX_KNOTS).contains(&values.len()),
        ReadError::Unsupported,
    )?;
    let mut out = Vec::new();
    out.try_reserve_exact(values.len())
        .map_err(|_| ReadError::Limit)?;
    for v in values {
        b.charge(128)?;
        let q = rational(source, id, v, b)?;
        if let Some(previous) = out.last() {
            check(compare(previous, &q, b)?.is_lt(), ReadError::Unsupported)?;
        }
        out.push(q);
    }
    Ok(out)
}
/// Sorted union of two knot lists over the SAME exact domain.
fn merge(a: &[BigRational], c: &[BigRational], b: &mut Budget<'_>) -> Result<Vec<BigRational>> {
    let (Some(a_end), Some(c_end)) = (a.last(), c.last()) else {
        return Err(ReadError::Unsupported);
    };
    check(
        compare(&a[0], &c[0], b)?.is_eq() && compare(a_end, c_end, b)?.is_eq(),
        ReadError::Unsupported,
    )?;
    let mut out = Vec::new();
    out.try_reserve_exact(a.len() + c.len())
        .map_err(|_| ReadError::Limit)?;
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < c.len() {
        b.charge(256)?;
        let next = match (a.get(i), c.get(j)) {
            (Some(x), Some(y)) => match compare(x, y, b)? {
                std::cmp::Ordering::Less => {
                    i += 1;
                    x
                }
                std::cmp::Ordering::Greater => {
                    j += 1;
                    y
                }
                std::cmp::Ordering::Equal => {
                    i += 1;
                    j += 1;
                    x
                }
            },
            (Some(x), None) => {
                i += 1;
                x
            }
            (None, Some(y)) => {
                j += 1;
                y
            }
            (None, None) => break,
        };
        out.push(next.clone());
    }
    Ok(out)
}
fn casteljau(points: &mut [V3], t: I) -> V3 {
    let n = points.len();
    for level in 1..n {
        for i in 0..n - level {
            let next = points[i + 1];
            for (x, y) in points[i].iter_mut().zip(next) {
                *x = (I::point(1.) - t) * *x + t * y;
            }
        }
    }
    points[0]
}
/// Controls of a Bezier curve restricted to `[t0, t1]` (blossoms): exact
/// polynomial reparameterisation, valid for any interval endpoints.
fn restrict_bezier(points: &[V3], t0: I, t1: I, b: &mut Budget<'_>) -> Result<Vec<V3>> {
    let n = points.len();
    check((1..=MAX_DEGREE + 1).contains(&n), ReadError::Unsupported)?;
    b.charge(n * n * n * 24 + 512)?;
    let mut out = Vec::new();
    out.try_reserve_exact(n).map_err(|_| ReadError::Limit)?;
    let mut work = [[I::point(0.); 3]; MAX_DEGREE + 1];
    for k in 0..n {
        work[..n].copy_from_slice(points);
        for level in 1..n {
            let t = if level < n - k { t0 } else { t1 };
            for i in 0..n - level {
                let next = work[i + 1];
                for (x, y) in work[i].iter_mut().zip(next) {
                    *x = (I::point(1.) - t) * *x + t * y;
                }
            }
        }
        check(work[0].iter().all(|v| v.valid()), ReadError::Unsupported)?;
        out.push(work[0]);
    }
    Ok(out)
}
fn value(p: &P, t: I, b: &mut Budget<'_>) -> Result<I> {
    b.charge(512)?;
    let [left, _] = p.split(t, b)?;
    let v = *left.controls().last().ok_or(ReadError::Unsupported)?;
    check(v.valid(), ReadError::Unsupported)?;
    Ok(v)
}
fn norm_hi(v: V3) -> f64 {
    norm(v).hi
}
fn polys(points: &[V3], b: &mut Budget<'_>) -> Result<[P; 3]> {
    let mut out = [P::constant(I::point(0.), b)?; 3];
    for (k, slot) in out.iter_mut().enumerate() {
        let values: Vec<I> = points.iter().map(|v| v[k]).collect();
        *slot = P::new(&values, b)?;
    }
    Ok(out)
}
/// Position and tangent bounds of an exact composite image against the
/// target on the whole cell, from Bernstein coefficient hulls.
fn exact_bounds(image: &[V3], target: &[P; 3], b: &mut Budget<'_>) -> Result<(f64, f64, f64)> {
    let degree = (image.len() - 1).max(target[0].degree());
    let image = polys(image, b)?;
    let mut difference = [P::constant(I::point(0.), b)?; 3];
    for (k, d) in difference.iter_mut().enumerate() {
        let a = image[k].elevate(degree, b)?;
        let c = target[k].elevate(degree, b)?;
        *d = a.add(&c.scale(I::point(-1.), b)?, b)?;
    }
    let maximum = |p: &[P; 3]| -> f64 {
        (0..=p[0].degree())
            .map(|j| norm(std::array::from_fn(|k| p[k].controls()[j])).hi)
            .fold(0_f64, f64::max)
    };
    let position = maximum(&difference);
    let derivative = maximum(&[
        difference[0].derivative(b)?,
        difference[1].derivative(b)?,
        difference[2].derivative(b)?,
    ]);
    let dc = [
        target[0].derivative(b)?,
        target[1].derivative(b)?,
        target[2].derivative(b)?,
    ];
    let ddc = [
        dc[0].derivative(b)?,
        dc[1].derivative(b)?,
        dc[2].derivative(b)?,
    ];
    let m = I::point(0.5);
    let dc_m: V3 = [
        value(&dc[0], m, b)?,
        value(&dc[1], m, b)?,
        value(&dc[2], m, b)?,
    ];
    let speed = norm(dc_m).lo - 0.5 * norm_hi([ddc[0].hull(), ddc[1].hull(), ddc[2].hull()]);
    Ok((position, derivative, speed))
}
/// One exact control row of the original net: a polynomial B-spline in the
/// other surface parameter, with its expanded knot enclosures.
struct Row {
    degree: usize,
    knots: Vec<I>,
    controls: Vec<V3>,
    /// Expanded index of the last copy of each distinct knot.
    ends: Vec<usize>,
}
impl Row {
    /// De Boor blossom of the span ending at expanded knot `s` (so
    /// `knots[s] <= t < knots[s+1]`), with `degree` arguments. Equal point
    /// arguments evaluate; mixed ends give the Bezier controls of a
    /// sub-interval of that span; interval arguments enclose.
    fn blossom(&self, s: usize, args: &[I], b: &mut Budget<'_>) -> Result<V3> {
        let p = self.degree;
        check(
            args.len() == p && s >= p && s < self.controls.len() && s + p < self.knots.len(),
            ReadError::Unsupported,
        )?;
        b.charge((p + 1) * (p + 1) * 48 + 256)?;
        let mut d = [[I::point(0.); 3]; MAX_DEGREE + 1];
        d[..=p].copy_from_slice(&self.controls[s - p..=s]);
        for (r, &t) in (1..=p).zip(args) {
            for i in (r..=p).rev() {
                let lo = self.knots[s - p + i];
                let hi = self.knots[s + 1 + i - r];
                let alpha = (t - lo) / (hi - lo);
                let previous = d[i - 1];
                for (x, y) in d[i].iter_mut().zip(previous) {
                    *x = (I::point(1.) - alpha) * y + alpha * *x;
                }
            }
        }
        check(d[p].iter().all(|v| v.valid()), ReadError::Unsupported)?;
        Ok(d[p])
    }
    /// Bezier controls of span `s` restricted to `[a, z]`.
    fn bezier(&self, s: usize, a: I, z: I, b: &mut Budget<'_>) -> Result<Vec<V3>> {
        let p = self.degree;
        let mut out = Vec::new();
        out.try_reserve_exact(p + 1).map_err(|_| ReadError::Limit)?;
        let mut args = [I::point(0.); MAX_DEGREE];
        for k in 0..=p {
            for (i, slot) in args.iter_mut().enumerate().take(p) {
                *slot = if i < p - k { a } else { z };
            }
            out.push(self.blossom(s, &args[..p], b)?);
        }
        Ok(out)
    }
    /// The exact hodograph row (degree minus one); span `s` becomes `s-1`.
    fn derivative(&self, b: &mut Budget<'_>) -> Result<Self> {
        let p = self.degree;
        check(p >= 1, ReadError::Unsupported)?;
        b.charge(self.controls.len() * 64 + 256)?;
        let mut controls = Vec::new();
        controls
            .try_reserve_exact(self.controls.len() - 1)
            .map_err(|_| ReadError::Limit)?;
        for i in 0..self.controls.len() - 1 {
            let width = self.knots[i + p + 1] - self.knots[i + 1];
            check(width.valid() && width.lo > 0., ReadError::Unsupported)?;
            let scale = I::point(p as f64) / width;
            controls.push(std::array::from_fn(|k| {
                scale * (self.controls[i + 1][k] - self.controls[i][k])
            }));
        }
        Ok(Self {
            degree: p - 1,
            knots: self.knots[1..self.knots.len() - 1].to_vec(),
            controls,
            ends: Vec::new(),
        })
    }
}
/// The original control row of the iso line `axis = constant`, if that
/// constant is a clamped end or an interior knot of multiplicity equal to
/// the degree (where the net interpolates a whole row exactly).
fn row(
    source: &Source<'_>,
    surface: u32,
    net: &crate::brep::NurbsSurface,
    axis: usize,
    constant: &BigRational,
    knots: &[Vec<BigRational>; 2],
    b: &mut Budget<'_>,
) -> Result<Option<Row>> {
    let field = fields(source, surface, 2, b)?;
    check(field.weights.is_none(), ReadError::Unsupported)?;
    let (_, multiplicities) = net.knots(axis);
    check(
        multiplicities.len() == knots[axis].len(),
        ReadError::Unsupported,
    )?;
    let p = net.degrees()[axis];
    let mut hit = None;
    for (k, q) in knots[axis].iter().enumerate() {
        b.charge(128)?;
        if compare(constant, q, b)?.is_eq() {
            hit = Some(k);
            break;
        }
    }
    let Some(k) = hit else {
        return Ok(None);
    };
    let last = knots[axis].len() - 1;
    let counts = net.counts();
    let first: usize = multiplicities[..k].iter().sum();
    let index = match (k, multiplicities[k]) {
        (0, m) if m == p + 1 => 0,
        (k, m) if k == last && m == p + 1 => counts[axis] - 1,
        (k, m) if k != 0 && k != last && m == p => first - 1,
        _ => return Ok(None),
    };
    let other = 1 - axis;
    let rows = list(field.controls)?;
    check(
        rows.len() == counts[0] && counts[other] <= MAX_ROW,
        ReadError::Unsupported,
    )?;
    let mut controls = Vec::new();
    controls
        .try_reserve_exact(counts[other])
        .map_err(|_| ReadError::Limit)?;
    for j in 0..counts[other] {
        b.charge(128)?;
        let (r, c) = if axis == 0 { (index, j) } else { (j, index) };
        let row = list(&rows[r])?;
        check(row.len() == counts[1], ReadError::Unsupported)?;
        controls.push(point(source, reference(&row[c])?, 3, b)?);
    }
    let (_, other_multiplicities) = net.knots(other);
    check(
        other_multiplicities.len() == knots[other].len(),
        ReadError::Unsupported,
    )?;
    let mut expanded = Vec::new();
    let mut ends = Vec::new();
    for (q, &m) in knots[other].iter().zip(other_multiplicities) {
        let value = enclose(q, b)?;
        for _ in 0..m {
            expanded.push(value);
        }
        ends.push(expanded.len() - 1);
    }
    check(
        expanded.len() == controls.len() + net.degrees()[other] + 1,
        ReadError::Unsupported,
    )?;
    Ok(Some(Row {
        degree: net.degrees()[other],
        knots: expanded,
        controls,
        ends,
    }))
}
enum Piece {
    Proved(f64),
    Refused,
    Exceeds,
}
/// One cell of row span `s` over the moving range `[ya, yb]` against the
/// target cell controls, bisected exactly when the bound is loose.
#[allow(
    clippy::too_many_arguments,
    reason = "Row, span, both ends, the target cell, the limit and the depth are independent proof inputs."
)]
fn prove_row(
    row: &Row,
    s: usize,
    ya: I,
    yb: I,
    target: &[V3],
    tol: f64,
    depth: usize,
    b: &mut Budget<'_>,
) -> Result<Piece> {
    b.charge(1024)?;
    let image = row.bezier(s, ya, yb, b)?;
    let (position, derivative, speed) = exact_bounds(&image, &polys(target, b)?, b)?;
    if position.is_finite() && position <= tol && derivative < speed {
        return Ok(Piece::Proved(position));
    }
    let last = target.len() - 1;
    for (e, t) in [
        (image[0], target[0]),
        (image[image.len() - 1], target[last]),
    ] {
        if norm(std::array::from_fn(|k| e[k] - t[k])).lo > tol {
            return Ok(Piece::Exceeds);
        }
    }
    if depth == MAX_DEPTH {
        return Ok(Piece::Refused);
    }
    let middle = (ya + yb) * I::point(0.5);
    let halves = [
        (
            ya,
            middle,
            restrict_bezier(target, I::point(0.), I::point(0.5), b)?,
        ),
        (
            middle,
            yb,
            restrict_bezier(target, I::point(0.5), I::point(1.), b)?,
        ),
    ];
    let mut worst = 0_f64;
    for (a, z, piece) in &halves {
        match prove_row(row, s, *a, *z, piece, tol, depth + 1, b)? {
            Piece::Proved(e) => worst = worst.max(e),
            other => return Ok(other),
        }
    }
    Ok(Piece::Proved(worst))
}
/// ORIGINAL UV control rationals. Every one must lie in the exact surface
/// domain, so the whole trace does (convex hull of the B-spline controls).
fn uv_controls(
    source: &Source<'_>,
    uv_id: u32,
    knots: &[Vec<BigRational>; 2],
    b: &mut Budget<'_>,
) -> Result<Vec<Q2>> {
    let field = fields(source, uv_id, 1, b)?;
    let controls = list(field.controls)?;
    let mut out = Vec::new();
    out.try_reserve_exact(controls.len())
        .map_err(|_| ReadError::Limit)?;
    for p in controls {
        b.charge(256)?;
        let id = reference(p)?;
        let [_, values] = source.args(id, "CARTESIAN_POINT")? else {
            return Err(ReadError::Invalid);
        };
        let values = list(values)?;
        check(values.len() == 2, ReadError::Unsupported)?;
        out.push([
            rational(source, id, &values[0], b)?,
            rational(source, id, &values[1], b)?,
        ]);
    }
    for q in &out {
        for (axis, k) in knots.iter().enumerate() {
            check(
                compare(&q[axis], &k[0], b)?.is_ge()
                    && compare(&q[axis], &k[k.len() - 1], b)?.is_le(),
                ReadError::Unsupported,
            )?;
        }
    }
    Ok(out)
}
fn polynomial_curve(curve: &crate::brep::NurbsCurve) -> bool {
    curve.weights().is_none_or(|w| w.iter().all(|&x| x == 1.))
}
fn has_c0_joint(net: &crate::brep::NurbsSurface) -> bool {
    (0..2).any(|axis| {
        let (_, m) = net.knots(axis);
        m.len() > 2 && m[1..m.len() - 1].contains(&net.degrees()[axis])
    })
}
/// The shared straight-ISO setup: the constant axis's original row, the
/// trace's exact breaks, and the row span each cell lies in.
struct Iso {
    row: Row,
    breaks: Vec<BigRational>,
    spans: Vec<usize>,
    /// The moving surface coordinate at the trace ends, as enclosures.
    moving: [I; 2],
    domain: [I; 2],
}
impl Iso {
    fn at(&self, t: I) -> I {
        self.moving[0]
            + (t - self.domain[0]) / (self.domain[1] - self.domain[0])
                * (self.moving[1] - self.moving[0])
    }
}
fn iso_setup(
    source: &Source<'_>,
    d: &super::super::uv::Definition,
    net: &crate::brep::NurbsSurface,
    uv: &crate::brep::NurbsCurve,
    extra: Option<&[BigRational]>,
    b: &mut Budget<'_>,
) -> Result<Option<Iso>> {
    if !has_c0_joint(net)
        || uv.degree() != 1
        || uv.points().len() != 2
        || !polynomial_curve(uv)
        || net.degrees().iter().any(|&p| p > MAX_DEGREE)
        || net.weights().is_some_and(|w| w.iter().any(|&x| x != 1.))
    {
        return Ok(None);
    }
    let knots = [
        distinct_knots(source, d.surface, 2, 0, b)?,
        distinct_knots(source, d.surface, 2, 1, b)?,
    ];
    let controls = uv_controls(source, d.curve_id, &knots, b)?;
    let ka = distinct_knots(source, d.curve_id, 1, 0, b)?;
    check(controls.len() == 2 && ka.len() == 2, ReadError::Unsupported)?;
    let mut axis = None;
    for (k, (first, last)) in controls[0].iter().zip(&controls[1]).enumerate() {
        if compare(first, last, b)?.is_eq() {
            axis = Some(k);
        }
    }
    let Some(axis) = axis else {
        return Ok(None);
    };
    let Some(row) = row(source, d.surface, net, axis, &controls[0][axis], &knots, b)? else {
        return Ok(None);
    };
    let moving = 1 - axis;
    let delta = operation(&controls[1][moving], &controls[0][moving], false, b)?;
    check(!delta.is_zero(), ReadError::Unsupported)?;
    let span = operation(&ka[1], &ka[0], false, b)?;
    // Exact trace parameter of every interior moving-axis knot it crosses.
    let mut crossings = Vec::new();
    for k in &knots[moving] {
        b.charge(256)?;
        let f = operation(
            &operation(k, &controls[0][moving], false, b)?,
            &delta,
            true,
            b,
        )?;
        if f > BigRational::from_integer(0.into()) && f < BigRational::from_integer(1.into()) {
            crossings.try_reserve(1).map_err(|_| ReadError::Limit)?;
            crossings.push(&ka[0] + &span * &f);
        }
    }
    crossings.sort();
    let mut all = vec![ka[0].clone()];
    all.extend(crossings);
    all.push(ka[1].clone());
    let breaks = match extra {
        Some(extra) => merge(&all, extra, b)?,
        None => all,
    };
    check(breaks.len() <= 2 * MAX_KNOTS, ReadError::Unsupported)?;
    // Each cell's row span, from the moving coordinate at both exact ends.
    let increasing = delta > BigRational::from_integer(0.into());
    let mk = &knots[moving];
    let mut spans = Vec::new();
    spans
        .try_reserve_exact(breaks.len())
        .map_err(|_| ReadError::Limit)?;
    let at = |t: &BigRational, b: &mut Budget<'_>| -> Result<BigRational> {
        let f = operation(&operation(t, &ka[0], false, b)?, &span, true, b)?;
        let step = &delta * &f;
        b.charge(256)?;
        Ok(&controls[0][moving] + step)
    };
    // Cells come in trace order, so the span index moves monotonically in
    // the direction of travel; each search continues from the last span.
    let mut i = if increasing { 0 } else { mk.len() - 2 };
    for pair in breaks.windows(2) {
        b.charge(512)?;
        let [a, z] = [at(&pair[0], b)?, at(&pair[1], b)?];
        let (lo, hi) = if increasing { (&a, &z) } else { (&z, &a) };
        loop {
            if compare(lo, &mk[i], b)?.is_ge() && compare(hi, &mk[i + 1], b)?.is_le() {
                break;
            }
            let next = if increasing { i + 1 } else { i.wrapping_sub(1) };
            check(next + 1 < mk.len(), ReadError::Unsupported)?;
            i = next;
        }
        spans.push(row.ends[i]);
    }
    Ok(Some(Iso {
        row,
        breaks,
        spans,
        moving: [
            enclose(&controls[0][moving], b)?,
            enclose(&controls[1][moving], b)?,
        ],
        domain: [enclose(&ka[0], b)?, enclose(&ka[1], b)?],
    }))
}
#[allow(
    clippy::too_many_arguments,
    reason = "Original owners, decoded curves, the finite edge and the limits are independent proof inputs."
)]
fn circle(
    source: &Source<'_>,
    curve_id: u32,
    d: &super::super::uv::Definition,
    net: &crate::brep::NurbsSurface,
    uv: &crate::brep::NurbsCurve,
    edge: &Edge,
    tol: f64,
    b: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    b.charge(8192)?;
    check(tol.is_finite() && tol > 0., ReadError::Invalid)?;
    let Some(iso) = iso_setup(source, d, net, uv, None, b)? else {
        return Ok(None);
    };
    check(super::super::si_radians(source, b)?, ReadError::Unsupported)?;
    let id = unwrapped(source, curve_id, b)?;
    let [_, placement, radius] = source.args(id, "CIRCLE")? else {
        return Err(ReadError::Unsupported);
    };
    let radius = scalar(source, id, radius, b)?;
    check(radius.valid() && radius.lo > 0., ReadError::Unsupported)?;
    let support = Support {
        frame: frame(source, reference(placement)?, 3, b)?,
        radius,
        slope: I::point(0.),
    };
    // The trace parameter is the circle angle. A vertex-derived edge range
    // may overshoot the exact domain by a few ulps; arc points beyond it lie
    // within `R * overshoot` of the proved end, added below.
    let [t0, t1] = iso.domain;
    check(edge.range[0] < edge.range[1], ReadError::Unsupported)?;
    let overshoot = (t0.lo - edge.range[0]).max(edge.range[1] - t1.hi).max(0.);
    check(overshoot <= 1e-9, ReadError::Unsupported)?;
    let mut maximum = 0_f64;
    for (pair, &s) in iso.breaks.windows(2).zip(&iso.spans) {
        b.charge(2048)?;
        let (ta, tb) = (enclose(&pair[0], b)?, enclose(&pair[1], b)?);
        let image = iso.row.bezier(s, iso.at(ta), iso.at(tb), b)?;
        let target = polys(&image, b)?;
        let theta = P::new(&[ta, tb], b)?;
        let zero = P::new(&[I::point(0.), I::point(0.)], b)?;
        let mut accept = |p: &CellProof, _: &mut Budget<'_>| -> Result<bool> {
            if p.position <= tol && p.derivative < p.speed {
                maximum = maximum.max(p.position);
                return Ok(true);
            }
            Ok(false)
        };
        if !analytic_polynomial::certify(&support, &[theta, zero], &target, tol, &mut accept, b)? {
            return Ok(Some(Proof::Unsupported));
        }
    }
    let total = I::point(maximum) + radius * I::point(overshoot);
    check(total.valid(), ReadError::Unsupported)?;
    if total.hi > tol {
        return Ok(Some(Proof::Unsupported));
    }
    b.charge(64)?;
    Ok(Some(Proof::Certified(total.hi)))
}
pub(super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    d: &super::super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    b: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    let (
        Surface::Nurbs { surface: net, .. },
        Curve::Nurbs {
            curve: uv,
            closed: false,
        },
    ) = (surface, &d.geometry)
    else {
        return Ok(None);
    };
    if !has_c0_joint(net) || uv.degree() != 1 || uv.points().len() != 2 {
        return Ok(None);
    }
    let target = match &edge.curve {
        Curve::Circle { .. } => return circle(source, curve_id, d, net, uv, edge, tol, b),
        Curve::Nurbs {
            curve: target,
            closed: false,
        } if target.degree() <= MAX_DEGREE && polynomial_curve(target) => target,
        _ => return Ok(None),
    };
    b.charge(8192)?;
    check(tol.is_finite() && tol > 0., ReadError::Invalid)?;
    let id = unwrapped(source, curve_id, b)?;
    let kc = distinct_knots(source, id, 1, 0, b)?;
    let Some(iso) = iso_setup(source, d, net, uv, Some(&kc), b)? else {
        return Ok(None);
    };
    let c: CurveCells<'_, '_, '_> =
        bezier::original_curve(source, id, target, Coordinates::Raw(3), b)?;
    check(
        c.matches(source, id) && matches!(c.coordinates(), Coordinates::Raw(3)),
        ReadError::Invalid,
    )?;
    let mut knots_c = Vec::new();
    knots_c
        .try_reserve_exact(kc.len())
        .map_err(|_| ReadError::Limit)?;
    for q in &kc {
        knots_c.push(enclose(q, b)?);
    }
    let derivative = iso.row.derivative(b)?;
    let span_t = iso.domain[1] - iso.domain[0];
    check(span_t.valid() && span_t.lo > 0., ReadError::Unsupported)?;
    let mut maximum = 0_f64;
    let mut cursor = 0usize;
    for (pair, &s) in iso.breaks.windows(2).zip(&iso.spans) {
        b.charge(2048)?;
        let (ta, tb) = (enclose(&pair[0], b)?, enclose(&pair[1], b)?);
        while cursor + 1 < c.cells().len()
            && tb.lo > knots_c[c.cells()[cursor].knot_indices()[1]].hi
        {
            cursor += 1;
        }
        let cell = &c.cells()[cursor];
        let [k0, k1] = cell.knot_indices();
        check(
            ta.lo >= knots_c[k0].lo && tb.hi <= knots_c[k1].hi,
            ReadError::Unsupported,
        )?;
        check(
            cell.controls()
                .iter()
                .all(|h| h[3].lo == 1. && h[3].hi == 1.),
            ReadError::Unsupported,
        )?;
        let width = knots_c[k1] - knots_c[k0];
        let f = [(ta - knots_c[k0]) / width, (tb - knots_c[k0]) / width];
        let controls_c: Vec<V3> = cell.controls().iter().map(|h| [h[0], h[1], h[2]]).collect();
        let (ya, yb) = (iso.at(ta), iso.at(tb));
        // A cell shorter than its own end uncertainty: bound it directly
        // by enclosures over the whole cell in GLOBAL parameter units.
        if (tb - ta).hi <= 1e-9 * span_t.lo {
            let ybox = I::new(ya.lo.min(yb.lo), ya.hi.max(yb.hi));
            let fbox = I::new(f[0].lo.min(f[1].lo), f[0].hi.max(f[1].hi));
            let args = [ybox; MAX_DEGREE];
            let point = iso.row.blossom(s, &args[..iso.row.degree], b)?;
            let partial = derivative.blossom(s - 1, &args[..derivative.degree], b)?;
            let rate = (iso.moving[1] - iso.moving[0]) / span_t;
            let mut row = controls_c.clone();
            let target_point = casteljau(&mut row, fbox);
            let degree = I::point((controls_c.len() - 1) as f64);
            let mut hodograph: Vec<V3> = controls_c
                .windows(2)
                .map(|w| std::array::from_fn(|k| degree * (w[1][k] - w[0][k]) / width))
                .collect();
            b.charge(controls_c.len() * controls_c.len() * 48 + 512)?;
            let target_tangent = casteljau(&mut hodograph, fbox);
            let position = norm(std::array::from_fn(|k| point[k] - target_point[k]));
            let tangent = norm(std::array::from_fn(|k| {
                partial[k] * rate - target_tangent[k]
            }));
            let speed = norm(target_tangent);
            if position.valid()
                && tangent.valid()
                && speed.valid()
                && position.hi <= tol
                && tangent.hi < speed.lo
            {
                maximum = maximum.max(position.hi);
                continue;
            }
            return Ok(Some(Proof::Unsupported));
        }
        let piece = restrict_bezier(&controls_c, f[0], f[1], b)?;
        match prove_row(&iso.row, s, ya, yb, &piece, tol, 0, b)? {
            Piece::Proved(e) => maximum = maximum.max(e),
            Piece::Refused | Piece::Exceeds => return Ok(Some(Proof::Unsupported)),
        }
    }
    b.charge(64)?;
    Ok(Some(Proof::Certified(maximum)))
}

#[cfg(test)]
mod tests;
