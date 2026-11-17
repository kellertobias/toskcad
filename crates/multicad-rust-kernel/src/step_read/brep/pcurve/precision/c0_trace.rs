//! ORIGINAL polynomial UV traces on a polynomial tensor surface whose
//! interior knots may be only C0 (multiplicity equal to the degree), against
//! an original polynomial 3D curve on the SAME parameter.
//!
//! Each exact common cell of the two curves is bounded by a second-order
//! interval model of `r(s) = S(uv(s)) - C(s)`: the value and first
//! derivative at the cell centre plus a whole-cell second-derivative hull.
//! Every surface patch whose CLOSED domain meets the UV hull box is bounded
//! with its own polynomial extended over that box, so a trace crossing a
//! joint is covered on both sides without any C1 assumption. A cell that
//! cannot earn both the position and the tangent bound is bisected exactly.
//! A centre value already beyond `tol` refuses the pair; nothing is repaired.
use super::super::Proof;
use super::analytic_polynomial::{self, Support};
use super::bezier::{self, Coordinates, CurveCells};
use super::cylinder_polynomial::{CellProof, compare, operation};
use super::polynomial_trace::Polynomial as P;
use super::surface_bezier::{self, SurfaceCells};
use super::*;
use crate::step_read::brep::check;
const MAX_DEGREE: usize = 11;
const MAX_KNOTS: usize = 1025;
const MAX_DEPTH: usize = 40;
const NET_BYTES: usize = 16 * 1024 * 1024;
type V3 = [I; 3];

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
            (None, None) => unreachable!("loop condition"),
        };
        out.push(next.clone());
    }
    Ok(out)
}
/// A tensor Bezier net in patch-local coordinates, row-major in u.
#[derive(Clone)]
struct Net {
    shape: [usize; 2],
    values: Vec<V3>,
}
impl Net {
    fn at(&self, i: usize, j: usize) -> V3 {
        self.values[i * self.shape[1] + j]
    }
    /// Exact hodograph net along `axis` (degree times forward difference).
    fn derivative(&self, axis: usize, b: &mut Budget<'_>) -> Result<Self> {
        let mut shape = self.shape;
        check(shape[axis] >= 2, ReadError::Unsupported)?;
        shape[axis] -= 1;
        let degree = I::point((self.shape[axis] - 1) as f64);
        b.charge(shape[0] * shape[1] * 48 + 256)?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(shape[0] * shape[1])
            .map_err(|_| ReadError::Limit)?;
        for i in 0..shape[0] {
            for j in 0..shape[1] {
                let (a, c) = if axis == 0 {
                    (self.at(i + 1, j), self.at(i, j))
                } else {
                    (self.at(i, j + 1), self.at(i, j))
                };
                values.push(std::array::from_fn(|k| degree * (a[k] - c[k])));
            }
        }
        Ok(Self { shape, values })
    }
    /// Directed de Casteljau at (possibly interval) local coordinates. An
    /// interval argument yields an enclosure of the range over that box,
    /// including boxes extending beyond `[0,1]^2` (polynomial extension).
    fn evaluate(&self, x: I, y: I, b: &mut Budget<'_>) -> Result<V3> {
        let [n, m] = self.shape;
        b.charge(n * m * 3 * 24 + 512)?;
        let mut column = [[I::point(0.); 3]; MAX_DEGREE + 1];
        let mut row = [[I::point(0.); 3]; MAX_DEGREE + 1];
        for (i, slot) in column.iter_mut().enumerate().take(n) {
            for (j, r) in row.iter_mut().enumerate().take(m) {
                *r = self.at(i, j);
            }
            *slot = casteljau(&mut row[..m], y);
        }
        let value = casteljau(&mut column[..n], x);
        check(value.iter().all(|v| v.valid()), ReadError::Unsupported)?;
        Ok(value)
    }
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
/// The value, first and second partial nets of one patch.
struct Nets {
    domain: [[I; 2]; 2],
    s: Net,
    su: Net,
    sv: Net,
}
fn nets(patch: &surface_bezier::Patch, domain: [[I; 2]; 2], b: &mut Budget<'_>) -> Result<Nets> {
    let shape = patch.shape();
    check(
        shape.iter().all(|&n| (2..=MAX_DEGREE + 1).contains(&n))
            && patch
                .controls()
                .iter()
                .all(|h| h[3].lo == 1. && h[3].hi == 1.),
        ReadError::Unsupported,
    )?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(patch.controls().len())
        .map_err(|_| ReadError::Limit)?;
    values.extend(patch.controls().iter().map(|h| [h[0], h[1], h[2]]));
    let s = Net { shape, values };
    let su = s.derivative(0, b)?;
    let sv = s.derivative(1, b)?;
    Ok(Nets { domain, s, su, sv })
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
struct Surface2<'o, 's, 'i> {
    knots: [Vec<I>; 2],
    cells: Vec<Option<usize>>,
    original: &'o SurfaceCells<'s, 'i>,
    cache: Vec<Option<Nets>>,
    bytes: usize,
}
impl Surface2<'_, '_, '_> {
    fn patch(&mut self, iu: usize, iv: usize, b: &mut Budget<'_>) -> Result<&Nets> {
        let index = iu * (self.knots[1].len() - 1) + iv;
        let patch = self.cells[index].ok_or(ReadError::Unsupported)?;
        if self.cache[patch].is_none() {
            let domain = [
                [self.knots[0][iu], self.knots[0][iu + 1]],
                [self.knots[1][iv], self.knots[1][iv + 1]],
            ];
            let made = nets(&self.original.patches()[patch], domain, b)?;
            let size = [&made.s, &made.su, &made.sv]
                .iter()
                .map(|n| n.values.len() * std::mem::size_of::<V3>())
                .sum::<usize>();
            self.bytes = self.bytes.checked_add(size).ok_or(ReadError::Limit)?;
            check(self.bytes <= NET_BYTES, ReadError::Limit)?;
            self.cache[patch] = Some(made);
        }
        self.cache[patch].as_ref().ok_or(ReadError::Unsupported)
    }
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
            let t = if level <= n - 1 - k { t0 } else { t1 };
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
/// Exact composition of a straight ISO trace (one UV coordinate constant,
/// the other affine) with one patch: the patch row or column at the fixed
/// coordinate, restricted to the trace's exact range.
fn iso_image(
    n: &Nets,
    constant: usize,
    x: [I; 2],
    y: [I; 2],
    b: &mut Budget<'_>,
) -> Result<Option<Vec<V3>>> {
    let [nu, nv] = n.s.shape;
    let (along_v, at) = match constant {
        0 => (true, x[0]),
        _ => (false, y[0]),
    };
    let (count, other) = if along_v { (nv, nu) } else { (nu, nv) };
    b.charge(count * other * 3 * 24 + 512)?;
    let mut curve = [[I::point(0.); 3]; MAX_DEGREE + 1];
    let mut column = [[I::point(0.); 3]; MAX_DEGREE + 1];
    for (j, slot) in curve.iter_mut().enumerate().take(count) {
        for (i, c) in column.iter_mut().enumerate().take(other) {
            *c = if along_v { n.s.at(i, j) } else { n.s.at(j, i) };
        }
        *slot = casteljau(&mut column[..other], at);
    }
    let ends = if along_v { y } else { x };
    Ok(Some(restrict_bezier(&curve[..count], ends[0], ends[1], b)?))
}
/// Position and tangent bounds of an exact composite image against the
/// target on the whole cell, from Bernstein coefficient hulls.
fn exact_bounds(image: &[V3], target: &[P; 3], b: &mut Budget<'_>) -> Result<(f64, f64, f64)> {
    let degree = (image.len() - 1).max(target[0].degree());
    let mut difference = [P::constant(I::point(0.), b)?; 3];
    for (k, d) in difference.iter_mut().enumerate() {
        let values: Vec<I> = image.iter().map(|v| v[k]).collect();
        let a = P::new(&values, b)?.elevate(degree, b)?;
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
/// One straight ISO cell inside ONE closed patch, from interval enclosures
/// of its exact ends. Both images are blossom restrictions, so they enclose
/// the exact restricted controls on the same exact parameter range.
#[allow(
    clippy::too_many_arguments,
    reason = "Patch, fixed axis, both local ends, the target cell and the limits are independent proof inputs."
)]
fn prove_iso(
    n: &Nets,
    axis: usize,
    x: [I; 2],
    y: [I; 2],
    target: &[V3],
    tol: f64,
    depth: usize,
    b: &mut Budget<'_>,
) -> Result<Piece> {
    b.charge(1024)?;
    let image = iso_image(n, axis, x, y, b)?.ok_or(ReadError::Unsupported)?;
    let polys: [P; 3] = {
        let mut out = [P::constant(I::point(0.), b)?; 3];
        for (k, slot) in out.iter_mut().enumerate() {
            let values: Vec<I> = target.iter().map(|v| v[k]).collect();
            *slot = P::new(&values, b)?;
        }
        out
    };
    let (position, derivative, speed) = exact_bounds(&image, &polys, b)?;
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
    let half = |e: [I; 2]| [e[0], (e[0] + e[1]) * I::point(0.5), e[1]];
    let (hx, hy) = (half(x), half(y));
    let halves = [
        restrict_bezier(target, I::point(0.), I::point(0.5), b)?,
        restrict_bezier(target, I::point(0.5), I::point(1.), b)?,
    ];
    let mut worst = 0_f64;
    for (k, piece) in halves.iter().enumerate() {
        match prove_iso(
            n,
            axis,
            [hx[k], hx[k + 1]],
            [hy[k], hy[k + 1]],
            piece,
            tol,
            depth + 1,
            b,
        )? {
            Piece::Proved(e) => worst = worst.max(e),
            other => return Ok(other),
        }
    }
    Ok(Piece::Proved(worst))
}
enum Piece {
    Proved(f64),
    Refused,
    Exceeds,
}
type Q2 = [BigRational; 2];
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
        let q = [
            rational(source, id, &values[0], b)?,
            rational(source, id, &values[1], b)?,
        ];
        out.push(q);
    }
    let ends = [0, 1].map(|axis| (&knots[axis][0], &knots[axis][knots[axis].len() - 1]));
    for q in &out {
        for axis in 0..2 {
            check(
                compare(&q[axis], ends[axis].0, b)?.is_ge()
                    && compare(&q[axis], ends[axis].1, b)?.is_le(),
                ReadError::Unsupported,
            )?;
        }
    }
    Ok(out)
}
/// The patch whose CLOSED domain exactly contains every given UV point, if
/// one does. Exact rational comparison: a point on a joint belongs to both.
fn locate(
    points: &[Q2],
    knots: &[Vec<BigRational>; 2],
    b: &mut Budget<'_>,
) -> Result<Option<(usize, usize)>> {
    let mut found = [0_usize; 2];
    for axis in 0..2 {
        let k = &knots[axis];
        let mut hit = None;
        for i in 0..k.len() - 1 {
            b.charge(128)?;
            let mut all = true;
            for q in points {
                if compare(&q[axis], &k[i], b)?.is_lt() || compare(&q[axis], &k[i + 1], b)?.is_gt()
                {
                    all = false;
                    break;
                }
            }
            if all {
                hit = Some(i);
                break;
            }
        }
        let Some(i) = hit else {
            return Ok(None);
        };
        found[axis] = i;
    }
    Ok(Some((found[0], found[1])))
}
fn polynomial_surface(surface: &crate::brep::NurbsSurface) -> bool {
    surface.weights().is_none_or(|w| w.iter().all(|&x| x == 1.))
}
fn polynomial_curve(curve: &crate::brep::NurbsCurve) -> bool {
    curve.weights().is_none_or(|w| w.iter().all(|&x| x == 1.))
}
type PatchIndex = ([Vec<I>; 2], Vec<Option<usize>>, Vec<Option<Nets>>);
/// Enclosed distinct surface knots and the patch of every knot cell.
fn index_patches(
    original: &SurfaceCells<'_, '_>,
    surface_knots: &[Vec<BigRational>; 2],
    b: &mut Budget<'_>,
) -> Result<PatchIndex> {
    let mut knots = [Vec::new(), Vec::new()];
    for (axis, out) in knots.iter_mut().enumerate() {
        out.try_reserve_exact(surface_knots[axis].len())
            .map_err(|_| ReadError::Limit)?;
        for q in &surface_knots[axis] {
            out.push(enclose(q, b)?);
        }
    }
    let spans = [knots[0].len() - 1, knots[1].len() - 1];
    let mut cells = Vec::new();
    cells
        .try_reserve_exact(spans[0] * spans[1])
        .map_err(|_| ReadError::Limit)?;
    cells.resize(spans[0] * spans[1], None);
    for (index, patch) in original.patches().iter().enumerate() {
        b.charge(64)?;
        let [[u0, u1], [v0, v1]] = patch.knot_indices();
        check(u1 == u0 + 1 && v1 == v0 + 1, ReadError::Unsupported)?;
        let slot = cells
            .get_mut(u0 * spans[1] + v0)
            .ok_or(ReadError::Unsupported)?;
        check(slot.is_none(), ReadError::Invalid)?;
        *slot = Some(index);
    }
    let mut cache = Vec::new();
    cache
        .try_reserve_exact(original.patches().len())
        .map_err(|_| ReadError::Limit)?;
    cache.resize_with(original.patches().len(), || None);
    Ok((knots, cells, cache))
}
fn has_c0_joint(net: &crate::brep::NurbsSurface) -> bool {
    (0..2).any(|axis| {
        let (_, m) = net.knots(axis);
        m.len() > 2 && m[1..m.len() - 1].contains(&net.degrees()[axis])
    })
}
/// A straight ISO UV trace on a C0 tensor surface against an original 3D
/// CIRCLE whose angle is the trace parameter. Each cell between exact
/// breaks lies in one closed patch; its exact composite image is proved
/// against the circle with the analytic support proof (circle = cylinder
/// slice at axial coordinate zero).
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
    if !has_c0_joint(net)
        || net.degrees().iter().any(|&p| p > MAX_DEGREE)
        || uv.degree() != 1
        || uv.points().len() != 2
        || !polynomial_surface(net)
        || !polynomial_curve(uv)
    {
        return Ok(None);
    }
    b.charge(8192)?;
    check(tol.is_finite() && tol > 0., ReadError::Invalid)?;
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
    let surface_knots = [
        distinct_knots(source, d.surface, 2, 0, b)?,
        distinct_knots(source, d.surface, 2, 1, b)?,
    ];
    let controls = uv_controls(source, d.curve_id, &surface_knots, b)?;
    let ka = distinct_knots(source, d.curve_id, 1, 0, b)?;
    check(controls.len() == 2 && ka.len() == 2, ReadError::Unsupported)?;
    let mut axis = None;
    for (k, (first, last)) in controls[0].iter().zip(&controls[1]).enumerate() {
        if compare(first, last, b)?.is_eq() {
            axis = Some(k);
        }
    }
    let axis = axis.ok_or(ReadError::Unsupported)?;
    let moving = 1 - axis;
    // The trace parameter is the circle angle: the finite edge must lie in
    // the trace's exact domain.
    let (t0, t1) = (enclose(&ka[0], b)?, enclose(&ka[1], b)?);
    // A vertex-derived edge range may overshoot the exact domain by a few
    // ulps. Arc points beyond it lie within `R * overshoot` of the proved
    // end point, which is added to the certified deviation below.
    check(edge.range[0] < edge.range[1], ReadError::Unsupported)?;
    let overshoot = (t0.lo - edge.range[0]).max(edge.range[1] - t1.hi).max(0.);
    check(overshoot <= 1e-9, ReadError::Unsupported)?;
    let delta = operation(&controls[1][moving], &controls[0][moving], false, b)?;
    check(!delta.is_zero(), ReadError::Unsupported)?;
    let increasing = delta > BigRational::from_integer(0.into());
    let span = operation(&ka[1], &ka[0], false, b)?;
    let mut breaks = vec![ka[0].clone()];
    let mut crossed = std::collections::BTreeMap::new();
    for (index, k) in surface_knots[moving].iter().enumerate() {
        b.charge(256)?;
        let f = operation(
            &operation(k, &controls[0][moving], false, b)?,
            &delta,
            true,
            b,
        )?;
        if f > BigRational::from_integer(0.into()) && f < BigRational::from_integer(1.into()) {
            let t = &ka[0] + &span * &f;
            crossed.insert(t.clone(), index);
            breaks.push(t);
        }
    }
    breaks.sort();
    breaks.push(ka[1].clone());
    check(breaks.len() <= MAX_KNOTS, ReadError::Unsupported)?;
    let original = surface_bezier::original_surface(source, d.surface, net, b)?;
    check(original.matches(source, d.surface), ReadError::Invalid)?;
    check(original.peak_bytes() <= 32 * 1024 * 1024, ReadError::Limit)?;
    let (knots, cells, cache) = index_patches(&original, &surface_knots, b)?;
    let mut s2 = Surface2 {
        knots,
        cells,
        original: &original,
        cache,
        bytes: 0,
    };
    let constant = enclose(&controls[0][axis], b)?;
    let ends = [
        [enclose(&controls[0][0], b)?, enclose(&controls[1][0], b)?],
        [enclose(&controls[0][1], b)?, enclose(&controls[1][1], b)?],
    ];
    let at = |t: I, k: usize| ends[k][0] + (t - t0) / (t1 - t0) * (ends[k][1] - ends[k][0]);
    let first = {
        let f = operation(&operation(&breaks[1], &ka[0], false, b)?, &span, true, b)?;
        let mut second = controls[0].clone();
        second[moving] = &controls[0][moving] + &delta * &f;
        [controls[0].clone(), second]
    };
    let (mut su, mut sv) = locate(&first, &surface_knots, b)?.ok_or(ReadError::Unsupported)?;
    let mut maximum = 0_f64;
    for (index, pair) in breaks.windows(2).enumerate() {
        b.charge(2048)?;
        if index > 0
            && let Some(&knot) = crossed.get(&pair[0])
        {
            let next = if increasing {
                knot
            } else {
                knot.checked_sub(1).ok_or(ReadError::Unsupported)?
            };
            if moving == 0 {
                su = next;
            } else {
                sv = next;
            }
        }
        let (ta, tb) = (enclose(&pair[0], b)?, enclose(&pair[1], b)?);
        let n = s2.patch(su, sv, b)?;
        let [[u0, u1], [v0, v1]] = n.domain;
        let local = |x: I, origin: I, width: I| (x - origin) / width;
        let (x, y) = if axis == 0 {
            (
                [local(constant, u0, u1 - u0); 2],
                [local(at(ta, 1), v0, v1 - v0), local(at(tb, 1), v0, v1 - v0)],
            )
        } else {
            (
                [local(at(ta, 0), u0, u1 - u0), local(at(tb, 0), u0, u1 - u0)],
                [local(constant, v0, v1 - v0); 2],
            )
        };
        let image = iso_image(n, axis, x, y, b)?.ok_or(ReadError::Unsupported)?;
        let mut target = [P::constant(I::point(0.), b)?; 3];
        for (k, slot) in target.iter_mut().enumerate() {
            let values: Vec<I> = image.iter().map(|v| v[k]).collect();
            *slot = P::new(&values, b)?;
        }
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
    if let (
        Surface::Nurbs { surface: net, .. },
        Curve::Nurbs {
            curve: uv,
            closed: false,
        },
        Curve::Circle { .. },
    ) = (surface, &d.geometry, &edge.curve)
    {
        return circle(source, curve_id, d, net, uv, edge, tol, b);
    }
    let (
        Surface::Nurbs { surface: net, .. },
        Curve::Nurbs {
            curve: uv,
            closed: false,
        },
        Curve::Nurbs {
            curve: target,
            closed: false,
        },
    ) = (surface, &d.geometry, &edge.curve)
    else {
        return Ok(None);
    };
    // This route owns tensor surfaces with a C0 interior joint, which the
    // C1 tensor-trace proof refuses by contract.
    let c0 = (0..2).any(|axis| {
        let (_, m) = net.knots(axis);
        m.len() > 2 && m[1..m.len() - 1].contains(&net.degrees()[axis])
    });
    // Only straight UV traces are proved here: their knot-line crossings
    // are exact, so every cell has one exact patch and an exact composite.
    if !c0
        || uv.degree() != 1
        || uv.points().len() != 2
        || net.degrees().iter().any(|&p| p > MAX_DEGREE)
        || target.degree() > MAX_DEGREE
        || !polynomial_surface(net)
        || !polynomial_curve(uv)
        || !polynomial_curve(target)
    {
        return Ok(None);
    }
    b.charge(8192)?;
    check(tol.is_finite() && tol > 0., ReadError::Invalid)?;
    let id = unwrapped(source, curve_id, b)?;
    // Cheap exact obligations first: the trace stays in the surface domain
    // and both curves share one exact parameter domain.
    let surface_knots = [
        distinct_knots(source, d.surface, 2, 0, b)?,
        distinct_knots(source, d.surface, 2, 1, b)?,
    ];
    let controls = uv_controls(source, d.curve_id, &surface_knots, b)?;
    let ka = distinct_knots(source, d.curve_id, 1, 0, b)?;
    let kc = distinct_knots(source, id, 1, 0, b)?;
    let mut points = merge(&ka, &kc, b)?;
    let original = surface_bezier::original_surface(source, d.surface, net, b)?;
    check(original.matches(source, d.surface), ReadError::Invalid)?;
    check(original.peak_bytes() <= 32 * 1024 * 1024, ReadError::Limit)?;
    let (knots, cells, cache) = index_patches(&original, &surface_knots, b)?;
    let mut s2 = Surface2 {
        knots,
        cells,
        original: &original,
        cache,
        bytes: 0,
    };
    let c: CurveCells<'_, '_, '_> =
        bezier::original_curve(source, id, target, Coordinates::Raw(3), b)?;
    check(
        c.matches(source, id) && matches!(c.coordinates(), Coordinates::Raw(3)),
        ReadError::Invalid,
    )?;
    // A straight UV trace crosses each surface knot line at an exact
    // rational parameter; those become cell breaks so every cell lies in
    // one closed patch.
    let line = uv.degree() == 1 && controls.len() == 2 && ka.len() == 2;
    let mut iso = None;
    if line {
        for (axis, (first, last)) in controls[0].iter().zip(&controls[1]).enumerate() {
            if compare(first, last, b)?.is_eq() {
                iso = Some((axis, enclose(first, b)?));
            }
        }
    }
    let along = |t: &BigRational, b: &mut Budget<'_>| -> Result<Q2> {
        let span = operation(&ka[1], &ka[0], false, b)?;
        let f = operation(&operation(t, &ka[0], false, b)?, &span, true, b)?;
        let mut out = [controls[0][0].clone(), controls[0][1].clone()];
        for (axis, slot) in out.iter_mut().enumerate() {
            let delta = operation(&controls[1][axis], &controls[0][axis], false, b)?;
            let scaled = &delta * &f;
            *slot = &controls[0][axis] + scaled;
        }
        Ok(out)
    };
    let mut crossed = std::collections::BTreeMap::new();
    if line {
        let mut crossings = Vec::new();
        for axis in 0..2 {
            let delta = operation(&controls[1][axis], &controls[0][axis], false, b)?;
            if delta.is_zero() {
                continue;
            }
            for (index, k) in surface_knots[axis].iter().enumerate() {
                b.charge(256)?;
                let f = operation(
                    &operation(k, &controls[0][axis], false, b)?,
                    &delta,
                    true,
                    b,
                )?;
                let zero = BigRational::from_integer(0.into());
                let one = BigRational::from_integer(1.into());
                if compare(&f, &zero, b)?.is_gt() && compare(&f, &one, b)?.is_lt() {
                    let span = operation(&ka[1], &ka[0], false, b)?;
                    crossings.try_reserve(1).map_err(|_| ReadError::Limit)?;
                    let t = &ka[0] + &span * &f;
                    crossed.insert(
                        t.clone(),
                        (axis, index, delta > BigRational::from_integer(0.into())),
                    );
                    crossings.push(t);
                }
            }
        }
        crossings.sort();
        crossings.dedup();
        check(crossings.len() < MAX_KNOTS, ReadError::Unsupported)?;
        let mut all = vec![ka[0].clone()];
        all.extend(crossings);
        all.push(ka[1].clone());
        points = merge(&points, &all, b)?;
    }
    let mut maximum = 0_f64;
    if let (true, Some((axis, constant))) = (line, iso) {
        // Every cell lies between two consecutive EXACT breaks (knots and
        // knot-line crossings), so it lies in one closed patch.
        let t0 = enclose(&ka[0], b)?;
        let t1 = enclose(&ka[1], b)?;
        let ends = [0, 1].map(|axis| controls.iter().map(|q| q[axis].clone()).collect::<Vec<_>>());
        let mut uv_ends = [[I::point(0.); 2]; 2];
        for (axis, slot) in uv_ends.iter_mut().enumerate() {
            *slot = [enclose(&ends[axis][0], b)?, enclose(&ends[axis][1], b)?];
        }
        let span = t1 - t0;
        check(span.valid() && span.lo > 0., ReadError::Unsupported)?;
        let at = |t: I, axis: usize| {
            uv_ends[axis][0] + (t - t0) / span * (uv_ends[axis][1] - uv_ends[axis][0])
        };
        let mut knots_c = Vec::new();
        knots_c
            .try_reserve_exact(kc.len())
            .map_err(|_| ReadError::Limit)?;
        for q in &kc {
            knots_c.push(enclose(q, b)?);
        }
        let mut cursor = 0usize;
        // Exact patch of the first cell; every later change of patch happens
        // at an exact knot-line crossing break.
        let first = [along(&points[0], b)?, along(&points[1], b)?];
        let (mut su, mut sv) = locate(&first, &surface_knots, b)?.ok_or(ReadError::Unsupported)?;
        for (index, pair) in points.windows(2).enumerate() {
            b.charge(2048)?;
            if index > 0
                && let Some(&(axis, knot, increasing)) = crossed.get(&pair[0])
            {
                let next = if increasing {
                    knot
                } else {
                    knot.checked_sub(1).ok_or(ReadError::Unsupported)?
                };
                if axis == 0 {
                    su = next;
                } else {
                    sv = next;
                }
            }
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
            let width = knots_c[k1] - knots_c[k0];
            let f = [(ta - knots_c[k0]) / width, (tb - knots_c[k0]) / width];
            let controls_c: Vec<V3> = cell.controls().iter().map(|h| [h[0], h[1], h[2]]).collect();
            check(
                cell.controls()
                    .iter()
                    .all(|h| h[3].lo == 1. && h[3].hi == 1.),
                ReadError::Unsupported,
            )?;
            let piece = restrict_bezier(&controls_c, f[0], f[1], b)?;
            check(
                su < s2.knots[0].len() - 1 && sv < s2.knots[1].len() - 1,
                ReadError::Unsupported,
            )?;
            let n = s2.patch(su, sv, b)?;
            let [[u0, u1], [v0, v1]] = n.domain;
            let local = |x: I, origin: I, width: I| (x - origin) / width;
            let x = [local(at(ta, 0), u0, u1 - u0), local(at(tb, 0), u0, u1 - u0)];
            let y = [local(at(ta, 1), v0, v1 - v0), local(at(tb, 1), v0, v1 - v0)];
            let (x, y) = if axis == 0 {
                ([local(constant, u0, u1 - u0); 2], y)
            } else {
                (x, [local(constant, v0, v1 - v0); 2])
            };
            // A cell shorter than its own end uncertainty: bound it directly
            // by enclosures over the whole cell in GLOBAL parameter units.
            if (tb - ta).hi <= 1e-9 * span.lo {
                let moving = 1 - axis;
                let (fixed_local, moving_box, moving_width) = if axis == 0 {
                    (
                        x[0],
                        I::new(y[0].lo.min(y[1].lo), y[0].hi.max(y[1].hi)),
                        v1 - v0,
                    )
                } else {
                    (
                        y[0],
                        I::new(x[0].lo.min(x[1].lo), x[0].hi.max(x[1].hi)),
                        u1 - u0,
                    )
                };
                let rate = (uv_ends[moving][1] - uv_ends[moving][0]) / span / moving_width;
                let (point, partial) = if axis == 0 {
                    (
                        n.s.evaluate(fixed_local, moving_box, b)?,
                        n.sv.evaluate(fixed_local, moving_box, b)?,
                    )
                } else {
                    (
                        n.s.evaluate(moving_box, fixed_local, b)?,
                        n.su.evaluate(moving_box, fixed_local, b)?,
                    )
                };
                let fbox = I::new(f[0].lo.min(f[1].lo), f[0].hi.max(f[1].hi));
                let mut row = controls_c.clone();
                let target_point = casteljau(&mut row, fbox);
                let degree = I::point((controls_c.len() - 1) as f64);
                let mut hodograph: Vec<V3> = controls_c
                    .windows(2)
                    .map(|w| std::array::from_fn(|k| degree * (w[1][k] - w[0][k]) / width))
                    .collect();
                b.charge(controls_c.len() * controls_c.len() * 48 + 512)?;
                let target_tangent = casteljau(&mut hodograph, fbox);
                let image_tangent: V3 = std::array::from_fn(|k| partial[k] * rate);
                let position = norm(std::array::from_fn(|k| point[k] - target_point[k]));
                let tangent = norm(std::array::from_fn(|k| {
                    image_tangent[k] - target_tangent[k]
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
            match prove_iso(n, axis, x, y, &piece, tol, 0, b)? {
                Piece::Proved(e) => maximum = maximum.max(e),
                Piece::Refused | Piece::Exceeds => return Ok(Some(Proof::Unsupported)),
            }
        }
        b.charge(64)?;
        return Ok(Some(Proof::Certified(maximum)));
    }
    // A straight trace that is not an iso line is not proved here.
    Ok(Some(Proof::Unsupported))
}

#[cfg(test)]
mod tests;
