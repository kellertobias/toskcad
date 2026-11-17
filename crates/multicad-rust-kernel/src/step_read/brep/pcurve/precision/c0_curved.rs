//! Curved ORIGINAL polynomial UV traces (degree 2..=11) on an original
//! polynomial tensor surface whose interior knots may be only C0, against
//! an original polynomial 3D curve on the SAME parameter.
//!
//! Every exact common cell gets a FOURTH-order interval Taylor model of
//! `r(s) = S(uv(s)) - C(s)` about the cell midpoint: exact midpoint jets of
//! both curves and of the patch (partials of total order <= 3 at the
//! midpoint), plus a whole-cell enclosure of the fourth normalized
//! coefficient (Faa di Bruno with partial hulls of total order <= 4 over
//! the trace's hull box and UV derivative hulls). Position and tangent
//! bounds follow from Taylor's theorem on the half-width 1/2.
//!
//! A UV span whose ORIGINAL controls lie in one closed patch (exact
//! rational test) is bounded with that patch alone. Otherwise every patch
//! whose closed domain meets the hull box is bounded with its own
//! polynomial extended over the box, so joint crossings need no C1
//! assumption. Cells that cannot earn both bounds are bisected exactly;
//! a midpoint already beyond `tol` inside its own patch refuses the pair.
use super::super::Proof;
use super::bezier::{self, Coordinates, CurveCells};
use super::cylinder_polynomial::compare;
use super::surface_bezier::{self, SurfaceCells};
use super::*;
use crate::step_read::brep::check;
const MAX_DEGREE: usize = 11;
const MAX_DEPTH: usize = 24;
const ORDER: usize = 5;
const FACTORIAL: [f64; ORDER] = [1., 1., 2., 6., 24.];
const NET_BYTES: usize = 32 * 1024 * 1024;
type V3 = [I; 3];
type Jet = [I; ORDER];

fn zero() -> Jet {
    [I::point(0.); ORDER]
}
fn mul(a: Jet, c: Jet) -> Jet {
    std::array::from_fn(|n| (0..=n).fold(I::point(0.), |sum, i| sum + a[i] * c[n - i]))
}
fn hull(a: I, c: I) -> I {
    I::new(a.lo.min(c.lo), a.hi.max(c.hi))
}
fn derivative<const N: usize>(points: &[[I; N]]) -> Vec<[I; N]> {
    let n = I::point((points.len() - 1) as f64);
    points
        .windows(2)
        .map(|pair| std::array::from_fn(|axis| n * (pair[1][axis] - pair[0][axis])))
        .collect()
}
fn evaluate<const N: usize>(points: &[[I; N]], t: I, b: &mut Budget<'_>) -> Result<[I; N]> {
    b.charge(points.len() * points.len() * N * 8 + 64)?;
    let mut values = points.to_vec();
    for n in (1..values.len()).rev() {
        for i in 0..n {
            let next = values[i + 1];
            for (value, next) in values[i].iter_mut().zip(next) {
                *value = (I::point(1.) - t) * *value + t * next;
            }
        }
    }
    let v = values[0];
    check(v.iter().all(|x| x.valid()), ReadError::Unsupported)?;
    Ok(v)
}
/// Normalized Taylor coefficients `f^(n)(t)/n!` at `t`.
fn jets<const N: usize>(points: &[[I; N]], t: I, b: &mut Budget<'_>) -> Result<[Jet; N]> {
    let mut out = [zero(); N];
    let mut current = points.to_vec();
    for (n, factorial) in FACTORIAL.into_iter().enumerate() {
        if current.is_empty() {
            break;
        }
        let value = evaluate(&current, t, b)?;
        for axis in 0..N {
            out[axis][n] = value[axis] / I::point(factorial);
        }
        current = derivative(&current);
    }
    Ok(out)
}
/// Whole-interval enclosures of the normalized coefficients (control hulls).
fn jet_bounds<const N: usize>(points: &[[I; N]], b: &mut Budget<'_>) -> Result<[Jet; N]> {
    let mut out = [zero(); N];
    let mut current = points.to_vec();
    for (n, factorial) in FACTORIAL.into_iter().enumerate() {
        if current.is_empty() {
            break;
        }
        b.charge(current.len() * N * 8 + 64)?;
        for axis in 0..N {
            let range = current
                .iter()
                .fold(current[0][axis], |range, p| hull(range, p[axis]));
            out[axis][n] = range / I::point(factorial);
        }
        current = derivative(&current);
    }
    Ok(out)
}
/// A patch net and all its partial nets of total order <= 4, row-major in u.
struct Partials {
    domain: [[I; 2]; 2],
    nets: Vec<(usize, usize, Vec<V3>)>,
}
impl Partials {
    fn new(patch: &surface_bezier::Patch, domain: [[I; 2]; 2], b: &mut Budget<'_>) -> Result<Self> {
        let [width, height] = patch.shape();
        check(
            (2..=MAX_DEGREE + 1).contains(&width)
                && (2..=MAX_DEGREE + 1).contains(&height)
                && patch
                    .controls()
                    .iter()
                    .all(|h| h[3].lo == 1. && h[3].hi == 1.),
            ReadError::Unsupported,
        )?;
        let base: Vec<V3> = patch
            .controls()
            .iter()
            .map(|h| [h[0], h[1], h[2]])
            .collect();
        let mut nets = Vec::new();
        nets.try_reserve_exact(15).map_err(|_| ReadError::Limit)?;
        let mut u = (width, height, base);
        for i in 0..ORDER {
            let mut v = u.clone();
            for j in 0..ORDER - i {
                b.charge(v.2.len() * 48 + 64)?;
                nets.push(v.clone());
                if j + 1 < ORDER - i {
                    v = Self::differentiate(&v, false);
                }
            }
            if i + 1 < ORDER {
                u = Self::differentiate(&u, true);
            }
        }
        Ok(Self { domain, nets })
    }
    fn differentiate(
        (w, h, p): &(usize, usize, Vec<V3>),
        along_u: bool,
    ) -> (usize, usize, Vec<V3>) {
        let (w, h) = (*w, *h);
        let n = if along_u { w } else { h };
        if n == 1 {
            return (w, h, vec![[I::point(0.); 3]; p.len()]);
        }
        let (nw, nh) = if along_u { (w - 1, h) } else { (w, h - 1) };
        let scale = I::point((n - 1) as f64);
        let mut out = Vec::with_capacity(nw * nh);
        for i in 0..nw {
            for j in 0..nh {
                let a = p[i * h + j];
                let c = if along_u {
                    p[(i + 1) * h + j]
                } else {
                    p[i * h + j + 1]
                };
                out.push(std::array::from_fn(|k| scale * (c[k] - a[k])));
            }
        }
        (nw, nh, out)
    }
    fn index(i: usize, j: usize) -> usize {
        ORDER * i - i * i.saturating_sub(1) / 2 + j
    }
    /// The partial (i, j) at a point or, with interval arguments, enclosed
    /// over that box (including its polynomial extension beyond [0,1]^2).
    fn partial(&self, i: usize, j: usize, x: I, y: I, b: &mut Budget<'_>) -> Result<V3> {
        let (w, h, p) = &self.nets[Self::index(i, j)];
        let rows = p
            .chunks_exact(*h)
            .map(|row| evaluate(row, y, b))
            .collect::<Result<Vec<_>>>()?;
        check(rows.len() == *w, ReadError::Unsupported)?;
        evaluate(&rows, x, b)
    }
    fn bytes(&self) -> usize {
        self.nets
            .iter()
            .map(|n| n.2.len() * std::mem::size_of::<V3>())
            .sum()
    }
}
/// Normalized jets of `S(x(s), y(s))` in patch-local coordinates. With
/// `range`, the result's fourth coefficient encloses the composition's
/// fourth normalized derivative over the whole cell.
fn composition(
    patch: &Partials,
    uv: [Jet; 2],
    point: [I; 2],
    range: bool,
    b: &mut Budget<'_>,
) -> Result<[Jet; 3]> {
    let mut delta = uv;
    delta[0][0] = I::point(0.);
    delta[1][0] = I::point(0.);
    let mut one = zero();
    one[0] = I::point(1.);
    let mut powers = [[one; ORDER]; 2];
    for axis in 0..2 {
        for n in 1..ORDER {
            powers[axis][n] = mul(powers[axis][n - 1], delta[axis]);
        }
    }
    let mut out = [zero(); 3];
    for i in 0..ORDER {
        for j in 0..ORDER - i {
            if !range && i + j >= ORDER - 1 {
                continue;
            }
            b.charge(3 * ORDER * ORDER * 8)?;
            let value = patch.partial(i, j, point[0], point[1], b)?;
            let monomial = mul(powers[0][i], powers[1][j]);
            let factor = I::point(FACTORIAL[i] * FACTORIAL[j]);
            for (axis, jet) in out.iter_mut().enumerate() {
                for (n, slot) in jet.iter_mut().enumerate() {
                    *slot = *slot + monomial[n] * (value[axis] / factor);
                }
            }
        }
    }
    check(
        out.iter().flatten().all(|x| x.valid()),
        ReadError::Unsupported,
    )?;
    Ok(out)
}
enum Bound {
    Proved(f64),
    Loose,
    Exceeds,
}
/// Fourth-order bounds of one cell against one patch polynomial.
fn bound(
    patch: &Partials,
    uv: &[[I; 2]],
    target: &[V3],
    allowed: f64,
    tol: f64,
    b: &mut Budget<'_>,
) -> Result<Bound> {
    let [[u0, u1], [v0, v1]] = patch.domain;
    let (wu, wv) = (u1 - u0, v1 - v0);
    check(
        wu.valid() && wu.lo > 0. && wv.valid() && wv.lo > 0.,
        ReadError::Unsupported,
    )?;
    let local: Vec<[I; 2]> = uv
        .iter()
        .map(|p| [(p[0] - u0) / wu, (p[1] - v0) / wv])
        .collect();
    let m = I::point(0.5);
    let mid = jets(&local, m, b)?;
    let range = jet_bounds(&local, b)?;
    let point = [mid[0][0], mid[1][0]];
    let box_ = [range[0][0], range[1][0]];
    let expected = composition(patch, mid, point, false, b)?;
    let expected_range = composition(patch, range, box_, true, b)?;
    let actual = jets(target, m, b)?;
    let actual_range = jet_bounds(target, b)?;
    // r(s) = sum_{n<4} c_n d^n + c_4(xi) d^4 with |d| <= 1/2.
    let mut position = [I::point(0.); 3];
    let mut tangent = [I::point(0.); 3];
    for axis in 0..3 {
        let c = |n: usize| I::point((expected[axis][n] - actual[axis][n]).magnitude());
        let c4 = I::point((expected_range[axis][4] - actual_range[axis][4]).magnitude());
        position[axis] = c(0)
            + c(1) * I::point(0.5)
            + c(2) * I::point(0.25)
            + c(3) * I::point(0.125)
            + c4 * I::point(0.0625);
        tangent[axis] = c(1) + c(2) + c(3) * I::point(0.75) + c4 * I::point(0.5);
    }
    let inside = point.iter().all(|p| p.lo >= 0. && p.hi <= 1.);
    let at_mid = norm(std::array::from_fn(|k| expected[k][0] - actual[k][0]));
    if inside && at_mid.lo > tol {
        return Ok(Bound::Exceeds);
    }
    let position = norm(position);
    let tangent = norm(tangent);
    // Target speed lower bound: |C'(1/2)| - sup|C''| / 2.
    let speed = norm(std::array::from_fn(|k| actual[k][1])).lo
        - 0.5 * norm(std::array::from_fn(|k| actual_range[k][2] * I::point(2.))).hi;
    if position.valid() && tangent.valid() && position.hi <= allowed && tangent.hi < speed {
        Ok(Bound::Proved(position.hi))
    } else {
        Ok(Bound::Loose)
    }
}
enum Piece {
    Proved(f64),
    Refused,
}
struct Patches<'o, 's, 'i> {
    knots: [Vec<I>; 2],
    cells: Vec<Option<usize>>,
    original: &'o SurfaceCells<'s, 'i>,
    cache: Vec<Option<Partials>>,
    bytes: usize,
}
impl Patches<'_, '_, '_> {
    fn spans(&self, axis: usize, box_: I) -> std::ops::Range<usize> {
        let k = &self.knots[axis];
        let first = (0..k.len() - 1)
            .find(|&i| k[i + 1].hi >= box_.lo)
            .unwrap_or(k.len() - 1);
        let last = (0..k.len() - 1)
            .rev()
            .find(|&i| k[i].lo <= box_.hi)
            .map_or(0, |i| i + 1);
        first..last.max(first)
    }
    fn get(&mut self, iu: usize, iv: usize, b: &mut Budget<'_>) -> Result<&Partials> {
        let index = iu * (self.knots[1].len() - 1) + iv;
        let patch = self.cells[index].ok_or(ReadError::Unsupported)?;
        if self.cache[patch].is_none() {
            let domain = [
                [self.knots[0][iu], self.knots[0][iu + 1]],
                [self.knots[1][iv], self.knots[1][iv + 1]],
            ];
            let made = Partials::new(&self.original.patches()[patch], domain, b)?;
            self.bytes = self
                .bytes
                .checked_add(made.bytes())
                .ok_or(ReadError::Limit)?;
            check(self.bytes <= NET_BYTES, ReadError::Limit)?;
            self.cache[patch] = Some(made);
        }
        self.cache[patch].as_ref().ok_or(ReadError::Unsupported)
    }
}
#[allow(
    clippy::too_many_arguments,
    reason = "Patches, both cell curves, the exact span patch, the limits and the depth are independent proof inputs."
)]
fn prove(
    patches: &mut Patches<'_, '_, '_>,
    uv: Vec<[I; 2]>,
    target: Vec<V3>,
    fixed: Option<(usize, usize)>,
    tol: f64,
    depth: usize,
    b: &mut Budget<'_>,
) -> Result<Piece> {
    b.charge(1024)?;
    let (us, vs) = match fixed {
        Some((iu, iv)) => (iu..iu + 1, iv..iv + 1),
        None => {
            let range = jet_bounds(&uv, b)?;
            (patches.spans(0, range[0][0]), patches.spans(1, range[1][0]))
        }
    };
    check(!us.is_empty() && !vs.is_empty(), ReadError::Unsupported)?;
    let mut worst = 0_f64;
    let mut loose = false;
    for iu in us {
        for iv in vs.clone() {
            let patch = patches.get(iu, iv, b)?;
            match bound(patch, &uv, &target, tol, tol, b)? {
                Bound::Proved(e) => worst = worst.max(e),
                Bound::Loose => loose = true,
                Bound::Exceeds => return Ok(Piece::Refused),
            }
        }
    }
    if !loose {
        return Ok(Piece::Proved(worst));
    }
    if depth == MAX_DEPTH {
        return Ok(Piece::Refused);
    }
    let mut worst = 0_f64;
    for half in [[0., 0.5], [0.5, 1.]] {
        let [a, z] = half.map(I::point);
        let uv_half = restrict(&uv, a, z, b)?;
        let target_half = restrict(&target, a, z, b)?;
        match prove(patches, uv_half, target_half, fixed, tol, depth + 1, b)? {
            Piece::Proved(e) => worst = worst.max(e),
            Piece::Refused => return Ok(Piece::Refused),
        }
    }
    Ok(Piece::Proved(worst))
}
/// Bezier controls restricted to `[t0, t1]` by blossoming.
fn restrict<const N: usize>(
    points: &[[I; N]],
    t0: I,
    t1: I,
    b: &mut Budget<'_>,
) -> Result<Vec<[I; N]>> {
    let n = points.len();
    check((1..=MAX_DEGREE + 1).contains(&n), ReadError::Unsupported)?;
    b.charge(n * n * n * N * 8 + 512)?;
    let mut out = Vec::new();
    out.try_reserve_exact(n).map_err(|_| ReadError::Limit)?;
    for k in 0..n {
        let mut work = points.to_vec();
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
/// The patch (closed domain, exact rationals) containing all `points`.
fn locate(
    points: &[[BigRational; 2]],
    knots: &[Vec<BigRational>; 2],
    b: &mut Budget<'_>,
) -> Result<Option<(usize, usize)>> {
    let mut found = [0_usize; 2];
    for (axis, k) in knots.iter().enumerate() {
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
/// Shared exact inputs from the straight-trace route.
pub(super) struct Inputs<'a> {
    pub(super) surface_knots: &'a [Vec<BigRational>; 2],
    pub(super) controls: &'a [[BigRational; 2]],
    pub(super) ka: &'a [BigRational],
    pub(super) kc: &'a [BigRational],
    pub(super) breaks: &'a [BigRational],
}
#[allow(
    clippy::too_many_arguments,
    reason = "Original owners, decoded geometry, exact inputs and limits are independent proof inputs."
)]
pub(super) fn certify(
    source: &Source<'_>,
    d: &super::super::uv::Definition,
    target_id: u32,
    net: &crate::brep::NurbsSurface,
    uv: &crate::brep::NurbsCurve,
    target: &crate::brep::NurbsCurve,
    inputs: Inputs<'_>,
    tol: f64,
    b: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    check(
        (2..=MAX_DEGREE).contains(&uv.degree()),
        ReadError::Unsupported,
    )?;
    // Exact patch of every ORIGINAL UV span (span k uses controls k-p..k).
    let p = uv.degree();
    let mut span_patch = Vec::new();
    span_patch
        .try_reserve_exact(inputs.ka.len())
        .map_err(|_| ReadError::Limit)?;
    let mut expanded = 0_usize;
    for &m in uv.knots().1.iter().take(inputs.ka.len() - 1) {
        expanded += m;
        let k = expanded - 1;
        check(k >= p && k < inputs.controls.len(), ReadError::Unsupported)?;
        span_patch.push(locate(
            &inputs.controls[k - p..=k],
            inputs.surface_knots,
            b,
        )?);
    }
    let original = surface_bezier::original_surface(source, d.surface, net, b)?;
    check(original.matches(source, d.surface), ReadError::Invalid)?;
    check(original.peak_bytes() <= 32 * 1024 * 1024, ReadError::Limit)?;
    let mut knots = [Vec::new(), Vec::new()];
    for (axis, out) in knots.iter_mut().enumerate() {
        out.try_reserve_exact(inputs.surface_knots[axis].len())
            .map_err(|_| ReadError::Limit)?;
        for q in &inputs.surface_knots[axis] {
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
    let mut patches = Patches {
        knots,
        cells,
        original: &original,
        cache,
        bytes: 0,
    };
    let a: CurveCells<'_, '_, '_> =
        bezier::original_curve(source, d.curve_id, uv, Coordinates::Raw(2), b)?;
    let c: CurveCells<'_, '_, '_> =
        bezier::original_curve(source, target_id, target, Coordinates::Raw(3), b)?;
    check(
        a.matches(source, d.curve_id)
            && c.matches(source, target_id)
            && matches!(a.coordinates(), Coordinates::Raw(2))
            && matches!(c.coordinates(), Coordinates::Raw(3)),
        ReadError::Invalid,
    )?;
    let enclosures = |k: &[BigRational], b: &mut Budget<'_>| -> Result<Vec<I>> {
        let mut out = Vec::new();
        out.try_reserve_exact(k.len())
            .map_err(|_| ReadError::Limit)?;
        for q in k {
            out.push(enclose(q, b)?);
        }
        Ok(out)
    };
    let ka = enclosures(inputs.ka, b)?;
    let kc = enclosures(inputs.kc, b)?;
    let mut maximum = 0_f64;
    let (mut ac, mut cc) = (0usize, 0usize);
    for pair in inputs.breaks.windows(2) {
        b.charge(2048)?;
        let (ta, tb) = (enclose(&pair[0], b)?, enclose(&pair[1], b)?);
        let piece = |data: &CurveCells<'_, '_, '_>,
                     knots: &[I],
                     cursor: &mut usize,
                     b: &mut Budget<'_>|
         -> Result<(usize, Vec<V3>)> {
            while *cursor + 1 < data.cells().len()
                && tb.lo > knots[data.cells()[*cursor].knot_indices()[1]].hi
            {
                *cursor += 1;
            }
            let cell = &data.cells()[*cursor];
            let [k0, k1] = cell.knot_indices();
            check(
                ta.lo >= knots[k0].lo && tb.hi <= knots[k1].hi,
                ReadError::Unsupported,
            )?;
            check(
                cell.controls()
                    .iter()
                    .all(|h| h[3].lo == 1. && h[3].hi == 1.),
                ReadError::Unsupported,
            )?;
            let width = knots[k1] - knots[k0];
            let f = [(ta - knots[k0]) / width, (tb - knots[k0]) / width];
            let controls: Vec<V3> = cell.controls().iter().map(|h| [h[0], h[1], h[2]]).collect();
            Ok((k0, restrict(&controls, f[0], f[1], b)?))
        };
        let (span, uv_piece) = piece(&a, &ka, &mut ac, b)?;
        let (_, target_piece) = piece(&c, &kc, &mut cc, b)?;
        let uv_piece: Vec<[I; 2]> = uv_piece.iter().map(|p| [p[0], p[1]]).collect();
        let fixed = *span_patch.get(span).ok_or(ReadError::Unsupported)?;
        match prove(&mut patches, uv_piece, target_piece, fixed, tol, 0, b)? {
            Piece::Proved(e) => maximum = maximum.max(e),
            Piece::Refused => return Ok(Some(Proof::Unsupported)),
        }
    }
    b.charge(64)?;
    Ok(Some(Proof::Certified(maximum)))
}

#[cfg(test)]
mod tests;
