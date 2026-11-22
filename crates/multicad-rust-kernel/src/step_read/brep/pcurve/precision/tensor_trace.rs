//! Whole ORIGINAL polynomial traces through continuous tensor spline patches.
//! A selected polynomial extension is only a proposal. Actual source partial
//! bounds account for every excursion across a knot before C0 is published;
//! IP2 uses actual patch derivatives, never the derivative of that allowance.
use super::super::Proof;
use super::bernstein80::{self as poly, Poly};
use super::bezier::{Coordinates, CurveCells};
use super::surface_bezier::SurfaceCells;
use super::*;
use crate::step_read::brep::check;
type H = [I; 4];
const NUMERIC_BYTES: usize = 16 * 1024 * 1024;
const MAX_CUTS: usize = 514;
const MAX_DEPTH: usize = 24;
#[cfg(test)]
mod oracle_tests;
#[cfg(test)]
mod source_tests;

fn numeric_preflight(source: &Source<'_>, id: u32, work: &mut Budget<'_>) -> Result<()> {
    fn visit(value: &Value<'_>, count: &mut usize, work: &mut Budget<'_>) -> Result<()> {
        work.charge(16)?;
        match value {
            Value::Number(_) => {
                *count = count.checked_add(1).ok_or(ReadError::Limit)?;
                check(*count <= 2048, ReadError::Limit)?;
            }
            Value::List(values) => {
                for v in values {
                    visit(v, count, work)?
                }
            }
            Value::Typed(e) => {
                for v in &e.args {
                    visit(v, count, work)?
                }
            }
            _ => {}
        }
        Ok(())
    }
    for entity in source.components(id)? {
        work.charge(16)?;
        let mut count = 0;
        for value in &entity.args {
            visit(value, &mut count, work)?;
        }
    }
    Ok(())
}
fn curve_numeric_preflight(source: &Source<'_>, id: u32, work: &mut Budget<'_>) -> Result<()> {
    numeric_preflight(source, id, work)?;
    let field = fields(source, id, 1, work)?;
    let controls = list(field.controls)?;
    check(controls.len() <= 1024, ReadError::Limit)?;
    for p in controls {
        work.charge(32)?;
        numeric_preflight(source, reference(p)?, work)?;
    }
    Ok(())
}

fn q_bits(q: &BigRational) -> Result<usize> {
    let n = q.numer().bits().max(q.denom().bits());
    check(n <= 4096, ReadError::Unsupported)?;
    usize::try_from(n).map_err(|_| ReadError::Limit)
}
fn q_op(
    a: &BigRational,
    b: &BigRational,
    divide: bool,
    work: &mut Budget<'_>,
) -> Result<BigRational> {
    let n = q_bits(a)?.saturating_add(q_bits(b)?).saturating_add(2);
    work.charge(bigint_work(n))?;
    check(!divide || !b.is_zero(), ReadError::Unsupported)?;
    let value = if divide { a / b } else { a - b };
    q_bits(&value)?;
    Ok(value)
}
fn q_cmp(a: &BigRational, b: &BigRational, work: &mut Budget<'_>) -> Result<std::cmp::Ordering> {
    work.charge(bigint_work(
        q_bits(a)?.saturating_add(q_bits(b)?).saturating_add(2),
    ))?;
    Ok(a.cmp(b))
}
/// The positive width of `ends`, shared by every normalization against them.
fn width(ends: &[BigRational; 2], work: &mut Budget<'_>) -> Result<BigRational> {
    let w = q_op(&ends[1], &ends[0], false, work)?;
    check(
        w > BigRational::from_integer(0.into()),
        ReadError::Unsupported,
    )?;
    Ok(w)
}
/// `(q - start) / w` for a width from [`width`].
fn normalized_by(
    q: &BigRational,
    start: &BigRational,
    w: &BigRational,
    work: &mut Budget<'_>,
) -> Result<BigRational> {
    let d = q_op(q, start, false, work)?;
    q_op(&d, w, true, work)
}
fn normalized(
    q: &BigRational,
    ends: &[BigRational; 2],
    work: &mut Budget<'_>,
) -> Result<BigRational> {
    let w = width(ends, work)?;
    normalized_by(q, &ends[0], &w, work)
}
fn raw_knots(
    source: &Source<'_>,
    id: u32,
    dimension: usize,
    work: &mut Budget<'_>,
) -> Result<Vec<Vec<BigRational>>> {
    work.charge(1024)?;
    let f = fields(source, id, dimension, work)?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(dimension)
        .map_err(|_| ReadError::Limit)?;
    for axis in 0..dimension {
        let mut values = Vec::new();
        if let Some((_, k)) = f.knots[axis] {
            let k = list(k)?;
            check((2..=257).contains(&k.len()), ReadError::Unsupported)?;
            work.charge(k.len().saturating_mul(16).saturating_add(64))?;
            values
                .try_reserve_exact(k.len())
                .map_err(|_| ReadError::Limit)?;
            for v in k {
                values.push(rational(source, id, v, work)?);
            }
        } else {
            values.try_reserve_exact(2).map_err(|_| ReadError::Limit)?;
            values.extend([
                BigRational::from_integer(0.into()),
                BigRational::from_integer(1.into()),
            ]);
        }
        check(
            values.capacity().saturating_mul(4096) <= NUMERIC_BYTES,
            ReadError::Limit,
        )?;
        result.push(values);
    }
    Ok(result)
}
fn unit(values: &[H]) -> bool {
    values.iter().all(|p| p[3].lo == 1. && p[3].hi == 1.)
}
fn cv_copy(values: &[H], work: &mut Budget<'_>) -> Result<Vec<H>> {
    check(
        !values.is_empty() && values.len() <= 12,
        ReadError::Unsupported,
    )?;
    work.charge(values.len() * 32 + 128)?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(values.len())
        .map_err(|_| ReadError::Limit)?;
    check(
        result.capacity().saturating_mul(std::mem::size_of::<H>()) <= 32 * 1024,
        ReadError::Limit,
    )?;
    result.extend_from_slice(values);
    Ok(result)
}
fn cv_split(values: &[H], at: I, work: &mut Budget<'_>) -> Result<[Vec<H>; 2]> {
    check(
        at.valid() && at.lo >= 0. && at.hi <= 1.,
        ReadError::Unsupported,
    )?;
    let degree = values.len().checked_sub(1).ok_or(ReadError::Unsupported)?;
    check(degree <= 11 && unit(values), ReadError::Unsupported)?;
    work.charge(degree * (degree + 1) / 2 * 64 + 4096)?;
    let mut row = cv_copy(values, work)?;
    let mut left = cv_copy(values, work)?;
    let mut right = cv_copy(values, work)?;
    left[0] = row[0];
    right[degree] = row[degree];
    for level in 1..=degree {
        for i in 0..=degree - level {
            let next = row[i + 1];
            for (x, y) in row[i].iter_mut().zip(next).take(3) {
                *x = (I::point(1.) - at) * *x + at * y;
            }
            // Both actual endpoints have exact unit weight; their affine
            // blend has exact unit weight independently of interval width.
            row[i][3] = I::point(1.);
        }
        left[level] = row[0];
        right[degree - level] = row[degree - level];
    }
    Ok([left, right])
}
fn restricted(
    values: &[H],
    lo: &BigRational,
    hi: &BigRational,
    work: &mut Budget<'_>,
) -> Result<Vec<H>> {
    let zero = BigRational::from_integer(0.into());
    let one = BigRational::from_integer(1.into());
    check(
        q_cmp(lo, &zero, work)? != std::cmp::Ordering::Less
            && q_cmp(hi, &one, work)? != std::cmp::Ordering::Greater
            && q_cmp(lo, hi, work)? == std::cmp::Ordering::Less,
        ReadError::Unsupported,
    )?;
    let mut result = if hi == &one {
        cv_copy(values, work)?
    } else {
        let [left, _] = cv_split(values, enclose(hi, work)?, work)?;
        left
    };
    if lo != &zero {
        let relative = q_op(lo, hi, true, work)?;
        let [_, right] = cv_split(&result, enclose(&relative, work)?, work)?;
        result = right;
    }
    Ok(result)
}

struct Grid {
    shape: [usize; 2],
    values: Vec<V>,
}
impl Grid {
    fn derivative(&self, axis: usize, span: I, work: &mut Budget<'_>) -> Result<Self> {
        work.charge(self.values.len() * 128 + 512)?;
        check(
            axis < 2 && self.shape[axis] >= 2 && span.valid() && span.lo > 0.,
            ReadError::Unsupported,
        )?;
        let mut shape = self.shape;
        shape[axis] -= 1;
        let mut values = Vec::new();
        values
            .try_reserve_exact(shape[0] * shape[1])
            .map_err(|_| ReadError::Limit)?;
        check(
            values.capacity() * std::mem::size_of::<V>() <= 128 * 1024,
            ReadError::Limit,
        )?;
        let factor = I::point(shape[axis] as f64) / span;
        for u in 0..shape[0] {
            for v in 0..shape[1] {
                let (a, b) = if axis == 0 {
                    (u * self.shape[1] + v, (u + 1) * self.shape[1] + v)
                } else {
                    (u * self.shape[1] + v, u * self.shape[1] + v + 1)
                };
                values.push(std::array::from_fn(|j| {
                    (self.values[b][j] - self.values[a][j]) * factor
                }));
            }
        }
        Ok(Self { shape, values })
    }
    fn hull(&self, work: &mut Budget<'_>) -> Result<V> {
        work.charge(self.values.len() * 32 + 64)?;
        let mut low = [f64::INFINITY; 3];
        let mut high = [f64::NEG_INFINITY; 3];
        for p in &self.values {
            for j in 0..3 {
                check(p[j].valid(), ReadError::Unsupported)?;
                low[j] = low[j].min(p[j].lo);
                high[j] = high[j].max(p[j].hi);
            }
        }
        Ok(std::array::from_fn(|j| I::new(low[j], high[j])))
    }
    fn maximum_norm(&self, work: &mut Budget<'_>) -> Result<f64> {
        work.charge(self.values.len() * 64 + 64)?;
        let mut maximum = 0_f64;
        for value in &self.values {
            let length = norm(*value);
            check(length.valid(), ReadError::Unsupported)?;
            maximum = maximum.max(length.hi);
        }
        Ok(maximum)
    }
    fn second_v_norm(&self, work: &mut Budget<'_>) -> Result<f64> {
        work.charge(64)?;
        if self.shape[1] < 3 {
            return Ok(0.);
        }
        self.derivative(1, I::point(1.), work)?
            .derivative(1, I::point(1.), work)?
            .maximum_norm(work)
    }
    fn at(&self, point: [f64; 2], work: &mut Budget<'_>) -> Result<V> {
        work.charge(self.shape[0] * self.shape[1] * (self.shape[0] + self.shape[1]) * 32 + 1024)?;
        check(
            self.shape.iter().all(|n| *n > 0 && *n <= 12)
                && point
                    .iter()
                    .all(|x| x.is_finite() && (0. ..=1.).contains(x)),
            ReadError::Unsupported,
        )?;
        let mut rows = [[I::point(0.); 3]; 12];
        for (u, out) in rows.iter_mut().enumerate().take(self.shape[0]) {
            let mut row = [[I::point(0.); 3]; 12];
            row[..self.shape[1]]
                .copy_from_slice(&self.values[u * self.shape[1]..(u + 1) * self.shape[1]]);
            for n in (1..self.shape[1]).rev() {
                for i in 0..n {
                    let next = row[i + 1];
                    for (x, y) in row[i].iter_mut().zip(next).take(3) {
                        *x = (I::point(1.) - I::point(point[1])) * *x + I::point(point[1]) * y;
                    }
                }
            }
            *out = row[0];
        }
        for n in (1..self.shape[0]).rev() {
            for i in 0..n {
                let next = rows[i + 1];
                for (x, y) in rows[i].iter_mut().zip(next).take(3) {
                    *x = (I::point(1.) - I::point(point[0])) * *x + I::point(point[0]) * y;
                }
            }
        }
        check(rows[0].iter().all(|x| x.valid()), ReadError::Unsupported)?;
        Ok(rows[0])
    }
    fn mvt(&self, box_: [I; 2], work: &mut Budget<'_>) -> Result<V> {
        check(
            box_.iter().all(|x| x.valid() && x.lo >= 0. && x.hi <= 1.),
            ReadError::Unsupported,
        )?;
        let center = box_.map(I::midpoint);
        let mut value = self.at(center, work)?;
        for axis in 0..2 {
            if self.shape[axis] > 1 {
                let derivative = self.derivative(axis, I::point(1.), work)?.hull(work)?;
                let radius = box_[axis] - I::point(center[axis]);
                for j in 0..3 {
                    value[j] = value[j] + derivative[j] * radius;
                }
            }
        }
        check(value.iter().all(|x| x.valid()), ReadError::Unsupported)?;
        let global = self.hull(work)?;
        for j in 0..3 {
            value[j] = I::new(value[j].lo.max(global[j].lo), value[j].hi.min(global[j].hi));
            check(value[j].valid(), ReadError::Unsupported)?;
        }
        Ok(value)
    }
}
struct Prepared {
    domain: [[I; 2]; 2],
    grid: Grid,
    partial: [Grid; 2],
    maximum: [f64; 2],
    second_v: f64,
}
#[derive(Clone, Copy)]
struct SurfaceBounds {
    first: [f64; 2],
    second_v: f64,
}

fn prepare(
    data: &SurfaceCells<'_, '_>,
    knots: &[Vec<BigRational>],
    ends: &[[BigRational; 2]; 2],
    work: &mut Budget<'_>,
) -> Result<(Vec<Prepared>, SurfaceBounds)> {
    check(
        knots.len() == 2 && data.patches().len() <= 128,
        ReadError::Unsupported,
    )?;
    let mut patches = Vec::new();
    patches
        .try_reserve_exact(data.patches().len())
        .map_err(|_| ReadError::Limit)?;
    let mut maximum = [0_f64; 2];
    let mut second_max = 0_f64;
    let widths = [width(&ends[0], work)?, width(&ends[1], work)?];
    let mut normalized_knots_cache: [Vec<Option<I>>; 2] = [Vec::new(), Vec::new()];
    for (cache, values) in normalized_knots_cache.iter_mut().zip(knots) {
        check(values.len() <= 257, ReadError::Unsupported)?;
        work.charge(values.len() * 16 + 64)?;
        cache
            .try_reserve_exact(values.len())
            .map_err(|_| ReadError::Limit)?;
        cache.resize(values.len(), None);
    }
    let mut held = patches.capacity() * std::mem::size_of::<Prepared>()
        + normalized_knots_cache
            .iter()
            .map(|c| c.capacity() * std::mem::size_of::<Option<I>>())
            .sum::<usize>();
    for (ordinal, patch) in data.patches().iter().enumerate() {
        work.charge(patch.controls().len() * 64 + 1024)?;
        let indices = patch.knot_indices();
        let counts = data.counts();
        check(
            counts[1] > 0
                && indices
                    == [
                        [ordinal / counts[1], ordinal / counts[1] + 1],
                        [ordinal % counts[1], ordinal % counts[1] + 1],
                    ],
            ReadError::Unsupported,
        )?;
        check(
            unit(patch.controls()) && patch.shape().iter().all(|x| *x > 1 && *x <= 12),
            ReadError::Unsupported,
        )?;
        let mut domain = [[I::point(0.); 2]; 2];
        for axis in 0..2 {
            for end in 0..2 {
                let index = indices[axis][end];
                let q = knots[axis].get(index).ok_or(ReadError::Unsupported)?;
                // Adjacent patches share knots: each is normalized once.
                let slot = normalized_knots_cache[axis]
                    .get_mut(index)
                    .ok_or(ReadError::Unsupported)?;
                domain[axis][end] = match *slot {
                    Some(value) => value,
                    None => {
                        let value = enclose(
                            &normalized_by(q, &ends[axis][0], &widths[axis], work)?,
                            work,
                        )?;
                        *slot = Some(value);
                        value
                    }
                };
            }
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(patch.controls().len())
            .map_err(|_| ReadError::Limit)?;
        check(
            values.capacity() * std::mem::size_of::<V>() <= 128 * 1024,
            ReadError::Limit,
        )?;
        for p in patch.controls() {
            values.push([p[0], p[1], p[2]]);
        }
        let grid = Grid {
            shape: patch.shape(),
            values,
        };
        let partial = [
            grid.derivative(0, domain[0][1] - domain[0][0], work)?,
            grid.derivative(1, domain[1][1] - domain[1][0], work)?,
        ];
        let upper = [
            partial[0].maximum_norm(work)?,
            partial[1].maximum_norm(work)?,
        ];
        let second_v = grid.second_v_norm(work)?;
        work.charge(128)?;
        let global_second = I::point(second_v) / (domain[1][1] - domain[1][0]).square();
        check(global_second.valid(), ReadError::Unsupported)?;
        second_max = second_max.max(global_second.hi);
        for axis in 0..2 {
            maximum[axis] = maximum[axis].max(upper[axis]);
        }
        held = held
            .checked_add(grid.values.capacity() * std::mem::size_of::<V>())
            .and_then(|n| {
                n.checked_add(
                    partial
                        .iter()
                        .map(|g| g.values.capacity() * std::mem::size_of::<V>())
                        .sum::<usize>(),
                )
            })
            .ok_or(ReadError::Limit)?;
        check(held <= 4 * 1024 * 1024, ReadError::Limit)?;
        patches.push(Prepared {
            domain,
            grid,
            partial,
            maximum: upper,
            second_v,
        });
    }
    check(
        !patches.is_empty() && maximum.iter().all(|x| x.is_finite()),
        ReadError::Unsupported,
    )?;
    Ok((
        patches,
        SurfaceBounds {
            first: maximum,
            second_v: second_max,
        },
    ))
}
fn uv_hull(uv: &[H], work: &mut Budget<'_>) -> Result<[I; 2]> {
    work.charge(uv.len() * 32 + 128)?;
    check(unit(uv), ReadError::Unsupported)?;
    let mut low = [f64::INFINITY; 2];
    let mut high = [f64::NEG_INFINITY; 2];
    for p in uv {
        for axis in 0..2 {
            check(p[axis].valid(), ReadError::Unsupported)?;
            low[axis] = low[axis].min(p[axis].lo.max(0.));
            high[axis] = high[axis].max(p[axis].hi.min(1.));
        }
    }
    let value = std::array::from_fn(|axis| I::new(low[axis], high[axis]));
    check(
        value.iter().all(|x| x.valid() && x.lo >= 0. && x.hi <= 1.),
        ReadError::Unsupported,
    )?;
    // UnitBox2's exact original normalization and positive Boehm blends
    // prove every actual control lies inside this square, even if the
    // directed representation has a tiny negative/over-one endpoint.
    Ok(value)
}
fn power(x: I, n: usize) -> I {
    (0..n).fold(I::point(1.), |a, _| a * x)
}
fn extension_error(
    patch: &Prepared,
    box_: [I; 2],
    global: SurfaceBounds,
    work: &mut Budget<'_>,
) -> Result<f64> {
    work.charge(1024)?;
    let mut outside = [0_f64; 2];
    let mut epsilon = [I::point(0.); 2];
    for axis in 0..2 {
        let below = patch.domain[axis][0] - I::point(box_[axis].lo);
        let above = I::point(box_[axis].hi) - patch.domain[axis][1];
        check(below.valid() && above.valid(), ReadError::Unsupported)?;
        outside[axis] = below.hi.max(above.hi).max(0.);
        let span = patch.domain[axis][1] - patch.domain[axis][0];
        check(span.valid() && span.lo > 0., ReadError::Unsupported)?;
        epsilon[axis] = I::point(outside[axis]) / I::point(span.lo);
    }
    if patch.grid.shape[0] == 2
        && patch.domain[0][0].lo == 0.
        && patch.domain[0][0].hi == 0.
        && patch.domain[0][1].lo == 1.
        && patch.domain[0][1].hi == 1.
    {
        // The caller first earns ORIGINAL m<p from validated raw integer
        // multiplicities. Hence the source and its selected ORIGINAL
        // polynomial share value AND global-V first jet at the boundary.
        // U remains a true convex [0,1] row blend. Integral Taylor gives
        // ||f-g|| <=(sup||f''||+sup||g''||)*outsideV^2/2, even across
        // additional C1 joins. Both derivatives here use GLOBAL normalized V.
        let span = patch.domain[1][1] - patch.domain[1][0];
        let enlarged = power(
            I::point(1.) + I::point(2.) * epsilon[1],
            patch.grid.shape[1].saturating_sub(3),
        );
        let extension = I::point(patch.second_v) / span.square() * enlarged;
        let error =
            (I::point(global.second_v) + extension) * I::point(outside[1]).square() / I::point(2.);
        check(error.valid(), ReadError::Unsupported)?;
        return Ok(error.hi);
    }
    let mut error = I::point(0.);
    for (axis, &beyond) in outside.iter().enumerate() {
        let degree = patch.partial[axis].shape.map(|n| n - 1);
        let factor = power(I::point(1.) + I::point(2.) * epsilon[0], degree[0])
            * power(I::point(1.) + I::point(2.) * epsilon[1], degree[1]);
        let extension = I::point(patch.maximum[axis]) * factor;
        error = error + (I::point(global.first[axis]) + extension) * I::point(beyond);
    }
    check(
        error.valid() && error.hi.is_finite(),
        ReadError::Unsupported,
    )?;
    Ok(error.hi)
}
fn actual_partials(patches: &[Prepared], box_: [I; 2], work: &mut Budget<'_>) -> Result<[V; 2]> {
    let mut low = [[f64::INFINITY; 3]; 2];
    let mut high = [[f64::NEG_INFINITY; 3]; 2];
    let mut found = false;
    for patch in patches {
        work.charge(256)?;
        if (0..2).any(|axis| {
            box_[axis].hi < patch.domain[axis][0].lo || box_[axis].lo > patch.domain[axis][1].hi
        }) {
            continue;
        }
        let mut local = [I::point(0.); 2];
        for axis in 0..2 {
            let span = patch.domain[axis][1] - patch.domain[axis][0];
            check(span.valid() && span.lo > 0., ReadError::Unsupported)?;
            let coordinate = (box_[axis] - patch.domain[axis][0]) / span;
            check(coordinate.valid(), ReadError::Unsupported)?;
            // The union includes this patch only at actual points inside
            // its own domain, so intersecting with [0,1] is a theorem.
            local[axis] = I::new(coordinate.lo.max(0.), coordinate.hi.min(1.));
        }
        if local.iter().any(|x| !x.valid()) {
            continue;
        }
        for axis in 0..2 {
            let value = patch.partial[axis].mvt(local, work)?;
            for j in 0..3 {
                low[axis][j] = low[axis][j].min(value[j].lo);
                high[axis][j] = high[axis][j].max(value[j].hi);
            }
        }
        found = true;
    }
    check(found, ReadError::Unsupported)?;
    Ok(std::array::from_fn(|axis| {
        std::array::from_fn(|j| I::new(low[axis][j], high[axis][j]))
    }))
}
fn derivative_hull(values: &[H], axis: usize, work: &mut Budget<'_>) -> Result<I> {
    work.charge(values.len() * 32 + 128)?;
    check(
        values.len() >= 2 && values.len() <= 12 && unit(values),
        ReadError::Unsupported,
    )?;
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    for pair in values.windows(2) {
        let value = (pair[1][axis] - pair[0][axis]) * I::point((values.len() - 1) as f64);
        check(value.valid(), ReadError::Unsupported)?;
        low = low.min(value.lo);
        high = high.max(value.hi);
    }
    Ok(I::new(low, high))
}
fn basis(variable: Poly, degree: usize, work: &mut Budget<'_>) -> Result<Vec<Poly>> {
    check(degree <= 11, ReadError::Unsupported)?;
    work.charge((degree + 1) * 32 + 512)?;
    let one = Poly::constant(variable.degree, I::point(1.), work)?;
    let complement = one.add(variable, true, work)?;
    let mut forward = Vec::new();
    let mut backward = Vec::new();
    let mut result = Vec::new();
    for v in [&mut forward, &mut backward, &mut result] {
        v.try_reserve_exact(degree + 1)
            .map_err(|_| ReadError::Limit)?;
        check(
            v.capacity() * std::mem::size_of::<Poly>() <= 32 * 1024,
            ReadError::Limit,
        )?;
    }
    forward.push(Poly::constant(0, I::point(1.), work)?);
    backward.push(Poly::constant(0, I::point(1.), work)?);
    for i in 1..=degree {
        forward.push(poly::product(forward[i - 1], variable, work)?);
        backward.push(poly::product(backward[i - 1], complement, work)?);
    }
    let mut choose = 1_u64;
    for i in 0..=degree {
        work.charge(64)?;
        if i > 0 {
            choose = choose
                .checked_mul((degree - i + 1) as u64)
                .ok_or(ReadError::Unsupported)?
                / (i as u64);
        }
        result.push(
            poly::product(forward[i], backward[degree - i], work)?
                .scale(I::point(choose as f64), work)?,
        );
    }
    Ok(result)
}
fn compose_tensor_bernstein(
    patch: &Prepared,
    uv: &[H],
    work: &mut Budget<'_>,
) -> Result<[Poly; 3]> {
    work.charge(4096)?;
    check(unit(uv), ReadError::Unsupported)?;
    let mut coordinates = [Poly::zero(0)?; 2];
    for (axis, out) in coordinates.iter_mut().enumerate() {
        let degree = uv.len().checked_sub(1).ok_or(ReadError::Unsupported)?;
        let mut control = Poly::from_controls(uv, axis)?;
        // Actual UnitBox2 coordinates and all positive de Casteljau cuts
        // preserve control membership. This intersection only tightens I.
        for c in &mut control.c[..=degree] {
            *c = I::new(c.lo.max(0.), c.hi.min(1.));
            check(c.valid(), ReadError::Unsupported)?;
        }
        let span = patch.domain[axis][1] - patch.domain[axis][0];
        check(span.valid() && span.lo > 0., ReadError::Unsupported)?;
        *out = control
            .add(
                Poly::constant(degree, patch.domain[axis][0], work)?,
                true,
                work,
            )?
            .scale(I::point(1.) / span, work)?;
    }
    let u = basis(coordinates[0], patch.grid.shape[0] - 1, work)?;
    let v = basis(coordinates[1], patch.grid.shape[1] - 1, work)?;
    let degree = u[0]
        .degree
        .checked_add(v[0].degree)
        .ok_or(ReadError::Limit)?;
    let mut result = [Poly::zero(degree)?; 3];
    for (i, a) in u.iter().enumerate() {
        for (j, b) in v.iter().enumerate() {
            let weight = poly::product(*a, *b, work)?;
            for (axis, out) in result.iter_mut().enumerate() {
                *out = out.add(
                    weight.scale(patch.grid.values[i * v.len() + j][axis], work)?,
                    false,
                    work,
                )?;
            }
        }
    }
    Ok(result)
}
fn row_power(grid: &Grid, row: usize, center: f64, work: &mut Budget<'_>) -> Result<[V; 12]> {
    check(
        row < grid.shape[0]
            && grid.shape[1] > 0
            && grid.shape[1] <= 12
            && center.is_finite()
            && (0. ..=1.).contains(&center),
        ReadError::Unsupported,
    )?;
    let n = grid.shape[1] - 1;
    // At n=11: 286 vector blends (2574 Iops),66 vector differences,
    //36 scales, bounded row copies and scalar loops. The complete ceiling
    //includes fixed12-slot arrays and all conversion/validity bookkeeping.
    work.charge(4096 + 64 * (n + 1) * (n + 1))?;
    let mut differences = [[I::point(0.); 3]; 12];
    differences[..=n].copy_from_slice(&grid.values[row * grid.shape[1]..(row + 1) * grid.shape[1]]);
    let mut result = [[I::point(0.); 3]; 12];
    let mut choose = 1_u64;
    for (k, slot) in result.iter_mut().enumerate().take(n + 1) {
        let mut evaluation = differences;
        for length in (1..=n - k).rev() {
            for i in 0..length {
                let next = evaluation[i + 1];
                for (x, y) in evaluation[i].iter_mut().zip(next).take(3) {
                    *x = (I::point(1.) - I::point(center)) * *x + I::point(center) * y;
                }
            }
        }
        if k > 0 {
            choose = choose
                .checked_mul((n - k + 1) as u64)
                .ok_or(ReadError::Unsupported)?
                / (k as u64);
        }
        *slot = evaluation[0].map(|v| v * I::point(choose as f64));
        for i in 0..n - k {
            let next = differences[i + 1];
            for (x, y) in differences[i].iter_mut().zip(next).take(3) {
                *x = y - *x;
            }
        }
    }
    check(
        result[..=n].iter().flatten().all(|v| v.valid()),
        ReadError::Unsupported,
    )?;
    Ok(result)
}
fn power_contract<const ROWS: usize>(
    coefficients: &[[V; 12]; ROWS],
    length: usize,
    variable: Poly,
    work: &mut Budget<'_>,
) -> Result<[[Poly; 3]; ROWS]> {
    check(
        (1..=2).contains(&ROWS) && (1..=12).contains(&length),
        ReadError::Unsupported,
    )?;
    work.charge(4096)?;
    // Horner over the centered power basis: r = c[n-1], r = r*w + c[k].
    // All 3*ROWS coordinate chains multiply by the SAME scalar argument, so
    // each step is one shared-weight product; adding a constant to a
    // Bernstein form adds it to every control. No elevation is needed.
    let mut flat = [Poly::zero(0)?; 6];
    for row in 0..ROWS {
        for j in 0..3 {
            flat[row * 3 + j] = Poly::constant(0, coefficients[row][length - 1][j], work)?;
        }
    }
    for k in (0..length - 1).rev() {
        let mut next = if ROWS == 1 {
            let [a, b, c] = poly::product_many([flat[0], flat[1], flat[2]], variable, work)?;
            [a, b, c, a, b, c]
        } else {
            poly::product_many(flat, variable, work)?
        };
        let degree = next[0].degree;
        work.charge((degree + 1) * 4 * 3 * ROWS + 512)?;
        for row in 0..ROWS {
            for j in 0..3 {
                let p = &mut next[row * 3 + j];
                for c in &mut p.c[..=degree] {
                    *c = *c + coefficients[row][k][j];
                }
            }
        }
        flat = next;
    }
    let result: [[Poly; 3]; ROWS] =
        std::array::from_fn(|row| std::array::from_fn(|j| flat[row * 3 + j]));
    check(
        result.iter().flatten().all(|p| p.valid()),
        ReadError::Unsupported,
    )?;
    Ok(result)
}
fn local_arguments(
    patch: &Prepared,
    uv: &[H],
    work: &mut Budget<'_>,
) -> Result<([Poly; 2], [I; 2])> {
    work.charge(4096)?;
    check(unit(uv), ReadError::Unsupported)?;
    let mut variables = [Poly::from_controls(uv, 0)?, Poly::from_controls(uv, 1)?];
    for (axis, variable) in variables.iter_mut().enumerate() {
        let degree = variable.degree;
        for c in &mut variable.c[..=degree] {
            check(c.valid(), ReadError::Unsupported)?;
            // Exact UnitBox membership belongs to the source-bound trace packet;
            // it remains true under every positive restriction in this caller.
            *c = (I::new(c.lo.max(0.), c.hi.min(1.)) - patch.domain[axis][0])
                / (patch.domain[axis][1] - patch.domain[axis][0]);
        }
        check(variable.valid(), ReadError::Unsupported)?;
    }
    let box_ = uv_hull(uv, work)?;
    let local = std::array::from_fn(|axis| {
        (box_[axis] - patch.domain[axis][0]) / (patch.domain[axis][1] - patch.domain[axis][0])
    });
    check(local.iter().all(|v| v.valid()), ReadError::Unsupported)?;
    Ok((variables, local))
}
fn compose_local(
    grid: &Grid,
    mut variables: [Poly; 2],
    center: f64,
    work: &mut Budget<'_>,
) -> Result<[Poly; 3]> {
    for c in &mut variables[1].c[..=variables[1].degree] {
        *c = *c - I::point(center);
    }
    let a = row_power(grid, 0, center, work)?;
    let b = row_power(grid, 1, center, work)?;
    let mut delta = b;
    work.charge(grid.shape[1] * 32 + 128)?;
    for i in 0..grid.shape[1] {
        for j in 0..3 {
            delta[i][j] = b[i][j] - a[i][j];
        }
    }
    let [a, b] = power_contract(&[a, delta], grid.shape[1], variables[1], work)?;
    let mut result = [Poly::zero(a[0].degree + variables[0].degree)?; 3];
    // The three coordinates share the u argument and the elevation, so each
    // is one shared-weight product.
    let contribution = poly::product_many(b, variables[0], work)?;
    let elevated = poly::elevate_many(a, a[0].degree + variables[0].degree, work)?;
    for j in 0..3 {
        result[j] = elevated[j].add(contribution[j], false, work)?;
    }
    // This is the exact identity A(v)+u*(B(v)-A(v)). Centered power
    // conversion uses all original coefficients; no Taylor term is dropped.
    Ok(result)
}
fn compose(patch: &Prepared, uv: &[H], work: &mut Budget<'_>) -> Result<[Poly; 3]> {
    if patch.grid.shape[0] != 2 {
        return compose_tensor_bernstein(patch, uv, work);
    }
    let (variables, local) = local_arguments(patch, uv, work)?;
    compose_local(
        &patch.grid,
        variables,
        local[1].midpoint().clamp(0., 1.),
        work,
    )
}
fn restricted_rows(grid: &Grid, range: [f64; 2], work: &mut Budget<'_>) -> Result<Grid> {
    let degree = grid.shape[1] - 1;
    check(
        grid.shape[0] == 2
            && degree <= 11
            && range.iter().all(|v| v.is_finite())
            && range[0] < range[1],
        ReadError::Unsupported,
    )?;
    // Each blossom has n(n+1)/2 three-coordinate affine blends. Arbitrary
    // finite literal endpoints are legal polynomial identities, including
    // excursions beyond the ORIGINAL patch; no positive-weight claim about
    // the old domain is used here.
    work.charge(2 * (degree + 1) * degree * (degree + 1) / 2 * 32 + 4096)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(2 * (degree + 1))
        .map_err(|_| ReadError::Limit)?;
    check(
        values.capacity() * std::mem::size_of::<V>() <= 16 * 1024,
        ReadError::Limit,
    )?;
    for row in 0..2 {
        for coefficient in 0..=degree {
            let mut data = [[I::point(0.); 3]; 12];
            data[..=degree]
                .copy_from_slice(&grid.values[row * (degree + 1)..(row + 1) * (degree + 1)]);
            for step in 0..degree {
                let t = I::point(if step < degree - coefficient {
                    range[0]
                } else {
                    range[1]
                });
                for i in 0..degree - step {
                    let next = data[i + 1];
                    for (x, y) in data[i].iter_mut().zip(next).take(3) {
                        *x = (I::point(1.) - t) * *x + t * y;
                    }
                }
            }
            check(data[0].iter().all(|v| v.valid()), ReadError::Unsupported)?;
            values.push(data[0]);
        }
    }
    Ok(Grid {
        shape: grid.shape,
        values,
    })
}
fn small_choose(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    let mut result = 1_u64;
    for i in 1..=k {
        result = result * (n - i + 1) as u64 / i as u64;
    }
    result
}
fn proposal(grid: &Grid, degree: usize, work: &mut Budget<'_>) -> Result<Grid> {
    check(
        grid.shape == [2, 10] && [3, 5, 7].contains(&degree),
        ReadError::Unsupported,
    )?;
    work.charge(8192)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(2 * (degree + 1))
        .map_err(|_| ReadError::Limit)?;
    for row in 0..2 {
        let mut candidates = [[0_f64; 3]; 8];
        for end in 0..2 {
            let mut differences = [[0_f64; 3]; 10];
            for (i, p) in differences.iter_mut().enumerate() {
                let index = if end == 0 { i } else { 9 - i };
                *p = grid.values[row * 10 + index].map(I::midpoint);
            }
            let mut front = [[0_f64; 3]; 4];
            for order in 0..=(degree - 1) / 2 {
                let factor = small_choose(9, order) as f64 / small_choose(degree, order) as f64;
                for axis in 0..3 {
                    let mut value = factor * differences[0][axis];
                    for (i, previous) in front.iter().enumerate().take(order) {
                        let sign = if (order - i) % 2 == 0 { 1. } else { -1. };
                        value -= sign * small_choose(order, i) as f64 * previous[axis];
                    }
                    front[order][axis] = value;
                }
                let index = if end == 0 { order } else { degree - order };
                candidates[index] = front[order];
                for i in 0..9 - order {
                    let next = differences[i + 1];
                    for (x, y) in differences[i].iter_mut().zip(next).take(3) {
                        *x = y - *x;
                    }
                }
            }
        }
        for c in &candidates[..=degree] {
            check(c.iter().all(|x| x.is_finite()), ReadError::Unsupported)?;
            values.push(c.map(I::point));
        }
    }
    Ok(Grid {
        shape: [2, degree + 1],
        values,
    })
}
fn elevation_error(original: &Grid, auxiliary: &Grid, work: &mut Budget<'_>) -> Result<f64> {
    let n = original.shape[1] - 1;
    let m = auxiliary.shape[1] - 1;
    check(
        original.shape[0] == 2 && auxiliary.shape[0] == 2 && m <= n && n <= 11,
        ReadError::Unsupported,
    )?;
    work.charge(2 * (n + 1) * (m + 1) * 64 + 4096)?;
    let mut error = 0_f64;
    for row in 0..2 {
        for k in 0..=n {
            let mut raised = [I::point(0.); 3];
            for i in k.saturating_sub(n - m)..=k.min(m) {
                let weight = I::point(small_choose(m, i) as f64)
                    * I::point(small_choose(n - m, k - i) as f64)
                    / I::point(small_choose(n, k) as f64);
                for (r, a) in raised.iter_mut().zip(auxiliary.values[row * (m + 1) + i]) {
                    *r = *r + weight * a;
                }
            }
            let difference =
                std::array::from_fn(|j| original.values[row * (n + 1) + k][j] - raised[j]);
            let value = norm(difference);
            check(value.valid(), ReadError::Unsupported)?;
            error = error.max(value.hi);
        }
    }
    Ok(error)
}
fn row_evaluate(grid: &Grid, row: usize, parameter: f64, work: &mut Budget<'_>) -> Result<V> {
    let n = grid.shape[1] - 1;
    check(
        row < grid.shape[0] && n <= 11 && parameter.is_finite(),
        ReadError::Unsupported,
    )?;
    work.charge(n * (n + 1) / 2 * 32 + 1024)?;
    let mut data = [[I::point(0.); 3]; 12];
    data[..=n].copy_from_slice(&grid.values[row * (n + 1)..(row + 1) * (n + 1)]);
    let t = I::point(parameter);
    for length in (1..=n).rev() {
        for i in 0..length {
            let next = data[i + 1];
            for (x, y) in data[i].iter_mut().zip(next).take(3) {
                *x = (I::point(1.) - t) * *x + t * y;
            }
        }
    }
    check(data[0].iter().all(|v| v.valid()), ReadError::Unsupported)?;
    Ok(data[0])
}
fn chord_model(
    grid: &Grid,
    range: [f64; 2],
    second_v: f64,
    allowed: f64,
    work: &mut Budget<'_>,
) -> Result<Option<(Grid, f64)>> {
    work.charge(1024)?;
    check(
        (1..=2).contains(&grid.shape[0])
            && grid.shape[1] >= 2
            && grid.shape[1] <= 12
            && range.iter().all(|v| v.is_finite())
            && range[0] < range[1]
            && second_v.is_finite()
            && second_v >= 0.,
        ReadError::Unsupported,
    )?;
    let below = I::point(0.) - I::point(range[0]);
    let above = I::point(range[1]) - I::point(1.);
    check(below.valid() && above.valid(), ReadError::Unsupported)?;
    let epsilon = below.hi.max(above.hi).max(0.);
    let amplitude = power(
        I::point(1.) + I::point(2.) * I::point(epsilon),
        grid.shape[1].saturating_sub(3),
    );
    let width = I::point(range[1]) - I::point(range[0]);
    let interpolation = I::point(second_v) * amplitude * width.square() / I::point(8.);
    check(interpolation.valid(), ReadError::Unsupported)?;
    let limit = (I::point(allowed) * I::point(0.125)).lo;
    if interpolation.hi > limit {
        return Ok(None);
    }
    work.charge(2048)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(grid.shape[0] * 2)
        .map_err(|_| ReadError::Limit)?;
    check(
        values.capacity() * std::mem::size_of::<V>() <= 16 * 1024,
        ReadError::Limit,
    )?;
    let mut endpoint = 0_f64;
    for row in 0..grid.shape[0] {
        for t in range {
            let actual = row_evaluate(grid, row, t, work)?;
            let literal = actual.map(|v| I::point(v.midpoint()));
            let error = norm(std::array::from_fn(|j| actual[j] - literal[j]));
            check(error.valid(), ReadError::Unsupported)?;
            endpoint = endpoint.max(error.hi);
            values.push(literal);
        }
    }
    // Both M2 and width use the SAME ORIGINAL polynomial's local chart.
    // The absolute Bernstein partition bound covers its entire extended
    // [a,b], not merely the old[0,1]. U is a true convex combination.
    let error = interpolation + I::point(endpoint);
    check(error.valid(), ReadError::Unsupported)?;
    if error.hi > limit {
        return Ok(None);
    }
    Ok(Some((
        Grid {
            shape: [grid.shape[0], 2],
            values,
        },
        error.hi,
    )))
}
fn position_composition(
    patch: &Prepared,
    uv: &[H],
    allowed: f64,
    work: &mut Budget<'_>,
) -> Result<([Poly; 3], f64)> {
    if patch.grid.shape != [2, 10]
        || patch.domain[0][0].lo != 0.
        || patch.domain[0][0].hi != 0.
        || patch.domain[0][1].lo != 1.
        || patch.domain[0][1].hi != 1.
    {
        return Ok((compose(patch, uv, work)?, 0.));
    }
    let (variables, local) = local_arguments(patch, uv, work)?;
    let range = [local[1].lo.next_down(), local[1].hi.next_up()];
    if !range.iter().all(|v| v.is_finite()) || range[0] >= range[1] {
        return Ok((compose(patch, uv, work)?, 0.));
    }
    if let Some((auxiliary, error)) =
        chord_model(&patch.grid, range, patch.second_v, allowed, work)?
    {
        let mut arguments = variables;
        let span = I::point(range[1]) - I::point(range[0]);
        check(span.valid() && span.lo > 0., ReadError::Unsupported)?;
        for c in &mut arguments[1].c[..=arguments[1].degree] {
            let q = (*c - I::point(range[0])) / span;
            check(q.valid(), ReadError::Unsupported)?;
            *c = I::new(q.lo.max(0.), q.hi.min(1.));
        }
        return Ok((compose_local(&auxiliary, arguments, 0.5, work)?, error));
    }
    let restricted = restricted_rows(&patch.grid, range, work)?;
    for degree in [3, 5, 7] {
        let auxiliary = proposal(&restricted, degree, work)?;
        let error = elevation_error(&restricted, &auxiliary, work)?;
        if error > (I::point(allowed) * I::point(0.125)).lo {
            continue;
        }
        let mut arguments = variables;
        let span = I::point(range[1]) - I::point(range[0]);
        check(span.valid() && span.lo > 0., ReadError::Unsupported)?;
        for coefficient in &mut arguments[1].c[..=arguments[1].degree] {
            let coordinate = (*coefficient - I::point(range[0])) / span;
            check(coordinate.valid(), ReadError::Unsupported)?;
            // The literal restriction interval contains the ENTIRE active
            // ORIGINAL argument/control box, even beyond its old patch.
            // Its NEW argument is therefore inside[0,1], where the complete
            // ORIGINAL-minus-elevated coefficient residual is convex.
            *coefficient = I::new(coordinate.lo.max(0.), coordinate.hi.min(1.));
        }
        return Ok((compose_local(&auxiliary, arguments, 0.5, work)?, error));
    }
    Ok((compose(patch, uv, work)?, 0.))
}
struct Node {
    uv: Vec<H>,
    target: Vec<H>,
    depth: usize,
}
fn normalized_knots(values: &[BigRational], work: &mut Budget<'_>) -> Result<Vec<BigRational>> {
    check(
        values.len() >= 2 && values.len() <= 257,
        ReadError::Unsupported,
    )?;
    work.charge(values.len() * 16 + 128)?;
    let ends = [values[0].clone(), values[values.len() - 1].clone()];
    let mut result = Vec::new();
    result
        .try_reserve_exact(values.len())
        .map_err(|_| ReadError::Limit)?;
    check(result.capacity() * 4096 <= NUMERIC_BYTES, ReadError::Limit)?;
    let w = width(&ends, work)?;
    for value in values {
        result.push(normalized_by(value, &ends[0], &w, work)?);
    }
    Ok(result)
}
fn union(a: &[BigRational], b: &[BigRational], work: &mut Budget<'_>) -> Result<Vec<BigRational>> {
    let bound = a.len().checked_add(b.len()).ok_or(ReadError::Limit)?;
    check(bound <= MAX_CUTS, ReadError::Unsupported)?;
    work.charge(bound * 32 + 128)?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(bound)
        .map_err(|_| ReadError::Limit)?;
    check(result.capacity() * 4096 <= NUMERIC_BYTES, ReadError::Limit)?;
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < b.len() {
        work.charge(32)?;
        let value = if i == a.len() {
            let q = &b[j];
            j += 1;
            q
        } else if j == b.len() {
            let q = &a[i];
            i += 1;
            q
        } else {
            match q_cmp(&a[i], &b[j], work)? {
                std::cmp::Ordering::Less => {
                    let q = &a[i];
                    i += 1;
                    q
                }
                std::cmp::Ordering::Greater => {
                    let q = &b[j];
                    j += 1;
                    q
                }
                std::cmp::Ordering::Equal => {
                    let q = &a[i];
                    i += 1;
                    j += 1;
                    q
                }
            }
        };
        work.charge(bigint_work(q_bits(value)?))?;
        result.push(value.clone());
    }
    check(result.len() >= 2, ReadError::Unsupported)?;
    Ok(result)
}
fn restriction_indices(
    values: &[BigRational],
    low: &BigRational,
    high: &BigRational,
    index: &mut usize,
    work: &mut Budget<'_>,
) -> Result<[BigRational; 2]> {
    while *index + 1 < values.len() - 1
        && q_cmp(low, &values[*index + 1], work)? != std::cmp::Ordering::Less
    {
        work.charge(32)?;
        *index += 1;
    }
    check(
        *index + 1 < values.len()
            && q_cmp(low, &values[*index], work)? != std::cmp::Ordering::Less
            && q_cmp(high, &values[*index + 1], work)? != std::cmp::Ordering::Greater,
        ReadError::Unsupported,
    )?;
    let ends = [values[*index].clone(), values[*index + 1].clone()];
    let w = width(&ends, work)?;
    Ok([
        normalized_by(low, &ends[0], &w, work)?,
        normalized_by(high, &ends[0], &w, work)?,
    ])
}
fn affine_q(
    value: &BigRational,
    from: &[BigRational; 2],
    to: &[BigRational; 2],
    work: &mut Budget<'_>,
) -> Result<BigRational> {
    let t = normalized(value, from, work)?;
    let width = q_op(&to[1], &to[0], false, work)?;
    let bits = q_bits(&t)?
        .saturating_add(q_bits(&width)?)
        .saturating_add(q_bits(&to[0])?)
        .saturating_add(4);
    work.charge(bigint_work(bits).saturating_mul(2))?;
    let result = &to[0] + width * t;
    q_bits(&result)?;
    Ok(result)
}
fn propose_parameter(
    q: f64,
    uv: &crate::brep::NurbsCurve,
    surface: &crate::brep::NurbsSurface,
    target: &crate::brep::NurbsCurve,
    window: [f64; 2],
    work: &mut Budget<'_>,
) -> Result<f64> {
    work.charge(1024)?;
    check(
        q.is_finite() && (0. ..=1.).contains(&q),
        ReadError::Unsupported,
    )?;
    if q == 0. {
        return Ok(window[0]);
    }
    if q == 1. {
        return Ok(window[1]);
    }
    let ur = uv.range();
    let tr = target.range();
    work.charge(16_384 + 8 * (uv.points().len() + surface.counts().iter().sum::<usize>()))?;
    let (p, _) = uv.evaluate(ur[0] + q * (ur[1] - ur[0]));
    let (point, _) = surface.evaluate(p[0], p[1]);
    check(point.iter().all(|v| v.is_finite()), ReadError::Unsupported)?;
    let mut proposal = window[0] + q * (window[1] - window[0]);
    for _ in 0..6 {
        work.charge(4096 + 8 * target.points().len())?;
        let (p, d) = target.evaluate(tr[0] + proposal * (tr[1] - tr[0]));
        let d = d.map(|v| v * (tr[1] - tr[0]));
        let den = d.iter().map(|v| v * v).sum::<f64>();
        check(
            den.is_finite() && den > 0. && p.iter().all(|v| v.is_finite()),
            ReadError::Unsupported,
        )?;
        let step = (0..3).map(|j| (p[j] - point[j]) * d[j]).sum::<f64>() / den;
        check(step.is_finite(), ReadError::Unsupported)?;
        let next = (proposal - step.clamp(-0.01, 0.01)).clamp(0., 1.);
        let converged = next == proposal || step.abs() <= 32. * f64::EPSILON;
        proposal = next;
        if converged {
            break;
        }
    }
    check(proposal.is_finite(), ReadError::Unsupported)?;
    Ok(proposal)
}
struct TargetPatch {
    domain: [I; 2],
    grid: Grid,
    derivative: Grid,
    maximum: f64,
    second_v: f64,
}
fn target_patches(
    packet: &CurveCells<'_, '_, '_>,
    knots: &[BigRational],
    work: &mut Budget<'_>,
) -> Result<(Vec<TargetPatch>, f64)> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(packet.cells().len())
        .map_err(|_| ReadError::Limit)?;
    let mut maximum = 0_f64;
    let mut held = result.capacity() * std::mem::size_of::<TargetPatch>();
    for (i, cell) in packet.cells().iter().enumerate() {
        work.charge(cell.controls().len() * 64 + 1024)?;
        check(
            cell.knot_indices() == [i, i + 1] && unit(cell.controls()),
            ReadError::Unsupported,
        )?;
        let domain = [enclose(&knots[i], work)?, enclose(&knots[i + 1], work)?];
        let mut values = Vec::new();
        values
            .try_reserve_exact(cell.controls().len())
            .map_err(|_| ReadError::Limit)?;
        values.extend(cell.controls().iter().map(|h| [h[0], h[1], h[2]]));
        let grid = Grid {
            shape: [1, values.len()],
            values,
        };
        let derivative = grid.derivative(1, domain[1] - domain[0], work)?;
        let upper = derivative.maximum_norm(work)?;
        maximum = maximum.max(upper);
        let second_v = grid.second_v_norm(work)?;
        held = held
            .checked_add(
                (grid.values.capacity() + derivative.values.capacity()) * std::mem::size_of::<V>(),
            )
            .ok_or(ReadError::Limit)?;
        check(held <= 1024 * 1024, ReadError::Limit)?;
        result.push(TargetPatch {
            domain,
            grid,
            derivative,
            maximum: upper,
            second_v,
        });
    }
    Ok((result, maximum))
}
fn target_box(phi: &[H], work: &mut Budget<'_>) -> Result<I> {
    work.charge(phi.len() * 32 + 128)?;
    check(unit(phi), ReadError::Unsupported)?;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in phi {
        check(p[0].valid(), ReadError::Unsupported)?;
        lo = lo.min(p[0].lo.max(0.));
        hi = hi.max(p[0].hi.min(1.));
    }
    let x = I::new(lo, hi);
    check(x.valid(), ReadError::Unsupported)?;
    Ok(x)
}
fn target_collar(p: &TargetPatch, box_: I, global: f64, work: &mut Budget<'_>) -> Result<f64> {
    work.charge(512)?;
    let below = p.domain[0] - I::point(box_.lo);
    let above = I::point(box_.hi) - p.domain[1];
    check(below.valid() && above.valid(), ReadError::Unsupported)?;
    let outside = below.hi.max(above.hi).max(0.);
    let span = p.domain[1] - p.domain[0];
    check(span.valid() && span.lo > 0., ReadError::Unsupported)?;
    let epsilon = I::point(outside) / I::point(span.lo);
    let extension =
        I::point(p.maximum) * power(I::point(1.) + I::point(2.) * epsilon, p.grid.shape[1] - 2);
    let error = (I::point(global) + extension) * I::point(outside);
    check(error.valid(), ReadError::Unsupported)?;
    Ok(error.hi)
}
fn target_derivative(patches: &[TargetPatch], box_: I, work: &mut Budget<'_>) -> Result<V> {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    let mut found = false;
    for p in patches {
        work.charge(128)?;
        if box_.hi < p.domain[0].lo || box_.lo > p.domain[1].hi {
            continue;
        }
        let local = (box_ - p.domain[0]) / (p.domain[1] - p.domain[0]);
        check(local.valid(), ReadError::Unsupported)?;
        let local = I::new(local.lo.max(0.), local.hi.min(1.));
        if !local.valid() {
            continue;
        }
        let value = p.derivative.mvt([I::point(0.), local], work)?;
        for j in 0..3 {
            lo[j] = lo[j].min(value[j].lo);
            hi[j] = hi[j].max(value[j].hi);
        }
        found = true;
    }
    check(found, ReadError::Unsupported)?;
    Ok(std::array::from_fn(|j| I::new(lo[j], hi[j])))
}
fn compose_target(p: &TargetPatch, phi: &[H], work: &mut Budget<'_>) -> Result<[Poly; 3]> {
    let mut argument = Poly::from_controls(phi, 0)?;
    for coefficient in &mut argument.c[..=argument.degree] {
        *coefficient = (*coefficient - p.domain[0]) / (p.domain[1] - p.domain[0]);
    }
    let box_ = target_box(phi, work)?;
    let local = (box_ - p.domain[0]) / (p.domain[1] - p.domain[0]);
    check(local.valid(), ReadError::Unsupported)?;
    let center = local.midpoint().clamp(0., 1.);
    for coefficient in &mut argument.c[..=argument.degree] {
        *coefficient = *coefficient - I::point(center);
    }
    let coefficients = row_power(&p.grid, 0, center, work)?;
    Ok(power_contract(&[coefficients], p.grid.shape[1], argument, work)?[0])
}
fn target_position(
    p: &TargetPatch,
    phi: &[H],
    allowed: f64,
    work: &mut Budget<'_>,
) -> Result<([Poly; 3], f64)> {
    let box_ = target_box(phi, work)?;
    let local = (box_ - p.domain[0]) / (p.domain[1] - p.domain[0]);
    check(local.valid(), ReadError::Unsupported)?;
    let range = [local.lo.next_down(), local.hi.next_up()];
    if let Some((model, error)) = chord_model(&p.grid, range, p.second_v, allowed, work)? {
        let span = I::point(range[1]) - I::point(range[0]);
        check(span.valid() && span.lo > 0., ReadError::Unsupported)?;
        let mut argument = Poly::from_controls(phi, 0)?;
        for c in &mut argument.c[..=argument.degree] {
            let old = (*c - p.domain[0]) / (p.domain[1] - p.domain[0]);
            let q = (old - I::point(range[0])) / span;
            check(q.valid(), ReadError::Unsupported)?;
            *c = I::new(q.lo.max(0.), q.hi.min(1.));
        }
        let coefficients = row_power(&model, 0, 0.5, work)?;
        for c in &mut argument.c[..=argument.degree] {
            *c = *c - I::point(0.5);
        }
        return Ok((
            power_contract(&[coefficients], model.shape[1], argument, work)?[0],
            error,
        ));
    }
    Ok((compose_target(p, phi, work)?, 0.))
}
struct TargetProof<'a> {
    patches: &'a [TargetPatch],
    maximum: f64,
    phase: f64,
}
fn mapped_cell(
    patches: &[Prepared],
    global: SurfaceBounds,
    target_proof: TargetProof<'_>,
    uv: Vec<H>,
    phi: Vec<H>,
    tol: f64,
    work: &mut Budget<'_>,
) -> Result<f64> {
    let TargetProof {
        patches: target,
        maximum: target_max,
        phase,
    } = target_proof;
    let mut stack: [Option<Node>; 26] = std::array::from_fn(|_| None);
    let (mut len, mut current, mut maximum) = (
        0_usize,
        Node {
            uv,
            target: phi,
            depth: 0,
        },
        0_f64,
    );
    loop {
        work.charge(1024)?;
        let box_ = uv_hull(&current.uv, work)?;
        let tbox = target_box(&current.target, work)?;
        let (mut si, mut se) = (0_usize, f64::INFINITY);
        let (mut ti, mut te) = (0_usize, f64::INFINITY);
        for (i, p) in patches.iter().enumerate() {
            let e = extension_error(p, box_, global, work)?;
            if e < se {
                si = i;
                se = e
            }
        }
        for (i, p) in target.iter().enumerate() {
            let e = target_collar(p, tbox, target_max, work)?;
            if e < te {
                ti = i;
                te = e
            }
        }
        let collar = I::point(se) + I::point(te);
        check(collar.valid(), ReadError::Unsupported)?;
        let allowed = (I::point(tol) - collar).lo;
        let mut accepted = None;
        if allowed > 0. {
            let partial = actual_partials(patches, box_, work)?;
            let du = derivative_hull(&current.uv, 0, work)?;
            let dv = derivative_hull(&current.uv, 1, work)?;
            let source: V = std::array::from_fn(|j| partial[0][j] * du + partial[1][j] * dv);
            // A partial finite Window supplies a positive affine correction A
            // of the nominal target map. The entire actual derivative image
            // of A(phi) lies inside this outward, natural-domain-clipped box.
            let expanded = tbox + I::new(-phase, phase);
            check(expanded.valid(), ReadError::Unsupported)?;
            let actual_box = I::new(expanded.lo.max(0.), expanded.hi.min(1.));
            check(actual_box.valid(), ReadError::Unsupported)?;
            let derivative = target_derivative(target, actual_box, work)?;
            let flow = (0..3).fold(I::point(0.), |sum, j| sum + source[j] * derivative[j]);
            check(flow.valid(), ReadError::Unsupported)?;
            // The private ORIGINAL proposal constructor earned phi'>0.
            // Every child restriction preserves it by a positive Jacobian.
            // Consequently this actual derivative-union dot, without a
            // polynomial degree101 product, earns the SAME IP2 direction.
            if flow.lo > 0. {
                let (a, proxy_error) =
                    position_composition(&patches[si], &current.uv, allowed, work)?;
                let (b, target_proxy_error) =
                    target_position(&target[ti], &current.target, allowed, work)?;
                let degree = a[0].degree.max(b[0].degree);
                let mut difference = [Poly::zero(degree)?; 3];
                let a = poly::elevate_many(a, degree, work)?;
                let b = poly::elevate_many(b, degree, work)?;
                for j in 0..3 {
                    difference[j] = a[j].add(b[j], true, work)?;
                }
                let proxy_allowance =
                    (I::point(allowed) - I::point(proxy_error) - I::point(target_proxy_error)).lo;
                check(proxy_allowance > 0., ReadError::Unsupported)?;
                let result = poly::bound(difference, 1., proxy_allowance, work);
                match result {
                    Ok(error) => {
                        let total = I::point(error)
                            + collar
                            + I::point(proxy_error)
                            + I::point(target_proxy_error);
                        if total.valid() && total.hi <= tol {
                            accepted = Some(total.hi)
                        }
                    }
                    Err(ReadError::Unsupported) => {}
                    Err(e) => return Err(e),
                }
            }
        }
        if let Some(error) = accepted {
            maximum = maximum.max(error);
            if len == 0 {
                return Ok(maximum);
            }
            len -= 1;
            current = stack[len].take().ok_or(ReadError::Invalid)?;
            continue;
        }
        check(current.depth < MAX_DEPTH, ReadError::Unsupported)?;
        let [u0, u1] = cv_split(&current.uv, I::point(0.5), work)?;
        let [t0, t1] = cv_split(&current.target, I::point(0.5), work)?;
        check(len < stack.len(), ReadError::Limit)?;
        stack[len] = Some(Node {
            uv: u1,
            target: t1,
            depth: current.depth + 1,
        });
        len += 1;
        current = Node {
            uv: u0,
            target: t0,
            depth: current.depth + 1,
        };
    }
}
#[allow(
    clippy::too_many_arguments,
    reason = "The two immutable Original packets, charts and separate actual derivative bounds define the coupled proof."
)]
fn spline(
    source: &Source<'_>,
    target_id: u32,
    target: &crate::brep::NurbsCurve,
    finite: (&super::super::uv::Definition, &Edge),
    decoded: (&crate::brep::NurbsCurve, &crate::brep::NurbsSurface),
    uv: &CurveCells<'_, '_, '_>,
    uv_knots: &[BigRational],
    patches: &[Prepared],
    maximum: SurfaceBounds,
    tol: f64,
    work: &mut Budget<'_>,
) -> Result<f64> {
    let target_packet =
        bezier::original_curve(source, target_id, target, Coordinates::Raw(3), work)?;
    check(
        target_packet.matches(source, target_id)
            && matches!(target_packet.coordinates(), Coordinates::Raw(3)),
        ReadError::Invalid,
    )?;
    let a = normalized_knots(uv_knots, work)?;
    let target_knots = raw_knots(source, target_id, 1, work)?;
    let b = normalized_knots(&target_knots[0], work)?;
    // Full-domain maps cover the whole bounded Original target. A proper
    // subwindow instead needs Original vertex/closest-point authority. The
    // present partial class has exact Original target natural domain [0,1].
    let zero = BigRational::from_integer(0.into());
    let one = BigRational::from_integer(1.into());
    let partial = uv_knots.first() != Some(&zero) || uv_knots.last() != Some(&one);
    let (mut map_window, mut phase, mut window_error) = ([0., 1.], 0_f64, 0_f64);
    if partial {
        check(
            target_knots[0].first() == Some(&zero) && target_knots[0].last() == Some(&one),
            ReadError::Unsupported,
        )?;
        let window = super::spline_window::prepare(
            source,
            target_id,
            finite.0,
            finite.1,
            &target_packet,
            work,
        )?;
        check(
            window.matches(source, target_id, finite.0, finite.1),
            ReadError::Invalid,
        )?;
        let ends = [uv_knots[0].clone(), uv_knots[uv_knots.len() - 1].clone()];
        let mut rounding = 0_f64;
        for i in 0..2 {
            work.charge(bigint_work(q_bits(&ends[i])?))?;
            map_window[i] = ends[i].to_f64().ok_or(ReadError::Unsupported)?;
            let delta = enclose(&ends[i], work)? - I::point(map_window[i]);
            check(delta.valid(), ReadError::Unsupported)?;
            rounding = rounding.max(delta.magnitude());
        }
        check(
            map_window[0] >= 0. && map_window[1] <= 1. && map_window[0] < map_window[1],
            ReadError::Unsupported,
        )?;
        let positive = window.raw_uv_width()
            - I::point(window.closest_delta_sum())
            - I::point(2.) * I::point(rounding);
        check(positive.valid() && positive.lo > 0., ReadError::Unsupported)?;
        let delta = I::point(window.closest_delta())
            + I::point(window.cached_window_delta())
            + I::point(rounding);
        check(delta.valid(), ReadError::Unsupported)?;
        phase = delta.hi;
        let error = I::point(window.speed()) * delta;
        check(error.valid(), ReadError::Unsupported)?;
        window_error = error.hi;
        // Window's normalized target numeric reserve is already included
        // once by the caller's curve_nurbs floor; do not add it again here.
    }
    let actual_tolerance = (I::point(tol) - I::point(window_error)).lo;
    check(actual_tolerance > 0., ReadError::Unsupported)?;
    // This initial class keeps every simultaneously live conservative Q
    // payload, both prepared grids and all polynomial temporaries inside the
    // declared 16MiB numerical arena. Broader knot counts require a separate
    // measured residency packet, rather than an assumed allocation bound.
    check(a.len() <= 65 && b.len() <= 65, ReadError::Unsupported)?;
    let (target_patches, target_max) = target_patches(&target_packet, &b, work)?;
    work.charge(a.len() * 64 + 4096)?;
    let mut midpoints = Vec::new();
    midpoints
        .try_reserve_exact(a.len() - 1)
        .map_err(|_| ReadError::Limit)?;
    for pair in a.windows(2) {
        midpoints.push(affine_q(
            &BigRational::new(1.into(), 2.into()),
            &[
                BigRational::from_integer(0.into()),
                BigRational::from_integer(1.into()),
            ],
            &[pair[0].clone(), pair[1].clone()],
            work,
        )?);
    }
    let cuts = union(&a, &midpoints, work)?;
    // The cubic proof handles crossings with actual target-knot collars;
    // no shared-knot snapping or artificial vanishing pullback is needed.
    work.charge(cuts.len() * 64 + 1024)?;
    let mut source_cuts = Vec::new();
    source_cuts
        .try_reserve_exact(cuts.len())
        .map_err(|_| ReadError::Limit)?;
    let mut mapped = Vec::new();
    mapped
        .try_reserve_exact(cuts.len())
        .map_err(|_| ReadError::Limit)?;
    check(
        (source_cuts.capacity() + mapped.capacity()) * std::mem::size_of::<f64>() <= 16 * 1024,
        ReadError::Limit,
    )?;
    for cut in &cuts {
        work.charge(bigint_work(q_bits(cut)?))?;
        let q = cut.to_f64().ok_or(ReadError::Unsupported)?;
        source_cuts.push(q);
        mapped.push(propose_parameter(
            q, decoded.0, decoded.1, target, map_window, work,
        )?);
    }
    check(
        mapped.first() == Some(&map_window[0])
            && mapped.last() == Some(&map_window[1])
            && mapped.windows(2).all(|p| p[0] < p[1]),
        ReadError::Unsupported,
    )?;
    let (mut ui, mut error) = (0_usize, 0_f64);
    for ((pair, map), samples) in cuts
        .windows(2)
        .zip(mapped.windows(2))
        .zip(source_cuts.windows(2))
    {
        let mut values = [0_f64; 2];
        for (i, fraction) in [1. / 3., 2. / 3.].into_iter().enumerate() {
            let q = samples[0] + fraction * (samples[1] - samples[0]);
            values[i] = propose_parameter(q, decoded.0, decoded.1, target, map_window, work)?;
        }
        work.charge(2048)?;
        let p0 = I::point(map[0]);
        let p3 = I::point(map[1]);
        let c0 = map[0];
        let c3 = map[1];
        let p1 = I::point(3. * values[0] - 1.5 * values[1] - (5. / 6.) * c0 + c3 / 3.);
        let p2 = I::point(3. * values[1] - 1.5 * values[0] + c0 / 3. - (5. / 6.) * c3);
        let parameters = [p0, p1, p2, p3];
        check(
            parameters
                .iter()
                .all(|p| p.valid() && p.lo >= 0. && p.hi <= 1.)
                && parameters.windows(2).all(|p| (p[1] - p[0]).lo > 0.),
            ReadError::Unsupported,
        )?;
        let mut phi = Vec::new();
        phi.try_reserve_exact(4).map_err(|_| ReadError::Limit)?;
        phi.extend(parameters.map(|q| [q, I::point(0.), I::point(0.), I::point(1.)]));
        let interval = restriction_indices(&a, &pair[0], &pair[1], &mut ui, work)?;
        let parent = uv.cells().get(ui).ok_or(ReadError::Unsupported)?;
        check(parent.knot_indices() == [ui, ui + 1], ReadError::Invalid)?;
        let controls = restricted(parent.controls(), &interval[0], &interval[1], work)?;
        let result = mapped_cell(
            patches,
            maximum,
            TargetProof {
                patches: &target_patches,
                maximum: target_max,
                phase,
            },
            controls,
            phi,
            actual_tolerance,
            work,
        );
        error = error.max(result?);
    }
    let total = I::point(error) + I::point(window_error);
    check(total.valid() && total.hi <= tol, ReadError::Unsupported)?;
    Ok(total.hi)
}

fn fixed_coordinate(
    source: &Source<'_>,
    id: u32,
    axis: usize,
    work: &mut Budget<'_>,
) -> Result<Option<BigRational>> {
    let field = fields(source, id, 1, work)?;
    let mut fixed = None;
    for p in list(field.controls)? {
        work.charge(128)?;
        let pid = reference(p)?;
        let [_, coordinates] = source.args(pid, "CARTESIAN_POINT")? else {
            return Err(ReadError::Invalid);
        };
        let coordinates = list(coordinates)?;
        check(coordinates.len() == 2, ReadError::Unsupported)?;
        let q = rational(source, pid, &coordinates[axis], work)?;
        if let Some(previous) = &fixed {
            if !q_cmp(previous, &q, work)?.is_eq() {
                return Ok(None);
            }
        } else {
            fixed = Some(q)
        }
    }
    Ok(fixed)
}
fn fixed_line_image(
    source: &Source<'_>,
    definition: &super::super::uv::Definition,
    original: &SurfaceCells<'_, '_>,
    knots: &[Vec<BigRational>],
    trace: &CurveCells<'_, '_, '_>,
    work: &mut Budget<'_>,
) -> Result<Option<[V; 2]>> {
    for (fixed, axis_knots) in knots.iter().enumerate() {
        let moving = 1 - fixed;
        if original.degrees()[moving] != 1 || original.counts()[moving] != 1 {
            continue;
        }
        let Some(value) = fixed_coordinate(source, definition.curve_id, fixed, work)? else {
            continue;
        };
        let mut at = None;
        for (i, knot) in axis_knots.iter().enumerate() {
            if q_cmp(knot, &value, work)?.is_eq() {
                at = Some(i);
                break;
            }
        }
        let Some(at) = at else { continue };
        let Some(patch) = original
            .patches()
            .iter()
            .find(|p| p.knot_indices()[fixed].contains(&at))
        else {
            continue;
        };
        let fixed_index = if patch.knot_indices()[fixed][0] == at {
            0
        } else {
            patch.shape()[fixed] - 1
        };
        work.charge(2048)?;
        check(
            unit(patch.controls()) && patch.shape()[moving] == 2,
            ReadError::Unsupported,
        )?;
        let row: [V; 2] = std::array::from_fn(|i| {
            let index = if moving == 0 {
                i * patch.shape()[1] + fixed_index
            } else {
                fixed_index * patch.shape()[1] + i
            };
            std::array::from_fn(|j| patch.controls()[index][j])
        });
        let mut sign = 0;
        for cell in trace.cells() {
            let d = derivative_hull(cell.controls(), moving, work)?;
            let current = if d.lo > 0. {
                1
            } else if d.hi < 0. {
                -1
            } else {
                return Err(ReadError::Unsupported);
            };
            check(sign == 0 || sign == current, ReadError::Unsupported)?;
            sign = current;
        }
        let first = trace
            .cells()
            .first()
            .ok_or(ReadError::Unsupported)?
            .controls()[0][moving];
        let last = *trace
            .cells()
            .last()
            .ok_or(ReadError::Unsupported)?
            .controls()
            .last()
            .ok_or(ReadError::Unsupported)?;
        let coordinates = [first, last[moving]];
        let points: [V; 2] =
            coordinates.map(|u| std::array::from_fn(|j| row[0][j] + (row[1][j] - row[0][j]) * u));
        return Ok(Some(points));
    }
    Ok(None)
}
#[allow(
    clippy::too_many_arguments,
    reason = "Original support, complete trace and finite target are independently bound."
)]
fn boundary_line(
    source: &Source<'_>,
    id: u32,
    definition: &super::super::uv::Definition,
    edge: &Edge,
    original: &SurfaceCells<'_, '_>,
    knots: &[Vec<BigRational>],
    trace: &CurveCells<'_, '_, '_>,
    tol: f64,
    work: &mut Budget<'_>,
) -> Result<f64> {
    let points = fixed_line_image(source, definition, original, knots, trace, work)?
        .ok_or(ReadError::Unsupported)?;
    let raw_id = unwrapped(source, id, work)?;
    let [_, _, direction_id] = source.args(raw_id, "LINE")? else {
        return Err(ReadError::Unsupported);
    };
    let direction = vector(source, reference(direction_id)?, 3, work)?;
    let flow = (0..3).fold(I::point(0.), |v, j| {
        v + (points[1][j] - points[0][j]) * direction[j]
    });
    check(flow.valid() && flow.lo > 0., ReadError::Unsupported)?;
    let target = I::point(analytic::curve_reserve(source, id, edge, work)?)
        + I::point(super::super::finite_axis::endpoint_reserve(
            source, id, edge, work,
        )?);
    check(target.valid() && target.hi < tol, ReadError::Unsupported)?;
    let available = (I::point(tol) - target).lo;
    let result = super::super::affine::nurbs(
        1,
        (&[0., 1.], &[2, 2]),
        None,
        &points,
        true,
        edge,
        available,
        work,
    )?;
    if let Proof::Certified(error) = result {
        let total = I::point(error) + target;
        check(total.valid() && total.hi <= tol, ReadError::Unsupported)?;
        return Ok(total.hi);
    }
    Err(ReadError::Unsupported)
}
fn raw_vertex(source: &Source<'_>, id: usize, work: &mut Budget<'_>) -> Result<V> {
    let id = u32::try_from(id).map_err(|_| ReadError::Unsupported)?;
    let [_, point_id] = source.args(id, "VERTEX_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let point_id = reference(point_id)?;
    numeric_preflight(source, point_id, work)?;
    point(source, point_id, 3, work)
}
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |sum, j| sum + a[j] * b[j])
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[allow(
    clippy::too_many_arguments,
    reason = "The actual Source support, UV natural image and finite master remain independently bound."
)]
fn bounded_line_target(
    source: &Source<'_>,
    id: u32,
    definition: &super::super::uv::Definition,
    edge: &Edge,
    original: &SurfaceCells<'_, '_>,
    knots: &[Vec<BigRational>],
    trace: &CurveCells<'_, '_, '_>,
    tol: f64,
    work: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    let Curve::Nurbs {
        curve,
        closed: false,
    } = &edge.curve
    else {
        return Ok(None);
    };
    if curve.degree() != 1 || curve.points().len() != 2 {
        return Ok(None);
    }
    let Some(image) = fixed_line_image(source, definition, original, knots, trace, work)? else {
        return Ok(None);
    };
    let data = bezier::original_curve(source, id, curve, Coordinates::Raw(3), work)?;
    check(
        data.matches(source, id) && matches!(data.coordinates(), Coordinates::Raw(3)),
        ReadError::Invalid,
    )?;
    check(
        data.cells().len() == 1 && data.degree() == 1,
        ReadError::Unsupported,
    )?;
    let controls = data.cells()[0].controls();
    check(
        controls.len() == 2 && unit(controls),
        ReadError::Unsupported,
    )?;
    work.charge(8192)?;
    let target: [V; 2] = std::array::from_fn(|i| [controls[i][0], controls[i][1], controls[i][2]]);
    let direction: V = std::array::from_fn(|j| target[1][j] - target[0][j]);
    let source_direction: V = std::array::from_fn(|j| image[1][j] - image[0][j]);
    let den = dot(direction, direction);
    let source_den = dot(source_direction, source_direction);
    check(
        den.valid() && den.lo > 0. && source_den.valid() && source_den.lo > 0.,
        ReadError::Unsupported,
    )?;
    let vertices = [
        raw_vertex(source, edge.vertices[0], work)?,
        raw_vertex(source, edge.vertices[1], work)?,
    ];
    let mut parameters = [I::point(0.); 2];
    for i in 0..2 {
        let delta = std::array::from_fn(|j| vertices[i][j] - target[0][j]);
        let q = dot(delta, direction) / den;
        check(q.valid(), ReadError::Unsupported)?;
        // A clamped unit-weight degree1 curve is an actual finite segment.
        // Its unique ORIGINAL closest parameter is the projected coordinate
        // clamped to[0,1], not the decoder's proposed edge.range.
        parameters[i] = I::new(q.lo.clamp(0., 1.), q.hi.clamp(0., 1.));
    }
    check(
        (parameters[1] - parameters[0]).lo > 0.,
        ReadError::Unsupported,
    )?;
    let used: [V; 2] = parameters.map(|q| std::array::from_fn(|j| target[0][j] + q * direction[j]));
    let mut selected = [0_f64; 2];
    let mut error = 0_f64;
    let mut lower = 0_f64;
    let source_length = norm(source_direction);
    check(
        source_length.valid() && source_length.lo > 0.,
        ReadError::Unsupported,
    )?;
    for i in 0..2 {
        let delta: V = std::array::from_fn(|j| used[i][j] - image[0][j]);
        let q = dot(delta, source_direction) / source_den;
        check(q.valid(), ReadError::Unsupported)?;
        selected[i] = q.midpoint().clamp(0., 1.);
        let witness = if q.hi < 0. {
            norm(delta)
        } else if q.lo > 1. {
            norm(std::array::from_fn(|j| used[i][j] - image[1][j]))
        } else {
            norm(cross(delta, source_direction)) / I::point(source_length.hi)
        };
        check(witness.valid(), ReadError::Unsupported)?;
        lower = lower.max(witness.lo);
        let matched: V =
            std::array::from_fn(|j| image[0][j] + I::point(selected[i]) * source_direction[j]);
        let bound = norm(std::array::from_fn(|j| used[i][j] - matched[j]));
        check(bound.valid(), ReadError::Unsupported)?;
        error = error.max(bound.hi);
    }
    if lower > tol {
        return Ok(Some(Proof::Disagreement(lower)));
    }
    let flow = dot(direction, source_direction);
    check(
        flow.valid()
            && flow.lo > 0.
            && selected.iter().all(|q| q.is_finite())
            && selected[0] < selected[1],
        ReadError::Unsupported,
    )?;
    let range = curve.range();
    check(
        range[0] < range[1] && edge.range[0] >= range[0] && edge.range[1] <= range[1],
        ReadError::Unsupported,
    )?;
    let span = I::point(range[1]) - I::point(range[0]);
    check(span.valid() && span.lo > 0., ReadError::Unsupported)?;
    let mut discrepancy = 0_f64;
    for (&end, &parameter) in edge.range.iter().zip(&parameters) {
        let cached = (I::point(end) - I::point(range[0])) / span;
        check(cached.valid(), ReadError::Unsupported)?;
        discrepancy = discrepancy.max((parameter - cached).magnitude());
    }
    let transfer = norm(direction) * I::point(discrepancy);
    let total = I::point(error) + transfer;
    check(total.valid(), ReadError::Unsupported)?;
    // Both used 3D and selected UV images are straight finite segments.
    // Their increasing endpoint interpolation covers EVERY point. The UV
    // natural curve may over-cover; only this earned subsegment is used.
    work.charge(64)?;
    Ok(Some(if total.hi <= tol {
        Proof::Certified(total.hi)
    } else {
        Proof::Unsupported
    }))
}

/// The proposed common parameter is a positive affine normalization of each
/// COMPLETE original natural domain. Coupled residuals earn that proposal;
/// endpoint equality or a cached finite-edge range never substitutes for it.
pub(super) fn certify(
    source: &Source<'_>,
    id: u32,
    definition: &super::super::uv::Definition,
    support: &Surface,
    edge: &Edge,
    tol: f64,
    work: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    match super::finite_boundary::certify(source, id, definition, support, edge, tol, work) {
        Ok(Some(Proof::Certified(error))) => return Ok(Some(Proof::Certified(error))),
        Ok(Some(Proof::Disagreement(lower))) => return Ok(Some(Proof::Disagreement(lower))),
        Ok(_) | Err(ReadError::Unsupported) => {}
        Err(e) => return Err(e),
    }
    let (
        Surface::Nurbs { surface, closed },
        Curve::Nurbs {
            curve,
            closed: false,
        },
    ) = (support, &definition.geometry)
    else {
        return Ok(None);
    };
    match &edge.curve {
        Curve::Nurbs {
            curve,
            closed: false,
        } if curve.degree() <= 11 && curve.points().len() <= 1024 => {}
        Curve::Line { .. } => {}
        _ => return Ok(None),
    }
    if closed.iter().any(|v| *v)
        || curve.degree() > 11
        || curve.points().len() > 1024
        || surface.degrees().iter().sum::<usize>() * curve.degree() > 80
    {
        return Ok(None);
    }
    // Three bounded extraction arenas, decoded arguments and all polynomial,
    // rational-cut and adaptive-stack temporaries coexist within this envelope.
    const ENVELOPE: usize = 3 * 8 * 1024 * 1024 + 4 * 1024 * 1024 + NUMERIC_BYTES + 256 * 1024;
    check(ENVELOPE <= 64 * 1024 * 1024, ReadError::Limit)?;
    work.charge(4096)?;
    for axis in 0..2 {
        let (_, mult) = surface.knots(axis);
        let p = surface.degrees()[axis];
        // C1 across every original interior knot makes the actual derivative
        // well-defined at a join. Its union still includes both adjacent jets.
        check(
            mult.len() >= 2 && mult[1..mult.len() - 1].iter().all(|m| *m < p),
            ReadError::Unsupported,
        )?;
    }
    let original = surface_bezier::original_surface(source, definition.surface, surface, work)?;
    check(
        original.matches(source, definition.surface) && original.degrees() == surface.degrees(),
        ReadError::Invalid,
    )?;
    check(
        original.peak_bytes() + NUMERIC_BYTES + 4 * 1024 * 1024 <= 64 * 1024 * 1024
            && original.retained_bytes() <= 8 * 1024 * 1024,
        ReadError::Limit,
    )?;
    let knots = raw_knots(source, definition.surface, 2, work)?;
    let min = [knots[0][0].clone(), knots[1][0].clone()];
    let max = [
        knots[0][knots[0].len() - 1].clone(),
        knots[1][knots[1].len() - 1].clone(),
    ];
    let ends = [
        [min[0].clone(), max[0].clone()],
        [min[1].clone(), max[1].clone()],
    ];
    let (patches, maximum) = prepare(&original, &knots, &ends, work)?;
    let trace = bezier::original_curve(
        source,
        definition.curve_id,
        curve,
        Coordinates::UnitBox2 {
            min: &min,
            max: &max,
        },
        work,
    )?;
    check(
        trace.matches(source, definition.curve_id),
        ReadError::Invalid,
    )?;
    match trace.coordinates() {
        Coordinates::UnitBox2 { min: a, max: b }
            if std::ptr::eq(a, &min) && std::ptr::eq(b, &max) => {}
        _ => return Err(ReadError::Invalid),
    }
    check(
        trace.retained_bytes() <= 8 * 1024 * 1024
            && trace.peak_bytes() + original.retained_bytes() + NUMERIC_BYTES + 4 * 1024 * 1024
                <= 64 * 1024 * 1024,
        ReadError::Limit,
    )?;
    let uvknots = raw_knots(source, definition.curve_id, 1, work)?;
    let bounded = bounded_line_target(
        source, id, definition, edge, &original, &knots, &trace, tol, work,
    )?;
    if let Some(Proof::Disagreement(lower)) = bounded {
        return Ok(Some(Proof::Disagreement(lower)));
    }
    // Original->decoded Body transfer is distinct from the RAW composition.
    // Its normalized-domain bound is added once, without treating decoded
    // coefficient identity as equality of the original decimal geometry.
    let field = fields(source, definition.surface, 2, work)?;
    let body = nurbs_error(
        source,
        definition.surface,
        &field,
        &surface.degrees(),
        &[surface.knots(0), surface.knots(1)],
        surface.points(),
        surface.weights(),
        3,
        None,
        work,
    )?
    .error;
    let target_floor = if matches!(edge.curve, Curve::Nurbs { .. }) {
        curve_numeric_preflight(source, id, work)?;
        curve_nurbs(source, id, &edge.curve, 3, work)?.error
    } else {
        0.
    };
    let floor = I::point(body) + I::point(target_floor);
    check(floor.valid() && floor.hi < tol, ReadError::Unsupported)?;
    let available = (I::point(tol) - floor).lo;
    if let Some(proof) = bounded {
        return match proof {
            Proof::Certified(error) => {
                let total = I::point(error) + floor;
                check(total.valid() && total.hi <= tol, ReadError::Unsupported)?;
                work.charge(64)?;
                Ok(Some(Proof::Certified(total.hi)))
            }
            other => Ok(Some(other)),
        };
    }
    let result = match &edge.curve {
        Curve::Nurbs { curve: target, .. } => spline(
            source,
            id,
            target,
            (definition, edge),
            (curve, surface),
            &trace,
            &uvknots[0],
            &patches,
            maximum,
            available,
            work,
        )?,
        Curve::Line { .. } => boundary_line(
            source, id, definition, edge, &original, &knots, &trace, available, work,
        )?,
        _ => return Ok(None),
    };
    let total = I::point(result) + floor;
    check(total.valid() && total.hi <= tol, ReadError::Unsupported)?;
    work.charge(64)?;
    Ok(Some(Proof::Certified(total.hi)))
}
