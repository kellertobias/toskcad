//! Original common-parameter subdivision and source-bound finite admission.
use super::super::super::Proof;
use super::super::bezier::{self, Coordinates, CurveCells};
use super::super::spline_window;
use super::*;
const MAX_KNOTS: usize = 129;
const MAX_UNION: usize = 256;
const SCRATCH: usize = 128 * 1024;
fn bits(q: &BigRational) -> Result<usize> {
    let n = q.numer().bits().max(q.denom().bits());
    check(n <= 4096, ReadError::Unsupported)?;
    usize::try_from(n).map_err(|_| ReadError::Limit)
}
fn compare(a: &BigRational, c: &BigRational, b: &mut Budget<'_>) -> Result<std::cmp::Ordering> {
    let size = bits(a)?.saturating_add(bits(c)?).saturating_add(2);
    // Structural equality alone proves equality even without relying on
    // canonical reduction. BigInt digit equality is linear, not a cross-product.
    b.charge(size.saturating_mul(4).saturating_add(256))?;
    if a.numer() == c.numer() && a.denom() == c.denom() {
        return Ok(std::cmp::Ordering::Equal);
    }
    b.charge(bigint_work(
        bits(a)?.saturating_add(bits(c)?).saturating_add(2),
    ))?;
    Ok(a.cmp(c))
}
fn operation(
    a: &BigRational,
    c: &BigRational,
    divide: bool,
    b: &mut Budget<'_>,
) -> Result<BigRational> {
    let n = bits(a)?.saturating_add(bits(c)?).saturating_add(2);
    check(n <= 4096, ReadError::Unsupported)?;
    b.charge(bigint_work(n))?;
    check(!divide || !c.is_zero(), ReadError::Unsupported)?;
    let value = if divide { a / c } else { a - c };
    bits(&value)?;
    Ok(value)
}
fn knots(source: &Source<'_>, id: u32, b: &mut Budget<'_>) -> Result<Vec<BigRational>> {
    let field = fields(source, id, 1, b)?;
    let Some((_, values)) = field.knots[0] else {
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
    check(out.capacity() <= MAX_KNOTS, ReadError::Limit)?;
    for v in values {
        b.charge(128)?;
        let q = rational(source, id, v, b)?;
        bits(&q)?;
        if let Some(previous) = out.last() {
            check(compare(previous, &q, b)?.is_lt(), ReadError::Unsupported)?;
        }
        out.push(q);
    }
    Ok(out)
}
fn merged_knots(
    a: &[BigRational],
    c: &[BigRational],
    b: &mut Budget<'_>,
) -> Result<Vec<BigRational>> {
    let mut out = Vec::new();
    out.try_reserve_exact(MAX_UNION)
        .map_err(|_| ReadError::Limit)?;
    check(out.capacity() <= MAX_UNION, ReadError::Limit)?;
    check(
        compare(&a[0], &c[0], b)?.is_ge()
            && compare(
                a.last().ok_or(ReadError::Unsupported)?,
                c.last().ok_or(ReadError::Unsupported)?,
                b,
            )?
            .is_le(),
        ReadError::Unsupported,
    )?;
    for value in a.iter().chain(c) {
        if compare(value, &a[0], b)?.is_lt()
            || compare(value, a.last().ok_or(ReadError::Unsupported)?, b)?.is_gt()
        {
            continue;
        }
        let mut low = 0;
        let mut high = out.len();
        let mut duplicate = false;
        while low < high {
            let mid = low + (high - low) / 2;
            match compare(&out[mid], value, b)? {
                std::cmp::Ordering::Less => low = mid + 1,
                std::cmp::Ordering::Greater => high = mid,
                std::cmp::Ordering::Equal => {
                    duplicate = true;
                    break;
                }
            }
        }
        if !duplicate {
            check(out.len() < MAX_UNION, ReadError::Unsupported)?;
            b.charge(
                out.len()
                    .saturating_sub(low)
                    .saturating_mul(256)
                    .saturating_add(4096),
            )?;
            out.insert(low, value.clone());
        }
    }
    Ok(out)
}
fn fraction(
    q: &BigRational,
    a: &BigRational,
    c: &BigRational,
    b: &mut Budget<'_>,
) -> Result<BigRational> {
    let top = operation(q, a, false, b)?;
    let bottom = operation(c, a, false, b)?;
    operation(&top, &bottom, true, b)
}
fn inside(q: &BigRational, b: &mut Budget<'_>) -> Result<I> {
    let zero = BigRational::from_integer(0.into());
    let one = BigRational::from_integer(1.into());
    check(
        compare(q, &zero, b)?.is_ge() && compare(q, &one, b)?.is_le(),
        ReadError::Unsupported,
    )?;
    let i = enclose(q, b)?;
    let i = I::new(i.lo.max(0.), i.hi.min(1.));
    check(i.valid(), ReadError::Unsupported)?;
    Ok(i)
}
fn restrict(
    data: &CurveCells<'_, '_, '_>,
    knots: &[BigRational],
    from: &BigRational,
    to: &BigRational,
    cursor: &mut usize,
    b: &mut Budget<'_>,
) -> Result<[P; 3]> {
    while *cursor < data.cells().len() {
        let cell = &data.cells()[*cursor];
        b.charge(64)?;
        let [i, j] = cell.knot_indices();
        let a = knots.get(i).ok_or(ReadError::Unsupported)?;
        let c = knots.get(j).ok_or(ReadError::Unsupported)?;
        if compare(from, c, b)?.is_ge() {
            *cursor += 1;
            continue;
        }
        check(
            compare(from, a, b)?.is_ge() && compare(to, c, b)?.is_le(),
            ReadError::Unsupported,
        )?;
        let left = fraction(from, a, c, b)?;
        let right = fraction(to, a, c, b)?;
        check(compare(&left, &right, b)?.is_lt(), ReadError::Unsupported)?;
        let right_i = inside(&right, b)?;
        let left_ratio = if left.is_zero() {
            BigRational::from_integer(0.into())
        } else {
            operation(&left, &right, true, b)?
        };
        let left_i = inside(&left_ratio, b)?;
        let mut out = [P::constant(I::point(0.), b)?; 3];
        check(
            cell.controls().len() == data.degree() + 1
                && cell
                    .controls()
                    .iter()
                    .all(|h| h[3].lo == 1. && h[3].hi == 1.),
            ReadError::Unsupported,
        )?;
        for axis in 0..3 {
            let mut values = [I::point(0.); 17];
            for (i, h) in cell.controls().iter().enumerate() {
                b.charge(16)?;
                values[i] = h[axis];
            }
            let p = P::new(&values[..=data.degree()], b)?;
            let p = if compare(to, c, b)?.is_eq() {
                p
            } else {
                p.split(right_i, b)?[0]
            };
            out[axis] = if left.is_zero() {
                p
            } else {
                p.split(left_i, b)?[1]
            };
        }
        return Ok(out);
    }
    Err(ReadError::Unsupported)
}

pub(in super::super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    d: &super::super::super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    b: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    let (
        Surface::Cylinder { .. },
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
    if uv.degree() > 3
        || target.degree() > 4
        || uv.points().len() > 1024
        || target.points().len() > 1024
    {
        return Ok(None);
    }
    b.charge(8192)?;
    check(tol.is_finite() && tol > 0., ReadError::Invalid)?;
    check(
        super::super::super::si_radians(source, b)?,
        ReadError::Unsupported,
    )?;
    let id = unwrapped(source, curve_id, b)?;
    spline_window::numeric_preflight(source, d.curve_id, b)?;
    spline_window::numeric_preflight(source, id, b)?;
    // Raw p<=4, <=1024 input controls and <=128 original spans have
    // <=513 decomposed controls/640 final cell controls. Old/new capacities,
    // tags, original numeric-token temporaries and fixed scratch fit 1MiB
    // per constructor, including conservative spare capacities.
    let maximum_q = (2 * MAX_KNOTS + MAX_UNION)
        .checked_mul(std::mem::size_of::<BigRational>() + 4096)
        .ok_or(ReadError::Limit)?;
    check(
        maximum_q
            .saturating_add(2 * 1024 * 1024)
            .saturating_add(SCRATCH)
            .saturating_add(spline_window::EXTRA_SCRATCH_BYTES)
            <= 8 * 1024 * 1024,
        ReadError::Limit,
    )?;
    let ka = knots(source, d.curve_id, b)?;
    let kc = knots(source, id, b)?;

    let a = bezier::original_curve(source, d.curve_id, uv, Coordinates::Raw(2), b)?;
    let c = bezier::original_curve(source, id, target, Coordinates::Raw(3), b)?;
    check(
        a.matches(source, d.curve_id)
            && c.matches(source, id)
            && matches!(a.coordinates(), Coordinates::Raw(2))
            && matches!(c.coordinates(), Coordinates::Raw(3)),
        ReadError::Unsupported,
    )?;

    check(
        a.peak_bytes() <= 1024 * 1024 && c.peak_bytes() <= 1024 * 1024,
        ReadError::Limit,
    )?;
    let points = merged_knots(&ka, &kc, b)?;
    let qbytes = (ka.capacity() + kc.capacity() + points.capacity())
        .checked_mul(std::mem::size_of::<BigRational>() + 4096)
        .ok_or(ReadError::Limit)?;
    let construction = a
        .peak_bytes()
        .max(a.retained_bytes().saturating_add(c.peak_bytes()));
    let resident = a
        .retained_bytes()
        .saturating_add(c.retained_bytes())
        .saturating_add(qbytes)
        .saturating_add(SCRATCH)
        .saturating_add(spline_window::EXTRA_SCRATCH_BYTES);
    check(
        construction.max(resident) <= 8 * 1024 * 1024,
        ReadError::Limit,
    )?;
    let window = spline_window::prepare(source, id, d, edge, &c, b)?;
    check(window.matches(source, id, d, edge), ReadError::Invalid)?;
    let domain = window.raw_uv_domain();
    for (i, q) in [&ka[0], ka.last().ok_or(ReadError::Unsupported)?]
        .into_iter()
        .enumerate()
    {
        let original = enclose(q, b)?;
        check(
            domain[i].valid() && domain[i].lo <= original.lo && domain[i].hi >= original.hi,
            ReadError::Invalid,
        )?;
    }
    let [_, placement, r] = source.args(d.surface, "CYLINDRICAL_SURFACE")? else {
        return Err(ReadError::Invalid);
    };
    let frame = frame(source, reference(placement)?, 3, b)?;
    let radius = scalar(source, d.surface, r, b)?;
    let phase = I::point(window.closest_delta()) + I::point(window.cached_window_delta());
    let width_error = I::point(window.closest_delta_sum())
        + I::point(2.) * I::point(window.cached_window_delta());
    let width = window.raw_uv_width();
    check(width.valid() && width.lo > 0., ReadError::Unsupported)?;
    let ratio = width_error / width;
    let speed_factor = I::point(1.) - ratio;
    check(
        ratio.valid() && speed_factor.valid() && speed_factor.lo > 0.,
        ReadError::Unsupported,
    )?;
    let reserve = I::point(window.speed()) * phase + I::point(window.numeric_reserve());
    let mut maximum = 0f64;
    let (mut source_cursor, mut target_cursor) = (0usize, 0usize);
    for pair in points.windows(2) {
        b.charge(4096)?;
        let p = restrict(&a, &ka, &pair[0], &pair[1], &mut source_cursor, b)?;
        let q = restrict(&c, &kc, &pair[0], &pair[1], &mut target_cursor, b)?;
        let proof = cell(&frame, radius, &[p[0], p[1]], &q, b)?;
        let span = enclose(&operation(&pair[1], &pair[0], false, b)?, b)?;
        let shifted = I::point(window.acceleration()) * phase * span
            + I::point(window.speed()) * span * ratio;
        let first = I::point(proof.derivative) + shifted;
        let global_speed = I::point(window.gamma()) * span * speed_factor;
        let local_speed = I::point(proof.speed) - shifted;
        check(
            global_speed.valid() && local_speed.valid(),
            ReadError::Unsupported,
        )?;
        let speed = I::point(global_speed.lo.max(local_speed.lo));
        let position = I::point(proof.position) + reserve;
        if !first.valid()
            || !speed.valid()
            || !position.valid()
            || first.hi >= speed.lo
            || position.hi > tol
        {
            return Ok(Some(Proof::Unsupported));
        }
        maximum = maximum.max(position.hi);
    }
    b.charge(64)?;
    Ok(Some(Proof::Certified(maximum)))
}
