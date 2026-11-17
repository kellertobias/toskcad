//! A finite original isocurve on a genuine clamped tensor boundary.
//! Whole original coefficient/jet bounds and an earned finite window, rather
//! than a cached endpoint or a nominal surface role, authorize this proof.
use super::bezier::{self, Coordinates};
use super::spline_window;
use super::*;
use crate::step_read::brep::check;
const MAX_CONTROLS: usize = 1024;
const LOCAL_BYTES: usize = 512 * 1024;
const EXTRACTION_BYTES: usize = 2 * 1024 * 1024;
// This envelope is input-independent. Enforce it when compiling rather than
// carrying an unreachable resource-failure branch through every query.
const _: () =
    assert!(EXTRACTION_BYTES + LOCAL_BYTES + spline_window::EXTRA_SCRATCH_BYTES <= 8 * 1024 * 1024);

fn equal(
    source: &Source<'_>,
    a: (u32, &Value<'_>),
    b: (u32, &Value<'_>),
    work: &mut Budget<'_>,
) -> Result<()> {
    check(exact_same(source, a, b, work)?, ReadError::Unsupported)
}
fn unit_weights(field: &Fields<'_, '_>) -> Result<()> {
    // The initial finite-boundary route is polynomial. No rational-speed or
    // rational tensor restriction is silently substituted for this contract.
    check(field.weights.is_none(), ReadError::Unsupported)
}
fn original_uv(
    source: &Source<'_>,
    d: &super::super::uv::Definition,
    work: &mut Budget<'_>,
) -> Result<Iso> {
    let iso = isoline(source, d, work)?;
    check(
        iso.moving == 1 && iso.rate.valid() && iso.rate.lo > 0.,
        ReadError::Unsupported,
    )?;
    let f = fields(source, d.curve_id, 1, work)?;
    unit_weights(&f)?;
    let Some((m, k)) = f.knots[0] else {
        return Err(ReadError::Unsupported);
    };
    let m = list(m)?;
    let k = list(k)?;
    check(m.len() == 2 && k.len() == 2, ReadError::Unsupported)?;
    for (i, id) in [
        iso.origin_id,
        iso.endpoint_id.ok_or(ReadError::Unsupported)?,
    ]
    .into_iter()
    .enumerate()
    {
        let [_, v] = source.args(id, "CARTESIAN_POINT")? else {
            return Err(ReadError::Invalid);
        };
        let v = list(v)?;
        check(v.len() == 2, ReadError::Unsupported)?;
        // Unit-weight p1 and these ORIGINAL endpoint equalities give V(t)=t,
        // not merely the same finite segment with an unknown parameter speed.
        equal(source, (id, &v[1]), (d.curve_id, &k[i]), work)?;
    }
    Ok(iso)
}
fn coefficient_errors(
    source: &Source<'_>,
    surface_id: u32,
    s: &crate::brep::NurbsSurface,
    id: u32,
    target: &crate::brep::NurbsCurve,
    iso: &Iso,
    work: &mut Budget<'_>,
) -> Result<(f64, f64)> {
    let sf = fields(source, surface_id, 2, work)?;
    let tf = fields(source, id, 1, work)?;
    unit_weights(&sf)?;
    unit_weights(&tf)?;
    let rows = list(sf.controls)?;
    let tc = list(tf.controls)?;
    check(
        rows.len() == 2 && tc.len() == target.points().len() && tc.len() <= MAX_CONTROLS,
        ReadError::Unsupported,
    )?;
    let range = s.range(0);
    let end = if iso.origin[0].lo <= range[0] && iso.origin[0].hi >= range[0] {
        0
    } else if iso.origin[0].lo <= range[1] && iso.origin[0].hi >= range[1] {
        1
    } else {
        return Err(ReadError::Unsupported);
    };
    // Cache overlap only proposes a boundary. Both actual UV endpoints must
    // equal the actual fixed surface knot before any row is selected.
    boundary(source, surface_id, &sf, 0, end, iso, work)?;
    let row = list(&rows[end])?;
    check(row.len() == tc.len(), ReadError::Unsupported)?;
    let Some((sm, sk)) = sf.knots[1] else {
        return Err(ReadError::Unsupported);
    };
    let Some((tm, tk)) = tf.knots[0] else {
        return Err(ReadError::Unsupported);
    };
    let (sm, sk, tm, tk) = (list(sm)?, list(sk)?, list(tm)?, list(tk)?);
    check(
        sk.len() == tk.len() && sm.len() == tm.len() && sm.len() == sk.len(),
        ReadError::Unsupported,
    )?;
    // The finite UV segment must belong to the ORIGINAL natural surface
    // domain, including a sub-ULP endpoint shadow. A small extrapolation
    // error is not permission to use an undefined surface parameter.
    for (i, point_id) in [
        iso.origin_id,
        iso.endpoint_id.ok_or(ReadError::Unsupported)?,
    ]
    .into_iter()
    .enumerate()
    {
        let [_, values] = source.args(point_id, "CARTESIAN_POINT")? else {
            return Err(ReadError::Invalid);
        };
        let values = list(values)?;
        let endpoint = rational(source, point_id, &values[1], work)?;
        let bound_index = if i == 0 { 0 } else { sk.len() - 1 };
        let bound = rational(source, surface_id, &sk[bound_index], work)?;
        let bits = endpoint
            .numer()
            .bits()
            .max(endpoint.denom().bits())
            .saturating_add(bound.numer().bits().max(bound.denom().bits()))
            .saturating_add(2);
        work.charge(bigint_work(
            usize::try_from(bits).map_err(|_| ReadError::Limit)?,
        ))?;
        check(
            if i == 0 {
                endpoint >= bound
            } else {
                endpoint <= bound
            },
            ReadError::Unsupported,
        )?;
    }
    for i in 0..sk.len() {
        equal(source, (surface_id, &sk[i]), (id, &tk[i]), work)?;
        equal(source, (surface_id, &sm[i]), (id, &tm[i]), work)?;
    }
    let p = target.degree();
    check(
        s.degrees()[1] == p && s.knots(0).1 == [2, 2] && (1..=11).contains(&p),
        ReadError::Unsupported,
    )?;
    let mut expanded = [I::point(0.); MAX_CONTROLS + 12];
    let mut count = 0_usize;
    for (i, &mult) in target.knots().1.iter().enumerate() {
        work.charge(mult.saturating_mul(8).saturating_add(128))?;
        integer(source, id, &tm[i], mult, work)?;
        let knot = scalar(source, id, &tk[i], work)?;
        check(
            count.checked_add(mult).is_some_and(|n| n <= expanded.len()),
            ReadError::Limit,
        )?;
        expanded[count..count + mult].fill(knot);
        count += mult;
    }
    check(count == tc.len() + p + 1, ReadError::Unsupported)?;
    let bytes = tc
        .len()
        .checked_mul(std::mem::size_of::<V>())
        .ok_or(ReadError::Limit)?;
    check(bytes <= 64 * 1024, ReadError::Limit)?;
    work.charge(tc.len().saturating_mul(32).saturating_add(512))?;
    let mut differences = Vec::new();
    differences
        .try_reserve_exact(tc.len())
        .map_err(|_| ReadError::Limit)?;
    check(
        differences.capacity() * std::mem::size_of::<V>() <= 64 * 1024,
        ReadError::Limit,
    )?;
    let mut position = 0_f64;
    for (a, b) in row.iter().zip(tc) {
        work.charge(512)?;
        let (a, b) = (
            point(source, reference(a)?, 3, work)?,
            point(source, reference(b)?, 3, work)?,
        );
        let delta: V = std::array::from_fn(|j| a[j] - b[j]);
        let length = norm(delta);
        check(length.valid(), ReadError::Unsupported)?;
        position = position.max(length.hi);
        differences.push(delta);
    }
    let mut derivative = 0_f64;
    for (i, ab) in differences.windows(2).enumerate() {
        work.charge(512)?;
        let h = expanded[i + p + 1] - expanded[i + 1];
        check(h.valid() && h.lo > 0., ReadError::Unsupported)?;
        let delta = std::array::from_fn(|j| I::point(p as f64) * (ab[1][j] - ab[0][j]) / h);
        let length = norm(delta);
        check(length.valid(), ReadError::Unsupported)?;
        derivative = derivative.max(length.hi);
    }
    Ok((position, derivative))
}
fn prove(
    source: &Source<'_>,
    curve_id: u32,
    d: &super::super::uv::Definition,
    support: &Surface,
    edge: &Edge,
    tol: f64,
    work: &mut Budget<'_>,
) -> Result<Option<super::super::Proof>> {
    let (
        Surface::Nurbs { surface: s, closed },
        Curve::Nurbs {
            curve: uv,
            closed: false,
        },
        Curve::Nurbs {
            curve: target,
            closed: false,
        },
    ) = (support, &d.geometry, &edge.curve)
    else {
        return Ok(None);
    };
    if s.weights().is_some()
        || uv.weights().is_some()
        || target.weights().is_some()
        || closed.iter().any(|x| *x)
        || s.degrees()[0] != 1
        || s.counts()[0] != 2
        || uv.degree() != 1
        || uv.points().len() != 2
        || target.degree() != s.degrees()[1]
        || target.points().len() != s.counts()[1]
        || target.points().len() > MAX_CONTROLS
    {
        return Ok(None);
    }
    // Eligibility only: already-complete natural rows belong to the older
    // full-boundary provider. Cache equality here can only decline this
    // extra route; it never grants original parameter membership.
    if uv.range() == target.range()
        || uv.points()[0][0] != uv.points()[1][0]
        || uv.points()[0][1] >= uv.points()[1][1]
    {
        return Ok(None);
    }
    work.charge(4096)?;
    check(tol.is_finite() && tol > 0., ReadError::Invalid)?;
    let id = unwrapped(source, curve_id, work)?;
    for record in [id, d.curve_id, d.surface] {
        spline_window::numeric_preflight(source, record, work)?;
    }
    // Bound all extraction overlap before its first allocation. With p<=11,
    // <=1024 input poles and <=128 spans, original_curve's fixed caps use
    // <2MiB (source net/growing net, tags, <=1536 final controls and scratch).
    check(
        (1..=11).contains(&target.degree()) && target.knots().0.len() <= 129,
        ReadError::Unsupported,
    )?;
    let iso = original_uv(source, d, work)?;
    let sf = fields(source, d.surface, 2, work)?;
    let body = nurbs_error(
        source,
        d.surface,
        &sf,
        &s.degrees(),
        &[s.knots(0), s.knots(1)],
        s.points(),
        s.weights(),
        3,
        None,
        work,
    )?
    .error;
    let data = bezier::original_curve(source, id, target, Coordinates::Raw(3), work)?;
    check(
        data.matches(source, id) && matches!(data.coordinates(), Coordinates::Raw(3)),
        ReadError::Invalid,
    )?;
    check(
        data.peak_bytes() <= EXTRACTION_BYTES
            && data.retained_bytes() + LOCAL_BYTES + spline_window::EXTRA_SCRATCH_BYTES
                <= 8 * 1024 * 1024,
        ReadError::Limit,
    )?;
    let window = spline_window::prepare_normalized(source, id, d, edge, &data, work)?;
    check(window.matches(source, id, d, edge), ReadError::Invalid)?;
    let (position, first) = coefficient_errors(source, d.surface, s, id, target, &iso, work)?;
    let phase = I::point(window.closest_delta()) + I::point(window.cached_window_delta());
    let ratio = (I::point(window.closest_delta_sum())
        + I::point(2.) * I::point(window.cached_window_delta()))
        / window.raw_uv_width();
    let speed = I::point(window.gamma()) * (I::point(1.) - ratio);
    let error = I::point(first)
        + I::point(window.acceleration()) * phase
        + I::point(window.speed()) * ratio;
    let upper = I::point(position)
        + I::point(window.speed()) * phase
        + I::point(body)
        + I::point(window.numeric_reserve());
    check(
        ratio.valid()
            && ratio.lo >= 0.
            && speed.valid()
            && speed.lo > 0.
            && error.valid()
            && upper.valid(),
        ReadError::Unsupported,
    )?;
    // Both actual jets are in the same finite affine chart. A difference
    // smaller than the actual target speed proves strictly positive IP2.
    if error.hi >= speed.lo || upper.hi > tol {
        return Ok(None);
    }
    work.charge(64)?;
    Ok(Some(super::super::Proof::Certified(upper.hi)))
}
pub(super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    d: &super::super::uv::Definition,
    support: &Surface,
    edge: &Edge,
    tol: f64,
    work: &mut Budget<'_>,
) -> Result<Option<super::super::Proof>> {
    match prove(source, curve_id, d, support, edge, tol, work) {
        Err(ReadError::Unsupported) => Ok(None),
        result => result,
    }
}
#[cfg(test)]
mod tests;
