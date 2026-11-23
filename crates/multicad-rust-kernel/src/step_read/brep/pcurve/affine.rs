//! Actual affine lifts and positive-basis finite NURBS correspondence.
use super::*;
use num_bigint::BigInt;
type V = [I; 3];
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn scale(a: V, b: I) -> V {
    a.map(|x| x * b)
}
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn valid(p: V) -> bool {
    p.into_iter().all(I::valid)
}

/// Actual value in the original natural domain; no clamping/extrapolation.
pub(super) fn at(curve: &Curve, t: f64, budget: &mut Budget<'_>) -> Result<V> {
    budget.charge(1024)?;
    if !t.is_finite() {
        return Err(ReadError::Unsupported);
    }
    let t = I::point(t);
    let p = match curve {
        Curve::Line { origin, direction } => add(vector(*origin), scale(vector(*direction), t)),
        Curve::Circle { frame, radius } => add(
            vector(frame.origin),
            add(
                scale(vector(frame.x), t.cos() * I::point(*radius)),
                scale(vector(frame.y), t.sin() * I::point(*radius)),
            ),
        ),
        Curve::Ellipse {
            frame,
            semiaxes: [a, b],
        } => add(
            vector(frame.origin),
            add(
                scale(vector(frame.x), t.cos() * I::point(*a)),
                scale(vector(frame.y), t.sin() * I::point(*b)),
            ),
        ),
        Curve::Parabola { frame, focal } => add(
            vector(frame.origin),
            add(
                scale(vector(frame.x), I::point(*focal) * t.square()),
                scale(vector(frame.y), I::point(2.) * I::point(*focal) * t),
            ),
        ),
        Curve::Nurbs {
            curve,
            closed: false,
        } => return nurbs_at(curve, t.lo, budget),
        _ => return Err(ReadError::Unsupported),
    };
    if valid(p) {
        Ok(p)
    } else {
        Err(ReadError::Unsupported)
    }
}
fn nurbs_at(curve: &NurbsCurve, t: f64, budget: &mut Budget<'_>) -> Result<V> {
    let [a, b] = curve.range();
    if t < a || t > b {
        return Err(ReadError::Unsupported);
    }
    let p = curve.degree();
    let (distinct, multiplicities) = curve.knots();
    if (t == a && multiplicities[0] == p + 1)
        || (t == b && multiplicities[multiplicities.len() - 1] == p + 1)
    {
        return Ok(vector(if t == a {
            curve.points()[0]
        } else {
            *curve.points().last().ok_or(ReadError::Invalid)?
        }));
    }
    budget.charge(
        curve
            .points()
            .len()
            .saturating_add(p + 1)
            .saturating_mul(32)
            .saturating_add(p.saturating_mul(p + 1).saturating_mul(256)),
    )?;
    let mut knots = Vec::new();
    knots
        .try_reserve_exact(curve.points().len() + p + 1)
        .map_err(|_| ReadError::Limit)?;
    for (&k, &m) in distinct.iter().zip(multiplicities) {
        knots.extend(std::iter::repeat_n(k, m));
    }
    let n = curve.points().len();
    let span = if t == b {
        n - 1
    } else {
        knots
            .partition_point(|&k| k <= t)
            .saturating_sub(1)
            .clamp(p, n - 1)
    };
    let maximum = curve
        .weights()
        .map_or(1., |w| w.iter().copied().fold(0_f64, f64::max));
    let mut d = [[I::point(0.); 4]; 33];
    for (j, slot) in d.iter_mut().enumerate().take(p + 1) {
        let index = span - p + j;
        let w = I::point(curve.weights().map_or(1., |w| w[index])) / I::point(maximum);
        *slot = [
            I::point(curve.points()[index][0]) * w,
            I::point(curve.points()[index][1]) * w,
            I::point(curve.points()[index][2]) * w,
            w,
        ];
    }
    for r in 1..=p {
        for j in (r..=p).rev() {
            let index = span - p + j;
            let den = I::point(knots[index + p - r + 1]) - I::point(knots[index]);
            if !den.valid() || den.lo <= 0. {
                return Err(ReadError::Unsupported);
            }
            let alpha = (I::point(t) - I::point(knots[index])) / den;
            let previous = d[j - 1];
            for (value, before) in d[j].iter_mut().zip(previous) {
                *value = (I::point(1.) - alpha) * before + alpha * *value;
            }
        }
    }
    if !d[p][3].valid() || d[p][3].lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let value = std::array::from_fn(|i| d[p][i] / d[p][3]);
    if valid(value) {
        Ok(value)
    } else {
        Err(ReadError::Unsupported)
    }
}

/// A positive affine map from the actual 3D LINE parameter into the unbounded
/// lifted UV line. A negative/unknown map does not earn IP2 same-sense.
pub(super) fn line(
    origin: V,
    direction: V,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    budget.charge(2048)?;
    let Curve::Line {
        origin: actual,
        direction: d,
    } = edge.curve
    else {
        return Ok(Proof::Unsupported);
    };
    let den = dot(direction, direction);
    if !den.valid() || den.lo <= 0. {
        return Ok(Proof::Unsupported);
    }
    let slope = dot(vector(d), direction) / den;
    if !slope.valid() || slope.lo <= 0. {
        return Ok(Proof::Unsupported);
    }
    let delta = difference(vector(actual), origin);
    let offset = dot(delta, direction) / den;
    let residual_origin = difference(delta, scale(direction, offset));
    let residual_direction = difference(vector(d), scale(direction, slope));
    let mut upper = 0_f64;
    for t in edge.range {
        let n = norm(add(residual_origin, scale(residual_direction, I::point(t))));
        if !n.valid() {
            return Ok(Proof::Unsupported);
        }
        upper = upper.max(n.hi);
    }
    Ok(if upper <= tol {
        Proof::Certified(upper)
    } else {
        Proof::Unsupported
    })
}
// Scale exact binary64 dyadics to a common power-of-two denominator. At most
// six values, <=2098 bits each; products <=4196 bits. No approximate ratio or
// floating equality of cross-products can authorize a knot map.
fn dyadics(values: &[f64]) -> Result<Vec<BigInt>> {
    if values.len() > 6 || values.iter().any(|x| !x.is_finite()) {
        return Err(ReadError::Invalid);
    }
    let parts: Vec<_> = values
        .iter()
        .map(|&x| {
            let bits = x.to_bits();
            let e = ((bits >> 52) & 2047) as i32;
            let m = (bits & ((1_u64 << 52) - 1)) | if e == 0 { 0 } else { 1_u64 << 52 };
            (
                if bits >> 63 == 0 {
                    m as i64
                } else {
                    -(m as i64)
                },
                if e == 0 { -1074 } else { e - 1023 - 52 },
            )
        })
        .collect();
    let minimum = parts.iter().map(|x| x.1).min().unwrap_or(0);
    Ok(parts
        .into_iter()
        .map(|(m, e)| BigInt::from(m) << ((e - minimum) as usize))
        .collect())
}
fn knots_match(
    a: (&[f64], &[usize]),
    b: (&[f64], &[usize]),
    forward: bool,
    budget: &mut Budget<'_>,
) -> Result<bool> {
    if a.0.len() != b.0.len() || a.1.len() != b.1.len() {
        return Ok(false);
    }
    budget.charge(a.0.len().saturating_mul(16))?;
    if forward && a == b {
        return Ok(true);
    }
    if a.0.len() < 2 {
        return Ok(false);
    }
    for i in 0..a.0.len() {
        budget.charge(32_768)?;
        let j = if forward { i } else { b.0.len() - 1 - i };
        if a.1[i] != b.1[j] {
            return Ok(false);
        }
        let v = dyadics(&[
            a.0[i],
            a.0[0],
            *a.0.last().unwrap(),
            b.0[j],
            b.0[0],
            *b.0.last().unwrap(),
        ])?;
        let left = (&v[0] - &v[1]) * (&v[5] - &v[4]);
        let right = if forward {
            (&v[3] - &v[4]) * (&v[2] - &v[1])
        } else {
            (&v[5] - &v[3]) * (&v[2] - &v[1])
        };
        if left != right {
            return Ok(false);
        }
    }
    Ok(true)
}
fn weights_match(
    a: Option<&[f64]>,
    b: Option<&[f64]>,
    count: usize,
    forward: bool,
    budget: &mut Budget<'_>,
) -> Result<bool> {
    if a.is_some_and(|w| w.len() != count) || b.is_some_and(|w| w.len() != count) {
        return Ok(false);
    }
    let av = |i| a.map_or(1., |w| w[i]);
    let bv = |i| b.map_or(1., |w| w[if forward { i } else { count - 1 - i }]);
    budget.charge(count.saturating_mul(16))?;
    if (0..count).all(|i| av(i) == bv(i)) {
        return Ok(true);
    }
    for i in 0..count {
        budget.charge(32_768)?;
        let v = dyadics(&[av(i), bv(0), bv(i), av(0)])?;
        if &v[0] * &v[1] != &v[2] * &v[3] {
            return Ok(false);
        }
    }
    Ok(true)
}
fn separating_witness(controls: &[V], edge: &Edge, budget: &mut Budget<'_>) -> Result<f64> {
    let mut lower = 0_f64;
    for t in edge.range {
        let point = match at(&edge.curve, t, budget) {
            Ok(p) => p,
            Err(ReadError::Unsupported) => continue,
            Err(e) => return Err(e),
        };
        let mut axes = [[I::point(0.); 3]; 8];
        for (slot, control) in axes
            .iter_mut()
            .zip([controls.first(), controls.last()].into_iter().flatten())
        {
            // Direction is a proposal only. ALL original positive-basis
            // control coefficients must prove the same separating functional.
            *slot = vector(std::array::from_fn(|i| {
                control[i].midpoint() - point[i].midpoint()
            }));
        }
        for axis in 0..3 {
            axes[2 + 2 * axis][axis] = I::point(1.);
            axes[3 + 2 * axis][axis] = I::point(-1.);
        }
        for axis in axes {
            let length = norm(axis);
            if !length.valid() || length.lo <= 0. {
                continue;
            }
            let mut projection = f64::INFINITY;
            for &p in controls {
                budget.charge(128)?;
                let d = dot(axis, difference(p, point));
                if !d.valid() {
                    projection = f64::NEG_INFINITY;
                    break;
                }
                projection = projection.min(d.lo);
            }
            if projection > 0. {
                let d = I::point(projection) / I::point(length.hi);
                if d.valid() {
                    lower = lower.max(d.lo);
                }
            }
        }
    }
    Ok(lower)
}
fn segment(
    controls: &[V],
    forward: bool,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    budget.charge(2048)?;
    let Curve::Line { origin, direction } = &edge.curve else {
        return Ok(Proof::Unsupported);
    };
    let den = dot(vector(*direction), vector(*direction));
    if !den.valid() || den.lo <= 0. {
        return Ok(Proof::Unsupported);
    }
    let points = if forward {
        [controls[0], controls[1]]
    } else {
        [controls[1], controls[0]]
    };
    let parameters = points.map(|p| dot(difference(p, vector(*origin)), vector(*direction)) / den);
    if parameters.into_iter().any(|p| !p.valid()) {
        return Ok(Proof::Unsupported);
    }
    let travel = parameters[1] - parameters[0];
    if travel.valid() && travel.lo > 0. {
        let mut perpendicular = 0_f64;
        for (p, t) in points.into_iter().zip(parameters) {
            let n = norm(difference(
                p,
                add(vector(*origin), scale(vector(*direction), t)),
            ));
            if !n.valid() {
                return Ok(Proof::Unsupported);
            }
            perpendicular = perpendicular.max(n.hi);
        }
        let length = den.sqrt();
        let gap = (I::point(parameters[0].hi) - I::point(edge.range[0]))
            .hi
            .max((I::point(edge.range[1]) - I::point(parameters[1].lo)).hi)
            .max(0.);
        let upper = (I::point(perpendicular) + I::point(gap) * length).hi;
        // Nonempty interior overlap permits an increasing finite affine
        // projected correspondence. Positive two-control rational weights
        // merely reparameterize this SAME segment by an increasing map.
        if edge.range[0].max(parameters[0].hi) < edge.range[1].min(parameters[1].lo)
            && upper.is_finite()
            && upper <= tol
        {
            return Ok(Proof::Certified(upper));
        }
    }
    let lower = separating_witness(controls, edge, budget)?;
    Ok(if lower > tol {
        Proof::Disagreement(lower)
    } else {
        Proof::Unsupported
    })
}
/// The source coefficients come from the ACTUAL affine lift. Identical positive
/// rational bases after an exact affine knot map prove the entire curve, not a
/// finite collection of evaluation points. `forward` is the UV-parameter sense
/// relative to the source's original increasing knot parameter.
#[allow(
    clippy::too_many_arguments,
    reason = "The original rational basis, orientation, finite target and proof budget are independently bound."
)]
pub(super) fn nurbs(
    degree: usize,
    knots: (&[f64], &[usize]),
    weights: Option<&[f64]>,
    controls: &[V],
    forward: bool,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    budget.charge(1024)?;
    if controls.is_empty() || controls.iter().any(|p| !valid(*p)) {
        return Ok(Proof::Unsupported);
    }
    budget.charge(controls.len().saturating_mul(16))?;
    if weights
        .is_some_and(|w| w.len() != controls.len() || w.iter().any(|x| !x.is_finite() || *x <= 0.))
    {
        return Err(ReadError::Invalid);
    }
    if degree == 1
        && controls.len() == 2
        && knots.0.len() == 2
        && knots.1 == [2, 2]
        && matches!(edge.curve, Curve::Line { .. })
    {
        return segment(controls, forward, edge, tol, budget);
    }
    let Curve::Nurbs {
        curve,
        closed: false,
    } = &edge.curve
    else {
        return Ok(Proof::Unsupported);
    };
    if controls.len() == curve.points().len()
        && degree == curve.degree()
        && knots_match(knots, curve.knots(), forward, budget)?
        && weights_match(weights, curve.weights(), controls.len(), forward, budget)?
    {
        let mut upper = 0_f64;
        for (i, &p) in controls.iter().enumerate() {
            budget.charge(128)?;
            let j = if forward { i } else { controls.len() - 1 - i };
            let n = norm(difference(p, vector(curve.points()[j])));
            if !n.valid() {
                return Ok(Proof::Unsupported);
            }
            upper = upper.max(n.hi);
        }
        if upper <= tol {
            return Ok(Proof::Certified(upper));
        }
    }
    // A failed parametric upper estimate is not Invalid. This separate point
    // versus ENTIRE positive-basis image proof is parameter-independent.
    let lower = separating_witness(controls, edge, budget)?;
    Ok(if lower > tol {
        Proof::Disagreement(lower)
    } else {
        Proof::Unsupported
    })
}

/// A finite degree-one positive rational spline is exactly its control
/// segment, with a strictly monotone (possibly non-affine) parameterization.
/// Its finite extent is retained; the caller must earn full natural coverage.
struct BoundaryIsoline {
    origin: Point,
    direction: Point,
    finite: Option<[Point; 2]>,
}
fn boundary_isoline(curve: &Curve, budget: &mut Budget<'_>) -> Result<Option<BoundaryIsoline>> {
    budget.charge(512)?;
    match curve {
        Curve::Line { origin, direction } => Ok(Some(BoundaryIsoline {
            origin: *origin,
            direction: *direction,
            finite: None,
        })),
        Curve::Nurbs {
            curve,
            closed: false,
        } if curve.degree() == 1 && curve.points().len() == 2 => {
            let (knots, m) = curve.knots();
            if knots.len() != 2 || m != [2, 2] {
                return Ok(None);
            }
            let a = curve.points()[0];
            let b = curve.points()[1];
            let direction = std::array::from_fn(|i| b[i] - a[i]);
            if direction.iter().any(|x| !x.is_finite()) {
                return Ok(None);
            }
            Ok(Some(BoundaryIsoline {
                origin: a,
                direction,
                finite: Some([a, b]),
            }))
        }
        _ => Ok(None),
    }
}
pub(super) fn certify(
    surface: &Surface,
    uv: &Curve,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Proof> {
    match (surface, uv) {
        (
            Surface::Plane(frame),
            Curve::Nurbs {
                curve,
                closed: false,
            },
        ) => {
            let mut controls = Vec::new();
            controls
                .try_reserve_exact(curve.points().len())
                .map_err(|_| ReadError::Limit)?;
            check(
                controls.capacity() * std::mem::size_of::<V>() <= MAX_WORKSPACE / 8,
                ReadError::Limit,
            )?;
            for &p in curve.points() {
                budget.charge(256)?;
                controls.push(add(
                    vector(frame.origin),
                    add(
                        scale(vector(frame.x), I::point(p[0])),
                        scale(vector(frame.y), I::point(p[1])),
                    ),
                ));
            }
            nurbs(
                curve.degree(),
                curve.knots(),
                curve.weights(),
                &controls,
                true,
                edge,
                tol,
                budget,
            )
        }
        (Surface::Nurbs { surface, closed }, uv) => {
            let Some(BoundaryIsoline {
                origin,
                direction,
                finite,
            }) = boundary_isoline(uv, budget)?
            else {
                return Ok(Proof::Unsupported);
            };
            let moving = if direction[0] != 0. && direction[1] == 0. {
                0
            } else if direction[0] == 0. && direction[1] != 0. {
                1
            } else {
                return Ok(Proof::Unsupported);
            };
            let fixed = 1 - moving;
            if let Some([a, b]) = finite {
                let range = surface.range(moving);
                let actual = if direction[moving] > 0. {
                    [a[moving], b[moving]]
                } else {
                    [b[moving], a[moving]]
                };
                if actual != range || a[fixed] != b[fixed] {
                    return Ok(Proof::Unsupported);
                }
            }
            // A certified closure still wraps parameters across the seam, and
            // this isoline proof has no seam convention for it: refuse, never
            // report a disagreement. No endpoint convention is substituted.
            if closed[moving] || closed[fixed] {
                return Ok(Proof::Unsupported);
            }
            let range = surface.range(fixed);
            let (_, fm) = surface.knots(fixed);
            let degrees = surface.degrees();
            if fm[0] != degrees[fixed] + 1 || fm[fm.len() - 1] != degrees[fixed] + 1 {
                return Ok(Proof::Unsupported);
            }
            let counts = surface.counts();
            let fixed_index = if origin[fixed] == range[0] {
                0
            } else if origin[fixed] == range[1] {
                counts[fixed] - 1
            } else {
                return Ok(Proof::Unsupported);
            };
            let mut controls = Vec::new();
            let mut weights = Vec::new();
            controls
                .try_reserve_exact(counts[moving])
                .map_err(|_| ReadError::Limit)?;
            if surface.weights().is_some() {
                weights
                    .try_reserve_exact(counts[moving])
                    .map_err(|_| ReadError::Limit)?;
            }
            check(
                controls.capacity() * std::mem::size_of::<V>() + weights.capacity() * 8
                    <= MAX_WORKSPACE / 8,
                ReadError::Limit,
            )?;
            for i in 0..counts[moving] {
                budget.charge(256)?;
                let index = if moving == 0 {
                    i * counts[1] + fixed_index
                } else {
                    fixed_index * counts[1] + i
                };
                controls.push(vector(surface.points()[index]));
                if let Some(w) = surface.weights() {
                    weights.push(w[index]);
                }
            }
            nurbs(
                degrees[moving],
                surface.knots(moving),
                (!weights.is_empty()).then_some(weights.as_slice()),
                &controls,
                direction[moving] > 0.,
                edge,
                tol,
                budget,
            )
        }
        _ => Ok(Proof::Unsupported),
    }
}
