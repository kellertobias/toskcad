//! Transport the ORIGINAL positive-dot TANGENT cone, not a unit-normal chord.
use super::*;
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn length(v: [I; 3]) -> I {
    v.into_iter()
        .fold(I::point(0.), |s, x| s + x.square())
        .sqrt()
}
fn facet(p: [[f64; 3]; 3]) -> [I; 3] {
    let [a, b, c] = p.map(|p| p.map(I::point));
    cross(
        std::array::from_fn(|k| b[k] - a[k]),
        std::array::from_fn(|k| c[k] - a[k]),
    )
}
fn cone(source: [I; 3], actual: [I; 3]) -> Result<f64> {
    let d = dot(source, actual);
    if !d.valid() || d.lo <= 0. {
        return Err(Error::Unsupported);
    }
    let r = I::point(length(cross(source, actual)).hi) / I::point(d.lo);
    if !r.valid() {
        return Err(Error::Unsupported);
    }
    Ok(r.hi)
}
fn spectrum(scale: I) -> Result<(f64, f64)> {
    let square = scale.square();
    if !scale.valid() || scale.lo <= 0. || !square.valid() || square.lo <= 0. {
        return Err(Error::Unsupported);
    }
    Ok((square.lo.min(1.), square.hi.max(1.)))
}
pub(super) fn child_angle(scale: I, angle: f64, w: &mut Work<'_>) -> Result<f64> {
    w.step(64)?;
    let (lo, hi) = spectrum(scale)?;
    let wanted = I::point(angle) - I::point(angle) / I::point(4096.);
    let half_gap = (I::point(hi) - I::point(lo)) / I::point(2.);
    let r = (I::point(lo) * wanted / (I::point(hi) + half_gap * wanted)).lo;
    if !r.is_finite() || r <= 0. {
        return Err(Error::Unsupported);
    }
    // The existing literal child accepts angles>=.01. At its lower frontier
    // keep that valid request and let the actual final proof refuse if needed.
    Ok(r.max(0.01).min(angle))
}
pub(super) fn transported(scale: I, r: f64, w: &mut Work<'_>) -> Result<f64> {
    w.step(64)?;
    if !r.is_finite() || r < 0. {
        return Err(Error::Invalid);
    }
    let (lo, hi) = spectrum(scale)?;
    // Existing literal cells earn N.F>0 and |N×F|<=r(N.F).
    // For unit f=F/|F|, N=a(f+v), v⊥f, |v|<=r. With K=diag(s,s,1)
    // and M=KᵀK, |vᵀMf|<=((hi-lo)/2)r and |Kf×Kv|<=hi*r.
    // Thus this SAME transformed cone remains positive-dot and tan(theta)
    // obeys the quotient below. No libm angle/chord conversion is assumed.
    transported_bounds(lo, hi, r)
}
fn transported_bounds(lo: f64, hi: f64, r: f64) -> Result<f64> {
    let denominator = I::point(lo) - (I::point(hi) - I::point(lo)) * I::point(r) / I::point(2.);
    if !denominator.valid() || denominator.lo <= 0. {
        return Err(Error::Unsupported);
    }
    let result = I::point(hi) * I::point(r) / denominator;
    if !result.valid() {
        return Err(Error::Unsupported);
    }
    Ok(result.hi)
}
pub(super) fn transported_metric(delta: f64, r: f64, w: &mut Work<'_>) -> Result<f64> {
    w.step(64)?;
    if !delta.is_finite() || !(0.0..1.0).contains(&delta) || !r.is_finite() || r < 0. {
        return Err(Error::Unsupported);
    }
    // delta bounds ||AᵀA-I||. Cofactor-Gram eigenvalues therefore lie
    // in[(1-delta)^2,(1+delta)^2]. Proper uniform file B only rotates/scales
    // this spectrum and does not change the quotient's condition.
    let lo = (I::point(1.) - I::point(delta)).square().lo;
    let hi = (I::point(1.) + I::point(delta)).square().hi;
    transported_bounds(lo, hi, r)
}
pub(super) fn world_rounding(
    raw: [[f64; 3]; 3],
    actual: [[f64; 3]; 3],
    axes: [[[f64; 2]; 3]; 3],
    placement: Option<&super::super::super::LiteralOccurrencePlacement>,
    w: &mut Work<'_>,
) -> Result<(f64, [I; 3])> {
    w.step(1024)?;
    let old = facet(raw);
    // B=sR, s>0, detR=+1. Its positive cofactor scale s² cancels from
    // normal direction; retain actual EXPRESS R and then actual Native A.
    let mut ideal = std::array::from_fn(|i| {
        (0..3).fold(I::point(0.), |v, j| {
            v + I::new(axes[j][i][0], axes[j][i][1]) * old[j]
        })
    });
    if let Some(a) = placement {
        ideal = a
            .normal_bounds(ideal.map(|x| [x.lo, x.hi]))?
            .map(|x| I::new(x[0], x[1]));
    }
    let actual = facet(actual);
    Ok((cone(ideal, actual)?, actual))
}
pub(super) fn rounding(
    raw: [[f64; 3]; 3],
    actual: [[f64; 3]; 3],
    scale: I,
    w: &mut Work<'_>,
) -> Result<f64> {
    w.step(512)?;
    let mut ideal = facet(raw);
    ideal[0] = ideal[0] * scale;
    ideal[1] = ideal[1] * scale;
    cone(ideal, facet(actual))
}
pub(super) fn planar(source: [f64; 3], triangle: [[f64; 3]; 3], w: &mut Work<'_>) -> Result<f64> {
    w.step(384)?;
    cone(source.map(I::point), facet(triangle))
}
pub(super) fn combine(a: f64, b: f64, w: &mut Work<'_>) -> Result<f64> {
    w.step(32)?;
    let denominator = I::point(1.) - I::point(a) * I::point(b);
    if !denominator.valid() || denominator.lo <= 0. {
        return Err(Error::Unsupported);
    }
    // atan(a)+atan(b)<pi/2 and tan addition retains a positive-dot cone.
    let r = (I::point(a) + I::point(b)) / denominator;
    if !r.valid() {
        return Err(Error::Unsupported);
    }
    Ok(r.hi)
}
