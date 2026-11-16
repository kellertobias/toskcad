//! Whole-disk plane support comparison from exact normalized numeric families.
use super::*;
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::point(0.), |s, k| s + a[k] * b[k])
}
fn abs(x: I) -> f64 {
    x.lo.abs().max(x.hi.abs())
}
pub(super) fn disk_transfer(
    origin: [I; 3],
    normal: [I; 3],
    reference: [I; 3],
    chart: CommonMouthChart,
    [radius, z]: [f64; 2],
    work: &mut Work<'_>,
) -> FResult<f64> {
    work.step(2048)?;
    if !radius.is_finite()
        || radius <= 0.
        || !z.is_finite()
        || chart.columns.into_iter().flatten().any(|x| !x.is_finite())
        || !chart.outer_top_z.is_finite()
    {
        return Err(PhillipsReadError::Unsupported);
    }
    let raw = support::Frame::new(origin, normal, reference).map_err(source_error)?;
    let [o, u, v, n] = chart.columns.map(|p| p.map(I::point));
    if dot(raw.axes[2], n).lo <= 0. {
        return Err(PhillipsReadError::Unsupported);
    }
    let delta = std::array::from_fn(|k| {
        o[k] + n[k] * (I::point(z) - I::point(chart.outer_top_z)) - raw.origin[k]
    });
    // |cos|,|sin|<=1 over the entire filled disk, including its collapsed
    // shared-crest boundary. No sampled face positions enter this bound.
    let result = (I::point(abs(dot(delta, raw.axes[2])))
        + I::point(radius)
            * (I::point(abs(dot(u, raw.axes[2]))) + I::point(abs(dot(v, raw.axes[2])))))
    .hi;
    if !result.is_finite() {
        return Err(PhillipsReadError::Unsupported);
    }
    Ok(result)
}
