//! Compare actual defining affine expressions, including the symbolic top-Z
//! subtraction. No independently rounded translated origin defines either map.
use super::*;
#[derive(Clone, Copy)]
pub(super) struct WorldDifference {
    pub translation: I,
    pub linear: I,
    pub target_norm: I,
    pub occurrence_determinant: I,
    /// Determinant of the target's own file operator (one at identity).
    pub target_determinant: I,
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn upper_norm(v: impl IntoIterator<Item = I>) -> I {
    let square = v.into_iter().fold(I::point(0.), |s, x| s + x.square());
    I::new(0., square.hi).sqrt()
}
fn columns(chart: CommonMouthChart, occurrence: [[I; 3]; 4]) -> [[I; 3]; 4] {
    let c = chart.columns.map(|p| p.map(I::point));
    let linear = |p: [I; 3]| {
        std::array::from_fn(|i| (0..3).fold(I::point(0.), |s, k| s + occurrence[k + 1][i] * p[k]))
    };
    let offset = linear(std::array::from_fn(|i| {
        c[0][i] - c[3][i] * I::point(chart.outer_top_z)
    }));
    [
        std::array::from_fn(|i| occurrence[0][i] + offset[i]),
        linear(c[1]),
        linear(c[2]),
        linear(c[3]),
    ]
}
pub(super) fn compare(
    source: CommonMouthChart,
    target: CommonMouthChart,
    occurrence: [[f64; 3]; 4],
) -> std::result::Result<WorldDifference, PhillipsReadError> {
    let identity =
        [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].map(|p| p.map(I::point));
    compare_through(source, target, occurrence, identity)
}
fn proper(frame: [[I; 3]; 4]) -> std::result::Result<I, PhillipsReadError> {
    for i in 0..3 {
        for j in 0..3 {
            let error = dot(frame[i + 1], frame[j + 1]) - I::point(if i == j { 1. } else { 0. });
            if !error.valid() || error.magnitude() > 1e-6 {
                return Err(PhillipsReadError::Invalid);
            }
        }
    }
    let determinant = dot(frame[1], cross(frame[2], frame[3]));
    if !determinant.valid() || determinant.lo <= 0. {
        return Err(PhillipsReadError::Invalid);
    }
    Ok(determinant)
}
/// As [compare], with the target's own decoded file operator (an enclosure
/// of a mapped item's exact map) applied to its chart instead of identity.
pub(super) fn compare_through(
    source: CommonMouthChart,
    target: CommonMouthChart,
    occurrence: [[f64; 3]; 4],
    target_operator: [[I; 3]; 4],
) -> std::result::Result<WorldDifference, PhillipsReadError> {
    if occurrence
        .into_iter()
        .flatten()
        .any(|v| !v.is_finite() || v.abs() > 1e7)
        || target_operator
            .into_iter()
            .flatten()
            .any(|v| !v.valid() || v.magnitude() > 1e7)
    {
        return Err(PhillipsReadError::Invalid);
    }
    let occurrence = occurrence.map(|p| p.map(I::point));
    let determinant = proper(occurrence)?;
    let target_determinant = proper(target_operator)?;
    let a = columns(source, occurrence);
    let b = columns(target, target_operator);
    let delta: [[I; 3]; 4] = std::array::from_fn(|j| std::array::from_fn(|i| a[j][i] - b[j][i]));
    let out = WorldDifference {
        translation: upper_norm(delta[0]),
        linear: upper_norm(delta[1..].iter().flatten().copied()),
        target_norm: upper_norm(b[1..].iter().flatten().copied()),
        occurrence_determinant: determinant,
        target_determinant,
    };
    if ![out.translation, out.linear, out.target_norm]
        .into_iter()
        .all(|v| v.valid() && v.hi.is_finite())
    {
        return Err(PhillipsReadError::Unsupported);
    }
    Ok(out)
}
