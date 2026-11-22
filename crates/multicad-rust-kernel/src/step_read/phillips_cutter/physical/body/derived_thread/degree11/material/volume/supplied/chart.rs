//! Actual affine expressions and complete-curve normalized-ray bounds.
use super::*;
type Matrix = [[I; 3]; 3]; // Columns, never rows.
pub(super) struct ChartBounds {
    pub translation: I,
    pub linear: I,
    pub author_norm: I,
    pub inverse_author_norm: I,
    pub author_determinant: I,
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
pub(super) fn norm(values: impl IntoIterator<Item = I>) -> I {
    let square = values.into_iter().fold(I::point(0.), |s, x| s + x.square());
    I::new(0., square.hi).sqrt()
}
// The largest eigenvalue of the actual symmetric Gram matrix is at most
// its maximum absolute row sum. This remains valid for arbitrary affine
// columns and is tight near a proper rotation, unlike a Frobenius sqrt(3).
fn operator(matrix: Matrix) -> I {
    if matrix.iter().flatten().any(|x| !x.valid()) {
        return I::new(f64::NEG_INFINITY, f64::INFINITY);
    }
    let maximum = (0..3)
        .map(|i| {
            (0..3).fold(I::point(0.), |s, j| {
                s + I::point(dot(matrix[i], matrix[j]).magnitude())
            })
        })
        .map(|x| x.hi)
        .fold(0., f64::max);
    I::new(0., maximum).sqrt()
}
pub(super) fn compare(frame: [Point; 3], chart: CommonMouthChart) -> PResult<ChartBounds> {
    let [origin, u, v] = frame;
    // EXACT retained binary64 cross-product columns, matching the original
    // writer. No mathematical reorthogonalization or rigid-frame assumption.
    let n =
        std::array::from_fn(|i| u[(i + 1) % 3] * v[(i + 2) % 3] - u[(i + 2) % 3] * v[(i + 1) % 3]);
    let a = [u, v, n].map(|p| p.map(I::point));
    let b = [chart.columns[1], chart.columns[2], chart.columns[3]].map(|p| p.map(I::point));
    let cofactors = [cross(a[1], a[2]), cross(a[2], a[0]), cross(a[0], a[1])];
    let determinant = dot(a[0], cofactors[0]);
    compatible(determinant.valid() && determinant.lo > 0. && determinant.hi.is_finite())?;
    let inverse = std::array::from_fn(|j| std::array::from_fn(|i| cofactors[i][j] / determinant));
    let delta = std::array::from_fn(|j| std::array::from_fn(|i| a[j][i] - b[j][i]));
    let translation = norm((0..3).map(|i| {
        I::point(origin[i])
            - (I::point(chart.columns[0][i])
                - I::point(chart.columns[3][i]) * I::point(chart.outer_top_z))
    }));
    let result = ChartBounds {
        translation,
        linear: operator(delta),
        author_norm: operator(a),
        inverse_author_norm: operator(inverse),
        author_determinant: determinant,
    };
    compatible(
        [
            result.translation,
            result.linear,
            result.author_norm,
            result.inverse_author_norm,
        ]
        .into_iter()
        .all(|x| x.valid() && x.hi.is_finite()),
    )?;
    Ok(result)
}
// BOTH curves must have norm>=rho over the complete domain. The joining
// segment need not stay outside that ball. Dnormalize=(Id-qq^T)/r gives
// ||Dnormalize(a)-Dnormalize(b)|| <=5||a-b||/rho² from the two endpoint norms.
// Derivatives are with respect to the SAME original normalized curve cell.
pub(super) fn ray(delta: I, derivative_delta: I, target_speed: I, rho: I) -> PResult<[I; 2]> {
    compatible(rho.valid() && rho.lo > 0.)?;
    let position = I::point(2.) * delta / rho;
    let derivative = derivative_delta / rho + I::point(5.) * delta * target_speed / rho.square();
    compatible(
        [position, derivative]
            .into_iter()
            .all(|x| x.valid() && x.hi.is_finite()),
    )?;
    Ok([position, derivative])
}
