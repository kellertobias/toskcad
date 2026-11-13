//! Any finite numeric preconditioner is qualified by the later directed test.
use super::*;
pub(super) type Matrix = [[I; 2]; 2];
pub(super) fn preconditioner(j: Matrix) -> BranchResult<[[f64; 2]; 2]> {
    let a = j.map(|row| row.map(|x| x.lo.midpoint(x.hi)));
    let scale = a.into_iter().flatten().map(f64::abs).fold(0., f64::max);
    if !scale.is_finite() || scale == 0. {
        return Err(PhillipsReadError::Unsupported);
    }
    let b = a.map(|row| row.map(|v| v / scale));
    let determinant = b[0][0] * b[1][1] - b[0][1] * b[1][0];
    if !determinant.is_finite() || determinant == 0. {
        return Err(PhillipsReadError::Unsupported);
    }
    let inverse = [[b[1][1], -b[0][1]], [-b[1][0], b[0][0]]]
        .map(|row| row.map(|v| (v / determinant) / scale));
    if inverse.into_iter().flatten().any(|v| !v.is_finite()) {
        return Err(PhillipsReadError::Unsupported);
    }
    Ok(inverse)
}
pub(super) fn contraction(m: [[f64; 2]; 2], j: Matrix) -> Matrix {
    std::array::from_fn(|r| {
        std::array::from_fn(|c| {
            I::point(if r == c { 1. } else { 0. })
                - (0..2).fold(I::point(0.), |s, k| s + I::point(m[r][k]) * j[k][c])
        })
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn directed_contraction_qualifies_scaled_nonsymmetric_preconditioners() {
        for scale in [1e-100, 1., 1e100] {
            let j = [[2., 0.25], [-0.5, 3.]].map(|row| row.map(|x| I::point(x * scale)));
            let m = preconditioner(j).unwrap();
            let q = contraction(m, j);
            assert!(
                q.into_iter()
                    .flatten()
                    .all(|x| x.valid() && x.magnitude() < 1e-12)
            );
        }
    }
    #[test]
    fn zero_singular_overflowing_and_unrepresentable_inverse_cannot_precondition_a_proof() {
        for j in [
            [[0.; 2]; 2],
            [[1., 2.], [2., 4.]],
            [[f64::INFINITY, 0.], [0., 1.]],
            [[f64::from_bits(1), 0.], [0., f64::from_bits(1)]],
        ] {
            assert!(preconditioner(j.map(|r| r.map(I::point))).is_err());
        }
    }
}
