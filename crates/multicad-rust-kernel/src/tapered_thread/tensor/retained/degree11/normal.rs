//! Exact degree-(21,21) polynomial P_u cross P_v of a degree-eleven cell.
//! Calculation only: a polynomial hull does not certify occupied material.
use super::*;
<<<<<<< HEAD
mod squared;
pub(super) use squared::{SQUARED11_WORK, SquaredNormal11};
=======
mod projected;
mod reduced;
pub use reduced::{ReducedNormal11, ReducedNormalCell11};
>>>>>>> d1743a74 (feat: enclose actual normals with certified elevation residuals)
/// Covers132*132 coefficient pairs, all three directed cross components,
/// product weights, differentiation, allocation/copies and setup.
pub const NORMAL11_WORK: usize = 400_000;
/// Root polynomial construction workspace and retained capacities. Callers
/// preflight this before construction; allocator bookkeeping/RSS are separate.
pub const NORMAL11_BYTES: usize = 64 * 1024;
pub struct NormalPolynomial11 {
    polynomial: NormalPolynomial,
    derivative: [[I; 3]; 2],
    derivative_hulls: [[[I; 3]; 9]; 2],
    hull: [I; 3],
}
fn choose(n: usize, k: usize) -> u64 {
    (0..k.min(n - k)).fold(1u64, |v, i| v * (n - i) as u64 / (i + 1) as u64)
}
impl NormalPolynomial11 {
    /// Coefficients use this supplied cell's unit-square coordinates. A caller
    /// owns aggregate work/storage charging; no source degree is inferred.
    pub fn new(cell: &TensorCell11) -> Result<Self, Error> {
        let du = derivative(&cell.0.controls, 0);
        let dv = derivative(&cell.0.controls, 1);
        let weights = std::array::from_fn::<_, 11, _>(|i| {
            std::array::from_fn::<_, 12, _>(|j| {
                I::point((choose(10, i) * choose(11, j)) as f64)
                    / I::point(choose(21, i + j) as f64)
            })
        });
        let mut controls = vec![vec![[I::point(0.); 3]; 22]; 22];
        for (i, row) in du.iter().enumerate() {
            for (j, u) in row.iter().enumerate() {
                for (a, row) in dv.iter().enumerate() {
                    for (b, v) in row.iter().enumerate() {
                        let weight = weights[i][a] * weights[b][j];
                        for k in 0..3 {
                            let cross =
                                u[(k + 1) % 3] * v[(k + 2) % 3] - u[(k + 2) % 3] * v[(k + 1) % 3];
                            controls[i + a][j + b][k] = controls[i + a][j + b][k] + weight * cross;
                        }
                    }
                }
            }
        }
        if controls.iter().flatten().flatten().any(|x| !x.valid()) {
            return Err(Error::Unresolved);
        }
        // Source tangent grids have finished their last use. Release them
        // before cloning the normal grid for its derivative-cache pass.
        drop(du);
        drop(dv);
        let hull = coefficient_hull(&controls);
        let derivative_hulls = std::array::from_fn(|axis| {
            let mut grid = controls.clone();
            std::array::from_fn(|order| {
                let value = coefficient_hull(&grid);
                if order < 8 {
                    let degree = if axis == 0 {
                        grid.len() - 1
                    } else {
                        grid[0].len() - 1
                    };
                    // Ascending updates read the next unmodified coefficient.
                    // No second derivative grid overlaps this single copy.
                    for i in 0..grid.len() - usize::from(axis == 0) {
                        for j in 0..grid[0].len() - usize::from(axis == 1) {
                            let next = if axis == 0 {
                                grid[i + 1][j]
                            } else {
                                grid[i][j + 1]
                            };
                            let current = grid[i][j];
                            grid[i][j] = std::array::from_fn(|k| {
                                (next[k] - current[k]) * I::point(degree as f64)
                            });
                        }
                    }
                    if axis == 0 {
                        grid.pop();
                    } else {
                        for row in &mut grid {
                            row.pop();
                        }
                    }
                }
                value
            })
        });
        let derivative = derivative_hulls.map(|rows| rows[1]);
        let result = Self {
            polynomial: NormalPolynomial { controls },
            derivative,
            derivative_hulls,
            hull,
        };
        if result.retained_budget_bytes() > NORMAL11_BYTES {
            return Err(Error::Limit);
        }
        Ok(result)
    }
    pub(crate) fn coefficients(&self) -> &[Vec<[I; 3]>] {
        &self.polynomial.controls
    }
    /// Fixed capacities requested by this constructor, checked again against
    /// actual capacities before a caller publishes a persistent cache.
    pub(crate) const fn retained_capacity_bound_bytes() -> usize {
        std::mem::size_of::<Self>()
            + 22 * std::mem::size_of::<Vec<BoxPoint>>()
            + 484 * std::mem::size_of::<BoxPoint>()
    }
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.polynomial.controls.capacity() * std::mem::size_of::<Vec<BoxPoint>>()
            + self
                .polynomial
                .controls
                .iter()
                .map(|r| r.capacity() * std::mem::size_of::<BoxPoint>())
                .sum::<usize>()
    }
    fn valid(uv: [[f64; 2]; 2]) -> bool {
        uv.iter().all(|d| {
            d[0].is_finite() && d[1].is_finite() && 0. <= d[0] && d[0] <= d[1] && d[1] <= 1.
        })
    }
    /// Original root-unit normal, without a child Jacobian factor.
    pub fn bounds(&self, uv: [[f64; 2]; 2]) -> Result<[[f64; 2]; 3], Error> {
        if !Self::valid(uv) {
            return Err(Error::Invalid);
        }
        let controls = restricted_grid(self.polynomial.controls.clone(), uv)?;
        let result = coefficient_hull(&controls);
        if result.iter().any(|x| !x.valid()) {
            return Err(Error::Unresolved);
        }
        Ok(result.map(|x| [x.lo, x.hi]))
    }
    /// Root-unit normal at whole parameter intervals. Factored Bernstein
    /// evaluation retains the anchor once, including source correlations.
    fn at_interval(&self, uv: [[f64; 2]; 2]) -> Result<[I; 3], Error> {
        if !Self::valid(uv) {
            return Err(Error::Invalid);
        }
        let bu = bernstein(21, I::new(uv[0][0], uv[0][1]));
        let bv = bernstein(21, I::new(uv[1][0], uv[1][1]));
        let anchor = self.polynomial.controls[0][0];
        let mut value = [I::point(0.); 3];
        for (i, row) in self.polynomial.controls.iter().enumerate() {
            for (j, p) in row.iter().enumerate() {
                if i == 0 && j == 0 {
                    continue;
                }
                let weight = bu[i] * bv[j];
                for k in 0..3 {
                    value[k] = value[k] + weight * (p[k] - anchor[k]);
                }
            }
        }
        let result = std::array::from_fn(|k| anchor[k] + value[k]);
        if result.iter().any(|x| !x.valid()) {
            return Err(Error::Unresolved);
        }
        Ok(result)
    }
    /// Normal and its first eight derivatives on THIS child unit square,
    /// including h_u*h_v and the selected axis chain powers. Root derivatives
    /// are differentiated before hulls to retain shared-control correlation.
    pub fn taylor8(
        &self,
        domain: [[f64; 2]; 2],
        along_v: bool,
    ) -> Result<[[[f64; 2]; 3]; 9], Error> {
        if !Self::valid(domain) {
            return Err(Error::Invalid);
        }
        self.polynomial
            .taylor8(domain, along_v)
            .map(|rows| rows.map(|p| p.map(|v| [v.lo, v.hi])))
    }
    pub fn at(&self, uv: [[f64; 2]; 2]) -> Result<[[f64; 2]; 3], Error> {
        self.at_interval(uv).map(|p| p.map(|v| [v.lo, v.hi]))
    }
    /// Exact polynomial centre plus whole-root derivative hulls. Every point
    /// of the requested box is covered by the mean-value theorem. This is a
    /// correlated normal-field enclosure, not a cross of independent tangents.
    pub fn lipschitz_bounds(&self, domain: [[f64; 2]; 2]) -> Result<[[f64; 2]; 3], Error> {
        if !Self::valid(domain) {
            return Err(Error::Invalid);
        }
        let mid = domain.map(|d| d[0] + (d[1] - d[0]) / 2.);
        let center = self.at_interval(mid.map(|m| [m, m]))?;
        let delta = std::array::from_fn::<_, 2, _>(|i| {
            I::new(domain[i][0], domain[i][1]) - I::point(mid[i])
        });
        let result = std::array::from_fn::<_, 3, _>(|k| {
            let b = center[k] + self.derivative[0][k] * delta[0] + self.derivative[1][k] * delta[1];
            I::new(b.lo.max(self.hull[k].lo), b.hi.min(self.hull[k].hi))
        });
        if result.iter().any(|x| !x.valid()) {
            return Err(Error::Unresolved);
        }
        Ok(result.map(|v| [v.lo, v.hi]))
    }
    /// Cached whole-root high derivative hulls with the cheaper proved
    /// Lipschitz enclosure for order zero. Higher coefficients use the same
    /// root derivatives and child Jacobian convention as taylor8.
    pub fn fast_taylor8(
        &self,
        domain: [[f64; 2]; 2],
        along_v: bool,
    ) -> Result<[[[f64; 2]; 3]; 9], Error> {
        let position = self.lipschitz_bounds(domain)?.map(|b| I::new(b[0], b[1]));
        let width = domain.map(|d| I::point(d[1]) - I::point(d[0]));
        let mut scale = width[0] * width[1];
        let mut factorial = 1.;
        let axis = usize::from(along_v);
        Ok(std::array::from_fn(|order| {
            let hull = if order == 0 {
                position
            } else {
                self.derivative_hulls[axis][order]
            };
            let value = hull.map(|v| v * scale / I::point(factorial));
            scale = scale * width[axis];
            factorial *= (order + 1) as f64;
            value.map(|v| [v.lo, v.hi])
        }))
    }
}

#[cfg(test)]
mod tests;
