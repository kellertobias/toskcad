//! Exact scalar N dot N, retaining component correlations before derivatives.
//! Numerical calculation only; no occupied-material admission is granted.
use super::*;
// 234256 pairs ×8 interval operations=1874048; cached differences/hulls
// <114k, binomial weights/integer setup <40k, whole-grid copies/checks and
// fixed metadata within the remaining 471k. Credits, not CPU instructions.
pub(crate) const SQUARED11_WORK: usize = 2_500_000;
const CONSTRUCTOR_BYTES: usize = 128 * 1024;
const GRID_BYTES: usize = 43 * 43 * std::mem::size_of::<I>() + 43 * std::mem::size_of::<Vec<I>>();
pub(crate) struct SquaredNormal11 {
    source: Arc<NormalPolynomial11>,
    controls: Vec<Vec<I>>,
    derivative: [[I; 9]; 2],
    hull: I,
}
fn hull(rows: &[Vec<I>]) -> I {
    rows.iter()
        .flatten()
        .fold(I::new(f64::INFINITY, f64::NEG_INFINITY), |b, x| {
            I::new(b.lo.min(x.lo), b.hi.max(x.hi))
        })
}
impl SquaredNormal11 {
    pub(crate) fn new(source: Arc<NormalPolynomial11>) -> Result<Self, Error> {
        // Preflight before allocating the scalar grid or derivative clone.
        let peak = source
            .retained_budget_bytes()
            .checked_add(2 * GRID_BYTES + 8192 + 4096 + std::mem::size_of::<Self>())
            .ok_or(Error::Limit)?;
        if peak > CONSTRUCTOR_BYTES {
            return Err(Error::Limit);
        }
        // All binomial integers and numerator products are <=2^53, so their
        // f64 conversion is exact; directed division encloses the weights.
        let weights = std::array::from_fn::<_, 22, _>(|i| {
            std::array::from_fn::<_, 22, _>(|a| {
                I::point((choose(21, i) * choose(21, a)) as f64)
                    / I::point(choose(42, i + a) as f64)
            })
        });
        let mut controls = vec![vec![I::point(0.); 43]; 43];
        for (i, row) in source.polynomial.controls.iter().enumerate() {
            for (j, p) in row.iter().enumerate() {
                for (a, row) in source.polynomial.controls.iter().enumerate() {
                    for (b, q) in row.iter().enumerate() {
                        let dot = if i == a && j == b {
                            p[0].square() + p[1].square() + p[2].square()
                        } else {
                            p[0] * q[0] + p[1] * q[1] + p[2] * q[2]
                        };
                        controls[i + a][j + b] =
                            controls[i + a][j + b] + weights[i][a] * weights[j][b] * dot;
                    }
                }
            }
        }
        if controls.iter().flatten().any(|x| !x.valid()) {
            return Err(Error::Unresolved);
        }
        let whole = hull(&controls);
        let derivative = std::array::from_fn(|axis| {
            let mut grid = controls.clone();
            std::array::from_fn(|order| {
                let value = hull(&grid);
                if order < 8 {
                    let degree = if axis == 0 {
                        grid.len() - 1
                    } else {
                        grid[0].len() - 1
                    };
                    for i in 0..grid.len() - usize::from(axis == 0) {
                        for j in 0..grid[0].len() - usize::from(axis == 1) {
                            let next = if axis == 0 {
                                grid[i + 1][j]
                            } else {
                                grid[i][j + 1]
                            };
                            grid[i][j] = (next - grid[i][j]) * I::point(degree as f64);
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
        let result = Self {
            source,
            controls,
            derivative,
            hull: whole,
        };
        if result.retained_budget_bytes() > NORMAL11_BYTES {
            return Err(Error::Limit);
        }
        Ok(result)
    }
    pub(crate) fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + 64
            + self.source().retained_budget_bytes()
            + self.controls.capacity() * std::mem::size_of::<Vec<I>>()
            + self
                .controls
                .iter()
                .map(|r| r.capacity() * std::mem::size_of::<I>())
                .sum::<usize>()
    }
    pub(crate) fn source(&self) -> &NormalPolynomial11 {
        &self.source
    }
    fn at(&self, uv: [[f64; 2]; 2]) -> Result<I, Error> {
        if !NormalPolynomial11::valid(uv) {
            return Err(Error::Invalid);
        }
        let bu = bernstein(42, I::new(uv[0][0], uv[0][1]));
        let bv = bernstein(42, I::new(uv[1][0], uv[1][1]));
        let anchor = self.controls[0][0];
        let mut delta = I::point(0.);
        for (i, row) in self.controls.iter().enumerate() {
            for (j, p) in row.iter().enumerate() {
                if i != 0 || j != 0 {
                    delta = delta + bu[i] * bv[j] * (*p - anchor);
                }
            }
        }
        (anchor + delta).checked()
    }
    /// Root-unit scalar field, intersected with its earned nonnegative range.
    pub(crate) fn bounds(&self, domain: [[f64; 2]; 2]) -> Result<I, Error> {
        if !NormalPolynomial11::valid(domain) {
            return Err(Error::Invalid);
        }
        let mid = domain.map(|d| d[0] + (d[1] - d[0]) * 0.5);
        let at = self.at(mid.map(|x| [x, x]))?;
        let delta = domain.map(|d| I::new(d[0], d[1]));
        let bound = at
            + self.derivative[0][1] * (delta[0] - I::point(mid[0]))
            + self.derivative[1][1] * (delta[1] - I::point(mid[1]));
        // N dot N >=0 is an algebraic identity, not an inferred tolerance or
        // epsilon. A zero lower bound still prevents a finite sqrt jet.
        I::new(
            bound.lo.max(self.hull.lo).max(0.),
            bound.hi.min(self.hull.hi),
        )
        .checked()
    }
    /// Coefficients of squared CHILD-unit normal length: (hu*hv)^2 and the
    /// selected axis chain powers are applied to original root derivatives.
    pub(crate) fn taylor8(&self, domain: [[f64; 2]; 2], along_v: bool) -> Result<[I; 9], Error> {
        let zero = self.bounds(domain)?;
        let width = domain.map(|d| I::point(d[1]) - I::point(d[0]));
        let mut scale = (width[0] * width[1]).square();
        let mut factorial = 1.;
        let axis = usize::from(along_v);
        Ok(std::array::from_fn(|order| {
            let root = if order == 0 {
                zero
            } else {
                self.derivative[axis][order]
            };
            let value = root * scale / I::point(factorial);
            scale = scale * width[axis];
            factorial *= (order + 1) as f64;
            value
        }))
    }
}
#[cfg(test)]
mod tests;
