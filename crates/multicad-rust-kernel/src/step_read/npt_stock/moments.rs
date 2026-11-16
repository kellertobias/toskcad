//! Fixed Bernstein product integrals. These weights have no source recipe,
//! geometry, closure or occupancy authority of their own.
use super::*;
use crate::tapered_thread::interval::I;
#[derive(Clone, Copy)]
pub(super) struct Moments {
    pub pair: [[I; 12]; 12],
    pub skew: [[I; 12]; 12],
}
fn binomials<const N: usize>() -> [f64; N] {
    let mut values = [0.; N];
    values[0] = 1.;
    for i in 1..N {
        values[i] = values[i - 1] * (N - i) as f64 / i as f64;
    }
    values
}
impl Moments {
    pub fn new(budget: &mut Budget<'_>) -> Result<Self> {
        //144 rational product weights and66 independent skew entries,
        // bounded small-integer binomials and publication fit4096 units.
        budget.step(4096)?;
        // All table integers and intermediate products are below2^24;
        // their binary64 arithmetic is exact. Only rational division and
        // subsequent matrix arithmetic need outward intervals.
        let c11 = binomials::<12>();
        let c10 = binomials::<11>();
        let c21 = binomials::<22>();
        let c22 = binomials::<23>();
        let product = |i: usize, j: usize| I::point(c11[i] * c10[j]) / I::point(22. * c21[i + j]);
        let derivative = |i: usize, j: usize| {
            let lower = if j == 0 {
                I::point(0.)
            } else {
                product(i, j - 1)
            };
            let upper = if j == 11 { I::point(0.) } else { product(i, j) };
            I::point(11.) * (lower - upper)
        };
        let mut pair = [[I::point(0.); 12]; 12];
        let mut skew = pair;
        for i in 0..12 {
            for j in i..12 {
                let value = I::point(c11[i] * c11[j]) / I::point(23. * c22[i + j]);
                pair[i][j] = value;
                pair[j][i] = value;
                if i != j {
                    let value = derivative(i, j) - derivative(j, i);
                    skew[i][j] = value;
                    skew[j][i] = -value;
                }
            }
        }
        Ok(Self { pair, skew })
    }
    /// Exact product integral ∫(x*y_u-y*x_u) du dv of the supplied
    /// literal tensor coefficients, enclosed by directed arithmetic.
    pub fn tensor_angular(&self, poles: &[[[I; 3]; 12]; 12], budget: &mut Budget<'_>) -> Result<I> {
        budget.step(8192)?;
        let mut derivative = [[I::point(0.); 12]; 12];
        for (i, row) in derivative.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                for (k, pole) in poles.iter().enumerate() {
                    *value = *value + self.skew[i][k] * pole[j][1];
                }
            }
        }
        let mut weighted = [[I::point(0.); 12]; 12];
        for (i, row) in weighted.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                for (k, weight) in self.pair[j].iter().enumerate() {
                    *value = *value + derivative[i][k] * *weight;
                }
            }
        }
        let mut integral = I::point(0.);
        for (i, row) in weighted.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                integral = integral + poles[i][j][0] * *value;
            }
        }
        Ok(integral)
    }
    /// Oriented planar boundary integral ∫(x*y'-y*x') dt.
    pub fn curve_angular(&self, poles: &[[I; 3]; 12], budget: &mut Budget<'_>) -> Result<I> {
        budget.step(1024)?;
        let mut integral = I::point(0.);
        for (i, a) in poles.iter().enumerate() {
            for (j, b) in poles.iter().enumerate() {
                integral = integral + a[0] * self.skew[i][j] * b[1];
            }
        }
        Ok(integral)
    }
}
#[cfg(test)]
mod tests;
