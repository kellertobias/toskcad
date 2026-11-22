//! Linear-time complex recurrence for the SAME affine original edge law.
//! This additive query reports borrowed preparation separately; the existing
//! general Series query and its published work semantics remain unchanged.
use super::*;
type Complex = [I; 2];
fn mul(a: Complex, b: Complex) -> Complex {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
#[derive(Clone, Copy, Debug)]
pub struct StockClipAffineCurveTaylor {
    coefficients: [[[f64; 2]; 3]; 13],
    norms: [f64; 13],
    work: usize,
    preparation_work: usize,
    bytes: usize,
}
impl StockClipAffineCurveTaylor {
    pub fn coefficients_mm(&self) -> [[[f64; 2]; 3]; 13] {
        self.coefficients
    }
    pub fn coefficient_l1_bounds_mm(&self) -> &[f64; 13] {
        &self.norms
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn preparation_work_used(&self) -> usize {
        self.preparation_work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl StockClipCurveLaw<'_> {
    pub fn affine_taylor(
        &self,
        parameter: [f64; 2],
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipAffineCurveTaylor, Error> {
        let retained = self.bytes - SCRATCH + size_of::<StockClipAffineCurveTaylor>();
        let mut meter = Meter::new(options, retained, interrupt)?;
        let t = sample::unit_parameter(parameter)?;
        // The affine phase recurrence uses12complex products, each6scalar
        // interval operations plus2divisions; coordinate Leibniz products,
        // final unit rotations and L1bounds use<400 more. Field evaluation,
        // publication and fixed arrays fit the remaining>500 units.
        meter.charge(1024)?;
        let delta: [I; 3] =
            std::array::from_fn(|axis| self.endpoints[1][axis] - self.endpoints[0][axis]);
        let at: [I; 3] = std::array::from_fn(|axis| self.endpoints[0][axis] + delta[axis] * t);
        let omega = I::point(2.) * pi() / I::point(self.plan.source.spec.pitch_mm);
        let theta = (at[2] - at[1]) * omega;
        let rate = (delta[2] - delta[1]) * omega;
        let zero = I::point(0.);
        let rotation = [theta.cos(), theta.sin()];
        let mut exponential = [[zero; 2]; 13];
        exponential[0] = [I::point(1.), zero];
        for n in 1..13 {
            exponential[n] = mul([zero, rate], exponential[n - 1]).map(|x| x / I::point(n as f64));
        }
        let sqrt_two = I::point(2.).sqrt();
        let mut norms = [0.; 13];
        let coefficients = std::array::from_fn(|n| {
            let mut xy = exponential[n].map(|x| x * at[0]);
            if n > 0 {
                for (axis, x) in xy.iter_mut().enumerate() {
                    *x = *x + delta[0] * exponential[n - 1][axis];
                }
            }
            let z = match n {
                0 => at[2],
                1 => delta[2],
                _ => zero,
            };
            norms[n] = (sqrt_two * (I::point(xy[0].magnitude()) + I::point(xy[1].magnitude()))
                + I::point(z.magnitude()))
            .hi;
            let xy = mul(rotation, xy);
            [bounds(xy[0]), bounds(xy[1]), bounds(z)]
        });
        // The sealed original edge has |theta'|<900 and r,z<=1e7.
        // All coefficients/intermediates through12 remain below1e50.
        meter.charge(1)?;
        Ok(StockClipAffineCurveTaylor {
            coefficients,
            norms,
            work: meter.used,
            preparation_work: self.work,
            bytes: retained + SCRATCH,
        })
    }
}
#[cfg(test)]
mod tests;
