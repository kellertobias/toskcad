//! Rectangular mixed Taylor coefficients of the original bilinear atlas.
//! Coefficients enclose source derivatives, never fitted surface derivatives.
use super::*;
const N: usize = 13;
type C = [I; 2];
fn add(a: C, b: C) -> C {
    [a[0] + b[0], a[1] + b[1]]
}
fn mul(a: C, b: C) -> C {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn scale(a: C, b: I) -> C {
    a.map(|x| x * b)
}
/// Entry[i][j] encloses d^(i+j)X / (du^i dv^j i!j!) at every
/// original chart center in the requested parameter rectangle.
#[derive(Clone, Copy, Debug)]
pub struct StockClipPatchTaylor {
    coefficients: [[[[f64; 2]; 3]; N]; N],
    work: usize,
    preparation_work: usize,
    bytes: usize,
}
impl StockClipPatchTaylor {
    pub fn coefficients_mm(&self) -> &[[[[f64; 2]; 3]; N]; N] {
        &self.coefficients
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Previously prepared borrowed law. Aggregate callers charge this once,
    /// then charge work_used() for EACH actual coefficient query.
    pub fn preparation_work_used(&self) -> usize {
        self.preparation_work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
// [value,du,dv,duv] at the original parameter center. The interval
// center is a family of centers; these are not divided differences.
fn scalar(c: [I; 4], u: I, v: I) -> [I; 4] {
    let a = c[1] - c[0];
    let b = c[2] - c[0];
    let d = c[3] - c[2] - c[1] + c[0];
    [c[0] + a * u + b * v + d * u * v, a + d * v, b + d * u, d]
}
impl StockClipPatchLaw<'_> {
    /// Exact source-field corner caches in (0,0),(1,0),(0,1),(1,1)
    /// order, each [radius,unwrapped phase,z]. They do not establish owner
    /// equality and retain the original symbolic contact enclosure.
    pub fn bilinear_fields_mm(&self) -> [[[f64; 2]; 3]; 4] {
        self.corners.map(|p| p.map(bounds))
    }
    pub fn taylor(
        &self,
        u: [f64; 2],
        v: [f64; 2],
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipPatchTaylor, Error> {
        // In addition to the returned array, the complex exponential and
        // radius/phase temporaries fit below16KiB of fixed stack workspace.
        let retained = self.bytes - SCRATCH + size_of::<StockClipPatchTaylor>();
        let mut meter = Meter::new(options, retained, interrupt)?;
        let u = sample::unit_parameter(u)?;
        let v = sample::unit_parameter(v)?;
        // Fixed recurrence: complex coefficients use <=2664 interval
        // operations, radius products/additions <=2162, scalar/phase setup
        // <128. Trig enclosure, bound publication and fixed loop overhead
        // together fit the remaining >3000 units. Prepared-law work is
        // reported separately; this call does not recompute its corners.
        meter.charge(8192)?;
        let [r, q, z]: [[I; 4]; 3] =
            std::array::from_fn(|j| scalar(self.corners.map(|p| p[j]), u, v));
        let omega = I::point(2.) * pi() / I::point(self.plan.source.spec.pitch_mm);
        let theta: [I; 4] = std::array::from_fn(|j| (z[j] - q[j]) * omega);
        let zero = I::point(0.);
        let mut exponential = [[[zero; 2]; N]; N];
        exponential[0][0] = [theta[0].cos(), theta[0].sin()];
        for j in 1..N {
            exponential[0][j] = scale(
                mul([zero, theta[2]], exponential[0][j - 1]),
                I::point(1.) / I::point(j as f64),
            );
        }
        // E_u = i(theta_u + theta_uv * delta_v) E. The mixed term
        // is retained explicitly; separable angular fits cannot omit it.
        for i in 1..N {
            for j in 0..N {
                let mut value = mul([zero, theta[1]], exponential[i - 1][j]);
                if j != 0 {
                    value = add(value, mul([zero, theta[3]], exponential[i - 1][j - 1]));
                }
                exponential[i][j] = scale(value, I::point(1.) / I::point(i as f64));
            }
        }
        let coefficients = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let mut xy = scale(exponential[i][j], r[0]);
                if i != 0 {
                    xy = add(xy, scale(exponential[i - 1][j], r[1]));
                }
                if j != 0 {
                    xy = add(xy, scale(exponential[i][j - 1], r[2]));
                }
                if i != 0 && j != 0 {
                    xy = add(xy, scale(exponential[i - 1][j - 1], r[3]));
                }
                let axial = match (i, j) {
                    (0, 0) => z[0],
                    (1, 0) => z[1],
                    (0, 1) => z[2],
                    (1, 1) => z[3],
                    _ => zero,
                };
                [bounds(xy[0]), bounds(xy[1]), bounds(axial)]
            })
        });
        // Sealed source r,z<=1e7, L/p<128, |q|<3p. Each normalized
        // scalar derivative is below1e9. Through(12,12), at most2^24
        // recurrence paths times(1e9)^24 remain below1e225. Interval
        // arithmetic and publication therefore cannot overflow on this
        // sealed domain. Requested geometric precision is a caller gate.
        meter.charge(1)?;
        Ok(StockClipPatchTaylor {
            coefficients,
            work: meter.used,
            preparation_work: self.work,
            bytes: retained + SCRATCH,
        })
    }
}
#[cfg(test)]
mod tests;
