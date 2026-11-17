//! Exact Bernstein Green integrals of actual planar cap boundaries.
use super::*;
pub(super) const SETUP_WORK: usize = 262144;
pub(super) const CURVE_WORK: usize = 32768;
pub(super) const SCRATCH_BYTES: usize = 128 * 1024;
pub(super) struct Weights {
    area: [[I; 11]; 12],
    moment: [[[I; 11]; 12]; 12],
}
fn choose(n: usize, k: usize) -> u64 {
    (0..k.min(n - k)).fold(1u64, |v, i| v * (n - i) as u64 / (i + 1) as u64)
}
impl Weights {
    pub(super) fn new() -> Self {
        Self {
            area: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    // Integral B_i^11 B_j^10 on [0,1]. Both integers are exact.
                    I::point((choose(11, i) * choose(10, j)) as f64)
                        / I::point((22 * choose(21, i + j)) as f64)
                })
            }),
            moment: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    std::array::from_fn(|k| {
                        // Integral B_i^11 B_j^11 B_k^10. All numerator/denominator
                        // integers are below2^53 before outward interval division.
                        I::point((choose(11, i) * choose(11, j) * choose(10, k)) as f64)
                            / I::point((33 * choose(32, i + j + k)) as f64)
                    })
                })
            }),
        }
    }
    /// Signed area and signed first moments in X/Y. Closure and orientation
    /// are supplied only by the caller's earned complete actual cap wire.
    pub(super) fn curve(&self, p: [[I; 2]; 12]) -> [I; 3] {
        let d = std::array::from_fn::<_, 11, _>(|j| {
            std::array::from_fn::<_, 2, _>(|k| (p[j + 1][k] - p[j][k]) * I::point(11.))
        });
        let mut result = [I::point(0.); 3];
        for (i, a) in p.iter().enumerate() {
            for (j, b) in d.iter().enumerate() {
                result[0] =
                    result[0] + self.area[i][j] * (a[0] * b[1] - a[1] * b[0]) / I::point(2.);
            }
            for (j, b) in p.iter().enumerate() {
                for (k, d) in d.iter().enumerate() {
                    let w = self.moment[i][j][k] / I::point(2.);
                    result[1] = result[1] + w * a[0] * b[0] * d[1];
                    result[2] = result[2] - w * a[1] * b[1] * d[0];
                }
            }
        }
        result
    }
}
