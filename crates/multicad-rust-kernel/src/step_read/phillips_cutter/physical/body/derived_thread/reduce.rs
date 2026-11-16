//! Directed degree elevation for measured repaired polynomial proposals.
use super::*;
fn choose(n: usize, k: usize) -> f64 {
    (0..k).fold(1., |a, i| a * (n - i) as f64 / (i + 1) as f64)
}
// Proposal weights are fixed rational endpoint-jet ratios. Their binary64
// values propose controls only; the separate directed elevation proves error.
const PROPOSE: [[f64; 5]; 5] = [
    [1., 0., 0., 0., 0.],
    [1., 11. / 9., 0., 0., 0.],
    [1., 22. / 9., 55. / 36., 0., 0.],
    [1., 11. / 3., 55. / 12., 55. / 28., 0.],
    [1., 44. / 9., 55. / 6., 55. / 7., 55. / 21.],
];
pub(super) fn reduce(p: [[f64; 3]; 12]) -> [[f64; 3]; 10] {
    let mut q = [[0.; 3]; 10];
    for reverse in [false, true] {
        let mut d = std::array::from_fn::<_, 12, _>(|i| p[if reverse { 11 - i } else { i }]);
        let mut jets = [[0.; 3]; 5];
        for (k, jet) in jets.iter_mut().enumerate() {
            *jet = d[0];
            for i in 0..11 - k {
                d[i] = std::array::from_fn(|a| d[i + 1][a] - d[i][a]);
            }
        }
        for i in 0..5 {
            q[if reverse { 9 - i } else { i }] =
                std::array::from_fn(|a| (0..=i).map(|k| PROPOSE[i][k] * jets[k][a]).sum());
        }
    }
    q
}
// Source-independent exact degree-elevation ratios, computed once per query.
// All binomials/products here are integers <2^53. Division stays directed.
pub(super) struct Elevation([[I; 10]; 12]);
impl Elevation {
    pub fn new() -> Self {
        Self(std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                if j <= i && i - j <= 2 {
                    I::point(choose(9, j) * choose(2, i - j)) / I::point(choose(11, i))
                } else {
                    I::point(0.)
                }
            })
        }))
    }
    pub fn apply(&self, p: [[I; 3]; 10]) -> [[I; 3]; 12] {
        std::array::from_fn(|i| {
            std::array::from_fn(|a| {
                (i.saturating_sub(2)..=9.min(i))
                    .fold(I::point(0.), |s, j| s + p[j][a] * self.0[i][j])
            })
        })
    }
}
pub(super) fn split(mut p: [[I; 3]; 12]) -> ([[I; 3]; 12], [[I; 3]; 12]) {
    let mut a = [[I::point(0.); 3]; 12];
    let mut b = a;
    for k in 0..12 {
        a[k] = p[0];
        b[11 - k] = p[11 - k];
        for i in 0..11 - k {
            p[i] = std::array::from_fn(|axis| (p[i][axis] + p[i + 1][axis]) * I::point(0.5));
        }
    }
    (a, b)
}
pub(super) fn l1(p: [I; 3]) -> f64 {
    p.into_iter()
        .fold(I::point(0.), |s, x| s + I::point(x.magnitude()))
        .hi
}
