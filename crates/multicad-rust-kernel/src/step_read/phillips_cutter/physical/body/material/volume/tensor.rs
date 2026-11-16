//! Exact directed Bernstein integral of the ACTUAL new polynomial patch.
//! For F=(x/2,y/2,0), side flux is(x*Nx+y*Ny)/2. N has degree17 in
//! each variable; multiplying degree9 P gives degree26. The separable weights
//! are integral(B_i^9 B_j^17)=C(9,i)C(17,j)/(27*C(26,i+j)). Integer
//! numerators<=3063060 and denominators<=280816200 are exact in binary64.
//! Each cell uses its OWN unit-square derivatives; all physical spans are
//! traversed once, without an extra span-count/Jacobian multiplier.
use super::*;
fn choose(n: usize, k: usize) -> u64 {
    let mut v = 1u64;
    for i in 0..k {
        v = v * (n - i) as u64 / (i + 1) as u64;
    }
    v
}
pub(super) fn weights() -> [[I; 18]; 10] {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let n = choose(9, i) * choose(17, j);
            let d = 27 * choose(26, i + j);
            I::point(n as f64) / I::point(d as f64)
        })
    })
}
pub(super) fn cell(p: [[f64; 3]; 100], w: &[[I; 18]; 10]) -> I {
    let geometry = crate::tapered_thread::calculation_cell(p.map(|p| p.map(I::point)));
    let normal = crate::tapered_thread::normal_polynomial(&geometry);
    let c = normal.coefficients();
    let mut volume = I::point(0.);
    for axis in 0..2 {
        let a = std::array::from_fn::<_, 18, _>(|u| {
            std::array::from_fn::<_, 10, _>(|v| {
                (0..10).fold(I::point(0.), |s, i| {
                    s + I::point(p[10 * i + v][axis]) * w[i][u]
                })
            })
        });
        for u in 0..18 {
            for v in 0..18 {
                let integrated = (0..10).fold(I::point(0.), |s, j| s + a[u][j] * w[j][v]);
                volume = volume + integrated * c[u][v][axis];
            }
        }
    }
    volume / I::point(2.)
}
