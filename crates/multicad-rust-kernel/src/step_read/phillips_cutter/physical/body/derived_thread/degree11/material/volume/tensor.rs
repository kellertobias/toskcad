//! Directed exact Bernstein flux of the original degree11 unit-square cell.
//! The derivative nets already include every source-to-cell chain factor.
//! No additional global span multiplier belongs in this integral.
use super::*;
type Weights = [[I; 22]; 12];
fn choose(n: usize, k: usize) -> u64 {
    (0..k.min(n - k)).fold(1u64, |v, i| v * (n - i) as u64 / (i + 1) as u64)
}
pub(super) struct Quadrature {
    #[cfg(test)]
    product: [[I; 12]; 11],
    integral: Weights,
}
impl Quadrature {
    pub(super) fn new() -> Self {
        Self {
            #[cfg(test)]
            product: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    I::point((choose(10, i) * choose(11, j)) as f64)
                        / I::point(choose(21, i + j) as f64)
                })
            }),
            // integral B_i^11 B_j^21 = C(11,i)C(21,j)/(33*C(32,i+j)).
            // Every integer numerator/denominator is below2^53 and exact.
            integral: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    I::point((choose(11, i) * choose(21, j)) as f64)
                        / I::point((33 * choose(32, i + j)) as f64)
                })
            }),
        }
    }
}
#[cfg(test)]
pub(super) fn cell(p: [[I; 3]; 144], q: &Quadrature) -> I {
    let du: [[[I; 3]; 12]; 11] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            std::array::from_fn(|k| (p[12 * (i + 1) + j][k] - p[12 * i + j][k]) * I::point(11.))
        })
    });
    let dv: [[[I; 3]; 11]; 12] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            std::array::from_fn(|k| (p[12 * i + j + 1][k] - p[12 * i + j][k]) * I::point(11.))
        })
    });
    let mut normal = [[[I::point(0.); 3]; 22]; 22];
    for (i, row) in du.iter().enumerate() {
        for (j, u) in row.iter().enumerate() {
            for (a, row) in dv.iter().enumerate() {
                for (b, v) in row.iter().enumerate() {
                    let weight = q.product[i][a] * q.product[b][j];
                    for k in 0..3 {
                        let cross =
                            u[(k + 1) % 3] * v[(k + 2) % 3] - u[(k + 2) % 3] * v[(k + 1) % 3];
                        normal[i + a][j + b][k] = normal[i + a][j + b][k] + weight * cross;
                    }
                }
            }
        }
    }
    contract(p, q, &normal.map(|row| row.to_vec()))
}
pub(super) fn cached_cell(
    p: [[I; 3]; 144],
    q: &Quadrature,
    normal: &crate::tapered_thread::NormalPolynomial11,
) -> I {
    contract(p, q, normal.coefficients())
}
fn contract(p: [[I; 3]; 144], q: &Quadrature, normal: &[Vec<[I; 3]>]) -> I {
    let mut value = I::point(0.);
    for k in 0..2 {
        let a: [[I; 12]; 22] = std::array::from_fn(|u| {
            std::array::from_fn(|v| {
                (0..12).fold(I::point(0.), |sum, i| {
                    sum + p[12 * i + v][k] * q.integral[i][u]
                })
            })
        });
        for (u, normal_row) in normal.iter().enumerate() {
            for (v, normal) in normal_row.iter().enumerate() {
                let integrated =
                    (0..12).fold(I::point(0.), |sum, j| sum + a[u][j] * q.integral[j][v]);
                value = value + integrated * normal[k];
            }
        }
    }
    value / I::point(2.)
}
