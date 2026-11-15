//! Distinct degree21 normal coefficients of ORIGINAL degree11 tensors.
//! Products are formed before restriction, preserving shared helical terms.
use super::*;
use std::sync::Arc;
#[derive(Clone)]
pub(super) struct Normal {
    coefficients: Arc<[[I; 3]]>,
    pub position_hessian: hessian::Hessian,
    pub hessian: hessian::Hessian,
    pub steering: [f64; 2],
    pub normal_aligned: [[I; 3]; 3],
    pub position_aligned: [[I; 3]; 3],
}
fn choose(n: usize, k: usize) -> f64 {
    (0..k.min(n - k)).fold(1u64, |value, i| value * (n - i) as u64 / (i + 1) as u64) as f64
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
impl Normal {
    pub fn coefficients(&self) -> &[[I; 3]] {
        &self.coefficients
    }
    #[cfg(test)]
    pub fn point(&self, uv: [f64; 2], budget: &mut Budget<'_>) -> Result<[I; 3]> {
        point_values(&self.coefficients, uv, budget)
    }
    pub fn new(
        source: &OwnedMaterial,
        tensor: &Tensor,
        u: usize,
        v: usize,
        budget: &mut Budget<'_>,
    ) -> Result<Self> {
        //17,424 paired cross products and degree21 Bernstein product weights,
        // then484 cofactor maps. This preparation is charged once per root.
        budget.step(524_288)?;
        let height = 11 * tensor.spans[1] + 1;
        let p = std::array::from_fn::<_, 144, _>(|i| {
            interval_point(
                source,
                tensor.poles[(11 * u + i / 12) * height + 11 * v + i % 12],
            )
        });
        let du = std::array::from_fn::<_, 132, _>(|i| {
            std::array::from_fn::<_, 3, _>(|k| {
                (p[(i / 12 + 1) * 12 + i % 12][k] - p[i][k]) * I::point(11.)
            })
        });
        let dv = std::array::from_fn::<_, 132, _>(|i| {
            std::array::from_fn::<_, 3, _>(|k| {
                (p[(i / 11) * 12 + i % 11 + 1][k] - p[(i / 11) * 12 + i % 11][k]) * I::point(11.)
            })
        });
        let wu = std::array::from_fn::<_, 11, _>(|i| {
            std::array::from_fn::<_, 12, _>(|a| {
                I::point(choose(10, i)) * I::point(choose(11, a)) / I::point(choose(21, i + a))
            })
        });
        let wv = std::array::from_fn::<_, 12, _>(|j| {
            std::array::from_fn::<_, 11, _>(|b| {
                I::point(choose(11, j)) * I::point(choose(10, b)) / I::point(choose(21, j + b))
            })
        });
        let mut coefficients = budget.allocate(484)?;
        coefficients.resize(484, [I::point(0.); 3]);
        for (i, row) in du.as_chunks::<12>().0.iter().enumerate() {
            for (j, &a) in row.iter().enumerate() {
                for (c, column) in dv.as_chunks::<11>().0.iter().enumerate() {
                    for (d, &b) in column.iter().enumerate() {
                        let n = cross(b, a);
                        let weight = wu[i][c] * wv[j][d];
                        let slot = &mut coefficients[(i + c) * 22 + j + d];
                        for k in 0..3 {
                            slot[k] = slot[k] + weight * n[k];
                        }
                    }
                }
            }
        }
        let a = source.0.linear;
        let cofactor = [cross(a[1], a[2]), cross(a[2], a[0]), cross(a[0], a[1])];
        for p in &mut coefficients {
            *p = std::array::from_fn(|k| {
                (0..3).fold(I::point(0.), |sum, j| sum + cofactor[j][k] * p[j])
            });
        }
        let geometry = crate::tapered_thread::calculation_cell11(p)
            .expect("sealed literal normal source cell");
        budget.step(4096)?;
        let midpoint = geometry
            .derivative_bounds_at([[0.5, 0.5]; 2])
            .expect("sealed source midpoint");
        let du = cells::vector(source, midpoint.du).map(I::midpoint);
        let dv = cells::vector(source, midpoint.dv).map(I::midpoint);
        let dot = |a: [f64; 3], b: [f64; 3]| (0..3).map(|k| a[k] * b[k]).sum::<f64>();
        let a = dot(du, du).sqrt();
        let b = dot(du, dv) / a;
        let c = (dot(dv, dv) - b * b).max(0.).sqrt();
        if !a.is_finite() || !b.is_finite() || !c.is_finite() || a <= 0. || c <= 0. {
            return Err(Error::Unsupported);
        }
        let steering = [
            ((b / a).clamp(-16., 16.) * 256.).round() / 256.,
            2f64.powi((c / a).log2().round().clamp(-32., 32.) as i32),
        ];
        let position_hessian = hessian::Hessian::new(11, &p, Some(source.0.linear), budget)?;
        let hessian = hessian::Hessian::new(21, &coefficients, None, budget)?;
        let normal_aligned = hessian.aligned(steering[0], budget)?;
        let position_aligned = position_hessian.aligned(steering[0], budget)?;
        hold(budget, 32)?;
        Ok(Self {
            coefficients: coefficients.into(),
            position_hessian,
            hessian,
            steering,
            normal_aligned,
            position_aligned,
        })
    }
}
pub(super) fn point_values(
    coefficients: &[[I; 3]],
    uv: [f64; 2],
    budget: &mut Budget<'_>,
) -> Result<[I; 3]> {
    budget.step(12_288)?;
    let weights = uv.map(basis::<22>);
    // Partition of unity retains the source anchor once, reducing
    // independent weight-rounding inflation around a large normal.
    let anchor = coefficients[0];
    let mut point = anchor;
    for i in 0..22 {
        for j in 0..22 {
            let w = weights[0][i] * weights[1][j];
            for k in 0..3 {
                point[k] = point[k] + w * (coefficients[i * 22 + j][k] - anchor[k]);
            }
        }
    }
    Ok(point)
}
pub(super) fn prepare(source: &OwnedMaterial, budget: &mut Budget<'_>) -> Result<Vec<Vec<Normal>>> {
    let mut cache = budget.allocate(source.tensors().len())?;
    for tensor in source.tensors() {
        let mut roots = budget.allocate(tensor.spans[0] * tensor.spans[1])?;
        for u in 0..tensor.spans[0] {
            for v in 0..tensor.spans[1] {
                roots.push(Normal::new(source, tensor, u, v, budget)?);
            }
        }
        cache.push(roots);
    }
    Ok(cache)
}
