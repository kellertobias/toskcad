//! Coupled directional derivatives of unchanged tensor polynomials.
//! Elevation aligns the three original partial fields before their linear
//! combination. Taking separate partial hulls would lose shear cancellation.
use super::*;
use std::sync::Arc;
#[derive(Clone)]
pub(super) struct Hessian {
    pub coefficients: [Arc<[[I; 3]]>; 3],
}
fn derivative(p: &[[I; 3]], dims: [usize; 2], axis: usize) -> (Vec<[I; 3]>, [usize; 2]) {
    let mut next = dims;
    next[axis] -= 1;
    let mut out = Vec::with_capacity(next[0] * next[1]);
    for u in 0..next[0] {
        for v in 0..next[1] {
            let a = p[u * dims[1] + v];
            let b = p[(u + usize::from(axis == 0)) * dims[1] + v + usize::from(axis == 1)];
            out.push(std::array::from_fn(|k| {
                (b[k] - a[k]) * I::point((dims[axis] - 1) as f64)
            }));
        }
    }
    (out, next)
}
fn elevate(p: &[[I; 3]], dims: [usize; 2], axis: usize) -> (Vec<[I; 3]>, [usize; 2]) {
    let mut next = dims;
    next[axis] += 1;
    let mut out = Vec::with_capacity(next[0] * next[1]);
    for u in 0..next[0] {
        for v in 0..next[1] {
            let i = if axis == 0 { u } else { v };
            let at = |i| {
                if axis == 0 {
                    p[i * dims[1] + v]
                } else {
                    p[u * dims[1] + i]
                }
            };
            let point = if i == 0 {
                at(0)
            } else if i == next[axis] - 1 {
                at(i - 1)
            } else {
                let a = at(i);
                let b = at(i - 1);
                let weight = I::point(i as f64) / I::point((next[axis] - 1) as f64);
                std::array::from_fn(|k| a[k] + weight * (b[k] - a[k]))
            };
            out.push(point);
        }
    }
    (out, next)
}
impl Hessian {
    pub fn new(
        degree: usize,
        p: &[[I; 3]],
        linear: Option<[[I; 3]; 3]>,
        budget: &mut Budget<'_>,
    ) -> Result<Self> {
        assert_eq!(
            p.len(),
            (degree + 1) * (degree + 1),
            "sealed original tensor polynomial"
        );
        let mut fields = Vec::with_capacity(3);
        for axes in [[0, 0], [0, 1], [1, 1]] {
            let (first, dims) = derivative(p, [degree + 1; 2], axes[0]);
            let (mut field, mut dims) = derivative(&first, dims, axes[1]);
            for axis in 0..2 {
                while dims[axis] < degree + 1 {
                    (field, dims) = elevate(&field, dims, axis);
                }
            }
            if let Some(a) = linear {
                for point in &mut field {
                    *point = std::array::from_fn(|k| {
                        (0..3).fold(I::point(0.), |sum, j| sum + a[j][k] * point[j])
                    });
                }
            }
            let mut retained = budget.allocate(field.len())?;
            retained.extend_from_slice(&field);
            hold(budget, 32)?;
            fields.push(Arc::from(retained));
        }
        Ok(Self {
            coefficients: fields
                .try_into()
                .expect("three actual Hessian partial fields"),
        })
    }
    pub fn weights(edges: [[I; 2]; 2]) -> [[I; 3]; 3] {
        let [a, b] = edges;
        [
            [a[0].square(), I::point(2.) * a[0] * a[1], a[1].square()],
            [a[0] * b[0], a[0] * b[1] + a[1] * b[0], a[1] * b[1]],
            [b[0].square(), I::point(2.) * b[0] * b[1], b[1].square()],
        ]
    }
    pub fn point(&self, index: usize, weight: [I; 3]) -> [I; 3] {
        std::array::from_fn(|k| {
            (0..3).fold(I::point(0.), |sum, j| {
                sum + weight[j] * self.coefficients[j][index][k]
            })
        })
    }
    pub fn aligned(&self, shear: f64, budget: &mut Budget<'_>) -> Result<[[I; 3]; 3]> {
        budget.step(self.coefficients[0].len() * 96)?;
        let s = I::point(shear);
        let weights = [
            [I::point(1.), I::point(0.), I::point(0.)],
            [-s, I::point(1.), I::point(0.)],
            [s.square(), -I::point(2.) * s, I::point(1.)],
        ];
        Ok(weights.map(|weight| {
            let mut hull = self.point(0, weight);
            for index in 1..self.coefficients[0].len() {
                let p = self.point(index, weight);
                for k in 0..3 {
                    hull[k] = I::new(hull[k].lo.min(p[k].lo), hull[k].hi.max(p[k].hi));
                }
            }
            hull
        }))
    }
    pub fn position_error(&self, weights: [[I; 3]; 3], budget: &mut Budget<'_>) -> Result<f64> {
        budget.step(self.coefficients[0].len() * 96)?;
        let maximum = weights.map(|weight| {
            (0..self.coefficients[0].len())
                .map(|i| cells::l1(self.point(i, weight)))
                .fold(0f64, f64::max)
        });
        Ok(
            ((I::point(maximum[0]) + I::point(2.) * I::point(maximum[1]) + I::point(maximum[2]))
                / I::point(8.))
            .hi,
        )
    }
}
