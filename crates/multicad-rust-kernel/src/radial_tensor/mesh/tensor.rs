//! Canonical source cells with indexed, knot-aware patch lookup.
use super::*;
use crate::tapered_thread::{calculation_cell, normal_polynomial};
#[cfg(test)]
mod tests;
type Grid = Vec<Vec<[I; 3]>>;
pub(super) struct Tensor {
    pub face: usize,
    pub source_index: usize,
    pub knots: [Vec<f64>; 2],
    cells: Vec<Patch>,
}
pub(super) struct Patch {
    pub controls: Grid,
    pub normal: Grid,
    pub domain: [[f64; 2]; 2],
    normal_derivative: [[I; 3]; 2],
    normal_hull: [I; 3],
    curvature: [f64; 2],
    normal_jet: std::sync::OnceLock<Option<normal_jet::NormalJet>>,
}
impl Tensor {
    pub fn new(
        source: &axial::ChartTensor,
        source_index: usize,
        chart: &axial::Chart,
        work: &mut Work<'_>,
    ) -> Result<Self, Error> {
        let knots = source
            .spans
            .map(|n| (0..=n).map(|i| i as f64 / n as f64).collect::<Vec<_>>());
        let height = 9 * source.spans[1] + 1;
        let mut cells = Vec::with_capacity(source.spans[0] * source.spans[1]);
        let bottom = std::array::from_fn(|k| {
            I::point(chart.bottom_point[k]) - I::point(chart.axis_origin[k])
        });
        let station = dot(bottom, chart.axes[2]);
        for u in 0..source.spans[0] {
            for v in 0..source.spans[1] {
                work.charge(16384)?;
                let local = std::array::from_fn(|i| {
                    source.controls[(9 * u + i / 10) * height + 9 * v + i % 10]
                });
                let controls = local.map(|p| {
                    std::array::from_fn(|k| {
                        I::point(chart.axis_origin[k])
                            + chart.axes[0][k] * p[0]
                            + chart.axes[1][k] * p[1]
                            + chart.axes[2][k] * (station + p[2])
                    })
                });
                // Translation never enters the derivative calculation.
                // Proper orthonormal placement rotates the exact local
                // normal field, preserving its polynomial degree/direction.
                let cell = calculation_cell(local);
                let normal: Grid = normal_polynomial(&cell)
                    .coefficients()
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|p| {
                                std::array::from_fn(|k| {
                                    chart.axes[0][k] * p[0]
                                        + chart.axes[1][k] * p[1]
                                        + chart.axes[2][k] * p[2]
                                })
                            })
                            .collect()
                    })
                    .collect();
                let controls = controls
                    .as_chunks::<10>()
                    .0
                    .iter()
                    .map(|row| row.to_vec())
                    .collect::<Vec<_>>();
                let domain = [
                    [knots[0][u], knots[0][u + 1]],
                    [knots[1][v], knots[1][v + 1]],
                ];
                let mut normal_hull = [I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
                let mut normal_derivative = [[I::new(f64::INFINITY, f64::NEG_INFINITY); 3]; 2];
                let mut curvature = [0_f64; 2];
                for (u, row) in normal.iter().enumerate() {
                    for (v, &point) in row.iter().enumerate() {
                        for k in 0..3 {
                            normal_hull[k].lo = normal_hull[k].lo.min(point[k].lo);
                            normal_hull[k].hi = normal_hull[k].hi.max(point[k].hi);
                        }
                        for axis in 0..2 {
                            if (axis == 0 && u == 17) || (axis == 1 && v == 17) {
                                continue;
                            }
                            let next = if axis == 0 {
                                normal[u + 1][v]
                            } else {
                                normal[u][v + 1]
                            };
                            let width = (I::point(domain[axis][1]) - I::point(domain[axis][0]))
                                .positive()?;
                            for k in 0..3 {
                                let d = (next[k] - point[k]) * I::point(17.) / width;
                                normal_derivative[axis][k].lo =
                                    normal_derivative[axis][k].lo.min(d.lo);
                                normal_derivative[axis][k].hi =
                                    normal_derivative[axis][k].hi.max(d.hi);
                            }
                        }
                    }
                }
                for axis in 0..2 {
                    let width =
                        (I::point(domain[axis][1]) - I::point(domain[axis][0])).positive()?;
                    for u in 0..10 {
                        for v in 0..10 {
                            if (axis == 0 && u >= 8) || (axis == 1 && v >= 8) {
                                continue;
                            }
                            let at = |d: usize| {
                                if axis == 0 {
                                    controls[u + d][v]
                                } else {
                                    controls[u][v + d]
                                }
                            };
                            let (a, b, c) = (at(0), at(1), at(2));
                            let value = (0..3)
                                .fold(I::point(0.), |sum, k| {
                                    sum + I::point(
                                        ((a[k] - I::point(2.) * b[k] + c[k]) * I::point(72.)
                                            / width.square())
                                        .magnitude(),
                                    )
                                })
                                .checked()?
                                .hi;
                            curvature[axis] = curvature[axis].max(value);
                        }
                    }
                }
                for value in normal_derivative.into_iter().flatten().chain(normal_hull) {
                    value.checked()?;
                }
                cells.push(Patch {
                    controls,
                    normal,
                    domain,
                    normal_derivative,
                    normal_hull,
                    curvature,
                    normal_jet: std::sync::OnceLock::new(),
                });
            }
        }
        Ok(Self {
            face: source.face,
            source_index,
            knots,
            cells,
        })
    }
    pub fn patches(&self, range: [[f64; 2]; 2], work: &mut Work<'_>) -> Result<Vec<&Patch>, Error> {
        let indices = std::array::from_fn::<_, 2, _>(|axis| {
            let first = self.knots[axis]
                .partition_point(|&k| k < range[axis][0])
                .saturating_sub(1);
            let last = self.knots[axis]
                .partition_point(|&k| k <= range[axis][1])
                .min(self.knots[axis].len() - 1);
            first..last
        });
        let count = indices[0].len() * indices[1].len();
        work.charge(count)?;
        let mut found = Vec::with_capacity(count);
        for u in indices[0].clone() {
            for v in indices[1].clone() {
                let patch = &self.cells[u * (self.knots[1].len() - 1) + v];
                // A C0 edge is already included by each incident closed
                // patch. Do not demand the opposite one-sided normal at
                // a zero-area intersection across that same crease.
                if (0..2).all(|axis| {
                    range[axis][0] < patch.domain[axis][1] && range[axis][1] > patch.domain[axis][0]
                }) {
                    found.push(patch);
                }
            }
        }
        Ok(found)
    }
    /// Whole-source cached derivative estimates only steer refinement;
    /// final acceptance uses the source/facet certificates below.
    pub fn hint(
        &self,
        range: [[f64; 2]; 2],
        deflection: f64,
        angle: f64,
        work: &mut Work<'_>,
    ) -> Result<[f64; 2], Error> {
        let mut score = [0_f64; 2];
        for patch in self.patches(range, work)? {
            work.charge(16)?;
            let q = patch.normal_hull;
            let qscale = q
                .into_iter()
                .map(|v| v.magnitude())
                .fold(0_f64, f64::max)
                .powi(2)
                .max(f64::MIN_POSITIVE);
            for axis in 0..2 {
                let d = patch.normal_derivative[axis];
                let turn = [
                    q[1] * d[2] - q[2] * d[1],
                    q[2] * d[0] - q[0] * d[2],
                    q[0] * d[1] - q[1] * d[0],
                ]
                .into_iter()
                .map(|v| v.magnitude())
                .sum::<f64>();
                let width = range[axis][1] - range[axis][0];
                score[axis] = score[axis].max(
                    turn * width / (qscale * angle)
                        + patch.curvature[axis] * width.powi(2) / deflection,
                );
            }
        }
        Ok(score)
    }
    pub fn at_rectified(
        &self,
        body: &RadialBoundary,
        uv: [f64; 2],
        work: &mut Work<'_>,
    ) -> Result<[I; 3], Error> {
        work.charge(512)?;
        let axial = &body.ordered.caps.axial;
        let chart = &axial.chart;
        let source = &axial.tensors[self.source_index];
        let indices = std::array::from_fn::<_, 2, _>(|axis| {
            self.knots[axis]
                .partition_point(|&k| k <= uv[axis])
                .saturating_sub(1)
                .min(self.knots[axis].len() - 2)
        });
        let cell = &self.cells[indices[0] * (self.knots[1].len() - 1) + indices[1]];
        let t: [I; 2] = std::array::from_fn(|axis| {
            (I::point(uv[axis]) - I::point(cell.domain[axis][0]))
                / (I::point(cell.domain[axis][1]) - I::point(cell.domain[axis][0]))
        });
        let height = 9 * source.spans[1] + 1;
        let rows = (0..=9)
            .map(|u| {
                let row = (0..=9)
                    .map(|v| source.controls[(9 * indices[0] + u) * height + 9 * indices[1] + v])
                    .collect::<Vec<_>>();
                evaluate(&row, t[1])
            })
            .collect::<Vec<_>>();
        let p = evaluate(&rows, t[0]);
        let bottom = std::array::from_fn(|k| {
            I::point(chart.bottom_point[k]) - I::point(chart.axis_origin[k])
        });
        let station = dot(bottom, chart.axes[2]);
        let z = I::point(chart.nominal_height) * I::point(uv[1]);
        let world = std::array::from_fn(|k| {
            I::point(chart.axis_origin[k])
                + chart.axes[0][k] * p[0]
                + chart.axes[1][k] * p[1]
                + chart.axes[2][k] * (station + z)
        });
        for value in world {
            value.checked()?;
        }
        Ok(world)
    }
}
fn evaluate(points: &[[I; 3]], t: I) -> [I; 3] {
    let mut row = points.to_vec();
    let one = I::point(1.) - t;
    while row.len() > 1 {
        row = row
            .windows(2)
            .map(|p| std::array::from_fn(|k| one * p[0][k] + t * p[1][k]))
            .collect();
    }
    row[0]
}

impl Patch {
    pub fn taylor_normal_box(
        &self,
        range: [[f64; 2]; 2],
        work: &mut Work<'_>,
    ) -> Result<[I; 3], Error> {
        if self.normal_jet.get().is_none() {
            let candidate = match normal_jet::NormalJet::new(&self.normal, self.domain, work) {
                Ok(jet) => Some(jet),
                Err(Error::Unresolved) => None,
                Err(error) => return Err(error),
            };
            let _ = self.normal_jet.set(candidate);
        }
        self.normal_jet
            .get()
            .and_then(Option::as_ref)
            .ok_or(Error::Unresolved)?
            .enclose(range, work)
    }

    /// Interval point evaluation plus a whole-polynomial derivative bound
    /// is an authority-preserving mean-value enclosure, not a sampled
    /// normal approximation. The Bernstein cone remains the fallback.
    pub fn normal_box(&self, range: [[f64; 2]; 2], work: &mut Work<'_>) -> Result<[I; 3], Error> {
        work.charge(2048)?;
        let range = std::array::from_fn::<_, 2, _>(|axis| {
            [
                range[axis][0].max(self.domain[axis][0]),
                range[axis][1].min(self.domain[axis][1]),
            ]
        });
        let middle = range.map(|d| d[0] + (d[1] - d[0]) * 0.5);
        let local = std::array::from_fn::<_, 2, _>(|axis| {
            let t = (I::point(middle[axis]) - I::point(self.domain[axis][0]))
                / (I::point(self.domain[axis][1]) - I::point(self.domain[axis][0]));
            I::new(t.lo.max(0.), t.hi.min(1.))
        });
        let weights = local.map(weights);
        let anchor = self.normal[0][0];
        let mut center = anchor;
        for (u, row) in self.normal.iter().enumerate() {
            for (v, &point) in row.iter().enumerate() {
                if u == 0 && v == 0 {
                    continue;
                }
                let weight = weights[0][u] * weights[1][v];
                for k in 0..3 {
                    center[k] = center[k] + weight * (point[k] - anchor[k]);
                }
            }
        }
        let delta = std::array::from_fn::<_, 2, _>(|axis| {
            I::new(range[axis][0], range[axis][1]) - I::point(middle[axis])
        });
        let result = std::array::from_fn(|k| {
            center[k]
                + self.normal_derivative[0][k] * delta[0]
                + self.normal_derivative[1][k] * delta[1]
        });
        for value in result {
            value.checked()?;
        }
        Ok(result)
    }
}
fn weights(t: I) -> [I; 18] {
    let one = I::point(1.) - t;
    let one = I::new(one.lo.max(0.), one.hi.min(1.));
    let power = |x: I, n: usize| (0..n).fold(I::point(1.), |v, _| v * x);
    std::array::from_fn(|i| I::point(trims::binomial(17, i)) * power(t, i) * power(one, 17 - i))
}
