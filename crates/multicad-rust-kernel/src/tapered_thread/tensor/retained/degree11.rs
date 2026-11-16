//! Fixed degree11 numerical cells, distinct from the degree9 source machinery.
//! No radial chart, occupied boundary, material, normal or mass is certified.
use super::*;
mod normal;
<<<<<<< HEAD
pub use normal::{NORMAL11_BYTES, NORMAL11_WORK, NormalPolynomial11};
mod surface;
pub use surface::{SurfaceFailure11, SurfaceOptions11, SurfaceProperties11, SurfaceStage11};
=======
pub use normal::{
    NORMAL11_BYTES, NORMAL11_WORK, NormalPolynomial11, ReducedNormal11, ReducedNormalCell11,
};
>>>>>>> d1743a74 (feat: enclose actual normals with certified elevation residuals)
const DEGREE: usize = 11;
#[derive(Clone, Debug, PartialEq)]
pub struct RetainedTensor11 {
    spans: [usize; 2],
    poles: Arc<[Point]>,
}
/// Degree11 only. Its private wrapper cannot enter degree9 normal helpers.
#[derive(Clone, Debug)]
pub struct TensorCell11(TensorCell);
/// Whole-leaf position and first derivatives in THIS leaf's unit square.
/// No second-derivative grids or arbitrary-domain restrictions are computed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorFirstDerivatives11 {
    pub position: [[f64; 2]; 3],
    pub du: [[f64; 2]; 3],
    pub dv: [[f64; 2]; 3],
}
fn grid(spans: [usize; 2], actual: usize) -> Result<[usize; 2], Error> {
    if spans.contains(&0) {
        return Err(Error::Invalid);
    }
    let dimensions = spans.map(|n| DEGREE.checked_mul(n).and_then(|x| x.checked_add(1)));
    let [Some(width), Some(height)] = dimensions else {
        return Err(Error::Limit);
    };
    let count = width.checked_mul(height).ok_or(Error::Limit)?;
    if count > MAX_POLES {
        return Err(Error::Limit);
    }
    if count != actual {
        return Err(Error::Invalid);
    }
    Ok([width, height])
}
impl RetainedTensor11 {
    /// Retain defining source-mm coefficients. Callers still own aggregate
    /// source/work/storage policy; this is the same524288 per-net cap as D9.
    pub fn new(spans: [usize; 2], poles: Vec<Point>) -> Result<Self, Error> {
        grid(spans, poles.len())?;
        if poles
            .iter()
            .flatten()
            .any(|x| !x.is_finite() || x.abs() > 1e7)
        {
            return Err(Error::Invalid);
        }
        Ok(Self {
            spans,
            poles: poles.into(),
        })
    }
    pub fn spans(&self) -> [usize; 2] {
        self.spans
    }
    pub fn dimensions(&self) -> [usize; 2] {
        self.spans.map(|n| DEGREE * n + 1)
    }
    pub fn poles(&self) -> &[Point] {
        &self.poles
    }
    pub fn cell(&self, index: [usize; 2]) -> Result<TensorCell11, Error> {
        if index.into_iter().zip(self.spans).any(|(i, n)| i >= n) {
            return Err(Error::Invalid);
        }
        let height = self.dimensions()[1];
        let controls = std::array::from_fn(|i| {
            self.poles[(11 * index[0] + i / 12) * height + 11 * index[1] + i % 12].map(I::point)
        });
        crate::tapered_thread::calculation_cell11(controls)
    }
}
impl TensorCell11 {
    pub fn domain(&self) -> [[f64; 2]; 2] {
        self.0.domain()
    }
    /// Actual interval Bernstein corner coefficients, counterclockwise in
    /// the cell's (u,v) chart: 00,10,11,01. Subdivision's directed transfers
    /// remain included; this is calculation data, not a solid admission.
    pub fn corner_bounds(&self) -> [[[f64; 2]; 3]; 4] {
        let last = self.0.controls.len() - 1;
        let last_v = self.0.controls[0].len() - 1;
        [[0, 0], [last, 0], [last, last_v], [0, last_v]]
            .map(|[u, v]| self.0.controls[u][v].map(|x| [x.lo, x.hi]))
    }
    /// Actual unit-cell endpoint derivatives crossed at00,10,11,01.
    /// These four values alone are NOT a whole-cell normal certificate.
    pub fn corner_normal_bounds(&self) -> [[[f64; 2]; 3]; 4] {
        [[0usize, 0usize], [11, 0], [11, 11], [0, 11]].map(|[u, v]| {
            let a = u.min(10);
            let b = v.min(10);
            let du: [I; 3] = std::array::from_fn(|k| {
                I::point(11.) * (self.0.controls[a + 1][v][k] - self.0.controls[a][v][k])
            });
            let dv: [I; 3] = std::array::from_fn(|k| {
                I::point(11.) * (self.0.controls[u][b + 1][k] - self.0.controls[u][b][k])
            });
            std::array::from_fn(|k| {
                let n = du[(k + 1) % 3] * dv[(k + 2) % 3] - du[(k + 2) % 3] * dv[(k + 1) % 3];
                [n.lo, n.hi]
            })
        })
    }
    pub fn bounds(&self) -> Result<[[f64; 2]; 3], Error> {
        self.0.bounds()
    }
    pub fn point_bounds(&self, uv: [[f64; 2]; 2]) -> Result<[[f64; 2]; 3], Error> {
        self.0.point_bounds(uv)
    }
    /// Derivatives use THIS leaf's unit square, with subdivision chain factors.
    /// Whole actual Pu and Pv-k*Pu coefficient hulls. Both derivative nets
    /// are elevated exactly to degree(11,11) BEFORE subtraction, retaining
    /// their shared-source dependence. For every finite k,
    /// Pu cross (Pv-k*Pu) equals the SAME source normal Pu cross Pv.
    /// The caller owns aggregate work/storage/cancellation policy.
    pub fn conditioned_derivative_bounds(&self, k: f64) -> Result<[[[f64; 2]; 3]; 2], Error> {
        let [u, r] = conditioned_controls(&self.0.controls, k)?;
        Ok([hull(u.into_iter())?, hull(r.into_iter())?])
    }
    /// Convex coefficient bound of this actual source minus its exact
    /// bilinear corner interpolant. The four corner coefficients are elevated
    /// to degree(11,11) before subtraction; no derivative or sampling estimate.
    pub fn bilinear_error_bound_mm(&self) -> Result<f64, Error> {
        let corners = [
            self.0.controls[0][0],
            self.0.controls[11][0],
            self.0.controls[11][11],
            self.0.controls[0][11],
        ];
        let weights: [I; 12] = std::array::from_fn(|i| I::point(i as f64) / I::point(11.));
        let mut error = 0f64;
        for (i, row) in self.0.controls.iter().enumerate() {
            let a = weights[i];
            let b = I::point(1.) - a;
            let top: [I; 3] = std::array::from_fn(|k| b * corners[0][k] + a * corners[1][k]);
            let bottom: [I; 3] = std::array::from_fn(|k| b * corners[3][k] + a * corners[2][k]);
            for (j, p) in row.iter().enumerate() {
                let a = weights[j];
                let b = I::point(1.) - a;
                let norm = (0..3).fold(I::point(0.), |sum, k| {
                    sum + I::point((p[k] - b * top[k] - a * bottom[k]).magnitude())
                });
                norm.checked()?;
                error = error.max(norm.hi);
            }
        }
        Ok(error)
    }
    pub fn derivative_bounds(&self) -> Result<TensorDerivatives, Error> {
        self.0.derivative_bounds()
    }
    pub fn first_derivative_bounds(&self) -> Result<TensorFirstDerivatives11, Error> {
        Ok(TensorFirstDerivatives11 {
            position: hull(self.0.controls.iter().flatten().copied())?,
            du: hull(derivative(&self.0.controls, 0).into_iter().flatten())?,
            dv: hull(derivative(&self.0.controls, 1).into_iter().flatten())?,
        })
    }
    /// Whole-interval first derivatives of this actual degree-eleven leaf.
    /// These are calculation bounds, without material or normal admission.
    pub fn first_derivative_bounds_at(
        &self,
        uv: [[f64; 2]; 2],
    ) -> Result<TensorFirstDerivatives11, Error> {
        if uv.iter().any(|d| {
            !d[0].is_finite() || !d[1].is_finite() || d[0] < 0. || d[0] > d[1] || d[1] > 1.
        }) {
            return Err(Error::Invalid);
        }
        let [position, du, dv] = tensor_first_bounds(&self.0, uv.map(|d| I::new(d[0], d[1])));
        let pair = |p: [I; 3]| p.map(|x| [x.lo, x.hi]);
        Ok(TensorFirstDerivatives11 {
            position: pair(position),
            du: pair(du),
            dv: pair(dv),
        })
    }
    /// Actual position Taylor coefficients through order eight on this leaf's
    /// unit square. Differentiate this degree-eleven net, never a D9 proxy.
    pub fn position_taylor8_bounds(&self, along_v: bool) -> [[[f64; 2]; 3]; 9] {
        let mut controls = self.0.controls.clone();
        let mut factorial = 1.;
        std::array::from_fn(|order| {
            let coefficient = coefficient_hull(&controls).map(|x| x / I::point(factorial));
            if order < 8 {
                controls = derivative(&controls, usize::from(along_v));
                factorial *= (order + 1) as f64;
            }
            coefficient.map(|x| [x.lo, x.hi])
        })
    }
    pub fn derivative_bounds_at(&self, uv: [[f64; 2]; 2]) -> Result<TensorDerivatives, Error> {
        self.0.derivative_bounds_at(uv)
    }
    pub fn subdivide(&self, axis: usize) -> Result<[Self; 2], Error> {
        self.0.subdivide(axis).map(|children| children.map(Self))
    }
}
// Complete actual coefficient nets, retained separately for the independent
// rational coefficient oracle. The public calculation operation returns hulls.
fn conditioned_controls(p: &[Vec<BoxPoint>], k: f64) -> Result<[Vec<BoxPoint>; 2], Error> {
    if !k.is_finite() {
        return Err(Error::Invalid);
    }
    let du = derivative(p, 0);
    let dv = derivative(p, 1);
    let mut elevated_u = Vec::with_capacity(144);
    let mut residual = Vec::with_capacity(144);
    for i in 0..12 {
        for j in 0..12 {
            let elevated = |p: &[Vec<BoxPoint>], axis: usize| {
                let n = if axis == 0 { i } else { j };
                let at = |n: usize| if axis == 0 { p[n][j] } else { p[i][n] };
                if n == 0 {
                    return at(0);
                }
                if n == 11 {
                    return at(10);
                }
                let weight = I::point(n as f64) / I::point(11.);
                let a = at(n - 1);
                let b = at(n);
                std::array::from_fn(|q| weight * a[q] + (I::point(1.) - weight) * b[q])
            };
            let u = elevated(&du, 0);
            let v = elevated(&dv, 1);
            elevated_u.push(u);
            residual.push(std::array::from_fn(|q| v[q] - I::point(k) * u[q]));
        }
    }
    Ok([elevated_u, residual])
}
/// Calculation-only interval coefficients in the same bounded source-mm chart;
/// this makes no source/geometry admission or manufacturing claim.
pub(crate) fn calculation_cell11(controls: [[I; 3]; 144]) -> Result<TensorCell11, Error> {
    if controls
        .iter()
        .flatten()
        .any(|x| !x.valid() || !x.lo.is_finite() || !x.hi.is_finite() || x.magnitude() > 1e7)
    {
        return Err(Error::Invalid);
    }
    Ok(TensorCell11(TensorCell {
        controls: controls
            .as_chunks::<12>()
            .0
            .iter()
            .map(|r| r.to_vec())
            .collect(),
        domain: [[0., 1.]; 2],
        depth: 0,
    }))
}
#[cfg(test)]
mod tests;
