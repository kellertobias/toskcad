//! Fixed degree11 numerical cells, distinct from the degree9 source machinery.
//! No radial chart, occupied boundary, material, normal or mass is certified.
use super::*;
const DEGREE: usize = 11;
#[derive(Clone, Debug, PartialEq)]
pub struct RetainedTensor11 {
    spans: [usize; 2],
    poles: Arc<[Point]>,
}
/// Degree11 only. Its private wrapper cannot enter degree9 normal helpers.
#[derive(Clone, Debug)]
pub struct TensorCell11(TensorCell);
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
    pub fn bounds(&self) -> Result<[[f64; 2]; 3], Error> {
        self.0.bounds()
    }
    pub fn point_bounds(&self, uv: [[f64; 2]; 2]) -> Result<[[f64; 2]; 3], Error> {
        self.0.point_bounds(uv)
    }
    /// Derivatives use THIS leaf's unit square, with subdivision chain factors.
    pub fn derivative_bounds(&self) -> Result<TensorDerivatives, Error> {
        self.0.derivative_bounds()
    }
    pub fn derivative_bounds_at(&self, uv: [[f64; 2]; 2]) -> Result<TensorDerivatives, Error> {
        self.0.derivative_bounds_at(uv)
    }
    pub fn subdivide(&self, axis: usize) -> Result<[Self; 2], Error> {
        self.0.subdivide(axis).map(|children| children.map(Self))
    }
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
