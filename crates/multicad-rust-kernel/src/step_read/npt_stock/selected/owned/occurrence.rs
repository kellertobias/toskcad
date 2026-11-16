//! A distinct exact Native occurrence chain; original STEP geometry is unchanged.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct LiteralTransform {
    matrix: [[f64; 3]; 3],
    translation: [f64; 3],
}
impl LiteralTransform {
    pub fn matrix(self) -> [[f64; 3]; 3] {
        self.matrix
    }
    pub fn translation_mm(self) -> [f64; 3] {
        self.translation
    }
}
#[derive(Clone)]
pub struct LiteralOccurrencePlacement {
    chain: Arc<[LiteralTransform]>,
    linear: [[I; 3]; 3],
    origin: [I; 3],
    cofactor: [[I; 3]; 3],
    determinant: I,
    norm1: f64,
    delta: f64,
    condition: I,
    work: usize,
    retained: usize,
    peak: usize,
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |sum, (a, b)| sum + a * b)
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn magnitude(value: I) -> I {
    I::new(0., value.magnitude())
}
fn metric(matrix: [[I; 3]; 3]) -> I {
    let columns: [[I; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| matrix[j][i]));
    let mut delta = I::point(0.);
    for i in 0..3 {
        for j in 0..3 {
            delta = delta
                + magnitude(dot(columns[i], columns[j]) - I::point(if i == j { 1. } else { 0. }));
        }
    }
    delta
}
fn checked(input: [[f64; 2]; 3]) -> Result<[I; 3]> {
    let value = input.map(interval);
    if value.iter().any(|v| !v.valid()) {
        return Err(Error::Invalid);
    }
    Ok(value)
}
fn finite(value: [I; 3]) -> Result<[[f64; 2]; 3]> {
    if value.iter().any(|v| !v.valid()) {
        return Err(Error::Unsupported);
    }
    Ok(value.map(bounds))
}
impl LiteralOccurrencePlacement {
    pub fn new(
        matrix: [[f64; 3]; 3],
        translation_mm: [f64; 3],
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        Self::append(None, matrix, translation_mm, options, interrupt)
    }
    pub fn transformed(
        &self,
        matrix: [[f64; 3]; 3],
        translation_mm: [f64; 3],
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        Self::append(Some(self), matrix, translation_mm, options, interrupt)
    }
    fn append(
        previous: Option<&Self>,
        matrix: [[f64; 3]; 3],
        translation: [f64; 3],
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        let mut budget = Budget::new(options, interrupt)?;
        budget.step(8192)?; // Both fixed matrix products, metrics, cofactors and publication.
        if matrix
            .iter()
            .flatten()
            .chain(&translation)
            .any(|v| !v.is_finite())
        {
            return Err(Error::Invalid);
        }
        if matrix.iter().flatten().any(|v| v.abs() > 2.)
            || translation.iter().any(|v| v.abs() > 1e7)
        {
            return Err(Error::Limit);
        }
        let outer = matrix.map(|r| r.map(I::point));
        if metric(outer).hi > 1e-6 {
            return Err(Error::Unsupported);
        }
        let identity = std::array::from_fn(|i| {
            std::array::from_fn(|j| I::point(if i == j { 1. } else { 0. }))
        });
        let inner = previous.map_or(identity, |v| v.linear);
        let offset = previous.map_or([I::point(0.); 3], |v| v.origin);
        let linear = std::array::from_fn(|i| {
            std::array::from_fn(|j| dot(outer[i], std::array::from_fn(|k| inner[k][j])))
        });
        let origin = std::array::from_fn(|i| dot(outer[i], offset) + I::point(translation[i]));
        let delta = metric(linear).hi;
        let cofactor = [
            cross(linear[1], linear[2]),
            cross(linear[2], linear[0]),
            cross(linear[0], linear[1]),
        ];
        let determinant = dot(linear[0], cofactor[0]);
        if !delta.is_finite()
            || delta > 1e-6
            || !determinant.valid()
            || determinant.lo <= 0.
            || origin.iter().any(|v| !v.valid())
        {
            return Err(Error::Unsupported);
        }
        let norm1 = (0..3)
            .map(|j| {
                (0..3)
                    .fold(I::point(0.), |s, i| s + magnitude(linear[i][j]))
                    .hi
            })
            .fold(0f64, f64::max);
        let old = previous.map_or(&[][..], |p| p.chain.as_ref());
        hold(&mut budget, previous.map_or(0, |p| p.retained))?;
        let mut chain = budget.allocate(old.len().checked_add(1).ok_or(Error::Limit)?)?;
        budget.step(old.len().checked_mul(12).ok_or(Error::Limit)?)?;
        chain.extend_from_slice(old);
        chain.push(LiteralTransform {
            matrix,
            translation,
        });
        let retained = size_of::<Self>()
            + 2 * size_of::<usize>()
            + chain.len() * size_of::<LiteralTransform>();
        hold(&mut budget, retained)?; // Arc conversion can coexist with old/new Vec storage.
        let chain = Arc::from(chain);
        let peak = budget.peak;
        Ok(Self {
            chain,
            linear,
            origin,
            cofactor,
            determinant,
            norm1,
            delta,
            condition: (I::point(1.) + I::point(delta)) / (I::point(1.) - I::point(delta)),
            work: budget.used,
            retained,
            peak,
        })
    }
    pub fn transforms(&self) -> &[LiteralTransform] {
        &self.chain
    }
    pub fn linear_rows_bounds(&self) -> [[[f64; 2]; 3]; 3] {
        self.linear.map(|r| r.map(bounds))
    }
    pub fn translation_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.origin.map(bounds)
    }
    pub fn determinant_bounds(&self) -> [f64; 2] {
        bounds(self.determinant)
    }
    pub fn metric_error(&self) -> f64 {
        self.delta
    }
    /// Spectral norm bound from the complete retained-chain Gram enclosure.
    pub fn linear_l2_norm_upper(&self) -> f64 {
        (I::point(1.) + I::point(self.delta)).sqrt().hi
    }
    pub fn linear_l1_norm(&self) -> f64 {
        self.norm1
    }
    pub fn area_scale_bounds(&self) -> [f64; 2] {
        [(1. - self.delta).next_down(), (1. + self.delta).next_up()]
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
    pub fn normal_condition_upper(&self) -> f64 {
        self.condition.hi
    }
    pub fn normal_bounds(&self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let input = checked(input)?;
        let normal = self.cofactor.map(|r| dot(r, input));
        let length = normal
            .into_iter()
            .fold(I::point(0.), |sum, v| sum + v.square())
            .sqrt();
        if !length.valid() || length.lo <= 0. {
            return Err(Error::Unsupported);
        }
        finite(normal.map(|v| v / length))
    }
    pub fn vector_bounds(&self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let input = checked(input)?;
        finite(self.linear.map(|r| dot(r, input)))
    }
    pub fn map_box_mm(&self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let input = checked(input)?;
        let result = std::array::from_fn(|i| dot(self.linear[i], input) + self.origin[i]);
        let result = finite(result)?;
        if result.iter().flatten().any(|v| v.abs() > 1e7) {
            return Err(Error::Limit);
        }
        Ok(result)
    }
    pub fn inverse_box_mm(&self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let input = checked(input)?;
        let shifted = std::array::from_fn(|i| input[i] - self.origin[i]);
        finite(std::array::from_fn(|i| {
            dot(std::array::from_fn(|j| self.cofactor[j][i]), shifted) / self.determinant
        }))
    }
}

#[cfg(test)]
mod tests;
