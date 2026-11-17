//! Source-bound Taylor3 comparison to ORIGINAL cubic coefficient enclosures.
//! This is a whole-curve C0/C1 calculation, never raw incidence/face filling or
//! occupied material. Actual decimal/control ownership belongs to its caller.
use super::*;
/// Complete Series5 branch/copy/validation upper bound, paid before evaluation.
/// Series5 add96/sub128/mul512/div1024/sqrt512/sincos1280/checked64;
/// the rounded-inner+cone branch costs<=14000, evaluated twice. The remaining
///4768 covers Taylor3 conversion, input/source validation, all span/cache
/// copies/control and checkpoints (the source event graph is borrowed).
pub const TRIM_SPAN_SOURCE_WORK: usize = 32_768;
/// Actual affine/control difference, all Bernstein hulls and publication.
pub const TRIM_SPAN_COMPARE_WORK: usize = 2_048;
const MAX_WORK: usize = 16_000_000;
const MAX_BYTES: usize = 64 * 1024 * 1024;
const SOURCE_SCRATCH: usize = 32 * 1024;
const COMPARE_SCRATCH: usize = 8 * 1024;
/// Cached jets of this SAME actual source/index/domain. No raw curve, geometry
/// recipe or sampled fit defines the source, and callers cannot inject jets.
pub struct TrimSpanCorrespondence<'events, 'source> {
    events: &'events BoundaryEvents<'source>,
    index: usize,
    parameters: [[f64; 2]; 2],
    polynomial: [[I; 3]; 4],
    fourth: [I; 3],
    h4: I,
}
pub struct TrimSpanComparison<'query, 'events, 'source> {
    source: &'query TrimSpanCorrespondence<'events, 'source>,
    controls: &'query [[[f64; 2]; 3]; 4],
    origin: [f64; 3],
    columns: [[f64; 3]; 3],
    position: f64,
    derivative: f64,
    peak: usize,
}
impl<'events, 'source> TrimSpanCorrespondence<'events, 'source> {
    pub fn events(&self) -> &'events BoundaryEvents<'source> {
        self.events
    }
    pub fn trim_index(&self) -> usize {
        self.index
    }
    pub fn source_parameter_bounds(&self) -> [[f64; 2]; 2] {
        self.parameters
    }
    pub fn construction_work_used(&self) -> usize {
        TRIM_SPAN_SOURCE_WORK
    }
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
    }
    pub fn construction_peak_bytes(&self) -> usize {
        self.retained_budget_bytes() + SOURCE_SCRATCH
    }
    /// Original controls must enclose the ACTUAL decimal source coefficients.
    /// The frame is an exact defining binary64 affine map, not a rounded
    /// approximation silently substituted for an actual STEP operator.
    pub fn compare_controls<'q>(
        &'q self,
        controls: &'q [[[f64; 2]; 3]; 4],
        origin: [f64; 3],
        columns: [[f64; 3]; 3],
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<TrimSpanComparison<'q, 'events, 'source>, Error> {
        checkpoint(interrupt)?;
        let peak = self.retained_budget_bytes()
            + std::mem::size_of::<TrimSpanComparison<'_, '_, '_>>()
            + std::mem::size_of_val(controls)
            + COMPARE_SCRATCH;
        admit(maximum_work, maximum_bytes, TRIM_SPAN_COMPARE_WORK, peak)?;
        let literal = controls.map(|p| p.map(|[lo, hi]| I::new(lo, hi)));
        if literal.iter().flatten().any(|x| !x.valid())
            || origin
                .into_iter()
                .chain(columns.into_iter().flatten())
                .any(|x| !x.is_finite())
        {
            return Err(Error::Invalid);
        }
        if literal.iter().flatten().any(|x| x.magnitude() > 1e12)
            || origin.into_iter().any(|x| x.abs() > 1e12)
            || columns.into_iter().flatten().any(|x| x.abs() > 1e9)
        {
            return Err(Error::Limit);
        }
        let [position, derivative] = difference_bounds(
            self.polynomial,
            self.fourth,
            self.h4,
            literal,
            origin,
            columns,
        )?;
        checkpoint(interrupt)?;
        Ok(TrimSpanComparison {
            source: self,
            controls,
            origin,
            columns,
            position,
            derivative,
            peak,
        })
    }
}
impl<'q, 'events, 'source> TrimSpanComparison<'q, 'events, 'source> {
    pub fn source(&self) -> &'q TrimSpanCorrespondence<'events, 'source> {
        self.source
    }
    pub fn original_control_bounds(&self) -> &'q [[[f64; 2]; 3]; 4] {
        self.controls
    }
    /// The exact defining source affine used in this numerical comparison.
    pub fn source_origin_mm(&self) -> [f64; 3] {
        self.origin
    }
    pub fn source_columns(&self) -> [[f64; 3]; 3] {
        self.columns
    }
    pub fn maximum_position_mm(&self) -> f64 {
        self.position
    }
    /// Derivative in the ORIGINAL cubic span's normalized unit parameter.
    /// Source trim factors and the actual endpoint-span h are included once.
    pub fn maximum_derivative_mm(&self) -> f64 {
        self.derivative
    }
    pub fn new_work_used(&self) -> usize {
        TRIM_SPAN_COMPARE_WORK
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
// Pure numerical enclosure, separated so independent polynomial oracles do
// not fabricate a source-bound cache/receipt or a qualified geometric source.
fn difference_bounds(
    polynomial: [[I; 3]; 4],
    fourth: [I; 3],
    h4: I,
    literal: [[I; 3]; 4],
    origin: [f64; 3],
    columns: [[f64; 3]; 3],
) -> Result<[f64; 2], Error> {
    let vector = |p: [I; 3]| {
        std::array::from_fn::<_, 3, _>(|i| {
            (0..3).fold(I::point(0.), |s, j| s + I::point(columns[j][i]) * p[j])
        })
    };
    let difference = polynomial.map(|p| {
        let v = vector(p);
        std::array::from_fn::<_, 3, _>(|i| I::point(origin[i]) + v[i])
    });
    let delta: [[I; 3]; 4] =
        std::array::from_fn(|i| std::array::from_fn(|j| difference[i][j] - literal[i][j]));
    let valid_norm = |value: [I; 3]| -> Result<I, Error> {
        if value.iter().any(|x| !x.valid()) {
            return Err(Error::Unresolved);
        }
        let n = norm(value);
        if !n.valid() {
            return Err(Error::Unresolved);
        }
        Ok(n)
    };
    let mut position = 0_f64;
    for value in delta {
        position = position.max(valid_norm(value)?.hi);
    }
    let mut derivative = 0_f64;
    for i in 0..3 {
        derivative = derivative.max(
            valid_norm(std::array::from_fn(|j| {
                I::point(3.) * (delta[i + 1][j] - delta[i][j])
            }))?
            .hi,
        );
    }
    // Taylor of the derivative, NOT differentiation of a Lagrange ξ:
    // h*f4*(h/2)^3/6 = f4*h4/48 = Series_c4*h4/2.
    // Each coordinate's remainder may use a different intermediate point;
    // whole-vector coordinate boxes enclose all of those choices.
    let fourth = valid_norm(vector(fourth))?;
    let position = I::point(position) + fourth * h4 / I::point(16.);
    let derivative = I::point(derivative) + fourth * h4 / I::point(2.);
    if !position.valid() || !derivative.valid() {
        return Err(Error::Unresolved);
    }
    Ok([position.hi, derivative.hi])
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn admit(work: usize, bytes: usize, required: usize, peak: usize) -> Result<(), Error> {
    if work > MAX_WORK || bytes > MAX_BYTES {
        return Err(Error::Invalid);
    }
    if work < required || bytes < peak {
        return Err(Error::Limit);
    }
    Ok(())
}
// Convert Σ q_k (s−1/2)^k to degree3 Bernstein coefficients. All factorials
// are ALREADY included in Series.c[k]; multiplying h^k gives the sole original
// source-normalized→actual-span chain factor. No extra span width is applied.
fn polynomial(midpoint: [Series<5>; 3], h: I) -> [[I; 3]; 4] {
    let powers = [I::point(1.), h, h * h, h * h * h];
    std::array::from_fn(|i| {
        std::array::from_fn(|axis| {
            let q = std::array::from_fn::<_, 4, _>(|k| midpoint[axis].c[k] * powers[k]);
            let a = [
                q[0] - q[1] / I::point(2.) + q[2] / I::point(4.) - q[3] / I::point(8.),
                q[1] - q[2] + I::point(0.75) * q[3],
                q[2] - I::point(1.5) * q[3],
                q[3],
            ];
            match i {
                0 => a[0],
                1 => a[0] + a[1] / I::point(3.),
                2 => a[0] + I::point(2.) * a[1] / I::point(3.) + a[2] / I::point(3.),
                _ => a.into_iter().fold(I::point(0.), |s, x| s + x),
            }
        })
    })
}
impl<'source> BoundaryEvents<'source> {
    /// Prepare a genuine whole-span jet cache. Endpoint intervals enclose the
    /// actual source-normalized knots; overlap/singular C4 charts refuse.
    /// The borrowed source/event graph is a separate caller-owned byte scope.
    pub fn prepare_trim_span_correspondence<'events>(
        &'events self,
        index: usize,
        parameters: [[f64; 2]; 2],
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<TrimSpanCorrespondence<'events, 'source>, Error> {
        checkpoint(interrupt)?;
        let peak = std::mem::size_of::<TrimSpanCorrespondence<'_, '_>>() + SOURCE_SCRATCH;
        admit(maximum_work, maximum_bytes, TRIM_SPAN_SOURCE_WORK, peak)?;
        let [a, b] = parameters.map(|[lo, hi]| I::new(lo, hi));
        if index >= 4 || !a.valid() || !b.valid() || a.lo < 0. || b.hi > 1. || a.lo >= b.hi {
            return Err(Error::Invalid);
        }
        let h = b - a;
        if a.hi >= b.lo || h.lo <= 0. {
            return Err(Error::Unresolved);
        }
        let m = (a + b) / I::point(2.);
        let midpoint = evaluate::trim_state::<Series<5>>(self, index, m)?;
        checkpoint(interrupt)?;
        let whole = evaluate::trim_state::<Series<5>>(self, index, I::new(a.lo, b.hi))?;
        let polynomial = polynomial(midpoint, h);
        let fourth = whole.map(|s| s.c[4]);
        let h4 = h * h * h * h;
        if polynomial
            .iter()
            .flatten()
            .chain(&fourth)
            .chain([&h4])
            .any(|x| !x.valid())
        {
            return Err(Error::Unresolved);
        }
        checkpoint(interrupt)?;
        Ok(TrimSpanCorrespondence {
            events: self,
            index,
            parameters,
            polynomial,
            fourth,
            h4,
        })
    }
}

#[cfg(test)]
mod tests;
