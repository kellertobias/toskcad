//! Source-bound Taylor4 comparison to ORIGINAL cubic coefficient enclosures.
//! This is a whole-curve C0/C1 calculation, never raw incidence/face filling or
//! occupied material. Actual decimal/control ownership belongs to its caller.
use super::*;
/// Complete Series5 branch/copy/validation upper bound, paid before evaluation.
/// Series5 add96/sub128/mul512/div1024/sqrt512/sincos1280/checked64;
/// each rounded-inner+cone evaluation costs<=14000. Two half-span midpoint
/// and two whole-half evaluations cost<=56000; remaining9536 covers directed
/// domains, two Taylor4 conversions, validation, copies and checkpoints.
pub const TRIM_SPAN_SOURCE_WORK: usize = 65_536;
/// Exact original-cubic bisection and both complete Taylor/affine comparisons.
pub const TRIM_SPAN_COMPARE_WORK: usize = 4_096;
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
    pieces: [SpanPiece; 2],
}
#[derive(Clone, Copy)]
struct SpanPiece {
    polynomial: [[I; 3]; 5],
    fourth_variation: [I; 3],
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
        let halves = bisect(literal);
        let mut position = 0_f64;
        let mut derivative = 0_f64;
        for (piece, controls) in self.pieces.into_iter().zip(halves) {
            let bounds = difference_bounds(
                piece.polynomial,
                piece.fourth_variation,
                piece.h4,
                controls,
                origin,
                columns,
            )?;
            position = position.max(bounds[0]);
            // Return ORIGINAL span normalized derivative, not half-span units.
            derivative = derivative.max((I::point(2.) * I::point(bounds[1])).hi);
        }
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
fn bisect(p: [[I; 3]; 4]) -> [[[I; 3]; 4]; 2] {
    let average = |a: [I; 3], b: [I; 3]| std::array::from_fn(|i| (a[i] + b[i]) / I::point(2.));
    let row = [
        average(p[0], p[1]),
        average(p[1], p[2]),
        average(p[2], p[3]),
    ];
    let next = [average(row[0], row[1]), average(row[1], row[2])];
    let middle = average(next[0], next[1]);
    [
        [p[0], row[0], next[0], middle],
        [middle, next[1], row[2], p[3]],
    ]
}
// Pure numerical enclosure, separated so independent polynomial oracles do
// not fabricate a source-bound cache/receipt or a qualified geometric source.
fn difference_bounds(
    polynomial: [[I; 3]; 5],
    fourth_variation: [I; 3],
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
    // Exact degree elevation of the ORIGINAL cubic, no refit or new samples.
    let elevated = std::array::from_fn::<_, 5, _>(|i| {
        std::array::from_fn::<_, 3, _>(|j| match i {
            0 => literal[0][j],
            4 => literal[3][j],
            _ => {
                I::point(i as f64 / 4.) * literal[i - 1][j]
                    + I::point((4 - i) as f64 / 4.) * literal[i][j]
            }
        })
    });
    let delta: [[I; 3]; 5] =
        std::array::from_fn(|i| std::array::from_fn(|j| difference[i][j] - elevated[i][j]));
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
    let derivative_controls = std::array::from_fn::<_, 4, _>(|i| {
        std::array::from_fn(|j| I::point(4.) * (delta[i + 1][j] - delta[i][j]))
    });
    let mut derivative = 0_f64;
    // Restrict the DERIVATIVE polynomial, not the source position curve.
    // Its coefficient units remain d/ds; restriction introduces no factor4.
    // Three cubic bisections plus16 fixed3D norms add <512 scalar actions,
    // within the original2048 allowance for each complete comparison.
    for half in bisect(derivative_controls) {
        for quarter in bisect(half) {
            for control in quarter {
                derivative = derivative.max(valid_norm(control)?.hi);
            }
        }
    }
    // Taylor3 for f and Taylor2 for f' have fourth derivative remainders.
    // We explicitly include the midpoint quartic term in the polynomial;
    // subtracting it leaves (c4(xi)-c4(midpoint)). This uses a whole-domain
    // enclosure, never differentiation of a Lagrange xi or an assumed C5.
    // In normalized span s, factors are h^4/16 and h^4/2 respectively.
    let remainder = valid_norm(vector(fourth_variation))?;
    let position = I::point(position) + remainder * h4 / I::point(16.);
    let derivative = I::point(derivative) + remainder * h4 / I::point(2.);
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
// Convert midpoint Taylor4 to degree4 Bernstein. Series already includes
// factorials; h^k applies each original source-to-span chain factor once.
fn polynomial(midpoint: [Series<5>; 3], h: I) -> [[I; 3]; 5] {
    let h2 = h * h;
    let powers = [I::point(1.), h, h2, h2 * h, h2 * h2];
    std::array::from_fn(|i| {
        std::array::from_fn(|axis| {
            let q = std::array::from_fn::<_, 5, _>(|k| midpoint[axis].c[k] * powers[k]);
            let a = [
                q[0] - q[1] / I::point(2.) + q[2] / I::point(4.) - q[3] / I::point(8.)
                    + q[4] / I::point(16.),
                q[1] - q[2] + I::point(0.75) * q[3] - q[4] / I::point(2.),
                q[2] - I::point(1.5) * q[3] + I::point(1.5) * q[4],
                q[3] - I::point(2.) * q[4],
                q[4],
            ];
            match i {
                0 => a[0],
                1 => a[0] + a[1] / I::point(4.),
                2 => a[0] + a[1] / I::point(2.) + a[2] / I::point(6.),
                3 => a[0] + I::point(0.75) * a[1] + a[2] / I::point(2.) + a[3] / I::point(4.),
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
        let half = h / I::point(2.);
        let h4 = half * half * half * half;
        let mut pieces = [SpanPiece {
            polynomial: [[I::point(0.); 3]; 5],
            fourth_variation: [I::point(0.); 3],
            h4,
        }; 2];
        for (piece, [left, right]) in pieces.iter_mut().zip([[a, m], [m, b]]) {
            let center = (left + right) / I::point(2.);
            let midpoint = evaluate::trim_state::<Series<5>>(self, index, center)?;
            checkpoint(interrupt)?;
            let whole = evaluate::trim_state::<Series<5>>(self, index, I::new(left.lo, right.hi))?;
            let polynomial = polynomial(midpoint, half);
            let fourth_variation = std::array::from_fn(|i| whole[i].c[4] - midpoint[i].c[4]);
            if polynomial
                .iter()
                .flatten()
                .chain(&fourth_variation)
                .chain([&h4])
                .any(|x| !x.valid())
            {
                return Err(Error::Unresolved);
            }
            *piece = SpanPiece {
                polynomial,
                fourth_variation,
                h4,
            };
        }
        checkpoint(interrupt)?;
        Ok(TrimSpanCorrespondence {
            events: self,
            index,
            parameters,
            pieces,
        })
    }
}

#[cfg(test)]
mod tests;
