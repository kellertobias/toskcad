//! Rational B-spline curves and surfaces with first derivatives.
use super::{Error, Point, Result, ensure, scale, sub};

/// The largest admitted control net of one curve or surface.
pub const MAX_CONTROL_POINTS: usize = 65_536;
pub const MAX_DEGREE: usize = 32;

type Homogeneous = [f64; 4];

fn knots(
    distinct: &[f64],
    multiplicities: &[usize],
    degree: usize,
    count: usize,
) -> Result<Vec<f64>> {
    ensure(
        (1..=MAX_DEGREE).contains(&degree) && distinct.len() == multiplicities.len(),
        Error::Invalid,
    )?;
    ensure(count > degree && count <= MAX_CONTROL_POINTS, Error::Limit)?;
    let mut expanded = Vec::with_capacity(count + degree + 1);
    for (&knot, &multiplicity) in distinct.iter().zip(multiplicities) {
        ensure(
            knot.is_finite()
                && (1..=degree + 1).contains(&multiplicity)
                && expanded.last().is_none_or(|&last| last < knot),
            Error::Invalid,
        )?;
        expanded.extend(std::iter::repeat_n(knot, multiplicity));
    }
    ensure(
        expanded.len() == count + degree + 1 && expanded[degree] < expanded[count],
        Error::Invalid,
    )?;
    Ok(expanded)
}

fn homogeneous(points: &[Point], weights: Option<&[f64]>) -> Result<Vec<Homogeneous>> {
    let unit = vec![1.; points.len()];
    let weights = weights.unwrap_or(&unit);
    ensure(weights.len() == points.len(), Error::Invalid)?;
    points
        .iter()
        .zip(weights)
        .map(|(p, &w)| {
            ensure(
                w.is_finite() && w > 0. && p.iter().all(|x| x.is_finite()),
                Error::Invalid,
            )?;
            Ok([p[0] * w, p[1] * w, p[2] * w, w])
        })
        .collect()
}

/// Room for the nonzero basis functions of the largest admitted degree.
const BASIS: usize = MAX_DEGREE + 1;

/// Span index, basis values and basis derivatives at `t`; the first
/// `degree + 1` entries of each array are the nonzero ones.
///
/// The Cox-de Boor triangle is built once on the stack. Its row of
/// `degree - 1` gives the lower-degree values the derivatives need, and
/// every row is formed by the same operations in the same order as a
/// separate evaluation of that degree, so both are bit-identical to it.
fn basis(knots: &[f64], degree: usize, t: f64) -> (usize, [f64; BASIS], [f64; BASIS]) {
    let count = knots.len() - degree - 1;
    let t = t.clamp(knots[degree], knots[count]);
    // The last span in `degree..count` whose knot does not exceed `t`, as
    // a scan up from `degree` finds it; knots never decrease.
    let above = degree + 1 + knots[degree + 1..count].partition_point(|&knot| t >= knot);
    let span = above - 1;
    let mut value = [0.; BASIS];
    value[0] = 1.;
    let mut lower = value;
    let mut next = [0.; BASIS];
    for j in 1..=degree {
        if j == degree {
            lower = value;
        }
        next[..=j].fill(0.);
        for r in 0..j {
            let i = span + 1 + r - j;
            let w = value[r] / (knots[i + j] - knots[i]);
            next[r + 1] += (t - knots[i]) * w;
            next[r] += (knots[i + j] - t) * w;
        }
        value[..=j].copy_from_slice(&next[..=j]);
    }
    let mut derivative = [0.; BASIS];
    for (m, slot) in derivative[..=degree].iter_mut().enumerate() {
        let i = span + m - degree;
        if m >= 1 {
            *slot += lower[m - 1] / (knots[i + degree] - knots[i]);
        }
        if m < degree {
            *slot -= lower[m] / (knots[i + degree + 1] - knots[i + 1]);
        }
        *slot *= degree as f64;
    }
    (span, value, derivative)
}

fn project(a: Homogeneous, derivatives: [Homogeneous; 2]) -> (Point, [Point; 2]) {
    let point = [a[0] / a[3], a[1] / a[3], a[2] / a[3]];
    let derivative =
        derivatives.map(|d| scale(sub([d[0], d[1], d[2]], scale(point, d[3])), 1. / a[3]));
    (point, derivative)
}

fn accumulate(into: &mut Homogeneous, control: Homogeneous, weight: f64) {
    for (slot, value) in into.iter_mut().zip(control) {
        *slot += value * weight;
    }
}

/// Distinct knot values strictly inside `(a, b)`, each span cut `degree` times.
fn breaks(knots: &[f64], degree: usize, period: Option<f64>, a: f64, b: f64) -> Vec<f64> {
    let (low, high) = (knots[degree], knots[knots.len() - degree - 1]);
    let mut distinct: Vec<f64> = knots[degree..knots.len() - degree].to_vec();
    distinct.dedup();
    let mut cuts = vec![a];
    let turns = period.map_or(0..=0, |p| {
        (((a - high) / p).floor() as i64)..=(((b - low) / p).ceil() as i64)
    });
    for turn in turns {
        let shift = turn as f64 * period.unwrap_or(0.);
        for pair in distinct.windows(2) {
            for piece in 0..degree {
                let cut = shift + pair[0] + (pair[1] - pair[0]) * piece as f64 / degree as f64;
                if cut > *cuts.last().unwrap_or(&a) + 1e-12 && cut < b - 1e-12 {
                    cuts.push(cut);
                }
            }
        }
    }
    cuts.push(b);
    cuts
}

#[derive(Clone, Debug, PartialEq)]
pub struct NurbsCurve {
    degree: usize,
    knots: Vec<f64>,
    control: Vec<Homogeneous>,
    distinct: Vec<f64>,
    multiplicities: Vec<usize>,
    points: Vec<Point>,
    weights: Option<Vec<f64>>,
}

impl NurbsCurve {
    pub fn new(
        degree: usize,
        points: Vec<Point>,
        weights: Option<Vec<f64>>,
        distinct: Vec<f64>,
        multiplicities: Vec<usize>,
    ) -> Result<Self> {
        let knots = knots(&distinct, &multiplicities, degree, points.len())?;
        let control = homogeneous(&points, weights.as_deref())?;
        Ok(Self {
            degree,
            knots,
            control,
            distinct,
            multiplicities,
            points,
            weights,
        })
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn points(&self) -> &[Point] {
        &self.points
    }

    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    pub fn knots(&self) -> (&[f64], &[usize]) {
        (&self.distinct, &self.multiplicities)
    }

    pub fn range(&self) -> [f64; 2] {
        [self.knots[self.degree], self.knots[self.control.len()]]
    }

    pub fn evaluate(&self, t: f64) -> (Point, Point) {
        let (span, value, derivative) = basis(&self.knots, self.degree, t);
        let mut sums = [[0.; 4]; 2];
        for k in 0..=self.degree {
            let control = self.control[span + k - self.degree];
            accumulate(&mut sums[0], control, value[k]);
            accumulate(&mut sums[1], control, derivative[k]);
        }
        let (point, [derivative, _]) = project(sums[0], [sums[1], sums[1]]);
        (point, derivative)
    }

    pub(crate) fn breaks(&self, period: Option<f64>, a: f64, b: f64) -> Vec<f64> {
        breaks(&self.knots, self.degree, period, a, b)
    }

    pub(crate) fn mapped(&self, map: &dyn Fn(Point) -> Point) -> Self {
        let points: Vec<_> = self.points.iter().map(|&p| map(p)).collect();
        let control = homogeneous(&points, self.weights.as_deref()).unwrap_or_default();
        Self {
            points,
            control,
            ..self.clone()
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct NurbsSurface {
    degrees: [usize; 2],
    knots: [Vec<f64>; 2],
    counts: [usize; 2],
    control: Vec<Homogeneous>,
    distinct: [Vec<f64>; 2],
    multiplicities: [Vec<usize>; 2],
    points: Vec<Point>,
    weights: Option<Vec<f64>>,
}

impl NurbsSurface {
    /// `points` is row-major: `counts[0]` rows along u, `counts[1]` along v.
    pub fn new(
        degrees: [usize; 2],
        counts: [usize; 2],
        points: Vec<Point>,
        weights: Option<Vec<f64>>,
        distinct: [Vec<f64>; 2],
        multiplicities: [Vec<usize>; 2],
    ) -> Result<Self> {
        ensure(
            counts[0].checked_mul(counts[1]) == Some(points.len()),
            Error::Invalid,
        )?;
        ensure(points.len() <= MAX_CONTROL_POINTS, Error::Limit)?;
        let knots = [
            knots(&distinct[0], &multiplicities[0], degrees[0], counts[0])?,
            knots(&distinct[1], &multiplicities[1], degrees[1], counts[1])?,
        ];
        let control = homogeneous(&points, weights.as_deref())?;
        Ok(Self {
            degrees,
            knots,
            counts,
            control,
            distinct,
            multiplicities,
            points,
            weights,
        })
    }

    pub fn degrees(&self) -> [usize; 2] {
        self.degrees
    }

    pub fn counts(&self) -> [usize; 2] {
        self.counts
    }

    pub fn points(&self) -> &[Point] {
        &self.points
    }

    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    pub fn knots(&self, axis: usize) -> (&[f64], &[usize]) {
        (&self.distinct[axis], &self.multiplicities[axis])
    }

    pub fn range(&self, axis: usize) -> [f64; 2] {
        [
            self.knots[axis][self.degrees[axis]],
            self.knots[axis][self.counts[axis]],
        ]
    }

    pub fn evaluate(&self, u: f64, v: f64) -> (Point, [Point; 2]) {
        let (su, nu, du) = basis(&self.knots[0], self.degrees[0], u);
        let (sv, nv, dv) = basis(&self.knots[1], self.degrees[1], v);
        let mut sums = [[0.; 4]; 3];
        for i in 0..=self.degrees[0] {
            let row = (su + i - self.degrees[0]) * self.counts[1];
            for j in 0..=self.degrees[1] {
                let control = self.control[row + sv + j - self.degrees[1]];
                accumulate(&mut sums[0], control, nu[i] * nv[j]);
                accumulate(&mut sums[1], control, du[i] * nv[j]);
                accumulate(&mut sums[2], control, nu[i] * dv[j]);
            }
        }
        project(sums[0], [sums[1], sums[2]])
    }

    pub(crate) fn breaks(&self, axis: usize, period: Option<f64>, a: f64, b: f64) -> Vec<f64> {
        breaks(&self.knots[axis], self.degrees[axis], period, a, b)
    }

    /// An upper bound of the distance between the two boundary curves at
    /// the ends of `axis`, over the WHOLE other parameter range.
    ///
    /// Each boundary is an exact rational curve in the other parameter whose
    /// homogeneous controls are the end-basis combinations of the net rows.
    /// Both share one basis, so with positive weights each boundary point is
    /// a convex combination of its controls. For weight ratios
    /// `r_j = w_b,j / w_a,j` in `[r_min, r_max]` the two convex weightings
    /// differ by at most `(r_max - r_min) / r_min` in total, giving
    /// `max_j |a_j - b_j| + (r_max - r_min) / r_min * max_j |b_j - b_0|`.
    /// A floating-point allowance covers the end-basis and division
    /// roundings. `None` means no bound was established (not closed).
    pub(crate) fn boundary_gap(&self, axis: usize) -> Option<f64> {
        let other = 1 - axis;
        let [low, high] = self.range(axis);
        let ends = [low, high].map(|t| basis(&self.knots[axis], self.degrees[axis], t));
        let count = self.counts[other];
        let control = |i: usize, j: usize| {
            if axis == 0 {
                self.control[i * self.counts[1] + j]
            } else {
                self.control[j * self.counts[1] + i]
            }
        };
        let row = |end: usize, j: usize| -> Homogeneous {
            let (span, values, _) = &ends[end];
            let mut sum = [0.; 4];
            for (k, &value) in values[..=self.degrees[axis]].iter().enumerate() {
                accumulate(&mut sum, control(span + k - self.degrees[axis], j), value);
            }
            sum
        };
        let (mut gap, mut spread, mut ratio, mut scale) =
            (0_f64, 0_f64, [f64::INFINITY, 0_f64], 0_f64);
        let mut first = None;
        for j in 0..count {
            let [a, b] = [row(0, j), row(1, j)];
            if !(a[3] > 0. && b[3] > 0. && a.iter().chain(&b).all(|x| x.is_finite())) {
                return None;
            }
            let pa = [a[0] / a[3], a[1] / a[3], a[2] / a[3]];
            let pb = [b[0] / b[3], b[1] / b[3], b[2] / b[3]];
            gap = gap.max(super::distance(pa, pb));
            let origin = *first.get_or_insert(pb);
            spread = spread.max(super::distance(pb, origin));
            let r = b[3] / a[3];
            ratio = [ratio[0].min(r), ratio[1].max(r)];
            scale = pa.iter().chain(&pb).fold(scale, |m, x| m.max(x.abs()));
        }
        let weights = if ratio[1] == ratio[0] {
            0.
        } else {
            (ratio[1] - ratio[0]) / ratio[0]
        };
        let rounding = 64. * f64::EPSILON * (self.degrees[axis] + 2) as f64 * scale;
        let bound = gap + weights * spread + rounding;
        bound.is_finite().then_some(bound)
    }

    pub(crate) fn mapped(&self, map: &dyn Fn(Point) -> Point) -> Self {
        let points: Vec<_> = self.points.iter().map(|&p| map(p)).collect();
        let control = homogeneous(&points, self.weights.as_deref()).unwrap_or_default();
        Self {
            points,
            control,
            ..self.clone()
        }
    }
}

/// One homogeneous Bézier patch of a surface: `rows` controls along u by
/// `columns` along v, row-major.
#[derive(Clone, Debug)]
pub(crate) struct Patch {
    pub rows: usize,
    pub columns: usize,
    pub control: Vec<[f64; 4]>,
}

fn lerp(a: Homogeneous, b: Homogeneous, alpha: f64) -> Homogeneous {
    std::array::from_fn(|i| a[i] + (b[i] - a[i]) * alpha)
}

/// Boehm insertion of one knot `t` into a clamped sequence.
fn insert(knots: &mut Vec<f64>, control: &mut Vec<Homogeneous>, degree: usize, t: f64) {
    let k = knots
        .iter()
        .rposition(|&x| x <= t)
        .expect("t lies in the clamped range");
    let mut next = Vec::with_capacity(control.len() + 1);
    for i in 0..=control.len() {
        if i + degree <= k {
            next.push(control[i]);
        } else if i > k {
            next.push(control[i - 1]);
        } else {
            let alpha = (t - knots[i]) / (knots[i + degree] - knots[i]);
            next.push(lerp(control[i - 1], control[i], alpha));
        }
    }
    knots.insert(k + 1, t);
    *control = next;
}

/// The Bézier segments of a clamped B-spline sequence: every interior knot
/// is raised to multiplicity `degree`, which leaves the curve unchanged.
/// `None` for unclamped ends.
pub(crate) fn bezier_segments(
    knots: &[f64],
    degree: usize,
    control: &[Homogeneous],
) -> Option<Vec<Vec<Homogeneous>>> {
    let n = control.len();
    if degree == 0 || knots.len() != n + degree + 1 {
        return None;
    }
    let (a, b) = (knots[degree], knots[n]);
    if knots[..=degree].iter().any(|&k| k != a) || knots[n..].iter().any(|&k| k != b) {
        return None;
    }
    let mut knots = knots.to_vec();
    let mut control = control.to_vec();
    let mut i = degree + 1;
    while i + degree + 1 < knots.len() {
        let t = knots[i];
        let multiplicity = knots.iter().filter(|&&k| k == t).count();
        for _ in multiplicity..degree {
            insert(&mut knots, &mut control, degree, t);
        }
        while i < knots.len() && knots[i] == t {
            i += 1;
        }
    }
    let segments = (control.len() - 1) / degree;
    Some(
        (0..segments)
            .map(|s| control[s * degree..=s * degree + degree].to_vec())
            .collect(),
    )
}

impl NurbsCurve {
    /// Homogeneous Bézier segments of the whole curve.
    pub(crate) fn bezier(&self) -> Option<Vec<Vec<Homogeneous>>> {
        bezier_segments(&self.knots, self.degree, &self.control)
    }
}

impl NurbsSurface {
    /// Homogeneous Bézier patches of the whole surface.
    pub(crate) fn bezier(&self) -> Option<Vec<Patch>> {
        let [p, q] = self.degrees;
        let [nu, nv] = self.counts;
        // Along u for every column of the net.
        let columns: Vec<Vec<Vec<Homogeneous>>> = (0..nv)
            .map(|j| {
                let column: Vec<_> = (0..nu).map(|i| self.control[i * nv + j]).collect();
                bezier_segments(&self.knots[0], p, &column)
            })
            .collect::<Option<_>>()?;
        let strips = columns.first()?.len();
        let mut patches = Vec::new();
        for s in 0..strips {
            // Each of the p + 1 rows of this strip, along v.
            let rows: Vec<Vec<Vec<Homogeneous>>> = (0..=p)
                .map(|r| {
                    let row: Vec<_> = columns.iter().map(|column| column[s][r]).collect();
                    bezier_segments(&self.knots[1], q, &row)
                })
                .collect::<Option<_>>()?;
            for t in 0..rows[0].len() {
                patches.push(Patch {
                    rows: p + 1,
                    columns: q + 1,
                    control: rows.iter().flat_map(|row| row[t].clone()).collect(),
                });
            }
        }
        Some(patches)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The allocating Cox-de Boor evaluation the stack triangle replaced.
    fn reference_values(knots: &[f64], span: usize, degree: usize, t: f64) -> Vec<f64> {
        let mut current = vec![1.];
        for j in 1..=degree {
            let mut next = vec![0.; j + 1];
            for r in 0..j {
                let i = span + 1 + r - j;
                let w = current[r] / (knots[i + j] - knots[i]);
                next[r + 1] += (t - knots[i]) * w;
                next[r] += (knots[i + j] - t) * w;
            }
            current = next;
        }
        current
    }

    fn reference_basis(knots: &[f64], degree: usize, t: f64) -> (usize, Vec<f64>, Vec<f64>) {
        let count = knots.len() - degree - 1;
        let t = t.clamp(knots[degree], knots[count]);
        let mut span = degree;
        while span + 1 < count && t >= knots[span + 1] {
            span += 1;
        }
        let value = reference_values(knots, span, degree, t);
        let lower = reference_values(knots, span, degree - 1, t);
        let mut derivative = vec![0.; degree + 1];
        for (m, slot) in derivative.iter_mut().enumerate() {
            let i = span + m - degree;
            if m >= 1 {
                *slot += lower[m - 1] / (knots[i + degree] - knots[i]);
            }
            if m < degree {
                *slot -= lower[m] / (knots[i + degree + 1] - knots[i + 1]);
            }
            *slot *= degree as f64;
        }
        (span, value, derivative)
    }

    /// Every degree, clamped and unclamped knots with repeated interior
    /// knots, parameters on and between knots and outside the range: the
    /// stack basis is bit-identical to the allocating one.
    #[test]
    fn stack_basis_is_bit_identical_to_the_allocating_evaluation() {
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut unit = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1_u64 << 53) as f64
        };
        let mut checked = 0;
        for degree in 1..=MAX_DEGREE {
            for clamped in [true, false] {
                let count = degree + 1 + (degree % 5);
                let mut knots = Vec::new();
                let mut knot = -1.7;
                while knots.len() < count + degree + 1 {
                    let repeat = if clamped && (knots.is_empty() || knots.len() >= count) {
                        degree + 1
                    } else if unit() < 0.2 {
                        2
                    } else {
                        1
                    };
                    for _ in 0..repeat.min(count + degree + 1 - knots.len()) {
                        knots.push(knot);
                    }
                    knot += 0.1 + unit();
                }
                let [low, high] = [knots[degree], knots[count]];
                let mut parameters: Vec<f64> = knots.clone();
                parameters.extend([low - 1., high + 1., low, high]);
                parameters.extend((0..40).map(|_| low + (high - low) * unit()));
                for t in parameters {
                    let (span, value, derivative) = basis(&knots, degree, t);
                    let (expected, values, derivatives) = reference_basis(&knots, degree, t);
                    assert_eq!(span, expected);
                    for (got, want) in [(&value, &values), (&derivative, &derivatives)] {
                        assert!(got[degree + 1..].iter().all(|&v| v.to_bits() == 0));
                        for (a, b) in got[..=degree].iter().zip(want.iter()) {
                            assert_eq!(a.to_bits(), b.to_bits(), "degree {degree} at {t}");
                        }
                    }
                    checked += 1;
                }
            }
        }
        assert!(checked > 2_000);
    }
}
