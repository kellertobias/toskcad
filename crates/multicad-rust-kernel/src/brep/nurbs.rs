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

/// Nonzero basis functions of `degree` on the knot span `span`.
fn values(knots: &[f64], span: usize, degree: usize, t: f64) -> Vec<f64> {
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

/// Span index, basis values and basis derivatives at `t`.
fn basis(knots: &[f64], degree: usize, t: f64) -> (usize, Vec<f64>, Vec<f64>) {
    let count = knots.len() - degree - 1;
    let t = t.clamp(knots[degree], knots[count]);
    let mut span = degree;
    while span + 1 < count && t >= knots[span + 1] {
        span += 1;
    }
    let value = values(knots, span, degree, t);
    let lower = values(knots, span, degree - 1, t);
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
            for (k, &value) in values.iter().enumerate() {
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
