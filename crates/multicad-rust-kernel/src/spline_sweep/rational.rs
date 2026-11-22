//! Exact single-span rational Bézier source curves in a planar sweep frame.
//!
//! The authored poles and positive weights are retained bit for bit; no
//! basis conversion, degree elevation or polynomial substitution happens.
//! Conics (degree two) and rational cubics share one homogeneous algorithm.

use super::{Error, Point};
use crate::revolution::{Sum, cross, distance_to_segment, dot, length, sub};
use multicad_precision::PortableMath as _;

const MAX_DEPTH: usize = 32;
const MAX_LENGTH_SPLITS: usize = 100_000;

/// Homogeneous control `[w x, w y, w]` of a rational Bézier span.
pub(crate) type Homogeneous = [f64; 3];

/// One rational Bézier span of degree two or three with positive weights.
#[derive(Clone, Debug, PartialEq)]
pub struct Rational {
    poles: Vec<Point>,
    weights: Vec<f64>,
    reversed: bool,
}

/// Fixed-capacity homogeneous span; `degree + 1` leading controls are used.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Span {
    pub(crate) controls: [Homogeneous; 4],
    pub(crate) degree: usize,
}

impl Rational {
    pub fn new(poles: Vec<Point>, weights: Vec<f64>, reversed: bool) -> Result<Self, Error> {
        if !(3..=4).contains(&poles.len())
            || weights.len() != poles.len()
            || !poles
                .iter()
                .flatten()
                .all(|value| value.is_finite() && value.abs() <= 1e8)
            || !weights
                .iter()
                .all(|weight| weight.is_finite() && *weight > 0. && *weight <= 1e8)
            || weights
                .iter()
                .fold(f64::INFINITY, |low: f64, w| low.min(*w))
                < 1e-8 * weights.iter().fold(0_f64, |high, w| high.max(*w))
            || poles.first() == poles.last()
        {
            return Err(Error::Invalid);
        }
        Ok(Self {
            poles,
            weights,
            reversed,
        })
    }

    pub fn degree(&self) -> usize {
        self.poles.len() - 1
    }

    pub fn poles(&self) -> &[Point] {
        &self.poles
    }

    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    pub fn reversed(&self) -> bool {
        self.reversed
    }

    pub fn reverse(&self) -> Self {
        Self {
            reversed: !self.reversed,
            ..self.clone()
        }
    }

    pub fn endpoints(&self) -> [Point; 2] {
        let ends = [self.poles[0], self.poles[self.poles.len() - 1]];
        if self.reversed {
            [ends[1], ends[0]]
        } else {
            ends
        }
    }

    /// The span in traversal order.
    pub(crate) fn span(&self) -> Span {
        let degree = self.degree();
        let mut controls = [[0.; 3]; 4];
        for (i, control) in controls.iter_mut().take(degree + 1).enumerate() {
            let source = if self.reversed { degree - i } else { i };
            let [x, y] = self.poles[source];
            let w = self.weights[source];
            *control = [x * w, y * w, w];
        }
        Span { controls, degree }
    }

    pub fn point(&self, fraction: f64) -> Point {
        if fraction <= 0. {
            return self.endpoints()[0];
        }
        if fraction >= 1. {
            return self.endpoints()[1];
        }
        self.span().point(fraction)
    }

    pub fn support(&self, direction: Point) -> f64 {
        let span = self.span();
        let ends = span.endpoints().map(|p| dot(direction, p));
        let mut low = ends[0].max(ends[1]);
        let scale = span
            .points()
            .iter()
            .map(|p| dot(direction, *p).abs())
            .fold(1., f64::max);
        let mut high = f64::NEG_INFINITY;
        let mut stack = vec![(span, 0)];
        while let Some((piece, depth)) = stack.pop() {
            let upper = piece
                .points()
                .iter()
                .map(|p| dot(direction, *p))
                .fold(f64::NEG_INFINITY, f64::max);
            let [a, b] = piece.endpoints().map(|p| dot(direction, p));
            low = low.max(a).max(b);
            if upper <= low || upper - a.max(b) <= 1e-14 * scale || depth == MAX_DEPTH {
                high = high.max(upper);
                continue;
            }
            let [first, second] = piece.split(0.5);
            stack.push((first, depth + 1));
            stack.push((second, depth + 1));
        }
        // The retained hull bound encloses the curve; it exceeds the exact
        // support by at most the refinement allowance.
        high.max(low)
    }

    /// `[area, x, y]` Green integrals about `origin`, as for polynomial spans.
    pub(super) fn moments(&self, origin: Point) -> [f64; 3] {
        // With homogeneous X(t) and weight W(t), P × P' = (X × X') / W².
        // The rational integrand is analytic on the span; composite
        // Gauss-Legendre on 32 equal parts reaches rounding level for any
        // admitted weight ratio below 1e8.
        let span = self.span().translated(origin);
        let mut sum = [Sum::default(); 3];
        const PARTS: usize = 32;
        for part in 0..PARTS {
            for (node, weight) in GAUSS_8 {
                let t = (part as f64 + node) / PARTS as f64;
                let [x, y, w] = span.homogeneous(t);
                let [dx, dy, _] = span.derivative(t);
                let flux = weight / PARTS as f64 * (x * dy - y * dx) / (w * w);
                let [px, py] = [x / w, y / w];
                for (sum, term) in sum
                    .iter_mut()
                    .zip([0.5 * flux, px * flux / 3., py * flux / 3.])
                {
                    sum.add(term);
                }
            }
        }
        sum.map(Sum::value)
    }

    /// Chord and control-polygon bounds for the exact rational length.
    /// Homogeneous de Casteljau steps are corner cuts of the projected
    /// control polygon, so its length never increases under subdivision.
    pub fn length_bounds_mm(&self, error_mm: f64) -> Result<[f64; 2], Error> {
        if !error_mm.is_finite() || error_mm <= 0. {
            return Err(Error::Invalid);
        }
        let mut budget = MAX_LENGTH_SPLITS;
        self.length_bounds_with_budget(error_mm, &mut budget)
    }

    pub(super) fn length_bounds_with_budget(
        &self,
        error_mm: f64,
        budget: &mut usize,
    ) -> Result<[f64; 2], Error> {
        measure(self.span(), error_mm, 0, budget)
    }

    pub(super) fn sample(
        &self,
        deflection: f64,
        angle: f64,
        limit: usize,
    ) -> Result<Vec<Point>, crate::display::MeshError> {
        let mut points = vec![self.endpoints()[0]];
        flatten(
            self.span(),
            deflection,
            angle.portable_tan(),
            0,
            limit,
            &mut points,
        )?;
        Ok(points)
    }
}

impl Span {
    pub(crate) fn controls(&self) -> &[Homogeneous] {
        &self.controls[..=self.degree]
    }

    /// Projected Euclidean control points; positive weights keep the curve
    /// in their convex hull and make it variation diminishing.
    pub(crate) fn points(&self) -> Vec<Point> {
        self.controls()
            .iter()
            .map(|&[x, y, w]| [x / w, y / w])
            .collect()
    }

    pub(crate) fn endpoints(&self) -> [Point; 2] {
        let [x, y, w] = self.controls[0];
        let [u, v, s] = self.controls[self.degree];
        [[x / w, y / w], [u / s, v / s]]
    }

    fn translated(mut self, [ox, oy]: Point) -> Self {
        for control in self.controls.iter_mut().take(self.degree + 1) {
            let w = control[2];
            control[0] -= ox * w;
            control[1] -= oy * w;
        }
        self
    }

    fn homogeneous(&self, t: f64) -> Homogeneous {
        let mut work = self.controls;
        for level in (1..=self.degree).rev() {
            for i in 0..level {
                work[i] = lerp3(work[i], work[i + 1], t);
            }
        }
        work[0]
    }

    fn derivative(&self, t: f64) -> Homogeneous {
        let degree = self.degree;
        let mut work = [[0.; 3]; 4];
        for (slot, pair) in work.iter_mut().zip(self.controls[..=degree].windows(2)) {
            *slot = std::array::from_fn(|axis| degree as f64 * (pair[1][axis] - pair[0][axis]));
        }
        for level in (1..degree).rev() {
            for i in 0..level {
                work[i] = lerp3(work[i], work[i + 1], t);
            }
        }
        work[0]
    }

    pub(crate) fn point(&self, t: f64) -> Point {
        let [x, y, w] = self.homogeneous(t);
        [x / w, y / w]
    }

    /// Homogeneous de Casteljau subdivision.
    pub(crate) fn split(&self, t: f64) -> [Self; 2] {
        let degree = self.degree;
        let mut left = [[0.; 3]; 4];
        let mut right = [[0.; 3]; 4];
        let mut work = self.controls;
        for level in 0..=degree {
            left[level] = work[0];
            right[degree - level] = work[degree - level];
            for i in 0..degree - level {
                work[i] = lerp3(work[i], work[i + 1], t);
            }
        }
        [
            Self {
                controls: left,
                degree,
            },
            Self {
                controls: right,
                degree,
            },
        ]
    }
}

fn lerp3(a: Homogeneous, b: Homogeneous, t: f64) -> Homogeneous {
    std::array::from_fn(|i| (b[i] - a[i]).mul_add(t, a[i]))
}

fn polygon_length(points: &[Point]) -> f64 {
    points
        .windows(2)
        .map(|pair| length(sub(pair[1], pair[0])))
        .sum()
}

fn measure(span: Span, error: f64, depth: usize, budget: &mut usize) -> Result<[f64; 2], Error> {
    *budget = budget.checked_sub(1).ok_or(Error::Limit)?;
    let points = span.points();
    let [a, b] = span.endpoints();
    let low = length(sub(b, a));
    let high = polygon_length(&points).max(low);
    if high - low <= error {
        return Ok([low, high]);
    }
    if depth == MAX_DEPTH {
        return Err(Error::Limit);
    }
    let [first, second] = span.split(0.5);
    let left = measure(first, error / 2., depth + 1, budget)?;
    let right = measure(second, error / 2., depth + 1, budget)?;
    Ok([left[0] + right[0], left[1] + right[1]])
}

fn flatten(
    span: Span,
    deflection: f64,
    cone: f64,
    depth: usize,
    limit: usize,
    points: &mut Vec<Point>,
) -> Result<(), crate::display::MeshError> {
    use crate::display::MeshError;

    let controls = span.points();
    let [start, end] = span.endpoints();
    let chord = sub(end, start);
    let resolved = length(chord) > 0.
        && controls[1..span.degree]
            .iter()
            .all(|&point| distance_to_segment(point, [start, end]) <= deflection)
        && controls.windows(2).all(|pair| {
            let derivative = sub(pair[1], pair[0]);
            derivative == [0.; 2] || dot(chord, derivative) * cone >= cross(chord, derivative).abs()
        });
    if resolved {
        if points.len() > limit {
            return Err(MeshError::BudgetExceeded);
        }
        points.push(end);
        return Ok(());
    }
    if depth == MAX_DEPTH {
        return Err(MeshError::UnresolvedBoundary);
    }
    for half in span.split(0.5) {
        flatten(half, deflection, cone, depth + 1, limit, points)?;
    }
    Ok(())
}

/// Eight-point Gauss-Legendre nodes and weights on [0, 1].
const GAUSS_8: [(f64, f64); 8] = [
    (0.019_855_071_751_231_856, 0.050_614_268_145_188_13),
    (0.101_666_761_293_186_63, 0.111_190_517_226_687_24),
    (0.237_233_795_041_835_5, 0.156_853_322_938_943_64),
    (0.408_282_678_752_175, 0.181_341_891_689_181),
    (0.591_717_321_247_825, 0.181_341_891_689_181),
    (0.762_766_204_958_164_5, 0.156_853_322_938_943_64),
    (0.898_333_238_706_813_4, 0.111_190_517_226_687_24),
    (0.980_144_928_248_768_1, 0.050_614_268_145_188_13),
];

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_1_SQRT_2, PI};

    fn quarter() -> Rational {
        Rational::new(
            vec![[10., 0.], [10., 10.], [0., 10.]],
            vec![1., FRAC_1_SQRT_2, 1.],
            false,
        )
        .unwrap()
    }

    #[test]
    fn quarter_circle_is_exact_in_points_area_length_and_support() {
        let arc = quarter();
        for i in 0..=16 {
            let [x, y] = arc.point(i as f64 / 16.);
            assert!((x.portable_hypot(y) - 10.).abs() < 1e-12);
        }
        // Quarter disk: the arc plus the two radii through the origin.
        let [area, mx, my] = arc.moments([0., 0.]);
        assert!((area - 25. * PI).abs() < 1e-11, "{area}");
        // ∫∫ x dA over the quarter disk is r³/3.
        assert!(
            (mx - 1000. / 3.).abs() < 1e-10 && (my - 1000. / 3.).abs() < 1e-10,
            "{mx} {my}"
        );
        let [low, high] = arc.length_bounds_mm(1e-6).unwrap();
        assert!(low <= 5. * PI && 5. * PI <= high && high - low <= 1e-6);
        let diagonal = arc.support([FRAC_1_SQRT_2, FRAC_1_SQRT_2]);
        assert!(
            (10. - 1e-12..=10. + 1e-11).contains(&diagonal),
            "{diagonal}"
        );
        assert_eq!(arc.support([1., 0.]), 10.);
    }

    #[test]
    fn reversal_and_samples_follow_traversal_and_bad_weights_refuse() {
        let arc = quarter().reverse();
        assert_eq!(arc.point(0.), [0., 10.]);
        assert_eq!(arc.point(1.), [10., 0.]);
        let samples = arc.sample(1e-3, 0.2, 1000).unwrap();
        assert_eq!(samples.first(), Some(&[0., 10.]));
        assert_eq!(samples.last(), Some(&[10., 0.]));
        assert!(
            samples
                .iter()
                .all(|[x, y]| (x.portable_hypot(*y) - 10.).abs() < 1e-3)
        );
        // Reversed area flips sign.
        assert!((arc.moments([0., 0.])[0] + 25. * PI).abs() < 1e-11);
        for weights in [vec![1., 0., 1.], vec![1., -1., 1.], vec![1., f64::NAN, 1.]] {
            assert_eq!(
                Rational::new(vec![[10., 0.], [10., 10.], [0., 10.]], weights, false),
                Err(Error::Invalid)
            );
        }
        assert_eq!(
            Rational::new(vec![[0., 0.], [1., 1.]], vec![1., 1.], false),
            Err(Error::Invalid)
        );
    }

    #[test]
    fn homogeneous_split_reproduces_the_curve() {
        let span = quarter().span();
        let [left, right] = span.split(0.5);
        for i in 0..=8 {
            let t = i as f64 / 8.;
            let a = left.point(t);
            let b = span.point(t / 2.);
            assert!(length(sub(a, b)) < 1e-12);
            let c = right.point(t);
            let d = span.point(0.5 + t / 2.);
            assert!(length(sub(c, d)) < 1e-12);
        }
    }
}
