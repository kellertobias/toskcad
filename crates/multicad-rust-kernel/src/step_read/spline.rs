//! Bounded nonrational cubic STEP curves with clamped ends and double interior
//! knots, and single-span rational Bézier curves of degree two or three with
//! positive weights. Decode the transmitted geometry independently of the
//! STEP writer.

use super::*;

const MAX_CUBIC_POLES: usize = 256;

#[derive(Clone, Debug, PartialEq)]
pub struct SplineEdge {
    poles_mm: Vec<Point>,
    knots: Vec<f64>,
    multiplicities: Vec<u8>,
    /// Positive weights of a single rational span; `None` for a polynomial
    /// cubic. A rational curve has exactly one Bézier span.
    weights: Option<Vec<f64>>,
}

impl SplineEdge {
    pub fn poles_mm(&self) -> &[Point] {
        &self.poles_mm
    }

    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    pub fn multiplicities(&self) -> &[u8] {
        &self.multiplicities
    }

    /// Weights of a single rational Bézier span, in pole order.
    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    /// Homogeneous `[w x, w y, w z, w]` controls of a rational span.
    fn homogeneous(&self) -> Option<Vec<[f64; 4]>> {
        let weights = self.weights.as_ref()?;
        Some(
            self.poles_mm
                .iter()
                .zip(weights)
                .map(|(&[x, y, z], &w)| [x * w, y * w, z * w, w])
                .collect(),
        )
    }

    /// Resolve the double interior knots into cubic Bezier spans. The shared
    /// endpoint follows from C1 continuity under unequal knot intervals.
    /// Polynomial curves only: a rational span has no polynomial spans.
    pub fn bezier_spans(
        &self,
    ) -> impl DoubleEndedIterator<Item = [Point; 4]> + ExactSizeIterator + '_ {
        debug_assert!(
            self.weights.is_none(),
            "polynomial spans of a rational curve"
        );
        let count = self.knots.len() - 1;
        let join = move |index: usize| {
            if index == 0 {
                return self.poles_mm[0];
            }
            if index == count {
                return self.poles_mm[2 * count + 1];
            }
            let [a, b, c] = [index - 1, index, index + 1].map(|i| self.knots[i]);
            let t = (b - a) / (c - a);
            lerp(self.poles_mm[2 * index], self.poles_mm[2 * index + 1], t)
        };
        (0..count).map(move |i| {
            [
                join(i),
                self.poles_mm[2 * i + 1],
                self.poles_mm[2 * i + 2],
                join(i + 1),
            ]
        })
    }

    pub(super) fn translated(&self, top: &Self, delta: Point, tolerance: f64) -> bool {
        self.knots == top.knots
            && self.multiplicities == top.multiplicities
            && self.weights == top.weights
            && self
                .poles_mm
                .iter()
                .zip(&top.poles_mm)
                .all(|(&low, &high)| length(sub(sub(high, low), delta)) <= tolerance)
    }

    pub(super) fn bounds(&self) -> [Point; 2] {
        if let Some(controls) = self.homogeneous() {
            return rational_bounds(&controls);
        }
        let mut bounds = [self.poles_mm[0]; 2];
        for span in self.bezier_spans() {
            include_point(&mut bounds, span[3]);
            for axis in 0..3 {
                let derivative = std::array::from_fn(|i| span[i + 1][axis] - span[i][axis]);
                for t in quadratic_roots(derivative).into_iter().flatten() {
                    if t > 0. && t < 1. {
                        include_point(&mut bounds, evaluate(span, t));
                    }
                }
            }
        }
        bounds
    }

    pub(super) fn signed_area(&self, origin: Point, normal: Point) -> f64 {
        if let Some(mut controls) = self.homogeneous() {
            for control in &mut controls {
                for axis in 0..3 {
                    control[axis] -= origin[axis] * control[3];
                }
            }
            return rational_signed_area(&controls, normal);
        }
        // P × P' has polynomial degree at most five on each cubic span.
        // Three-point Gauss integration is exact for that degree.
        let d = (3_f64 / 5.).sqrt() * 0.5;
        self.bezier_spans()
            .map(|span| {
                let span = span.map(|p| sub(p, origin));
                let derivative: [Point; 3] = std::array::from_fn(|i| sub(span[i + 1], span[i]));
                [(0.5 - d, 5. / 18.), (0.5, 4. / 9.), (0.5 + d, 5. / 18.)]
                    .into_iter()
                    .map(|(t, weight)| {
                        let tangent = lerp(
                            lerp(derivative[0], derivative[1], t),
                            lerp(derivative[1], derivative[2], t),
                            t,
                        );
                        1.5 * weight * dot(cross(evaluate(span, t), tangent), normal)
                    })
                    .sum::<f64>()
            })
            .sum()
    }
}

fn lerp4(a: [f64; 4], b: [f64; 4], t: f64) -> [f64; 4] {
    std::array::from_fn(|i| (b[i] - a[i]).mul_add(t, a[i]))
}

fn rational_at(controls: &[[f64; 4]], t: f64) -> ([f64; 4], [f64; 4]) {
    let degree = controls.len() - 1;
    let mut work = controls.to_vec();
    for level in (2..=degree).rev() {
        for i in 0..level {
            work[i] = lerp4(work[i], work[i + 1], t);
        }
    }
    let derivative = std::array::from_fn(|i| degree as f64 * (work[1][i] - work[0][i]));
    (lerp4(work[0], work[1], t), derivative)
}

fn rational_split(controls: &[[f64; 4]]) -> [Vec<[f64; 4]>; 2] {
    let degree = controls.len() - 1;
    let mut left = Vec::with_capacity(degree + 1);
    let mut right = vec![[0.; 4]; degree + 1];
    let mut work = controls.to_vec();
    for level in 0..=degree {
        left.push(work[0]);
        right[degree - level] = work[degree - level];
        for i in 0..degree - level {
            work[i] = lerp4(work[i], work[i + 1], 0.5);
        }
    }
    [left, right]
}

fn projected(control: [f64; 4]) -> Point {
    [0, 1, 2].map(|axis| control[axis] / control[3])
}

/// Axis extents by hull branch and bound: positive weights keep the
/// curve inside its projected control hull, which converges under halving.
fn rational_bounds(controls: &[[f64; 4]]) -> [Point; 2] {
    let first = projected(controls[0]);
    let last = projected(controls[controls.len() - 1]);
    let mut bounds = [first; 2];
    include_point(&mut bounds, last);
    let scale = controls
        .iter()
        .flat_map(|&c| projected(c))
        .map(f64::abs)
        .fold(1., f64::max);
    let ends = bounds;
    for axis in 0..3 {
        for sign in [-1., 1.] {
            let mut best = sign * ends[usize::from(sign > 0.)][axis];
            let mut high = f64::NEG_INFINITY;
            let mut stack = vec![(controls.to_vec(), 0)];
            while let Some((piece, depth)) = stack.pop() {
                let values: Vec<_> = piece.iter().map(|&c| sign * projected(c)[axis]).collect();
                let upper = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let ends = values[0].max(values[values.len() - 1]);
                best = best.max(ends);
                if upper <= best || upper - ends <= 1e-14 * scale || depth == 48 {
                    high = high.max(upper);
                    continue;
                }
                for half in rational_split(&piece) {
                    stack.push((half, depth + 1));
                }
            }
            bounds[usize::from(sign > 0.)][axis] = sign * high.max(best);
        }
    }
    bounds
}

/// `½∫ (P × P') · n` with `P × P' = (X × X') / W²` for homogeneous `X, W`.
/// The analytic integrand is resolved by composite eight-point Gauss rules.
fn rational_signed_area(controls: &[[f64; 4]], normal: Point) -> f64 {
    const PARTS: usize = 32;
    const GAUSS: [(f64, f64); 4] = [
        (0.183_434_642_495_649_8, 0.362_683_783_378_362),
        (0.525_532_409_916_329, 0.313_706_645_877_887_3),
        (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
        (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    ];
    let mut total = 0.;
    for part in 0..PARTS {
        for (node, weight) in GAUSS {
            for t in [(1. - node) / 2., (1. + node) / 2.] {
                let t = (part as f64 + t) / PARTS as f64;
                let ([x, y, z, w], [dx, dy, dz, _]) = rational_at(controls, t);
                let flux = dot(cross([x, y, z], [dx, dy, dz]), normal) / (w * w);
                total += 0.25 * weight / PARTS as f64 * flux;
            }
        }
    }
    total
}

fn lerp(a: Point, b: Point, t: f64) -> Point {
    std::array::from_fn(|i| (b[i] - a[i]).mul_add(t, a[i]))
}

fn evaluate([a, b, c, d]: [Point; 4], t: f64) -> Point {
    let [a, b, c] = [lerp(a, b, t), lerp(b, c, t), lerp(c, d, t)];
    lerp(lerp(a, b, t), lerp(b, c, t), t)
}

/// Roots of a quadratic supplied by its three Bernstein coefficients.
fn quadratic_roots([a, b, c]: [f64; 3]) -> [Option<f64>; 2] {
    let [a, b, c] = [a - 2. * b + c, 2. * (b - a), a];
    if a == 0. {
        return [(b != 0.).then(|| -c / b), None];
    }
    let discriminant = b.mul_add(b, -4. * a * c);
    if discriminant < 0. {
        return [None; 2];
    }
    let q = -0.5 * (b + discriminant.sqrt().copysign(b));
    if q == 0. {
        return [Some(-b / (2. * a)), None];
    }
    [Some(q / a), Some(c / q)]
}

impl Shell<'_, '_> {
    pub(super) fn spline_curve(&self, id: u32) -> Result<SplineEdge> {
        let [
            Value::String(_),
            Value::Number(3.),
            poles,
            Value::Enum("UNSPECIFIED"),
            Value::Enum("F"),
            Value::Enum("F"),
            Value::List(multiplicities),
            Value::List(knots),
            Value::Enum("UNSPECIFIED" | "PIECEWISE_BEZIER_KNOTS"),
        ] = self.args(id, "B_SPLINE_CURVE_WITH_KNOTS")?
        else {
            return Err(ReadError::Unsupported);
        };
        let poles = refs(poles)?;
        if poles.len() > MAX_CUBIC_POLES {
            return Err(ReadError::Limit);
        }
        if knots.len() < 2 || multiplicities.len() != knots.len() || poles.len() != 2 * knots.len()
        {
            return Err(ReadError::Invalid);
        }
        let multiplicities = multiplicities
            .iter()
            .enumerate()
            .map(|(i, value)| {
                let expected = if i == 0 || i + 1 == knots.len() { 4 } else { 2 };
                if !matches!(value, Value::Number(n) if *n == f64::from(expected)) {
                    return Err(ReadError::Unsupported);
                }
                Ok(expected)
            })
            .collect::<Result<Vec<_>>>()?;
        let knots = knots
            .iter()
            .map(|value| match value {
                Value::Number(knot) if knot.is_finite() && knot.abs() <= 1e12 => Ok(*knot),
                _ => Err(ReadError::Invalid),
            })
            .collect::<Result<Vec<_>>>()?;
        if knots.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(ReadError::Invalid);
        }
        let poles_mm = poles
            .into_iter()
            .map(|id| self.point(id))
            .collect::<Result<Vec<_>>>()?;
        Ok(SplineEdge {
            poles_mm,
            knots,
            multiplicities,
            weights: None,
        })
    }

    /// A complex single-span rational Bézier curve of degree two or three:
    /// `B_SPLINE_CURVE`, `B_SPLINE_CURVE_WITH_KNOTS` with end multiplicities
    /// `degree + 1`, and positive `RATIONAL_B_SPLINE_CURVE` weights.
    pub(super) fn rational_curve(&self, id: u32) -> Result<SplineEdge> {
        let parts = self.components(id)?;
        let component = |name: &str| {
            parts
                .iter()
                .find(|entity| entity.name == name)
                .map(|entity| entity.args.as_slice())
                .ok_or(ReadError::Unsupported)
        };
        let [
            Value::Number(degree),
            poles,
            Value::Enum("UNSPECIFIED" | "CIRCULAR_ARC" | "ELLIPTIC_ARC" | "CONIC_ARC"),
            Value::Enum("F"),
            Value::Enum("F"),
        ] = component("B_SPLINE_CURVE")?
        else {
            return Err(ReadError::Unsupported);
        };
        let [
            Value::List(multiplicities),
            Value::List(knots),
            Value::Enum("UNSPECIFIED" | "PIECEWISE_BEZIER_KNOTS" | "QUASI_UNIFORM_KNOTS"),
        ] = component("B_SPLINE_CURVE_WITH_KNOTS")?
        else {
            return Err(ReadError::Unsupported);
        };
        let [Value::List(weights)] = component("RATIONAL_B_SPLINE_CURVE")? else {
            return Err(ReadError::Invalid);
        };
        if *degree != 2. && *degree != 3. {
            return Err(ReadError::Unsupported);
        }
        let order = *degree as u8 + 1;
        let poles = refs(poles)?;
        if poles.len() != usize::from(order) || weights.len() != poles.len() {
            return Err(ReadError::Unsupported);
        }
        let numbers = |values: &[Value<'_>]| {
            values
                .iter()
                .map(|value| match value {
                    Value::Number(n) if n.is_finite() && n.abs() <= 1e12 => Ok(*n),
                    _ => Err(ReadError::Invalid),
                })
                .collect::<Result<Vec<_>>>()
        };
        let multiplicities = numbers(multiplicities)?;
        let knots = numbers(knots)?;
        let weights = numbers(weights)?;
        if multiplicities != [f64::from(order); 2] {
            return Err(ReadError::Unsupported);
        }
        if knots.len() != 2 || knots[0] >= knots[1] || weights.iter().any(|w| *w <= 0.) {
            return Err(ReadError::Invalid);
        }
        let poles_mm = poles
            .into_iter()
            .map(|id| self.point(id))
            .collect::<Result<Vec<_>>>()?;
        Ok(SplineEdge {
            poles_mm,
            knots,
            multiplicities: vec![order; 2],
            weights: Some(weights),
        })
    }

    /// A polynomial cubic or single-span rational sweep rim curve.
    pub(super) fn sweep_curve(&self, id: u32) -> Result<SplineEdge> {
        if self.components(id)?.len() > 1 {
            self.rational_curve(id)
        } else {
            self.spline_curve(id)
        }
    }

    pub(super) fn spline_extrusion(&self, id: u32) -> Result<(u32, Point)> {
        self.spline_extrusion_vector(id)
            .map(|(curve, direction, _)| (curve, direction))
    }

    /// The swept curve, unit direction and parameter speed of `S(u, v) = C(u) + v V`.
    pub(super) fn spline_extrusion_vector(&self, id: u32) -> Result<(u32, Point, f64)> {
        let [Value::String(_), curve, vector] = self.args(id, "SURFACE_OF_LINEAR_EXTRUSION")?
        else {
            return Err(ReadError::Invalid);
        };
        let [Value::String(_), direction, Value::Number(scale)] =
            self.args(reference(vector)?, "VECTOR")?
        else {
            return Err(ReadError::Invalid);
        };
        if !scale.is_finite() || *scale <= 0. {
            return Err(ReadError::Invalid);
        }
        Ok((
            reference(curve)?,
            self.direction(reference(direction)?)?,
            *scale,
        ))
    }
}

#[cfg(test)]
mod tests;
