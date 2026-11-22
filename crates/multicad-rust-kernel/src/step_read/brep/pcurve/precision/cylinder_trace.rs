//! Original polynomial degree-eleven rails against a finite cylinder trace.
//! Hermite C0 and the independent Taylor bound for the first derivative are
//! separate proofs; no derivative of an unknown remainder is taken.
use super::super::Proof;
use super::bezier::{self, Coordinates};
use super::*;
use crate::step_read::brep::check;

const SCRATCH: usize = 16 * 1024;
const MAX_BYTES: usize = 8 * 1024 * 1024;
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V, s: I) -> V {
    a.map(|x| x * s)
}
fn power(a: I, n: usize) -> I {
    (0..n).fold(I::point(1.), |v, _| v * a)
}
fn factorial(n: usize) -> f64 {
    (1..=n).fold(1., |v, i| v * i as f64)
}
fn choose(n: usize, k: usize) -> f64 {
    (0..k).fold(1., |v, i| v * (n - i) as f64 / (i + 1) as f64)
}
fn trig(a: I, cosine: bool) -> Result<I> {
    check(a.valid() && a.magnitude() <= 1e6, ReadError::Unsupported)?;
    let m = a.midpoint();
    let radius = (a - I::point(m)).magnitude();
    let at = if cosine {
        I::point(m).cos()
    } else {
        I::point(m).sin()
    };
    Ok(I::new(
        (at.lo - radius).next_down().max(-1.),
        (at.hi + radius).next_up().min(1.),
    ))
}
struct Wave {
    origin: V,
    axes: [V; 3],
    radius: I,
    theta: I,
    height: I,
    angle_step: I,
    height_step: I,
}
impl Wave {
    fn derivative(&self, order: usize, t: f64) -> Result<V> {
        let theta = self.theta + self.angle_step * I::point(t);
        let (c, s) = (trig(theta, true)?, trig(theta, false)?);
        let (c, s) = match order % 4 {
            0 => (c, s),
            1 => (-s, c),
            2 => (-c, -s),
            _ => (s, -c),
        };
        let wave = scale(
            add(scale(self.axes[0], c), scale(self.axes[1], s)),
            self.radius * power(self.angle_step, order),
        );
        Ok(match order {
            0 => add(
                self.origin,
                add(
                    wave,
                    scale(self.axes[2], self.height + self.height_step * I::point(t)),
                ),
            ),
            1 => add(wave, scale(self.axes[2], self.height_step)),
            _ => wave,
        })
    }
    fn derivative12(&self) -> V {
        // Each component has its own legitimate scalar remainder. Combining
        // their bounds does not assume one common mean-value point.
        let amplitude =
            I::point(self.radius.magnitude()) * power(I::point(self.angle_step.magnitude()), 12);
        std::array::from_fn(|i| {
            amplitude
                * (I::point(self.axes[0][i].magnitude()) + I::point(self.axes[1][i].magnitude()))
        })
    }
    fn speed_lower(&self) -> f64 {
        ((self.radius * self.angle_step).square() + self.height_step.square())
            .sqrt()
            .lo
    }
}

fn hermite11(wave: &Wave, budget: &mut Budget<'_>) -> Result<[V; 12]> {
    budget.charge(16_384)?;
    let mut out = [[I::point(0.); 3]; 12];
    for end in 0..2 {
        let mut derivatives = [[I::point(0.); 3]; 6];
        for (k, d) in derivatives.iter_mut().enumerate() {
            *d = wave.derivative(k, end as f64)?;
        }
        for j in 0..6 {
            let index = if end == 0 { j } else { 11 - j };
            for (k, d) in derivatives.iter().enumerate().take(j + 1) {
                let falling = (0..k).fold(1., |v, i| v * (11 - i) as f64);
                let coefficient = I::point(choose(j, k)) / I::point(falling);
                let sign = if end == 1 && k % 2 != 0 { -1. } else { 1. };
                out[index] = add(out[index], scale(*d, coefficient * I::point(sign)));
            }
        }
    }
    Ok(out)
}

fn derivative_taylor10(wave: &Wave, budget: &mut Budget<'_>) -> Result<[V; 11]> {
    budget.charge(32_768)?;
    let mut centered = [[I::point(0.); 3]; 11];
    for (k, d) in centered.iter_mut().enumerate() {
        *d = scale(
            wave.derivative(k + 1, 0.5)?,
            I::point(1.) / I::point(factorial(k)),
        );
    }
    let mut monomials = [[I::point(0.); 3]; 11];
    for (j, value) in monomials.iter_mut().enumerate() {
        for (k, d) in centered.iter().enumerate().skip(j) {
            let coefficient = I::point(choose(k, j)) * power(I::point(-0.5), k - j);
            *value = add(*value, scale(*d, coefficient));
        }
    }
    let mut bernstein = [[I::point(0.); 3]; 11];
    for (j, value) in bernstein.iter_mut().enumerate() {
        for (k, d) in monomials.iter().enumerate().take(j + 1) {
            *value = add(
                *value,
                scale(*d, I::point(choose(j, k)) / I::point(choose(10, k))),
            );
        }
    }
    Ok(bernstein)
}

fn cell(wave: &Wave, controls: &[[I; 4]], tol: f64, budget: &mut Budget<'_>) -> Result<Proof> {
    budget.charge(4096)?;
    if controls.len() != 12 || controls.iter().any(|p| p[3].lo != 1. || p[3].hi != 1.) {
        return Err(ReadError::Unsupported);
    }
    let reference = hermite11(wave, budget)?;
    let mut position: f64 = 0.;
    for (actual, ideal) in controls.iter().zip(reference) {
        let error = norm(difference([actual[0], actual[1], actual[2]], ideal));
        check(error.valid(), ReadError::Unsupported)?;
        position = position.max(error.hi);
    }
    let m12 = wave.derivative12();
    let remainder = norm(m12.map(|v| v / (I::point(factorial(12)) * I::point(4096.))));
    let upper = I::point(position) + remainder;
    if !upper.valid() || upper.hi > tol {
        return Ok(Proof::Unsupported);
    }
    // Compare the actual polynomial derivative to a separately constructed
    // Taylor polynomial of the analytic derivative. Its remainder uses the
    // twelfth analytic derivative and is NEVER a differentiated Hermite error.
    let derivative = derivative_taylor10(wave, budget)?;
    let mut first: f64 = 0.;
    for (pair, ideal) in controls.windows(2).zip(derivative) {
        let actual = std::array::from_fn(|i| I::point(11.) * (pair[1][i] - pair[0][i]));
        let error = norm(difference(actual, ideal));
        check(error.valid(), ReadError::Unsupported)?;
        first = first.max(error.hi);
    }
    let tail = norm(m12.map(|v| v / (I::point(factorial(11)) * I::point(2048.))));
    let derivative_error = I::point(first) + tail;
    let speed = wave.speed_lower();
    if !derivative_error.valid()
        || !speed.is_finite()
        || speed <= 0.
        || derivative_error.hi >= speed
    {
        return Ok(Proof::Unsupported);
    }
    budget.charge(64)?;
    Ok(Proof::Certified(upper.hi))
}

/// Final Original-REAL C0/IP2 certificate. The entire actual bounded target
/// spline is covered by the full finite positive UV segment, so any genuine
/// finite EDGE_CURVE restriction is covered without a cached endpoint oracle.
pub(super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    let (
        Surface::Cylinder { .. },
        Curve::Nurbs {
            curve: uv_curve,
            closed: false,
        },
        Curve::Nurbs {
            curve: target,
            closed: false,
        },
    ) = (surface, &definition.geometry, &edge.curve)
    else {
        return Ok(None);
    };
    if uv_curve.degree() != 1
        || uv_curve.points().len() != 2
        || uv_curve.knots().1 != [2, 2]
        || target.degree() != 11
    {
        return Ok(None);
    }
    budget.charge(4096)?;
    check(
        super::super::si_radians(source, budget)?,
        ReadError::Unsupported,
    )?;
    curve_nurbs(source, definition.curve_id, &definition.geometry, 2, budget)?;
    let uv_fields = fields(source, definition.curve_id, 1, budget)?;
    let uv_controls = list(uv_fields.controls)?;
    check(uv_controls.len() == 2, ReadError::Unsupported)?;
    let a = point(source, reference(&uv_controls[0])?, 2, budget)?;
    let b = point(source, reference(&uv_controls[1])?, 2, budget)?;
    let [_, placement, r] = source.args(definition.surface, "CYLINDRICAL_SURFACE")? else {
        return Err(ReadError::Invalid);
    };
    let f = frame(source, reference(placement)?, 3, budget)?;
    let radius = scalar(source, definition.surface, r, budget)?;
    check(radius.valid() && radius.lo > 0., ReadError::Unsupported)?;
    let curve_id = unwrapped(source, curve_id, budget)?;
    let data = bezier::original_curve(source, curve_id, target, Coordinates::Raw(3), budget)?;
    check(
        data.matches(source, curve_id) && matches!(data.coordinates(), Coordinates::Raw(3)),
        ReadError::Unsupported,
    )?;
    check(
        data.degree() == 11
            && !data.cells().is_empty()
            && data.peak_bytes() <= MAX_BYTES
            && data
                .retained_bytes()
                .checked_add(SCRATCH)
                .is_some_and(|n| n <= MAX_BYTES),
        ReadError::Limit,
    )?;
    // The packet's original ordinal provenance binds its cells back to the
    // actual source knot records. No inferred or decoded knot chart is used.
    // Child construction and the retained-packet/caller phase do not overlap;
    // both independent peaks above remain inside the same 8MiB child ceiling.
    let raw_fields = fields(source, curve_id, 1, budget)?;
    let raw_knots = raw_fields.knots[0]
        .map(|(_, knots)| list(knots))
        .transpose()?;
    let knot_count = raw_knots.map_or(2, |values| values.len());
    check(
        data.cells().len().checked_add(1) == Some(knot_count),
        ReadError::Unsupported,
    )?;
    let original_knot = |index: usize, budget: &mut Budget<'_>| -> Result<I> {
        match raw_knots {
            Some(values) => scalar(
                source,
                curve_id,
                values.get(index).ok_or(ReadError::Unsupported)?,
                budget,
            ),
            // Only an actual BEZIER_CURVE lacking a knot component earns
            // its schema-defined natural chart [0,1].
            None if index <= 1 => Ok(I::point(index as f64)),
            None => Err(ReadError::Unsupported),
        }
    };
    let first = original_knot(0, budget)?;
    let last = original_knot(knot_count - 1, budget)?;
    let width = last - first;
    check(width.valid() && width.lo > 0., ReadError::Unsupported)?;
    let mut maximum: f64 = 0.;
    for (ordinal, current) in data.cells().iter().enumerate() {
        budget.charge(256)?;
        check(
            current.knot_indices() == [ordinal, ordinal + 1],
            ReadError::Unsupported,
        )?;
        let domain = [
            original_knot(ordinal, budget)?,
            original_knot(ordinal + 1, budget)?,
        ];
        check(
            domain
                .iter()
                .zip(current.domain())
                .all(|(raw, stored)| raw.lo == stored.lo && raw.hi == stored.hi),
            ReadError::Unsupported,
        )?;
        let mut normalized = [I::point(0.); 2];
        for (i, v) in domain.into_iter().enumerate() {
            let q = (v - first) / width;
            check(q.valid(), ReadError::Unsupported)?;
            // Extractor membership proves the true normalized knot lies in
            // [0,1]; intersect the enclosure, never move an input parameter.
            normalized[i] = I::new(q.lo.max(0.), q.hi.min(1.));
        }
        check(normalized.into_iter().all(I::valid), ReadError::Unsupported)?;
        let delta = normalized[1] - normalized[0];
        check(delta.valid() && delta.lo > 0., ReadError::Unsupported)?;
        let angle_step = (b[0] - a[0]) * delta;
        let height_step = (b[1] - a[1]) * delta;
        let wave = Wave {
            origin: f.origin,
            axes: f.axes,
            radius,
            theta: a[0] + (b[0] - a[0]) * normalized[0],
            height: a[1] + (b[1] - a[1]) * normalized[0],
            angle_step,
            height_step,
        };
        match cell(&wave, current.controls(), tol, budget)? {
            Proof::Certified(error) => maximum = maximum.max(error),
            _ => return Ok(Some(Proof::Unsupported)),
        }
    }
    budget.charge(64)?;
    Ok(Some(Proof::Certified(maximum)))
}

#[cfg(test)]
mod tests;
