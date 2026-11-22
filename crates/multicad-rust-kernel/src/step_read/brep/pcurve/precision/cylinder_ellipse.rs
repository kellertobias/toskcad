//! Whole original polynomial UV traces against an independent ellipse.
//! The ellipse supplies the geometric coefficients; no cut recipe is inferred.
use super::super::Proof;
use super::bezier::{self, Coordinates};
use super::*;
use crate::step_read::brep::check;

const SCRATCH: usize = 64 * 1024;
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V, s: I) -> V {
    a.map(|v| v * s)
}
fn power(x: I, n: usize) -> I {
    (0..n).fold(I::point(1.), |v, _| v * x)
}
fn factorial(n: usize) -> f64 {
    (1..=n).fold(1., |v, i| v * i as f64)
}
fn choose(n: usize, k: usize) -> f64 {
    (0..k).fold(1., |v, i| v * (n - i) as f64 / (i + 1) as f64)
}
fn trig(x: I, cosine: bool) -> Result<I> {
    check(x.valid() && x.magnitude() <= 1e6, ReadError::Unsupported)?;
    let m = x.midpoint();
    let r = (x - I::point(m)).magnitude();
    let y = if cosine {
        I::point(m).cos()
    } else {
        I::point(m).sin()
    };
    Ok(I::new(
        (y.lo - r).next_down().max(-1.),
        (y.hi + r).next_up().min(1.),
    ))
}
struct Height {
    center: I,
    cosine: I,
    sine: I,
    start: I,
    width: I,
}
impl Height {
    fn derivative(&self, n: usize, t: f64) -> Result<I> {
        let theta = self.start + self.width * I::point(t);
        let (c, s) = (trig(theta, true)?, trig(theta, false)?);
        let (c, s) = match n % 4 {
            0 => (c, s),
            1 => (-s, c),
            2 => (-c, -s),
            _ => (s, -c),
        };
        Ok((self.cosine * c + self.sine * s) * power(self.width, n)
            + if n == 0 { self.center } else { I::point(0.) })
    }
    fn m12(&self) -> I {
        (I::point(self.cosine.magnitude()) + I::point(self.sine.magnitude()))
            * power(I::point(self.width.magnitude()), 12)
    }
    fn hermite(&self) -> Result<[I; 12]> {
        let mut out = [I::point(0.); 12];
        for end in 0..2 {
            for j in 0..6 {
                let index = if end == 0 { j } else { 11 - j };
                for k in 0..=j {
                    let falling = (0..k).fold(1., |v, i| v * (11 - i) as f64);
                    let sign = if end == 1 && k % 2 != 0 { -1. } else { 1. };
                    out[index] = out[index]
                        + self.derivative(k, end as f64)?
                            * (I::point(choose(j, k) * sign) / I::point(falling));
                }
            }
        }
        Ok(out)
    }
    fn derivative_taylor(&self) -> Result<[I; 11]> {
        let mut powers = [I::point(0.); 11];
        for (j, value) in powers.iter_mut().enumerate() {
            for k in j..=10 {
                *value = *value
                    + self.derivative(k + 1, 0.5)? / I::point(factorial(k))
                        * I::point(choose(k, j))
                        * power(I::point(-0.5), k - j);
            }
        }
        let mut b = [I::point(0.); 11];
        for (j, value) in b.iter_mut().enumerate() {
            for (k, p) in powers.iter().enumerate().take(j + 1) {
                *value = *value + *p * (I::point(choose(j, k)) / I::point(choose(10, k)));
            }
        }
        Ok(b)
    }
}

pub(super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    d: &super::super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    let (
        Surface::Cylinder { .. },
        Curve::Nurbs {
            curve: uv,
            closed: false,
        },
        Curve::Ellipse { .. },
    ) = (surface, &d.geometry, &edge.curve)
    else {
        return Ok(None);
    };
    if uv.degree() != 11 || uv.points().len() != 12 {
        return Ok(None);
    }
    budget.charge(131_072)?;
    check(
        super::super::si_radians(source, budget)?,
        ReadError::Unsupported,
    )?;
    let packet = bezier::original_curve(source, d.curve_id, uv, Coordinates::Raw(2), budget)?;
    check(
        packet.matches(source, d.curve_id)
            && matches!(packet.coordinates(), Coordinates::Raw(2))
            && packet.degree() == 11
            && packet.cells().len() == 1,
        ReadError::Unsupported,
    )?;
    let peak = packet
        .peak_bytes()
        .max(packet.retained_bytes().saturating_add(SCRATCH));
    check(peak <= 8 * 1024 * 1024, ReadError::Limit)?;
    let cell = &packet.cells()[0];
    let controls = cell.controls();
    check(
        controls.len() == 12 && controls.iter().all(|p| p[3].lo == 1. && p[3].hi == 1.),
        ReadError::Unsupported,
    )?;
    // Rebind both domain endpoints to actual Original knot records.
    let fields = fields(source, d.curve_id, 1, budget)?;
    let (_, raw_knots) = fields.knots[0].ok_or(ReadError::Unsupported)?;
    let knots = list(raw_knots)?;
    let ord = cell.knot_indices();
    let start = scalar(
        source,
        d.curve_id,
        knots.get(ord[0]).ok_or(ReadError::Unsupported)?,
        budget,
    )?;
    let end = scalar(
        source,
        d.curve_id,
        knots.get(ord[1]).ok_or(ReadError::Unsupported)?,
        budget,
    )?;
    let domain = cell.domain();
    check(
        start.lo == domain[0].lo
            && start.hi == domain[0].hi
            && end.lo == domain[1].lo
            && end.hi == domain[1].hi,
        ReadError::Unsupported,
    )?;
    let width = end - start;
    check(width.valid() && width.lo > 0., ReadError::Unsupported)?;
    let [_, placement, r] = source.args(d.surface, "CYLINDRICAL_SURFACE")? else {
        return Err(ReadError::Invalid);
    };
    let cylinder = frame(source, reference(placement)?, 3, budget)?;
    let r = scalar(source, d.surface, r, budget)?;
    let id = unwrapped(source, curve_id, budget)?;
    let [_, placement, a, b] = source.args(id, "ELLIPSE")? else {
        return Err(ReadError::Invalid);
    };
    let ellipse = frame(source, reference(placement)?, 3, budget)?;
    let a = scalar(source, id, a, budget)?;
    let b = scalar(source, id, b, budget)?;
    check(
        r.valid() && r.lo > 0. && a.valid() && a.lo > 0. && b.valid() && b.lo > 0.,
        ReadError::Unsupported,
    )?;
    let ea = scale(ellipse.axes[0], a);
    let eb = scale(ellipse.axes[1], b);
    let center = difference(ellipse.origin, cylinder.origin);
    let phase = controls[0][0] - start;
    let (c, s) = (trig(phase, true)?, trig(phase, false)?);
    let radial_cos = scale(
        add(scale(cylinder.axes[0], c), scale(cylinder.axes[1], s)),
        r,
    );
    let radial_sin = scale(
        add(scale(cylinder.axes[0], -s), scale(cylinder.axes[1], c)),
        r,
    );
    let height = Height {
        center: dot(center, cylinder.axes[2]),
        cosine: dot(ea, cylinder.axes[2]),
        sine: dot(eb, cylinder.axes[2]),
        start,
        width,
    };
    let center_error = norm(difference(center, scale(cylinder.axes[2], height.center)));
    let cos_error = norm(difference(
        ea,
        add(radial_cos, scale(cylinder.axes[2], height.cosine)),
    ));
    let sin_error = norm(difference(
        eb,
        add(radial_sin, scale(cylinder.axes[2], height.sine)),
    ));
    let coefficient_error = center_error + cos_error + sin_error;
    let ideal = height.hermite()?;
    let derivative = height.derivative_taylor()?;
    let (mut angle_error, mut angle_derivative, mut z_error, mut z_derivative) =
        (0f64, 0f64, 0f64, 0f64);
    for (j, p) in controls.iter().enumerate() {
        let u = controls[0][0] + width * (I::point(j as f64) / I::point(11.));
        let ue = p[0] - u;
        let ze = p[1] - ideal[j];
        check(ue.valid() && ze.valid(), ReadError::Unsupported)?;
        angle_error = angle_error.max(ue.magnitude());
        z_error = z_error.max(ze.magnitude());
    }
    for (j, p) in controls.windows(2).enumerate() {
        let ue = I::point(11.) * (p[1][0] - p[0][0]) - width;
        let ze = I::point(11.) * (p[1][1] - p[0][1]) - derivative[j];
        check(ue.valid() && ze.valid(), ReadError::Unsupported)?;
        angle_derivative = angle_derivative.max(ue.magnitude());
        z_derivative = z_derivative.max(ze.magnitude());
    }
    let z_error = I::point(z_error) + height.m12() / (I::point(factorial(12)) * I::point(4096.));
    let z_derivative =
        I::point(z_derivative) + height.m12() / (I::point(factorial(11)) * I::point(2048.));
    let z_norm = norm(cylinder.axes[2]);
    let xy_norm = norm(cylinder.axes[0]) + norm(cylinder.axes[1]);
    let position = coefficient_error + r * xy_norm * I::point(angle_error) + z_norm * z_error;
    // Independent first-derivative enclosure, never differentiated C0 error.
    let first = (cos_error + sin_error) * width
        + r * xy_norm
            * (I::point(angle_derivative) + I::point(width.magnitude()) * I::point(angle_error))
        + z_norm * z_derivative;
    let speed = I::point(a.lo.min(b.lo)) * I::point(width.lo);
    if !first.valid() || !speed.valid() || first.hi >= speed.lo {
        return Ok(Some(Proof::Unsupported));
    }
    let finite = super::super::finite_axis::endpoint_reserve(source, curve_id, edge, budget)?;
    // Natural ellipse knots and the retained edge can use different period
    // representatives. Bound the actual first-ray angle, without identifying
    // binary64 TAU with the true period; retain the complete positive width.
    let phase_delta = start - I::point(edge.range[0]);
    let forward = trig(phase_delta, true)?;
    let sideways = trig(phase_delta, false)?;
    check(forward.lo > 0., ReadError::Unsupported)?;
    let first_delta = I::point(sideways.magnitude()) / I::point(forward.lo);
    let width_delta = width - (I::point(edge.range[1]) - I::point(edge.range[0]));
    check(
        first_delta.valid() && width_delta.valid(),
        ReadError::Unsupported,
    )?;
    let domain_error = (first_delta + I::point(width_delta.magnitude())).hi;
    // Compare tangents at the final paired finite-window points too. This
    // ellipse reserve bounds (a+b)*each original endpoint angle error.
    // Both |T'| and |T''| are <= a+b, independently of the Hermite remainder.
    let endpoint_angle = I::point(finite) / I::point((a + b).lo);
    let paired_phase = I::point(domain_error) + endpoint_angle;
    let paired_width = I::point(width_delta.magnitude()) + I::point(2.) * endpoint_angle;
    let paired_first =
        first + (a + b) * (I::point(width.magnitude()) * paired_phase + paired_width);
    let actual_span =
        I::point(edge.range[1]) - I::point(edge.range[0]) - I::point(2.) * endpoint_angle;
    let paired_speed = I::point(a.lo.min(b.lo)) * I::point(actual_span.lo);
    if !endpoint_angle.valid()
        || !paired_first.valid()
        || !actual_span.valid()
        || actual_span.lo <= 0.
        || !paired_speed.valid()
        || paired_first.hi >= paired_speed.lo
    {
        return Ok(Some(Proof::Unsupported));
    }
    let upper = position + I::point(finite) + (a + b) * I::point(domain_error);
    budget.charge(64)?;
    Ok(Some(if upper.valid() && upper.hi <= tol {
        Proof::Certified(upper.hi)
    } else {
        Proof::Unsupported
    }))
}

#[cfg(test)]
mod tests;
