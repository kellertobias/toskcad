//! Exact v integrals of the original charts, not fitted carrier area.
use super::*;
type S = Series<5>;
fn c(v: I) -> S {
    S::interval(v)
}
fn n(v: f64) -> S {
    S::constant(v)
}
fn value(s: Series<6>) -> S {
    S {
        c: std::array::from_fn(|i| s.c[i]),
    }
}
fn derivative(s: Series<6>) -> S {
    S {
        c: std::array::from_fn(|i| s.c[i + 1] * I::point((i + 1) as f64)),
    }
}
fn atan(x: S) -> Result<S, Error> {
    // d(atan x)/du=x'/(1+x²). Only derivative coefficients0..3
    // contribute to atan coefficients1..4; the unused fifth coefficient
    // of the derivative is not asserted as a source derivative bound.
    let dx = S {
        c: std::array::from_fn(|i| {
            if i < 4 {
                x.c[i + 1] * I::point((i + 1) as f64)
            } else {
                I::point(0.)
            }
        }),
    };
    let derivative = dx / (n(1.) + x * x);
    let mut result = S::interval(x.c[0].atan());
    for i in 1..5 {
        result.c[i] = derivative.c[i - 1] / I::point(i as f64);
    }
    checked(result)
}
fn checked(v: S) -> Result<S, Error> {
    v.checked().map_err(|_| Error::Unresolved)
}
fn frame(mut raw: [S; 4], quarter: u8, lower: bool) -> [S; 4] {
    if lower {
        raw[2] = -raw[2];
    }
    let [x, y] = match quarter {
        0 => [raw[1], raw[2]],
        1 => [-raw[2], raw[1]],
        2 => [-raw[1], -raw[2]],
        _ => [raw[2], -raw[1]],
    };
    [raw[0], x, y, raw[3]]
}
enum TrimKind {
    Envelope,
    Wall,
    Trace,
}
fn trim_density(
    events: &BoundaryEvents<'_>,
    index: usize,
    u: I,
    kind: TrimKind,
) -> Result<[S; 4], Error> {
    let p = evaluate::trim_state::<Series<6>>(events, index, u)?;
    let [x, y, z] = p.map(value);
    let [dx, dy, dz] = p.map(derivative);
    let beta = c(events.beta_tangent);
    Ok(match kind {
        TrimKind::Envelope => {
            let radius = checked((x * x + y * y).sqrt())?;
            let phi = atan(y / x)?;
            let width = c(pi() / I::point(2.)) - n(2.) * phi;
            let speed = if events.trims[index].envelope == EnvelopeRole::EntryRoundover {
                // Exact circular meridian speed, including either endpoint.
                c(I::point(events.source().spec.mouth_radius_mm)
                    * (events.trims[index].parameters[1] - events.trims[index].parameters[0]))
            } else {
                let i = usize::from(events.trims[index].envelope == EnvelopeRole::Floor);
                let slope = events.cone_slopes[i];
                // In the sealed source y decreases as depth increases.
                // x*offset_y+y>=0 and r-q*x>0 imply z_u>=0, including
                // H0's zero derivative at y=0. Keep that analytic factor.
                dz * c((I::point(1.) + slope.square()).sqrt())
            };
            let area = radius * speed * width;
            let xy = radius * speed * (x - y);
            [area, xy, xy, z * area]
        }
        TrimKind::Wall => {
            let along = dx + beta * dz;
            let speed = checked(((n(1.) + beta * beta) * dy * dy + along * along).sqrt())?;
            let area = z * speed;
            [
                area,
                area * (x + beta * z / n(2.)),
                area * y,
                area * z / n(2.),
            ]
        }
        TrimKind::Trace => {
            let tx = x + beta * z;
            let dtx = dx + beta * dz;
            // Cavity admission proves the positive front intercept, including
            // H0. Source y_u<0 fixes this determinant's indicated sign.
            let jacobian = y * dtx - tx * dy;
            [
                jacobian / n(2.),
                jacobian * tx / n(3.),
                jacobian * y / n(3.),
                n(0.),
            ]
        }
    })
}
pub(super) fn evaluate(atlas: &PatchAtlas<'_>, index: usize, u: I) -> Result<[S; 4], Error> {
    let events = atlas.topology.events();
    let cutter = events.source();
    let corner = &cutter.corners[0];
    let envelope = &cutter.envelope;
    let (raw, q, lower) = match atlas.patches[index].law {
        Law::Envelope(index, q) => (
            trim_density(events, index, u, TrimKind::Envelope)?,
            q,
            false,
        ),
        Law::Wall(index, q, lower) => (trim_density(events, index, u, TrimKind::Wall)?, q, lower),
        Law::Trace(index, q, lower) => (trim_density(events, index, u, TrimKind::Trace)?, q, lower),
        Law::Gap(q) => {
            let p = events.event(EventRole::MouthEnvelope).point;
            let phi = (p[1] / p[0]).atan();
            let width = pi() / I::point(2.) - I::point(2.) * phi;
            let theta = c(phi) + S::variable(u) * c(width);
            let [sin, cos] = theta.sin_cos();
            let radius = c(envelope.entry);
            let area = radius * radius * c(width) / n(2.);
            let moment = radius * radius * radius * c(width) / n(3.);
            ([area, moment * cos, moment * sin, n(0.)], q, false)
        }
        Law::Floor(q) => {
            let radius = c(corner.b_half) * (n(1.) - S::variable(u));
            let depth = c(corner.b_depth) + S::variable(u) * c(envelope.tip - corner.b_depth);
            let speed =
                c((corner.b_half.square() + (envelope.tip - corner.b_depth).square()).sqrt());
            let area = radius * speed * c(pi() / I::point(2.));
            let xy = radius * radius * speed;
            ([area, xy, xy, depth * area], q, false)
        }
    };
    let raw = frame(raw, q, lower);
    Ok([
        checked(raw[0])?,
        checked(raw[1])?,
        checked(raw[2])?,
        checked(raw[3])?,
    ])
}
