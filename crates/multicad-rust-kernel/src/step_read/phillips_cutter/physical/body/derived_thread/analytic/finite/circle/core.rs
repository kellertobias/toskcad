//! Source-free bounded circle correspondence. RawCircle defines the family
//! obtained by EXACT normalization of the supplied normal/reference variables;
//! the Frame is an interval evaluation cache of that coupled construction.
//! File decoders must bind those variables to original numeric leaves/owners.
use super::*;
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::point(0.), |s, k| s + a[k] * b[k])
}
fn delta(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|k| a[k] - b[k])
}

pub(super) struct RawCircle {
    frame: support::Frame,
    radius: I,
}
impl RawCircle {
    pub(super) fn new(
        origin: [I; 3],
        normal: [I; 3],
        reference: [I; 3],
        radius: I,
        work: &mut Work<'_>,
    ) -> FResult<Self> {
        work.step(1024)?;
        if !radius.valid()
            || radius.lo <= 0.
            || origin
                .into_iter()
                .chain(normal)
                .chain(reference)
                .any(|x| !x.valid())
        {
            return Err(PhillipsReadError::Unsupported);
        }
        let frame = support::Frame::new(origin, normal, reference).map_err(source_error)?;
        Ok(Self { frame, radius })
    }
}
#[derive(Clone, Copy)]
pub(super) struct Input {
    pub target: ArcTarget,
    pub chart: CommonMouthChart,
    pub vertices: [[I; 3]; 2],
    pub same_sense: bool,
}
pub(super) struct Correspondence {
    pub distance: f64,
    pub derivative: f64,
    pub incidence: f64,
}
pub(super) fn certify(
    raw: &RawCircle,
    input: Input,
    work: &mut Work<'_>,
) -> FResult<Correspondence> {
    let Input {
        target,
        chart,
        vertices,
        same_sense,
    } = input;
    work.step(2048)?;
    if vertices.into_iter().flatten().any(|x| !x.valid())
        || chart.columns.into_iter().flatten().any(|x| !x.is_finite())
        || !chart.outer_top_z.is_finite()
    {
        return Err(PhillipsReadError::Unsupported);
    }
    let r = raw.radius;
    let [o, u, v, n] = chart.columns.map(|p| p.map(I::point));
    let alignment = dot(raw.frame.axes[2], n);
    let sign = if alignment.lo > 0. {
        1.
    } else if alignment.hi < 0. {
        -1.
    } else {
        return Err(PhillipsReadError::Unsupported);
    };
    let source_sign = if target.left_hand() { 1. } else { -1. };
    if (sign * source_sign > 0.) != same_sense {
        return Err(PhillipsReadError::Unsupported);
    }
    let center =
        std::array::from_fn(|k| o[k] + n[k] * (I::point(target.z()) - I::point(chart.outer_top_z)));
    let radius = I::point(target.radius());
    let project = |actual: [I; 3]| -> FResult<([I; 3], [I; 3], f64)> {
        let d = delta(actual, raw.frame.origin);
        let xy = [dot(raw.frame.axes[0], d), dot(raw.frame.axes[1], d)];
        let length = (xy[0].square() + xy[1].square()).sqrt();
        if !length.valid() || length.lo <= 0. {
            return Err(PhillipsReadError::Unsupported);
        }
        let direction = std::array::from_fn(|k| {
            raw.frame.axes[0][k] * xy[0] / length + raw.frame.axes[1][k] * xy[1] / length
        });
        let projected = std::array::from_fn(|k| raw.frame.origin[k] + r * direction[k]);
        Ok((projected, direction, reduce::l1(delta(actual, projected))))
    };
    let actual_start = vertices[0];
    let (projected_start, start_direction, start_incidence) = project(actual_start)?;
    let source_start = target.endpoint(false, work)?;
    let [c0, s0] = [source_start[0] / radius, source_start[1] / radius];
    let side = support::cross(raw.frame.axes[2], start_direction).map(|x| I::point(sign) * x);
    // Exact parameter-origin rotation of the unchanged actual circle. Its
    // original U is arbitrary (the real writer uses the crest start direction).
    // Source sin/cos at the proposed start are evaluation caches of the same
    // real phase expression; the original projected start fixes this raw phase.
    let raw_u = std::array::from_fn(|k| r * (start_direction[k] * c0 - side[k] * s0));
    let raw_v = std::array::from_fn(|k| r * (start_direction[k] * s0 + side[k] * c0));
    let a = delta(raw_u, u.map(|x| radius * x));
    let b = delta(raw_v, v.map(|x| radius * x));
    let basis_error = (I::point(reduce::l1(a)) + I::point(reduce::l1(b))).hi;
    let support_error =
        (I::point(reduce::l1(delta(raw.frame.origin, center))) + I::point(basis_error)).hi;
    let mut chord = 0_f64;
    let mut incidence = start_incidence;
    for (end, &actual) in vertices.iter().enumerate() {
        let (projected, proposed) = if end == 0 {
            (projected_start, source_start)
        } else {
            let (projected, _, reserve) = project(actual)?;
            incidence = incidence.max(reserve);
            (projected, target.endpoint(true, work)?)
        };
        let matched = std::array::from_fn(|k| {
            raw.frame.origin[k]
                + raw_u[k] * (proposed[0] / radius)
                + raw_v[k] * (proposed[1] / radius)
        });
        chord = chord.max(reduce::l1(delta(projected, matched)));
    }
    if (I::point(chord) / r).hi >= 0.5 {
        return Err(PhillipsReadError::Unsupported);
    }
    // For chord/r<1/2, minor angular distance <=2*chord/r. No sampled angle or
    // libm inverse-trig estimate enters either endpoint or branch association.
    let angular = I::point(2.) * I::point(chord) / r;
    let rotation = rotation::Rotation::new(work)?;
    let pi = rotation.pi(work)?;
    let span = target.angular_span(pi);
    if span.lo <= (I::point(4.) * angular).hi
        || (I::point(2.) * pi - span).lo <= (I::point(4.) * angular).hi
    {
        return Err(PhillipsReadError::Unsupported);
    }

    let result = Correspondence {
        distance: (I::point(support_error) + r * angular + I::point(incidence)).hi,
        derivative: (I::point(basis_error) * span
            + r * (I::point(2.) * angular + angular * (span + I::point(2.) * angular)))
            .hi,
        incidence,
    };
    if ![result.distance, result.derivative, result.incidence]
        .into_iter()
        .all(|x| x.is_finite() && x >= 0.)
    {
        return Err(PhillipsReadError::Unsupported);
    }
    Ok(result)
}
#[cfg(test)]
mod tests;
