//! Original-REAL cone isocurves in explicitly converted angular units.
//! The returned bound is final: neither decoded coefficients nor a unit-name
//! approximation supplies authority for these certificates.
use super::*;
use precision::{frame, point, scalar, vector as raw_vector};
type V = [I; 3];

fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V, b: I) -> V {
    a.map(|x| x * b)
}
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn unit(a: V) -> Result<V> {
    let n = norm(a);
    check(n.valid() && n.lo > 0., ReadError::Unsupported)?;
    Ok(a.map(|x| x / n))
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

fn angle_unit(source: &Source<'_>, budget: &mut Budget<'_>) -> Result<I> {
    let mut found = None;
    for parts in source.entities.values() {
        budget.charge(16)?;
        for e in parts {
            budget.charge(16)?;
            if e.name != "GLOBAL_UNIT_ASSIGNED_CONTEXT" {
                continue;
            }
            let [units] = e.args.as_slice() else {
                return Err(ReadError::Invalid);
            };
            let units = list(units)?;
            check(units.len() <= 64, ReadError::Limit)?;
            for value in units {
                budget.charge(64)?;
                let id = reference(value)?;
                if source.component(id, "PLANE_ANGLE_UNIT")?.is_none() {
                    continue;
                }
                if let Some(previous) = found {
                    check(previous == id, ReadError::Unsupported)?;
                }
                found = Some(id);
            }
        }
    }
    let id = found.ok_or(ReadError::Unsupported)?;
    if matches!(
        source.component(id, "SI_UNIT")?,
        Some([Value::Omitted, Value::Enum("RADIAN")])
    ) {
        return Ok(I::point(1.));
    }
    let [_, measure] = source
        .component(id, "CONVERSION_BASED_UNIT")?
        .ok_or(ReadError::Unsupported)?
    else {
        return Err(ReadError::Invalid);
    };
    let measure_id = reference(measure)?;
    let parts = source.entities.get(&measure_id).ok_or(ReadError::Invalid)?;
    let [e] = parts.as_slice() else {
        return Err(ReadError::Unsupported);
    };
    check(
        matches!(
            e.name,
            "PLANE_ANGLE_MEASURE_WITH_UNIT" | "MEASURE_WITH_UNIT"
        ),
        ReadError::Unsupported,
    )?;
    let [
        Value::Typed(Entity {
            name: "PLANE_ANGLE_MEASURE",
            args,
            ..
        }),
        base,
    ] = e.args.as_slice()
    else {
        return Err(ReadError::Unsupported);
    };
    let [value] = args.as_slice() else {
        return Err(ReadError::Invalid);
    };
    let base = reference(base)?;
    check(
        source.component(base, "PLANE_ANGLE_UNIT")?.is_some(),
        ReadError::Unsupported,
    )?;
    check(
        matches!(
            source.component(base, "SI_UNIT")?,
            Some([Value::Omitted, Value::Enum("RADIAN")])
        ),
        ReadError::Unsupported,
    )?;
    let factor = scalar(source, measure_id, value, budget)?;
    check(factor.valid() && factor.lo > 0., ReadError::Unsupported)?;
    Ok(factor)
}

fn line(source: &Source<'_>, id: u32, dim: usize, budget: &mut Budget<'_>) -> Result<(V, V)> {
    let [_, p, d] = source.args(id, "LINE")? else {
        return Err(ReadError::Invalid);
    };
    let vector_id = reference(d)?;
    let mut direction = raw_vector(source, vector_id, dim, budget)?;
    let [_, axis, _] = source.args(vector_id, "VECTOR")? else {
        return Err(ReadError::Invalid);
    };
    let axis_id = reference(axis)?;
    let [_, coordinates] = source.args(axis_id, "DIRECTION")? else {
        return Err(ReadError::Invalid);
    };
    let coordinates = list(coordinates)?;
    check(coordinates.len() == dim, ReadError::Unsupported)?;
    for i in 0..dim {
        // Division by a positive norm and multiplication by a finite vector
        // magnitude preserve an ORIGINAL exact zero. Never infer this from
        // a decoded zero or from a small interval.
        if zero(scalar(source, axis_id, &coordinates[i], budget)?) {
            direction[i] = I::point(0.);
        }
    }
    Ok((point(source, reference(p)?, dim, budget)?, direction))
}
fn zero(a: I) -> bool {
    a.lo == 0. && a.hi == 0.
}

pub(super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    definition: &uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    if !matches!(
        (surface, &definition.geometry),
        (Surface::Cone { .. }, Curve::Line { .. })
    ) {
        return Ok(None);
    }
    budget.charge(16384)?;
    let factor = angle_unit(source, budget)?;
    let [_, placement, radius, angle] = source.args(definition.surface, "CONICAL_SURFACE")? else {
        return Err(ReadError::Invalid);
    };
    let f = frame(source, reference(placement)?, 3, budget)?;
    let radius = scalar(source, definition.surface, radius, budget)?;
    let angle = scalar(source, definition.surface, angle, budget)? * factor;
    let cosine = trig(angle, true)?;
    check(
        angle.lo > 0. && angle.hi < std::f64::consts::FRAC_PI_2 && cosine.lo > 0.,
        ReadError::Unsupported,
    )?;
    let tangent = trig(angle, false)? / cosine;
    let (p, d) = line(source, definition.curve_id, 2, budget)?;
    let curve_id = precision::unwrapped(source, curve_id, budget)?;
    if zero(d[0]) && matches!(edge.curve, Curve::Line { .. }) {
        let theta = p[0] * factor;
        let radial = add(
            scale(f.axes[0], trig(theta, true)?),
            scale(f.axes[1], trig(theta, false)?),
        );
        let a = add(
            f.origin,
            add(
                scale(radial, radius + p[1] * tangent),
                scale(f.axes[2], p[1]),
            ),
        );
        let b = scale(add(f.axes[2], scale(radial, tangent)), d[1]);
        let (target_a, target_b) = line(source, curve_id, 3, budget)?;
        let nn = dot(f.axes[2], f.axes[2]);
        let bb = dot(b, b);
        let dd = dot(target_b, target_b);
        check(
            nn.lo > 0. && bb.lo > 0. && dd.lo > 0. && (d[1].lo > 0. || d[1].hi < 0.),
            ReadError::Unsupported,
        )?;
        let shift = (dot(difference(target_a, f.origin), f.axes[2]) / nn - p[1]) / d[1];
        let rate = dot(target_b, f.axes[2]) / (nn * d[1]);
        check(rate.lo > 0., ReadError::Unsupported)?;
        let mut ends = [I::point(0.); 2];
        for (i, &vertex) in edge.vertices.iter().enumerate() {
            budget.charge(512)?;
            let [_, vertex_point] = source.args(
                u32::try_from(vertex).map_err(|_| ReadError::Limit)?,
                "VERTEX_POINT",
            )?
            else {
                return Err(ReadError::Invalid);
            };
            let v = point(source, reference(vertex_point)?, 3, budget)?;
            ends[i] = dot(difference(v, target_a), target_b) / dd;
        }
        let range = I::new(ends[0].lo.min(ends[1].lo), ends[0].hi.max(ends[1].hi));
        let delta_a = difference(add(a, scale(b, shift)), target_a);
        let delta_b = difference(scale(b, rate), target_b);
        let upper = norm(add(delta_a, scale(delta_b, range)));
        if upper.valid() && upper.hi <= tolerance {
            return Ok(Some(Proof::Certified(upper.hi)));
        }
        // Distance to the complete lifted supporting line is a lower bound
        // on distance to any finite use of that line, independent of a map.
        let mut lower: f64 = 0.;
        for end in ends {
            let delta = difference(add(target_a, scale(target_b, end)), a);
            let perpendicular = difference(delta, scale(b, dot(delta, b) / bb));
            let distance = norm(perpendicular);
            if distance.valid() {
                lower = lower.max(distance.lo);
            }
        }
        return Ok(Some(if lower > tolerance {
            Proof::Disagreement(lower)
        } else {
            Proof::Unsupported
        }));
    }
    if zero(d[1]) && matches!(edge.curve, Curve::Circle { .. }) {
        let [_, target_placement, target_radius] = source.args(curve_id, "CIRCLE")? else {
            return Err(ReadError::Invalid);
        };
        let target = frame(source, reference(target_placement)?, 3, budget)?;
        let target_radius = scalar(source, curve_id, target_radius, budget)?;
        let orientation = dot(f.axes[2], target.axes[2]);
        let sign = if orientation.lo > 0. {
            1.
        } else if orientation.hi < 0. {
            -1.
        } else {
            return Ok(Some(Proof::Unsupported));
        };
        check(
            (d[0] * factor * I::point(sign)).lo > 0.,
            ReadError::Unsupported,
        )?;
        let r = radius + p[1] * tangent;
        let center = add(f.origin, scale(f.axes[2], p[1]));
        let x = scale(f.axes[0], r);
        let y = scale(f.axes[1], r);
        let phase = unit([dot(x, target.axes[0]), dot(x, target.axes[1]), I::point(0.)])?;
        let target_x = scale(
            add(
                scale(target.axes[0], phase[0]),
                scale(target.axes[1], phase[1]),
            ),
            target_radius,
        );
        let target_y = scale(
            add(
                scale(target.axes[0], -phase[1]),
                scale(target.axes[1], phase[0]),
            ),
            target_radius * I::point(sign),
        );
        let upper = norm(difference(center, target.origin))
            + norm(difference(x, target_x))
            + norm(difference(y, target_y));
        // The original UV LINE has unbounded angular domain. A positive
        // nonzero angular rate covers every finite master circle window.
        return Ok(Some(if upper.valid() && upper.hi <= tolerance {
            Proof::Certified(upper.hi)
        } else {
            Proof::Unsupported
        }));
    }
    Ok(None)
}

#[cfg(test)]
mod tests;
