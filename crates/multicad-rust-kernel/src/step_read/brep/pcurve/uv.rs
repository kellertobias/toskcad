//! Dedicated dimension-two geometry parser. The ordinary 3D parser is untouched.
use super::*;
use crate::brep::nurbs::{MAX_CONTROL_POINTS, MAX_DEGREE};

pub(super) struct Definition {
    pub surface: u32,
    pub curve_id: u32,
    pub geometry: Curve,
}
pub(super) fn definition(
    source: &Source<'_>,
    id: u32,
    budget: &mut Budget<'_>,
) -> Result<Definition> {
    budget.charge(256)?;
    let [_, surface, representation] = source.args(id, "PCURVE")? else {
        return Err(ReadError::Invalid);
    };
    let surface = reference(surface)?;
    let [_, items, context] =
        source.args(reference(representation)?, "DEFINITIONAL_REPRESENTATION")?
    else {
        return Err(ReadError::Invalid);
    };
    let context = reference(context)?;
    let Some([dimension]) = source.component(context, "GEOMETRIC_REPRESENTATION_CONTEXT")? else {
        return Err(ReadError::Unsupported);
    };
    let dimension = super::precision::scalar(source, context, dimension, budget)?;
    check(
        dimension.valid() && dimension.lo == 2. && dimension.hi == 2.,
        ReadError::Unsupported,
    )?;
    let items = refs(items)?;
    let [curve_id] = items.as_slice() else {
        return Err(ReadError::Unsupported);
    };
    Ok(Definition {
        surface,
        curve_id: *curve_id,
        geometry: curve(source, *curve_id, budget)?,
    })
}
fn pair(source: &Source<'_>, id: u32, kind: &str, budget: &mut Budget<'_>) -> Result<Point> {
    budget.charge(256)?;
    let [_, coordinates] = source.args(id, kind)? else {
        return Err(ReadError::Invalid);
    };
    let [a, b] = list(coordinates)? else {
        return Err(ReadError::Invalid);
    };
    // These are parameter-space coordinates, not bounded 3D lengths. A large
    // legal planar gauge can cancel exactly against the support's origin.
    Ok([number(a)?, number(b)?, 0.])
}
fn direction(source: &Source<'_>, id: u32, budget: &mut Budget<'_>) -> Result<Point> {
    let p = pair(source, id, "DIRECTION", budget)?;
    let scale = p.into_iter().map(f64::abs).fold(0_f64, f64::max);
    check(scale > 0., ReadError::Invalid)?;
    normalized(p.map(|x| x / scale))
}
fn frame(source: &Source<'_>, id: u32, budget: &mut Budget<'_>) -> Result<Frame> {
    budget.charge(256)?;
    let [_, origin, reference_axis] = source.args(id, "AXIS2_PLACEMENT_2D")? else {
        return Err(ReadError::Invalid);
    };
    let x = match reference_axis {
        Value::Omitted => [1., 0., 0.],
        value => direction(source, reference(value)?, budget)?,
    };
    Ok(Frame {
        origin: pair(source, reference(origin)?, "CARTESIAN_POINT", budget)?,
        x,
        y: [-x[1], x[0], 0.],
        z: [0., 0., 1.],
    })
}
fn curve(source: &Source<'_>, id: u32, budget: &mut Budget<'_>) -> Result<Curve> {
    budget.charge(256)?;
    let parts = source.components(id)?;
    let [entity] = parts else {
        return spline(source, id, budget);
    };
    match (entity.name, entity.args.as_slice()) {
        ("LINE", [_, origin, vector_id]) => {
            let [_, d, magnitude] = source.args(reference(vector_id)?, "VECTOR")? else {
                return Err(ReadError::Invalid);
            };
            let magnitude = number(magnitude)?;
            check(magnitude > 0., ReadError::Invalid)?;
            Ok(Curve::Line {
                origin: pair(source, reference(origin)?, "CARTESIAN_POINT", budget)?,
                direction: direction(source, reference(d)?, budget)?.map(|x| x * magnitude),
            })
        }
        ("CIRCLE", [_, axis, r]) => Ok(Curve::Circle {
            frame: frame(source, reference(axis)?, budget)?,
            radius: positive(r)?,
        }),
        ("ELLIPSE", [_, axis, a, b]) => Ok(Curve::Ellipse {
            frame: frame(source, reference(axis)?, budget)?,
            semiaxes: [positive(a)?, positive(b)?],
        }),
        ("HYPERBOLA", [_, axis, a, b]) => Ok(Curve::Hyperbola {
            frame: frame(source, reference(axis)?, budget)?,
            semiaxes: [positive(a)?, positive(b)?],
        }),
        ("PARABOLA", [_, axis, f]) => Ok(Curve::Parabola {
            frame: frame(source, reference(axis)?, budget)?,
            focal: positive(f)?,
        }),
        ("B_SPLINE_CURVE_WITH_KNOTS" | "BEZIER_CURVE", _) => spline(source, id, budget),
        // A trimmed curve carries an additional finite-domain contract, never
        // silently stripped down to its unbounded basis support.
        _ => Err(ReadError::Unsupported),
    }
}
fn positive(value: &Value<'_>) -> Result<f64> {
    let value = number(value)?;
    check(value > 0., ReadError::Invalid)?;
    Ok(value)
}
fn spline(source: &Source<'_>, id: u32, budget: &mut Budget<'_>) -> Result<Curve> {
    let parts = source.components(id)?;
    let (base, knots, weights) = match parts {
        [entity] => match (entity.name, entity.args.as_slice()) {
            ("B_SPLINE_CURVE_WITH_KNOTS", [_, base @ .., multiplicities, knots, _])
                if base.len() == 5 =>
            {
                (base, Some((multiplicities, knots)), None)
            }
            ("BEZIER_CURVE", [_, base @ ..]) if base.len() == 5 => (base, None, None),
            _ => return Err(ReadError::Unsupported),
        },
        _ => {
            let base = source
                .component(id, "B_SPLINE_CURVE")?
                .ok_or(ReadError::Unsupported)?;
            let knots = match source.component(id, "B_SPLINE_CURVE_WITH_KNOTS")? {
                Some([m, k, _]) => Some((m, k)),
                Some(_) => return Err(ReadError::Invalid),
                None if source.component(id, "BEZIER_CURVE")?.is_some() => None,
                None => return Err(ReadError::Unsupported),
            };
            let weights = match source.component(id, "RATIONAL_B_SPLINE_CURVE")? {
                Some([w]) => Some(w),
                Some(_) => return Err(ReadError::Invalid),
                None => None,
            };
            (base, knots, weights)
        }
    };
    let [degree, controls, ..] = base else {
        return Err(ReadError::Invalid);
    };
    let degree = count(degree)?;
    check((1..=MAX_DEGREE).contains(&degree), ReadError::Invalid)?;
    let controls = refs(controls)?;
    check(
        controls.len() > degree && controls.len() <= MAX_CONTROL_POINTS,
        ReadError::Limit,
    )?;
    budget.charge(controls.len().saturating_mul(128))?;
    let mut points = Vec::new();
    points
        .try_reserve_exact(controls.len())
        .map_err(|_| ReadError::Limit)?;
    check(
        points.capacity() * std::mem::size_of::<Point>() <= MAX_WORKSPACE / 8,
        ReadError::Limit,
    )?;
    for id in controls {
        points.push(pair(source, id, "CARTESIAN_POINT", budget)?);
    }
    let weights = match weights {
        Some(w) => {
            check(list(w)?.len() == points.len(), ReadError::Invalid)?;
            budget.charge(points.len().saturating_mul(64))?;
            Some(reals(w)?)
        }
        None => None,
    };
    let (knots, multiplicities) = match knots {
        Some((m, k)) => {
            check(
                list(m)?.len() <= MAX_CONTROL_POINTS + MAX_DEGREE + 1
                    && list(k)?.len() == list(m)?.len(),
                ReadError::Limit,
            )?;
            budget.charge(list(k)?.len().saturating_mul(128))?;
            (reals(k)?, counts(m)?)
        }
        None => (vec![0., 1.], vec![degree + 1; 2]),
    };
    let curve = NurbsCurve::new(degree, points, weights, knots, multiplicities)?;
    // Natural-domain supplementary values only. Do not fabricate a closure or
    // wrap the last endpoint onto the first based on an approximate sample.
    Ok(Curve::Nurbs {
        curve,
        closed: false,
    })
}
