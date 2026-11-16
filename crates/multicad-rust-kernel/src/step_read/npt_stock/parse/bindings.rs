//! Actual pcurve/support equality is earned from owning records and original
//! decimal leaves. Overlapping caches never establish geometric identity.
use super::*;
fn pcurve<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    id: u32,
    budget: &mut Budget<'_>,
) -> Result<(u32, u32)> {
    let entity = kind(records, id, "PCURVE", budget)?;
    let [Value::String(_), support, representation] = entity.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let rep = kind(
        records,
        reference(representation)?,
        "DEFINITIONAL_REPRESENTATION",
        budget,
    )?;
    let [Value::String(_), items, context] = rep.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let [item] = list(items)? else {
        return Err(Error::Unsupported);
    };
    let context = records.get(&reference(context)?).ok_or(Error::Invalid)?;
    let [geometry, parametric, representation] = context.as_slice() else {
        return Err(Error::Unsupported);
    };
    if geometry.name != "GEOMETRIC_REPRESENTATION_CONTEXT"
        || parametric.name != "PARAMETRIC_REPRESENTATION_CONTEXT"
        || representation.name != "REPRESENTATION_CONTEXT"
    {
        return Err(Error::Unsupported);
    }
    let [dimension] = geometry.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let mut numbers = numeric::Numbers::new(geometry, budget)?;
    if numbers.integer(dimension, budget)? != 2 {
        return Err(Error::Unsupported);
    }
    numbers.finish()?;
    Ok((reference(support)?, reference(item)?))
}
fn line<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    id: u32,
    budget: &mut Budget<'_>,
) -> Result<usize> {
    let entity = kind(records, id, "B_SPLINE_CURVE_WITH_KNOTS", budget)?;
    let [
        Value::String(_),
        degree,
        poles,
        Value::Enum("UNSPECIFIED"),
        Value::Enum("F"),
        Value::Enum("F"),
        mult,
        knots,
        Value::Enum("UNSPECIFIED"),
    ] = entity.args.as_slice()
    else {
        return Err(Error::Unsupported);
    };
    let mut numbers = numeric::Numbers::new(entity, budget)?;
    if numbers.integer(degree, budget)? != 1 {
        return Err(Error::Unsupported);
    }
    let mult = list(mult)?;
    let knots = list(knots)?;
    if mult.len() != 2 || knots.len() != 2 {
        return Err(Error::Unsupported);
    }
    for value in mult {
        if numbers.integer(value, budget)? != 2 {
            return Err(Error::Unsupported);
        }
    }
    for (value, expected) in knots.iter().zip([0., 1.]) {
        if !numeric::expected(numbers.scalar(value, budget)?, expected, budget)? {
            return Err(Error::Unsupported);
        }
    }
    numbers.finish()?;
    let poles = list(poles)?;
    if poles.len() != 2 {
        return Err(Error::Unsupported);
    }
    let mut coordinates = [[false; 2]; 2];
    let mut zero = [[false; 2]; 2];
    for (i, pole) in poles.iter().enumerate() {
        let p = point(
            kind(records, reference(pole)?, "CARTESIAN_POINT", budget)?,
            2,
            budget,
        )?;
        for (axis, value) in p.iter().enumerate() {
            coordinates[i][axis] = numeric::expected(*value, 1., budget)?;
            zero[i][axis] = numeric::expected(*value, 0., budget)?;
        }
        budget.release(&p);
        drop(p);
    }
    for (side, expected) in [
        [[0usize, 0], [0, 1]],
        [[0, 1], [1, 1]],
        [[1, 0], [1, 1]],
        [[0, 0], [1, 0]],
    ]
    .into_iter()
    .enumerate()
    {
        if (0..2).all(|i| {
            (0..2).all(|axis| {
                if expected[i][axis] == 0 {
                    zero[i][axis]
                } else {
                    coordinates[i][axis]
                }
            })
        }) {
            return Ok(side);
        }
    }
    Err(Error::Unsupported)
}
fn cap_curve<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    points: &[Point<'input>],
    id: u32,
    curve: &Curve,
    height: usize,
    budget: &mut Budget<'_>,
) -> Result<()> {
    let entity = kind(records, id, "B_SPLINE_CURVE_WITH_KNOTS", budget)?;
    let [
        Value::String(_),
        degree,
        poles,
        Value::Enum("UNSPECIFIED"),
        Value::Enum("F"),
        Value::Enum("F"),
        mult,
        knots,
        Value::Enum("UNSPECIFIED"),
    ] = entity.args.as_slice()
    else {
        return Err(Error::Unsupported);
    };
    let mut numbers = numeric::Numbers::new(entity, budget)?;
    if numbers.integer(degree, budget)? != 11
        || uniform_knots(&mut numbers, mult, knots, budget)? != curve.spans
    {
        return Err(Error::Unsupported);
    }
    numbers.finish()?;
    let poles = list(poles)?;
    if poles.len() != curve.poles.len() {
        return Err(Error::Invalid);
    }
    for (value, index) in poles.iter().zip(&curve.poles) {
        let p = point(
            kind(records, reference(value)?, "CARTESIAN_POINT", budget)?,
            2,
            budget,
        )?;
        for (a, b) in p.iter().zip(points[*index].coordinates) {
            if !numeric::equal(a.spelling(), b.spelling(), budget)? {
                return Err(Error::Unsupported);
            }
        }
        if !numeric::equal(
            points[*index].coordinates[2].spelling(),
            points[height].coordinates[2].spelling(),
            budget,
        )? {
            return Err(Error::Unsupported);
        }
        budget.release(&p);
        drop(p);
    }
    Ok(())
}

#[cfg(test)]
mod tests;
pub(super) struct BindingInput<'a, 'input> {
    pub records: &'a BTreeMap<u32, Vec<Entity<'input>>>,
    pub points: &'a [Point<'input>],
    pub curves: &'a [Curve],
    pub tensors: &'a [Tensor],
    pub edges: &'a [Edge],
}
pub(super) fn face(
    input: BindingInput<'_, '_>,
    face: &Face,
    support_id: u32,
    budget: &mut Budget<'_>,
) -> Result<Option<[usize; 4]>> {
    let BindingInput {
        records,
        points,
        curves,
        tensors,
        edges,
    } = input;
    let mut sides = [0usize; 4];
    for (index, coedge) in face.coedges.iter().enumerate() {
        let edge = &edges[coedge.edge];
        let mut selected = None;
        for id in edge.pcurves {
            let (support, geometry) = pcurve(records, id, budget)?;
            if support == support_id {
                if selected.is_some() {
                    return Err(Error::Invalid);
                }
                selected = Some(geometry);
            }
        }
        let geometry = selected.ok_or(Error::Unsupported)?;
        let curve = &curves[edge.curve];
        match face.support {
            Support::Tensor(tensor) => {
                let side = line(records, geometry, budget)?;
                if sides[..index].contains(&side) || coedge.forward != (side < 2) {
                    return Err(Error::Unsupported);
                }
                sides[index] = side;
                let tensor = &tensors[tensor];
                let dims = tensor.spans.map(|n| 11 * n + 1);
                let length = if side % 2 == 0 { dims[1] } else { dims[0] };
                if curve.poles.len() != length {
                    return Err(Error::Unsupported);
                }
                for (i, &point) in curve.poles.iter().enumerate() {
                    budget.step(1)?;
                    let [u, v] = match side {
                        0 => [0, i],
                        1 => [i, dims[1] - 1],
                        2 => [dims[0] - 1, i],
                        _ => [i, 0],
                    };
                    let surface_point = tensor.poles[u * dims[1] + v];
                    if point != surface_point {
                        // Distinct original owners may define the SAME exact
                        // boundary coefficient. Earn equality from the actual
                        // decimal leaves, never from overlapping float caches.
                        // Keep both point records and all topology identities.
                        for (a, b) in points[point]
                            .coordinates
                            .iter()
                            .zip(points[surface_point].coordinates)
                        {
                            if !numeric::equal(a.spelling(), b.spelling(), budget)? {
                                return Err(Error::Unsupported);
                            }
                        }
                    }
                }
            }
            Support::Plane { height } => {
                cap_curve(records, points, geometry, curve, height, budget)?
            }
        }
    }
    if matches!(face.support, Support::Tensor(_))
        && (0..4).any(|i| sides[(i + 1) % 4] != (sides[i] + 1) % 4)
    {
        return Err(Error::Unsupported);
    }
    Ok(matches!(face.support, Support::Tensor(_)).then_some(sides))
}
