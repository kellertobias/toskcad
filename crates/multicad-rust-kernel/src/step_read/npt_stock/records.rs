use super::*;
pub(super) fn kind<'arena, 'input>(
    records: &'arena BTreeMap<u32, Vec<Entity<'input>>>,
    id: u32,
    name: &str,
    budget: &mut Budget<'_>,
) -> Result<&'arena Entity<'input>> {
    budget.step(32)?;
    let records = records.get(&id).ok_or(Error::Invalid)?;
    match records.as_slice() {
        [entity] if entity.name == name => Ok(entity),
        _ => Err(Error::Unsupported),
    }
}
pub(super) fn reference(value: &Value<'_>) -> Result<u32> {
    if let Value::Reference(id) = value {
        Ok(*id)
    } else {
        Err(Error::Invalid)
    }
}
pub(super) fn list<'a, 'input>(value: &'a Value<'input>) -> Result<&'a [Value<'input>]> {
    if let Value::List(values) = value {
        Ok(values)
    } else {
        Err(Error::Invalid)
    }
}
pub(super) fn boolean(value: &Value<'_>) -> Result<bool> {
    match value {
        Value::Enum("T") => Ok(true),
        Value::Enum("F") => Ok(false),
        _ => Err(Error::Invalid),
    }
}
pub(super) fn point<'input>(
    entity: &Entity<'input>,
    dimensions: usize,
    budget: &mut Budget<'_>,
) -> Result<Vec<Scalar<'input>>> {
    let [Value::String(_), values] = entity.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let values = list(values)?;
    if values.len() != dimensions {
        return Err(Error::Unsupported);
    }
    let mut numbers = numeric::Numbers::new(entity, budget)?;
    let mut point = budget.allocate(dimensions)?;
    for value in values {
        let scalar = numbers.scalar(value, budget)?;
        if scalar
            .bounds()
            .into_iter()
            .any(|x| !x.is_finite() || x.abs() > 1e7)
        {
            return Err(Error::Limit);
        }
        point.push(scalar);
    }
    numbers.finish()?;
    Ok(point)
}
pub(super) fn point_index(points: &[Point<'_>], id: u32, budget: &mut Budget<'_>) -> Result<usize> {
    budget.step(32)?;
    points
        .binary_search_by_key(&id, |p| p.id)
        .map_err(|_| Error::Invalid)
}
pub(super) fn uniform_knots(
    numbers: &mut numeric::Numbers<'_>,
    multiplicities: &Value<'_>,
    knots: &Value<'_>,
    budget: &mut Budget<'_>,
) -> Result<usize> {
    let mult = list(multiplicities)?;
    let knots = list(knots)?;
    if mult.len() < 2 || mult.len() != knots.len() {
        return Err(Error::Invalid);
    }
    let spans = knots.len() - 1;
    if !spans.is_power_of_two() || spans > 256 {
        return Err(Error::Unsupported);
    }
    for (index, value) in mult.iter().enumerate() {
        if numbers.integer(value, budget)? != if index == 0 || index == spans { 12 } else { 11 } {
            return Err(Error::Unsupported);
        }
    }
    for (index, value) in knots.iter().enumerate() {
        if !numeric::expected(
            numbers.scalar(value, budget)?,
            index as f64 / spans as f64,
            budget,
        )? {
            return Err(Error::Unsupported);
        }
    }
    Ok(spans)
}
pub(super) fn multiplicities(
    numbers: &mut numeric::Numbers<'_>,
    values: &Value<'_>,
    budget: &mut Budget<'_>,
) -> Result<usize> {
    let values = list(values)?;
    if values.len() < 2 {
        return Err(Error::Invalid);
    }
    let spans = values.len() - 1;
    if !spans.is_power_of_two() || spans > 256 {
        return Err(Error::Unsupported);
    }
    for (index, value) in values.iter().enumerate() {
        if numbers.integer(value, budget)? != if index == 0 || index == spans { 12 } else { 11 } {
            return Err(Error::Unsupported);
        }
    }
    Ok(spans)
}
pub(super) fn knots(
    numbers: &mut numeric::Numbers<'_>,
    values: &Value<'_>,
    spans: usize,
    budget: &mut Budget<'_>,
) -> Result<()> {
    let values = list(values)?;
    if values.len() != spans + 1 {
        return Err(Error::Invalid);
    }
    for (index, value) in values.iter().enumerate() {
        if !numeric::expected(
            numbers.scalar(value, budget)?,
            index as f64 / spans as f64,
            budget,
        )? {
            return Err(Error::Unsupported);
        }
    }
    Ok(())
}
