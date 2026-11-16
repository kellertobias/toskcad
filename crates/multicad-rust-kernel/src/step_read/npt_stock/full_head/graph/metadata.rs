//! Typed metadata edges cannot authorize geometry dependencies.
use super::*;
fn target(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    value: &Value<'_>,
    name: &str,
    b: &mut Budget<'_>,
) -> Result<()> {
    kind(records, reference(value)?, name, b)?;
    Ok(())
}

#[cfg(test)]
mod tests;
fn integer(e: &Entity<'_>, v: &Value<'_>, expected: usize, b: &mut Budget<'_>) -> Result<()> {
    let mut n = numeric::Numbers::new(e, b)?;
    if n.integer(v, b)? != expected {
        return Err(Error::Unsupported);
    }
    n.finish()
}
fn strings2(e: &Entity<'_>) -> Result<()> {
    if matches!(e.args.as_slice(), [Value::String(_), Value::String(_)]) {
        Ok(())
    } else {
        Err(Error::Unsupported)
    }
}
fn context(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    value: &Value<'_>,
    dimension: usize,
    b: &mut Budget<'_>,
) -> Result<()> {
    b.step(128)?;
    let parts = records.get(&reference(value)?).ok_or(Error::Invalid)?;
    match parts.as_slice() {
        [g, p, r]
            if dimension == 2
                && g.name == "GEOMETRIC_REPRESENTATION_CONTEXT"
                && p.name == "PARAMETRIC_REPRESENTATION_CONTEXT"
                && r.name == "REPRESENTATION_CONTEXT" =>
        {
            Ok(())
        }
        [g, u, units, r]
            if dimension == 3
                && g.name == "GEOMETRIC_REPRESENTATION_CONTEXT"
                && u.name == "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT"
                && units.name == "GLOBAL_UNIT_ASSIGNED_CONTEXT"
                && r.name == "REPRESENTATION_CONTEXT" =>
        {
            Ok(())
        }
        _ => Err(Error::Unsupported),
    }
}
fn unit(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    v: &Value<'_>,
    role: &str,
    b: &mut Budget<'_>,
) -> Result<()> {
    b.step(128)?;
    let p = records.get(&reference(v)?).ok_or(Error::Invalid)?;
    if p.iter().any(|e| e.name == role) {
        Ok(())
    } else {
        Err(Error::Unsupported)
    }
}
fn identity_axis(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    value: &Value<'_>,
    b: &mut Budget<'_>,
) -> Result<()> {
    let axis = kind(records, reference(value)?, "AXIS2_PLACEMENT_3D", b)?;
    let [Value::String(_), origin, normal, reference_dir] = axis.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    for (index, value) in [origin, normal, reference_dir].into_iter().enumerate() {
        let e = kind(
            records,
            reference(value)?,
            if index == 0 {
                "CARTESIAN_POINT"
            } else {
                "DIRECTION"
            },
            b,
        )?;
        let p = point(e, 3, b)?;
        for (coordinate, scalar) in p.iter().enumerate() {
            b.step(32)?;
            let must_zero = match index {
                0 => true,
                1 => coordinate != 2,
                _ => coordinate == 1,
            };
            if must_zero && !numeric::expected(*scalar, 0., b)? {
                return Err(Error::Unsupported);
            }
            if ((index == 1 && coordinate == 2) || (index == 2 && coordinate == 0))
                && (!scalar.bounds()[0].is_finite() || scalar.bounds()[0] <= 0.)
            {
                return Err(Error::Unsupported);
            }
        }
        b.release(&p);
    }
    // Original CompoundDocument::Shell::placement already checked finite,
    // nondegenerate directions. Exact zero coordinates and positive scales
    // prove +Z normal and +X projected reference algebraically, including a
    // reference's permitted Z component; no rounded normalization is trusted.
    Ok(())
}
pub(super) fn check(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    parts: &[Entity<'_>],
    b: &mut Budget<'_>,
) -> Result<()> {
    b.step(128)?;
    if parts.len() > 1 {
        match parts {
            [g, p, r]
                if g.name == "GEOMETRIC_REPRESENTATION_CONTEXT"
                    && p.name == "PARAMETRIC_REPRESENTATION_CONTEXT"
                    && r.name == "REPRESENTATION_CONTEXT" =>
            {
                let [dimension] = g.args.as_slice() else {
                    return Err(Error::Unsupported);
                };
                integer(g, dimension, 2, b)?;
                if !p.args.is_empty() {
                    return Err(Error::Unsupported);
                }
                strings2(r)?;
            }
            [g, u, units, r]
                if g.name == "GEOMETRIC_REPRESENTATION_CONTEXT"
                    && u.name == "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT"
                    && units.name == "GLOBAL_UNIT_ASSIGNED_CONTEXT"
                    && r.name == "REPRESENTATION_CONTEXT" =>
            {
                let [dimension] = g.args.as_slice() else {
                    return Err(Error::Unsupported);
                };
                integer(g, dimension, 3, b)?;
                strings2(r)?;
                let [uncertainty] = u.args.as_slice() else {
                    return Err(Error::Unsupported);
                };
                let [value] = list(uncertainty)? else {
                    return Err(Error::Unsupported);
                };
                target(records, value, "UNCERTAINTY_MEASURE_WITH_UNIT", b)?;
                let [values] = units.args.as_slice() else {
                    return Err(Error::Unsupported);
                };
                let [length, angle, solid] = list(values)? else {
                    return Err(Error::Unsupported);
                };
                unit(records, length, "LENGTH_UNIT", b)?;
                unit(records, angle, "PLANE_ANGLE_UNIT", b)?;
                unit(records, solid, "SOLID_ANGLE_UNIT", b)?;
            }
            [l, n, s]
                if l.name == "LENGTH_UNIT" && n.name == "NAMED_UNIT" && s.name == "SI_UNIT" =>
            {
                if !l.args.is_empty()
                    || !matches!(n.args.as_slice(), [Value::Derived])
                    || !matches!(
                        s.args.as_slice(),
                        [Value::Enum("MILLI"), Value::Enum("METRE")]
                    )
                {
                    return Err(Error::Unsupported);
                }
            }
            [n, a, s]
                if n.name == "NAMED_UNIT"
                    && a.name == "PLANE_ANGLE_UNIT"
                    && s.name == "SI_UNIT" =>
            {
                if !a.args.is_empty()
                    || !matches!(n.args.as_slice(), [Value::Derived])
                    || !matches!(s.args.as_slice(), [Value::Omitted, Value::Enum("RADIAN")])
                {
                    return Err(Error::Unsupported);
                }
            }
            [n, s, a]
                if n.name == "NAMED_UNIT"
                    && s.name == "SI_UNIT"
                    && a.name == "SOLID_ANGLE_UNIT" =>
            {
                if !a.args.is_empty()
                    || !matches!(n.args.as_slice(), [Value::Derived])
                    || !matches!(
                        s.args.as_slice(),
                        [Value::Omitted, Value::Enum("STERADIAN")]
                    )
                {
                    return Err(Error::Unsupported);
                }
            }
            _ => return Err(Error::Unsupported),
        }
        return Ok(());
    }
    let [e] = parts else {
        return Err(Error::Invalid);
    };
    match e.name {
        "APPLICATION_CONTEXT" => {
            if !matches!(e.args.as_slice(), [Value::String(_)]) {
                return Err(Error::Unsupported);
            }
        }
        "PRODUCT_CONTEXT" | "PRODUCT_DEFINITION_CONTEXT" => {
            let [Value::String(_), context, Value::String(_)] = e.args.as_slice() else {
                return Err(Error::Unsupported);
            };
            target(records, context, "APPLICATION_CONTEXT", b)?;
        }
        "PRODUCT" => {
            let [
                Value::String(_),
                Value::String(_),
                Value::String(_),
                contexts,
            ] = e.args.as_slice()
            else {
                return Err(Error::Unsupported);
            };
            let [c] = list(contexts)? else {
                return Err(Error::Unsupported);
            };
            target(records, c, "PRODUCT_CONTEXT", b)?;
        }
        "PRODUCT_DEFINITION_FORMATION" => {
            let [Value::String(_), Value::String(_), product] = e.args.as_slice() else {
                return Err(Error::Unsupported);
            };
            target(records, product, "PRODUCT", b)?;
        }
        "PRODUCT_DEFINITION" => {
            let [Value::String(_), Value::String(_), formation, context] = e.args.as_slice() else {
                return Err(Error::Unsupported);
            };
            target(records, formation, "PRODUCT_DEFINITION_FORMATION", b)?;
            target(records, context, "PRODUCT_DEFINITION_CONTEXT", b)?;
        }
        "PRODUCT_DEFINITION_SHAPE" => {
            let [Value::String(_), Value::String(_), definition] = e.args.as_slice() else {
                return Err(Error::Unsupported);
            };
            target(records, definition, "PRODUCT_DEFINITION", b)?;
        }
        "SHAPE_DEFINITION_REPRESENTATION" => {
            let [shape, representation] = e.args.as_slice() else {
                return Err(Error::Unsupported);
            };
            target(records, shape, "PRODUCT_DEFINITION_SHAPE", b)?;
            target(
                records,
                representation,
                "ADVANCED_BREP_SHAPE_REPRESENTATION",
                b,
            )?;
        }
        "ADVANCED_BREP_SHAPE_REPRESENTATION" => {
            let [Value::String(_), items, c] = e.args.as_slice() else {
                return Err(Error::Unsupported);
            };
            context(records, c, 3, b)?;
            let items = list(items)?;
            if items.len() != 2 {
                return Err(Error::Unsupported);
            }
            let mut solids = 0;
            let mut axes = 0;
            for item in items {
                b.step(128)?;
                let [target] = records
                    .get(&reference(item)?)
                    .ok_or(Error::Invalid)?
                    .as_slice()
                else {
                    return Err(Error::Unsupported);
                };
                match target.name {
                    "MANIFOLD_SOLID_BREP" => solids += 1,
                    "AXIS2_PLACEMENT_3D" => {
                        axes += 1;
                        identity_axis(records, item, b)?;
                    }
                    _ => return Err(Error::Unsupported),
                }
            }
            if solids != 1 || axes != 1 {
                return Err(Error::Unsupported);
            }
        }
        "DEFINITIONAL_REPRESENTATION" => {
            let [Value::String(_), items, c] = e.args.as_slice() else {
                return Err(Error::Unsupported);
            };
            context(records, c, 2, b)?;
            let [curve] = list(items)? else {
                return Err(Error::Unsupported);
            };
            target(records, curve, "B_SPLINE_CURVE_WITH_KNOTS", b)?;
        }
        "UNCERTAINTY_MEASURE_WITH_UNIT" => {
            let [
                Value::Typed(measure),
                length,
                Value::String(_),
                Value::String(_),
            ] = e.args.as_slice()
            else {
                return Err(Error::Unsupported);
            };
            if measure.name != "LENGTH_MEASURE"
                || !matches!(measure.args.as_slice(),[Value::Number(n)] if n.is_finite()&&*n>0.)
            {
                return Err(Error::Unsupported);
            }
            unit(records, length, "LENGTH_UNIT", b)?;
        }
        "APPLICATION_PROTOCOL_DEFINITION" => {
            let [
                Value::String("international standard"),
                Value::String("automotive_design"),
                year,
                c,
            ] = e.args.as_slice()
            else {
                return Err(Error::Unsupported);
            };
            integer(e, year, 2000, b)?;
            target(records, c, "APPLICATION_CONTEXT", b)?;
        }
        "PRODUCT_RELATED_PRODUCT_CATEGORY" => {
            let [Value::String("part"), Value::Omitted, products] = e.args.as_slice() else {
                return Err(Error::Unsupported);
            };
            let [p] = list(products)? else {
                return Err(Error::Unsupported);
            };
            target(records, p, "PRODUCT", b)?;
        }
        "CARTESIAN_POINT"
        | "PLANE"
        | "DIRECTION"
        | "AXIS2_PLACEMENT_3D"
        | "VECTOR"
        | "LINE"
        | "B_SPLINE_CURVE_WITH_KNOTS"
        | "B_SPLINE_SURFACE_WITH_KNOTS"
        | "SURFACE_CURVE"
        | "PCURVE"
        | "VERTEX_POINT"
        | "EDGE_CURVE"
        | "ORIENTED_EDGE"
        | "EDGE_LOOP"
        | "FACE_BOUND"
        | "FACE_OUTER_BOUND"
        | "ADVANCED_FACE"
        | "CLOSED_SHELL"
        | "MANIFOLD_SOLID_BREP" => {}
        _ => return Err(Error::Unsupported),
    }
    Ok(())
}
