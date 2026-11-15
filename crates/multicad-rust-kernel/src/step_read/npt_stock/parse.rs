use super::*;
use records::*;
mod bindings;
mod topology;
fn points<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    budget: &mut Budget<'_>,
) -> Result<Vec<Point<'input>>> {
    budget.step(records.len())?;
    let count=records.values().filter(|values|matches!(values.as_slice(),[Entity {name:"CARTESIAN_POINT",args,..}]if matches!(args.as_slice(),[_,Value::List(p)]if p.len()==3))).count();
    if count > 262_144 + 4096 {
        return Err(Error::Limit);
    }
    let mut points = budget.allocate(count)?;
    for (&id, values) in records {
        let [entity] = values.as_slice() else {
            continue;
        };
        if entity.name != "CARTESIAN_POINT"
            || !matches!(entity.args.as_slice(),[_,Value::List(p)]if p.len()==3)
        {
            continue;
        }
        let coordinates = point(entity, 3, budget)?;
        let fixed = coordinates
            .as_slice()
            .try_into()
            .expect("three actual coordinates");
        budget.release(&coordinates);
        drop(coordinates);
        points.push(Point {
            id,
            coordinates: fixed,
        });
    }
    Ok(points)
}
fn curve<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    points: &[Point<'input>],
    id: u32,
    budget: &mut Budget<'_>,
) -> Result<Curve> {
    let entity = kind(records, id, "B_SPLINE_CURVE_WITH_KNOTS", budget)?;
    let [
        Value::String(_),
        degree,
        poles,
        Value::Enum("UNSPECIFIED"),
        Value::Enum("F"),
        Value::Enum("F"),
        mult,
        values,
        Value::Enum("UNSPECIFIED"),
    ] = entity.args.as_slice()
    else {
        return Err(Error::Unsupported);
    };
    let mut numbers = numeric::Numbers::new(entity, budget)?;
    if numbers.integer(degree, budget)? != 11 {
        return Err(Error::Unsupported);
    }
    let spans = uniform_knots(&mut numbers, mult, values, budget)?;
    numbers.finish()?;
    let refs = list(poles)?;
    if refs.len() != 11 * spans + 1 {
        return Err(Error::Invalid);
    }
    let mut poles = budget.allocate(refs.len())?;
    for value in refs {
        poles.push(point_index(points, reference(value)?, budget)?);
    }
    Ok(Curve { id, spans, poles })
}
fn tensor<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    points: &[Point<'input>],
    id: u32,
    budget: &mut Budget<'_>,
    total: &mut usize,
) -> Result<Tensor> {
    let entity = kind(records, id, "B_SPLINE_SURFACE_WITH_KNOTS", budget)?;
    let [
        Value::String(_),
        du,
        dv,
        rows,
        Value::Enum("UNSPECIFIED"),
        Value::Enum("F"),
        Value::Enum("F"),
        Value::Enum("F"),
        um,
        vm,
        uk,
        vk,
        Value::Enum("UNSPECIFIED"),
    ] = entity.args.as_slice()
    else {
        return Err(Error::Unsupported);
    };
    let mut numbers = numeric::Numbers::new(entity, budget)?;
    if numbers.integer(du, budget)? != 11 || numbers.integer(dv, budget)? != 11 {
        return Err(Error::Unsupported);
    }
    let spans = [
        multiplicities(&mut numbers, um, budget)?,
        multiplicities(&mut numbers, vm, budget)?,
    ];
    knots(&mut numbers, uk, spans[0], budget)?;
    knots(&mut numbers, vk, spans[1], budget)?;
    numbers.finish()?;
    let dims = spans.map(|n| 11 * n + 1);
    let rows = list(rows)?;
    if rows.len() != dims[0] {
        return Err(Error::Invalid);
    }
    *total += dims[0] * dims[1];
    if *total > 262_144 {
        return Err(Error::Limit);
    }
    let mut poles = budget.allocate(dims[0] * dims[1])?;
    for row in rows {
        let row = list(row)?;
        if row.len() != dims[1] {
            return Err(Error::Invalid);
        }
        for value in row {
            poles.push(point_index(points, reference(value)?, budget)?);
        }
    }
    Ok(Tensor { id, spans, poles })
}
fn plane<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    points: &[Point<'input>],
    id: u32,
    budget: &mut Budget<'_>,
) -> Result<Support> {
    let entity = kind(records, id, "PLANE", budget)?;
    let [Value::String(_), axis] = entity.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let axis = kind(records, reference(axis)?, "AXIS2_PLACEMENT_3D", budget)?;
    let [Value::String(_), origin, normal, tangent] = axis.args.as_slice() else {
        return Err(Error::Invalid);
    };
    for (direction, expected) in [(normal, [0., 0., 1.]), (tangent, [1., 0., 0.])] {
        let entity = kind(records, reference(direction)?, "DIRECTION", budget)?;
        let [Value::String(_), values] = entity.args.as_slice() else {
            return Err(Error::Invalid);
        };
        let values = list(values)?;
        if values.len() != 3 {
            return Err(Error::Invalid);
        }
        let mut numbers = numeric::Numbers::new(entity, budget)?;
        for (value, expected) in values.iter().zip(expected) {
            if !numeric::expected(numbers.scalar(value, budget)?, expected, budget)? {
                return Err(Error::Unsupported);
            }
        }
        numbers.finish()?;
    }
    let height = point_index(points, reference(origin)?, budget)?;
    for scalar in &points[height].coordinates[..2] {
        if !numeric::expected(*scalar, 0., budget)? {
            return Err(Error::Unsupported);
        }
    }
    Ok(Support::Plane { height })
}
pub(super) fn body<'document, 'input>(
    document: &'document CompoundDocument<'input>,
    solid: u32,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Parsed<'document, 'input>> {
    let mut budget = Budget::new(options, interrupt)?;
    budget.step(1)?;
    let [member] = document.represented_members() else {
        return Err(Error::Unsupported);
    };
    if member.item_entity() != solid || member.mapping_link().is_some() {
        return Err(Error::Unsupported);
    }
    let records = document.records();
    let entity = kind(records, solid, "MANIFOLD_SOLID_BREP", &mut budget)?;
    let [Value::String(_), shell] = entity.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let shell = kind(records, reference(shell)?, "CLOSED_SHELL", &mut budget)?;
    let [Value::String(_), faces] = shell.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let face_ids = list(faces)?;
    if face_ids.len() != 10 {
        return Err(Error::Unsupported);
    }
    let points = points(records, &mut budget)?;
    let mut tensors = budget.allocate(8)?;
    let mut vertices = budget.allocate(12)?;
    let mut curves = budget.allocate(20)?;
    let mut edges = budget.allocate(20)?;
    let mut faces = budget.allocate(10)?;
    let mut total = 0;
    for value in face_ids {
        let id = reference(value)?;
        if faces.iter().any(|face: &Face| face.id == id) {
            return Err(Error::Invalid);
        }
        let entity = kind(records, id, "ADVANCED_FACE", &mut budget)?;
        let [Value::String(_), bounds, surface, forward] = entity.args.as_slice() else {
            return Err(Error::Invalid);
        };
        let bounds = list(bounds)?;
        let [bound] = bounds else {
            return Err(Error::Unsupported);
        };
        let forward = boolean(forward)?;
        let surface = reference(surface)?;
        let support = match records.get(&surface).map(Vec::as_slice) {
            Some(
                [
                    Entity {
                        name: "B_SPLINE_SURFACE_WITH_KNOTS",
                        ..
                    },
                ],
            ) => {
                if forward || tensors.len() == 8 {
                    return Err(Error::Unsupported);
                }
                if tensors.iter().any(|t: &Tensor| t.id == surface) {
                    return Err(Error::Unsupported);
                }
                let value = tensor(records, &points, surface, &mut budget, &mut total)?;
                let index = tensors.len();
                tensors.push(value);
                Support::Tensor(index)
            }
            Some([Entity { name: "PLANE", .. }]) => plane(records, &points, surface, &mut budget)?,
            _ => return Err(Error::Unsupported),
        };
        let coedges = topology::wire(
            records,
            reference(bound)?,
            &points,
            &mut vertices,
            &mut curves,
            &mut edges,
            &mut budget,
        )?;
        let face = Face {
            id,
            support,
            coedges,
            forward,
        };
        bindings::face(
            bindings::BindingInput {
                records,
                points: &points,
                curves: &curves,
                tensors: &tensors,
                edges: &edges,
            },
            &face,
            surface,
            &mut budget,
        )?;
        faces.push(face);
    }
    topology::shell(&vertices, &edges, &faces, &mut budget)?;
    if vertices.len() != 12 || edges.len() != 20 || curves.len() != 20 || tensors.len() != 8 {
        return Err(Error::Unsupported);
    }
    budget.step(1)?;
    Ok(Parsed {
        document,
        solid,
        points,
        vertices,
        curves,
        tensors,
        edges,
        faces,
        work: budget.used,
        bytes: budget.peak,
    })
}
