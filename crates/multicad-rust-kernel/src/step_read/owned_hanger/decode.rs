use super::*;
type Arena<'i> = BTreeMap<u32, Vec<Entity<'i>>>;
fn one<'a, 'i>(arena: &'a Arena<'i>, id: u32, name: &str) -> Result<&'a Entity<'i>, Error> {
    let e = arena.get(&id).ok_or(Error::Invalid)?;
    if e.len() != 1 || e[0].name != name {
        return Err(Error::Unsupported);
    }
    Ok(&e[0])
}
fn component<'a, 'i>(entities: &'a [Entity<'i>], name: &str) -> Result<&'a Entity<'i>, Error> {
    // Every arena comes directly from part21::parse, whose record parser
    // rejects duplicate component names. Do not rescan each sealed record.
    entities
        .iter()
        .find(|e| e.name == name)
        .ok_or(Error::Unsupported)
}
fn components(entities: &[Entity<'_>], names: &[&str]) -> Result<(), Error> {
    if entities.len() != names.len() {
        return Err(Error::Unsupported);
    }
    for name in names {
        component(entities, name)?;
    }
    Ok(())
}
fn number_list<'a, 'i>(v: &'a Value<'i>) -> Result<&'a [Value<'i>], Error> {
    let values = list(v)?;
    if values.iter().any(|v| !matches!(v, Value::Number(_))) {
        return Err(Error::Invalid);
    }
    Ok(values)
}
fn reference(v: &Value<'_>) -> Result<u32, Error> {
    match v {
        Value::Reference(id) => Ok(*id),
        _ => Err(Error::Invalid),
    }
}
fn list<'a, 'i>(v: &'a Value<'i>) -> Result<&'a [Value<'i>], Error> {
    match v {
        Value::List(v) => Ok(v),
        _ => Err(Error::Invalid),
    }
}
fn refs(v: &Value<'_>) -> Result<Vec<u32>, Error> {
    list(v)?.iter().map(reference).collect()
}
fn boolean(v: &Value<'_>) -> Result<bool, Error> {
    match v {
        Value::Enum("T") => Ok(true),
        Value::Enum("F") => Ok(false),
        _ => Err(Error::Invalid),
    }
}
fn nat(s: Scalar<'_>) -> Result<usize, Error> {
    s.text.parse().map_err(|_| Error::Invalid)
}
fn numeric<'i>(e: &Entity<'i>, b: &mut Budget<'_>) -> Result<Vec<Scalar<'i>>, Error> {
    b.charge(e.source.ok_or(Error::Invalid)?.len())?;
    let tokens = part21::numeric_spellings(e).ok_or(Error::Invalid)?;
    let out = tokens
        .into_iter()
        .map(|text| {
            let value = text.parse::<f64>().map_err(|_| Error::Invalid)?;
            if !value.is_finite() {
                return Err(Error::Invalid);
            }
            b.charge(text.len() * (1 + usize::from(value == 0. || value == 1.)))?;
            Ok(Scalar::new(text, value))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(out)
}
fn sizes(values: &[Value<'_>]) -> Result<usize, Error> {
    let mut total = 0_usize;
    for value in values {
        if let Value::List(v) = value {
            total = total
                .checked_add(v.capacity() * size_of::<Value<'_>>() + sizes(v)?)
                .ok_or(Error::Limit)?;
        }
        if let Value::Typed(e) = value {
            total = total
                .checked_add(e.args.capacity() * size_of::<Value<'_>>() + sizes(&e.args)?)
                .ok_or(Error::Limit)?;
        }
    }
    Ok(total)
}
fn validate_knots(values: &[Scalar<'_>]) -> Result<(), Error> {
    if values.len() < 2
        || values
            .windows(2)
            .any(|p| p[0].interval().hi >= p[1].interval().lo)
    {
        return Err(Error::Invalid);
    }
    Ok(())
}
fn point<'i>(id: u32, e: &Entity<'i>, b: &mut Budget<'_>) -> Result<Point<'i>, Error> {
    if e.args.len() != 2 {
        return Err(Error::Invalid);
    }
    let dim = number_list(&e.args[1])?.len();
    if !matches!(dim, 2 | 3) {
        return Err(Error::Unsupported);
    }
    let coordinates = numeric(e, b)?;
    if coordinates.len() != dim {
        return Err(Error::Invalid);
    }
    b.retain(128 + coordinates.capacity() * size_of::<Scalar<'_>>())?;
    Ok(Point {
        entity: id,
        coordinates,
    })
}
fn curve<'i>(id: u32, entities: &[Entity<'i>], b: &mut Budget<'_>) -> Result<Curve<'i>, Error> {
    let simple = entities.len() == 1;
    if !simple {
        components(
            entities,
            &[
                "BOUNDED_CURVE",
                "B_SPLINE_CURVE",
                "B_SPLINE_CURVE_WITH_KNOTS",
                "CURVE",
                "GEOMETRIC_REPRESENTATION_ITEM",
                "RATIONAL_B_SPLINE_CURVE",
                "REPRESENTATION_ITEM",
            ],
        )?;
    }
    let base = component(
        entities,
        if simple {
            "B_SPLINE_CURVE_WITH_KNOTS"
        } else {
            "B_SPLINE_CURVE"
        },
    )?;
    if base.args.len() != if simple { 9 } else { 5 } {
        return Err(Error::Invalid);
    }
    let tail = if simple { 3 } else { 2 };
    if !matches!(
        &base.args[tail..tail + 3],
        [
            Value::Enum("UNSPECIFIED"),
            Value::Enum("F"),
            Value::Enum("F")
        ]
    ) {
        return Err(Error::Unsupported);
    }
    let numbers = numeric(base, b)?;
    let degree = nat(*numbers.first().ok_or(Error::Invalid)?)?;
    if degree == 0 || degree > 32 {
        return Err(Error::Unsupported);
    }
    let points = refs(
        base.args
            .get(if simple { 2 } else { 1 })
            .ok_or(Error::Invalid)?,
    )?;
    let with = component(entities, "B_SPLINE_CURVE_WITH_KNOTS")?;
    let offset = if simple { 6 } else { 0 };
    let multis = list(with.args.get(offset).ok_or(Error::Invalid)?)?;
    let nk = list(with.args.get(offset + 1).ok_or(Error::Invalid)?)?.len();
    if multis.len() != nk {
        return Err(Error::Invalid);
    }
    let all = numeric(with, b)?;
    let skip = usize::from(simple);
    if all.len() != skip + 2 * nk {
        return Err(Error::Invalid);
    }
    for (i, &m) in all[skip..skip + nk].iter().enumerate() {
        if nat(m)? != degree + usize::from(i == 0 || i + 1 == nk) {
            return Err(Error::Unsupported);
        }
    }
    let knots = all[skip + nk..].to_vec();
    validate_knots(&knots)?;
    if points.len() != degree.checked_mul(nk - 1).ok_or(Error::Limit)? + 1 {
        return Err(Error::Invalid);
    }
    let weights = if entities.iter().any(|e| e.name == "RATIONAL_B_SPLINE_CURVE") {
        let rational = component(entities, "RATIONAL_B_SPLINE_CURVE")?;
        if rational.args.len() != 1 || number_list(&rational.args[0])?.len() != points.len() {
            return Err(Error::Invalid);
        }
        let ws = numeric(rational, b)?;
        if ws.len() != points.len() || ws.iter().any(|w| w.interval().lo <= 0.) {
            return Err(Error::Invalid);
        }
        ws
    } else {
        Vec::new()
    };
    b.retain(
        256 + points.capacity() * size_of::<u32>()
            + (weights.capacity() + knots.capacity()) * size_of::<Scalar<'_>>(),
    )?;
    Ok(Curve {
        entity: id,
        degree,
        points,
        weights,
        knots,
    })
}
fn tensor<'i>(id: u32, entities: &[Entity<'i>], b: &mut Budget<'_>) -> Result<Tensor<'i>, Error> {
    components(
        entities,
        &[
            "BOUNDED_SURFACE",
            "B_SPLINE_SURFACE",
            "B_SPLINE_SURFACE_WITH_KNOTS",
            "GEOMETRIC_REPRESENTATION_ITEM",
            "RATIONAL_B_SPLINE_SURFACE",
            "REPRESENTATION_ITEM",
            "SURFACE",
        ],
    )?;
    let base = component(entities, "B_SPLINE_SURFACE")?;
    if base.args.len() != 7
        || !matches!(
            &base.args[3..],
            [
                Value::Enum("UNSPECIFIED"),
                Value::Enum("F"),
                Value::Enum("F"),
                Value::Enum("F")
            ]
        )
    {
        return Err(Error::Unsupported);
    }
    let nums = numeric(base, b)?;
    if nums.len() != 2 {
        return Err(Error::Invalid);
    }
    let degrees = [nat(nums[0])?, nat(nums[1])?];
    if !matches!(degrees, [1 | 2, 10]) {
        return Err(Error::Unsupported);
    }
    let rows = list(base.args.get(2).ok_or(Error::Invalid)?)?;
    if rows.len() != degrees[0] + 1 {
        return Err(Error::Invalid);
    }
    let dimensions = [rows.len(), list(rows.first().ok_or(Error::Invalid)?)?.len()];
    let mut points = Vec::with_capacity(dimensions[0] * dimensions[1]);
    for row in rows {
        let row = refs(row)?;
        if row.len() != dimensions[1] {
            return Err(Error::Invalid);
        }
        points.extend(row);
    }
    let with = component(entities, "B_SPLINE_SURFACE_WITH_KNOTS")?;
    if with.args.len() != 5 {
        return Err(Error::Invalid);
    }
    let nu = list(&with.args[0])?.len();
    let nv = list(&with.args[1])?.len();
    if nu != 2 || list(&with.args[2])?.len() != 2 || list(&with.args[3])?.len() != nv {
        return Err(Error::Unsupported);
    }
    let nums = numeric(with, b)?;
    if nums.len() != nu + 2 * nv + 2 {
        return Err(Error::Invalid);
    }
    if nat(nums[0])? != degrees[0] + 1 || nat(nums[1])? != degrees[0] + 1 {
        return Err(Error::Unsupported);
    }
    for (i, &m) in nums[nu..nu + nv].iter().enumerate() {
        if nat(m)? != 10 + usize::from(i == 0 || i + 1 == nv) {
            return Err(Error::Unsupported);
        }
    }
    if nums[nu + nv].text != "0." || nums[nu + nv + 1].text != "1." {
        return Err(Error::Unsupported);
    }
    let knots_v = nums[nu + nv + 2..].to_vec();
    validate_knots(&knots_v)?;
    if dimensions[1] != 10 * (nv - 1) + 1 {
        return Err(Error::Invalid);
    }
    let rational = component(entities, "RATIONAL_B_SPLINE_SURFACE")?;
    let weights = numeric(rational, b)?;
    if weights.len() != points.len() || weights.iter().any(|w| w.interval().lo <= 0.) {
        return Err(Error::Invalid);
    }
    let weight_rows = list(rational.args.first().ok_or(Error::Invalid)?)?;
    if weight_rows.len() != dimensions[0]
        || weight_rows.iter().any(|row| {
            list(row).map_or(true, |r| {
                r.len() != dimensions[1] || r.iter().any(|v| !matches!(v, Value::Number(_)))
            })
        })
    {
        return Err(Error::Invalid);
    }
    b.retain(
        size_of::<Tensor<'_>>()
            + points.capacity() * size_of::<u32>()
            + (weights.capacity() + knots_v.capacity()) * size_of::<Scalar<'_>>(),
    )?;
    Ok(Tensor {
        entity: id,
        degrees,
        dimensions,
        points,
        weights,
        knots_v,
    })
}
fn plane(arena: &Arena<'_>, id: u32) -> Result<Plane, Error> {
    let e = one(arena, id, "PLANE")?;
    let axis = reference(e.args.get(1).ok_or(Error::Invalid)?)?;
    let a = one(arena, axis, "AXIS2_PLACEMENT_3D")?;
    if a.args.len() != 4 {
        return Err(Error::Invalid);
    }
    Ok(Plane {
        entity: id,
        axis,
        origin: reference(&a.args[1])?,
        normal: reference(&a.args[2])?,
        reference: reference(&a.args[3])?,
    })
}
fn direction(
    arena: &Arena<'_>,
    id: u32,
    expected: [bool; 3],
    b: &mut Budget<'_>,
) -> Result<(), Error> {
    let d = one(arena, id, "DIRECTION")?;
    if d.args.len() != 2 || number_list(&d.args[1])?.len() != 3 {
        return Err(Error::Invalid);
    }
    let coordinates = numeric(d, b)?;
    if coordinates.len() != 3
        || coordinates
            .iter()
            .zip(expected)
            .any(|(v, one)| if one { !v.is_one() } else { !v.is_zero() })
    {
        return Err(Error::Unsupported);
    }
    Ok(())
}
fn active(arena: &Arena<'_>, b: &mut Budget<'_>) -> Result<(u32, u32, Vec<u32>), Error> {
    let mut shape = arena
        .values()
        .flat_map(|v| v.iter())
        .filter(|e| e.name == "SHAPE_DEFINITION_REPRESENTATION");
    let link = shape.next().ok_or(Error::Unsupported)?;
    if shape.next().is_some() || link.args.len() != 2 {
        return Err(Error::Unsupported);
    }
    let definition = one(arena, reference(&link.args[0])?, "PRODUCT_DEFINITION_SHAPE")?;
    if definition.args.len() != 3 {
        return Err(Error::Invalid);
    }
    one(arena, reference(&definition.args[2])?, "PRODUCT_DEFINITION")?;
    let rep = one(
        arena,
        reference(&link.args[1])?,
        "ADVANCED_BREP_SHAPE_REPRESENTATION",
    )?;
    if rep.args.len() != 3 {
        return Err(Error::Invalid);
    }
    let items = refs(&rep.args[1])?;
    let mut solids = Vec::new();
    for id in items {
        let entity = arena.get(&id).ok_or(Error::Invalid)?;
        if entity.len() == 1 && entity[0].name == "MANIFOLD_SOLID_BREP" {
            solids.push(id);
            continue;
        }
        let axis = one(arena, id, "AXIS2_PLACEMENT_3D")?;
        if axis.args.len() != 4 {
            return Err(Error::Unsupported);
        }
        let p = one(arena, reference(&axis.args[1])?, "CARTESIAN_POINT")?;
        let coordinates = numeric(p, b)?;
        if coordinates.len() != 3 || coordinates.iter().any(|v| !v.is_zero()) {
            return Err(Error::Unsupported);
        }
        for (position, expected) in [(2, [false, false, true]), (3, [true, false, false])] {
            direction(arena, reference(&axis.args[position])?, expected, b)?;
        }
    }
    if solids.len() != 1 {
        return Err(Error::Unsupported);
    }
    let context = arena.get(&reference(&rep.args[2])?).ok_or(Error::Invalid)?;
    let dim = numeric(component(context, "GEOMETRIC_REPRESENTATION_CONTEXT")?, b)?;
    if dim.len() != 1 || nat(dim[0])? != 3 {
        return Err(Error::Unsupported);
    }
    let units = component(context, "GLOBAL_UNIT_ASSIGNED_CONTEXT")?;
    let ids = refs(units.args.first().ok_or(Error::Invalid)?)?;
    let mut millimetre = 0;
    for id in ids {
        let unit = arena.get(&id).ok_or(Error::Invalid)?;
        if component(unit, "LENGTH_UNIT").is_ok() {
            let si = component(unit, "SI_UNIT")?;
            if !matches!(
                si.args.as_slice(),
                [Value::Enum("MILLI"), Value::Enum("METRE")]
            ) {
                return Err(Error::Unsupported);
            }
            millimetre += 1;
        }
    }
    if millimetre != 1 {
        return Err(Error::Invalid);
    }
    let solid = solids[0];
    let body = one(arena, solid, "MANIFOLD_SOLID_BREP")?;
    let shell = reference(body.args.get(1).ok_or(Error::Invalid)?)?;
    let body = one(arena, shell, "CLOSED_SHELL")?;
    let faces = refs(body.args.get(1).ok_or(Error::Invalid)?)?;
    if faces.is_empty() {
        return Err(Error::Invalid);
    }
    Ok((solid, shell, faces))
}
fn pcurve(arena: &Arena<'_>, id: u32, b: &mut Budget<'_>) -> Result<Pcurve, Error> {
    let e = one(arena, id, "PCURVE")?;
    if e.args.len() != 3 {
        return Err(Error::Invalid);
    }
    let surface = reference(&e.args[1])?;
    let rep = one(arena, reference(&e.args[2])?, "DEFINITIONAL_REPRESENTATION")?;
    if rep.args.len() != 3 {
        return Err(Error::Invalid);
    }
    let curves = refs(&rep.args[1])?;
    if curves.len() != 1 {
        return Err(Error::Unsupported);
    }
    let context = arena.get(&reference(&rep.args[2])?).ok_or(Error::Invalid)?;
    component(context, "PARAMETRIC_REPRESENTATION_CONTEXT")?;
    let dimension = component(context, "GEOMETRIC_REPRESENTATION_CONTEXT")?;
    let dim = numeric(dimension, b)?;
    if dim.len() != 1 || nat(dim[0])? != 2 {
        return Err(Error::Unsupported);
    }
    Ok(Pcurve {
        surface,
        curve: curves[0],
    })
}
fn edge(arena: &Arena<'_>, id: u32, b: &mut Budget<'_>) -> Result<Edge, Error> {
    let e = one(arena, id, "EDGE_CURVE")?;
    if e.args.len() != 5 {
        return Err(Error::Invalid);
    }
    let vertices = [reference(&e.args[1])?, reference(&e.args[2])?];
    if vertices[0] == vertices[1] {
        return Err(Error::Unsupported);
    }
    let s = one(arena, reference(&e.args[3])?, "SURFACE_CURVE")?;
    if s.args.len() != 4 || !matches!(s.args[3], Value::Enum("CURVE_3D")) {
        return Err(Error::Unsupported);
    }
    let ids = refs(&s.args[2])?;
    if ids.len() != 2 {
        return Err(Error::Unsupported);
    }
    Ok(Edge {
        entity: id,
        vertices,
        curve: reference(&s.args[1])?,
        pcurves: [pcurve(arena, ids[0], b)?, pcurve(arena, ids[1], b)?],
        same_sense: boolean(&e.args[4])?,
    })
}
pub(super) fn read<'i>(
    input: &'i str,
    options: ReadOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Decoded<'i>, Error> {
    if options.maximum_work == 0
        || options.maximum_work > 256_000_000
        || options.retained_budget_bytes == 0
        || options.retained_budget_bytes > 64 * 1024 * 1024
    {
        return Err(Error::Invalid);
    }
    let mut b = Budget {
        work: 0,
        remaining: options.maximum_work,
        bytes: 0,
        limit: options.retained_budget_bytes,
        interrupt,
    };
    b.charge(input.len())?;
    b.retain(input.len() + size_of::<Decoded<'_>>())?;
    let arena = part21::parse(input).map_err(|e| match e {
        super::super::ReadError::Invalid => Error::Invalid,
        super::super::ReadError::Limit => Error::Limit,
        _ => Error::Unsupported,
    })?;
    let before_arena = b.bytes;
    b.retain(448)?;
    // Pinned BTree B=6: u32 + Vec24 nodes fit448 bytes, Point32
    // nodes fit576, Curve88 nodes fit1280 (including12 internal links).
    // Every non-root node owns >=5 keys: nodes<=1+(N-1)/5.
    // Charge root plus96/128/256 per AST/point/curve entry respectively.
    for entities in arena.values() {
        b.charge(entities.len())?;
        b.retain(96 + entities.capacity() * size_of::<Entity<'_>>())?;
        for e in entities {
            b.retain(e.args.capacity() * size_of::<Value<'_>>() + sizes(&e.args)?)?;
        }
    }
    let arena_bytes = b.bytes - before_arena;
    let (solid, shell, face_ids) = active(&arena, &mut b)?;
    b.retain(576 + 1280 + 384 + 192 + 448 + 224 + 576)?;
    let mut points = BTreeMap::new();
    let mut curves = BTreeMap::new();
    let mut tensors = Vec::new();
    let mut planes = Vec::new();
    let mut surface_ids = BTreeMap::new();
    for (&id, entities) in &arena {
        b.charge(1)?;
        if entities.len() == 1 && entities[0].name == "CARTESIAN_POINT" {
            points.insert(id, point(id, &entities[0], &mut b)?);
        }
        if entities
            .iter()
            .any(|e| e.name == "B_SPLINE_CURVE" || e.name == "B_SPLINE_CURVE_WITH_KNOTS")
        {
            curves.insert(id, curve(id, entities, &mut b)?);
        }
        if entities.iter().any(|e| e.name == "B_SPLINE_SURFACE") {
            surface_ids.insert(id, Surface::Tensor(tensors.len()));
            tensors.push(tensor(id, entities, &mut b)?);
        }
        if entities.len() == 1 && entities[0].name == "PLANE" {
            surface_ids.insert(id, Surface::Plane(planes.len()));
            planes.push(plane(&arena, id)?);
        }
    }
    for curve in curves.values() {
        let dim = points
            .get(&curve.points[0])
            .ok_or(Error::Invalid)?
            .coordinates
            .len();
        for id in &curve.points {
            b.charge(1)?;
            if points.get(id).is_none_or(|p| p.coordinates.len() != dim) {
                return Err(Error::Invalid);
            }
        }
    }
    for tensor in &tensors {
        for id in &tensor.points {
            b.charge(1)?;
            if points.get(id).is_none_or(|p| p.coordinates.len() != 3) {
                return Err(Error::Invalid);
            }
        }
    }
    for plane in &planes {
        let p = points.get(&plane.origin).ok_or(Error::Invalid)?;
        if p.coordinates.len() != 3
            || p.coordinates[0].text != "0."
            || p.coordinates[1].text != "0."
        {
            return Err(Error::Unsupported);
        }
        for (id, expected) in [
            (plane.normal, ["0.", "0.", "1."]),
            (plane.reference, ["1.", "0.", "0."]),
        ] {
            direction(&arena, id, expected.map(|v| v == "1."), &mut b)?;
        }
    }
    let mut faces = Vec::new();
    let mut edges = BTreeMap::new();
    let mut vertices = BTreeMap::new();
    let mut uses: BTreeMap<u32, Vec<(usize, bool)>> = BTreeMap::new();
    let mut unique_faces = BTreeSet::new();
    for id in face_ids {
        b.charge(1)?;
        if !unique_faces.insert(id) {
            return Err(Error::Invalid);
        }
        let e = one(&arena, id, "ADVANCED_FACE")?;
        if e.args.len() != 4 {
            return Err(Error::Invalid);
        }
        let surface = reference(&e.args[2])?;
        let kind = *surface_ids.get(&surface).ok_or(Error::Unsupported)?;
        let mut loops = Vec::new();
        let bounds = refs(&e.args[1])?;
        if bounds.is_empty() || bounds.len() > 2 {
            return Err(Error::Unsupported);
        }
        for (ring, bound) in bounds.into_iter().enumerate() {
            let bound = one(
                &arena,
                bound,
                if ring == 0 {
                    "FACE_OUTER_BOUND"
                } else {
                    "FACE_BOUND"
                },
            )?;
            if bound.args.len() != 3 || !boolean(&bound.args[2])? {
                return Err(Error::Unsupported);
            }
            let wire = one(&arena, reference(&bound.args[1])?, "EDGE_LOOP")?;
            let ids = refs(wire.args.get(1).ok_or(Error::Invalid)?)?;
            if ids.len() < 3 {
                return Err(Error::Unsupported);
            }
            let mut coedges = Vec::with_capacity(ids.len());
            for id in ids {
                b.charge(1)?;
                let e = one(&arena, id, "ORIENTED_EDGE")?;
                if e.args.len() != 5
                    || !matches!(e.args[1], Value::Derived)
                    || !matches!(e.args[2], Value::Derived)
                {
                    return Err(Error::Invalid);
                }
                let id = reference(&e.args[3])?;
                let forward = boolean(&e.args[4])?;
                if let std::collections::btree_map::Entry::Vacant(slot) = edges.entry(id) {
                    slot.insert(edge(&arena, id, &mut b)?);
                }
                let edge = edges.get(&id).ok_or(Error::Invalid)?;
                if !edge.pcurves.iter().any(|pc| pc.surface == surface) {
                    return Err(Error::Invalid);
                }
                uses.entry(id).or_default().push((faces.len(), forward));
                coedges.push(Coedge { edge: id, forward });
            }
            for (i, c) in coedges.iter().enumerate() {
                let edge = &edges[&c.edge];
                let next = &coedges[(i + 1) % coedges.len()];
                let next_edge = &edges[&next.edge];
                if edge.vertices[usize::from(c.forward)]
                    != next_edge.vertices[usize::from(!next.forward)]
                {
                    return Err(Error::Invalid);
                }
            }
            loops.push(coedges);
        }
        faces.push(Face {
            entity: id,
            surface,
            kind,
            same_sense: boolean(&e.args[3])?,
            loops,
        });
    }
    for (&id, edge) in &edges {
        b.charge(1)?;
        let owners = &uses[&id];
        if owners.len() != 2 || owners[0].0 == owners[1].0 || owners[0].1 == owners[1].1 {
            return Err(Error::Invalid);
        }
        let mut paired = edge.pcurves.iter().map(|pc| pc.surface).collect::<Vec<_>>();
        paired.sort_unstable();
        let mut actual = owners
            .iter()
            .map(|(f, _)| faces[*f].surface)
            .collect::<Vec<_>>();
        actual.sort_unstable();
        if paired != actual {
            return Err(Error::Invalid);
        }
        for &id in &edge.vertices {
            let e = one(&arena, id, "VERTEX_POINT")?;
            let point = reference(e.args.get(1).ok_or(Error::Invalid)?)?;
            if points.get(&point).is_none_or(|p| p.coordinates.len() != 3) {
                return Err(Error::Invalid);
            }
            vertices.insert(id, point);
        }
        if curves
            .get(&edge.curve)
            .is_none_or(|c| points[&c.points[0]].coordinates.len() != 3)
        {
            return Err(Error::Invalid);
        }
        for pc in &edge.pcurves {
            if curves
                .get(&pc.curve)
                .is_none_or(|c| points[&c.points[0]].coordinates.len() != 2)
            {
                return Err(Error::Invalid);
            }
        }
    }
    super::links::validate(&faces, &edges, vertices.len(), &mut b)?;
    let contribution = faces
        .iter()
        .map(|f| 2 - f.loops.len() as isize)
        .sum::<isize>();
    if vertices.len() as isize - edges.len() as isize + contribution != 2 {
        return Err(Error::Unsupported);
    }
    let mut reached = BTreeSet::from([0_usize]);
    loop {
        let before = reached.len();
        for owners in uses.values() {
            if owners.iter().any(|(f, _)| reached.contains(f)) {
                reached.extend(owners.iter().map(|(f, _)| *f));
            }
        }
        b.charge(uses.len() * 2)?;
        if reached.len() == before {
            break;
        }
    }
    if reached.len() != faces.len() {
        return Err(Error::Invalid);
    }
    b.retain(
        tensors.capacity() * size_of::<Tensor<'_>>()
            + planes.capacity() * size_of::<Plane>()
            + faces.capacity() * size_of::<Face>()
            + vertices.len() * 1024
            + edges.len() * (1024 + size_of::<Edge>())
            + faces
                .iter()
                .map(|f| {
                    f.loops.capacity() * size_of::<Vec<Coedge>>()
                        + f.loops
                            .iter()
                            .map(|l| l.capacity() * size_of::<Coedge>())
                            .sum::<usize>()
                })
                .sum::<usize>(),
    )?;
    drop(arena);
    b.bytes -= arena_bytes;
    Ok(Decoded {
        input,
        points,
        curves,
        tensors,
        planes,
        faces,
        edges,
        vertices,
        solid,
        shell,
        work: b.work,
        bytes: b.bytes,
    })
}
