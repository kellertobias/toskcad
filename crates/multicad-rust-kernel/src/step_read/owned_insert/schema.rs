//! Exact schema metadata and original finite carrier records.
use super::*;
pub(super) fn point<'i>(id: u32, e: &Entity<'i>, w: &mut Ledger<'_>) -> Result<PointLeaf<'i>> {
    let [Value::String(_), coordinates] = e.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let values = list(coordinates)?;
    if !(2..=3).contains(&values.len()) {
        return Err(Error::Unsupported);
    }
    let mut n = Numeric::new(e, w)?;
    let coordinates = values
        .iter()
        .map(|v| n.scalar(v))
        .collect::<Result<Vec<_>>>()?;
    n.finish()?;
    Ok(PointLeaf {
        entity: id,
        coordinates,
    })
}
fn integers(values: &Value<'_>, n: &mut Numeric<'_>) -> Result<Vec<u32>> {
    list(values)?.iter().map(|v| n.integer(v)).collect()
}
fn scalars<'i>(values: &Value<'_>, n: &mut Numeric<'i>) -> Result<Vec<Scalar<'i>>> {
    list(values)?.iter().map(|v| n.scalar(v)).collect()
}
fn knot_domain(degree: u32, poles: usize, mults: &[u32], knots: &[Scalar<'_>]) -> Result<()> {
    if !(2..=1025).contains(&knots.len()) {
        return Err(Error::Limit);
    }
    if knots.len() != mults.len()
        || poles != (knots.len() - 1) * degree as usize + 1
        || !knots[0].is_zero()
        || !knots.last().unwrap().is_one()
    {
        return Err(Error::Invalid);
    }
    for (i, &m) in mults.iter().enumerate() {
        let expected = if i == 0 || i + 1 == mults.len() {
            degree + 1
        } else {
            degree
        };
        if m != expected {
            return Err(Error::Unsupported);
        }
    }
    if knots.windows(2).any(|p| p[0].value() >= p[1].value()) {
        return Err(Error::Invalid);
    }
    Ok(())
}
pub(super) fn curve<'i>(id: u32, e: &Entity<'i>, w: &mut Ledger<'_>) -> Result<SplineCurve<'i>> {
    let [
        Value::String(_),
        degree,
        points,
        Value::Enum("UNSPECIFIED"),
        Value::Enum("F"),
        Value::Enum("F"),
        multiplicities,
        knots,
        Value::Enum("UNSPECIFIED"),
    ] = e.args.as_slice()
    else {
        return Err(Error::Unsupported);
    };
    if !(2..=1025).contains(&list(knots)?.len()) || list(points)?.len() > 11 * 1024 + 1 {
        return Err(Error::Limit);
    }
    let mut n = Numeric::new(e, w)?;
    let degree = n.integer(degree)?;
    if !matches!(degree, 1 | 11) {
        return Err(Error::Unsupported);
    }
    let poles = refs(points)?;
    let multiplicities = integers(multiplicities, &mut n)?;
    let knots = scalars(knots, &mut n)?;
    n.finish()?;
    knot_domain(degree, poles.len(), &multiplicities, &knots)?;
    Ok(SplineCurve {
        entity: id,
        degree,
        poles,
        multiplicities,
        knots,
    })
}
pub(super) fn tensor<'i>(id: u32, e: &Entity<'i>, w: &mut Ledger<'_>) -> Result<TensorSurface<'i>> {
    let [
        Value::String(_),
        du,
        dv,
        net,
        Value::Enum("UNSPECIFIED"),
        Value::Enum("F"),
        Value::Enum("F"),
        Value::Enum("F"),
        mu,
        mv,
        ku,
        kv,
        Value::Enum("UNSPECIFIED"),
    ] = e.args.as_slice()
    else {
        return Err(Error::Unsupported);
    };
    if !(2..=1025).contains(&list(kv)?.len())
        || list(mv)?.len() > 1025
        || list(mu)?.len() > 2
        || list(ku)?.len() > 2
    {
        return Err(Error::Limit);
    }
    let mut n = Numeric::new(e, w)?;
    let degrees = [n.integer(du)?, n.integer(dv)?];
    if degrees != [1, 11] {
        return Err(Error::Unsupported);
    }
    let [a, b] = list(net)? else {
        return Err(Error::Unsupported);
    };
    let rows = [refs(a)?, refs(b)?];
    if rows[0].len() != rows[1].len() {
        return Err(Error::Invalid);
    }
    let mu = integers(mu, &mut n)?;
    let multiplicities_v = integers(mv, &mut n)?;
    let ku = scalars(ku, &mut n)?;
    let knots_v = scalars(kv, &mut n)?;
    n.finish()?;
    knot_domain(1, 2, &mu, &ku)?;
    // knot_domain(1, 2, ...) already proves ku.len()==2 from its exact
    // poles==(knots.len()-1)*degree+1 invariant.
    knot_domain(11, rows[0].len(), &multiplicities_v, &knots_v)?;
    Ok(TensorSurface {
        entity: id,
        rows,
        multiplicities_v,
        knots_v,
    })
}
pub(super) fn plane(
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    id: u32,
    e: &Entity<'_>,
) -> Result<PlaneSurface> {
    let [Value::String(_), axis] = e.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let axis = reference(axis)?;
    let [Value::String(_), origin, normal, reference_axis] =
        kind(entities, axis, "AXIS2_PLACEMENT_3D")?.args.as_slice()
    else {
        return Err(Error::Unsupported);
    };
    Ok(PlaneSurface {
        entity: id,
        axis,
        origin: reference(origin)?,
        normal: reference(normal)?,
        reference: reference(reference_axis)?,
    })
}
pub(super) fn face(
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    id: u32,
    tensors: &[TensorSurface<'_>],
    planes: &[PlaneSurface],
    w: &mut Ledger<'_>,
) -> Result<Face> {
    let [Value::String(_), bounds, support, same_sense] =
        kind(entities, id, "ADVANCED_FACE")?.args.as_slice()
    else {
        return Err(Error::Invalid);
    };
    let surface_entity = reference(support)?;
    let surface = if let Ok(i) = tensors.binary_search_by_key(&surface_entity, |s| s.entity) {
        SurfaceIndex::Tensor(i)
    } else if let Ok(i) = planes.binary_search_by_key(&surface_entity, |s| s.entity) {
        SurfaceIndex::Plane(i)
    } else {
        return Err(Error::Unsupported);
    };
    let mut loops = Vec::new();
    for bound in refs(bounds)? {
        w.step(1)?;
        let e = simple(entities, bound)?;
        let outer = match e.name {
            "FACE_OUTER_BOUND" => true,
            "FACE_BOUND" => false,
            _ => return Err(Error::Unsupported),
        };
        let [Value::String(_), wire, Value::Enum("T")] = e.args.as_slice() else {
            return Err(Error::Unsupported);
        };
        let wire = reference(wire)?;
        let [Value::String(_), uses] = kind(entities, wire, "EDGE_LOOP")?.args.as_slice() else {
            return Err(Error::Invalid);
        };
        let mut coedges = Vec::new();
        for co in refs(uses)? {
            w.step(1)?;
            let [
                Value::String(_),
                Value::Derived,
                Value::Derived,
                edge,
                forward,
            ] = kind(entities, co, "ORIENTED_EDGE")?.args.as_slice()
            else {
                return Err(Error::Unsupported);
            };
            coedges.push(Coedge {
                edge: reference(edge)?,
                forward: sense(forward)?,
            });
        }
        if coedges.len() < 3 {
            return Err(Error::Invalid);
        }
        loops.push(Loop {
            bound,
            entity: wire,
            outer,
            coedges,
        });
    }
    let expected = match surface {
        SurfaceIndex::Tensor(_) => 1,
        SurfaceIndex::Plane(_) => 2,
    };
    if loops.len() != expected || !loops[0].outer || loops.iter().skip(1).any(|l| l.outer) {
        return Err(Error::Unsupported);
    }
    Ok(Face {
        entity: id,
        surface_entity,
        surface,
        same_sense: sense(same_sense)?,
        loops,
    })
}
pub(super) fn edge(
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    id: u32,
    w: &mut Ledger<'_>,
) -> Result<Edge> {
    w.step(1)?;
    let [Value::String(_), a, b, curve, Value::Enum("T")] =
        kind(entities, id, "EDGE_CURVE")?.args.as_slice()
    else {
        return Err(Error::Unsupported);
    };
    let ends = [reference(a)?, reference(b)?];
    if ends[0] == ends[1] {
        return Err(Error::Unsupported);
    }
    let surface_curve = reference(curve)?;
    let [Value::String(_), curve, associated, Value::Enum("CURVE_3D")] =
        kind(entities, surface_curve, "SURFACE_CURVE")?
            .args
            .as_slice()
    else {
        return Err(Error::Unsupported);
    };
    let associated = refs(associated)?;
    let [a, b] = associated.as_slice() else {
        return Err(Error::Unsupported);
    };
    let pcurves = [pcurve(entities, *a)?, pcurve(entities, *b)?];
    Ok(Edge {
        entity: id,
        ends,
        curve_entity: reference(curve)?,
        surface_curve,
        pcurves,
    })
}
fn pcurve(entities: &BTreeMap<u32, Vec<Entity<'_>>>, id: u32) -> Result<Pcurve> {
    let [Value::String(_), support, rep] = kind(entities, id, "PCURVE")?.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let representation = reference(rep)?;
    let [Value::String(_), curves, context] =
        kind(entities, representation, "DEFINITIONAL_REPRESENTATION")?
            .args
            .as_slice()
    else {
        return Err(Error::Invalid);
    };
    let curves = refs(curves)?;
    let [curve] = curves.as_slice() else {
        return Err(Error::Unsupported);
    };
    Ok(Pcurve {
        entity: id,
        support_entity: reference(support)?,
        curve_entity: *curve,
        representation,
        context: reference(context)?,
    })
}
