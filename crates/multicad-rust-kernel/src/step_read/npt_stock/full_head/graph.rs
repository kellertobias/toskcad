use super::super::records::{boolean, kind, list, point, reference};
use super::*;
use crate::step::pipe_hardware::numbers as exact;
use crate::tapered_thread::interval::I;
use std::collections::BTreeSet;
pub(super) struct Plan<'i> {
    pub head: HeadInput,
    pub shell: u32,
    pub thread_faces: Vec<u32>,
    pub inner_wire: u32,
    pub physical_plane: u32,
    pub cap_height: Scalar<'i>,
    pub rehomes: Vec<(u32, u32)>,
}
struct Face {
    id: u32,
    support: u32,
    outer: Vec<(u32, bool)>,
    inner: Vec<(u32, bool)>,
    inner_wire: u32,
    forward: bool,
    tensor: bool,
}
fn wire(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    id: u32,
    outer: bool,
    b: &mut Budget<'_>,
) -> Result<(u32, Vec<(u32, bool)>)> {
    let bound = kind(
        records,
        id,
        if outer {
            "FACE_OUTER_BOUND"
        } else {
            "FACE_BOUND"
        },
        b,
    )?;
    let [Value::String(_), wire, Value::Enum("T")] = bound.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    let id = reference(wire)?;
    let wire = kind(records, id, "EDGE_LOOP", b)?;
    let [Value::String(_), uses] = wire.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let uses = list(uses)?;
    if !matches!(uses.len(), 4 | 6) {
        return Err(Error::Unsupported);
    }
    let mut out = b.allocate(uses.len())?;
    for value in uses {
        let e = kind(records, reference(value)?, "ORIENTED_EDGE", b)?;
        let [
            Value::String(_),
            Value::Derived,
            Value::Derived,
            edge,
            forward,
        ] = e.args.as_slice()
        else {
            return Err(Error::Invalid);
        };
        out.push((reference(edge)?, boolean(forward)?));
    }
    Ok((id, out))
}
fn xyz<'i>(
    records: &BTreeMap<u32, Vec<Entity<'i>>>,
    id: u32,
    name: &str,
    require_exact: bool,
    b: &mut Budget<'_>,
) -> Result<[Scalar<'i>; 3]> {
    let e = kind(records, id, name, b)?;
    let p = point(e, 3, b)?;
    let out: [Scalar<'i>; 3] = p.as_slice().try_into().map_err(|_| Error::Invalid)?;
    if require_exact {
        for s in out {
            let (cost, bytes) = exact::cost(s.value()).map_err(geometry_error)?;
            b.step(cost.checked_mul(64).ok_or(Error::Limit)?)?;
            if bytes > 4096
                || !exact::exact_leaf(s.spelling(), s.value()).map_err(geometry_error)?
            {
                return Err(Error::Unsupported);
            }
        }
    }
    b.release(&p);
    Ok(out)
}
fn vertex<'i>(
    records: &BTreeMap<u32, Vec<Entity<'i>>>,
    id: u32,
    exact: bool,
    b: &mut Budget<'_>,
) -> Result<[Scalar<'i>; 3]> {
    let e = kind(records, id, "VERTEX_POINT", b)?;
    let [Value::String(_), point] = e.args.as_slice() else {
        return Err(Error::Invalid);
    };
    xyz(records, reference(point)?, "CARTESIAN_POINT", exact, b)
}
fn plane(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    id: u32,
    b: &mut Budget<'_>,
) -> Result<head::Plane> {
    let e = kind(records, id, "PLANE", b)?;
    let [Value::String(_), axis] = e.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let axis = kind(records, reference(axis)?, "AXIS2_PLACEMENT_3D", b)?;
    let [Value::String(_), p, n, r] = axis.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    Ok(head::Plane {
        origin: xyz(records, reference(p)?, "CARTESIAN_POINT", true, b)?.map(|s| s.value()),
        normal: xyz(records, reference(n)?, "DIRECTION", true, b)?.map(|s| s.value()),
        reference: xyz(records, reference(r)?, "DIRECTION", true, b)?.map(|s| s.value()),
    })
}
fn edge<'a, 'i>(
    records: &'a BTreeMap<u32, Vec<Entity<'i>>>,
    id: u32,
    b: &mut Budget<'_>,
) -> Result<(u32, u32, u32)> {
    let e = kind(records, id, "EDGE_CURVE", b)?;
    let [Value::String(_), a, z, c, Value::Enum("T")] = e.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    Ok((reference(a)?, reference(z)?, reference(c)?))
}
fn references(
    values: &[Value<'_>],
    out: &mut Vec<u32>,
    seen: &mut BTreeSet<u32>,
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    b: &mut Budget<'_>,
) -> Result<()> {
    for v in values {
        b.step(1)?;
        match v {
            Value::Reference(id) => {
                b.step(128)?;
                if !records.contains_key(id) {
                    return Err(Error::Invalid);
                }
                if seen.insert(*id) {
                    out.push(*id);
                }
            }
            Value::List(v) => references(v, out, seen, records, b)?,
            Value::Typed(e) => references(&e.args, out, seen, records, b)?,
            _ => {}
        }
    }
    Ok(())
}
fn consume(records: &BTreeMap<u32, Vec<Entity<'_>>>, b: &mut Budget<'_>) -> Result<()> {
    let n = records.len();
    let tree = (1 + n / 5)
        .checked_mul(20 * size_of::<usize>() + 11 * size_of::<u32>())
        .ok_or(Error::Limit)?;
    b.step(n.checked_mul(32).ok_or(Error::Limit)?)?;
    if tree > b.limit - b.bytes {
        return Err(Error::Limit);
    }
    b.bytes += tree;
    b.peak = b.peak.max(b.bytes);
    let mut pending = b.allocate(n)?;
    let mut seen = BTreeSet::new();
    let mut roots = 0;
    for (&id, parts) in records {
        for e in parts {
            b.step(512 + 4 * e.name.len())?;
            if e.name.contains("TRANSFORM") || e.name.contains("MAPPED") {
                return Err(Error::Unsupported);
            }
            if matches!(
                e.name,
                "MAPPED_ITEM"
                    | "REPRESENTATION_MAP"
                    | "ITEM_DEFINED_TRANSFORMATION"
                    | "CARTESIAN_TRANSFORMATION_OPERATOR_3D"
                    | "CARTESIAN_TRANSFORMATION_OPERATOR_3DNON_UNIFORM"
                    | "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"
            ) {
                return Err(Error::Unsupported);
            }
            if e.name == "SHAPE_DEFINITION_REPRESENTATION" {
                roots += 1;
                seen.insert(id);
                pending.push(id);
            }
        }
    }
    if roots != 1 {
        return Err(Error::Unsupported);
    }
    while let Some(id) = pending.pop() {
        b.step(128)?;
        for e in &records[&id] {
            references(&e.args, &mut pending, &mut seen, records, b)?;
        }
    }
    // These standard metadata records point INTO the product graph; they
    // cannot be reached by following its forward references. They may only
    // name the already consumed context/product, never grant new geometry.
    let mut protocols = 0;
    let mut categories = 0;
    for (&id, parts) in records {
        b.step(128)?;
        if seen.contains(&id) {
            continue;
        }
        let [e] = parts.as_slice() else {
            return Err(Error::Unsupported);
        };
        match e.name {
            "APPLICATION_PROTOCOL_DEFINITION" => {
                protocols += 1;
                let [
                    Value::String("international standard"),
                    Value::String("automotive_design"),
                    year,
                    context,
                ] = e.args.as_slice()
                else {
                    return Err(Error::Unsupported);
                };
                let mut numbers = numeric::Numbers::new(e, b)?;
                if numbers.integer(year, b)? != 2000 {
                    return Err(Error::Unsupported);
                }
                numbers.finish()?;
                let context = reference(context)?;
                kind(records, context, "APPLICATION_CONTEXT", b)?;
                if !seen.contains(&context) {
                    return Err(Error::Unsupported);
                }
            }
            "PRODUCT_RELATED_PRODUCT_CATEGORY" => {
                categories += 1;
                let [Value::String("part"), Value::Omitted, products] = e.args.as_slice() else {
                    return Err(Error::Unsupported);
                };
                let [product] = list(products)? else {
                    return Err(Error::Unsupported);
                };
                let product = reference(product)?;
                kind(records, product, "PRODUCT", b)?;
                if !seen.contains(&product) {
                    return Err(Error::Unsupported);
                }
            }
            _ => return Err(Error::Unsupported),
        }
        seen.insert(id);
    }
    if protocols != 1 || categories != 1 || seen.len() != n {
        return Err(Error::Unsupported);
    }
    b.release(&pending);
    drop(pending);
    drop(seen);
    b.bytes -= tree;
    Ok(())
}
pub(super) fn decode<'i>(doc: &CompoundDocument<'i>, b: &mut Budget<'_>) -> Result<Plan<'i>> {
    let records = doc.records();
    let [member] = doc.represented_members() else {
        return Err(Error::Unsupported);
    };
    if member.mapping_link().is_some() || member.item_entity() != member.local_solid_entity() {
        return Err(Error::Unsupported);
    }
    consume(records, b)?;
    let solid = member.item_entity();
    let e = kind(records, solid, "MANIFOLD_SOLID_BREP", b)?;
    let [Value::String(_), shell] = e.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let shell = reference(shell)?;
    let e = kind(records, shell, "CLOSED_SHELL", b)?;
    let [Value::String(_), faces] = e.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let faces = list(faces)?;
    if faces.len() != 17 {
        return Err(Error::Unsupported);
    }
    let mut decoded = b.allocate(17)?;
    for v in faces {
        let id = reference(v)?;
        if decoded.iter().any(|f: &Face| f.id == id) {
            return Err(Error::Invalid);
        }
        let e = kind(records, id, "ADVANCED_FACE", b)?;
        let [Value::String(_), bounds, support, forward] = e.args.as_slice() else {
            return Err(Error::Invalid);
        };
        let bounds = list(bounds)?;
        if bounds.is_empty() || bounds.len() > 2 {
            return Err(Error::Unsupported);
        }
        let (_, outer) = wire(records, reference(&bounds[0])?, true, b)?;
        let (inner_wire, inner) = if bounds.len() == 2 {
            wire(records, reference(&bounds[1])?, false, b)?
        } else {
            (0, b.allocate(0)?)
        };
        let support = reference(support)?;
        b.step(128)?;
        let tensor = match records.get(&support).map(Vec::as_slice) {
            Some([e]) if e.name == "PLANE" => false,
            Some([e]) if e.name == "B_SPLINE_SURFACE_WITH_KNOTS" => true,
            _ => return Err(Error::Unsupported),
        };
        decoded.push(Face {
            id,
            support,
            outer,
            inner,
            inner_wire,
            forward: boolean(forward)?,
            tensor,
        });
    }
    let mut annulus = None;
    let mut bottom = None;
    let mut thread_faces = b.allocate(9)?;
    for (i, f) in decoded.iter().enumerate() {
        if !f.inner.is_empty() {
            if annulus.replace(i).is_some() || f.tensor || f.inner.len() != 4 || f.outer.len() != 6
            {
                return Err(Error::Unsupported);
            }
        } else if f.tensor {
            thread_faces.push(f.id);
        } else {
            let mut all_spline = true;
            for &(id, _) in &f.outer {
                let (_, _, c) = edge(records, id, b)?;
                b.step(128)?;
                all_spline &= matches!(records.get(&c).map(Vec::as_slice),Some([e]) if e.name=="SURFACE_CURVE");
            }
            if all_spline && bottom.replace(i).is_some() {
                return Err(Error::Unsupported);
            }
        }
    }
    if thread_faces.len() != 8 {
        return Err(Error::Unsupported);
    }
    let annulus = annulus.ok_or(Error::Unsupported)?;
    let bottom = bottom.ok_or(Error::Unsupported)?;
    if decoded[bottom].outer.len() != 4 {
        return Err(Error::Unsupported);
    }
    thread_faces.push(decoded[bottom].id);
    let mut vertex_ids = b.allocate(12)?;
    let mut vertices = b.allocate(12)?;
    let mut edge_ids = b.allocate(18)?;
    let mut edges = b.allocate(18)?;
    let mut head_faces = b.allocate(8)?;
    for (i, f) in decoded.iter().enumerate() {
        if f.tensor || i == bottom {
            continue;
        }
        let p = plane(records, f.support, b)?;
        let mut uses = b.allocate(f.outer.len())?;
        for &(id, forward) in &f.outer {
            let index = if let Some(j) = edge_ids.iter().position(|v| *v == id) {
                j
            } else {
                if edges.len() == 18 {
                    return Err(Error::Unsupported);
                }
                let (a, z, c) = edge(records, id, b)?;
                let mut endpoint = [0usize; 2];
                for (out, id) in endpoint.iter_mut().zip([a, z]) {
                    *out = if let Some(j) = vertex_ids.iter().position(|v| *v == id) {
                        j
                    } else {
                        if vertices.len() == 12 {
                            return Err(Error::Unsupported);
                        }
                        let v = vertex(records, id, true, b)?.map(|s| s.value());
                        let j = vertices.len();
                        vertex_ids.push(id);
                        vertices.push(v);
                        j
                    };
                }
                let e = kind(records, c, "LINE", b)?;
                let [Value::String(_), origin, vector] = e.args.as_slice() else {
                    return Err(Error::Invalid);
                };
                let origin = xyz(records, reference(origin)?, "CARTESIAN_POINT", true, b)?
                    .map(|s| s.value());
                let vector = kind(records, reference(vector)?, "VECTOR", b)?;
                let [Value::String(_), direction, magnitude] = vector.args.as_slice() else {
                    return Err(Error::Invalid);
                };
                let mut n = numeric::Numbers::new(vector, b)?;
                let magnitude = n.scalar(magnitude, b)?;
                n.finish()?;
                b.step(8192)?;
                if !exact::exact_leaf(magnitude.spelling(), 1.).map_err(geometry_error)? {
                    return Err(Error::Unsupported);
                }
                let delta =
                    xyz(records, reference(direction)?, "DIRECTION", true, b)?.map(|s| s.value());
                if origin != vertices[endpoint[0]] {
                    return Err(Error::Unsupported);
                }
                let mut norm = I::point(0.);
                for axis in 0..3 {
                    b.step(8192)?;
                    if !exact::exact_difference(
                        vertices[endpoint[1]][axis],
                        vertices[endpoint[0]][axis],
                        delta[axis],
                    )
                    .map_err(geometry_error)?
                    {
                        return Err(Error::Unsupported);
                    }
                    norm = norm + I::point(delta[axis]) * I::point(delta[axis]);
                }
                let length = norm.sqrt();
                if length.lo <= 0. || !length.hi.is_finite() {
                    return Err(Error::Unsupported);
                }
                let j = edges.len();
                edge_ids.push(id);
                edges.push(endpoint);
                j
            };
            uses.push((index, forward));
        }
        let mut inner = b.allocate(f.inner.len())?;
        inner.extend_from_slice(&f.inner);
        head_faces.push(head::HeadFace {
            plane: p,
            outer: uses,
            inner,
            forward: f.forward,
        });
    }
    if vertices.len() != 12 || edges.len() != 18 || head_faces.len() != 8 {
        return Err(Error::Unsupported);
    }
    let a = &decoded[annulus];
    let mut rehomes = b.allocate(4)?;
    let mut height = None;
    for &(id, _) in &a.inner {
        let (vertex_id, _, c) = edge(records, id, b)?;
        let p = vertex(records, vertex_id, false, b)?;
        if height.is_none() {
            height = Some(p[2]);
        }
        let curve = kind(records, c, "SURFACE_CURVE", b)?;
        let [Value::String(_), _, pcurves, Value::Enum("PCURVE_S1")] = curve.args.as_slice() else {
            return Err(Error::Unsupported);
        };
        let pcurves = list(pcurves)?;
        if pcurves.len() != 2 {
            return Err(Error::Unsupported);
        }
        let mut found = 0;
        for value in pcurves {
            let id = reference(value)?;
            let p = kind(records, id, "PCURVE", b)?;
            let [Value::String(_), support, representation] = p.args.as_slice() else {
                return Err(Error::Invalid);
            };
            if reference(support)? == a.support {
                found += 1;
                rehomes.push((id, reference(representation)?));
            }
        }
        if found != 1 {
            return Err(Error::Unsupported);
        }
    }
    let head = HeadInput {
        vertices,
        edges,
        faces: head_faces,
        auxiliary_cap_entity: 0,
    };
    // Fixed scratch owners no longer survive the decoded packet.
    b.release(&vertex_ids);
    b.release(&edge_ids);
    for f in &decoded {
        b.release(&f.outer);
        b.release(&f.inner);
    }
    b.release(&decoded);
    Ok(Plan {
        head,
        shell,
        thread_faces,
        inner_wire: a.inner_wire,
        physical_plane: a.support,
        cap_height: height.ok_or(Error::Invalid)?,
        rehomes,
    })
}
