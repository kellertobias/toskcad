//! Bounded non-transform representation ownership shared by file readers.
use super::*;
// A plain relationship has no placement operator. We retain ownership only;
// datum-only wrappers cannot silently contribute omitted geometry. The EXPRESS
// graph treats representation relationships as undirected links (ISO10303-41).
const ROOT_LINK_NODES: usize = 64;
const ROOT_DATUM_ITEMS: usize = 64;
pub(super) fn product_shape<E: From<ReadError>>(
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    top: u32,
    charge: &mut dyn FnMut(usize) -> std::result::Result<(), E>,
) -> std::result::Result<u32, E> {
    fn owner<E: From<ReadError>>(
        entities: &BTreeMap<u32, Vec<Entity<'_>>>,
        nodes: &[(u32, Option<u32>)],
        charge: &mut dyn FnMut(usize) -> std::result::Result<(), E>,
    ) -> std::result::Result<Option<u32>, E> {
        let mut owner = None;
        for records in entities.values() {
            charge(1)?;
            let [e] = records.as_slice() else { continue };
            if e.name != "SHAPE_DEFINITION_REPRESENTATION" {
                continue;
            }
            let [Value::Reference(shape), Value::Reference(rep)] = e.args.as_slice() else {
                return Err(ReadError::Invalid.into());
            };
            charge(nodes.len())?;
            if nodes.iter().any(|(id, _)| id == rep)
                && Shell::new(entities, 1e-7, false)
                    .args(*shape, "PRODUCT_DEFINITION_SHAPE")
                    .is_ok()
            {
                if owner.is_some_and(|old| old != *shape) {
                    return Err(ReadError::Unsupported.into());
                }
                owner = Some(*shape);
            }
        }
        Ok(owner)
    }
    let mut nodes = [(0, None); ROOT_LINK_NODES];
    nodes[0] = (top, None);
    // Preserve ordinary direct ownership without expanding unrelated metadata.
    if let Some(shape) = owner(entities, &nodes[..1], charge)? {
        return Ok(shape);
    }
    let mut count = 1;
    let mut next = 0;
    let shell = Shell::new(entities, 1e-7, false);
    while next < count {
        let (current, incoming) = nodes[next];
        for (&id, parts) in entities {
            charge(1)?;
            let [e] = parts.as_slice() else {
                charge(parts.len().saturating_mul(2))?;
                if parts
                    .iter()
                    .any(|e| e.name == "SHAPE_REPRESENTATION_RELATIONSHIP")
                    && parts.iter().any(|e| {
                        matches!(e.name, "REPRESENTATION_RELATIONSHIP" | "SHAPE_REPRESENTATION_RELATIONSHIP")
                            && matches!(e.args.as_slice(), [_, _, Value::Reference(a), Value::Reference(b)]
                                if *a == current || *b == current)
                    })
                {
                    return Err(ReadError::Unsupported.into());
                }
                continue;
            };
            if e.name != "SHAPE_REPRESENTATION_RELATIONSHIP" {
                continue;
            }
            let [
                Value::String(_),
                description,
                Value::Reference(a),
                Value::Reference(b),
            ] = e.args.as_slice()
            else {
                return Err(ReadError::Invalid.into());
            };
            if !matches!(description, Value::String(_) | Value::Omitted) {
                return Err(ReadError::Invalid.into());
            }
            if incoming == Some(id) || (*a != current && *b != current) {
                continue;
            }
            let other = if *a == current { *b } else { *a };
            // Cyclic/parallel relationship paths have no unique ownership
            // route in this bounded view. They are unsupported, not repaired.
            charge(count)?;
            if nodes[..count].iter().any(|(node, _)| *node == other) {
                return Err(ReadError::Unsupported.into());
            }
            let [Value::String(_), items, context] = shell.args(other, "SHAPE_REPRESENTATION")?
            else {
                return Err(ReadError::Invalid.into());
            };
            shell.check_units(reference(context)?)?;
            // Borrow the parser arena: no temporary Vec or tree allocation is
            // hidden outside the caller's retained-arena byte accounting.
            let datums = list(items)?;
            if datums.len() > ROOT_DATUM_ITEMS {
                return Err(ReadError::Limit.into());
            }
            for (i, value) in datums.iter().enumerate() {
                charge(i + 1)?;
                let datum = reference(value)?;
                if datums[..i]
                    .iter()
                    .any(|old| matches!(old, Value::Reference(id) if *id == datum))
                {
                    return Err(ReadError::Invalid.into());
                }
                shell.args(datum, "AXIS2_PLACEMENT_3D")?;
                shell.placement(datum)?;
            }
            if count == ROOT_LINK_NODES {
                return Err(ReadError::Limit.into());
            }
            nodes[count] = (other, Some(id));
            count += 1;
        }
        next += 1;
    }
    owner(entities, &nodes[..count], charge)?.ok_or_else(|| ReadError::Invalid.into())
}
