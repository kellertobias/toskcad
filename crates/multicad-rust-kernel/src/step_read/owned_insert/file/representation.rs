//! Explicit represented-instance ownership. No operator/body geometry admission.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MappingLink {
    map: u32,
    representation: u32,
    origin: u32,
    operator: u32,
}
impl MappingLink {
    pub fn map_entity(self) -> u32 {
        self.map
    }
    pub fn source_representation_entity(self) -> u32 {
        self.representation
    }
    pub fn source_origin_entity(self) -> u32 {
        self.origin
    }
    pub fn operator_entity(self) -> u32 {
        self.operator
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Member {
    item: u32,
    solid: u32,
    mapping: Option<MappingLink>,
}
impl Member {
    pub fn item_entity(self) -> u32 {
        self.item
    }
    pub fn local_solid_entity(self) -> u32 {
        self.solid
    }
    pub fn mapping_link(self) -> Option<MappingLink> {
        self.mapping
    }
}
pub(super) fn read(
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    w: &mut Ledger<'_>,
) -> Result<Vec<Member>> {
    let mut reps = entities.iter().filter_map(|(&id, r)| {
        matches!(r.as_slice(),[e]if e.name=="ADVANCED_BREP_SHAPE_REPRESENTATION").then_some(id)
    });
    let rep = reps.next().ok_or(Error::Unsupported)?;
    if reps.next().is_some() {
        return Err(Error::Unsupported);
    }
    let [Value::String(_), items, context] =
        kind(entities, rep, "ADVANCED_BREP_SHAPE_REPRESENTATION")?
            .args
            .as_slice()
    else {
        return Err(Error::Invalid);
    };
    let linked=entities.values().any(|r|matches!(r.as_slice(),[Entity{name:"SHAPE_DEFINITION_REPRESENTATION",args,..}]if matches!(args.as_slice(),[Value::Reference(shape),Value::Reference(target)]if *target==rep&&kind(entities,*shape,"PRODUCT_DEFINITION_SHAPE").is_ok())));
    if !linked {
        return Err(Error::Invalid);
    }
    let shell = Shell::new(entities, 1e-7, false);
    shell
        .check_units(reference(context)?)
        .map_err(Error::from)?;
    let mut seen = BTreeSet::new();
    let mut members = Vec::new();
    for id in refs(items)? {
        w.step(1)?;
        if !seen.insert(id) {
            return Err(Error::Invalid);
        }
        match simple(entities, id)?.name {
            "AXIS2_PLACEMENT_3D" => {
                shell.placement(id).map_err(Error::from)?;
            }
            "MANIFOLD_SOLID_BREP" | "BREP_WITH_VOIDS" => members.push(Member {
                item: id,
                solid: id,
                mapping: None,
            }),
            "MAPPED_ITEM" => members.push(mapped(entities, id, w)?),
            _ => return Err(Error::Unsupported),
        }
    }
    if members.is_empty() {
        return Err(Error::Unsupported);
    }
    if members.len() > crate::step::MAX_COMPOUND_SOLIDS {
        return Err(Error::Limit);
    }
    Ok(members)
}
fn mapped(
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    item: u32,
    w: &mut Ledger<'_>,
) -> Result<Member> {
    let [Value::String(_), source, target] = kind(entities, item, "MAPPED_ITEM")?.args.as_slice()
    else {
        return Err(Error::Invalid);
    };
    let map = reference(source)?;
    let operator = reference(target)?;
    let [origin, rep] = kind(entities, map, "REPRESENTATION_MAP")?.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let origin = reference(origin)?;
    let rep = reference(rep)?;
    let [Value::String(_), items, context] =
        kind(entities, rep, "SHAPE_REPRESENTATION")?.args.as_slice()
    else {
        return Err(Error::Invalid);
    };
    Shell::new(entities, 1e-7, false)
        .check_units(reference(context)?)
        .map_err(Error::from)?;
    let ids = refs(items)?;
    if ids.len() != 2 || ids[0] == ids[1] || !ids.contains(&origin) {
        return Err(Error::Unsupported);
    }
    let solid = ids
        .into_iter()
        .find(|id| *id != origin)
        .expect("two distinct items include origin");
    kind(entities, solid, "MANIFOLD_SOLID_BREP")?;
    // Initial supported map source is exactly the declared identity datum.
    // Never fit a source placement or borrow file tolerance as datum uncertainty.
    identity(entities, origin, w)?;
    let target = entities.get(&operator).ok_or(Error::Invalid)?;
    if !target
        .iter()
        .any(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR_3D")
    {
        return Err(Error::Unsupported);
    }
    Ok(Member {
        item,
        solid,
        mapping: Some(MappingLink {
            map,
            representation: rep,
            origin,
            operator,
        }),
    })
}
fn identity(entities: &BTreeMap<u32, Vec<Entity<'_>>>, id: u32, w: &mut Ledger<'_>) -> Result<()> {
    let [Value::String(_), point, normal, reference_axis] =
        kind(entities, id, "AXIS2_PLACEMENT_3D")?.args.as_slice()
    else {
        return Err(Error::Invalid);
    };
    let point = reference(point)?;
    let p = schema::point(point, kind(entities, point, "CARTESIAN_POINT")?, w)?;
    if p.coordinates().len() != 3 || p.coordinates().iter().any(|v| !v.is_zero()) {
        return Err(Error::Unsupported);
    }
    for (value, axis) in [(normal, 2), (reference_axis, 0)] {
        match value {
            Value::Omitted => {}
            Value::Reference(id) => {
                let e = kind(entities, *id, "DIRECTION")?;
                let [Value::String(_), values] = e.args.as_slice() else {
                    return Err(Error::Invalid);
                };
                let values = list(values)?;
                if values.len() != 3 {
                    return Err(Error::Unsupported);
                }
                let mut n = Numeric::new(e, w)?;
                for (i, value) in values.iter().enumerate() {
                    let v = n.scalar(value)?;
                    if (i == axis && v.value() <= 0.) || (i != axis && !v.is_zero()) {
                        return Err(Error::Unsupported);
                    }
                }
                n.finish()?;
            }
            _ => return Err(Error::Unsupported),
        }
    }
    Ok(())
}
