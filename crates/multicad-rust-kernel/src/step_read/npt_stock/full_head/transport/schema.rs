//! Only explicit represented-member ownership is removed from the local view.
use super::*;
use crate::step_read::npt_stock::records::{kind, list, reference};
pub(super) struct Envelope {
    pub top: u32,
    pub solid: u32,
    pub axis: u32,
    pub context: u32,
    pub item: u32,
    pub operator: u32,
    pub removable: [u32; 16],
    pub count: usize,
}
impl Envelope {
    fn allow(&mut self, id: u32) -> Result<()> {
        if self.removable[..self.count].contains(&id) {
            return Ok(());
        }
        let slot = self
            .removable
            .get_mut(self.count)
            .ok_or(Error::Unsupported)?;
        *slot = id;
        self.count += 1;
        Ok(())
    }
    fn axis(
        &mut self,
        records: &BTreeMap<u32, Vec<Entity<'_>>>,
        id: u32,
        b: &mut Budget<'_>,
    ) -> Result<()> {
        let e = kind(records, id, "AXIS2_PLACEMENT_3D", b)?;
        let [Value::String(_), p, z, x] = e.args.as_slice() else {
            return Err(Error::Unsupported);
        };
        self.allow(id)?;
        for (value, name) in [(p, "CARTESIAN_POINT"), (z, "DIRECTION"), (x, "DIRECTION")] {
            let id = reference(value)?;
            leaf(records, id, name, b)?;
            self.allow(id)?;
        }
        Ok(())
    }
}
fn leaf(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    id: u32,
    name: &str,
    b: &mut Budget<'_>,
) -> Result<()> {
    let e = kind(records, id, name, b)?;
    let p = crate::step_read::npt_stock::records::point(e, 3, b)?;
    b.release(&p);
    Ok(())
}
pub(super) fn read(doc: &CompoundDocument<'_>, b: &mut Budget<'_>) -> Result<Envelope> {
    b.step(128)?;
    let [member] = doc.represented_members() else {
        return Err(Error::Unsupported);
    };
    let link = member.mapping_link().ok_or(Error::Unsupported)?;
    let records = doc.records();
    let mut top = None;
    for (&id, parts) in records {
        b.step(1)?;
        if matches!(parts.as_slice(),[e] if e.name=="ADVANCED_BREP_SHAPE_REPRESENTATION")
            && top.replace(id).is_some()
        {
            return Err(Error::Unsupported);
        }
    }
    let top = top.ok_or(Error::Unsupported)?;
    let e = kind(records, top, "ADVANCED_BREP_SHAPE_REPRESENTATION", b)?;
    let [Value::String(_), items, context] = e.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    let items = list(items)?;
    if items.len() != 2 {
        return Err(Error::Unsupported);
    }
    let ids = [reference(&items[0])?, reference(&items[1])?];
    if ids[0] == ids[1] || !ids.contains(&member.item_entity()) {
        return Err(Error::Unsupported);
    }
    let axis = ids
        .into_iter()
        .find(|id| *id != member.item_entity())
        .ok_or(Error::Unsupported)?;
    kind(records, axis, "AXIS2_PLACEMENT_3D", b)?;
    let context = reference(context)?;
    let mut out = Envelope {
        top,
        solid: member.local_solid_entity(),
        axis,
        context,
        item: member.item_entity(),
        operator: link.operator_entity(),
        removable: [0; 16],
        count: 0,
    };
    let item = kind(records, out.item, "MAPPED_ITEM", b)?;
    if !matches!(item.args.as_slice(),[Value::String(_),Value::Reference(m),Value::Reference(o)] if *m==link.map_entity()&&*o==out.operator)
    {
        return Err(Error::Unsupported);
    }
    let map = kind(records, link.map_entity(), "REPRESENTATION_MAP", b)?;
    if !matches!(map.args.as_slice(),[Value::Reference(a),Value::Reference(r)] if *a==link.source_origin_entity()&&*r==link.source_representation_entity())
    {
        return Err(Error::Unsupported);
    }
    let rep = kind(
        records,
        link.source_representation_entity(),
        "SHAPE_REPRESENTATION",
        b,
    )?;
    let [Value::String(_), items, ctx] = rep.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    let items = list(items)?;
    if reference(ctx)? != context || items.len() != 2 {
        return Err(Error::Unsupported);
    }
    let ids = [reference(&items[0])?, reference(&items[1])?];
    if ids[0] == ids[1] || !ids.contains(&out.solid) || !ids.contains(&link.source_origin_entity())
    {
        return Err(Error::Unsupported);
    }
    // The neutral Member parser already proves this exact source datum is
    // identity. The unchanged top-axis is independently checked in local view.
    out.axis(records, link.source_origin_entity(), b)?;
    for id in [
        out.item,
        link.map_entity(),
        link.source_representation_entity(),
        out.operator,
    ] {
        out.allow(id)?;
    }
    let operator = records.get(&out.operator).ok_or(Error::Invalid)?;
    // The actual neutral operator evaluator validates complete component
    // names/arity, all scalar leaves, EXPRESS projection and handedness.
    // Here only its explicitly typed coordinate dependencies may disappear.
    let base = operator
        .iter()
        .find(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
        .ok_or(Error::Unsupported)?;
    let dim = operator
        .iter()
        .find(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR_3D")
        .ok_or(Error::Unsupported)?;
    let [x, y, p, _] = base.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    let [z] = dim.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    for (value, name) in [
        (p, "CARTESIAN_POINT"),
        (x, "DIRECTION"),
        (y, "DIRECTION"),
        (z, "DIRECTION"),
    ] {
        if matches!(value, Value::Omitted) && name == "DIRECTION" {
            continue;
        }
        let id = reference(value)?;
        leaf(records, id, name, b)?;
        out.allow(id)?;
    }
    Ok(out)
}
