//! Selected dependency traversal borrows the neutral arena; no AST is copied.
use super::*;
use std::collections::BTreeSet;
const MAX_SELECTED: usize = 327_680;
fn tree_bytes(count: usize) -> usize {
    // Same pinned BTree minimum occupancy bound as CompoundDocument packets.
    (1 + count / 5) * (20 * size_of::<usize>() + 11 * size_of::<u32>())
}
fn insert(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    id: u32,
    selected: &mut BTreeSet<u32>,
    pending: &mut Vec<u32>,
    budget: &mut Budget<'_>,
) -> Result<()> {
    budget.step(96)?;
    if selected.contains(&id) {
        return Ok(());
    }
    // Budget::new seals maximum_work<=16M. EVERY unique dependency, including
    // the seed, passed the96-unit charge above on this same ledger. Thus the
    // count is below floor(16M/96), which is below this storage frontier.
    assert!(
        selected.len() < MAX_SELECTED,
        "continuing selected work bound"
    );
    if !records.contains_key(&id) {
        return Err(Error::Invalid);
    }
    let old = tree_bytes(selected.len());
    let next = tree_bytes(selected.len() + 1);
    if next - old > budget.limit - budget.bytes {
        return Err(Error::Limit);
    }
    budget.bytes += next - old;
    budget.peak = budget.peak.max(budget.bytes);
    selected.insert(id);
    pending.push(id);
    Ok(())
}
fn references(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    values: &[Value<'_>],
    selected: &mut BTreeSet<u32>,
    pending: &mut Vec<u32>,
    budget: &mut Budget<'_>,
) -> Result<()> {
    for value in values {
        budget.step(1)?;
        match value {
            Value::Reference(id) => insert(records, *id, selected, pending, budget)?,
            Value::List(values) => references(records, values, selected, pending, budget)?,
            Value::Typed(entity) => references(records, &entity.args, selected, pending, budget)?,
            _ => {}
        }
    }
    Ok(())
}
pub(super) fn points<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    item: u32,
    budget: &mut Budget<'_>,
) -> Result<(Vec<Point<'input>>, Vec<u32>)> {
    let mut pending = budget.allocate(MAX_SELECTED)?;
    let tree = tree_bytes(0);
    if tree > budget.limit - budget.bytes {
        return Err(Error::Limit);
    }
    budget.bytes += tree;
    budget.peak = budget.peak.max(budget.bytes);
    let mut selected = BTreeSet::new();
    insert(records, item, &mut selected, &mut pending, budget)?;
    while let Some(id) = pending.pop() {
        budget.step(32)?;
        for entity in &records[&id] {
            references(records, &entity.args, &mut selected, &mut pending, budget)?;
        }
    }
    budget.release(&pending);
    drop(pending);
    budget.step(selected.len() * 64)?;
    let count=selected.iter().filter(|id|matches!(records[id].as_slice(),[Entity{name:"CARTESIAN_POINT",args,..}] if matches!(args.as_slice(),[_,Value::List(p)]if p.len()==3))).count();
    // Each counted point is one already charged unique dependency. The same
    // floor(16M/96) bound is strictly below the Cartesian registry cap too.
    assert!(
        count <= 262_144 + 4096,
        "continuing selected point work bound"
    );
    let mut points = budget.allocate(count)?;
    for id in &selected {
        append_point(*id, &records[id], &mut points, budget)?;
    }
    let mut dependencies = budget.allocate(selected.len())?;
    dependencies.extend(selected.iter().copied());
    budget.bytes -= tree_bytes(selected.len());
    drop(selected);
    Ok((points, dependencies))
}
