//! One original represented file arena shared by explicit member queries.
use super::*;
use std::collections::BTreeSet;
mod representation;
mod view;
pub use representation::{MappingLink, Member};
pub use view::MemberDocument;
pub struct File<'input> {
    pub(super) input: &'input str,
    pub(super) entities: Arc<BTreeMap<u32, Vec<Entity<'input>>>>,
    solids: Vec<u32>,
    members: Vec<Member>,
    work: usize,
    bytes: usize,
}
fn valid(options: Options) -> Result<()> {
    if options.maximum_work == 0
        || options.maximum_work > 16_000_000
        || options.retained_budget_bytes == 0
        || options.retained_budget_bytes > 64 * 1024 * 1024
    {
        return Err(Error::Invalid);
    }
    Ok(())
}
impl<'input> File<'input> {
    /// Retain the actual millimetre represented local solids once. This is not
    /// geometric admission: each requested member earns its own finite packet.
    pub fn parse(
        input: &'input str,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        valid(options)?;
        let mut w = Ledger {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        w.step(input.len())?;
        let entities = part21::parse(input).map_err(Error::from)?;
        let mut bytes = budget::arena(input, &entities, &mut w)?;
        // The Arc allocation owns the BTreeMap header and two reference counts.
        bytes = bytes.saturating_add(
            size_of::<Self>()
                + size_of::<BTreeMap<u32, Vec<Entity<'_>>>>()
                + 2 * size_of::<usize>(),
        );
        if bytes > options.retained_budget_bytes {
            return Err(Error::Limit);
        }
        for records in entities.values() {
            w.step(1)?;
            for e in records {
                if e.name == "GEOMETRIC_REPRESENTATION_CONTEXT" {
                    let [dimension] = e.args.as_slice() else {
                        return Err(Error::Invalid);
                    };
                    let mut n = Numeric::new(e, &mut w)?;
                    if !matches!(n.integer(dimension)?, 2 | 3) {
                        return Err(Error::Unsupported);
                    }
                    n.finish()?;
                }
            }
        }
        let members = representation::read(&entities, &mut w)?;
        let solids = members
            .iter()
            .filter(|m| m.mapping_link().is_none())
            .map(|m| m.local_solid_entity())
            .collect::<Vec<_>>();
        w.step(solids.len())?;
        bytes = bytes
            .saturating_add(solids.capacity() * size_of::<u32>())
            .saturating_add(members.capacity() * size_of::<Member>());
        if bytes > options.retained_budget_bytes {
            return Err(Error::Limit);
        }
        w.step(1)?;
        Ok(Self {
            input,
            entities: Arc::new(entities),
            solids,
            members,
            work: w.used,
            bytes,
        })
    }
    pub fn represented_members(&self) -> &[Member] {
        &self.members
    }
    /// Unmapped top-level solid items only; mapped instances are explicit members.
    pub fn represented_solids(&self) -> &[u32] {
        &self.solids
    }
    pub fn literal_input(&self) -> &'input str {
        self.input
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
    pub fn shares_arena(&self, member: &Parsed<'input>) -> bool {
        Arc::ptr_eq(&self.entities, &member.entities)
    }
    /// Follow only this represented member's actual finite face reference graph.
    /// Unrelated member carriers never select its planes or redefine its domain.
    /// The full immutable original arena remains retained; query scratch and its
    /// retained source charge are reserved once, without reparsing or AST clones.
    pub fn member(
        &self,
        solid: u32,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Parsed<'input>> {
        valid(options)?;
        let mut w = Ledger {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        w.step(self.work)?;
        w.step(self.members.len())?;
        let member = *self
            .members
            .iter()
            .find(|m| m.item_entity() == solid)
            .ok_or(Error::Unsupported)?;
        // Same pinned B=6 insertion-only BTree node bound as budget::arena,
        // now with u32 keys and no value payload. All pending IDs are unique.
        let set_bytes = (1 + self.entities.len() / 5)
            .saturating_mul(20 * size_of::<usize>() + 11 * size_of::<u32>());
        let pending_bytes = self.entities.len().saturating_mul(size_of::<u32>());
        let scratch = set_bytes.saturating_add(2 * pending_bytes);
        if self.bytes.saturating_add(scratch) > options.retained_budget_bytes {
            return Err(Error::Limit);
        }
        let shell = Shell::new(&self.entities, 1e-7, false);
        let faces = shell
            .member_faces(member.local_solid_entity())
            .map_err(Error::from)?;
        let mut selected = BTreeSet::new();
        let mut pending = Vec::new();
        pending
            .try_reserve_exact(self.entities.len())
            .map_err(|_| Error::Limit)?;
        for id in faces {
            w.step(1)?;
            if selected.insert(id) {
                pending.push(id);
            }
        }
        while let Some(id) = pending.pop() {
            w.step(1)?;
            let records = self.entities.get(&id).ok_or(Error::Invalid)?;
            for e in records {
                walk(&self.entities, &e.args, &mut selected, &mut pending, &mut w)?;
            }
        }
        drop(pending);
        let child = Options {
            retained_budget_bytes: options
                .retained_budget_bytes
                .checked_sub(set_bytes)
                .ok_or(Error::Limit)?,
            ..options
        };
        decode(self, member, Some(&selected), child, interrupt, w.used)
    }
}
fn walk(
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    values: &[Value<'_>],
    selected: &mut BTreeSet<u32>,
    pending: &mut Vec<u32>,
    w: &mut Ledger<'_>,
) -> Result<()> {
    for v in values {
        w.step(1)?;
        match v {
            Value::Reference(id) => {
                // Membership is checked BEFORE growing either collection, so
                // their total unique IDs never exceed the preflight arena size.
                if !entities.contains_key(id) {
                    return Err(Error::Invalid);
                }
                if selected.insert(*id) {
                    pending.push(*id);
                }
            }
            Value::List(v) => walk(entities, v, selected, pending, w)?,
            Value::Typed(e) => walk(entities, &e.args, selected, pending, w)?,
            _ => {}
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
