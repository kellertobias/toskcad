//! Neutral represented topology ownership. Material, supports, finite trims
//! and occurrence operators remain independent per-member prerequisites.
use super::*;
use crate::step_read::owned_insert::CompoundDocument;
const WORK: usize = 16_000_000;
const BYTES: usize = 64 * 1024 * 1024;
#[derive(Clone, Copy)]
pub struct DocumentOwnershipOptions {
    pub maximum_work: usize,
    pub maximum_query_bytes: usize,
    pub maximum_members: usize,
}
impl Default for DocumentOwnershipOptions {
    fn default() -> Self {
        Self {
            maximum_work: WORK,
            maximum_query_bytes: BYTES,
            maximum_members: 1024,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DocumentOwnershipError {
    Invalid,
    Unsupported,
    Limit,
    Cancelled,
}
impl From<ReadError> for DocumentOwnershipError {
    fn from(e: ReadError) -> Self {
        match e {
            ReadError::Invalid => Self::Invalid,
            ReadError::Unsupported => Self::Unsupported,
            ReadError::Limit => Self::Limit,
        }
    }
}
pub struct DocumentOwnership<'d, 'i> {
    document: &'d CompoundDocument<'i>,
    topology: usize,
    work: usize,
    peak: usize,
}
impl<'d, 'i> DocumentOwnership<'d, 'i> {
    pub fn document(&self) -> &'d CompoundDocument<'i> {
        self.document
    }
    pub fn represented_member_count(&self) -> usize {
        self.document.represented_members().len()
    }
    pub fn topology_entity_count(&self) -> usize {
        self.topology
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_peak_bytes(&self) -> usize {
        self.peak
    }
}
struct Budget<'a> {
    used: usize,
    maximum: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Budget<'_> {
    fn step(&mut self, n: usize) -> std::result::Result<(), DocumentOwnershipError> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(DocumentOwnershipError::Cancelled);
        }
        self.used = self
            .used
            .checked_add(n)
            .filter(|n| *n <= self.maximum)
            .ok_or(DocumentOwnershipError::Limit)?;
        Ok(())
    }
}
/// Every actual local BREP has one exclusive topological owner, with no orphan
/// topology. Repeated properly represented maps of that SAME local BREP share
/// its local owner; this receipt does not validate any mapping or occupied
/// material. Input/neutral arena remain a separate borrowed resident scope.
/// Work units cover component/tag visits, direct topology references, and
/// bounded64-unit map lookups; they are not CPU instruction counts.
pub fn validate_document_topology_ownership<'d, 'i>(
    document: &'d CompoundDocument<'i>,
    tolerance: f64,
    options: DocumentOwnershipOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> std::result::Result<DocumentOwnership<'d, 'i>, DocumentOwnershipError> {
    let mut budget = Budget {
        used: 0,
        maximum: options.maximum_work,
        interrupt,
    };
    budget.step(1)?;
    if !tolerance.is_finite()
        || tolerance <= 0.
        || options.maximum_work > WORK
        || options.maximum_query_bytes == 0
        || options.maximum_query_bytes > BYTES
        || options.maximum_members == 0
        || options.maximum_members > 1024
    {
        return Err(DocumentOwnershipError::Invalid);
    }
    if document.represented_members().len() > options.maximum_members {
        return Err(DocumentOwnershipError::Limit);
    }
    for (i, member) in document.represented_members().iter().enumerate() {
        budget.step(i + 1)?;
        if document.represented_members()[..i]
            .iter()
            .any(|prior| prior.item_entity() == member.item_entity())
        {
            return Err(DocumentOwnershipError::Invalid);
        }
    }
    let mut references = 0usize;
    for parts in document.records().values() {
        for entity in parts {
            budget.step(32)?;
            if topology(entity.name) {
                for arg in &entity.args {
                    budget.step(1)?;
                    match arg {
                        Value::Reference(_) => {
                            references = references
                                .checked_add(1)
                                .ok_or(DocumentOwnershipError::Limit)?
                        }
                        Value::List(items) => {
                            for item in items {
                                budget.step(1)?;
                                references += usize::from(matches!(item, Value::Reference(_)));
                            }
                        }
                        _ => {}
                    }
                }
            }
        }
    }
    let (topology_count, _) = inspect_literal_entities(document.records().values().flatten())
        .map_err(DocumentOwnershipError::from)?;
    let slots = references
        .checked_add(1)
        .ok_or(DocumentOwnershipError::Limit)?;
    // A u32/u32 BTree node is <256bytes on the pinned64-bit toolchain, and
    // there are fewer than two nodes per inserted pair.1KiB/pair includes
    // sparse-root/split overlap. The exact-reserved queue bounds ALL pushes:
    // each owned topology node expands once; repeated owners stop immediately.
    let mut peak = 16 * 1024
        + topology_count
            .checked_mul(1024)
            .ok_or(DocumentOwnershipError::Limit)?
        + slots
            .checked_mul(std::mem::size_of::<u32>())
            .ok_or(DocumentOwnershipError::Limit)?;
    if peak > options.maximum_query_bytes {
        return Err(DocumentOwnershipError::Limit);
    }
    let mut pending = Vec::new();
    pending
        .try_reserve_exact(slots)
        .map_err(|_| DocumentOwnershipError::Limit)?;
    peak = peak
        .max(16 * 1024 + topology_count * 1024 + pending.capacity() * std::mem::size_of::<u32>());
    if peak > options.maximum_query_bytes {
        return Err(DocumentOwnershipError::Limit);
    }
    let mut owners = BTreeMap::new();
    let shell = Shell::new(document.records(), tolerance, true);
    for member in document.represented_members() {
        budget.step(1)?;
        claim_topology_by(
            &shell,
            member.local_solid_entity(),
            &mut owners,
            &mut pending,
            &mut || budget.step(128),
        )?;
    }
    if owners.len() != topology_count {
        return Err(DocumentOwnershipError::Invalid);
    }
    budget.step(1)?;
    Ok(DocumentOwnership {
        document,
        topology: topology_count,
        work: budget.used,
        peak,
    })
}
#[cfg(test)]
mod tests;
