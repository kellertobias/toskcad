//! Shared document storage and selected-member queries have separate budgets.
//! The original input (at most 64MiB) and full AST are charged ONCE to the
//! document. A live packet/view is additional storage/work: its receipt does
//! not include the document and must never be reported as the combined peak.
use super::*;

#[derive(Clone, Copy, Debug)]
pub struct DocumentOptions {
    pub maximum_work: usize,
    pub retained_budget_bytes: usize,
}
impl Default for DocumentOptions {
    fn default() -> Self {
        Self {
            maximum_work: 512_000_000,
            retained_budget_bytes: 512 * 1024 * 1024,
        }
    }
}
/// One validated representation arena. No solid, source recipe, or mass claim.
pub struct CompoundDocument<'input> {
    source: File<'input>,
}
/// Selected original records, with immutable original representation ownership.
/// Borrowing the document keeps its arena/input authority alive independently
/// of the bounded packet AST. No geometry or numeric leaf is reconstructed.
pub struct MemberPacket<'document, 'input> {
    document: &'document CompoundDocument<'input>,
    parsed: Parsed<'input>,
    bytes: usize,
}
/// The ACTUAL selected source-local token graph, with original represented
/// member and optional MappingLink retained. This is not a world-space solid.
/// Mapped coordinates are deliberately NOT applied or silently discarded;
/// callers must separately earn actual operator/source-origin transport.
pub struct SelectedSourceDocument<'document, 'input> {
    document: &'document CompoundDocument<'input>,
    member: Member,
    text: String,
    work: usize,
    bytes: usize,
}
impl<'document, 'input> SelectedSourceDocument<'document, 'input> {
    pub fn document(&self) -> &'document CompoundDocument<'input> {
        self.document
    }
    pub fn member(&self) -> Member {
        self.member
    }
    pub fn mapping_link(&self) -> Option<MappingLink> {
        self.member.mapping_link()
    }
    /// Ordinary singleton wrapper around the original LOCAL solid closure.
    /// Numeric/ref/component tokens and original IDs remain unchanged.
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Additional output/closure scratch/metadata peak; the borrowed document
    /// has its own live retained-byte scope and is NOT included again here.
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'document, 'input> MemberPacket<'document, 'input> {
    pub fn document(&self) -> &'document CompoundDocument<'input> {
        self.document
    }
    pub fn parsed(&self) -> &Parsed<'input> {
        &self.parsed
    }
    pub fn work_used(&self) -> usize {
        self.parsed.work_used()
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'input> CompoundDocument<'input> {
    pub fn parse(
        input: &'input str,
        options: DocumentOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        if options.maximum_work == 0
            || options.maximum_work > 512_000_000
            || options.retained_budget_bytes == 0
            || options.retained_budget_bytes > 512 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        Ok(Self {
            source: File::parse_scoped(
                input,
                options.maximum_work,
                options.retained_budget_bytes,
                interrupt,
            )?,
        })
    }
    pub fn represented_members(&self) -> &[Member] {
        self.source.represented_members()
    }
    pub fn literal_input(&self) -> &'input str {
        self.source.input
    }
    /// Immutable neutral records for another bounded decoder in this module.
    /// This supplies the original arena, never a geometry-admission receipt.
    pub(in crate::step_read) fn records(&self) -> &BTreeMap<u32, Vec<Entity<'input>>> {
        &self.source.entities
    }
    pub fn work_used(&self) -> usize {
        self.source.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.source.bytes
    }
    /// Local ordinary-member graph view. Its query ledger starts at zero;
    /// global document work/storage is separately reported by this document.
    pub fn local_member_document<'a>(
        &'a self,
        item: u32,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MemberDocument<'a, 'input>> {
        self.source
            .member_document_scoped(item, options, interrupt, 0, 0)
    }
    /// Select actual source-local geometry, including represented MAPPED_ITEM
    /// instances, without flattening or rebuilding any geometry. This typed
    /// result retains actual Member+MappingLink and borrows this document.
    /// The ordinary local_member_document API continues to refuse mappings.
    pub fn local_source_document<'a>(
        &'a self,
        item: u32,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<SelectedSourceDocument<'a, 'input>> {
        let data = self.source.source_document_data_scoped(
            item,
            options,
            interrupt,
            size_of::<SelectedSourceDocument<'_, '_>>(),
        )?;
        Ok(SelectedSourceDocument {
            document: self,
            member: data.member,
            text: data.text,
            work: data.work,
            bytes: data.bytes,
        })
    }
    /// Decode ONLY the selected actual dependency closure. Original IDs,
    /// component/numeric slices, source origin and map/operator are retained.
    /// The unchanged 16M/64MiB packet limits include closure scratch, copied AST
    /// and decoded arrays; document storage/work is an additional live scope.
    pub fn member_packet<'a>(
        &'a self,
        item: u32,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MemberPacket<'a, 'input>> {
        valid(options)?;
        let mut w = Ledger {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        w.step(self.source.members.len())?;
        let member = *self
            .source
            .members
            .iter()
            .find(|m| m.item_entity() == item)
            .ok_or(Error::Unsupported)?;
        let count = self.source.entities.len();
        let set_bytes =
            (1 + count / 5).saturating_mul(20 * size_of::<usize>() + 11 * size_of::<u32>());
        let pending_bytes = count.saturating_mul(size_of::<u32>());
        let fixed = size_of::<File<'_>>()
            + size_of::<MemberPacket<'_, '_>>()
            + size_of::<BTreeMap<u32, Vec<Entity<'_>>>>()
            + 2 * size_of::<usize>();
        admit_bytes(
            fixed
                .saturating_add(set_bytes)
                .saturating_add(pending_bytes),
            options.retained_budget_bytes,
        )?;
        let mut selected = BTreeSet::new();
        let mut pending = Vec::new();
        pending.try_reserve_exact(count).map_err(|_| Error::Limit)?;
        let pending_bytes = pending.capacity().saturating_mul(size_of::<u32>());
        // Pinned rust58cf8f9c raw_vec::grow_exact sets cap=len+additional,
        // and set_ptr_and_cap stores that requested count (including Global).
        // For this new non-ZST Vec, capacity=count: the preflight above already
        // admits the allocation. Keep measured capacity in the peak receipt.
        selected.insert(item);
        pending.push(item);
        while let Some(id) = pending.pop() {
            w.step(1)?;
            let records = self
                .source
                .entities
                .get(&id)
                .expect("validated member or checked dependency before enqueue");
            for e in records {
                walk(
                    &self.source.entities,
                    &e.args,
                    &mut selected,
                    &mut pending,
                    &mut w,
                )?;
            }
        }
        drop(pending);
        // Original nested Vec capacities are an upper bound for exact-length
        // copies. Preflight the selected AST before allocating any copied node.
        let mut bytes = fixed.saturating_add((1 + selected.len() / 5).saturating_mul(
            20 * size_of::<usize>() + 11 * (size_of::<u32>() + size_of::<Vec<Entity<'_>>>()),
        ));
        for id in &selected {
            bytes = bytes.saturating_add(record_bytes(&self.source.entities[id], &mut w)?);
        }
        admit_bytes(
            bytes.saturating_add(set_bytes),
            options.retained_budget_bytes,
        )?;
        let mut entities = BTreeMap::new();
        for id in &selected {
            entities.insert(*id, copy_records(&self.source.entities[id], &mut w)?);
        }
        // On the same pinned RawVec contract every copied Vec has capacity=len
        // <= the original capacity charged above. Map node count is unchanged.
        // Thus measured<=preflight bytes; no second fallible admission is needed.
        // See rust-lang/rust58cf8f9c35b355199727627d07ae34a4e2846d67
        // library/alloc/src/raw_vec/mod.rs::grow_exact/set_ptr_and_cap.
        let measured = budget::arena("", &entities, &mut w)?.saturating_add(fixed);
        let bytes = bytes.max(measured);
        drop(selected);
        let packet = File {
            input: self.source.input,
            entities: Arc::new(entities),
            solids: vec![],
            members: vec![],
            view: self.source.view,
            work: 0,
            bytes,
        };
        let parsed = decode(&packet, member, None, options, interrupt, w.used)?;
        let bytes = parsed
            .retained_budget_bytes()
            .max(bytes.saturating_add(set_bytes))
            .max(
                fixed
                    .saturating_add(set_bytes)
                    .saturating_add(pending_bytes),
            );
        Ok(MemberPacket {
            document: self,
            parsed,
            bytes,
        })
    }
}
fn admit_bytes(bytes: usize, maximum: usize) -> Result<()> {
    if bytes > maximum {
        return Err(Error::Limit);
    }
    Ok(())
}
fn record_bytes(records: &Vec<Entity<'_>>, w: &mut Ledger<'_>) -> Result<usize> {
    let mut n = records.capacity().saturating_mul(size_of::<Entity<'_>>());
    for e in records {
        n = n.saturating_add(value_bytes(&e.args, w)?);
    }
    Ok(n)
}
fn value_bytes(values: &Vec<Value<'_>>, w: &mut Ledger<'_>) -> Result<usize> {
    let mut n = values.capacity().saturating_mul(size_of::<Value<'_>>());
    for v in values {
        w.step(1)?;
        n = n.saturating_add(match v {
            Value::List(v) => value_bytes(v, w)?,
            Value::Typed(e) => value_bytes(&e.args, w)?,
            _ => 0,
        });
    }
    Ok(n)
}
fn copy_records<'a>(records: &[Entity<'a>], w: &mut Ledger<'_>) -> Result<Vec<Entity<'a>>> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(records.len())
        .map_err(|_| Error::Limit)?;
    for e in records {
        w.step(1)?;
        result.push(Entity {
            source: e.source,
            name: e.name,
            args: copy_values(&e.args, w)?,
        });
    }
    Ok(result)
}
fn copy_values<'a>(values: &[Value<'a>], w: &mut Ledger<'_>) -> Result<Vec<Value<'a>>> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(values.len())
        .map_err(|_| Error::Limit)?;
    for v in values {
        w.step(1)?;
        result.push(match v {
            Value::String(s) => Value::String(s),
            Value::Enum(s) => Value::Enum(s),
            Value::Reference(id) => Value::Reference(*id),
            Value::Number(n) => Value::Number(*n),
            Value::List(v) => Value::List(copy_values(v, w)?),
            Value::Omitted => Value::Omitted,
            Value::Derived => Value::Derived,
            Value::Typed(e) => Value::Typed(Entity {
                source: e.source,
                name: e.name,
                args: copy_values(&e.args, w)?,
            }),
        });
    }
    Ok(result)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod source_tests;
