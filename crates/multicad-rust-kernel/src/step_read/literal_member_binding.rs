//! Exact local-record correspondence to an independently qualified candidate.
//! This proves neither occupied material nor placement: callers earn those
//! separately. Original numeric spellings, sharing and every reachable owner
//! survive; only entity numbering may change.
use super::{owned_insert::CompoundDocument, part21::Entity};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Different,
    Limit,
    WorkLimit,
    StorageLimit,
    EntityLimit,
    Cancelled,
}
#[derive(Clone, Copy)]
pub struct Options {
    pub maximum_work: usize,
    /// Additional query storage; the two original document arenas are separate.
    pub maximum_bytes: usize,
    pub maximum_entities: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            maximum_entities: 65_536,
        }
    }
}
pub struct Binding {
    entities: usize,
    work: usize,
    peak_bytes: usize,
}
impl Binding {
    pub fn matched_entities(&self) -> usize {
        self.entities
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_peak_bytes(&self) -> usize {
        self.peak_bytes
    }
}
type Result<T> = std::result::Result<T, Error>;
// Conservative allocation reserve per distinct pair: two BTreeMap nodes and
// the Vec queue (including capacity growth). Even a one-entry root fits here.
const PAIR_BYTES: usize = 1024;
const FIXED_BYTES: usize = 4096;
struct State<'a> {
    options: Options,
    interrupt: Option<&'a dyn Fn() -> bool>,
    work: usize,
    peak: usize,
    forward: BTreeMap<u32, u32>,
    reverse: BTreeMap<u32, u32>,
    pending: Vec<(u32, u32)>,
}
impl State<'_> {
    fn charge(&mut self, amount: usize) -> Result<()> {
        if self.interrupt.is_some_and(|cancel| cancel()) {
            return Err(Error::Cancelled);
        }
        self.work = self.work.checked_add(amount).ok_or(Error::Limit)?;
        if self.work > self.options.maximum_work {
            return Err(Error::WorkLimit);
        }
        Ok(())
    }
    fn reserve(&mut self, pairs: usize, scratch: usize) -> Result<()> {
        let bytes = pairs
            .checked_mul(PAIR_BYTES)
            .and_then(|n| n.checked_add(FIXED_BYTES))
            .and_then(|n| n.checked_add(scratch))
            .ok_or(Error::Limit)?;
        if bytes > self.options.maximum_bytes {
            return Err(Error::StorageLimit);
        }
        self.peak = self.peak.max(bytes);
        Ok(())
    }
    fn reference(&mut self, a: u32, b: u32) -> Result<()> {
        self.charge(64 * (1 + (self.forward.len() + 1).ilog2() as usize))?;
        if let Some(&expected) = self.forward.get(&a) {
            return if expected == b {
                Ok(())
            } else {
                Err(Error::Different)
            };
        }
        if self.reverse.contains_key(&b) {
            return Err(Error::Different);
        }
        let count = self.forward.len() + 1;
        if count > self.options.maximum_entities {
            return Err(Error::EntityLimit);
        }
        self.reserve(count, 0)?;
        self.forward.insert(a, b);
        self.reverse.insert(b, a);
        self.pending.push((a, b));
        Ok(())
    }
    fn entity(&mut self, a: &Entity<'_>, b: &Entity<'_>) -> Result<()> {
        component(
            self,
            a.source.expect("validated original component"),
            b.source.expect("validated original component"),
        )
    }
}
trait Visitor {
    fn charge(&mut self, amount: usize) -> Result<()>;
    fn reference(&mut self, a: u32, b: u32) -> Result<()>;
}
impl Visitor for State<'_> {
    fn charge(&mut self, amount: usize) -> Result<()> {
        State::charge(self, amount)
    }
    fn reference(&mut self, a: u32, b: u32) -> Result<()> {
        State::reference(self, a, b)
    }
}
fn component(state: &mut impl Visitor, mut a: &str, mut b: &str) -> Result<()> {
    state.charge(a.len().checked_add(b.len()).ok_or(Error::Limit)?)?;
    loop {
        state.charge(1)?;
        match (token(&mut a), token(&mut b)) {
            (None, None) => return Ok(()),
            (Some(Token::Reference(a)), Some(Token::Reference(b))) => state.reference(a, b)?,
            (Some(Token::Literal(a)), Some(Token::Literal(b))) if a == b => {}
            _ => return Err(Error::Different),
        }
    }
}
pub(crate) enum Token<'a> {
    Literal(&'a str),
    Reference(u32),
}
// This is a zero-allocation comparison cursor over components ALREADY parsed
// by the bounded Part21 parser, not an alternative untrusted-input parser.
// Original numeric/quoted tokens are compared literally; only #IDs are mapped.
pub(crate) fn token<'a>(rest: &mut &'a str) -> Option<Token<'a>> {
    loop {
        *rest = rest.trim_start_matches(|c: char| c.is_ascii_whitespace());
        if let Some(comment) = rest.strip_prefix("/*") {
            *rest = &comment[comment.find("*/").expect("validated comment") + 2..];
        } else {
            break;
        }
    }
    let first = *rest.as_bytes().first()?;
    let bytes = rest.as_bytes();
    let end = match first {
        b'\'' => {
            let mut i = 1;
            loop {
                if bytes[i] == b'\'' {
                    if bytes.get(i + 1) == Some(&b'\'') {
                        i += 2;
                    } else {
                        break i + 1;
                    }
                } else {
                    i += 1;
                }
            }
        }
        b'#' => 1 + bytes[1..].iter().take_while(|c| c.is_ascii_digit()).count(),
        b'+' | b'-' | b'0'..=b'9' => bytes
            .iter()
            .take_while(|c| matches!(c, b'+' | b'-' | b'.' | b'E' | b'0'..=b'9'))
            .count(),
        b'A'..=b'Z' | b'_' => bytes
            .iter()
            .take_while(|c| c.is_ascii_alphanumeric() || **c == b'_')
            .count(),
        b'.' => 2 + rest[1..].find('.').expect("validated enumeration"),
        _ => 1,
    };
    let value = &rest[..end];
    *rest = &rest[end..];
    Some(if first == b'#' {
        Token::Reference(value[1..].parse().expect("validated reference"))
    } else {
        Token::Literal(value)
    })
}
/// Bind only the selected members' LOCAL reachable geometry records. Each item
/// must belong to its supplied validated representation document. Mapped
/// occurrences intentionally require a separate placement proof.
pub fn bind_local_members(
    expected: &CompoundDocument<'_>,
    expected_item: u32,
    actual: &CompoundDocument<'_>,
    actual_item: u32,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Binding> {
    if options.maximum_work == 0
        || options.maximum_work > 512_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
        || options.maximum_entities == 0
        || options.maximum_entities > 1_048_576
    {
        return Err(Error::Invalid);
    }
    let mut state = State {
        options,
        interrupt,
        work: 0,
        peak: 0,
        forward: BTreeMap::new(),
        reverse: BTreeMap::new(),
        pending: Vec::new(),
    };
    state.charge(expected.represented_members().len() + actual.represented_members().len() + 1)?;
    let a = expected
        .represented_members()
        .iter()
        .find(|m| m.item_entity() == expected_item)
        .ok_or(Error::Invalid)?;
    let b = actual
        .represented_members()
        .iter()
        .find(|m| m.item_entity() == actual_item)
        .ok_or(Error::Invalid)?;
    state.reference(a.local_solid_entity(), b.local_solid_entity())?;
    while let Some((a, b)) = state.pending.pop() {
        state.charge(
            128 * (1 + expected
                .records()
                .len()
                .max(actual.records().len())
                .max(1)
                .ilog2() as usize),
        )?;
        let a = expected.records().get(&a).ok_or(Error::Invalid)?;
        let b = actual.records().get(&b).ok_or(Error::Invalid)?;
        if a.len() != b.len() {
            return Err(Error::Different);
        }
        for (a, b) in a.iter().zip(b) {
            state.entity(a, b)?;
        }
    }
    state.charge(1)?;
    Ok(Binding {
        entities: state.forward.len(),
        work: state.work,
        peak_bytes: state.peak,
    })
}

#[cfg(test)]
mod tests;

/// Immutable components from an already qualified retained arena. This is an
/// internal exact-record view, never an alternative geometry admission API.
pub(crate) trait RetainedRecords {
    fn record_count(&self) -> usize;
    fn record_index(&self, id: u32) -> Option<usize>;
    fn component(&self, record: usize, index: usize) -> Option<&str>;
}
struct Dense<'a, A: RetainedRecords, B: RetainedRecords> {
    a: &'a A,
    b: &'a B,
    options: Options,
    interrupt: Option<&'a dyn Fn() -> bool>,
    work: usize,
    forward: Vec<u32>,
    reverse: Vec<u32>,
    pending: Vec<(usize, usize)>,
    matched: usize,
}
impl<A: RetainedRecords, B: RetainedRecords> Visitor for Dense<'_, A, B> {
    fn charge(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.work = self.work.checked_add(n).ok_or(Error::Limit)?;
        if self.work > self.options.maximum_work {
            return Err(Error::WorkLimit);
        }
        Ok(())
    }
    fn reference(&mut self, a: u32, b: u32) -> Result<()> {
        // Two sorted-index lookups; counts<=2^20 bound each binary search.
        // Use the same32-unit lookup convention as the selected NPT reader.
        self.charge(64)?;
        let a = self.a.record_index(a).ok_or(Error::Invalid)?;
        let b = self.b.record_index(b).ok_or(Error::Invalid)?;
        if self.forward[a] != u32::MAX {
            return if self.forward[a] as usize == b {
                Ok(())
            } else {
                Err(Error::Different)
            };
        }
        if self.reverse[b] != u32::MAX {
            return Err(Error::Different);
        }
        // Counts were bounded by maximum_entities<=2^20 before allocation.
        self.forward[a] = u32::try_from(b).expect("bounded retained record index");
        self.reverse[b] = u32::try_from(a).expect("bounded retained record index");
        self.pending.push((a, b));
        self.matched += 1;
        Ok(())
    }
}
fn allocate<T>(count: usize, used: &mut usize, limit: usize) -> Result<Vec<T>> {
    let requested = count
        .checked_mul(std::mem::size_of::<T>())
        .ok_or(Error::StorageLimit)?;
    if requested > limit.saturating_sub(*used) {
        return Err(Error::StorageLimit);
    }
    let mut out = Vec::new();
    out.try_reserve_exact(count)
        .map_err(|_| Error::StorageLimit)?;
    *used = used
        .checked_add(out.capacity() * std::mem::size_of::<T>())
        .ok_or(Error::StorageLimit)?;
    if *used > limit {
        return Err(Error::StorageLimit);
    }
    Ok(out)
}
pub(crate) fn bind_retained<A: RetainedRecords, B: RetainedRecords>(
    a: &A,
    a_root: u32,
    b: &B,
    b_root: u32,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Binding> {
    if options.maximum_work == 0
        || options.maximum_work > 512_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
        || options.maximum_entities == 0
        || options.maximum_entities > 1_048_576
    {
        return Err(Error::Invalid);
    }
    if a.record_count() > options.maximum_entities || b.record_count() > options.maximum_entities {
        return Err(Error::EntityLimit);
    }
    let mut bytes = FIXED_BYTES;
    let mut forward = allocate(a.record_count(), &mut bytes, options.maximum_bytes)?;
    forward.resize(a.record_count(), u32::MAX);
    let mut reverse = allocate(b.record_count(), &mut bytes, options.maximum_bytes)?;
    reverse.resize(b.record_count(), u32::MAX);
    let pending = allocate(
        a.record_count().min(b.record_count()),
        &mut bytes,
        options.maximum_bytes,
    )?;
    let mut state = Dense {
        a,
        b,
        options,
        interrupt,
        work: 0,
        forward,
        reverse,
        pending,
        matched: 0,
    };
    state.charge(a.record_count() + b.record_count() + 1)?;
    state.reference(a_root, b_root)?;
    while let Some((a, b)) = state.pending.pop() {
        state.charge(32)?;
        let mut index = 0;
        loop {
            match (state.a.component(a, index), state.b.component(b, index)) {
                (Some(a), Some(b)) => component(&mut state, a, b)?,
                (None, None) if index > 0 => break,
                _ => return Err(Error::Different),
            }
            index += 1;
        }
    }
    Ok(Binding {
        entities: state.matched,
        work: state.work,
        peak_bytes: bytes,
    })
}
struct DocumentRecords<'a, 'input> {
    ids: &'a [u32],
    records: Vec<&'a [Entity<'input>]>,
}
impl RetainedRecords for DocumentRecords<'_, '_> {
    fn record_count(&self) -> usize {
        self.ids.len()
    }
    fn record_index(&self, id: u32) -> Option<usize> {
        self.ids.binary_search(&id).ok()
    }
    fn component(&self, record: usize, index: usize) -> Option<&str> {
        self.records.get(record)?.get(index)?.source
    }
}
/// Bind an already retained original to independently parsed selected file
/// records. Selection IDs come from that file's actual dependency traversal.
/// The neutral document and original retained arena are separate resident scopes.
pub(crate) fn bind_owned_document<A: RetainedRecords>(
    a: &A,
    a_root: u32,
    document: &CompoundDocument<'_>,
    ids: &[u32],
    actual_root: u32,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Binding> {
    if ids.len() > options.maximum_entities {
        return Err(Error::EntityLimit);
    }
    if !ids.windows(2).all(|w| w[0] < w[1]) {
        return Err(Error::Invalid);
    }
    let mut bytes = 0;
    let mut records = allocate(ids.len(), &mut bytes, options.maximum_bytes)?;
    let mut work = ids.len().checked_mul(32).ok_or(Error::WorkLimit)?;
    if work >= options.maximum_work {
        return Err(Error::WorkLimit);
    }
    for id in ids {
        if interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        records.push(document.records().get(id).ok_or(Error::Invalid)?.as_slice());
    }
    let b = DocumentRecords { ids, records };
    let result = bind_retained(
        a,
        a_root,
        &b,
        actual_root,
        Options {
            maximum_work: options.maximum_work - work,
            maximum_bytes: options
                .maximum_bytes
                .checked_sub(bytes)
                .ok_or(Error::StorageLimit)?,
            ..options
        },
        interrupt,
    )?;
    work += result.work;
    Ok(Binding {
        entities: result.entities,
        work,
        peak_bytes: bytes + result.peak_bytes,
    })
}
