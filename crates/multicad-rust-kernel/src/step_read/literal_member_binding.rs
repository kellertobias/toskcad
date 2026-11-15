//! Exact local-record correspondence to an independently qualified candidate.
//! This proves neither occupied material nor placement: callers earn those
//! separately. Original numeric spellings, sharing and every reachable owner
//! survive; only entity numbering may change.
use super::{
    owned_insert::CompoundDocument,
    part21::{self, Entity, Value},
};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Different,
    Limit,
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
            return Err(Error::Limit);
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
            return Err(Error::Limit);
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
            return Err(Error::Limit);
        }
        self.reserve(count, 0)?;
        self.forward.insert(a, b);
        self.reverse.insert(b, a);
        self.pending.push((a, b));
        Ok(())
    }
    fn entity(&mut self, a: &Entity<'_>, b: &Entity<'_>, numbers: bool) -> Result<()> {
        self.charge(1 + a.name.len() + b.name.len())?;
        if a.name != b.name || a.args.len() != b.args.len() {
            return Err(Error::Different);
        }
        if numbers {
            let sources = a
                .source
                .ok_or(Error::Invalid)?
                .len()
                .checked_add(b.source.ok_or(Error::Invalid)?.len())
                .ok_or(Error::Limit)?;
            // A numeric-token slice occupies two machine words. Reserve twice
            // that per input byte before either temporary Vec is allocated.
            let scratch = sources.checked_mul(32).ok_or(Error::Limit)?;
            self.reserve(self.forward.len(), scratch)?;
            self.charge(sources.checked_mul(8).ok_or(Error::Limit)?)?;
            let x = part21::numeric_spellings(a).ok_or(Error::Invalid)?;
            let y = part21::numeric_spellings(b).ok_or(Error::Invalid)?;
            if x != y {
                return Err(Error::Different);
            }
        }
        for (a, b) in a.args.iter().zip(&b.args) {
            self.value(a, b)?;
        }
        Ok(())
    }
    fn value(&mut self, a: &Value<'_>, b: &Value<'_>) -> Result<()> {
        self.charge(1)?;
        match (a, b) {
            (Value::String(a), Value::String(b)) | (Value::Enum(a), Value::Enum(b)) => {
                self.charge(a.len() + b.len())?;
                if a != b {
                    return Err(Error::Different);
                }
            }
            (Value::Reference(a), Value::Reference(b)) => self.reference(*a, *b)?,
            // The complete original number-token sequence was compared above.
            // Comparing float caches here would accept distinct decimal leaves.
            (Value::Number(_), Value::Number(_))
            | (Value::Omitted, Value::Omitted)
            | (Value::Derived, Value::Derived) => {}
            (Value::List(a), Value::List(b)) => {
                if a.len() != b.len() {
                    return Err(Error::Different);
                }
                for (a, b) in a.iter().zip(b) {
                    self.value(a, b)?;
                }
            }
            (Value::Typed(a), Value::Typed(b)) => self.entity(a, b, false)?,
            _ => return Err(Error::Different),
        }
        Ok(())
    }
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
        || options.maximum_work > 64_000_000
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
            state.entity(a, b, true)?;
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
