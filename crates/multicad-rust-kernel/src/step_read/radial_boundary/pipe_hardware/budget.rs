//! Allocation-free grammar preflight for this codec's cold parser phase.
//! It grants no geometric or document admission. The real neutral parser and
//! represented-shell decoder still process the unchanged original file.
use crate::{
    axial_pipe_thread::Error,
    step_read::part21::{Entity, Value},
};
use std::mem::size_of;
type Result<T> = std::result::Result<T, Error>;
const TOKEN: usize = 64 * 1024;
const VALUES: usize = 8_388_608;
const NODES: usize = 1_048_576;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug)]
pub(crate) struct Receipt {
    pub input_bytes: usize,
    pub represented_record_count: usize,
    pub values: usize,
    pub references: usize,
    pub components: usize,
    pub maximum_argument_depth: usize,
    pub retained_arena_bytes_reserved: usize,
    pub parser_peak_bytes_reserved: usize,
    pub scan_work_used: usize,
    pub delegated_parser_work_reserved: usize,
    pub aggregate_work_reserved: usize,
}
pub(crate) fn reserve(
    input: &str,
    maximum_work: usize,
    maximum_bytes: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Receipt> {
    if maximum_work == 0
        || maximum_work > 64_000_000
        || maximum_bytes == 0
        || maximum_bytes > 64 * 1024 * 1024
    {
        return Err(Error::Invalid);
    }
    if input.len() > 64 * 1024 * 1024 || input.len() > maximum_bytes || input.len() > maximum_work {
        return Err(Error::Limit);
    }
    let mut p = Scan::new(input, maximum_work, interrupt);
    p.keyword("ISO-10303-21")?;
    p.keyword("HEADER")?;
    for required in ["FILE_DESCRIPTION", "FILE_NAME", "FILE_SCHEMA"] {
        if p.entity(0)? != required {
            return Err(Error::Invalid);
        }
        p.require(b';')?;
    }
    p.trivia()?;
    while !p.rest().starts_with("ENDSEC") {
        if matches!(
            p.entity(0)?,
            "FILE_DESCRIPTION" | "FILE_NAME" | "FILE_SCHEMA"
        ) {
            return Err(Error::Invalid);
        }
        p.require(b';')?;
        p.trivia()?;
    }
    p.keyword("ENDSEC")?;
    p.trivia()?;
    if !p.rest().starts_with("DATA") {
        return Err(Error::Invalid);
    }
    p.advance(4)?;
    if p.peek()? == Some(b'(') {
        return Err(Error::Unsupported);
    }
    p.require(b';')?;
    while p.take(b'#')? {
        p.records = p.records.checked_add(1).ok_or(Error::Limit)?;
        if p.records > NODES {
            return Err(Error::Limit);
        }
        p.id()?;
        p.require(b'=')?;
        if p.take(b'(')? {
            let mut names: [Option<&str>; 64] = [None; 64];
            let mut count = 0usize;
            loop {
                if count == 64 {
                    return Err(Error::Limit);
                }
                let name = p.entity(0)?;
                for old in names[..count].iter().flatten() {
                    let old: &str = old;
                    let comparisons = if old.len() == name.len() {
                        old.len()
                    } else {
                        1
                    };
                    p.component_comparisons = p
                        .component_comparisons
                        .checked_add(comparisons)
                        .ok_or(Error::Limit)?;
                    p.charge(comparisons)?;
                    if old == name {
                        return Err(Error::Invalid);
                    }
                }
                names[count] = Some(name);
                count += 1;
                if p.take(b')')? {
                    break;
                }
            }
            let allocation = capacity(count)?
                .checked_mul(size_of::<Entity<'_>>())
                .ok_or(Error::Limit)?;
            p.components = p.components.checked_add(count).ok_or(Error::Limit)?;
            p.component_bytes = p
                .component_bytes
                .checked_add(allocation)
                .ok_or(Error::Limit)?;
            p.largest_vector_bytes = p.largest_vector_bytes.max(allocation);
        } else {
            p.entity(0)?;
            p.components = p.components.checked_add(1).ok_or(Error::Limit)?;
            p.component_bytes = p
                .component_bytes
                .checked_add(size_of::<Entity<'_>>())
                .ok_or(Error::Limit)?;
            p.largest_vector_bytes = p.largest_vector_bytes.max(size_of::<Entity<'_>>());
        }
        p.require(b';')?;
    }
    p.keyword("ENDSEC")?;
    p.trivia()?;
    if p.rest().starts_with("DATA") {
        return Err(Error::Unsupported);
    }
    p.keyword("END-ISO-10303-21")?;
    p.trivia()?;
    if p.offset != input.len() {
        return Err(Error::Invalid);
    }
    // Pinned1.99.1 insertion-only BTree B=6:11KV slots/12children;
    // >=5keys/nonroot. Extra64nodes cover all simultaneous insertion splits.
    let node = 20 * size_of::<usize>() + 11 * (size_of::<u32>() + size_of::<Vec<Entity<'_>>>());
    let nodes = 1usize.checked_add(p.records / 5).ok_or(Error::Limit)?;
    let tree = nodes.checked_mul(node).ok_or(Error::Limit)?;
    let vectors = p
        .value_bytes
        .checked_add(p.component_bytes)
        .ok_or(Error::Limit)?;
    let retained = input
        .len()
        .checked_add(vectors)
        .and_then(|n| n.checked_add(tree))
        .and_then(|n| n.checked_add(size_of::<std::collections::BTreeMap<u32, Vec<Entity<'_>>>>()))
        .ok_or(Error::Limit)?;
    // One Vec reallocates at a time: a nested value finishes before its parent
    // push. Final-arena storage plus the largest single final vector encloses
    // old/new buffer overlap. Header vectors remain conservatively counted.
    let peak = retained
        .checked_add(p.largest_vector_bytes)
        .and_then(|n| n.checked_add(64 * node + 64 * 1024))
        .ok_or(Error::Limit)?;
    #[cfg(test)]
    eprintln!(
        "COLD_RESERVATION bytes={} records={} values={} refs={} arena={} peak={} maxbytes={}",
        input.len(),
        p.records,
        p.values,
        p.references,
        retained,
        peak,
        maximum_bytes
    );
    if peak > maximum_bytes {
        return Err(Error::Limit);
    }
    // Separate delegated upper reservation: lexical/number visits, bounded
    // BTree comparisons and complex-record name comparisons. Not measured CPU.
    let delegated = input
        .len()
        .checked_mul(8)
        .and_then(|n| n.checked_add(p.values.checked_mul(16)?))
        .and_then(|n| n.checked_add(p.records.checked_mul(64)?))
        .and_then(|n| n.checked_add(p.component_comparisons))
        .ok_or(Error::Limit)?;
    let used = maximum_work - p.remaining;
    let aggregate = used.checked_add(delegated).ok_or(Error::Limit)?;
    #[cfg(test)]
    eprintln!(
        "COLD_WORK_RESERVATION scan={used} delegated={delegated} aggregate={aggregate} maximum={maximum_work}"
    );
    if aggregate > maximum_work {
        return Err(Error::Limit);
    }
    Ok(Receipt {
        input_bytes: input.len(),
        represented_record_count: p.records,
        values: p.values,
        references: p.references,
        components: p.components,
        maximum_argument_depth: p.depth,
        retained_arena_bytes_reserved: retained,
        parser_peak_bytes_reserved: peak,
        scan_work_used: used,
        delegated_parser_work_reserved: delegated,
        aggregate_work_reserved: aggregate,
    })
}
fn capacity(len: usize) -> Result<usize> {
    if len == 0 {
        Ok(0)
    } else {
        len.checked_next_power_of_two()
            .map(|n| n.max(4))
            .ok_or(Error::Limit)
    }
}
struct Scan<'a, 'c> {
    input: &'a str,
    offset: usize,
    remaining: usize,
    interrupt: Option<&'c dyn Fn() -> bool>,
    values: usize,
    references: usize,
    components: usize,
    records: usize,
    depth: usize,
    value_bytes: usize,
    component_bytes: usize,
    component_comparisons: usize,
    largest_vector_bytes: usize,
}
impl<'a, 'c> Scan<'a, 'c> {
    fn new(input: &'a str, maximum_work: usize, interrupt: Option<&'c dyn Fn() -> bool>) -> Self {
        Self {
            input,
            offset: 0,
            remaining: maximum_work,
            interrupt,
            values: 0,
            references: 0,
            components: 0,
            records: 0,
            depth: 0,
            value_bytes: 0,
            component_bytes: 0,
            component_comparisons: 0,
            largest_vector_bytes: 0,
        }
    }
    fn rest(&self) -> &'a str {
        &self.input[self.offset..]
    }
    fn charge(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        Ok(())
    }
    fn advance(&mut self, n: usize) -> Result<()> {
        let end = self.offset.checked_add(n).ok_or(Error::Limit)?;
        if end > self.input.len() {
            return Err(Error::Invalid);
        }
        while self.offset < end {
            let n = (end - self.offset).min(4096);
            self.charge(n)?;
            self.offset += n;
        }
        Ok(())
    }
    fn trivia(&mut self) -> Result<()> {
        loop {
            loop {
                let n = self
                    .rest()
                    .bytes()
                    .take(4096)
                    .take_while(u8::is_ascii_whitespace)
                    .count();
                self.advance(n)?;
                if n < 4096 {
                    break;
                }
            }
            if self.rest().starts_with("/*") {
                self.advance(2)?;
                loop {
                    let end = self.offset.saturating_add(4096).min(self.input.len());
                    let look = end.saturating_add(1).min(self.input.len());
                    let found = self.input.as_bytes()[self.offset..look]
                        .windows(2)
                        .position(|p| p == b"*/");
                    if let Some(n) = found {
                        self.advance(n + 2)?;
                        break;
                    }
                    if end == self.input.len() {
                        return Err(Error::Invalid);
                    }
                    self.advance(end - self.offset)?;
                }
            } else {
                return Ok(());
            }
        }
    }
    fn peek(&mut self) -> Result<Option<u8>> {
        self.trivia()?;
        Ok(self.input.as_bytes().get(self.offset).copied())
    }
    fn take(&mut self, b: u8) -> Result<bool> {
        if self.peek()? == Some(b) {
            self.advance(1)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }
    fn require(&mut self, b: u8) -> Result<()> {
        if self.take(b)? {
            Ok(())
        } else {
            Err(Error::Invalid)
        }
    }
    fn keyword(&mut self, s: &str) -> Result<()> {
        self.trivia()?;
        if !self.rest().starts_with(s) {
            return Err(Error::Invalid);
        }
        self.advance(s.len())?;
        self.require(b';')
    }
    fn token(&mut self, accept: impl Fn(u8) -> bool) -> Result<&'a str> {
        let start = self.offset;
        let n = self
            .rest()
            .bytes()
            .take_while(|&b| accept(b))
            .take(TOKEN + 1)
            .count();
        if n == 0 {
            return Err(Error::Invalid);
        }
        if n > TOKEN {
            return Err(Error::Limit);
        }
        self.advance(n)?;
        Ok(&self.input[start..self.offset])
    }
    fn name(&mut self) -> Result<&'a str> {
        if !self
            .peek()?
            .is_some_and(|b| b.is_ascii_uppercase() || b == b'_')
        {
            return Err(Error::Invalid);
        }
        self.token(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_')
    }
    fn id(&mut self) -> Result<()> {
        let id = self
            .token(|b| b.is_ascii_digit())?
            .parse::<u32>()
            .map_err(|_| Error::Limit)?;
        if id == 0 { Err(Error::Invalid) } else { Ok(()) }
    }
    fn quoted(&mut self) -> Result<()> {
        let start = self.offset;
        loop {
            let remaining = TOKEN.checked_sub(self.offset - start).ok_or(Error::Limit)?;
            let next = self
                .rest()
                .as_bytes()
                .iter()
                .take(remaining + 1)
                .position(|b| *b == b'\'');
            let end = match next {
                Some(n) => n,
                None if self.rest().len() > remaining => return Err(Error::Limit),
                None => return Err(Error::Invalid),
            };
            self.advance(end + 1)?;
            if !self.rest().starts_with('\'') {
                return Ok(());
            }
            self.advance(1)?;
        }
    }
    fn node(&mut self) -> Result<()> {
        self.charge(1)?;
        self.values = self.values.checked_add(1).ok_or(Error::Limit)?;
        if self.values > VALUES {
            Err(Error::Limit)
        } else {
            Ok(())
        }
    }
    fn entity(&mut self, depth: usize) -> Result<&'a str> {
        self.node()?;
        let name = self.name()?;
        self.args(depth)?;
        Ok(name)
    }
    fn args(&mut self, depth: usize) -> Result<()> {
        if depth >= 24 {
            return Err(Error::Limit);
        }
        self.depth = self.depth.max(depth + 1);
        self.require(b'(')?;
        let mut count = 0usize;
        if !self.take(b')')? {
            loop {
                self.value(depth + 1)?;
                count = count.checked_add(1).ok_or(Error::Limit)?;
                if self.take(b')')? {
                    break;
                }
                self.require(b',')?;
            }
        }
        let allocation = capacity(count)?
            .checked_mul(size_of::<Value<'_>>())
            .ok_or(Error::Limit)?;
        self.value_bytes = self
            .value_bytes
            .checked_add(allocation)
            .ok_or(Error::Limit)?;
        self.largest_vector_bytes = self.largest_vector_bytes.max(allocation);
        Ok(())
    }
    fn value(&mut self, depth: usize) -> Result<()> {
        self.node()?;
        match self.peek()? {
            Some(b'\'') => {
                self.advance(1)?;
                self.quoted()
            }
            Some(b'.') => {
                self.advance(1)?;
                self.token(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_')?;
                if !self.rest().starts_with('.') {
                    return Err(Error::Invalid);
                }
                self.advance(1)
            }
            Some(b'#') => {
                self.references = self.references.checked_add(1).ok_or(Error::Limit)?;
                self.advance(1)?;
                self.id()
            }
            Some(b'+' | b'-' | b'0'..=b'9') => {
                let number = self
                    .token(|b| matches!(b, b'+' | b'-' | b'.' | b'E' | b'0'..=b'9'))?
                    .parse::<f64>()
                    .map_err(|_| Error::Invalid)?;
                if number.is_finite() {
                    Ok(())
                } else {
                    Err(Error::Invalid)
                }
            }
            Some(b'(') => self.args(depth),
            Some(b'$' | b'*') => self.advance(1),
            Some(b'A'..=b'Z' | b'_') => self.entity(depth).map(|_| ()),
            Some(b'"') => Err(Error::Unsupported),
            _ => Err(Error::Invalid),
        }
    }
}
