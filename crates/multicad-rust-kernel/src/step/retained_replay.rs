//! Bounded token replay for components already parsed and independently
//! qualified by their source-specific consumer. Entity renumbering preserves
//! numeric spelling, exact references and record/component incidence. This
//! core never recognizes occupied material or supplies geometry authority.
use super::*;
use crate::step_read::literal_member_binding::{RetainedRecords, Token, token};
use std::mem::size_of;
const FIXED: usize = 128 * 1024;
pub(crate) const OUTPUT: usize = 16 * 1024;
pub(crate) trait Budget {
    type Error;
    fn step(&mut self, n: usize) -> Result<(), Self::Error>;
    fn observe(&mut self, n: usize) -> Result<(), Self::Error>;
    fn invalid() -> Self::Error;
    fn limit() -> Self::Error;
}
// Every counted preflight/emit operation checks its child cap BEFORE doing
// work. The original NPT facade supplies sixteen million; new consumers may
// choose a separately named cap without resetting their shared outer ledger.
struct EncodingBudget<'a, W: Budget> {
    caller: &'a mut W,
    used: usize,
    maximum: usize,
}
impl<W: Budget> Budget for EncodingBudget<'_, W> {
    type Error = W::Error;
    fn step(&mut self, n: usize) -> Result<(), Self::Error> {
        let used = self
            .used
            .checked_add(n)
            .filter(|n| *n <= self.maximum)
            .ok_or(W::limit())?;
        self.caller.step(n)?;
        self.used = used;
        Ok(())
    }
    fn observe(&mut self, n: usize) -> Result<(), Self::Error> {
        self.caller.observe(n)
    }
    fn invalid() -> Self::Error {
        W::invalid()
    }
    fn limit() -> Self::Error {
        W::limit()
    }
}
#[derive(Clone)]
pub(crate) struct Node {
    pub(crate) id: u32,
    pub(crate) body: Range<usize>,
}
pub(crate) struct Trace {
    pub(crate) nodes: Vec<Node>,
    pub(crate) member: u32,
}
pub(crate) struct View<'a> {
    pub(crate) data: &'a str,
    pub(crate) trace: &'a Trace,
}
impl RetainedRecords for View<'_> {
    fn record_count(&self) -> usize {
        self.trace.nodes.len()
    }
    fn record_index(&self, id: u32) -> Option<usize> {
        self.trace.nodes.binary_search_by_key(&id, |n| n.id).ok()
    }
    fn component(&self, record: usize, index: usize) -> Option<&str> {
        if index != 0 {
            return None;
        }
        Some(&self.data[self.trace.nodes.get(record)?.body.clone()])
    }
}
pub(crate) fn digits(mut n: usize) -> usize {
    let mut digits = 1;
    while n >= 10 {
        n /= 10;
        digits += 1;
    }
    digits
}
fn rebased<W: Budget>(
    view: &impl RetainedRecords,
    start: usize,
    id: u32,
) -> Result<usize, W::Error> {
    start
        .checked_add(view.record_index(id).ok_or(W::invalid())?)
        .ok_or(W::limit())
}
pub(crate) fn encode<W: Budget>(
    view: &impl RetainedRecords,
    root: u32,
    file: &mut StepFile,
    base: usize,
    work: &mut W,
    maximum_encoding_work: usize,
) -> Result<Trace, W::Error> {
    let mut bounded = EncodingBudget {
        caller: work,
        used: 0,
        maximum: maximum_encoding_work,
    };
    let work = &mut bounded;
    let start = file.next;
    let end = start.checked_add(view.record_count()).ok_or(W::limit())?;
    if end - 1 > crate::step_read::MAX_STEP_ENTITIES {
        return Err(W::limit());
    }
    let mut output = 0usize;
    for index in 0..view.record_count() {
        work.step(32)?;
        let mut component = 0;
        while let Some(mut value) = view.component(index, component) {
            work.step(value.len())?;
            while let Some(token) = token(&mut value) {
                let length = match token {
                    Token::Literal(s) => s.len(),
                    Token::Reference(id) => {
                        work.step(32)?;
                        1 + digits(rebased::<W>(view, start, id)?)
                    }
                };
                output = output.checked_add(length).ok_or(W::limit())?;
            }
            component += 1;
        }
        if component == 0 {
            return Err(W::invalid());
        }
        output = output
            .checked_add(digits(start + index) + 8)
            .ok_or(W::limit())?;
    }
    let reserve = output.checked_add(OUTPUT).ok_or(W::limit())?;
    let requested = file.data.len().checked_add(reserve).ok_or(W::limit())?;
    let trace_bytes = view
        .record_count()
        .checked_mul(size_of::<Node>())
        .ok_or(W::limit())?;
    work.observe(base + requested + trace_bytes + FIXED)?;
    file.data
        .try_reserve_exact(reserve)
        .map_err(|_| W::limit())?;
    let mut nodes = Vec::new();
    nodes
        .try_reserve_exact(view.record_count())
        .map_err(|_| W::limit())?;
    work.observe(base + file.data.capacity() + nodes.capacity() * size_of::<Node>() + FIXED)?;
    for index in 0..view.record_count() {
        work.step(32)?;
        let id = file.next;
        file.next += 1;
        let prefix = file.data.len();
        write!(file.data, "#{id}=").expect("bounded String write");
        let body_start = file.data.len();
        let complex = view.component(index, 1).is_some();
        if complex {
            file.data.push('(');
        }
        let mut component = 0;
        while let Some(mut value) = view.component(index, component) {
            work.step(value.len())?;
            while let Some(value) = token(&mut value) {
                match value {
                    Token::Literal(value) => {
                        work.step(value.len() + 1)?;
                        file.data.push_str(value);
                    }
                    Token::Reference(old) => {
                        work.step(48)?;
                        write!(file.data, "#{}", rebased::<W>(view, start, old)?)
                            .expect("bounded reference write");
                    }
                }
            }
            component += 1;
        }
        if complex {
            file.data.push(')');
        }
        let body_end = file.data.len();
        file.data.push_str(";\n");
        // Token writes above already charge the body once; only the
        // record framing is additional (there is no temporary body copy).
        work.step(body_start - prefix + 2 + usize::from(complex) * 2)?;
        nodes.push(Node {
            id: u32::try_from(id).expect("bounded STEP entity range"),
            body: body_start..body_end,
        });
    }
    work.observe(base + file.data.capacity() + nodes.capacity() * size_of::<Node>() + FIXED)?;
    Ok(Trace {
        nodes,
        member: u32::try_from(rebased::<W>(view, start, root)?).expect("bounded STEP root"),
    })
}
