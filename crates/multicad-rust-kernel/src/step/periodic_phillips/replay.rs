//! Qualified repaired packet -> exact selected rebased token graph. Final
//! represented count/order/ownership is separately checked on the full file.
use super::*;
use crate::{
    step::{self, retained_replay},
    step_read::{literal_member_binding as binding, owned_insert as read},
};
const FIXED: usize = 128 * 1024;
const BYTES: usize = 64 * 1024 * 1024;
#[derive(Clone, Copy)]
pub struct AppendOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    /// Explicit new Phi encoding stage; old NPT retains its original16M cap.
    pub maximum_encoding_work: usize,
}
impl Default for AppendOptions {
    fn default() -> Self {
        Self {
            maximum_work: 96_000_000,
            maximum_bytes: BYTES,
            maximum_encoding_work: 32_000_000,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct AppendReceipt {
    pub matched_entities: usize,
    pub document_work: usize,
    pub selection_work: usize,
    pub encoding_work: usize,
    pub binding_work: usize,
    pub work_used: usize,
    pub accounted_peak_bytes: usize,
    /// Completed independent candidate qualification, separate from new append.
    pub candidate: Verification,
    pub proposal_work: usize,
}
/// Append-only resources of one exact singleton replay.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ReplayCore {
    pub(crate) matched_entities: usize,
    pub(crate) document_work: usize,
    pub(crate) selection_work: usize,
    pub(crate) encoding_work: usize,
    pub(crate) binding_work: usize,
    pub(crate) work_used: usize,
    pub(crate) accounted_peak_bytes: usize,
}
struct Ledger<'a> {
    used: usize,
    maximum: usize,
    peak: usize,
    bytes: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Ledger<'_> {
    fn step(&mut self, n: usize) -> Result<(), Error> {
        checkpoint(self.interrupt)?;
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        if self.used > self.maximum {
            return Err(Error::Limit);
        }
        Ok(())
    }
    fn observe(&mut self, n: usize) -> Result<(), Error> {
        if n > self.bytes {
            return Err(Error::Limit);
        }
        self.peak = self.peak.max(n);
        Ok(())
    }
    fn remaining(&self) -> Result<usize, Error> {
        self.maximum
            .checked_sub(self.used)
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)
    }
    fn room(&self, n: usize) -> Result<usize, Error> {
        self.bytes
            .checked_sub(n)
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)
    }
}
impl retained_replay::Budget for Ledger<'_> {
    type Error = Error;
    fn step(&mut self, n: usize) -> Result<(), Error> {
        self.step(n)
    }
    fn observe(&mut self, n: usize) -> Result<(), Error> {
        self.observe(n)
    }
    fn invalid() -> Error {
        Error::Invalid
    }
    fn limit() -> Error {
        Error::Limit
    }
}
impl From<step::StepError> for Error {
    fn from(e: step::StepError) -> Self {
        match e {
            step::StepError::Limit => Self::Limit,
            step::StepError::Unsupported => Self::Unsupported,
            _ => Self::Invalid,
        }
    }
}
fn document_error(e: read::Error) -> Error {
    match e {
        read::Error::Invalid => Error::Invalid,
        read::Error::Limit => Error::Limit,
        read::Error::Cancelled => Error::Cancelled,
        read::Error::Unsupported => Error::Unsupported,
    }
}
fn binding_error(e: binding::Error) -> Error {
    match e {
        binding::Error::Invalid => Error::Invalid,
        binding::Error::Different => Error::Unsupported,
        binding::Error::Cancelled => Error::Cancelled,
        _ => Error::Limit,
    }
}
fn parse<'a>(
    input: &'a str,
    resident: usize,
    work: &mut Ledger<'_>,
) -> Result<read::CompoundDocument<'a>, Error> {
    let document = read::CompoundDocument::parse(
        input,
        read::DocumentOptions {
            maximum_work: work.remaining()?.min(512_000_000),
            retained_budget_bytes: work.room(resident)?.min(512 * 1024 * 1024),
        },
        work.interrupt,
    )
    .map_err(document_error)?;
    work.step(document.work_used())?;
    work.observe(resident + document.retained_budget_bytes())?;
    Ok(document)
}
fn select<'i>(
    document: &read::CompoundDocument<'i>,
    root: u32,
    resident: usize,
    work: &mut Ledger<'_>,
) -> Result<binding::SelectedRecords<'i>, Error> {
    let records = binding::selected_records(
        document,
        root,
        binding::Options {
            maximum_work: work.remaining()?.min(16_000_000),
            maximum_bytes: work.room(resident)?.min(BYTES),
            maximum_entities: 1_048_576,
        },
        work.interrupt,
    )
    .map_err(binding_error)?;
    work.step(records.work_used())?;
    work.observe(resident + records.query_peak_bytes())?;
    Ok(records)
}
impl step::Writer {
    /// Atomically append a sealed independently qualified repaired candidate.
    /// Original decimal tokens and graph sharing survive exact entity rebasing.
    /// No material authority comes from this replay core or its wrapper header.
    pub fn append_qualified_periodic_phillips(
        &mut self,
        candidate: &QualifiedStep,
        options: AppendOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<AppendReceipt, Error> {
        let candidate_bytes = candidate.text.capacity()
            + candidate.qualification.material.retained_budget_bytes()
            + std::mem::size_of::<QualifiedStep>();
        let core = self.append_qualified_singleton(
            candidate.as_str(),
            candidate_bytes,
            options,
            None,
            interrupt,
        )?;
        Ok(AppendReceipt {
            matched_entities: core.matched_entities,
            document_work: core.document_work,
            selection_work: core.selection_work,
            encoding_work: core.encoding_work,
            binding_work: core.binding_work,
            work_used: core.work_used,
            accounted_peak_bytes: core.accounted_peak_bytes,
            candidate: candidate.qualification.verification(),
            proposal_work: candidate.proposal_work_used(),
        })
    }
    /// Shared exact rebased replay of one sealed, already independently
    /// qualified singleton text. `candidate_bytes` is the caller's complete
    /// resident candidate (text, owned material and wrapper). This core adds
    /// no material authority: callers must have qualified `text` itself.
    /// With `placement` (origin, u, v, u×v), the rebased local solid is held
    /// by one declared map and placed by one `MAPPED_ITEM`; the caller must
    /// have proved that frame proper and orthonormal.
    pub(crate) fn append_qualified_singleton(
        &mut self,
        text: &str,
        candidate_bytes: usize,
        options: AppendOptions,
        placement: Option<[[f64; 3]; 4]>,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ReplayCore, Error> {
        checkpoint(interrupt)?;
        if options.maximum_work == 0
            || options.maximum_work > 128_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > BYTES
            || options.maximum_encoding_work == 0
            || options.maximum_encoding_work > 32_000_000
        {
            return Err(Error::Invalid);
        }
        let mut work = Ledger {
            used: 0,
            maximum: options.maximum_work,
            peak: 0,
            bytes: options.maximum_bytes,
            interrupt,
        };
        let base = candidate_bytes.checked_add(FIXED).ok_or(Error::Limit)?;
        work.observe(base)?;
        let source_document = parse(text, base, &mut work)?;
        // Callers pass sealed qualified singletons; refuse anything else.
        if source_document.represented_members().len() != 1 {
            return Err(Error::Invalid);
        }
        let root = source_document.represented_members()[0].item_entity();
        let source_document_work = source_document.work_used();
        let source = select(
            &source_document,
            root,
            base + source_document.retained_budget_bytes(),
            &mut work,
        )?;
        let source_selection_work = source.work_used();
        // Original component slices borrow candidate.text, not its neutral AST.
        // Drop that arena before parsing the actual emitted payload below.
        drop(source_document);
        let mut receipt = None;
        self.append_typed(|file| {
            let resident = base + source.retained_budget_bytes() + file.data.capacity();
            work.observe(resident)?;
            let before = work.used;
            let start = file.data.len();
            let trace = retained_replay::encode(
                &source,
                root,
                file,
                // Include old String capacity while exact reserve may allocate
                // a new buffer. No allocator in-place realloc is assumed.
                base + source.retained_budget_bytes() + file.data.capacity(),
                &mut work,
                options.maximum_encoding_work,
            )?;
            let encoding_work = work.used - before;
            let trace_bytes = trace.nodes.capacity() * std::mem::size_of::<retained_replay::Node>();
            // Parse the ACTUAL emitted payload bytes, not a reconstructed recipe
            // or control net. A neutral wrapper merely supplies one represented
            // root to the existing parser. IDs and all payload tokens are intact.
            let mut actual = step::StepFile::new();
            let length = file.data.len() - start;
            work.step(length + 1)?;
            let reserve = length + 16 * 1024;
            work.observe(
                base + source.retained_budget_bytes()
                    + file.data.capacity()
                    + trace_bytes
                    + reserve
                    + FIXED,
            )?;
            actual
                .data
                .try_reserve_exact(reserve)
                .map_err(|_| Error::Limit)?;
            actual.data.push_str(&file.data[start..]);
            actual.next = file.next;
            work.step(16_384)?; // fixed singleton representation/unit metadata
            actual.finish(vec![step::Entity(trace.member as usize)]);
            let actual_resident = base
                + source.retained_budget_bytes()
                + file.data.capacity()
                + trace_bytes
                + actual.data.capacity();
            let actual_document = parse(&actual.data, actual_resident, &mut work)?;
            assert_eq!(
                actual_document.represented_members().len(),
                1,
                "one actual replay payload wrapper"
            );
            let actual_root = actual_document.represented_members()[0].item_entity();
            let actual_document_work = actual_document.work_used();
            let actual_records = select(
                &actual_document,
                actual_root,
                actual_resident + actual_document.retained_budget_bytes(),
                &mut work,
            )?;
            let actual_selection_work = actual_records.work_used();
            drop(actual_document);
            let bind_resident = actual_resident + actual_records.retained_budget_bytes();
            let bound = binding::bind_retained(
                &source,
                root,
                &actual_records,
                actual_root,
                binding::Options {
                    maximum_work: work.remaining()?.min(512_000_000),
                    maximum_bytes: work.room(bind_resident)?.min(BYTES),
                    maximum_entities: 1_048_576,
                },
                interrupt,
            )
            .map_err(binding_error)?;
            work.step(bound.work_used())?;
            work.observe(bind_resident + bound.query_peak_bytes())?;
            receipt = Some(ReplayCore {
                matched_entities: bound.matched_entities(),
                document_work: source_document_work + actual_document_work,
                selection_work: source_selection_work + actual_selection_work,
                encoding_work,
                binding_work: bound.work_used(),
                work_used: work.used,
                accounted_peak_bytes: work.peak,
            });
            let member = step::Entity(trace.member as usize);
            Ok::<_, Error>(match placement {
                Some(frame) => {
                    work.step(4096)?;
                    file.mapped_singleton(member, frame)
                }
                None => member,
            })
        })?;
        Ok(receipt.expect("successful typed atomic replay"))
    }
}
#[cfg(test)]
mod tests;
