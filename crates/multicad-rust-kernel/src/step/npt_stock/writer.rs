//! Closed source-known STEP candidate. Independent decoded material remains
//! a separate gate; no Native export/import authority is granted here.
use super::*;
use crate::step::{Coordinates, Entities, Entity, StepFile, Writer};
use std::fmt::{self, Write as _};
mod emit;
#[derive(Clone, Copy, Debug)]
pub struct EncodingReceipt {
    work: usize,
    bytes: usize,
}
impl EncodingReceipt {
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.bytes
    }
}
#[derive(Debug)]
pub struct CandidateFile {
    text: String,
    receipt: EncodingReceipt,
}
impl CandidateFile {
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn into_string(self) -> String {
        self.text
    }
    pub fn receipt(&self) -> EncodingReceipt {
        self.receipt
    }
}
pub fn write_source_candidate(
    plan: &SurfacePlan<'_>,
    maximum_work: usize,
    maximum_bytes: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<CandidateFile, Error> {
    let mut writer = Writer::new();
    let receipt =
        writer.append_source_npt_candidate(plan, maximum_work, maximum_bytes, interrupt)?;
    let text = writer.finish()?;
    Ok(CandidateFile { text, receipt })
}
impl Writer {
    /// Append one sealed source-known NPT packet using the common atomic
    /// record/ID transaction. This does not admit decoded material.
    pub fn append_source_npt_candidate(
        &mut self,
        plan: &SurfacePlan<'_>,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<EncodingReceipt, Error> {
        let mut receipt = None;
        self.append_typed(|file| {
            let mut encoder = Encoder::new(file, plan, maximum_work, maximum_bytes, interrupt)?;
            let solid = emit::body(&mut encoder, plan)?;
            encoder.reserve(4096)?;
            encoder.work.charge(4096)?;
            receipt = Some(EncodingReceipt {
                work: encoder.work.used,
                bytes: encoder.peak,
            });
            Ok::<_, Error>(solid)
        })?;
        Ok(receipt.expect("successful sealed NPT append"))
    }
}
// Only a sealed verified candidate may replay with a distinct actual-encoding
// phase. Its preparation was already charged by verified construction. The
// legacy public source-candidate API above retains its inherited-work contract.
pub(super) fn append_verified_replay(
    writer: &mut Writer,
    candidate: &super::verified::VerifiedCandidate<'_, '_>,
    maximum_work: usize,
    maximum_bytes: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<EncodingReceipt, Error> {
    let plan = candidate.source_plan();
    let mut receipt = None;
    writer.append_typed(|file| {
        let mut work = Work::new(maximum_work, interrupt)?;
        work.charge(0)?;
        let mut encoder = Encoder::with_prepared_work(file, plan, work, maximum_bytes)?;
        let solid = emit::body(&mut encoder, plan)?;
        encoder.reserve(4096)?;
        encoder.work.charge(4096)?;
        receipt = Some(EncodingReceipt {
            work: encoder.work.used,
            bytes: encoder.peak,
        });
        Ok::<_, Error>(solid)
    })?;
    Ok(receipt.expect("successful sealed verified replay"))
}
struct Counter(usize);
impl fmt::Write for Counter {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        self.0 = self.0.checked_add(s.len()).ok_or(fmt::Error)?;
        Ok(())
    }
}
struct Encoder<'file, 'cancel> {
    file: &'file mut StepFile,
    work: Work<'cancel>,
    memory: Memory,
    peak: usize,
}
impl<'file, 'cancel> Encoder<'file, 'cancel> {
    fn new(
        file: &'file mut StepFile,
        plan: &SurfacePlan<'_>,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&'cancel dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut work = Work::new(maximum_work, interrupt)?;
        work.charge(plan.work_used())?;
        Self::with_prepared_work(file, plan, work, maximum_bytes)
    }
    fn with_prepared_work(
        file: &'file mut StepFile,
        plan: &SurfacePlan<'_>,
        work: Work<'cancel>,
        maximum_bytes: usize,
    ) -> Result<Self, Error> {
        if maximum_bytes == 0 || maximum_bytes > 64 * 1024 * 1024 {
            return Err(Error::Invalid);
        }
        let used = plan.retained_budget_bytes() + file.data.capacity() + 4096;
        let memory = Memory::new(used, maximum_bytes)?;
        Ok(Self {
            file,
            work,
            memory,
            peak: used,
        })
    }
    fn allocate<T>(&mut self, count: usize) -> Result<Vec<T>, Error> {
        let out = self.memory.allocate(count, &mut self.work)?;
        self.peak = self.peak.max(self.memory.used);
        Ok(out)
    }
    fn entity(&mut self, args: fmt::Arguments<'_>) -> Result<Entity, Error> {
        let mut count = Counter(0);
        writeln!(&mut count, "{}={args};", Entity(self.file.next)).map_err(|_| Error::Limit)?;
        // One encoded-byte unit covers both bounded formatting visits. The
        // exact charge and allocation admission precede transaction mutation.
        self.work.charge(count.0 + 1)?;
        self.reserve(count.0)?;
        Ok(self.file.entity(args))
    }
    fn reserve(&mut self, additional: usize) -> Result<(), Error> {
        self.work.charge(1)?;
        let required = self
            .file
            .data
            .len()
            .checked_add(additional)
            .ok_or(Error::Limit)?;
        let old = self.file.data.capacity();
        if required <= old {
            return Ok(());
        }
        let desired = old.saturating_add((old / 2).max(64 * 1024));
        let target = desired.max(required).min(self.memory.remaining());
        if target < required {
            return Err(Error::Limit);
        }
        let before = self.memory.used;
        self.memory.add(target)?;
        self.peak = self.peak.max(self.memory.used);
        self.work.charge(
            1 + self.file.data.len().div_ceil(size_of::<usize>()) + target.div_ceil(4096),
        )?;
        self.file
            .data
            .try_reserve_exact(target - self.file.data.len())
            .map_err(|_| Error::Limit)?;
        // Pinned RawVec exact growth retains the requested capacity; allocator
        // RSS is separate. Keep and check the actual returned capacity.
        self.memory.used = before - old;
        self.memory.add(self.file.data.capacity())?;
        Ok(())
    }
    fn point(&mut self, p: [f64; 3]) -> Result<Entity, Error> {
        self.entity(format_args!("CARTESIAN_POINT('',{})", Coordinates(p)))
    }
}

#[cfg(test)]
mod tests;
