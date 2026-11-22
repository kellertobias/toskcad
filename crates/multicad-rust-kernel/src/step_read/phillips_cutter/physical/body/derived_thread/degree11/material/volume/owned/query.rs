//! Whole-box predicates of the SAME owned D11 material. The thread predicate
//! consumes actual tensor coefficients and its earned whole(theta,z) radial
//! graph; neither a thread recipe nor equal-v rails determine membership.
use super::*;
use crate::{
    head_profile::{RepairedHeadError, RepairedHeadLocation, RepairedHeadQueryOptions},
    phillips_drive::Classification as C,
};
mod thread;
const HEAD_WORK: usize = 16_000_000;
const THREAD_WORK: usize = 64_000_000;
const MAX_QUERY_WORK: usize = HEAD_WORK + THREAD_WORK + 16384;
const MAX_THREAD_CELLS: usize = 65_536;
const MAX_THREAD_DEPTH: usize = 64;
// Head query:65 fixed(range,depth) entries (<2KiB), an inline <=5-part
// source/view, and fixed two-coordinate normalized-ray/cubic evaluations.
// ArcSource::evaluate has no heap, recursion or adaptive jets.256KiB bounds
// those simultaneous scalar/array temporaries and call frames conservatively.
const HEAD_SCRATCH: usize = 256 * 1024;
const METADATA_SCRATCH: usize = 64 * 1024;
type R<T> = std::result::Result<T, Degree11MaterialQueryError>;
#[derive(Clone, Copy, Debug)]
pub struct Degree11MaterialQueryOptions {
    /// Fresh bounded query credits, not prior material construction work or
    /// measured CPU instructions. Head/thread retain separate16M/64M ceilings.
    pub maximum_work: usize,
    /// Complete reachable owned source plus conservative simultaneous scratch.
    pub maximum_bytes: usize,
    pub maximum_head_cells: usize,
    pub maximum_thread_cells: usize,
    pub maximum_thread_depth: usize,
}
impl Default for Degree11MaterialQueryOptions {
    fn default() -> Self {
        Self {
            maximum_work: MAX_QUERY_WORK,
            maximum_bytes: MAX_BYTES,
            maximum_head_cells: 4096,
            maximum_thread_cells: 4096,
            maximum_thread_depth: 48,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Degree11MaterialQueryError {
    Invalid,
    Limit,
    Cancelled,
    Head(RepairedHeadError),
    Thread(crate::tapered_thread::Error),
    Cavity(crate::phillips_drive::Error),
}
pub struct Degree11MaterialQuery<'s> {
    source: &'s OwnedDegree11Material,
    bounds: [[f64; 2]; 3],
    location: C,
    used: usize,
    remaining: usize,
    cells: usize,
    peak: usize,
}
impl Degree11MaterialQuery<'_> {
    pub fn source(&self) -> &OwnedDegree11Material {
        self.source
    }
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        self.bounds
    }
    pub fn location(&self) -> C {
        self.location
    }
    /// Debited conservative primitive-work credits; not measured executions.
    pub fn work_used(&self) -> usize {
        self.used
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
    pub fn visited_thread_cells(&self) -> usize {
        self.cells
    }
    /// Conservative source+scratch reservation, not allocator RSS.
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
struct QueryWork<'a> {
    remaining: usize,
    thread_remaining: usize,
    cells: usize,
    maximum_cells: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl QueryWork<'_> {
    fn charge(&mut self, n: usize) -> R<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Degree11MaterialQueryError::Cancelled);
        }
        self.remaining = self
            .remaining
            .checked_sub(n)
            .ok_or(Degree11MaterialQueryError::Limit)?;
        Ok(())
    }
    fn thread_charge(&mut self, n: usize) -> R<()> {
        self.charge(n)?;
        self.thread_remaining = self
            .thread_remaining
            .checked_sub(n)
            .ok_or(Degree11MaterialQueryError::Limit)?;
        Ok(())
    }
    fn cell(&mut self) -> R<()> {
        self.thread_charge(thread::NODE_WORK)?;
        if self.cells == self.maximum_cells {
            return Err(Degree11MaterialQueryError::Limit);
        }
        self.cells += 1;
        Ok(())
    }
}
fn valid(options: Degree11MaterialQueryOptions, bounds: [[f64; 2]; 3]) -> bool {
    options.maximum_work <= MAX_QUERY_WORK
        && options.maximum_bytes <= MAX_BYTES
        && options.maximum_head_cells > 0
        && options.maximum_head_cells <= 1_000_000
        && options.maximum_thread_cells > 0
        && options.maximum_thread_cells <= MAX_THREAD_CELLS
        && options.maximum_thread_depth <= MAX_THREAD_DEPTH
        && bounds.iter().all(|b| {
            b[0].is_finite() && b[1].is_finite() && b[0] <= b[1] && b.iter().all(|x| x.abs() <= 1e7)
        })
}
impl OwnedDegree11Material {
    /// Classifies EVERY point of this finite material-local coordinate box.
    /// Actual boundary/join crossings may remain unresolved. Typed refusal and
    /// cancellation publish no result and never mutate the immutable source.
    pub fn classify_local_box(
        &self,
        bounds: [[f64; 2]; 3],
        options: Degree11MaterialQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> R<Degree11MaterialQuery<'_>> {
        if !valid(options, bounds) {
            return Err(Degree11MaterialQueryError::Invalid);
        }
        let mut work = QueryWork {
            remaining: options.maximum_work,
            thread_remaining: THREAD_WORK,
            cells: 0,
            maximum_cells: options.maximum_thread_cells,
            interrupt,
        };
        // Covers the bounded source descriptor/capacity walk before allocation.
        work.charge(4096)?;
        let peak = self
            .retained_budget_bytes()
            .checked_add(thread::scratch_bytes(options.maximum_thread_depth))
            .and_then(|x| x.checked_add(HEAD_SCRATCH + METADATA_SCRATCH))
            .ok_or(Degree11MaterialQueryError::Limit)?;
        if peak > options.maximum_bytes {
            return Err(Degree11MaterialQueryError::Limit);
        }
        let location = if bounds
            .iter()
            .zip(self.local_enclosure_bounds_mm())
            .any(|(b, e)| b[1] < e[0] || b[0] > e[1])
        {
            C::Outside
        } else {
            self.locate_box(bounds.map(|b| I::new(b[0], b[1])), options, &mut work)?
        };
        work.charge(1)?;
        Ok(Degree11MaterialQuery {
            source: self,
            bounds,
            location,
            used: options.maximum_work - work.remaining,
            remaining: work.remaining,
            cells: work.cells,
            peak,
        })
    }
    fn locate_box(
        &self,
        p: [I; 3],
        options: Degree11MaterialQueryOptions,
        work: &mut QueryWork<'_>,
    ) -> R<C> {
        // Fixed classify::enclosure: one envelope and <=4 Corner::classify
        // calls. Per corner <=5 interval multiplies,3 squares,16 add/sub,
        // fixed rotations/selections; envelope <=4 multiplies,1 division,
        // 6 squares/1 sqrt and12 add/sub. At <=24 scalar arithmetic/selection
        // credits per multiply, <=12 per square and <=6 per add/sub, plus
        // bounded predicates/copies,4096 credits covers the complete call.
        // No trig evaluation, allocation, recursion or adaptive loop occurs.
        work.charge(4096)?;
        // SAME actual Cutter translation as the earned head/material source.
        let cavity_box = [p[0], p[1], p[2] - I::point(self.chart().outer_top_z)];
        let cavity = self
            .cavity()
            .classify_box_mm(cavity_box.map(|b| [b.lo, b.hi]))
            .map_err(Degree11MaterialQueryError::Cavity)?;
        if cavity == C::Inside {
            return Ok(C::Outside);
        }
        let head = self.head().view();
        let h = head
            .classify_box(
                p.map(|b| [b.lo, b.hi]),
                RepairedHeadQueryOptions {
                    maximum_work: work.remaining.min(HEAD_WORK),
                    maximum_cells: options.maximum_head_cells,
                },
                work.interrupt,
            )
            .map_err(Degree11MaterialQueryError::Head)?;
        work.charge(h.work_used())?;
        let union = if h.location() == RepairedHeadLocation::Inside {
            C::Inside
        } else {
            let t = thread::classify(self, p, options.maximum_thread_depth, work)?;
            match (h.location(), t) {
                (_, C::Inside) => C::Inside,
                (RepairedHeadLocation::Outside, C::Outside) => C::Outside,
                _ => C::UnresolvedBoundary,
            }
        };
        Ok(match (union, cavity) {
            (C::Outside, _) => C::Outside,
            (C::Inside, C::Outside) => C::Inside,
            _ => C::UnresolvedBoundary,
        })
    }
}
#[cfg(test)]
pub(crate) mod tests;
