//! Whole-box predicates of the same owned head/thread minus actual cavity.
use super::*;
use crate::{
    head_profile::{RepairedHeadError, RepairedHeadLocation, RepairedHeadQueryOptions},
    phillips_drive::Classification as C,
    radial_tensor::{BoxQueryOptions, Location},
};
const MAXIMUM_WORK: usize = 18_008_192;
type R<T> = std::result::Result<T, OwnedMaterialQueryError>;
#[derive(Clone, Copy)]
pub struct OwnedMaterialQueryOptions {
    /// One aggregate balance. Head and thread keep their original 16M/2M caps.
    pub maximum_work: usize,
    pub maximum_head_cells: usize,
}
impl Default for OwnedMaterialQueryOptions {
    fn default() -> Self {
        Self {
            maximum_work: MAXIMUM_WORK,
            maximum_head_cells: 4096,
        }
    }
}
/// Provider errors retain their origin; no uncertain result becomes material.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OwnedMaterialQueryError {
    Invalid,
    Limit,
    Cancelled,
    Head(RepairedHeadError),
    Thread(crate::tapered_thread::Error),
    Cavity(crate::phillips_drive::Error),
}
pub struct OwnedMaterialQuery<'s> {
    source: &'s OwnedDerivedMaterial,
    bounds: [[f64; 2]; 3],
    location: C,
    used: usize,
    remaining: usize,
}
impl OwnedMaterialQuery<'_> {
    pub fn source(&self) -> &OwnedDerivedMaterial {
        self.source
    }
    /// Material-local chart, before the common-mouth construction map.
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        self.bounds
    }
    pub fn location(&self) -> C {
        self.location
    }
    pub fn work_used(&self) -> usize {
        self.used
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
}
struct QueryWork<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl QueryWork<'_> {
    fn charge(&mut self, n: usize) -> R<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(OwnedMaterialQueryError::Cancelled);
        }
        self.remaining = self
            .remaining
            .checked_sub(n)
            .ok_or(OwnedMaterialQueryError::Limit)?;
        Ok(())
    }
}
impl OwnedDerivedMaterial {
    /// Classifies every point in this finite MATERIAL-LOCAL box. This uses
    /// actual retained sources, not raw STEP mass, nominal recipes or meshes.
    /// Virtual join disks and surface crossings may remain unresolved.
    pub fn classify_local_box(
        &self,
        bounds: [[f64; 2]; 3],
        options: OwnedMaterialQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> R<OwnedMaterialQuery<'_>> {
        let mut w = QueryWork {
            remaining: options.maximum_work,
            interrupt,
        };
        w.charge(64)?;
        let box_intervals = bounds.map(|b| I::new(b[0], b[1]));
        if options.maximum_work > MAXIMUM_WORK
            || options.maximum_head_cells == 0
            || options.maximum_head_cells > 1_000_000
            || box_intervals.into_iter().any(|b| !b.valid())
        {
            return Err(OwnedMaterialQueryError::Invalid);
        }
        let location = if bounds
            .iter()
            .zip(self.local_enclosure)
            .any(|(b, e)| b[1] < e[0] || b[0] > e[1])
        {
            C::Outside
        } else {
            self.locate_box(box_intervals, options.maximum_head_cells, &mut w)?
        };
        w.charge(1)?;
        Ok(OwnedMaterialQuery {
            source: self,
            bounds,
            location,
            used: options.maximum_work - w.remaining,
            remaining: w.remaining,
        })
    }
    fn locate_box(&self, p: [I; 3], cells: usize, w: &mut QueryWork<'_>) -> R<C> {
        // Cavity geometry starts at its mouth at z=0. The OWNED material
        // carries the same cavity translated by the exact common station.
        w.charge(4096)?;
        let cavity_box = [p[0], p[1], p[2] - I::point(self.chart.outer_top_z)];
        let cavity = self
            .cavity
            .classify_box_mm(cavity_box.map(|b| [b.lo, b.hi]))
            .map_err(OwnedMaterialQueryError::Cavity)?;
        if cavity == C::Inside {
            return Ok(C::Outside);
        }
        let bounds = p.map(|b| [b.lo, b.hi]);
        let head = self.head.view();
        let query = head
            .classify_box(
                bounds,
                RepairedHeadQueryOptions {
                    maximum_work: w.remaining.min(16_000_000),
                    maximum_cells: cells,
                },
                w.interrupt,
            )
            .map_err(OwnedMaterialQueryError::Head)?;
        w.charge(query.work_used())?;
        let union = if query.location() == RepairedHeadLocation::Inside {
            C::Inside
        } else {
            let thread = self
                .thread
                .classify_box(
                    bounds,
                    BoxQueryOptions {
                        maximum_work: w.remaining.min(2_000_000),
                    },
                    w.interrupt,
                )
                .map_err(OwnedMaterialQueryError::Thread)?;
            w.charge(thread.work_used())?;
            match (query.location(), thread.location()) {
                (_, Location::Inside) => C::Inside,
                (RepairedHeadLocation::Outside, Location::Outside) => C::Outside,
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
mod tests;
