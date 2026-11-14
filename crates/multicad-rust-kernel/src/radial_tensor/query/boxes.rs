//! Actual input enclosures, sharing one explicit core/sector query balance.
use super::*;
const MAX_WORK: usize = 2_000_000;
#[derive(Clone, Copy)]
pub struct BoxQueryOptions {
    pub maximum_work: usize,
}
impl Default for BoxQueryOptions {
    fn default() -> Self {
        Self {
            maximum_work: MAX_WORK,
        }
    }
}
pub struct BoxQuery<'s> {
    source: &'s RadialBoundary,
    bounds: [[f64; 2]; 3],
    location: Location,
    used: usize,
    remaining: usize,
}
impl BoxQuery<'_> {
    pub fn source(&self) -> &RadialBoundary {
        self.source
    }
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        self.bounds
    }
    pub fn location(&self) -> Location {
        self.location
    }
    pub fn work_used(&self) -> usize {
        self.used
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
}
impl RadialBoundary {
    /// Classify every point in a closed WORLD coordinate box, including input
    /// transform uncertainty. Source sewing collars remain in force. This
    /// predicate proves no Boolean fragment, mesh, or manufacturing boundary.
    pub fn classify_box(
        &self,
        bounds: [[f64; 2]; 3],
        options: BoxQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<BoxQuery<'_>, Error> {
        mesh::parallel::local_budget(|| {
            let mut work = Work {
                remaining: options.maximum_work,
                interrupt,
            };
            work.step()?;
            let p = bounds.map(|[lo, hi]| I::new(lo, hi));
            if options.maximum_work > MAX_WORK || p.into_iter().any(|p| !p.valid()) {
                return Err(Error::Invalid);
            }
            // Fixed affine projection, collar/core arithmetic and final
            // publication all consume this SAME balance as tensor refinement.
            work.charge(128)?;
            let q = self.ordered.caps.axial.chart.project_interval(p);
            let mut location = self.classify_core_projected(q)?;
            if location == Location::Unresolved {
                location = self.classify_projected(q, &mut work)?;
            }
            work.step()?;
            Ok(BoxQuery {
                source: self,
                bounds,
                location,
                used: options.maximum_work - work.remaining,
                remaining: work.remaining,
            })
        })
    }
}
#[cfg(test)]
mod tests;
