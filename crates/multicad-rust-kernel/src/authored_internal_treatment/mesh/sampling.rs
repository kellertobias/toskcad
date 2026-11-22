//! Computational partition proposals only. No proposal certifies source/facets.
use super::*;
pub(super) const MAX_ROWS: u32 = 65_536;
pub(super) const MAX_SUBDIVISIONS: usize = 1024;
const ALTERNATE_OFFSET: u32 = 3;
// Half the finest permitted angular cell: no legal dyadic row places the
// original cut on this shifted sampling seam. Source phase stays unchanged.
pub(super) const TURN_OFFSET: f64 = 1. / (2. * MAX_ROWS as f64);
pub(super) enum Attempt<T> {
    Complete(T),
    Refine,
    Repartition,
}
pub(super) struct Grid {
    pub subdivisions: usize,
    base_rows: u32,
    alternative: bool,
    pub shifted: bool,
}
impl Grid {
    pub fn new() -> Self {
        Self {
            subdivisions: 1,
            base_rows: 64,
            alternative: false,
            shifted: false,
        }
    }
    pub fn rows(&self) -> u32 {
        self.base_rows
            + if self.alternative {
                ALTERNATE_OFFSET
            } else {
                0
            }
    }
    pub fn advance(&mut self, repartition: bool) -> Result<(), Error> {
        if repartition && !self.alternative && self.base_rows <= MAX_ROWS - ALTERNATE_OFFSET {
            self.alternative = true;
        } else if repartition && !self.shifted && self.base_rows < MAX_ROWS {
            self.shifted = true;
            self.alternative = false;
        } else {
            if self.subdivisions >= MAX_SUBDIVISIONS || self.base_rows >= MAX_ROWS {
                return Err(Error::Limit);
            }
            // The preceding chart bounds and private state transitions prove
            // these products fit. Failed alternate returns to original dyadic
            // refinement; no partial candidate or ownership cache survives.
            self.subdivisions *= 2;
            self.base_rows *= 2;
            self.alternative = false;
        }
        Ok(())
    }
    pub fn refine_angular(&mut self) -> Result<(), Error> {
        if self.base_rows >= MAX_ROWS {
            return Err(Error::Limit);
        }
        self.base_rows *= 2;
        self.alternative = false;
        Ok(())
    }
}
/// A nominal short literal edge relative to its endpoint evaluation radii
/// proposes another partition. This does NOT classify source separation, erase
/// a positive domain, or certify any facet; every acceptance proof runs again.
pub(super) fn short_literal_edge(points: [[f64; 3]; 3], errors: [f64; 3]) -> bool {
    [(0, 1), (1, 2), (2, 0)].into_iter().any(|(a, b)| {
        let edge = (0..3).fold(0., |sum, k| sum + (points[a][k] - points[b][k]).abs());
        edge <= errors[a] + errors[b]
    })
}

/// A literal edge shorter than a fixed multiple of the binary64 spacing at
/// its coordinates: its incident facets cannot have resolvable orientation
/// or separation against facets of ordinary size. A partition proposal only.
pub(super) fn conditioning_edge(points: [[f64; 3]; 2]) -> bool {
    let scale = points.iter().flatten().fold(1f64, |m, x| m.max(x.abs()));
    let edge = (0..3).fold(0., |sum, k| sum + (points[0][k] - points[1][k]).abs());
    edge <= 1024. * f64::EPSILON * scale
}

#[cfg(test)]
mod tests;
