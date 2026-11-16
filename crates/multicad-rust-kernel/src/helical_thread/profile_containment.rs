//! Whole-period radial containment of two ACTUAL retained thread profiles.
//! Placement and finite material domains remain the consumer's responsibility.
use super::*;
use crate::tapered_thread::interval::I;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProfileContainmentError {
    Invalid,
    Unresolved,
    Limit,
    Cancelled,
}
#[derive(Clone, Copy, Debug)]
pub struct ProfileContainmentOptions {
    /// Maximum complete pair-profile enclosure queries. The pending-cell array
    /// has 65 entries, with at most 64 subdivision levels and no allocation.
    pub maximum_cells: usize,
}
impl Default for ProfileContainmentOptions {
    fn default() -> Self {
        Self {
            maximum_cells: 65_536,
        }
    }
}
#[derive(Debug)]
pub struct ProfileContainment<'a> {
    external: &'a Thread,
    cavity: &'a Thread,
    lower: f64,
    cells: usize,
}
impl ProfileContainment<'_> {
    pub fn external(&self) -> &Thread {
        self.external
    }
    pub fn cavity(&self) -> &Thread {
        self.cavity
    }
    /// Zero is earned only for identical retained profile coefficients. It
    /// permits coincident boundary contact, not positive solid intersection.
    pub fn minimum_gap_mm(&self) -> f64 {
        self.lower
    }
    pub fn cells_visited(&self) -> usize {
        self.cells
    }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), ProfileContainmentError> {
    if interrupt.is_some_and(|f| f()) {
        Err(ProfileContainmentError::Cancelled)
    } else {
        Ok(())
    }
}
#[derive(Clone, Copy)]
struct Cell {
    lo: f64,
    hi: f64,
    depth: u8,
}
impl Thread {
    /// Prove external radius <= cavity radius at EVERY common phase. No
    /// sampling, nominal-fit lookup or manufacturing tolerance closes a gap.
    /// Both profiles must use the exact same pitch/start/hand phase law.
    pub fn prove_profile_containment<'a>(
        &'a self,
        cavity: &'a Thread,
        options: ProfileContainmentOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ProfileContainment<'a>, ProfileContainmentError> {
        use ProfileContainmentError::*;
        checkpoint(interrupt)?;
        if options.maximum_cells == 0
            || options.maximum_cells > 1_000_000
            || self.spec.pitch_mm != cavity.spec.pitch_mm
            || self.spec.starts != cavity.spec.starts
            || self.spec.left_hand != cavity.spec.left_hand
        {
            return Err(Invalid);
        }
        let mut cavity_profile = cavity.spec;
        cavity_profile.length_mm = self.spec.length_mm;
        if self.spec == cavity_profile
            && self.breaks == cavity.breaks
            && self.flanks == cavity.flanks
            && self.root == cavity.root
        {
            checkpoint(interrupt)?;
            return Ok(ProfileContainment {
                external: self,
                cavity,
                lower: 0.,
                cells: 0,
            });
        }
        let initial = Cell {
            lo: 0.,
            hi: self.spec.pitch_mm,
            depth: 0,
        };
        let mut stack = [initial; 65];
        let mut len = 1;
        let mut cells = 0;
        let mut lower = f64::INFINITY;
        while len != 0 {
            checkpoint(interrupt)?;
            if cells == options.maximum_cells {
                return Err(Limit);
            }
            cells += 1;
            len -= 1;
            let cell = stack[len];
            let phase = I::new(cell.lo, cell.hi);
            let outer = self.radius_at_phase_interval(phase).ok_or(Unresolved)?;
            let inner = cavity.radius_at_phase_interval(phase).ok_or(Unresolved)?;
            let gap = inner - outer;
            if gap.lo > 0. {
                lower = lower.min(gap.lo);
                continue;
            }
            if gap.hi < 0. {
                return Err(Unresolved);
            }
            if cell.depth == 64 {
                return Err(Limit);
            }
            let mid = cell.lo + (cell.hi - cell.lo) * 0.5;
            if mid <= cell.lo || mid >= cell.hi {
                return Err(Unresolved);
            }
            // Visit lower phases first; the deferred siblings fit the fixed
            // array because a depth-64 cell is never subdivided further.
            stack[len] = Cell {
                lo: mid,
                hi: cell.hi,
                depth: cell.depth + 1,
            };
            stack[len + 1] = Cell {
                lo: cell.lo,
                hi: mid,
                depth: cell.depth + 1,
            };
            len += 2;
        }
        checkpoint(interrupt)?;
        Ok(ProfileContainment {
            external: self,
            cavity,
            lower,
            cells,
        })
    }
}

#[cfg(test)]
mod tests;
