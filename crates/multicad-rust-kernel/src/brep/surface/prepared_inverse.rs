//! Optional point-independent NURBS inverse proposal grid, owned by one borrow.
use super::*;
use std::cell::OnceCell;
// Optional cache construction/storage only, not the inverse's total workspace.
// Existing per-query ranking/sort/Newton storage retains its original behavior.
const MAX_BYTES: usize = 4 * 1024 * 1024;
// NurbsSurface constructors bound degree<=32. Evaluation retains at most five
// basis Vecs of <=33 doubles, their headers and small fixed arithmetic frames.
const EVALUATE_SCRATCH: usize = 8 * 1024;
#[derive(Clone, Copy)]
pub(super) struct Candidate {
    pub uv: [f64; 2],
    pub position: Point,
}
pub(crate) struct PreparedInverse<'a> {
    surface: &'a Surface,
    grid: OnceCell<Option<Vec<Candidate>>>,
    maximum_bytes: usize,
}
impl<'a> PreparedInverse<'a> {
    pub(crate) fn new(surface: &'a Surface) -> Self {
        Self {
            surface,
            grid: OnceCell::new(),
            maximum_bytes: MAX_BYTES,
        }
    }
    pub(crate) fn invert(&self, point: Point, seed: Option<[f64; 2]>) -> [f64; 2] {
        self.surface.invert_prepared(point, seed, Some(self))
    }
    pub(super) fn candidates(&self) -> Option<&[Candidate]> {
        self.grid.get_or_init(|| self.prepare()).as_deref()
    }
    fn prepare(&self) -> Option<Vec<Candidate>> {
        let Surface::Nurbs { surface, .. } = self.surface else {
            return None;
        };
        let degrees = surface.degrees();
        let mut sizes = [0usize; 2];
        for axis in 0..2 {
            let [a, b] = surface.range(axis);
            let n = surface
                .knots(axis)
                .0
                .iter()
                .filter(|&&k| k >= a && k <= b)
                .count();
            sizes[axis] = n
                .checked_sub(1)?
                .checked_mul(degrees[axis])?
                .checked_add(2)?;
        }
        let header =
            std::mem::size_of::<Self>() + 2 * std::mem::size_of::<Vec<f64>>() + EVALUATE_SCRATCH;
        let cuts_bytes = sizes[0].checked_add(sizes[1])?.checked_mul(8)?;
        let count = sizes[0].checked_mul(sizes[1])?;
        if header
            .checked_add(cuts_bytes)?
            .checked_add(count.checked_mul(std::mem::size_of::<Candidate>())?)?
            > self.maximum_bytes
        {
            return None;
        }
        let first = cuts(surface, 0, sizes[0])?;
        let first_bytes = first.capacity().checked_mul(8)?;
        // Check each actual capacity before the next optional allocation.
        if header
            .checked_add(first_bytes)?
            .checked_add(sizes[1].checked_mul(8)?)?
            .checked_add(count.checked_mul(std::mem::size_of::<Candidate>())?)?
            > self.maximum_bytes
        {
            return None;
        }
        let second = cuts(surface, 1, sizes[1])?;
        let cut_capacity = first_bytes.checked_add(second.capacity().checked_mul(8)?)?;
        let actual_count = first.len().checked_mul(second.len())?;
        if header
            .checked_add(cut_capacity)?
            .checked_add(actual_count.checked_mul(std::mem::size_of::<Candidate>())?)?
            > self.maximum_bytes
        {
            return None;
        }
        let cuts = [first, second];
        let mut grid = Vec::new();
        grid.try_reserve_exact(cuts[0].len().checked_mul(cuts[1].len())?)
            .ok()?;
        let actual = header.checked_add(cut_capacity)?.checked_add(
            grid.capacity()
                .checked_mul(std::mem::size_of::<Candidate>())?,
        )?;
        if actual > self.maximum_bytes {
            return None;
        }
        for &u in &cuts[0] {
            for &v in &cuts[1] {
                let uv = [u, v];
                grid.push(Candidate {
                    uv,
                    position: self.surface.evaluate(uv).0,
                });
            }
        }
        Some(grid)
    }
}
// Same open natural-knot cuts and arithmetic order as NurbsSurface::breaks
// with period=None. Use fallible, preflighted storage for the optional cache.
fn cuts(surface: &NurbsSurface, axis: usize, count: usize) -> Option<Vec<f64>> {
    let [a, b] = surface.range(axis);
    let degree = surface.degrees()[axis];
    let mut cuts = Vec::new();
    cuts.try_reserve_exact(count).ok()?;
    cuts.push(a);
    for pair in surface.knots(axis).0.windows(2) {
        if pair[0] < a || pair[1] > b {
            continue;
        }
        for piece in 0..degree {
            let shift = 0.;
            let cut = shift + pair[0] + (pair[1] - pair[0]) * piece as f64 / degree as f64;
            if cut > *cuts.last().unwrap_or(&a) + 1e-12 && cut < b - 1e-12 {
                cuts.push(cut);
            }
        }
    }
    cuts.push(b);
    Some(cuts)
}
#[cfg(test)]
mod tests;
