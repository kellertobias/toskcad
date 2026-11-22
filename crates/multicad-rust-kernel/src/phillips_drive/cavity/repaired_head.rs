//! Whole actual repaired-source containment, without a legacy profile inverse.
//! The local radial/floor margins are not Euclidean affine wall thicknesses.
use super::*;
use crate::head_profile::RepairedHeadMeridian;
const CELLS: usize = 65536;
const WORK: usize = 16_000_000;
const BYTES: usize = 64 * 1024 * 1024;
const SCRATCH: usize = 16 * 1024;
#[derive(Clone, Copy)]
pub struct RepairedHeadContainmentOptions {
    pub minimum_clearance_mm: f64,
    pub maximum_cells: usize,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for RepairedHeadContainmentOptions {
    fn default() -> Self {
        Self {
            minimum_clearance_mm: 1e-7,
            maximum_cells: CELLS,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
pub struct RepairedHeadContainment<'s, 'h, 'p> {
    source: &'s Cutter,
    head: &'h RepairedHeadMeridian<'p>,
    cells: usize,
    radial: f64,
    floor: f64,
    work: usize,
}
impl<'s, 'h, 'p> RepairedHeadContainment<'s, 'h, 'p> {
    pub fn source(&self) -> &'s Cutter {
        self.source
    }
    pub fn head(&self) -> &'h RepairedHeadMeridian<'p> {
        self.head
    }
    pub fn cells(&self) -> usize {
        self.cells
    }
    pub fn minimum_radial_clearance_mm(&self) -> f64 {
        self.radial
    }
    pub fn axial_floor_clearance_mm(&self) -> f64 {
        self.floor
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + std::mem::size_of::<Cutter>()
            + self.head.retained_budget_bytes()
    }
    /// Caller already validates cells<=CELLS. Shared-source aggregate queries
    /// count their own headers once and add this actual frontier/scratch bound.
    pub(crate) fn scratch_budget_bytes(cells: usize) -> usize {
        4 * (cells + 1) * std::mem::size_of::<Cell>() + SCRATCH
    }
}
#[derive(Clone, Copy)]
struct Cell {
    part: usize,
    parameter: [f64; 2],
    depth: [f64; 2],
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn disjoint(a: [f64; 2], b: I) -> bool {
    a[1] < b.lo || a[0] > b.hi
}
impl<'s> ClosedBoundaryCertificate<'s> {
    /// Certify containment of the entire closed cavity (except its common mouth
    /// disk) in this distinct repaired head. The mouth lies at the head's exact
    /// top station. Product cells cover every head parameter/depth pair; only
    /// certified axial disjointness removes cells. Otherwise a whole-source
    /// radius lower bound is compared to the entire cavity suffix upper bound.
    /// The head's positive monotone atlas supplies every intervening axial
    /// section; no old circle/rounded-station intersection is assumed.
    pub fn certify_repaired_head<'h, 'p>(
        &self,
        head: &'h RepairedHeadMeridian<'p>,
        options: RepairedHeadContainmentOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RepairedHeadContainment<'s, 'h, 'p>, Error> {
        checkpoint(interrupt)?;
        if !options.minimum_clearance_mm.is_finite()
            || options.minimum_clearance_mm <= 0.
            || options.maximum_cells == 0
            || options.maximum_cells > CELLS
            || options.maximum_work > WORK
            || options.maximum_bytes > BYTES
        {
            return Err(Error::Invalid);
        }
        let mut remaining = options
            .maximum_work
            .checked_sub(self.work_used())
            .and_then(|x| x.checked_sub(head.work_used()))
            .ok_or(Error::Limit)?;
        // Includes BOTH live prerequisite headers, their source values and the
        // immutable head Profile once. At most five initial roots and one extra
        // pending cell per binary split give <=4*(C+1) allocated Cell slots
        // under the pinned Global Vec geometric growth: for C=1,2 the at-most
        // seven pending values fit eight slots; C>=3 uses at most C+5 values,
        // hence capacity <=2*(C+5)<=4*(C+1). No custom allocator is used.
        let peak = std::mem::size_of::<RepairedHeadContainment<'_, '_, '_>>()
            + std::mem::size_of::<Self>()
            + std::mem::size_of::<Cutter>()
            + head.retained_budget_bytes()
            + RepairedHeadContainment::scratch_budget_bytes(options.maximum_cells);
        if peak > options.maximum_bytes {
            return Err(Error::Limit);
        }
        // Fixed station/option arithmetic and <=5 initial pending descriptors
        // are charged BEFORE constructing or allocating the query frontier.
        remaining = remaining.checked_sub(128).ok_or(Error::Limit)?;
        let source = self.source();
        let top = head.station(0).unwrap();
        let end = head.station(head.segment_count()).unwrap();
        let tip = source.envelope.tip;
        let floor = (I::point(end[1]) - I::point(top[1]) - tip).lo;
        if floor <= options.minimum_clearance_mm
            || (I::point(top[0]) - source.envelope.entry).lo <= options.minimum_clearance_mm
        {
            return Err(Error::Unresolved);
        }
        let mut pending = Vec::new();
        for part in 0..head.segment_count() {
            pending.push(Cell {
                part,
                parameter: [0., 1.],
                depth: [0., tip.hi],
            });
        }
        let mut cells = 0;
        let mut radial = f64::INFINITY;
        while let Some(cell) = pending.pop() {
            checkpoint(interrupt)?;
            if cells == options.maximum_cells {
                return Err(Error::Limit);
            }
            // Source evaluation <=257 fixed interval operations, suffix <=100,
            // and fixed bounds/splitting. Inherited head proof is counted once;
            // the sealed source evaluator inherits its constructor qualification.
            remaining = remaining.checked_sub(512).ok_or(Error::Limit)?;
            cells += 1;
            let z = I::point(top[1]) + I::new(cell.depth[0], cell.depth[1]);
            let stations = [
                head.station(cell.part).unwrap()[1],
                head.station(cell.part + 1).unwrap()[1],
            ];
            if disjoint(stations, z) {
                continue;
            }
            let radius = source.radial_suffix_bound_mm(cell.depth[0])?;
            let minimum = head
                .segment_certificate(cell.part)
                .unwrap()
                .minimum_radius_mm;
            let coarse = (I::point(minimum) - I::point(radius)).lo;
            if coarse > options.minimum_clearance_mm {
                radial = radial.min(coarse);
                continue;
            }
            // Root indices are built from this immutable head's actual count;
            // binary children retain0<=a<=b<=1. Its constructor already proves
            // whole-unit finite boxes; this typed handoff adds no new policy.
            let bounds = head.certified_interval_bounds(cell.part, cell.parameter);
            checkpoint(interrupt)?;
            if disjoint(bounds.position_mm[1], z) {
                continue;
            }
            let gap = (I::point(bounds.position_mm[0][0]) - I::point(radius)).lo;
            if gap > options.minimum_clearance_mm {
                radial = radial.min(gap);
                continue;
            }
            // Grid choice only. Both children retain the complete parameter ×
            // depth product; the same whole-source tests remain mandatory.
            let depth_width = (cell.depth[1] - cell.depth[0]) / tip.hi;
            let parameter_width = cell.parameter[1] - cell.parameter[0];
            let depth = depth_width >= parameter_width;
            let [a, b] = if depth { cell.depth } else { cell.parameter };
            let mid = a + (b - a) / 2.;
            if mid <= a || mid >= b {
                return Err(Error::Unresolved);
            }
            for range in [[mid, b], [a, mid]] {
                pending.push(if depth {
                    Cell {
                        depth: range,
                        ..cell
                    }
                } else {
                    Cell {
                        parameter: range,
                        ..cell
                    }
                });
            }
        }
        // At every cavity depth the monotone head atlas has an actual source
        // point, so at least one covering product cell cannot be discarded.
        // Successful exhaustion therefore took a nonempty minimum of finite
        // accepted gaps; radial is finite and strictly exceeds the request.
        checkpoint(interrupt)?;
        remaining = remaining.checked_sub(1).ok_or(Error::Limit)?;
        Ok(RepairedHeadContainment {
            source,
            head,
            cells,
            radial,
            floor,
            work: options.maximum_work - remaining,
        })
    }
}
#[cfg(test)]
mod tests;
