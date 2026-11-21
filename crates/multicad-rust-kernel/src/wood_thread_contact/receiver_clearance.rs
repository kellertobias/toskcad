//! Original convex nut stock against an actual retained normal-sweep section.
//! Exterior cells are earned from original boundary owners, never family names.
use super::*;
use crate::normal_sweep::{Edge, Sweep};
mod bored;
mod cells;
mod screw;
mod seated_stock;
mod stock;
pub use bored::{
    BoredReceiverClearance, HardwarePiece, ReceivingBore, ToolBand, certify_round_pieces_clearance,
    profiled_screw_pieces,
};
use cells::Cell;
pub use screw::{
    ProfiledScrewPlacement, ProfiledScrewReceiverClearance,
    certify_profiled_screw_receiver_clearance,
};
pub use seated_stock::{
    SeatedStockPlacement, SeatedStockReceiverClearance, certify_seated_stock_receiver_clearance,
};
use stock::StockEnvelope;
#[derive(Clone, Copy)]
pub struct ReceivingSweep<'a> {
    pub source: &'a Sweep,
    pub frame: Frame,
    pub occurrence: Frame,
    pub reverse: bool,
}
#[derive(Clone, Copy)]
pub struct ReceiverClearanceOptions {
    pub maximum_overlap_volume_mm3: f64,
    /// Bounded rational primitive evaluations, including rejected candidates.
    pub maximum_work: usize,
}
pub struct Line8ReceiverClearance<'a> {
    nut: ClockedLine8Placement<'a>,
    receiver: ReceivingSweep<'a>,
    pub overlap_volume_upper_mm3: f64,
    pub used_exterior_cell: bool,
    pub boundary_owners: Vec<usize>,
    pub work_used: usize,
    pub accounted_peak_bytes: usize,
}
impl<'a> Line8ReceiverClearance<'a> {
    pub fn nut(&self) -> ClockedLine8Placement<'a> {
        self.nut
    }
    pub fn receiver(&self) -> ReceivingSweep<'a> {
        self.receiver
    }
}
struct Work<'a> {
    used: usize,
    maximum: usize,
    interrupt: &'a dyn Fn() -> bool,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> R<()> {
        checkpoint(self.interrupt)?;
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        if self.used > self.maximum {
            return Err(Error::Limit);
        }
        Ok(())
    }
}
/// The receipt bounds occupied overlap, including measured positive placement
/// residuals. It does not certify screw mating or the other receiving member.
pub fn certify_line8_receiver_clearance<'a>(
    nut: ClockedLine8Placement<'a>,
    receiver: ReceivingSweep<'a>,
    options: ReceiverClearanceOptions,
    interrupt: &dyn Fn() -> bool,
) -> R<Line8ReceiverClearance<'a>> {
    checkpoint(interrupt)?;
    if !options.maximum_overlap_volume_mm3.is_finite()
        || options.maximum_overlap_volume_mm3 < 0.
        || options.maximum_overlap_volume_mm3 > 1.
        || options.maximum_work == 0
        || options.maximum_work > 16_000_000
    {
        return Err(Error::Invalid);
    }
    let rings = receiver.source.rings();
    let outer = rings.first().ok_or(Error::Unsupported)?;
    if outer.len() > 256 {
        return Err(Error::Limit);
    }
    let mut work = Work {
        used: 0,
        maximum: options.maximum_work,
        interrupt,
    };
    work.charge(std::mem::size_of_val(rings))?;
    // Source scans and actual retained capacities are separate from scratch.
    // At most 2048 live guarded scalars (16,384 combined bits each),64 raw
    // <=131,072-bit arithmetic values, pinned GCD scratch and bounded vectors
    // fit the 16MiB reserve. No old scalar guard or caller work limit is renewed.
    let retained = receiver.source.retained_bytes() + nut.source.retained_bytes();
    let bytes = retained + 16 * 1024 * 1024 + outer.len() * 1024;
    if bytes > 64 * 1024 * 1024 {
        return Err(Error::Limit);
    }
    work.charge(rings.iter().map(Vec::len).sum::<usize>() * std::mem::size_of::<Edge>())?;
    let envelope = StockEnvelope::new(nut, receiver, &mut work)?;
    let mut best = envelope.outer_cover_bound(receiver, &mut work)?;
    let mut owners = Vec::new();
    if best > q(options.maximum_overlap_volume_mm3)? {
        for start in 0..outer.len() {
            work.charge(1)?;
            if let Some(cell) = Cell::from_boundary(outer, start, &mut work)? {
                let upper = envelope.cell_bound(&cell, &mut work)?;
                if upper < best {
                    best = upper;
                    owners = cell.owners;
                }
                if best <= q(options.maximum_overlap_volume_mm3)? {
                    break;
                }
            }
        }
    }
    if best > q(options.maximum_overlap_volume_mm3)? {
        return Err(Error::Unsupported);
    }
    let upper = if best.is_zero() {
        0.
    } else {
        // The exact comparison above independently proves the caller's binary
        // threshold itself is an upper bound. Tighten the outward cache with
        // that earned bound; equality preserves the existing volume semantics.
        enclosure(&best)?.hi.min(options.maximum_overlap_volume_mm3)
    };
    checkpoint(interrupt)?;
    Ok(Line8ReceiverClearance {
        nut,
        receiver,
        overlap_volume_upper_mm3: upper,
        used_exterior_cell: !owners.is_empty(),
        boundary_owners: owners,
        work_used: work.used,
        accounted_peak_bytes: bytes,
    })
}
#[cfg(test)]
mod tests;
