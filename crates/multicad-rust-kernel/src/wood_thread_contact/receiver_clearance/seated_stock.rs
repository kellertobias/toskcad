//! Actual original cells enclose the seated stock material: cavities remove only.
use super::*;
use crate::orthogonal_bore_seats::Seats;
#[derive(Clone, Copy)]
pub struct SeatedStockPlacement<'a> {
    pub source: &'a Seats,
    pub frame: Frame,
}
pub struct SeatedStockReceiverClearance<'a> {
    stock: SeatedStockPlacement<'a>,
    receiver: ReceivingSweep<'a>,
    pub overlap_volume_upper_mm3: f64,
    /// Bounds retain the original source cell order; no cells are merged.
    pub source_cell_upper_mm3: Vec<f64>,
    pub work_used: usize,
    pub accounted_peak_bytes: usize,
}
impl<'a> SeatedStockReceiverClearance<'a> {
    pub fn stock(&self) -> SeatedStockPlacement<'a> {
        self.stock
    }
    pub fn receiver(&self) -> ReceivingSweep<'a> {
        self.receiver
    }
}
pub fn certify_seated_stock_receiver_clearance<'a>(
    stock: SeatedStockPlacement<'a>,
    receiver: ReceivingSweep<'a>,
    options: ReceiverClearanceOptions,
    interrupt: &dyn Fn() -> bool,
) -> R<SeatedStockReceiverClearance<'a>> {
    checkpoint(interrupt)?;
    if !options.maximum_overlap_volume_mm3.is_finite()
        || !(0. ..=1.).contains(&options.maximum_overlap_volume_mm3)
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
    let bytes = stock.source.retained_bytes()
        + receiver.source.retained_bytes()
        + 16 * 1024 * 1024
        + outer.len() * 1024;
    if bytes > 64 * 1024 * 1024 {
        return Err(Error::Limit);
    }
    let mut work = Work {
        used: 0,
        maximum: options.maximum_work,
        interrupt,
    };
    work.charge(
        rings.iter().map(Vec::len).sum::<usize>() * std::mem::size_of::<Edge>()
            + std::mem::size_of_val(stock.source.stock().cells()),
    )?;
    let threshold = q(options.maximum_overlap_volume_mm3)?;
    let mut total = Q::zero();
    let mut source_cell_upper_mm3 = Vec::with_capacity(stock.source.stock().cells().len());
    for cell in stock.source.stock().cells() {
        work.charge(64)?;
        let bounds = [
            [q(cell.min[0])?, q(cell.max[0])?],
            [q(cell.min[1])?, q(cell.max[1])?],
            [q(cell.min[2])?, q(cell.max[2])?],
        ];
        let envelope = StockEnvelope::source_box(stock.frame, bounds, receiver, &mut work)?;
        // The actual receiver's original support enclosure is sufficient for
        // this restricted external stock proof. No groove datum or catalog
        // designation supplies a void. More difficult arrangements refuse.
        let upper = envelope.outer_cover_bound(receiver, &mut work)?;
        total = bounded(total + &upper)?;
        if total > threshold {
            return Err(Error::Unsupported);
        }
        source_cell_upper_mm3.push(if upper.is_zero() {
            0.
        } else {
            enclosure(&upper)?
                .hi
                .min(options.maximum_overlap_volume_mm3)
        });
    }
    let upper = if total.is_zero() {
        0.
    } else {
        enclosure(&total)?
            .hi
            .min(options.maximum_overlap_volume_mm3)
    };
    checkpoint(interrupt)?;
    Ok(SeatedStockReceiverClearance {
        stock,
        receiver,
        overlap_volume_upper_mm3: upper,
        source_cell_upper_mm3,
        work_used: work.used,
        accounted_peak_bytes: bytes,
    })
}
