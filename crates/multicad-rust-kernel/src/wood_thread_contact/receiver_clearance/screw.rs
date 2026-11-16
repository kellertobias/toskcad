//! Complete original profiled screw envelope; recesses only remove material.
use super::*;
use crate::{head_profile::Curve, profiled_screw::Screw, screw_recess::Recess};
#[derive(Clone, Copy)]
pub struct ProfiledScrewPlacement<'a> {
    pub source: &'a Screw,
    pub frame: Frame,
}
pub struct ProfiledScrewReceiverClearance<'a> {
    screw: ProfiledScrewPlacement<'a>,
    receiver: ReceivingSweep<'a>,
    pub overlap_volume_upper_mm3: f64,
    /// One list per original meridian segment, followed by the full thread.
    pub boundary_owners: Vec<Vec<usize>>,
    pub work_used: usize,
    pub accounted_peak_bytes: usize,
}
impl<'a> ProfiledScrewReceiverClearance<'a> {
    pub fn screw(&self) -> ProfiledScrewPlacement<'a> {
        self.screw
    }
    pub fn receiver(&self) -> ReceivingSweep<'a> {
        self.receiver
    }
}
/// The full source material is contained in the union of original segment
/// envelopes and the complete original thread support. No visible/hidden band,
/// recess or catalog recipe is substituted for those original sources.
pub fn certify_profiled_screw_receiver_clearance<'a>(
    screw: ProfiledScrewPlacement<'a>,
    receiver: ReceivingSweep<'a>,
    options: ReceiverClearanceOptions,
    interrupt: &dyn Fn() -> bool,
) -> R<ProfiledScrewReceiverClearance<'a>> {
    checkpoint(interrupt)?;
    if !options.maximum_overlap_volume_mm3.is_finite()
        || !(0. ..=1.).contains(&options.maximum_overlap_volume_mm3)
        || options.maximum_work == 0
        || options.maximum_work > 16_000_000
    {
        return Err(Error::Invalid);
    }
    // This finite-stack implementation qualifies the actual flat-hex source.
    // Other recess families need their retained allocation ledger separately.
    if !matches!(screw.source.recess(), Recess::FlatHex(_)) {
        return Err(Error::Unsupported);
    }
    let rings = receiver.source.rings();
    let outer = rings.first().ok_or(Error::Unsupported)?;
    if outer.len() > 256 {
        return Err(Error::Limit);
    }
    let bytes = receiver.source.retained_bytes()
        + std::mem::size_of::<Screw>()
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
            + std::mem::size_of::<Screw>(),
    )?;
    let error = q(screw.source.profile().rounding_error_mm())?;
    let plain = q(screw.source.spec().shank_length_mm)?;
    let mut boxes = Vec::with_capacity(6);
    for segment in screw.source.profile().segments() {
        work.charge(64)?;
        let mut radius = q(segment.start[0].max(segment.end[0]))?;
        let mut lo = q(segment.start[1].min(segment.end[1]))?;
        let mut hi = q(segment.start[1].max(segment.end[1]))?;
        if let Curve::Arc {
            center, radius_mm, ..
        } = segment.curve
        {
            let r = q(radius_mm)?;
            radius = radius.max(bounded(q(center[0])? + &r)?);
            lo = lo.min(bounded(q(center[1])? - &r)?);
            hi = hi.max(bounded(q(center[1])? + &r)?);
        }
        radius = bounded(radius + &error)?;
        boxes.push([
            [-radius.clone(), radius],
            [-Q::zero(), Q::zero()],
            [
                bounded(lo - &plain - &error)?,
                bounded(hi - &plain + &error)?,
            ],
        ]);
        let last = boxes.last_mut().unwrap();
        last[1] = last[0].clone();
    }
    let thread = screw.source.thread().source_coordinate_bounds_mm();
    boxes.push([
        [q(thread[0][0])?, q(thread[0][1])?],
        [q(thread[1][0])?, q(thread[1][1])?],
        [q(thread[2][0])?, q(thread[2][1])?],
    ]);
    let mut total = Q::zero();
    let mut boundary_owners = Vec::with_capacity(boxes.len());
    let threshold = q(options.maximum_overlap_volume_mm3)?;
    let target = bounded(&threshold / Q::from_integer(boxes.len().into()))?;
    for bounds in boxes {
        let envelope = StockEnvelope::source_box(screw.frame, bounds, receiver, &mut work)?;
        let mut best = envelope.outer_cover_bound(receiver, &mut work)?;
        let mut owners = Vec::new();
        if best > threshold {
            for start in 0..outer.len() {
                if best <= target {
                    break;
                }
                work.charge(1)?;
                if let Some(cell) = Cell::from_boundary(outer, start, &mut work)? {
                    let bound = envelope.cell_bound(&cell, &mut work)?;
                    if bound < best {
                        best = bound;
                        owners = cell.owners;
                    }
                }
            }
            let ranges = envelope.projected_ranges(&mut work)?;
            for first in 0..outer.len() {
                if best <= target {
                    break;
                }
                for second in first + 1..outer.len() {
                    if best <= target {
                        break;
                    }
                    if let Some(cell) = Cell::mouth(outer, first, second, &ranges, &mut work)? {
                        let bound = envelope.cell_bound(&cell, &mut work)?;
                        if bound < best {
                            best = bound;
                            owners = cell.owners;
                        }
                    }
                }
            }
        }
        total = bounded(total + best)?;
        if total > threshold {
            return Err(Error::Unsupported);
        }
        boundary_owners.push(owners);
    }
    let upper = if total.is_zero() {
        0.
    } else {
        enclosure(&total)?
            .hi
            .min(options.maximum_overlap_volume_mm3)
    };
    checkpoint(interrupt)?;
    Ok(ProfiledScrewReceiverClearance {
        screw,
        receiver,
        overlap_volume_upper_mm3: upper,
        boundary_owners,
        work_used: work.used,
        accounted_peak_bytes: bytes,
    })
}
