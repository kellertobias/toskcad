//! Round hardware against a receiving sweep, optionally less a coaxial
//! stepped tool. The hardware is a union of axial cylinder pieces in its own
//! frame. Each piece's square envelope is subdivided until every box that
//! can hold hardware material is proved empty of receiving material: outside
//! the original stock, inside one convex exterior cell of its section, or
//! inside one convex band of the tool, which removed that material. A piece
//! coaxial with the tool is admitted exactly when it lies inside one band.
//! The receiving material is a subset of the original sweep, so every bound
//! of the original sweep also bounds it. Boxes left at the depth limit keep
//! their exact slab bounds in the total.
use super::*;
use crate::{head_profile::Curve, profiled_screw::Screw};

/// An axial cylinder of hardware material along the hardware frame's third
/// axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HardwarePiece {
    pub radius_mm: f64,
    pub span_mm: [f64; 2],
}

/// One band of the receiving tool along its axis from its origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ToolBand {
    pub start_mm: f64,
    pub end_mm: f64,
    pub radius_mm: f64,
}

/// A coaxial stepped tool that removed material from the receiving sweep,
/// in the receiver's body coordinates.
#[derive(Clone, Copy)]
pub struct ReceivingBore<'a> {
    pub origin_mm: [f64; 3],
    pub axis: [f64; 3],
    pub bands: &'a [ToolBand],
}

pub struct BoredReceiverClearance {
    pub overlap_volume_upper_mm3: f64,
    pub boxes: usize,
    pub work_used: usize,
}

const MAX_DEPTH: usize = 28;
const MAX_BOXES: usize = 200_000;

/// The original profiled screw as axial pieces: each meridian segment's
/// radial and axial bounds and the complete thread support.
pub fn profiled_screw_pieces(screw: &Screw) -> R<Vec<HardwarePiece>> {
    let error = screw.profile().rounding_error_mm();
    let plain = screw.spec().shank_length_mm;
    let mut pieces = Vec::new();
    for segment in screw.profile().segments() {
        let mut radius = segment.start[0].max(segment.end[0]);
        let mut lo = segment.start[1].min(segment.end[1]);
        let mut hi = segment.start[1].max(segment.end[1]);
        if let Curve::Arc {
            center,
            radius_mm,
            sweep,
        } = segment.curve
        {
            // Only the extremes the arc actually passes bound it.
            let start = (segment.start[1] - center[1]).atan2(segment.start[0] - center[0]);
            let passes = |angle: f64| {
                let t = ((angle - start) * sweep.signum()).rem_euclid(std::f64::consts::TAU);
                t <= sweep.abs()
            };
            if passes(0.) {
                radius = radius.max(center[0] + radius_mm);
            }
            if passes(-std::f64::consts::FRAC_PI_2) {
                lo = lo.min(center[1] - radius_mm);
            }
            if passes(std::f64::consts::FRAC_PI_2) {
                hi = hi.max(center[1] + radius_mm);
            }
        }
        // Outward rounding keeps the bounds conservative.
        pieces.push(HardwarePiece {
            radius_mm: (radius + error).next_up(),
            span_mm: [
                (lo - plain - error).next_down(),
                (hi - plain + error).next_up(),
            ],
        });
    }
    let thread = screw.thread().source_coordinate_bounds_mm();
    let reach = thread[0]
        .iter()
        .chain(&thread[1])
        .map(|v| v.abs())
        .fold(0., f64::max);
    // The thread box is radial: its half width bounds the major radius and
    // the retained root circle about the axis.
    pieces.push(HardwarePiece {
        radius_mm: reach.next_up(),
        span_mm: [thread[2][0], thread[2][1]],
    });
    if pieces
        .iter()
        .any(|p| !(p.radius_mm.is_finite() && p.radius_mm > 0. && p.span_mm[1] > p.span_mm[0]))
    {
        return Err(Error::Invalid);
    }
    Ok(pieces)
}

struct Band {
    range: [Q; 2],
    radius2: Q,
    radius: Q,
}

struct LocalBore {
    origin: [Q; 3],
    axis: [Q; 3],
    bands: Vec<Band>,
}

impl LocalBore {
    /// The band containing every corner, if any.
    fn contains(&self, corners: &[[Q; 3]], work: &mut Work<'_>) -> R<bool> {
        'band: for band in &self.bands {
            for p in corners {
                work.charge(24)?;
                let d: [Q; 3] = std::array::from_fn(|i| &p[i] - &self.origin[i]);
                let t = dot(&d, &self.axis)?;
                if t < band.range[0] || t > band.range[1] {
                    continue 'band;
                }
                let radial = bounded(dot(&d, &d)? - &t * &t)?;
                if radial > band.radius2 {
                    continue 'band;
                }
            }
            return Ok(true);
        }
        Ok(false)
    }
}

/// Bound hardware-material overlap with the receiving material.
pub fn certify_round_pieces_clearance(
    pieces: &[HardwarePiece],
    frame: Frame,
    receiver: ReceivingSweep<'_>,
    bore: Option<ReceivingBore<'_>>,
    options: ReceiverClearanceOptions,
    interrupt: &dyn Fn() -> bool,
) -> R<BoredReceiverClearance> {
    checkpoint(interrupt)?;
    if !options.maximum_overlap_volume_mm3.is_finite()
        || !(0. ..=1.).contains(&options.maximum_overlap_volume_mm3)
        || options.maximum_work == 0
        || options.maximum_work > 16_000_000
        || pieces.is_empty()
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
    work.charge(rings.iter().map(Vec::len).sum::<usize>() * std::mem::size_of::<Edge>())?;
    let placement = StockEnvelope::placement(frame, receiver, &mut work)?;
    // Exterior cells of the original section that do not depend on a box.
    let mut cells = Vec::new();
    for start in 0..outer.len() {
        if let Some(cell) = Cell::from_boundary(outer, start, &mut work)? {
            cells.push(cell);
        }
    }
    let own_extent = |edges: [usize; 2]| -> R<[[Q; 2]; 2]> {
        let mut extent: [[Q; 2]; 2] = std::array::from_fn(|_| [Q::zero(), Q::zero()]);
        let points: Vec<[f64; 2]> = edges
            .iter()
            .flat_map(|&i| [outer[i].start(), outer[i].end()])
            .collect();
        for axis in 0..2 {
            let values: Vec<Q> = points.iter().map(|p| q(p[axis])).collect::<R<_>>()?;
            extent[axis] = [
                values.iter().min().unwrap().clone(),
                values.iter().max().unwrap().clone(),
            ];
        }
        Ok(extent)
    };
    let mut mouths = Vec::new();
    for first in 0..outer.len() {
        for second in first + 1..outer.len() {
            if Cell::mouth(
                outer,
                first,
                second,
                &own_extent([first, second])?,
                &mut work,
            )?
            .is_some()
            {
                mouths.push((first, second));
            }
        }
    }
    let local = match bore {
        None => None,
        Some(bore) => {
            let origin = StockEnvelope::receiver_local(receiver, bore.origin_mm, true, &mut work)?;
            let axis = StockEnvelope::receiver_local(receiver, bore.axis, false, &mut work)?;
            if dot(&axis, &axis)? != Q::from_integer(1.into()) {
                return Err(Error::Unsupported);
            }
            let mut bands = Vec::new();
            for band in bore.bands {
                let radius = q(band.radius_mm)?;
                bands.push(Band {
                    range: [q(band.start_mm)?, q(band.end_mm)?],
                    radius2: bounded(&radius * &radius)?,
                    radius,
                });
            }
            Some(LocalBore {
                origin,
                axis,
                bands,
            })
        }
    };
    let threshold = q(options.maximum_overlap_volume_mm3)?;
    let mut total = Q::zero();
    let mut boxes = 0usize;
    let (origin, columns, _) = &placement;
    let piece_axis = columns[2].clone();
    for piece in pieces {
        let rho = q(piece.radius_mm)?;
        let span = [q(piece.span_mm[0])?, q(piece.span_mm[1])?];
        // A piece coaxial with the tool lies inside one band exactly.
        if let Some(bore) = &local
            && cross(&piece_axis, &bore.axis)?.iter().all(Q::is_zero)
        {
            let sign = dot(&piece_axis, &bore.axis)?;
            let d: [Q; 3] = std::array::from_fn(|i| &origin[i] - &bore.origin[i]);
            let base = dot(&d, &bore.axis)?;
            let offset2 = bounded(dot(&d, &d)? - &base * &base)?;
            let ends = [
                bounded(&base + &sign * &span[0])?,
                bounded(&base + &sign * &span[1])?,
            ];
            let [low, high] = [
                ends[0].clone().min(ends[1].clone()),
                ends[0].clone().max(ends[1].clone()),
            ];
            let inside = bore.bands.iter().any(|band| {
                let room = &band.radius - &rho;
                low >= band.range[0]
                    && high <= band.range[1]
                    && room >= Q::zero()
                    && offset2 <= &room * &room
            });
            if inside {
                continue;
            }
        }
        let mut stack = vec![(
            [
                [-rho.clone(), rho.clone()],
                [-rho.clone(), rho.clone()],
                span.clone(),
            ],
            0usize,
        )];
        let rho2 = bounded(&rho * &rho)?;
        while let Some((bounds, depth)) = stack.pop() {
            boxes += 1;
            if boxes > MAX_BOXES {
                return Err(Error::Limit);
            }
            work.charge(64)?;
            // The nearest point of the cross-section to the piece axis.
            let near = |range: &[Q; 2]| {
                if range[0] > Q::zero() {
                    range[0].clone()
                } else if range[1] < Q::zero() {
                    range[1].clone()
                } else {
                    Q::zero()
                }
            };
            let [nx, ny] = [near(&bounds[0]), near(&bounds[1])];
            if bounded(&nx * &nx + &ny * &ny)? >= rho2 {
                continue;
            }
            let envelope = StockEnvelope::placed_box(&placement, bounds.clone(), &mut work)?;
            let mut best = envelope.outer_cover_bound(receiver, &mut work)?;
            for cell in &cells {
                if best.is_zero() {
                    break;
                }
                best = best.min(envelope.cell_bound(cell, &mut work)?);
            }
            if !best.is_zero() && !mouths.is_empty() {
                let ranges = envelope.projected_ranges(&mut work)?;
                for &(first, second) in &mouths {
                    if best.is_zero() {
                        break;
                    }
                    if let Some(cell) = Cell::mouth(outer, first, second, &ranges, &mut work)? {
                        best = best.min(envelope.cell_bound(&cell, &mut work)?);
                    }
                }
            }
            if best.is_zero() {
                continue;
            }
            if let Some(bore) = &local
                && bore.contains(envelope.corners(), &mut work)?
            {
                continue;
            }
            if depth >= MAX_DEPTH {
                total = bounded(total + best)?;
                if total > threshold {
                    return Err(Error::Unsupported);
                }
                continue;
            }
            let extent: [Q; 3] = std::array::from_fn(|i| &bounds[i][1] - &bounds[i][0]);
            let cross_max = extent[0].clone().max(extent[1].clone());
            let axis = if extent[2] > Q::from_integer(2.into()) * &cross_max {
                2
            } else if extent[0] >= extent[1] {
                0
            } else {
                1
            };
            let middle =
                bounded((&bounds[axis][0] + &bounds[axis][1]) / Q::from_integer(2.into()))?;
            let mut low = bounds.clone();
            low[axis][1] = middle.clone();
            let mut high = bounds;
            high[axis][0] = middle;
            stack.push((low, depth + 1));
            stack.push((high, depth + 1));
        }
    }
    let upper = if total.is_zero() {
        0.
    } else {
        enclosure(&total)?
            .hi
            .min(options.maximum_overlap_volume_mm3)
    };
    checkpoint(interrupt)?;
    Ok(BoredReceiverClearance {
        overlap_volume_upper_mm3: upper,
        boxes,
        work_used: work.used,
    })
}
