//! Exact cap proof partition. Partition events do not force mesh vertices.
use super::*;
use num_bigint::{BigInt, Sign};
use num_rational::BigRational as Q;
use num_traits::{One, ToPrimitive, Zero};
use std::{cmp::Ordering, mem::size_of};
#[cfg(test)]
mod tests;
const MAX_BITS: u64 = 8192;
// Temporary bounded Q operands/conversion buffers, separate from retained data.
const SCRATCH_BYTES: usize = 256 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Outer,
    Inner,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CapEnd {
    Start,
    End,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CapKnot {
    pub side: Side,
    pub source_start: u32,
    pub source_boundary: u8,
    pub event: usize,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CapCell {
    pub events: [usize; 2],
    pub outer_band: u8,
    pub inner_band: u8,
    pub outer_period: i64,
    pub inner_period: i64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapCellBounds {
    pub turns: [[f64; 2]; 2],
    pub outer_radius_mm: [[f64; 2]; 2],
    pub inner_radius_mm: [[f64; 2]; 2],
}
struct Angle {
    n: Box<[u8]>,
    d: Box<[u8]>,
}
struct Raw {
    angle: Angle,
    knot: Option<(Side, u32, u8)>,
}
/// A source-owned exact partition of a cap's angular domain. Source knots are
/// retained separately from the union partition; no proximity welding or forced
/// cross-ring vertex insertion follows from this object.
pub struct CapDomains {
    owner: Insert,
    end: CapEnd,
    angles: Box<[Angle]>,
    knots: Box<[CapKnot]>,
    cells: Box<[CapCell]>,
    retained_bytes: usize,
    construction_work: usize,
}
struct Work<'a> {
    used: usize,
    limit: usize,
    bytes: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        checkpoint(self.interrupt)?;
        if n > self.limit - self.used {
            return Err(Error::Limit);
        }
        self.used += n;
        Ok(())
    }
    fn memory(&self, n: usize) -> Result<(), Error> {
        if n > self.bytes {
            Err(Error::Limit)
        } else {
            Ok(())
        }
    }
}
fn work(options: CapOptions, interrupt: Option<&dyn Fn() -> bool>) -> Result<Work<'_>, Error> {
    let mut w = Work {
        used: 0,
        limit: options.maximum_work,
        bytes: options.maximum_bytes,
        interrupt,
    };
    w.charge(1)?;
    w.memory(SCRATCH_BYTES)?;
    Ok(w)
}
fn binary(x: f64) -> Q {
    Q::from_float(x).expect("sealed finite authored dimension")
}
#[derive(Clone, Copy)]
enum Op {
    Add,
    Sub,
    Mul,
    Div,
}
fn op(a: &Q, b: &Q, kind: Op, w: &mut Work<'_>) -> Result<Q, Error> {
    let [an, ad, bn, bd] = [
        a.numer().bits(),
        a.denom().bits(),
        b.numer().bits(),
        b.denom().bits(),
    ];
    let upper = match kind {
        Op::Add | Op::Sub => (an + bd).max(bn + ad).max(ad + bd) + 1,
        Op::Mul => (an + bn).max(ad + bd),
        Op::Div => (an + bd).max(ad + bn),
    };
    if upper > MAX_BITS {
        return Err(Error::Limit);
    }
    let words = upper.div_ceil(64) as usize;
    w.charge(1 + words * words)?;
    Ok(match kind {
        Op::Add => a + b,
        Op::Sub => a - b,
        Op::Mul => a * b,
        Op::Div => a / b,
    })
}
fn cmp(a: &Q, b: &Q, w: &mut Work<'_>) -> Result<Ordering, Error> {
    let difference = op(a, b, Op::Sub, w)?;
    Ok(difference.cmp(&Q::zero()))
}
fn modulo(a: &Q, period: &Q, w: &mut Work<'_>) -> Result<(Q, i64), Error> {
    let quotient = op(a, period, Op::Div, w)?;
    // Floor performs another integer division even when its result is small.
    // Charge the actual operand limb bound before invoking library arithmetic.
    let words = (quotient.numer().bits() + quotient.denom().bits()).div_ceil(64) as usize;
    w.charge(1 + words * words)?;
    // Source turn bounds keep the result within128 periods plus8starts;
    // the checked conversion below also prevents silent saturation.
    let k = quotient.floor().to_integer().to_i64().ok_or(Error::Limit)?;
    let shift = op(period, &Q::from_integer(k.into()), Op::Mul, w)?;
    Ok((op(a, &shift, Op::Sub, w)?, k))
}
impl Angle {
    fn encode(value: &Q, w: &mut Work<'_>) -> Result<Self, Error> {
        w.charge(1 + (value.numer().bits() + value.denom().bits()).div_ceil(64) as usize)?;
        let (_, n) = value.numer().to_bytes_be();
        let (_, d) = value.denom().to_bytes_be();
        Ok(Self {
            n: n.into_boxed_slice(),
            d: d.into_boxed_slice(),
        })
    }
    fn decode(&self, w: &mut Work<'_>) -> Result<Q, Error> {
        let words = (self.n.len() + self.d.len()).div_ceil(8);
        w.charge(1 + words * words)?;
        Ok(Q::new(
            BigInt::from_bytes_be(Sign::Plus, &self.n),
            BigInt::from_bytes_be(Sign::Plus, &self.d),
        ))
    }
    fn bytes(&self) -> usize {
        self.n.len() + self.d.len()
    }
}
fn reserve<T>(n: usize, w: &mut Work<'_>) -> Result<Vec<T>, Error> {
    w.charge(1)?;
    let mut v = Vec::new();
    v.try_reserve_exact(n).map_err(|_| Error::Limit)?;
    Ok(v)
}
fn station(owner: &Insert, end: CapEnd) -> Q {
    binary(match end {
        CapEnd::Start => 0.,
        CapEnd::End => owner.full_length_mm(),
    })
}
fn lead(source: &AuthoredThread, w: &mut Work<'_>) -> Result<Q, Error> {
    op(
        &binary(source.spec().pitch_mm),
        &Q::from_integer(source.spec().starts.into()),
        Op::Mul,
        w,
    )
}
fn phase(source: &AuthoredThread, z: &Q, turn: &Q, w: &mut Work<'_>) -> Result<Q, Error> {
    let delta = op(&lead(source, w)?, turn, Op::Mul, w)?;
    op(
        z,
        &delta,
        if source.spec().left_hand {
            Op::Add
        } else {
            Op::Sub
        },
        w,
    )
}
fn band_at(source: &AuthoredThread, z: &Q, turn: &Q, w: &mut Work<'_>) -> Result<(u8, i64), Error> {
    let q = phase(source, z, turn, w)?;
    let (q, k) = modulo(&q, &binary(source.spec().pitch_mm), w)?;
    for (i, b) in source.bands().into_iter().take(4).enumerate() {
        if cmp(&q, &binary(b.phase_mm[1]), w)? != Ordering::Greater {
            return Ok((i as u8, k));
        }
    }
    // Exact modulo is in[0,pitch), and sealed positive source bands cover it.
    // Failure of the first four upper comparisons selects the final band.
    Ok((4, k))
}
fn radius(
    source: &AuthoredThread,
    z: &Q,
    turn: &Q,
    band: u8,
    k: i64,
    w: &mut Work<'_>,
) -> Result<Q, Error> {
    let b = source.bands()[usize::from(band)];
    let offset = op(
        &binary(source.spec().pitch_mm),
        &Q::from_integer(k.into()),
        Op::Mul,
        w,
    )?;
    let q = op(&phase(source, z, turn, w)?, &offset, Op::Sub, w)?;
    let width = op(&binary(b.phase_mm[1]), &binary(b.phase_mm[0]), Op::Sub, w)?;
    let u = op(
        &op(&q, &binary(b.phase_mm[0]), Op::Sub, w)?,
        &width,
        Op::Div,
        w,
    )?;
    let dr = op(&binary(b.radius_mm[1]), &binary(b.radius_mm[0]), Op::Sub, w)?;
    op(
        &binary(b.radius_mm[0]),
        &op(&dr, &u, Op::Mul, w)?,
        Op::Add,
        w,
    )
}
fn enclosure(q: &Q, w: &mut Work<'_>) -> Result<[f64; 2], Error> {
    w.charge(1)?;
    let x = q.to_f64().ok_or(Error::Unresolved)?;
    // Validate the actual conversion against the retained exact quotient;
    // no assumed rounding quality of a general-purpose conversion is authority.
    let pair = [x.next_down(), x.next_up()];
    if !pair.into_iter().all(f64::is_finite)
        || cmp(&binary(pair[0]), q, w)? == Ordering::Greater
        || cmp(&binary(pair[1]), q, w)? == Ordering::Less
    {
        return Err(Error::Unresolved);
    }
    Ok(pair)
}
impl Insert {
    pub fn cap_domains(
        &self,
        end: CapEnd,
        options: CapOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<CapDomains, Error> {
        let mut w = work(options, interrupt)?;
        let count = 5 * (self.outer.spec().starts + self.inner.spec().starts) as usize;
        // The closed construction formulas below perform at most20 finite
        // binary conversions per event. Each has<=1075bit operands; reserve
        // their limb/GCD cost before any conversion, separately from exact ops.
        w.charge(20 * (count + 2) * 300)?;
        let planned_headers = size_of::<CapDomains>()
            + (count + 2) * size_of::<Raw>()
            + count * size_of::<CapKnot>()
            + (count + 2) * size_of::<Angle>()
            + (count + 1) * size_of::<CapCell>();
        w.memory(SCRATCH_BYTES + planned_headers + (count + 2) * 2 * 1024)?;
        let mut raw = reserve::<Raw>(count + 2, &mut w)?;
        let mut knots = reserve::<CapKnot>(count, &mut w)?;
        let mut angles = reserve::<Angle>(count + 2, &mut w)?;
        let mut cells = reserve::<CapCell>(count + 1, &mut w)?;
        let headers = size_of::<CapDomains>()
            + raw.capacity() * size_of::<Raw>()
            + knots.capacity() * size_of::<CapKnot>()
            + angles.capacity() * size_of::<Angle>()
            + cells.capacity() * size_of::<CapCell>();
        // Bound all compact payload before filling any arena; no retained
        // BigInt exists and each bounded rational has at most2*1024 bytes.
        w.memory(SCRATCH_BYTES + headers + (count + 2) * 2 * 1024)?;
        let z = station(self, end);
        for v in [Q::zero(), Q::one()] {
            raw.push(Raw {
                angle: Angle::encode(&v, &mut w)?,
                knot: None,
            });
        }
        for (side, source) in [(Side::Outer, &self.outer), (Side::Inner, &self.inner)] {
            for start in 0..source.spec().starts {
                for (boundary, b) in source.bands().into_iter().enumerate() {
                    let shift = op(
                        &binary(source.spec().pitch_mm),
                        &Q::from_integer(start.into()),
                        Op::Mul,
                        &mut w,
                    )?;
                    let q = op(&shift, &binary(b.phase_mm[0]), Op::Add, &mut w)?;
                    let delta = if source.spec().left_hand {
                        op(&q, &z, Op::Sub, &mut w)?
                    } else {
                        op(&z, &q, Op::Sub, &mut w)?
                    };
                    let angle = op(&delta, &lead(source, &mut w)?, Op::Div, &mut w)?;
                    let (angle, _) = modulo(&angle, &Q::one(), &mut w)?;
                    raw.push(Raw {
                        angle: Angle::encode(&angle, &mut w)?,
                        knot: Some((side, start, boundary as u8)),
                    });
                }
            }
        }
        // Bounded insertion sort charges actual comparisons and arithmetic;
        // cancellation never hides inside an infallible library comparator.
        for i in 1..raw.len() {
            let mut j = i;
            while j > 0 {
                let a = raw[j - 1].angle.decode(&mut w)?;
                let b = raw[j].angle.decode(&mut w)?;
                if cmp(&a, &b, &mut w)? != Ordering::Greater {
                    break;
                }
                raw.swap(j - 1, j);
                j -= 1;
            }
        }
        let mut previous = None;
        for entry in raw {
            let value = entry.angle.decode(&mut w)?;
            let distinct = if let Some(prev) = &previous {
                cmp(prev, &value, &mut w)? != Ordering::Equal
            } else {
                true
            };
            if distinct {
                angles.push(entry.angle);
                previous = Some(value);
            }
            if let Some((side, start, boundary)) = entry.knot {
                knots.push(CapKnot {
                    side,
                    source_start: start,
                    source_boundary: boundary,
                    event: angles.len() - 1,
                });
            }
        }
        for event in 0..angles.len() - 1 {
            let a = angles[event].decode(&mut w)?;
            let b = angles[event + 1].decode(&mut w)?;
            let midpoint = op(
                &op(&a, &b, Op::Add, &mut w)?,
                &Q::from_integer(2.into()),
                Op::Div,
                &mut w,
            )?;
            let (outer_band, outer_period) = band_at(&self.outer, &z, &midpoint, &mut w)?;
            let (inner_band, inner_period) = band_at(&self.inner, &z, &midpoint, &mut w)?;
            cells.push(CapCell {
                events: [event, event + 1],
                outer_band,
                inner_band,
                outer_period,
                inner_period,
            });
        }
        let retained_bytes = size_of::<CapDomains>()
            + angles.len() * size_of::<Angle>()
            + knots.len() * size_of::<CapKnot>()
            + cells.len() * size_of::<CapCell>()
            + angles.iter().map(Angle::bytes).sum::<usize>();
        // Shrinking the bounded Vec headers into boxed slices can briefly
        // retain old descriptors; the separate256KiB scratch includes these.
        w.memory(SCRATCH_BYTES + retained_bytes)?;
        w.charge(1)?;
        Ok(CapDomains {
            owner: *self,
            end,
            angles: angles.into_boxed_slice(),
            knots: knots.into_boxed_slice(),
            cells: cells.into_boxed_slice(),
            retained_bytes,
            construction_work: w.used,
        })
    }
}
impl CapDomains {
    pub fn owner(&self) -> &Insert {
        &self.owner
    }
    pub fn end(&self) -> CapEnd {
        self.end
    }
    pub fn knots(&self) -> &[CapKnot] {
        &self.knots
    }
    pub fn cells(&self) -> &[CapCell] {
        &self.cells
    }
    pub fn retained_bytes(&self) -> usize {
        self.retained_bytes
    }
    pub fn construction_work(&self) -> usize {
        self.construction_work
    }
    pub fn cell_bounds(
        &self,
        index: usize,
        options: CapOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<CapCellBounds, Error> {
        let mut w = work(options, interrupt)?;
        w.memory(SCRATCH_BYTES + self.retained_bytes)?;
        let cell = self.cells.get(index).ok_or(Error::Invalid)?;
        // Endpoint radii, station and six checked enclosures use at most60 conversions.
        w.charge(60 * 300)?;
        let turns = [
            self.angles[cell.events[0]].decode(&mut w)?,
            self.angles[cell.events[1]].decode(&mut w)?,
        ];
        let z = station(&self.owner, self.end);
        let mut result = CapCellBounds {
            turns: [[0.; 2]; 2],
            outer_radius_mm: [[0.; 2]; 2],
            inner_radius_mm: [[0.; 2]; 2],
        };
        for (i, t) in turns.iter().enumerate() {
            result.turns[i] = enclosure(t, &mut w)?;
            result.outer_radius_mm[i] = enclosure(
                &radius(
                    &self.owner.outer,
                    &z,
                    t,
                    cell.outer_band,
                    cell.outer_period,
                    &mut w,
                )?,
                &mut w,
            )?;
            result.inner_radius_mm[i] = enclosure(
                &radius(
                    &self.owner.inner,
                    &z,
                    t,
                    cell.inner_band,
                    cell.inner_period,
                    &mut w,
                )?,
                &mut w,
            )?;
        }
        w.charge(1)?;
        Ok(result)
    }
}
