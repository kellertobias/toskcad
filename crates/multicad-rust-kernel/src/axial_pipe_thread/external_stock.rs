//! Distinct derived NPT stock-cut material. The original clipped band is
//! radially mapped onto the retained EXACT rational stock radius. Original
//! pitch/phase metadata and all owners stay; the measured map transfer is a
//! separate policy gate. Untouched stock datums remain authoritative.
use super::authored::{
    StockClipBand, StockClipCurve, StockClipFace, StockClipOptions, StockClipPlan, StockClipRail,
    StockClipStation, StockClipVertex,
};
use super::{Error, I};
use crate::{Rational, cylinder::AxialCylinder};
use num_rational::BigRational;
use num_traits::{Signed, ToPrimitive};
use std::mem::size_of;
mod measures;
mod mesh;
mod query;
mod topology;
pub use measures::{FaceMeasure, Volume};
pub use mesh::MeshCandidate;
pub use query::{Point, Query};
pub use topology::{Boundary, Coedge, Curve, Face, FaceLoop, Vertex};
#[derive(Clone, Copy, Debug)]
pub struct ExternalStockCut {
    thread: StockClipPlan,
    stock: AxialCylinder,
    radius: I,
    center: [I; 2],
    bottom: I,
    top: I,
    join: I,
    tail: I,
    radial: I,
    transfer: f64,
    work: usize,
}
struct Work<'a> {
    used: usize,
    left: usize,
    limit: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl<'a> Work<'a> {
    fn new(
        options: StockClipOptions,
        retained: usize,
        interrupt: Option<&'a dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        if options.maximum_work == 0
            || options.maximum_work > 16_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        if retained
            .checked_add(16 * 1024)
            .is_none_or(|v| v > options.maximum_bytes)
        {
            return Err(Error::Limit);
        }
        Ok(Self {
            used: 0,
            left: options.maximum_work,
            limit: options.maximum_bytes,
            interrupt,
        })
    }
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.left = self.left.checked_sub(n).ok_or(Error::Limit)?;
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
fn exact(value: Rational) -> BigRational {
    let (n, d) = value.ratio();
    BigRational::new(n.into(), d.into())
}
fn enclosure(value: &BigRational) -> Result<I, Error> {
    let x = value.to_f64().ok_or(Error::Unresolved)?;
    let rounded = BigRational::from_float(x).ok_or(Error::Unresolved)?;
    let bounds = match value.cmp(&rounded) {
        std::cmp::Ordering::Less => I::new(x.next_down(), x),
        std::cmp::Ordering::Equal => I::point(x),
        std::cmp::Ordering::Greater => I::new(x, x.next_up()),
    };
    let lo = BigRational::from_float(bounds.lo).ok_or(Error::Unresolved)?;
    let hi = BigRational::from_float(bounds.hi).ok_or(Error::Unresolved)?;
    if &lo > value || &hi < value {
        return Err(Error::Unresolved);
    }
    Ok(bounds)
}
impl ExternalStockCut {
    /// A distinct radially derived threaded band plus the exact original stock
    /// tail. λ=R_stock_exact/R_original_binary is an exact ratio, not equality
    /// of display caches. Native admission must check the measured transfer.
    pub fn new(
        thread: StockClipPlan,
        stock: AxialCylinder,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut work = Work::new(options, size_of::<Self>() + 64 * 1024, interrupt)?;
        work.charge(thread.work_used())?;
        let r = exact(stock.circle().radius());
        let old = BigRational::from_float(thread.stock_radius_limit_mm()).ok_or(Error::Invalid)?;
        let height = exact(stock.height());
        let cut =
            BigRational::from_float(thread.source().spec().length_mm).ok_or(Error::Invalid)?;
        // Rational components<=60bits and exact finite-binary datums<=1075bits
        // keep normalized values below1203bits and comparison products below
        // 4096bits.64liveBigInts plus pinned scratch fit reserved64KiB.
        let bits = (old.numer().bits() + old.denom().bits())
            .max(cut.numer().bits() + cut.denom().bits()) as usize
            + 128;
        let words = bits.div_ceil(32);
        work.charge(4096 * (1 + words * words))?;
        if height <= cut
            || height > BigRational::from_integer(10_000_000.into())
            || r > BigRational::from_integer(10_000_000.into())
        {
            return Err(Error::Invalid);
        }
        let center = stock
            .circle()
            .center()
            .map(|v| enclosure(&exact(v)))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        let bottom = exact(stock.bottom());
        let top = exact(stock.top());
        let join = &bottom + &cut;
        let tail = &height - &cut;
        let radius = enclosure(&r)?;
        let radial = enclosure(&(&r / &old))?;
        let tail = enclosure(&tail)?;
        if radius.lo <= 0. || radial.lo <= 0. || tail.lo <= 0. {
            return Err(Error::Unresolved);
        }
        let transfer = enclosure(&(&r - &old).abs())?.hi;
        work.charge(1)?;
        Ok(Self {
            thread,
            stock,
            radius,
            center: [center[0], center[1]],
            bottom: enclosure(&bottom)?,
            top: enclosure(&top)?,
            join: enclosure(&join)?,
            tail,
            radial,
            transfer,
            work: work.used,
        })
    }
    /// Unchanged original authoring metadata; actual thread material uses the
    /// separately retained radial map and stock translation.
    pub fn thread(&self) -> &StockClipPlan {
        &self.thread
    }
    pub fn stock(&self) -> &AxialCylinder {
        &self.stock
    }
    pub fn stock_axis_base_bounds_mm(&self) -> [[f64; 2]; 3] {
        let p = [self.center[0], self.center[1], self.bottom];
        p.map(|v| [v.lo, v.hi])
    }
    /// Rounded display cache only; exact datums are available through stock().
    pub fn stock_length_mm(&self) -> f64 {
        self.stock.height().as_f64()
    }
    pub fn stock_radius_mm(&self) -> f64 {
        self.stock.circle().radius().as_f64()
    }
    pub fn derived_major_radius_bounds_mm(&self) -> [f64; 2] {
        let r = I::point(self.thread.source().spec().major_radius_mm) * self.radial;
        [r.lo, r.hi]
    }
    pub fn stock_radius_bounds_mm(&self) -> [f64; 2] {
        [self.radius.lo, self.radius.hi]
    }
    pub fn stock_top_bounds_mm(&self) -> [f64; 2] {
        [self.top.lo, self.top.hi]
    }
    pub fn stock_bottom_bounds_mm(&self) -> [f64; 2] {
        [self.bottom.lo, self.bottom.hi]
    }
    pub fn thread_join_bounds_mm(&self) -> [f64; 2] {
        [self.join.lo, self.join.hi]
    }
    pub fn radial_scale_bounds(&self) -> [f64; 2] {
        [self.radial.lo, self.radial.hi]
    }
    pub fn original_source_position_transfer_mm(&self) -> f64 {
        self.transfer
    }
    /// All original normalized patch/coedge derivatives. Translation is the
    /// common stock placement; pitch and phase are unchanged by the radial map.
    pub fn original_source_derivative_transfer_mm(&self) -> f64 {
        let tau = I::new(
            std::f64::consts::TAU.next_down(),
            std::f64::consts::TAU.next_up(),
        );
        let turns = I::point(self.thread.source().spec().length_mm)
            / I::point(self.thread.source().spec().pitch_mm);
        (I::point(self.transfer) * (I::point(1.) + tau * (I::point(2.) + turns))).hi
    }
    pub fn tail_length_bounds_mm(&self) -> [f64; 2] {
        [self.tail.lo, self.tail.hi]
    }
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        let bounds = self.stock.bounds();
        std::array::from_fn(|axis| {
            [
                enclosure(&exact(bounds[0][axis]))
                    .expect("admitted exact stock")
                    .lo,
                enclosure(&exact(bounds[1][axis]))
                    .expect("admitted exact stock")
                    .hi,
            ]
        })
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        size_of::<Self>()
    }
}
#[cfg(test)]
mod tests;
