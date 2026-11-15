//! Finite conical seats on original, disjoint orthogonal through bores.
//! Protected full cylinders prove the whole cone against the unchanged stock.
use crate::{
    orthogonal::OrthogonalSolid,
    orthogonal_boundary::{self, Boundary},
    orthogonal_pockets::{self, OrthogonalThroughBore, PocketError},
    tapered_thread::interval::I,
};
use std::{f64::consts::PI, rc::Rc};
mod mesh;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Seat {
    pub bore: usize,
    pub end: bool,
    pub radius_mm: f64,
    pub depth_mm: f64,
}
#[derive(Clone, Copy)]
pub struct Options {
    pub clearance_mm: f64,
    pub minimum_feature_mm: f64,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Limit,
    Cancelled,
}
pub struct Seats {
    stock: Rc<OrthogonalSolid>,
    boundary: Boundary,
    bores: Rc<[OrthogonalThroughBore]>,
    seats: Vec<Seat>,
    volume: I,
    pub work_used: usize,
    pub accounted_peak_bytes: usize,
}
fn checkpoint(interrupt: &dyn Fn() -> bool) -> Result<(), Error> {
    if interrupt() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn failure(e: PocketError) -> Error {
    match e {
        PocketError::InvalidInput => Error::Invalid,
        PocketError::BudgetExceeded => Error::Limit,
        PocketError::UnresolvedBoundary => Error::Unsupported,
    }
}
fn depth_interval(bore: &OrthogonalThroughBore) -> I {
    let [a, b] = bore.centers_mm;
    let axis = (0..3)
        .find(|&i| a[i] != b[i])
        .expect("validated through bore");
    if a[axis] < b[axis] {
        I::point(b[axis]) - I::point(a[axis])
    } else {
        I::point(a[axis]) - I::point(b[axis])
    }
}
fn boundary_bytes(b: &Boundary) -> usize {
    std::mem::size_of::<Boundary>()
        + b.faces.capacity() * std::mem::size_of::<orthogonal_boundary::BoundaryFace>()
        + b.faces
            .iter()
            .map(|f| {
                f.rings.capacity() * std::mem::size_of::<Vec<[f64; 3]>>()
                    + f.rings
                        .iter()
                        .map(|r| r.capacity() * std::mem::size_of::<[f64; 3]>())
                        .sum::<usize>()
                    + f.quad_indices.capacity() * std::mem::size_of::<usize>()
            })
            .sum::<usize>()
}
impl Seats {
    pub fn new(
        stock: Rc<OrthogonalSolid>,
        bores: Rc<[OrthogonalThroughBore]>,
        seats: &[Seat],
        options: Options,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<Self, Error> {
        checkpoint(interrupt)?;
        if !options.clearance_mm.is_finite()
            || options.clearance_mm < 0.
            || !options.minimum_feature_mm.is_finite()
            || options.minimum_feature_mm <= 0.
            || options.maximum_work == 0
            || options.maximum_work > 16_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
            || seats.is_empty()
            || bores.is_empty()
        {
            return Err(Error::Invalid);
        }
        // A finite, deliberately narrow provider. Existing general Boolean and
        // original drilling APIs retain their independent admission limits.
        if stock.cells().len() > 64
            || stock.quads().len() > 256
            || bores.len() > 16
            || seats.len() > 2 * bores.len()
        {
            return Err(Error::Limit);
        }
        if stock.solid_count() != 1 {
            return Err(Error::Unsupported);
        }
        let q = stock.quads().len();
        // Conservative primitive bound for boundary construction and both full
        // original/protected cylinder validations, including rejected cylinders.
        let work = q * q * 32
            + (q + stock.cells().len() + bores.len()) * bores.len() * 2048
            + seats.len() * 128;
        if work > options.maximum_work {
            return Err(Error::Limit);
        }
        let initial = stock.retained_bytes()
            + std::mem::size_of_val(bores.as_ref())
            + 16 * 1024 * 1024
            + q * 8192;
        if initial > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let boundary = orthogonal_boundary::build(stock.quads()).map_err(|e| match e {
            orthogonal_boundary::BoundaryError::Limit => Error::Limit,
            _ => Error::Unsupported,
        })?;
        checkpoint(interrupt)?;
        orthogonal_pockets::validate_with_through_bores(
            &stock,
            &boundary,
            &[],
            &bores,
            options.clearance_mm,
            options.minimum_feature_mm,
        )
        .map_err(failure)?;
        let mut ordered = seats.to_vec();
        ordered.sort_by_key(|s| (s.bore, s.end));
        let mut protected = bores.to_vec();
        let mut depths = vec![[0.; 2]; bores.len()];
        for seat in &ordered {
            checkpoint(interrupt)?;
            let original = bores.get(seat.bore).ok_or(Error::Invalid)?;
            if !seat.radius_mm.is_finite()
                || !seat.depth_mm.is_finite()
                || (I::point(seat.radius_mm) - I::point(original.radius_mm)).lo
                    < options.minimum_feature_mm
                || seat.depth_mm < options.minimum_feature_mm
            {
                return Err(Error::Invalid);
            }
            let depth = &mut depths[seat.bore][usize::from(seat.end)];
            if *depth != 0. {
                return Err(Error::Invalid);
            }
            *depth = seat.depth_mm;
            protected[seat.bore].radius_mm = protected[seat.bore].radius_mm.max(seat.radius_mm);
        }
        for (bore, depths) in bores.iter().zip(&depths) {
            let depth = depth_interval(bore) - I::point(depths[0]) - I::point(depths[1]);
            if depth.lo <= options.minimum_feature_mm {
                return Err(Error::Unsupported);
            }
        }
        // Larger complete cylinders contain all cone points. Their original
        // entry/exit face ownership and pairwise separation remain mandatory.
        orthogonal_pockets::validate_with_through_bores(
            &stock,
            &boundary,
            &[],
            &protected,
            options.clearance_mm,
            options.minimum_feature_mm,
        )
        .map_err(failure)?;
        checkpoint(interrupt)?;
        let pi = I {
            lo: PI.next_down(),
            hi: PI.next_up(),
        };
        let mut volume = I::point(0.);
        for cell in stock.cells() {
            checkpoint(interrupt)?;
            volume = volume
                + (0..3).fold(I::point(1.), |v, i| {
                    v * (I::point(cell.max[i]) - I::point(cell.min[i]))
                });
        }
        for bore in bores.iter() {
            volume = volume
                - pi * I::point(bore.radius_mm) * I::point(bore.radius_mm) * depth_interval(bore);
        }
        for s in &ordered {
            let r = I::point(bores[s.bore].radius_mm);
            let outer = I::point(s.radius_mm);
            volume = volume
                - pi * I::point(s.depth_mm) * (outer - r) * (outer + r * I::point(2.))
                    / I::point(3.);
        }
        if !volume.lo.is_finite() || !volume.hi.is_finite() || volume.lo <= 0. {
            return Err(Error::Unsupported);
        }
        let bytes = stock.retained_bytes()
            + boundary_bytes(&boundary) * 2
            + std::mem::size_of_val(bores.as_ref())
            + ordered.capacity() * std::mem::size_of::<Seat>()
            + 16 * 1024 * 1024
            + q * 8192;
        if bytes > options.maximum_bytes {
            return Err(Error::Limit);
        }
        checkpoint(interrupt)?;
        Ok(Self {
            stock,
            boundary,
            bores,
            seats: ordered,
            volume,
            work_used: work,
            accounted_peak_bytes: bytes,
        })
    }
    pub fn stock(&self) -> &OrthogonalSolid {
        &self.stock
    }
    pub fn boundary(&self) -> &Boundary {
        &self.boundary
    }
    pub fn bores(&self) -> &[OrthogonalThroughBore] {
        &self.bores
    }
    pub fn seats(&self) -> &[Seat] {
        &self.seats
    }
    pub fn volume_bounds_mm3(&self) -> [f64; 2] {
        [self.volume.lo, self.volume.hi]
    }
}
