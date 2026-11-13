//! Paired authored thread material. This is not a mesh or exchange certificate.
//! Both laws retain their own true-period phase, hand, starts and finite length.
use crate::{authored_flat_thread::AuthoredThread, tapered_thread::interval::I};
use num_rational::BigRational;

mod caps;
mod classify;
pub use caps::{CapCell, CapCellBounds, CapDomains, CapEnd, CapKnot, CapOptions, Side};
#[cfg(test)]
mod tests;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Unresolved,
    Limit,
    Cancelled,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Classification {
    Inside,
    Outside,
    UnresolvedBoundary,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QueryOptions {
    pub volume_error_mm3: f64,
    /// Each actual linear source band consumes one cell; no sampling is used.
    pub maximum_cells: usize,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Volume {
    pub material_mm3: [f64; 2],
    pub outer_mm3: [f64; 2],
    pub cavity_mm3: [f64; 2],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Insert {
    outer: AuthoredThread,
    inner: AuthoredThread,
    minimum_wall_mm: f64,
}
impl Insert {
    /// Through cavity in the common source chart, with original phase origins.
    /// Unequal lengths need an independently qualified blind-floor graph.
    pub fn new(
        outer: AuthoredThread,
        inner: AuthoredThread,
        minimum_wall_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        checkpoint(interrupt)?;
        if !minimum_wall_mm.is_finite() || minimum_wall_mm <= 0. {
            return Err(Error::Invalid);
        }
        if outer.full_length_mm() != inner.full_length_mm() {
            return Err(Error::Unsupported);
        }
        // All operands are admitted finite literal source data. Binary-rational
        // comparison preserves exact equality, including nonbinary dimensions.
        let q = |x| BigRational::from_float(x).expect("sealed finite authored source");
        if q(outer.minimum_radius_mm()) - q(inner.maximum_radius_mm()) < q(minimum_wall_mm) {
            return Err(Error::Unsupported);
        }
        checkpoint(interrupt)?;
        Ok(Self {
            outer,
            inner,
            minimum_wall_mm,
        })
    }
    pub fn outer(&self) -> &AuthoredThread {
        &self.outer
    }
    pub fn inner(&self) -> &AuthoredThread {
        &self.inner
    }
    pub fn minimum_wall_mm(&self) -> f64 {
        self.minimum_wall_mm
    }
    pub fn full_length_mm(&self) -> f64 {
        self.outer.full_length_mm()
    }
    /// Exact complete source envelope, for separation only, never material.
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        let r = self.outer.maximum_radius_mm();
        [[-r, r], [-r, r], [0., self.full_length_mm()]]
    }
    pub fn separated_from_box(
        &self,
        other: [[f64; 2]; 3],
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<bool, Error> {
        checkpoint(interrupt)?;
        let other = valid_box(other)?;
        let own = self.bounds_mm();
        let result = (0..3).any(|a| own[a][1] < other[a].lo || other[a].hi < own[a][0]);
        checkpoint(interrupt)?;
        Ok(result)
    }
    /// Every angular cross-section spans a full true period. Integrating the
    /// actual normalized affine band yields width*(a²+ab+b²)/3. Starts cancel
    /// against their own lead, independently for the two retained sources.
    pub fn volume(
        &self,
        options: QueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Volume, Error> {
        checkpoint(interrupt)?;
        if !options.volume_error_mm3.is_finite() || options.volume_error_mm3 <= 0. {
            return Err(Error::Invalid);
        }
        if options.maximum_cells < 10 {
            return Err(Error::Limit);
        }
        let outer = integrated_volume(&self.outer, interrupt)?;
        let inner = integrated_volume(&self.inner, interrupt)?;
        let material = outer - inner;
        if !material.valid() || material.lo <= 0. || material.width() > options.volume_error_mm3 {
            return Err(Error::Unresolved);
        }
        let pair = |v: I| [v.lo, v.hi];
        let result = Volume {
            material_mm3: pair(material),
            outer_mm3: pair(outer),
            cavity_mm3: pair(inner),
        };
        checkpoint(interrupt)?;
        Ok(result)
    }
}
fn integrated_volume(source: &AuthoredThread, f: Option<&dyn Fn() -> bool>) -> Result<I, Error> {
    let mut sum = I::point(0.);
    for band in source.bands() {
        checkpoint(f)?;
        let [a, b] = band.radius_mm.map(I::point);
        let width = I::point(band.phase_mm[1]) - I::point(band.phase_mm[0]);
        sum = sum + width * (a.square() + a * b + b.square()) / I::point(3.);
    }
    let pi = I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    );
    Ok(pi * I::point(source.full_length_mm()) * sum / I::point(source.spec().pitch_mm))
}
fn valid_box(bounds: [[f64; 2]; 3]) -> Result<[I; 3], Error> {
    let p = bounds.map(|[lo, hi]| I::new(lo, hi));
    if p.into_iter().all(I::valid) {
        Ok(p)
    } else {
        Err(Error::Invalid)
    }
}
fn checkpoint(f: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if f.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
