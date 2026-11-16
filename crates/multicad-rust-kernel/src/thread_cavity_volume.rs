//! Shared numerical integration of explicit affine cavity source providers.
//! This module owns no geometry and never converts one source kind into another.
use crate::tapered_thread::interval::I;
use crate::thread_internal_treatment::{Enclosure, Error, QueryOptions};
pub(crate) struct Entry {
    pub radius: f64,
    pub depth: f64,
    pub slope: [f64; 2],
}
pub(crate) trait Source {
    fn phase_knots(&self) -> [f64; 6];
    fn profile(&self, band: usize, q: I) -> (I, I);
    fn pitch_mm(&self) -> f64;
    fn full_length_mm(&self) -> f64;
    fn full_length_interval(&self) -> I {
        I::point(self.full_length_mm())
    }
    fn prefix_length_mm(&self) -> f64;
    fn terminal_radius_mm(&self) -> Option<f64>;
    fn entry(&self) -> Option<Entry>;
}
struct Work<'a> {
    left: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn step(&mut self) -> Result<(), Error> {
        if self.interrupt.is_some_and(|check| check()) {
            return Err(Error::Cancelled);
        }
        if self.left == 0 {
            return Err(Error::Limit);
        }
        self.left -= 1;
        Ok(())
    }
}
/// Static source dispatch; reserve is a source-owned nonnegative volume bound.
/// The authored exact-turn provider uses zero, the literal rate provider its
/// independently derived angular-sliver reserve. It is not a caller tolerance.
pub(crate) fn integrate<S: Source>(
    source: &S,
    options: QueryOptions,
    reserve: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Enclosure, Error> {
    if !options.tolerance_mm3.is_finite() || options.tolerance_mm3 <= 0. {
        return Err(Error::Invalid);
    }
    if options.maximum_cells == 0 || options.maximum_cells > 1_000_000 {
        return Err(Error::Limit);
    }
    let mut work = Work {
        left: options.maximum_cells,
        interrupt,
    };
    work.step()?;
    let available = accuracy_budget(options, reserve)?;
    let pi = I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    );
    let factor = pi / I::point(source.pitch_mm());
    let mut sum = source.terminal_radius_mm().map_or(I::point(0.), |radius| {
        pi * I::point(radius).square()
            * (source.full_length_interval() - I::point(source.prefix_length_mm()))
    });
    // Both sealed providers qualify dimension, factor, derivative and slope
    // bounds before this entry point. Restricted cells preserve those finite
    // intervals; query options change subdivision only. Final certification
    // still requires a finite nonnegative enclosure within requested accuracy.
    let knots = source.phase_knots();
    let mut pending = Vec::with_capacity(64);
    for band in 0..5 {
        pending.push((knots[band], knots[band + 1], band, available / 10.));
    }
    while let Some((a, b, band, tolerance)) = pending.pop() {
        work.step()?;
        let value = cell(source, a, b, band, factor);
        if value.width() <= tolerance {
            sum = sum + value;
            continue;
        }
        let mid = a.midpoint(b);
        if mid == a || mid == b || tolerance == 0. {
            return Err(Error::Unresolved);
        }
        pending.push((a, mid, band, tolerance / 2.));
        pending.push((mid, b, band, tolerance / 2.));
    }
    let sum = I::new((sum.lo - reserve).next_down(), (sum.hi + reserve).next_up());
    if !sum.valid() || sum.lo < 0. || sum.width() > options.tolerance_mm3 {
        return Err(Error::Unresolved);
    }
    Ok(Enclosure {
        lower: sum.lo,
        upper: sum.hi,
    })
}
fn accuracy_budget(options: QueryOptions, reserve: f64) -> Result<f64, Error> {
    // Each private provider derives a finite nonnegative reserve from its
    // immutable source. Reserve consumes output width before any integration
    // cell is allocated; the caller's requested accuracy never becomes error.
    let available = if reserve == 0. {
        options.tolerance_mm3
    } else {
        (I::point(options.tolerance_mm3) - I::point(2.) * I::point(reserve)).lo
    };
    if available <= 0. {
        return Err(Error::Unresolved);
    }
    Ok(available)
}
fn column<S: Source>(source: &S, r: I) -> I {
    let f = I::point(source.prefix_length_mm());
    let Some(e) = source.entry() else {
        return r.square() * f;
    };
    let a = I::point(e.radius);
    let k = I::new(e.slope[0], e.slope[1]);
    let g = (a - r) / k;
    let g = I::new(g.lo.clamp(0., e.depth), g.hi.clamp(0., e.depth));
    r.square() * (f - g) + g * (a.square() - a * k * g + k.square() * g.square() / I::point(3.))
}
fn cell<S: Source>(source: &S, a: f64, b: f64, band: usize, factor: I) -> I {
    let width = I::point(b) - I::point(a);
    let midpoint = (I::point(a) + I::point(b)) / I::point(2.);
    let (r, _) = source.profile(band, midpoint);
    let (range, slope) = source.profile(band, I::new(a, b));
    let f = I::point(source.prefix_length_mm());
    let extra = source.entry().map_or(I::point(0.), |e| {
        I::point(range.magnitude()) / I::new(e.slope[0], e.slope[1])
    });
    // The integrated MAX column is C1 across g=0/c. Its derivative is
    // Lipschitz: H'=2r(F-g), a.e. H''<=2(F+r/k). This ordinary midpoint
    // bound is valid across a threshold; no high-order unsplit-kink rule.
    let second = I::point(2.) * slope.square() * (f + extra);
    let error = second * width.square() / I::point(24.);
    (column(source, r) + I::new(-error.hi, error.hi)) * width * factor
}
