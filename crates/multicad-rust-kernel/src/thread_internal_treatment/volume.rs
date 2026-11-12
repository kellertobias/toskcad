//! Finite source-column integration, independent of any tessellation.
use super::*;
use crate::tapered_thread::interval::I;
fn checked(value: I) -> Result<I, Error> {
    if value.valid() {
        Ok(value)
    } else {
        Err(Error::Unresolved)
    }
}
fn clamp(value: I, lower: f64, upper: f64) -> I {
    I::new(value.lo.clamp(lower, upper), value.hi.clamp(lower, upper))
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
impl Cavity {
    pub fn removed_volume_mm3(&self, options: QueryOptions) -> Result<Enclosure, Error> {
        self.removed_volume_with_interrupt(options, None)
    }
    pub fn removed_volume_with_interrupt(
        &self,
        options: QueryOptions,
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
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let factor = pi / I::point(self.thread.spec().pitch_mm);
        let mut sum = if self.intent.terminal.is_some() {
            pi * I::point(self.terminal_radius).square()
                * (I::point(self.full_length_mm()) - I::point(self.terminal_station))
        } else {
            I::point(0.)
        };
        let knots = self.thread.profile_breaks();
        let mut pending = Vec::with_capacity(64);
        for band in 0..5 {
            pending.push((
                knots[band],
                knots[band + 1],
                band,
                options.tolerance_mm3 / 10.,
            ));
        }
        while let Some((a, b, band, tolerance)) = pending.pop() {
            work.step()?;
            let value = self.volume_cell(a, b, band, factor)?;
            if value.width() <= tolerance {
                sum = checked(sum + value)?;
                continue;
            }
            let mid = a.midpoint(b);
            if mid == a || mid == b || tolerance == 0. {
                return Err(Error::Unresolved);
            }
            pending.push((a, mid, band, tolerance / 2.));
            pending.push((mid, b, band, tolerance / 2.));
        }
        let sum = checked(sum)?;
        if sum.lo < 0. || sum.width() > options.tolerance_mm3 {
            return Err(Error::Unresolved);
        }
        Ok(Enclosure {
            lower: sum.lo,
            upper: sum.hi,
        })
    }
    pub(super) fn profile(&self, band: usize, q: I) -> (I, I) {
        let spec = self.thread.spec();
        let knots = self.thread.profile_breaks();
        if band == 0 || band == 4 {
            return (I::point(spec.major_radius_mm), I::point(0.));
        }
        if band == 2 {
            return (I::point(self.thread.root_radius_mm()), I::point(0.));
        }
        let anchor = if band == 1 { knots[1] } else { knots[4] };
        let slope = I::point(self.thread.slope_at_phase(anchor, band));
        (
            I::point(spec.major_radius_mm) + slope * (q - I::point(anchor)),
            slope,
        )
    }
    fn column(&self, r: I) -> Result<I, Error> {
        let f = I::point(self.terminal_station);
        let Some(c) = self.intent.entry_depth_mm else {
            return checked(r.square() * f);
        };
        let a = I::point(self.entry_radius);
        let [lo, hi] = self.entry_slope.expect("validated enabled cone");
        let k = I::new(lo, hi);
        let g = clamp(checked((a - r) / k)?, 0., c);
        // Two positive physical height columns. This remains correct at both
        // clamp endpoints, including a last-bit source flank/root discrepancy.
        checked(
            r.square() * (f - g)
                + g * (a.square() - a * k * g + k.square() * g.square() / I::point(3.)),
        )
    }
    fn volume_cell(&self, a: f64, b: f64, band: usize, factor: I) -> Result<I, Error> {
        let width = I::point(b) - I::point(a);
        let midpoint = (I::point(a) + I::point(b)) / I::point(2.);
        let (r, _) = self.profile(band, midpoint);
        let (range, slope) = self.profile(band, I::new(a, b));
        let f = I::point(self.terminal_station);
        let extra = if let Some([lo, hi]) = self.entry_slope {
            I::point(range.magnitude()) / I::new(lo, hi)
        } else {
            I::point(0.)
        };
        // H(r)=integral max(r,cone(z))² dz is C1 across both clamp thresholds.
        // H'=2r(F-g); a.e. H''<=2(F+r/k). Thus the first derivative is
        // Lipschitz and the ordinary midpoint error is valid across a kink;
        // no high-derivative quadrature is applied across an unsplit MAX.
        let second = checked(I::point(2.) * slope.square() * (f + extra))?;
        let error = second * width.square() / I::point(24.);
        checked((self.column(r)? + I::new(-error.hi, error.hi)) * width * factor)
    }
}
