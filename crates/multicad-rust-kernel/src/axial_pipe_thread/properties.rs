use super::*;
use crate::tapered_thread::{Scalar, Series};
use std::{cmp::Ordering, collections::BinaryHeap};
type S = Series<13>;
fn c(x: f64) -> S {
    S::constant(x)
}
fn power(x: I, n: usize) -> I {
    (0..n).fold(I::point(1.), |a, _| a * x)
}
struct Work<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn step(&mut self) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(1).ok_or(Error::Limit)?;
        Ok(())
    }
}
fn positive_part(x: I) -> I {
    I::new(x.lo.max(0.), x.hi.max(0.))
}
impl Thread {
    fn a<const N: usize>(&self, band: usize, q: Series<N>) -> Series<N> {
        let k = self.spec.radial_taper;
        let shift = if k < 0. {
            Series::constant(k) * Series::constant(self.spec.length_mm)
        } else {
            Series::constant(0.)
        };
        match self.bands[band].curve {
            Curve::Line {
                phase_mm,
                radius_mm,
                slope,
            } => {
                let anchor = Series::constant(radius_mm)
                    - Series::constant(k) * Series::constant(phase_mm)
                    + shift;
                if slope == k {
                    anchor
                } else {
                    anchor
                        + (Series::constant(slope) - Series::constant(k))
                            * (q - Series::constant(phase_mm))
                }
            }
            Curve::Circle {
                phase_mm,
                radius_mm,
                circle_radius_mm: r,
                sign,
            } => {
                let c = Series::constant;
                let d = q - c(phase_mm);
                c(radius_mm) + c(sign) * (c(r) * c(r) - d * d).sqrt() - c(k) * q + shift
            }
        }
    }
    fn axial_integral(&self, band: usize, q: S, mode: u8) -> S {
        let a = self.a(band, q);
        let l = self.spec.length_mm;
        let k = self.spec.radial_taper.abs();
        match mode {
            0 => c(l) * a * a + c(k) * c(l) * c(l) * a + c(k) * c(k) * c(l) * c(l) * c(l) / c(3.),
            1 => {
                let limit = self.spec.crest_radius_limit_mm.expect("clipped mode");
                c(limit) * c(limit) * c(l)
            }
            _ => {
                let limit = self.spec.crest_radius_limit_mm.expect("clipped mode");
                c(limit) * c(limit) * c(l)
                    + (c(-2.) * c(limit) * c(limit) * c(limit) / c(3.) + c(limit) * c(limit) * a
                        - a * a * a / c(3.))
                        / c(k)
            }
        }
    }
    fn integrate(
        &self,
        band: usize,
        [lo, hi]: [I; 2],
        mode: u8,
        budget: f64,
        depth: u32,
        work: &mut Work<'_>,
    ) -> Result<I, Error> {
        work.step()?;
        let width = hi - lo;
        let mid = (lo + hi) * I::point(0.5);
        let half = width * I::point(0.5);
        let center = self
            .axial_integral(band, S::variable(mid), mode)
            .checked()?;
        let whole = self
            .axial_integral(band, S::variable(I::new(lo.lo, hi.hi)), mode)
            .checked()?;
        let mut result = I::point(0.);
        for n in (0..12).step_by(2) {
            result =
                result + center.c[n] * I::point(2.) * power(half, n + 1) / I::point((n + 1) as f64);
        }
        let remainder =
            I::point(2.) * I::point(whole.c[12].magnitude()) * power(positive_part(half), 13)
                / I::point(13.);
        result = result + I::new(-remainder.hi, remainder.hi);
        if result.width() <= budget {
            return result.checked();
        }
        if width.lo <= 0. || depth >= 12 {
            return Err(Error::Unresolved);
        }
        Ok(
            self.integrate(band, [lo, mid], mode, budget / 2., depth + 1, work)?
                + self.integrate(band, [mid, hi], mode, budget / 2., depth + 1, work)?,
        )
    }
    /// Volume of the complete finite radial material, including OD-clipped
    /// imperfect crests. Junction closure faces have zero volume, not omitted
    /// positive-width material. No decoded manufacturing mass is claimed.
    pub fn volume_bounds_mm3(
        &self,
        accuracy: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<[f64; 2], Error> {
        if !accuracy.is_finite() || accuracy <= 0. {
            return Err(Error::Invalid);
        }
        let pi = I::new(PI.next_down(), PI.next_up());
        let factor = pi / I::point(self.spec.pitch_mm);
        let budget = accuracy / (factor.hi * 32.);
        let mut work = Work {
            remaining: 8192,
            interrupt,
        };
        let mut sum = I::point(0.);
        for (band, b) in self.bands.iter().enumerate() {
            work.step()?;
            let [lo, hi] = b.phase_mm;
            let mut strips = Vec::<I>::new();
            if let (
                Some(limit),
                Curve::Line {
                    phase_mm: q,
                    radius_mm: r,
                    slope,
                },
            ) = (self.spec.crest_radius_limit_mm, b.curve)
                && slope != self.spec.radial_taper
            {
                let coefficient = I::point(slope) - I::point(self.spec.radial_taper);
                let origin = I::point(r) - I::point(self.spec.radial_taper) * I::point(q)
                    + if self.spec.radial_taper < 0. {
                        I::point(self.spec.radial_taper) * I::point(self.spec.length_mm)
                    } else {
                        I::point(0.)
                    };
                for threshold in [
                    I::point(limit),
                    I::point(limit)
                        - I::point(self.spec.radial_taper.abs()) * I::point(self.spec.length_mm),
                ] {
                    let root = I::point(q) + (threshold - origin) / coefficient;
                    let pad =
                        32. * f64::EPSILON * (1. + self.spec.pitch_mm + root.midpoint().abs());
                    let strip = I::new(
                        (root.lo - pad).next_down().max(lo),
                        (root.hi + pad).next_up().min(hi),
                    );
                    if strip.lo < strip.hi {
                        strips.push(strip);
                    }
                }
            }
            strips.sort_by(|a, b| a.lo.total_cmp(&b.lo));
            let mut cursor = lo;
            for strip in strips {
                if strip.lo > cursor {
                    sum = sum
                        + self.integral_piece(
                            band,
                            I::new(cursor, strip.lo),
                            false,
                            budget,
                            &mut work,
                        )?;
                }
                let end = cursor.max(strip.hi);
                sum = sum
                    + self.integral_piece(
                        band,
                        I::new(cursor.max(strip.lo), end),
                        true,
                        budget,
                        &mut work,
                    )?;
                cursor = end;
            }
            if cursor < hi {
                sum = sum
                    + self.integral_piece(band, I::new(cursor, hi), false, budget, &mut work)?;
            }
        }
        let result = (factor * sum).checked()?;
        if result.lo <= 0. || result.width() > accuracy {
            return Err(Error::Unresolved);
        }
        Ok(bounds(result))
    }
    fn integral_piece(
        &self,
        band: usize,
        domain: I,
        direct: bool,
        budget: f64,
        work: &mut Work<'_>,
    ) -> Result<I, Error> {
        work.step()?;
        let base = self.a(band, S::interval(domain)).c[0];
        let end = base + I::point(self.spec.radial_taper.abs()) * I::point(self.spec.length_mm);
        let mode = if let Some(limit) = self.spec.crest_radius_limit_mm {
            if base.lo >= limit {
                Some(1)
            } else if end.hi <= limit {
                Some(0)
            } else if base.hi <= limit && end.lo >= limit {
                Some(2)
            } else {
                None
            }
        } else {
            Some(0)
        };
        if let Some(mode) = mode.filter(|_| !direct) {
            self.integrate(
                band,
                [I::point(domain.lo), I::point(domain.hi)],
                mode,
                budget,
                0,
                work,
            )
        } else {
            // Root strips are numerical integration domains, never removed or
            // merged source geometry. Include their entire positive material.
            let radius = ceiling(
                base + I::point(self.spec.radial_taper.abs()) * I::new(0., self.spec.length_mm),
                self.spec.crest_radius_limit_mm,
            );
            (radius.square()
                * I::point(self.spec.length_mm)
                * (I::point(domain.hi) - I::point(domain.lo)))
            .checked()
        }
    }
    /// Bounded support of actual finite material, including both radial caps
    /// and the axis. Unresolved or exhausted searches never return an estimate.
    pub fn support_bounds(
        &self,
        direction: [f64; 3],
        accuracy: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<[f64; 2], Error> {
        if direction.into_iter().any(|x| !x.is_finite()) || !accuracy.is_finite() || accuracy <= 0.
        {
            return Err(Error::Invalid);
        }
        let mut work = Work {
            remaining: 32768,
            interrupt,
        };
        work.step()?;
        if direction == [0.; 3] {
            return Ok([0.; 2]);
        }
        if direction[0] == 0. && direction[1] == 0. {
            let cap = (I::point(direction[2]) * I::point(self.spec.length_mm)).checked()?;
            let result = I::new(cap.lo.max(0.), cap.hi.max(0.));
            if result.width() > accuracy {
                return Err(Error::Unresolved);
            }
            return Ok(bounds(result));
        }
        let mut lower = (I::point(direction[2]) * I::point(self.spec.length_mm))
            .lo
            .max(0.);
        let mut heap = BinaryHeap::new();
        let mut serial = 0;
        for (i, b) in self.bands.iter().enumerate() {
            heap.push(self.support_cell(
                i,
                I::new(b.phase_mm[0], b.phase_mm[1]),
                I::new(0., self.spec.length_mm),
                direction,
                &mut lower,
                serial,
            )?);
            serial += 1;
        }
        loop {
            let cell = heap
                .pop()
                .expect("five initial cells; every subdivision replaces one with two");
            work.step()?;
            if cell.upper - lower <= accuracy {
                return Ok([lower, cell.upper.max(lower)]);
            }
            let split_q = cell.q.width() >= cell.z.width();
            let halves = if split_q {
                cell.q.split()?
            } else {
                cell.z.split()?
            };
            for half in halves {
                let (q, z) = if split_q {
                    (half, cell.z)
                } else {
                    (cell.q, half)
                };
                let child = self.support_cell(cell.band, q, z, direction, &mut lower, serial)?;
                serial += 1;
                heap.push(child);
            }
        }
    }
    fn support_cell(
        &self,
        band: usize,
        q: I,
        z: I,
        d: [f64; 3],
        lower: &mut f64,
        serial: usize,
    ) -> Result<Cell, Error> {
        let dot = |p: [[f64; 2]; 3]| {
            (0..3).fold(I::point(0.), |v, i| {
                v + I::point(d[i]) * I::new(p[i][0], p[i][1])
            })
        };
        let upper = dot(self.surface_point_bounds(band, bounds(q), bounds(z))?)
            .checked()?
            .hi;
        let point = dot(self.surface_point_bounds(band, [q.midpoint(); 2], [z.midpoint(); 2])?)
            .checked()?;
        *lower = lower.max(point.lo);
        Ok(Cell {
            band,
            q,
            z,
            upper,
            serial,
        })
    }
}
struct Cell {
    band: usize,
    q: I,
    z: I,
    upper: f64,
    serial: usize,
}
impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Cell {}
impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Cell {
    fn cmp(&self, other: &Self) -> Ordering {
        self.upper
            .total_cmp(&other.upper)
            .then(self.serial.cmp(&other.serial))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn support_queue_orders_actual_bounds_and_exhaustion_refuses() {
        let t = Thread::new(
            Spec {
                major_radius_mm: 5.,
                pitch_mm: 1.,
                length_mm: 5.,
                radial_taper: 1. / 32.,
                section: Section::NptFlat {
                    depth_mm: 0.5,
                    flat_width_mm: 0.1,
                },
                crest_radius_limit_mm: None,
            },
            1e-7,
            1e-6,
        )
        .unwrap();
        let mut lower = 0.;
        let a = t
            .support_cell(
                0,
                I::new(0., 0.05),
                I::new(0., 1.),
                [0., 0., 1.],
                &mut lower,
                0,
            )
            .unwrap();
        let b = t
            .support_cell(
                0,
                I::new(0., 0.05),
                I::new(3., 4.),
                [0., 0., 1.],
                &mut lower,
                1,
            )
            .unwrap();
        assert!(a < b);
        assert!(a == a);
        assert!(a != b);
        let mut heap = BinaryHeap::from([a, b]);
        assert_eq!(heap.pop().unwrap().serial, 1);
        assert_eq!(
            Work {
                remaining: 0,
                interrupt: None
            }
            .step(),
            Err(Error::Limit)
        );
    }
}
