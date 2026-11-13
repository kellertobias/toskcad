//! Whole-box queries of the continuous authored law, not legacy membership.
use super::*;
impl AuthoredMetricThread {
    pub fn classify_box(
        &self,
        bounds: [[f64; 2]; 3],
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Classification, Error> {
        use Classification::*;
        checkpoint(interrupt)?;
        let p = bounds.map(interval);
        if p.into_iter().any(|x| !x.valid()) {
            return Err(Error::Invalid);
        }
        let [x, y, z] = p;
        let rho = (x.square() + y.square()).sqrt();
        let result =
            if z.hi < 0. || z.lo > self.spec.length_mm || rho.lo > self.spec.major_radius_mm {
                Outside
            } else if !rho.valid() {
                UnresolvedBoundary
            } else {
                let r = self.phase_radius(z - self.rate() * angle(x, y))?;
                if rho.lo > r.hi {
                    Outside
                } else if z.lo > 0. && z.hi < self.spec.length_mm && rho.hi < r.lo {
                    Inside
                } else {
                    UnresolvedBoundary
                }
            };
        checkpoint(interrupt)?;
        Ok(result)
    }
    pub(super) fn phase_radius(&self, q: I) -> Result<I, Error> {
        let quotient = q / I::point(self.spec.pitch_mm);
        if !q.valid() || !quotient.valid() {
            return Err(Error::Unresolved);
        }
        if quotient.width() >= 1. {
            return Ok(I::new(
                self.minimum_radius_mm()[0],
                self.spec.major_radius_mm,
            ));
        }
        let knots = self.knots();
        let mut result = I::new(f64::INFINITY, f64::NEG_INFINITY);
        // Axial overlap and the sealed finite turn count bound narrow indices.
        for period in quotient.lo.floor() as i64 - 1..=quotient.hi.floor() as i64 + 1 {
            let phase = q - I::point(period as f64) * I::point(self.spec.pitch_mm);
            for band in 0..5 {
                let p = I::new(
                    phase.lo.max(knots[band][0].lo),
                    phase.hi.min(knots[band + 1][0].hi),
                );
                if p.lo > p.hi {
                    continue;
                }
                let u = (p - knots[band][0]) / (knots[band + 1][0] - knots[band][0]);
                let u = [u.lo.max(0.), u.hi.min(1.)];
                let r = interval(self.profile_bounds(band, u)?.radius_mm);
                result.lo = result.lo.min(r.lo);
                result.hi = result.hi.max(r.hi);
            }
        }
        result.checked().map_err(|_| Error::Unresolved)
    }
}
fn angle(x: I, y: I) -> I {
    let pi = pi();
    if x.lo > 0. {
        (y / x).atan()
    } else if x.hi < 0. && y.lo >= 0. {
        (y / x).atan() + pi
    } else if x.hi < 0. && y.hi <= 0. {
        (y / x).atan() - pi
    } else if y.lo > 0. {
        pi / I::point(2.) - (x / y).atan()
    } else if y.hi < 0. {
        -pi / I::point(2.) - (x / y).atan()
    } else {
        I::new(-pi.hi, pi.hi)
    }
}
