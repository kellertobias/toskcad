//! Whole-box membership in the two actual shared-boundary phase laws.
use super::*;
impl Insert {
    pub fn classify_box(
        &self,
        bounds: [[f64; 2]; 3],
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Classification, Error> {
        use Classification::*;
        checkpoint(interrupt)?;
        let [x, y, z] = valid_box(bounds)?;
        let rho = (x.square() + y.square()).sqrt();
        let result = if z.hi < 0. || z.lo > self.full_length_mm() {
            Outside
        } else if !rho.valid() {
            UnresolvedBoundary
        } else {
            let angle = angle(x, y);
            let outer = radius(&self.outer, z, angle);
            let inner = radius(&self.inner, z, angle);
            match (outer, inner) {
                (Some(o), Some(i)) if rho.lo > o.hi || rho.hi < i.lo => Outside,
                (Some(o), Some(i))
                    if z.lo > 0.
                        && z.hi < self.full_length_mm()
                        && rho.hi < o.lo
                        && rho.lo > i.hi =>
                {
                    Inside
                }
                _ => UnresolvedBoundary,
            }
        };
        checkpoint(interrupt)?;
        Ok(result)
    }
}
fn angle(x: I, y: I) -> I {
    let pi = I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    );
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
fn radius(source: &AuthoredThread, z: I, theta: I) -> Option<I> {
    let q = z - theta / source.angular_gain();
    let pitch = source.spec().pitch_mm;
    let quotient = q / I::point(pitch);
    if !quotient.valid() {
        return None;
    }
    let bands = source.bands();
    let mut bound = I::new(f64::INFINITY, f64::NEG_INFINITY);
    if quotient.width() >= 1. {
        return Some(I::new(
            source.minimum_radius_mm(),
            source.maximum_radius_mm(),
        ));
    }
    // Axial overlap plus narrow quotient bounds its indices by sealed finite
    // source turn limits; at most four neighboring periods can contribute.
    for k in quotient.lo.floor() as i64 - 1..=quotient.hi.floor() as i64 + 1 {
        let phase = q - I::point(k as f64) * I::point(pitch);
        for band in bands {
            let p = I::new(
                phase.lo.max(band.phase_mm[0]),
                phase.hi.min(band.phase_mm[1]),
            );
            if p.lo > p.hi {
                continue;
            }
            let [a, b] = band.radius_mm;
            let u = (p - I::point(band.phase_mm[0]))
                / (I::point(band.phase_mm[1]) - I::point(band.phase_mm[0]));
            let r = I::point(a) + (I::point(b) - I::point(a)) * u;
            bound.lo = bound.lo.min(r.lo.max(a.min(b)));
            bound.hi = bound.hi.max(r.hi.min(a.max(b)));
        }
    }
    bound.valid().then_some(bound)
}
