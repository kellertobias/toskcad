//! Strict whole-box membership, preserving the two distinct source phase laws.
use super::*;
impl Hanger {
    pub fn classify_box(
        &self,
        bounds: [[f64; 2]; 3],
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Classification, Error> {
        use Classification::*;
        checkpoint(interrupt)?;
        let p = bounds.map(|[lo, hi]| I::new(lo, hi));
        if p.into_iter().any(|v| !v.valid()) {
            return Err(Error::Invalid);
        }
        let [x, y, z] = p;
        let rho = (x.square() + y.square()).sqrt();
        let plain = -self.spec.plain_length_mm;
        let result = if z.hi < self.metric_bottom[0]
            || z.lo > self.wood().full_length_mm()
            || rho.lo > self.wood().maximum_radius_mm()
        {
            Outside
        } else if !rho.valid() {
            UnresolvedBoundary
        } else if z.lo > self.metric_bottom[1]
            && z.hi < self.wood().full_length_mm()
            && rho.hi < self.core_radius_mm
        {
            Inside
        } else if z.lo > 0. && z.hi < self.wood().full_length_mm() {
            match wood_radius(self.wood(), z, angle(x, y)) {
                Some(r) if rho.hi < r.lo => Inside,
                Some(r) if rho.lo > r.hi => Outside,
                _ => UnresolvedBoundary,
            }
        } else if z.lo > plain && z.hi < 0. {
            if rho.hi < self.spec.plain_radius_mm {
                Inside
            } else if rho.lo > self.spec.plain_radius_mm {
                Outside
            } else {
                UnresolvedBoundary
            }
        } else if z.lo > self.metric_bottom[1] && z.hi < plain {
            checkpoint(interrupt)?;
            // The original local metric law is evaluated at z+P+Lm. Each
            // operation encloses the symbolic sum; no rounded origin is used.
            let local =
                z + I::point(self.spec.plain_length_mm) + I::point(self.metric().spec().length_mm);
            let translated = [bounds[0], bounds[1], pair(local)];
            match self
                .metric()
                .classify_box(translated)
                .map_err(|_| Error::Invalid)?
            {
                crate::helical_thread::BoxClassification::Inside => Inside,
                crate::helical_thread::BoxClassification::Outside => Outside,
                crate::helical_thread::BoxClassification::UnresolvedBoundary => UnresolvedBoundary,
            }
        } else {
            UnresolvedBoundary
        };
        checkpoint(interrupt)?;
        Ok(result)
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
fn wood_radius(source: &AuthoredThread, z: I, angle: I) -> Option<I> {
    let q = z - angle / source.angular_gain();
    let pitch = source.spec().pitch_mm;
    let quotient = q / I::point(pitch);
    if !quotient.valid() {
        return None;
    }
    if quotient.width() >= 1. {
        return Some(I::new(
            source.minimum_radius_mm(),
            source.maximum_radius_mm(),
        ));
    }
    let mut result = I::new(f64::INFINITY, f64::NEG_INFINITY);
    // Only strict finite wood-domain boxes reach this point. The sealed source
    // bounds L/p+2N<=128; narrow angle boxes require at most4 nearby periods.
    for k in quotient.lo.floor() as i64 - 1..=quotient.hi.floor() as i64 + 1 {
        let phase = q - I::point(k as f64) * I::point(pitch);
        for b in source.bands() {
            let p = I::new(phase.lo.max(b.phase_mm[0]), phase.hi.min(b.phase_mm[1]));
            if p.lo > p.hi {
                continue;
            }
            let [a, b_radius] = b.radius_mm;
            let u =
                (p - I::point(b.phase_mm[0])) / (I::point(b.phase_mm[1]) - I::point(b.phase_mm[0]));
            let r = I::point(a) + (I::point(b_radius) - I::point(a)) * u;
            result.lo = result.lo.min(r.lo.max(a.min(b_radius)));
            result.hi = result.hi.max(r.hi.min(a.max(b_radius)));
        }
    }
    result.valid().then_some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn overflowing_query_phase_never_chooses_a_nominal_band() {
        let spec = super::super::tests::fixture(false, false);
        assert!(wood_radius(&spec.wood, I::point(f64::MAX), I::point(f64::MIN)).is_none());
    }
}
