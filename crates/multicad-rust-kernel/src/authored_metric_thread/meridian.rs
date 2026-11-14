//! Exact symbolic circle chart of this authored source's finite rounded root.
//! Enclosed controls are not independently rounded physical inputs or a STEP
//! certificate. Literal serialization must prove its own transfer afterwards.
use super::*;

#[derive(Clone, Copy, Debug)]
pub struct RationalMeridian {
    source: AuthoredMetricThread,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeridianBounds {
    /// Coordinates are [original phase q, radius], in millimetres.
    pub point_mm: [[f64; 2]; 2],
    /// Derivatives use this chart's normalized parameter t, not phase q.
    pub first_derivative_mm: [[f64; 2]; 2],
    pub second_derivative_mm: [[f64; 2]; 2],
}
impl AuthoredMetricThread {
    pub fn root_meridian(&self) -> RationalMeridian {
        RationalMeridian { source: *self }
    }
}
impl RationalMeridian {
    pub fn source(&self) -> &AuthoredMetricThread {
        &self.source
    }
    /// Homogeneous quadratic Bernstein controls [w*q,w*r,w]. The middle
    /// weight is EXACTLY 1/2 in the authored expression, not a fitted float.
    pub fn homogeneous_controls(&self) -> [[[f64; 2]; 3]; 3] {
        let p = I::point(self.source.spec.pitch_mm);
        let r = self.source.root_radius();
        let c = self.source.center_radius();
        [
            [
                p * I::point(3.) / I::point(8.),
                c - r / I::point(2.),
                I::point(1.),
            ],
            [
                p / I::point(4.),
                (c - I::point(2.) * r) / I::point(2.),
                I::point(0.5),
            ],
            [
                p * I::point(5.) / I::point(8.),
                c - r / I::point(2.),
                I::point(1.),
            ],
        ]
        .map(|row| row.map(pair))
    }
    /// Whole closed parameter interval; no sampled or fitted bound. Fixed
    /// stack-only arithmetic (<256 scalar interval operations), no allocation.
    pub fn bounds(
        &self,
        parameter: [f64; 2],
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MeridianBounds, Error> {
        checkpoint(interrupt)?;
        let t = I::new(parameter[0], parameter[1]);
        if !t.valid() || t.lo < 0. || t.hi > 1. {
            return Err(Error::Invalid);
        }
        let one = I::point(1.);
        let p = I::point(self.source.spec.pitch_mm);
        let r = self.source.root_radius();
        let c = self.source.center_radius();
        let x = I::point(2.) * t - one;
        let x2 = x.square();
        // W=1-t+t²=(3+x²)/4 >=3/4. q' has strictly positive
        // numerator p*(3-x²); this proves the entire trim is traversed once.
        let w = (I::point(3.) + x2) / I::point(4.);
        let q = p / I::point(2.) + p * x / (I::point(8.) * w);
        let radius = c + r - I::point(1.5) * r / w;
        let dq = p * (I::point(3.) - x2) / (I::point(16.) * w.square());
        let dr = I::point(1.5) * r * x / w.square();
        let ddq = p * x * (x2 - I::point(9.)) / (I::point(16.) * w * w * w);
        let ddr = I::point(9.) * r * (one - x2) / (I::point(4.) * w * w * w);
        // Constructor bounds p,R,C by finite engineering dimensions; here
        // |x|<=1 and W>=3/4. No generic overflow/domain handoff is needed.
        // With u=1/2+x/(2W), X=sqrt(3)*(u-1/2) and
        // H=(3/2-W)/W>0. Thus q=3p/8+(p/4)u and radius=C-R*H
        // are the SAME original profile expressions, including both tangents.
        let result = MeridianBounds {
            point_mm: [q, radius].map(pair),
            first_derivative_mm: [dq, dr].map(pair),
            second_derivative_mm: [ddq, ddr].map(pair),
        };
        checkpoint(interrupt)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
