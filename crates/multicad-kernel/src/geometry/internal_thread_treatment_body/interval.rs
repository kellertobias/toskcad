//! Small outward enclosures for the actual native source-to-world map.
use multicad_rust_kernel::authored_internal_treatment::Enclosure;
use std::ops::{Add, Div, Mul, Neg, Sub};
#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::geometry) struct Interval {
    pub lo: f64,
    pub hi: f64,
}
impl Interval {
    pub const ZERO: Self = Self { lo: 0., hi: 0. };
    pub const ONE: Self = Self { lo: 1., hi: 1. };
    pub fn point(v: f64) -> Self {
        Self { lo: v, hi: v }
    }
    pub fn new(lo: f64, hi: f64) -> Self {
        Self { lo, hi }
    }
    pub fn from(v: Enclosure) -> Self {
        Self::new(v.lower, v.upper)
    }
    pub fn source(self) -> Enclosure {
        Enclosure {
            lower: self.lo,
            upper: self.hi,
        }
    }
    pub fn valid(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite() && self.lo <= self.hi
    }
    pub fn width(self) -> f64 {
        if self.lo == self.hi {
            0.
        } else {
            (self.hi - self.lo).next_up()
        }
    }
    pub fn magnitude(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }
    pub fn square(self) -> Self {
        if !self.valid() {
            return Self::new(f64::NEG_INFINITY, f64::INFINITY);
        }
        let lo = if self.lo <= 0. && self.hi >= 0. {
            0.
        } else {
            self.lo.abs().min(self.hi.abs()).powi(2).next_down().max(0.)
        };
        Self::new(lo, self.magnitude().powi(2).next_up())
    }
    pub fn sqrt(self) -> Self {
        if !self.valid() || self.lo < 0. {
            return Self::new(f64::NEG_INFINITY, f64::INFINITY);
        }
        Self::new(self.lo.sqrt().next_down().max(0.), self.hi.sqrt().next_up())
    }
    pub fn intersect(self, r: Self) -> Option<Self> {
        self.checked()?;
        r.checked()?;
        Self::new(self.lo.max(r.lo), self.hi.min(r.hi)).checked()
    }
    pub fn checked(self) -> Option<Self> {
        self.valid().then_some(self)
    }
}
impl Add for Interval {
    type Output = Self;
    fn add(self, r: Self) -> Self {
        if !self.valid() || !r.valid() {
            return Self::new(f64::NEG_INFINITY, f64::INFINITY);
        }
        if self == Self::ZERO {
            return r;
        }
        if r == Self::ZERO {
            return self;
        }
        Self::new((self.lo + r.lo).next_down(), (self.hi + r.hi).next_up())
    }
}
impl Neg for Interval {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.hi, -self.lo)
    }
}
impl Sub for Interval {
    type Output = Self;
    fn sub(self, r: Self) -> Self {
        self + (-r)
    }
}
impl Mul for Interval {
    type Output = Self;
    fn mul(self, r: Self) -> Self {
        if !self.valid() || !r.valid() {
            return Self::new(f64::NEG_INFINITY, f64::INFINITY);
        }
        if self == Self::ZERO || r == Self::ZERO {
            return Self::ZERO;
        }
        if self == Self::ONE {
            return r;
        }
        if r == Self::ONE {
            return self;
        }
        let values = [
            self.lo * r.lo,
            self.lo * r.hi,
            self.hi * r.lo,
            self.hi * r.hi,
        ];
        if values.into_iter().any(|v| !v.is_finite()) {
            return Self::new(f64::NEG_INFINITY, f64::INFINITY);
        }
        Self::new(
            values.into_iter().fold(f64::INFINITY, f64::min).next_down(),
            values
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        )
    }
}
impl Div for Interval {
    type Output = Self;
    fn div(self, r: Self) -> Self {
        if !self.valid() || !r.valid() {
            return Self::new(f64::NEG_INFINITY, f64::INFINITY);
        }
        if r.lo <= 0. && r.hi >= 0. {
            return Self::new(f64::NEG_INFINITY, f64::INFINITY);
        }
        self * Self::new((1. / r.hi).next_down(), (1. / r.lo).next_up())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn contains(i: Interval, x: f64) {
        assert!(i.valid() && i.lo <= x && x <= i.hi, "{i:?} misses {x}");
    }
    #[test]
    fn directed_scalar_algebra_and_exact_identities() {
        let a = Interval::new(-2., 3.);
        assert_eq!(Interval::ZERO + a, a);
        assert_eq!(a + Interval::ZERO, a);
        assert_eq!(Interval::ONE * a, a);
        assert_eq!(a * Interval::ONE, a);
        assert_eq!(Interval::ZERO * a, Interval::ZERO);
        assert_eq!(a * Interval::ZERO, Interval::ZERO);
        contains(a + Interval::new(4., 5.), 2.);
        contains(a + Interval::new(4., 5.), 8.);
        contains(a - Interval::point(1.), -3.);
        contains(a * Interval::new(-4., 5.), -10.);
        contains(a * Interval::new(-4., 5.), 15.);
        contains(a / Interval::new(2., 4.), -1.);
        contains(a / Interval::new(-4., -2.), -1.5);
        contains(a.square(), 0.);
        contains(a.square(), 9.);
        contains(Interval::new(-3., -2.).square(), 4.);
        contains(Interval::new(-3., -2.).square(), 9.);
        contains(Interval::new(4., 9.).sqrt(), 2.);
        contains(Interval::new(4., 9.).sqrt(), 3.);
        assert_eq!(Interval::point(2.).width(), 0.);
        assert!(a.width() >= 5.);
        assert_eq!(a.magnitude(), 3.);
        assert_eq!(Interval::from(a.source()), a);
        assert_eq!(
            a.intersect(Interval::new(1., 4.)),
            Some(Interval::new(1., 3.))
        );
        assert!(a.intersect(Interval::new(4., 5.)).is_none());
        // The rounded product is below its exact value by the FMA residual.
        let x = 1. + f64::EPSILON;
        let p = x * x;
        let residual = x.mul_add(x, -p);
        let result = Interval::point(x) * Interval::point(x);
        assert!(residual > 0. && result.lo <= p && result.hi > p);
    }
    #[test]
    fn invalid_overflow_and_zero_divisors_cannot_be_hidden() {
        for bad in [
            Interval::point(f64::NAN),
            Interval::point(f64::INFINITY),
            Interval::new(2., 1.),
        ] {
            assert!(bad.checked().is_none());
            assert!((bad + Interval::ZERO).checked().is_none());
            assert!((Interval::ZERO + bad).checked().is_none());
            assert!((bad * Interval::ZERO).checked().is_none());
            assert!((Interval::ONE * bad).checked().is_none());
            assert!((bad / Interval::ONE).checked().is_none());
            assert!((Interval::ONE / bad).checked().is_none());
            assert!(bad.square().checked().is_none());
            assert!(bad.sqrt().checked().is_none());
            assert!(bad.intersect(Interval::ONE).is_none());
            assert!(Interval::ONE.intersect(bad).is_none());
        }
        assert!((Interval::ONE / Interval::new(-1., 1.)).checked().is_none());
        assert!(
            (Interval::point(f64::MAX) * Interval::point(2.))
                .checked()
                .is_none()
        );
        assert!(
            (Interval::point(f64::MAX) + Interval::point(f64::MAX))
                .checked()
                .is_none()
        );
        assert!(Interval::new(-1., 4.).sqrt().checked().is_none());
        assert!(Interval::point(f64::MAX).square().checked().is_none());
    }
}
