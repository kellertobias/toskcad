//! Outward Taylor coefficients through order12. A range-valued expansion point
//! encloses derivatives throughout that range; it is not a sampled remainder.
use super::*;
pub(super) const ORDER: usize = 12;
#[derive(Clone, Copy)]
pub(super) struct S(pub [I; ORDER + 1]);
impl S {
    pub fn constant(v: I) -> Self {
        let mut c = [I::point(0.); ORDER + 1];
        c[0] = v;
        Self(c)
    }
    pub fn variable(v: I) -> Self {
        let mut s = Self::constant(v);
        s.0[1] = I::point(1.);
        s
    }
    pub fn add(self, r: Self) -> Self {
        Self(std::array::from_fn(|n| self.0[n] + r.0[n]))
    }
    pub fn sub(self, r: Self) -> Self {
        self.add(r.scale(I::point(-1.)))
    }
    pub fn scale(self, r: I) -> Self {
        Self(self.0.map(|v| v * r))
    }
    pub fn mul(self, r: Self) -> Self {
        Self(std::array::from_fn(|n| {
            (0..=n).fold(I::point(0.), |v, k| v + self.0[k] * r.0[n - k])
        }))
    }
    pub fn square(self) -> Self {
        self.mul(self)
    }
    pub fn sqrt(self) -> Result<Self, Error> {
        if !self.0[0].valid() || self.0[0].lo <= 0. {
            return Err(Error::Unresolved);
        }
        let mut c = [I::point(0.); ORDER + 1];
        c[0] = self.0[0].sqrt();
        let divisor = I::point(2.) * c[0];
        for n in 1..=ORDER {
            c[n] = (self.0[n] - (1..n).fold(I::point(0.), |v, k| v + c[k] * c[n - k])) / divisor;
        }
        Self(c).checked()
    }
    pub fn checked(self) -> Result<Self, Error> {
        if self.0.into_iter().all(I::valid) {
            Ok(self)
        } else {
            Err(Error::Unresolved)
        }
    }
}
/// Exact affine source angle; base and slope are independent outward bounds
/// for its actual coefficients, not derivative fits from evaluated points.
pub(super) fn sin_cos(base: I, slope: I, u: I) -> (S, S) {
    let a = base + slope * u;
    let (s, c) = (a.sin(), a.cos());
    let mut factor = I::point(1.);
    let (mut sine, mut cosine) = ([I::point(0.); ORDER + 1], [I::point(0.); ORDER + 1]);
    for n in 0..=ORDER {
        if n > 0 {
            factor = factor * slope / I::point(n as f64)
        }
        sine[n] = factor * [s, c, -s, -c][n % 4];
        cosine[n] = factor * [c, -s, -c, s][n % 4];
    }
    (S(sine), S(cosine))
}
