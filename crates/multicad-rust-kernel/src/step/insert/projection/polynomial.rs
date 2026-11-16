//! Bounded Bernstein algebra for this family's degree<=32 literal tensors.
use super::*;
#[derive(Clone, Copy)]
pub(crate) struct P {
    pub degree: usize,
    pub c: [I; 33],
}
const fn choices() -> [[u64; 33]; 33] {
    let mut c = [[0; 33]; 33];
    let mut n = 0;
    while n < 33 {
        c[n][0] = 1;
        c[n][n] = 1;
        let mut j = 1;
        while j < n {
            c[n][j] = c[n - 1][j - 1] + c[n - 1][j];
            j += 1;
        }
        n += 1;
    }
    c
}
const CHOOSE: [[u64; 33]; 33] = choices();
impl P {
    pub fn controls(values: &[I]) -> Self {
        let mut c = [I::point(0.); 33];
        c[..values.len()].copy_from_slice(values);
        Self {
            degree: values.len() - 1,
            c,
        }
    }
    pub fn add(self, b: Self, w: &mut Work<'_>) -> Result<Self, Error> {
        w.charge(self.degree + 1)?;
        let mut c = self;
        for i in 0..=self.degree {
            c.c[i] = self.c[i] + b.c[i];
        }
        Ok(c)
    }
    pub fn sub(self, b: Self, w: &mut Work<'_>) -> Result<Self, Error> {
        w.charge(self.degree + 1)?;
        let mut c = self;
        for i in 0..=self.degree {
            c.c[i] = self.c[i] - b.c[i];
        }
        Ok(c)
    }
    pub fn scale(self, f: f64, w: &mut Work<'_>) -> Result<Self, Error> {
        w.charge(self.degree + 1)?;
        let mut c = self;
        for p in &mut c.c[..=self.degree] {
            *p = *p * I::point(f);
        }
        Ok(c)
    }
    pub fn product(self, b: Self, w: &mut Work<'_>) -> Result<Self, Error> {
        let n = self.degree;
        let m = b.degree;
        w.charge((n + 1) * (m + 1))?;
        let mut out = Self {
            degree: n + m,
            c: [I::point(0.); 33],
        };
        for k in 0..=n + m {
            for i in k.saturating_sub(m)..=n.min(k) {
                let weight = I::point(CHOOSE[n][i] as f64) * I::point(CHOOSE[m][k - i] as f64)
                    / I::point(CHOOSE[n + m][k] as f64);
                out.c[k] = out.c[k] + self.c[i] * b.c[k - i] * weight;
            }
        }
        Ok(out)
    }
    pub fn derivative(self, w: &mut Work<'_>) -> Result<Self, Error> {
        w.charge(self.degree)?;
        let mut out = self;
        out.degree -= 1;
        for i in 0..=out.degree {
            out.c[i] = I::point(self.degree as f64) * (self.c[i + 1] - self.c[i]);
        }
        Ok(out)
    }
    pub fn times_t(self, w: &mut Work<'_>) -> Result<Self, Error> {
        w.charge(self.degree + 2)?;
        let mut out = Self {
            degree: self.degree + 1,
            c: [I::point(0.); 33],
        };
        for i in 1..=out.degree {
            out.c[i] = self.c[i - 1] * I::point(i as f64) / I::point(out.degree as f64);
        }
        Ok(out)
    }
    /// Called ONLY when the constant coefficient is the SAME owned point
    /// difference, hence exactly zero. No interval-overlap factor removal.
    pub fn remove_t(self, w: &mut Work<'_>) -> Result<Self, Error> {
        w.charge(self.degree)?;
        let mut out = self;
        out.degree -= 1;
        for i in 0..=out.degree {
            out.c[i] = self.c[i + 1] * I::point(self.degree as f64) / I::point((i + 1) as f64);
        }
        Ok(out)
    }
    pub fn split(self, w: &mut Work<'_>) -> Result<[Self; 2], Error> {
        w.charge((self.degree + 1) * (self.degree + 1))?;
        let mut q = self;
        let mut a = self;
        let mut b = self;
        for level in 1..=self.degree {
            for i in 0..=self.degree - level {
                q.c[i] = (q.c[i] + q.c[i + 1]) / I::point(2.);
            }
            a.c[level] = q.c[0];
            b.c[self.degree - level] = q.c[self.degree - level];
        }
        Ok([a, b])
    }
    pub fn hull(self) -> I {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for p in &self.c[..=self.degree] {
            if !p.valid() {
                return I::new(f64::NAN, f64::NAN);
            }
            lo = lo.min(p.lo);
            hi = hi.max(p.hi);
        }
        I::new(lo, hi)
    }
}
pub(super) fn dot_xy(a: [P; 3], b: [P; 3], w: &mut Work<'_>) -> Result<P, Error> {
    a[0].product(b[0], w)?.add(a[1].product(b[1], w)?, w)
}
pub(super) fn derivative(a: [P; 3], w: &mut Work<'_>) -> Result<[P; 3], Error> {
    Ok([
        a[0].derivative(w)?,
        a[1].derivative(w)?,
        a[2].derivative(w)?,
    ])
}
#[cfg(test)]
pub(super) fn test_choices() -> [[u64; 33]; 33] {
    choices()
}
