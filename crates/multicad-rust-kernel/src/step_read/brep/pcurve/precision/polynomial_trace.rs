//! Small directed Bernstein algebra for original analytic trim composition.
//! No polynomial operation grants a curve, parameter map or geometry receipt.
use super::*;
use crate::step_read::brep::check;
const MAX: usize = 16;
#[derive(Clone, Copy)]
pub(super) struct Polynomial {
    degree: usize,
    controls: [I; MAX + 1],
}
fn choose(n: usize, k: usize) -> f64 {
    // Every intermediate is an integer <= C(16,8), exactly representable.
    (0..k).fold(1usize, |v, i| v * (n - i) / (i + 1)) as f64
}
impl Polynomial {
    pub(super) fn new(values: &[I], b: &mut Budget<'_>) -> Result<Self> {
        b.charge(128)?;
        check(
            !values.is_empty() && values.len() <= MAX + 1,
            ReadError::Unsupported,
        )?;
        check(values.iter().all(|v| v.valid()), ReadError::Unsupported)?;
        let mut controls = [I::point(0.); MAX + 1];
        controls[..values.len()].copy_from_slice(values);
        Ok(Self {
            degree: values.len() - 1,
            controls,
        })
    }
    pub(super) fn constant(value: I, b: &mut Budget<'_>) -> Result<Self> {
        Self::new(&[value], b)
    }
    pub(super) fn degree(&self) -> usize {
        self.degree
    }
    pub(super) fn controls(&self) -> &[I] {
        &self.controls[..=self.degree]
    }
    pub(super) fn hull(&self) -> I {
        I::new(
            self.controls()
                .iter()
                .fold(f64::INFINITY, |v, x| v.min(x.lo)),
            self.controls()
                .iter()
                .fold(f64::NEG_INFINITY, |v, x| v.max(x.hi)),
        )
    }
    pub(super) fn elevate(&self, n: usize, b: &mut Budget<'_>) -> Result<Self> {
        check(n >= self.degree && n <= MAX, ReadError::Unsupported)?;
        if n == self.degree {
            b.charge(256)?;
            return Ok(*self);
        }
        if self.degree == 0 {
            b.charge((n + 1) * 64 + 256)?;
            let out = [self.controls[0]; MAX + 1];
            return Self::new(&out[..=n], b);
        }
        b.charge(
            (n + 1)
                .saturating_mul(self.degree + 1)
                .saturating_mul(256)
                .saturating_add(4096),
        )?;
        let mut out = [I::point(0.); MAX + 1];
        for (i, value) in out.iter_mut().enumerate().take(n + 1) {
            for j in i.saturating_sub(n - self.degree)..=i.min(self.degree) {
                let weight = I::point(choose(self.degree, j))
                    * I::point(choose(n - self.degree, i - j))
                    / I::point(choose(n, i));
                *value = *value + weight * self.controls[j];
            }
        }
        Self::new(&out[..=n], b)
    }
    pub(super) fn add(&self, other: &Self, b: &mut Budget<'_>) -> Result<Self> {
        let n = self.degree.max(other.degree);
        let a = self.elevate(n, b)?;
        let c = other.elevate(n, b)?;
        b.charge((n + 1) * 64)?;
        let out: [I; MAX + 1] = std::array::from_fn(|i| a.controls[i] + c.controls[i]);
        Self::new(&out[..=n], b)
    }
    pub(super) fn scale(&self, value: I, b: &mut Budget<'_>) -> Result<Self> {
        b.charge((self.degree + 1) * 64)?;
        let out = self.controls.map(|v| v * value);
        Self::new(&out[..=self.degree], b)
    }
    pub(super) fn multiply(&self, other: &Self, b: &mut Budget<'_>) -> Result<Self> {
        let n = self
            .degree
            .checked_add(other.degree)
            .ok_or(ReadError::Limit)?;
        check(n <= MAX, ReadError::Unsupported)?;
        b.charge(
            (self.degree + 1)
                .saturating_mul(other.degree + 1)
                .saturating_mul(256)
                .saturating_add(4096),
        )?;
        let mut out = [I::point(0.); MAX + 1];
        for (i, a) in self.controls().iter().enumerate() {
            for (j, c) in other.controls().iter().enumerate() {
                let weight = I::point(choose(self.degree, i)) * I::point(choose(other.degree, j))
                    / I::point(choose(n, i + j));
                out[i + j] = out[i + j] + weight * (*a) * (*c);
            }
        }
        Self::new(&out[..=n], b)
    }
    pub(super) fn derivative(&self, b: &mut Budget<'_>) -> Result<Self> {
        b.charge((self.degree + 1) * 64)?;
        if self.degree == 0 {
            return Self::constant(I::point(0.), b);
        }
        let mut out = [I::point(0.); MAX + 1];
        for (i, p) in self.controls().windows(2).enumerate() {
            out[i] = I::point(self.degree as f64) * (p[1] - p[0]);
        }
        Self::new(&out[..self.degree], b)
    }
    pub(super) fn split(&self, t: I, b: &mut Budget<'_>) -> Result<[Self; 2]> {
        check(
            t.valid() && t.lo >= 0. && t.hi <= 1.,
            ReadError::Unsupported,
        )?;
        b.charge(
            (self.degree + 1)
                .saturating_mul(self.degree + 1)
                .saturating_mul(128)
                .saturating_add(4096),
        )?;
        let mut net = self.controls;
        let mut left = [I::point(0.); MAX + 1];
        let mut right = left;
        left[0] = net[0];
        right[self.degree] = net[self.degree];
        for depth in 1..=self.degree {
            for i in 0..=self.degree - depth {
                net[i] = (I::point(1.) - t) * net[i] + t * net[i + 1];
            }
            left[depth] = net[0];
            right[self.degree - depth] = net[self.degree - depth];
        }
        Ok([
            Self::new(&left[..=self.degree], b)?,
            Self::new(&right[..=self.degree], b)?,
        ])
    }
}

#[cfg(test)]
mod tests;
