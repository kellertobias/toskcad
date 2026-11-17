//! Bounded directed Bernstein algebra for actual curve correspondences.
//! Degree<=80 uses exact u128 binomial recurrence (including the temporary
//! multiply: C(80,40)*80<2^84), then adjacent-float enclosures of those integers.
//! No polynomial or sampled proposal is a geometry permission by itself.
use super::*;
#[cfg(test)]
mod oracle_tests;
const N: usize = 81;
pub(super) const SCRATCH: usize = 2 * 1024 * 1024;
#[derive(Clone, Copy)]
pub(super) struct Poly {
    pub(super) c: [I; N],
    pub(super) degree: usize,
}
impl Poly {
    #[cfg(test)]
    pub(super) fn constant(degree: usize, value: I, budget: &mut Budget<'_>) -> Result<Self> {
        budget.charge((degree + 1).saturating_mul(8).saturating_add(512))?;
        let mut p = Self::zero(degree)?;
        p.c[..=degree].fill(value);
        if p.valid() {
            Ok(p)
        } else {
            Err(ReadError::Unsupported)
        }
    }
    #[cfg(test)]
    pub(super) fn elevate(self, degree: usize, budget: &mut Budget<'_>) -> Result<Self> {
        if degree < self.degree {
            return Err(ReadError::Unsupported);
        }
        if degree == self.degree {
            budget.charge(512)?;
            return Ok(self);
        }
        product(
            self,
            Self::constant(degree - self.degree, I::point(1.), budget)?,
            budget,
        )
    }
    pub(super) fn zero(degree: usize) -> Result<Self> {
        if degree >= N {
            return Err(ReadError::Unsupported);
        }
        Ok(Self {
            c: [I::point(0.); N],
            degree,
        })
    }
    pub(super) fn from_controls(values: &[[I; 4]], axis: usize) -> Result<Self> {
        let mut p = Self::zero(values.len().checked_sub(1).ok_or(ReadError::Invalid)?)?;
        for (i, v) in values.iter().enumerate() {
            p.c[i] = v[axis];
        }
        Ok(p)
    }
    pub(super) fn valid(self) -> bool {
        self.c[..=self.degree].iter().all(|x| x.valid())
    }
    pub(super) fn scale(self, value: I, budget: &mut Budget<'_>) -> Result<Self> {
        budget.charge((self.degree + 1) * 16 + 512)?;
        let mut p = self;
        for c in &mut p.c[..=p.degree] {
            *c = *c * value;
        }
        if p.valid() {
            Ok(p)
        } else {
            Err(ReadError::Unsupported)
        }
    }
    pub(super) fn add(self, other: Self, subtract: bool, budget: &mut Budget<'_>) -> Result<Self> {
        if self.degree != other.degree {
            return Err(ReadError::Unsupported);
        }
        budget.charge((self.degree + 1) * 16 + 512)?;
        let mut p = self;
        for (i, c) in p.c[..=p.degree].iter_mut().enumerate() {
            *c = if subtract {
                *c - other.c[i]
            } else {
                *c + other.c[i]
            };
        }
        if p.valid() {
            Ok(p)
        } else {
            Err(ReadError::Unsupported)
        }
    }
    pub(super) fn derivative(self, budget: &mut Budget<'_>) -> Result<Self> {
        budget.charge(512)?;
        if self.degree == 0 {
            return Self::zero(0);
        }
        budget.charge(self.degree * 24 + 64)?;
        let mut p = Self::zero(self.degree - 1)?;
        for (i, c) in p.c[..=p.degree].iter_mut().enumerate() {
            *c = I::point(self.degree as f64) * (self.c[i + 1] - self.c[i]);
        }
        if p.valid() {
            Ok(p)
        } else {
            Err(ReadError::Unsupported)
        }
    }
    pub(super) fn split(self, budget: &mut Budget<'_>) -> Result<[Self; 2]> {
        budget.charge(self.degree * (self.degree + 1) / 2 * 16 + (self.degree + 1) * 32 + 512)?;
        let (mut a, mut b) = (Self::zero(self.degree)?, Self::zero(self.degree)?);
        let mut row = self.c;
        a.c[0] = row[0];
        b.c[self.degree] = row[self.degree];
        for level in 1..=self.degree {
            for i in 0..=self.degree - level {
                row[i] = (row[i] + row[i + 1]) * I::point(0.5);
            }
            a.c[level] = row[0];
            b.c[self.degree - level] = row[self.degree - level];
        }
        if a.valid() && b.valid() {
            Ok([a, b])
        } else {
            Err(ReadError::Unsupported)
        }
    }
}
fn binomial_row(n: usize, budget: &mut Budget<'_>) -> Result<[I; N]> {
    if n >= N {
        return Err(ReadError::Unsupported);
    }
    budget.charge((n + 1) * 32 + 512)?;
    let mut row = [I::point(0.); N];
    let mut value = 1_u128;
    for (i, slot) in row[..=n].iter_mut().enumerate() {
        if i > 0 {
            value = value
                .checked_mul((n - i + 1) as u128)
                .ok_or(ReadError::Unsupported)?
                / (i as u128);
        }
        let f = value as f64;
        *slot = if f as u128 == value {
            I::point(f)
        } else {
            I::new(f.next_down(), f.next_up())
        };
    }
    Ok(row)
}
pub(super) fn product(a: Poly, b: Poly, budget: &mut Budget<'_>) -> Result<Poly> {
    let degree = a
        .degree
        .checked_add(b.degree)
        .ok_or(ReadError::Unsupported)?;
    budget.charge((a.degree + 1) * (b.degree + 1) * 16 + (degree + 1) * 8 + 512)?;
    let mut p = Poly::zero(degree)?;
    let ar = binomial_row(a.degree, budget)?;
    let br = binomial_row(b.degree, budget)?;
    let cr = binomial_row(degree, budget)?;
    for i in 0..=a.degree {
        for j in 0..=b.degree {
            let weight = ar[i] * br[j] / cr[i + j];
            p.c[i + j] = p.c[i + j] + weight * a.c[i] * b.c[j];
        }
    }
    if p.valid() {
        Ok(p)
    } else {
        Err(ReadError::Unsupported)
    }
}
pub(super) fn dot_polynomial(a: [Poly; 3], b: [Poly; 3], budget: &mut Budget<'_>) -> Result<Poly> {
    product(a[0], b[0], budget)?
        .add(product(a[1], b[1], budget)?, false, budget)?
        .add(product(a[2], b[2], budget)?, false, budget)
}
pub(super) fn derivative_numerator(h: [Poly; 4], budget: &mut Budget<'_>) -> Result<[Poly; 3]> {
    let dw = h[3].derivative(budget)?;
    let mut out = [Poly::zero(0)?; 3];
    for i in 0..3 {
        out[i] = product(h[i].derivative(budget)?, h[3], budget)?.add(
            product(h[i], dw, budget)?,
            true,
            budget,
        )?;
    }
    Ok(out)
}
pub(super) fn positive(p: Poly, budget: &mut Budget<'_>) -> Result<bool> {
    budget.charge((p.degree + 1) * 8 + 64)?;
    Ok(p.c[..=p.degree].iter().all(|x| x.valid() && x.lo > 0.))
}
#[derive(Clone, Copy)]
struct Node {
    p: [Poly; 3],
    depth: usize,
}
fn vector_upper(p: [Poly; 3], budget: &mut Budget<'_>) -> Result<f64> {
    budget.charge((p[0].degree + 1) * 64 + 64)?;
    let mut value = 0_f64;
    for i in 0..=p[0].degree {
        let n = norm([p[0].c[i], p[1].c[i], p[2].c[i]]);
        if !n.valid() {
            return Err(ReadError::Unsupported);
        }
        value = value.max(n.hi);
    }
    Ok(value)
}
pub(super) fn bound(
    mut p: [Poly; 3],
    denominator: f64,
    allowed: f64,
    budget: &mut Budget<'_>,
) -> Result<f64> {
    if !denominator.is_finite() || denominator <= 0. || !allowed.is_finite() || allowed <= 0. {
        return Err(ReadError::Unsupported);
    }
    // Fixed depth<=24 has at most25 pending siblings. The caller preflights
    // the entire 2MiB numerical stack before entering this implementation.
    let mut stack: [Option<Node>; 26] = std::array::from_fn(|_| None);
    let (mut len, mut depth, mut maximum) = (0_usize, 0_usize, 0_f64);
    loop {
        budget.charge(64)?;
        let upper = (I::point(vector_upper(p, budget)?) / I::point(denominator)).hi;
        if !upper.is_finite() {
            return Err(ReadError::Unsupported);
        }
        if upper <= allowed {
            maximum = maximum.max(upper);
            if len == 0 {
                return Ok(maximum);
            }
            len -= 1;
            let node = stack[len].take().ok_or(ReadError::Invalid)?;
            p = node.p;
            depth = node.depth;
            continue;
        }
        if depth == 24 {
            return Err(ReadError::Unsupported);
        }
        let x = p[0].split(budget)?;
        let y = p[1].split(budget)?;
        let z = p[2].split(budget)?;
        if len >= stack.len() {
            return Err(ReadError::Limit);
        }
        stack[len] = Some(Node {
            p: [x[1], y[1], z[1]],
            depth: depth + 1,
        });
        len += 1;
        p = [x[0], y[0], z[0]];
        depth += 1;
    }
}
