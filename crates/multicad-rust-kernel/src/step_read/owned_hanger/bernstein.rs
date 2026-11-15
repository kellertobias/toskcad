//! Whole-polynomial interval arithmetic for decoded curve/surface composition.
use super::*;
pub(super) type Polynomial = Vec<I>;
const MAX: usize = 160;
pub(super) struct Algebra {
    choose: Vec<Vec<I>>,
}
fn p(v: f64) -> I {
    I::point(v)
}
impl Algebra {
    pub(super) fn new(b: &mut Budget<'_>) -> Result<Self, Error> {
        b.charge((MAX + 1) * (MAX + 1) * 8)?;
        let mut choose = Vec::with_capacity(MAX + 1);
        choose.push(vec![p(1.)]);
        for n in 1..=MAX {
            let mut row = vec![p(1.); n + 1];
            for k in 1..n {
                row[k] = choose[n - 1][k - 1] + choose[n - 1][k];
            }
            choose.push(row);
        }
        Ok(Self { choose })
    }
    pub(super) fn multiply(
        &self,
        a: &[I],
        b: &[I],
        budget: &mut Budget<'_>,
    ) -> Result<Polynomial, Error> {
        let n = a.len().checked_sub(1).ok_or(Error::Invalid)?;
        let m = b.len().checked_sub(1).ok_or(Error::Invalid)?;
        let degree = n.checked_add(m).ok_or(Error::Limit)?;
        if degree > MAX {
            return Err(Error::Unsupported);
        }
        budget.charge(
            a.len()
                .checked_mul(b.len())
                .and_then(|v| v.checked_mul(80))
                .ok_or(Error::Limit)?,
        )?;
        let mut out = vec![p(0.); degree + 1];
        for (i, &a) in a.iter().enumerate() {
            for (j, &b) in b.iter().enumerate() {
                out[i + j] = out[i + j]
                    + a * b * self.choose[n][i] * self.choose[m][j] / self.choose[degree][i + j];
            }
        }
        if out.iter().any(|v| !v.valid()) {
            return Err(Error::Unsupported);
        }
        Ok(out)
    }
    pub(super) fn power(&self, a: &[I], n: usize, b: &mut Budget<'_>) -> Result<Polynomial, Error> {
        let mut out = vec![p(1.)];
        for _ in 0..n {
            out = self.multiply(&out, a, b)?;
        }
        Ok(out)
    }
    pub(super) fn basis(
        &self,
        x: &[I],
        complement: &[I],
        degree: usize,
        index: usize,
        b: &mut Budget<'_>,
    ) -> Result<Polynomial, Error> {
        if index > degree || degree > MAX {
            return Err(Error::Invalid);
        }
        let a = self.power(x, index, b)?;
        let c = self.power(complement, degree - index, b)?;
        let mut out = self.multiply(&a, &c, b)?;
        scale(&mut out, self.choose[degree][index], b)?;
        Ok(out)
    }
}
pub(super) fn add(a: &mut [I], b: &[I], budget: &mut Budget<'_>) -> Result<(), Error> {
    if a.len() != b.len() {
        return Err(Error::Invalid);
    }
    budget.charge(a.len() * 8)?;
    for (a, b) in a.iter_mut().zip(b) {
        *a = *a + *b;
    }
    Ok(())
}
pub(super) fn subtract(a: &[I], b: &[I], budget: &mut Budget<'_>) -> Result<Polynomial, Error> {
    if a.len() != b.len() {
        return Err(Error::Invalid);
    }
    budget.charge(a.len() * 8)?;
    Ok(a.iter().zip(b).map(|(&a, &b)| a - b).collect())
}
pub(super) fn scale(a: &mut [I], factor: I, b: &mut Budget<'_>) -> Result<(), Error> {
    b.charge(a.len() * 16)?;
    for a in a {
        *a = *a * factor;
    }
    Ok(())
}
pub(super) fn derivative(a: &[I], b: &mut Budget<'_>) -> Result<Polynomial, Error> {
    if a.len() < 2 {
        return Err(Error::Invalid);
    }
    b.charge(a.len() * 24)?;
    Ok(a.windows(2)
        .map(|w| (w[1] - w[0]) * p((a.len() - 1) as f64))
        .collect())
}
pub(super) fn split(a: &[I], t: I, b: &mut Budget<'_>) -> Result<[Polynomial; 2], Error> {
    if a.is_empty() || a.len() > MAX + 1 || !t.valid() {
        return Err(Error::Invalid);
    }
    b.charge(a.len() * a.len() * 40)?;
    let n = a.len() - 1;
    let mut work = a.to_vec();
    let mut left = vec![p(0.); a.len()];
    let mut right = left.clone();
    for level in 0..=n {
        left[level] = work[0];
        right[n - level] = work[n - level];
        for i in 0..n - level {
            work[i] = work[i] * (p(1.) - t) + work[i + 1] * t;
        }
    }
    Ok([left, right])
}
pub(super) fn restrict(a: &[I], lo: I, hi: I, b: &mut Budget<'_>) -> Result<Polynomial, Error> {
    if lo.lo < 0. || hi.hi > 1. || lo.hi >= hi.lo {
        return Err(Error::Unsupported);
    }
    let [left, _] = split(a, hi, b)?;
    let [_, out] = split(&left, lo / hi, b)?;
    Ok(out)
}
pub(super) fn evaluate(a: &[I], t: I, b: &mut Budget<'_>) -> Result<I, Error> {
    Ok(split(a, t, b)?[0][a.len() - 1])
}
pub(super) fn magnitude(a: &[I]) -> f64 {
    if a.iter().any(|v| !v.valid()) {
        return f64::INFINITY;
    }
    a.iter().map(|v| v.magnitude()).fold(0., f64::max)
}
pub(super) fn lower(a: &[I]) -> f64 {
    if a.iter().any(|v| !v.valid()) {
        return f64::NEG_INFINITY;
    }
    a.iter().map(|v| v.lo).fold(f64::INFINITY, f64::min)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bernstein_product_derivative_and_restriction_enclose_known_polynomials() {
        let mut b = Budget {
            work: 0,
            remaining: 16_000_000,
            bytes: 0,
            limit: 64 * 1024 * 1024,
            interrupt: None,
        };
        let a = Algebra::new(&mut b).unwrap();
        let x = vec![p(0.), p(1.)];
        let complement = vec![p(1.), p(0.)];
        let product = a.multiply(&x, &complement, &mut b).unwrap();
        for (actual, expected) in product.iter().zip([0., 0.5, 0.]) {
            assert!(actual.lo <= expected && actual.hi >= expected)
        }
        let dx = derivative(&x, &mut b).unwrap();
        assert!(dx[0].lo <= 1. && dx[0].hi >= 1.);
        let restricted = restrict(&product, p(0.25), p(0.75), &mut b).unwrap();
        for t in [0., 0.2, 0.5, 0.9, 1.] {
            let actual = evaluate(&restricted, p(t), &mut b).unwrap();
            let original = 0.25 + 0.5 * t;
            let exact = original * (1. - original);
            assert!(actual.lo <= exact && actual.hi >= exact);
        }
        let basis = a.basis(&x, &complement, 3, 1, &mut b).unwrap();
        for t in [0., 0.25, 0.5, 1.] {
            let actual = evaluate(&basis, p(t), &mut b).unwrap();
            let exact = 3. * t * (1. - t) * (1. - t);
            assert!(actual.lo <= exact && actual.hi >= exact)
        }
        assert!(magnitude(&product) >= 0.5);
        assert!(lower(&product) <= 0.);
    }
}
