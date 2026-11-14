use super::*;
pub(super) type Poly = Vec<Z>;
const BITS: u64 = 8192;
pub(super) struct Arithmetic<'w, 'i> {
    pub work: &'w mut Work<'i>,
}
fn trim(p: &mut Poly) {
    while p.len() > 1 && p.last().unwrap().is_zero() {
        p.pop();
    }
}
fn words(z: &Z) -> usize {
    z.bits().div_ceil(64).max(1) as usize
}
impl Arithmetic<'_, '_> {
    fn reserve(&mut self, bits: u64, cost: usize) -> Result<()> {
        self.work.step(cost)?;
        if bits > BITS {
            return Err(ReadError::Limit);
        }
        Ok(())
    }
    fn add(&mut self, a: &Z, b: &Z) -> Result<Z> {
        self.reserve(a.bits().max(b.bits()) + 1, 1 + words(a).max(words(b)))?;
        Ok(a + b)
    }
    fn sub(&mut self, a: &Z, b: &Z) -> Result<Z> {
        self.reserve(a.bits().max(b.bits()) + 1, 1 + words(a).max(words(b)))?;
        Ok(a - b)
    }
    pub(super) fn mul(&mut self, a: &Z, b: &Z) -> Result<Z> {
        self.reserve(a.bits() + b.bits(), 1 + words(a) * words(b))?;
        Ok(a * b)
    }
    pub(super) fn shift(&mut self, a: &Z, n: u32) -> Result<Z> {
        self.reserve(
            a.bits() + u64::from(n),
            1 + words(a) + (n as usize).div_ceil(64),
        )?;
        Ok(a << n)
    }
    fn division(&mut self, a: &Z, b: &Z) -> Result<(Z, Z)> {
        self.reserve(a.bits().max(b.bits()), 1 + words(a) * words(b))?;
        // Euclid calls only with a nonzero remainder; primitive division uses
        // positive content>1, and square-free division uses nonzero gcd(P,P').
        Ok((a / b, a % b))
    }
    fn content(&mut self, p: &Poly) -> Result<Z> {
        let mut gcd = Z::zero();
        for coefficient in p {
            let mut b = coefficient.abs();
            while !b.is_zero() {
                let (_, r) = self.division(&gcd, &b)?;
                gcd = b;
                b = r;
            }
        }
        Ok(gcd)
    }
    fn primitive(&mut self, mut p: Poly) -> Result<Poly> {
        trim(&mut p);
        let content = self.content(&p)?;
        if content > Z::one() {
            for coefficient in &mut p {
                *coefficient = self.division(coefficient, &content)?.0;
            }
        }
        Ok(p)
    }
    // Positive pseudo-division preserves the sign of the true rational
    // remainder even when the divisor's leading coefficient is negative.
    fn remainder(&mut self, mut p: Poly, q: &Poly) -> Result<Poly> {
        let lc = q.last().unwrap().abs();
        let sign = if q.last().unwrap().is_positive() {
            1
        } else {
            -1
        };
        while p.len() >= q.len() && !(p.len() == 1 && p[0].is_zero()) {
            let shift = p.len() - q.len();
            let lead = self.mul(p.last().unwrap(), &Z::from(sign))?;
            for c in &mut p {
                *c = self.mul(c, &lc)?;
            }
            for (i, c) in q.iter().enumerate() {
                let term = self.mul(&lead, c)?;
                p[i + shift] = self.sub(&p[i + shift], &term)?;
            }
            p = self.primitive(p)?;
        }
        Ok(p)
    }
    pub fn derivative(&mut self, p: &Poly) -> Result<Poly> {
        self.work.step(1)?;
        if p.len() == 1 {
            return Ok(vec![Z::zero()]);
        }
        p.iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| self.mul(c, &Z::from(i)))
            .collect()
    }
    pub fn gcd(&mut self, mut p: Poly, mut q: Poly) -> Result<Poly> {
        while !(q.len() == 1 && q[0].is_zero()) {
            let r = self.remainder(p, &q)?;
            p = q;
            q = r;
        }
        p = self.primitive(p)?;
        if p.last().unwrap().is_negative() {
            for c in &mut p {
                self.work.step(1 + words(c))?;
                *c = -c.clone();
            }
        }
        Ok(p)
    }
    pub fn divide(&mut self, mut p: Poly, q: &Poly) -> Result<Poly> {
        let mut result = vec![Z::zero(); p.len() - q.len() + 1];
        while p.len() >= q.len() && !(p.len() == 1 && p[0].is_zero()) {
            let index = p.len() - q.len();
            // The only caller divides primitive integerP by its primitive
            // gcd(P,P'). Gauss' lemma earns an integer quotient at every step.
            let (c, _) = self.division(p.last().unwrap(), q.last().unwrap())?;
            result[index] = c.clone();
            for (i, term) in q.iter().enumerate() {
                let t = self.mul(&c, term)?;
                p[i + index] = self.sub(&p[i + index], &t)?;
            }
            trim(&mut p);
        }
        // Exact Euclidean gcd divisibility earns the final zero remainder.
        trim(&mut result);
        Ok(result)
    }
    pub fn sturm(&mut self, p: Poly) -> Result<Vec<Poly>> {
        let derivative = self.derivative(&p)?;
        let mut result = vec![p, derivative];
        while result.last().unwrap().len() > 1 {
            let n = result.len();
            let mut r = self.remainder(result[n - 2].clone(), &result[n - 1])?;
            // InputP is P/gcd(P,P'), hence square-free. Its Sturm chain
            // reaches a nonzero constant before any zero remainder.
            for c in &mut r {
                self.work.step(1 + words(c))?;
                *c = -c.clone();
            }
            result.push(r);
        }
        Ok(result)
    }
    pub fn evaluate(&mut self, p: &Poly, x: DyadicParameter) -> Result<Z> {
        let mut value = p.last().unwrap().clone();
        for (i, c) in p.iter().rev().skip(1).enumerate() {
            let product = self.mul(&value, &Z::from(x.numerator))?;
            let scaled = self.shift(c, x.depth * (i as u32 + 1))?;
            value = self.add(&product, &scaled)?;
        }
        Ok(value)
    }

    pub fn bernstein(&mut self, p: &Poly) -> Result<Poly> {
        let n = p.len() - 1;
        let factorial = (1..=n as u64).product::<u64>();
        let choose = |n: usize, k: usize| -> u64 {
            (0..k).fold(1, |v, j| v * (n - j) as u64 / (j + 1) as u64)
        };
        let mut output = Vec::with_capacity(n + 1);
        for k in 0..=n {
            let mut coefficient = Z::zero();
            for (j, power) in p.iter().enumerate().take(k + 1) {
                let factor = Z::from(choose(k, j) * (factorial / choose(n, j)));
                let term = self.mul(power, &factor)?;
                coefficient = self.add(&coefficient, &term)?;
            }
            output.push(coefficient);
        }
        Ok(output)
    }
    pub fn split_bernstein(&mut self, mut p: Poly) -> Result<[Poly; 2]> {
        let n = p.len() - 1;
        let mut left = vec![Z::zero(); n + 1];
        let mut right = left.clone();
        // A common positive2^n scaling removes all de Casteljau denominators.
        // Keep that common scaling; no expensive primitive/GCD is needed.
        for level in 0..=n {
            left[level] = self.shift(&p[0], (n - level) as u32)?;
            right[n - level] = self.shift(&p[n - level], (n - level) as u32)?;
            for i in 0..n - level {
                p[i] = self.add(&p[i], &p[i + 1])?;
            }
        }
        Ok([left, right])
    }
    pub fn radius_polynomial(&mut self, poles: [Point; 10], radius: f64) -> Result<Poly> {
        // Twenty-one fixed finite-float decodes before any integer buffers.
        self.work.step(64)?;
        fn parts(v: f64) -> (Z, i32) {
            let bits = v.to_bits();
            let e = ((bits >> 52) & 2047) as i32;
            let mantissa = (bits & ((1u64 << 52) - 1)) | if e == 0 { 0 } else { 1u64 << 52 };
            let mut m = Z::from(mantissa);
            if bits >> 63 != 0 {
                m = -m;
            }
            (m, if e == 0 { -1074 } else { e - 1023 - 52 })
        }
        let raw = poles
            .iter()
            .flat_map(|p| [p[0], p[1]])
            .chain([radius])
            .map(parts)
            .collect::<Vec<_>>();
        // The earned repaired head circle has strictly positive radius, so
        // raw[20] is nonzero. Zero coordinate leaves need no exponent shift.
        let exponent = raw[..20]
            .iter()
            .filter(|(m, _)| !m.is_zero())
            .fold(raw[20].1, |minimum, (_, e)| minimum.min(*e));
        let mut integers = Vec::new();
        for (m, e) in raw {
            integers.push(if m.is_zero() {
                m
            } else {
                self.shift(&m, (e - exponent) as u32)?
            });
        }
        let choose = |n: usize, k: usize| -> u64 {
            (0..k).fold(1, |v, j| v * (n - j) as u64 / (j + 1) as u64)
        };
        let mut xy = [vec![Z::zero(); 10], vec![Z::zero(); 10]];
        for axis in 0..2 {
            for i in 0..10 {
                for j in 0..10 - i {
                    let factor =
                        Z::from(choose(9, i) * choose(9 - i, j)) * if j % 2 == 0 { 1 } else { -1 };
                    let term = self.mul(&integers[2 * i + axis], &factor)?;
                    xy[axis][i + j] = self.add(&xy[axis][i + j], &term)?;
                }
            }
        }
        let mut p = vec![Z::zero(); 19];
        for row in xy {
            for (i, a) in row.iter().enumerate() {
                for (j, b) in row.iter().enumerate() {
                    let t = self.mul(a, b)?;
                    p[i + j] = self.add(&p[i + j], &t)?;
                }
            }
        }
        let square = self.mul(&integers[20], &integers[20])?;
        p[0] = self.sub(&p[0], &square)?;
        self.primitive(p)
    }
}
