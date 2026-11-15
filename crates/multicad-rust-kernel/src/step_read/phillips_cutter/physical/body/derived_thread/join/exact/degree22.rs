//! Fixed degree22 radius comparison of actual degree11 rim coefficients.
//! This does not enlarge the existing degree18 interface or its scratch bound.
use super::*;
/// Only the fixed12-pole radius conversion can create this finite polynomial.
pub(in super::super) struct Polynomial22(pub(super) Poly);
impl Polynomial22 {
    pub(in super::super) fn degree(&self) -> usize {
        self.0.len() - 1
    }
    pub(in super::super) fn into_polynomial(self) -> Poly {
        self.0
    }
}
impl Arithmetic<'_, '_> {
    pub(in super::super) fn radius_polynomial11(
        &mut self,
        poles: [Point; 12],
        radius: f64,
    ) -> Result<Polynomial22> {
        self.work.step(80)?;
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
        // The sealed carrier admitted finite physical poles; the SAME earned
        // head has a strictly positive literal circle radius. No fit or phase
        // reconstruction enters these exact binary64 coefficient leaves.
        let raw = poles
            .iter()
            .flat_map(|p| [p[0], p[1]])
            .chain([radius])
            .map(parts)
            .collect::<Vec<_>>();
        let exponent = raw[..24]
            .iter()
            .filter(|(m, _)| !m.is_zero())
            .fold(raw[24].1, |minimum, (_, e)| minimum.min(*e));
        let mut integers = Vec::with_capacity(25);
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
        let mut xy = [vec![Z::zero(); 12], vec![Z::zero(); 12]];
        for axis in 0..2 {
            for i in 0..12 {
                for j in 0..12 - i {
                    let factor = Z::from(choose(11, i) * choose(11 - i, j))
                        * if j % 2 == 0 { 1 } else { -1 };
                    let term = self.mul(&integers[2 * i + axis], &factor)?;
                    xy[axis][i + j] = self.add(&xy[axis][i + j], &term)?;
                }
            }
        }
        let mut p = vec![Z::zero(); 23];
        for row in xy {
            for (i, a) in row.iter().enumerate() {
                for (j, b) in row.iter().enumerate() {
                    let term = self.mul(a, b)?;
                    p[i + j] = self.add(&p[i + j], &term)?;
                }
            }
        }
        let square = self.mul(&integers[24], &integers[24])?;
        p[0] = self.sub(&p[0], &square)?;
        self.primitive(p).map(Polynomial22)
    }
    /// Common positive n! scaling of the exact power-to-Bernstein transform.
    /// n<=22 by the opaque fixed12-pole producer. Factorials and all quotients
    /// remain BigInt;22! exceeds the old degree18 u64 representation.
    pub(in super::super) fn bernstein22(&mut self, p: &Polynomial22) -> Result<Poly> {
        let n = p.0.len() - 1;
        let mut factorial = Z::one();
        for j in 1..=n {
            factorial = self.mul(&factorial, &Z::from(j))?;
        }
        let choose = |n: usize, k: usize| -> u64 {
            (0..k).fold(1, |v, j| v * (n - j) as u64 / (j + 1) as u64)
        };
        let mut output = Vec::with_capacity(n + 1);
        for k in 0..=n {
            let mut coefficient = Z::zero();
            for (j, power) in p.0.iter().enumerate().take(k + 1) {
                // Binomial(n,j) divides n! exactly, including n=0.
                let factor = self.division(&factorial, &Z::from(choose(n, j)))?.0;
                let factor = self.mul(&factor, &Z::from(choose(k, j)))?;
                let term = self.mul(power, &factor)?;
                coefficient = self.add(&coefficient, &term)?;
            }
            output.push(coefficient);
        }
        Ok(output)
    }
}
