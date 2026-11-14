use super::*;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
pub(super) struct Grid(pub(super) BigInt);
pub(super) fn charge(bits: u64, w: &mut Work<'_>) -> Result<(), Error> {
    if bits > 16384 {
        return Err(Error::Limit);
    }
    let words = (bits as usize).div_ceil(64);
    w.charge(1 + words * words)
}
pub(super) fn rounded_ratio(n: BigInt, d: &BigInt) -> BigInt {
    let negative = n.is_negative();
    let absolute = n.abs();
    let q = &absolute / d;
    let r = &absolute % d;
    let twice = &r << 1;
    let q = if twice > *d || (twice == *d && (&q % BigInt::from(2)) != BigInt::zero()) {
        q + 1
    } else {
        q
    };
    if negative { -q } else { q }
}
impl Grid {
    pub(super) fn rational(q: &N, bits: usize, w: &mut Work<'_>) -> Result<Self, Error> {
        charge(q.0.numer().bits() + bits as u64 + q.0.denom().bits(), w)?;
        Ok(Self(rounded_ratio(q.0.numer() << bits, q.0.denom())))
    }
    pub(super) fn actual(&self, origin: f64, bits: usize) -> BigRational {
        BigRational::from_float(origin).unwrap()
            + BigRational::new(self.0.clone(), BigInt::one() << bits)
    }
}
pub(super) fn upper(q: &BigRational, w: &mut Work<'_>) -> Result<f64, Error> {
    Ok(N(q.clone()).enclosure(w)?.hi)
}
pub(super) fn dyadic_text(value: &Grid, bits: usize, w: &mut Work<'_>) -> Result<String, Error> {
    // This typed numerator has denominator2^bits, so reducing factors of2 is
    // exact and needs no generic-rational domain guess or float conversion.
    let trailing = value
        .0
        .trailing_zeros()
        .unwrap_or(bits as u64)
        .min(bits as u64) as usize;
    let k = bits - trailing;
    charge(value.0.bits() + 3 * k as u64 + 1, w)?;
    let n = (value.0.abs() >> trailing) * BigInt::from(5).pow(k as u32);
    let digits = n.to_str_radix(10);
    w.charge(digits.len() + k + 1)?;
    let sign = if value.0.is_negative() { "-" } else { "" };
    Ok(if k == 0 {
        format!("{sign}{digits}.")
    } else if digits.len() > k {
        format!(
            "{sign}{}.{}",
            &digits[..digits.len() - k],
            &digits[digits.len() - k..]
        )
    } else {
        format!("{sign}0.{}{digits}", "0".repeat(k - digits.len()))
    })
}
