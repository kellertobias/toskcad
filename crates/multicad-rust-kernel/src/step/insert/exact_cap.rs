//! Rare-span sign proof from actual decimal coefficients, not float-cache boxes.
//! All46 retained rationals and bounded bigint temporaries fit this conservative
//! scratch charge (16Kbit operand/cross-product preflight inherited from N).
use super::domains::N;
use super::*;
use num_bigint::BigInt;
use num_rational::BigRational;
use std::cmp::Ordering;
pub(crate) const SCRATCH: usize = 4 * 1024 * 1024;
fn decimal(text: &str, w: &mut Work<'_>) -> Result<N, Error> {
    w.charge(text.len() + 1)?;
    if text.len() > 512 {
        return Err(Error::Limit);
    }
    let (mantissa, exponent) = text.split_once('E').unwrap_or((text, "0"));
    let fraction = mantissa
        .split_once('.')
        .map_or(0, |(_, fraction)| fraction.len());
    let digits = mantissa.replace('.', "");
    let n = BigInt::parse_bytes(digits.as_bytes(), 10).ok_or(Error::Invalid)?;
    let exponent = exponent
        .parse::<i32>()
        .map_err(|_| Error::Limit)?
        .checked_sub(fraction as i32)
        .ok_or(Error::Limit)?;
    if exponent.unsigned_abs() > 4096 {
        return Err(Error::Limit);
    }
    let bits = (digits.len() + exponent.unsigned_abs() as usize).saturating_mul(4);
    if bits > 16384 {
        return Err(Error::Limit);
    }
    let words = bits.div_ceil(64);
    w.charge(1 + words * words)?;
    let ten = BigInt::from(10).pow(exponent.unsigned_abs());
    Ok(N(if exponent < 0 {
        BigRational::new(n, ten)
    } else {
        BigRational::from_integer(n * ten)
    }))
}
fn choose(n: usize, k: usize) -> i64 {
    let mut value = 1;
    for i in 0..k {
        value = value * (n - i) as i64 / (i + 1) as i64;
    }
    value
}
/// Exact degree21 Bernstein coefficients of sign*(X Y'−Y X'). A positive
/// coefficient minimum proves the WHOLE finite polynomial span strictly
/// angular-monotone. Failure remains a refusal, never a guessed direction.
pub(crate) fn angular_lower(
    tokens: [[&str; 2]; 12],
    reverse: bool,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<f64, Error> {
    memory.add(SCRATCH)?;
    let mut xy: [[_; 2]; 12] = std::array::from_fn(|_| std::array::from_fn(|_| N::integer(0)));
    for i in 0..12 {
        for axis in 0..2 {
            xy[i][axis] = decimal(tokens[i][axis], w)?;
        }
    }
    let mut coefficients: [N; 22] = std::array::from_fn(|_| N::integer(0));
    for i in 0..12 {
        for j in 0..11 {
            w.charge(64)?;
            let dy = xy[j + 1][1].sub(&xy[j][1], w)?;
            let dx = xy[j + 1][0].sub(&xy[j][0], w)?;
            let cross = xy[i][0].mul(&dy, w)?.sub(&xy[i][1].mul(&dx, w)?, w)?;
            let weight = N::integer(if reverse { -11 } else { 11 } * choose(11, i) * choose(10, j));
            let term = cross
                .mul(&weight, w)?
                .div(&N::integer(choose(21, i + j)), w)?;
            coefficients[i + j] = coefficients[i + j].add(&term, w)?;
        }
    }
    let zero = N::integer(0);
    let mut minimum = &coefficients[0];
    for c in &coefficients {
        if c.cmp(&zero, w)? != Ordering::Greater {
            return Err(Error::Unresolved);
        }
        if c.cmp(minimum, w)? == Ordering::Less {
            minimum = c;
        }
    }
    let lower = minimum.enclosure(w)?.lo;
    if lower <= 0. {
        return Err(Error::Unresolved);
    }
    memory.used -= SCRATCH;
    Ok(lower)
}

#[cfg(test)]
mod tests;
