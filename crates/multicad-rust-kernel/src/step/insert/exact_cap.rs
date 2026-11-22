//! Rare-span sign proof from actual decimal coefficients, not float-cache boxes.
//! Local tensor/arc coordinates and bounded bigint temporaries fit the4MiB
//! scratch charge (16Kbit operation preflight inherited from N). The retained
//! immutable-point memo is separately charged before storing fresh clones.
use super::domains::N;
use super::*;
use num_bigint::BigInt;
use num_rational::BigRational;
use std::cmp::Ordering;
mod cache;
mod enclosure;
pub(crate) use cache::{DecimalMemo, LiteralPoint};
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
/// Original finite decimal ordering when binary64 caches cannot distinguish
/// actual coefficients. This is a file-value proof, never a tolerance test.
pub(crate) fn compare_decimal(
    a: &str,
    b: &str,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Ordering, Error> {
    memory.add(SCRATCH)?;
    let a = decimal(a, w)?;
    let b = decimal(b, w)?;
    let result = a.cmp(&b, w)?;
    memory.used -= SCRATCH;
    Ok(result)
}
/// Exact ORIGINAL decimal value versus a source-owned dyadic. The selected
/// body graph and STEP mapping remain separate proof obligations.
pub(crate) fn equals_dyadic(
    text: &str,
    sign: i8,
    numerator: &[u8],
    fractional_bits: usize,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<bool, Error> {
    memory.add(SCRATCH)?;
    w.charge(numerator.len() + 1)?;
    if numerator.len() > 2048 || fractional_bits > 16383 {
        return Err(Error::Limit);
    }
    let sign = if sign < 0 {
        num_bigint::Sign::Minus
    } else if sign > 0 {
        num_bigint::Sign::Plus
    } else {
        num_bigint::Sign::NoSign
    };
    let intended = BigInt::from_bytes_le(sign, numerator);
    // The canonical candidate encodes INTEGER real tokens. Avoid building and
    // reducing a denominator when the actual grammar already gives one.
    let equal = if fractional_bits == 0 && text.ends_with('.') {
        w.charge(text.len() + 1)?;
        if text.len() > 512 {
            return Err(Error::Limit);
        }
        let actual =
            BigInt::parse_bytes(&text.as_bytes()[..text.len() - 1], 10).ok_or(Error::Invalid)?;
        w.charge(1 + (actual.bits().max(intended.bits()) as usize).div_ceil(32))?;
        actual == intended
    } else {
        let actual = decimal(text, w)?;
        let bits = intended.bits() + fractional_bits as u64 + 1;
        let words = (bits as usize).div_ceil(64);
        w.charge(1 + words * words)?;
        let expected = N(BigRational::new(
            intended,
            BigInt::from(1) << fractional_bits,
        ));
        actual.cmp(&expected, w)? == Ordering::Equal
    };
    memory.used -= SCRATCH;
    Ok(equal)
}
/// Differential coefficients of the actual encoded tensor. Subtraction and
/// owned-corner factor removal precede conversion to outward binary64 bounds.
pub(crate) fn tensor_differentials<'a>(
    tokens: [[LiteralPoint<'a>; 12]; 2],
    factor: bool,
    memo: &mut DecimalMemo<'a>,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<[[projection::polynomial::P; 3]; 3], Error> {
    use projection::polynomial::P;
    memory.add(SCRATCH)?;
    let mut rows: [[[N; 12]; 3]; 2] =
        std::array::from_fn(|_| std::array::from_fn(|_| std::array::from_fn(|_| N::integer(0))));
    for r in 0..2 {
        for i in 0..12 {
            for (axis, value) in memo.point(tokens[r][i], memory, w)?.into_iter().enumerate() {
                rows[r][axis][i] = value;
            }
        }
    }
    let zero = P::controls(&[I::point(0.)]);
    let mut out = [[zero; 3]; 3];
    for axis in 0..3 {
        let mut d: [N; 12] = std::array::from_fn(|_| N::integer(0));
        for (i, x) in d.iter_mut().enumerate() {
            *x = rows[1][axis][i].sub(&rows[0][axis][i], w)?;
        }
        let degree = if factor { 10 } else { 11 };
        if factor {
            if d[0].cmp(&N::integer(0), w)? != Ordering::Equal {
                return Err(Error::Invalid);
            }
            for i in 0..11 {
                d[i] = d[i + 1]
                    .mul(&N::integer(11), w)?
                    .div(&N::integer((i + 1) as i64), w)?;
            }
        }
        let mut reduced = [I::point(0.); 12];
        for i in 0..=degree {
            reduced[i] = enclosure::rational(&d[i], w)?;
        }
        out[0][axis] = P::controls(&reduced[..=degree]);
        let mut da = [I::point(0.); 11];
        for (i, x) in da.iter_mut().enumerate() {
            *x = enclosure::rational(
                &rows[0][axis][i + 1]
                    .sub(&rows[0][axis][i], w)?
                    .mul(&N::integer(11), w)?,
                w,
            )?;
        }
        out[1][axis] = P::controls(&da);
        let mut dd = [I::point(0.); 11];
        for i in 0..degree {
            dd[i] = enclosure::rational(
                &d[i + 1].sub(&d[i], w)?.mul(&N::integer(degree as i64), w)?,
                w,
            )?;
        }
        out[2][axis] = P::controls(&dd[..degree]);
    }
    memory.used -= SCRATCH;
    Ok(out)
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
