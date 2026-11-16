//! Directed bounds from exact rational integer prefixes. The arithmetic result
//! is never replaced by a binary64 proposal or an authored geometric recipe.
use super::*;
use num_bigint::{BigUint, Sign};
use num_traits::ToPrimitive;
fn prefix(n: &BigUint, w: &mut Work<'_>) -> Result<(I, i32), Error> {
    let bits = n.bits() as usize;
    let shift = bits.saturating_sub(53);
    w.charge(4 + bits.div_ceil(32))?;
    let head = (n >> shift)
        .to_u64()
        .expect("retained prefix has at most53 bits") as f64;
    Ok((
        I::new(head, if shift == 0 { head } else { head + 1. }),
        shift as i32,
    ))
}
pub(in crate::step::insert::exact_cap) fn rational(
    value: &N,
    w: &mut Work<'_>,
) -> Result<I, Error> {
    if value.0.numer().sign() == Sign::NoSign {
        w.charge(1)?;
        return Ok(I::point(0.));
    }
    let (numerator, ns) = prefix(value.0.numer().magnitude(), w)?;
    // Every N in this module comes from canonical BigRational::new or its
    // canonical arithmetic; the denominator is strictly positive.
    let (denominator, ds) = prefix(value.0.denom().magnitude(), w)?;
    w.charge(1)?;
    let mut out = numerator / denominator;
    let mut exponent = ns - ds;
    while exponent != 0 {
        w.charge(1)?;
        let step = exponent.clamp(-512, 512);
        // Every chunk lies in the normal binary64 exponent domain; construct
        // its exact power of two directly, without a libm precision premise.
        let scale = f64::from_bits(((step + 1023) as u64) << 52);
        out = out * I::point(scale);
        exponent -= step;
        if !out.valid() {
            return Err(Error::Unresolved);
        }
    }
    if value.0.numer().sign() == Sign::Minus {
        out = I::new(-out.hi, -out.lo);
    }
    Ok(out)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn directed_prefix_ratio_contains_independent_exact_rationals_across_numeric_ranges() {
        let mut w = Work::new(16_000_000, None).unwrap();
        for text in [
            "0.",
            "1.",
            "-1.",
            "4.0000000000000001",
            "1E-300",
            "-1E300",
            "1E-1000",
        ] {
            let value = decimal(text, &mut w).unwrap();
            let bound = rational(&value, &mut w).unwrap();
            assert!(BigRational::from_float(bound.lo).unwrap() <= value.0);
            assert!(BigRational::from_float(bound.hi).unwrap() >= value.0);
        }
        for bits in [1_usize, 52, 53, 54, 80, 256, 1024] {
            for negative in [false, true] {
                let n = (BigInt::from(1) << bits) + BigInt::from(37);
                let n = if negative { -n } else { n };
                let value = N(BigRational::new(n, BigInt::from(3) << (bits / 2)));
                let bound = rational(&value, &mut w).unwrap();
                assert!(BigRational::from_float(bound.lo).unwrap() <= value.0);
                assert!(BigRational::from_float(bound.hi).unwrap() >= value.0);
            }
        }
        // Joint exponent normalization must handle a ratio near one whose
        // numerator and denominator individually exceed binary64's range.
        for negative in [false, true] {
            let huge = BigInt::from(1) << 4096_usize;
            let n = &huge + BigInt::from(1);
            let value = N(BigRational::new(
                if negative { -n } else { n },
                &huge - BigInt::from(1),
            ));
            let bound = rational(&value, &mut w).unwrap();
            assert!(BigRational::from_float(bound.lo).unwrap() <= value.0);
            assert!(BigRational::from_float(bound.hi).unwrap() >= value.0);
        }
        let finite = decimal("1E308", &mut w).unwrap();
        let large_derivative = finite.mul(&N::integer(11), &mut w).unwrap();
        assert!(matches!(
            rational(&large_derivative, &mut w),
            Err(Error::Unresolved)
        ));
        let mut cancelled = Work::new(16_000_000, Some(&|| true)).unwrap();
        assert!(matches!(
            rational(&finite, &mut cancelled),
            Err(Error::Cancelled)
        ));
        let mut short = Work::new(1, None).unwrap();
        assert!(matches!(rational(&finite, &mut short), Err(Error::Limit)));
    }
}
