//! Exact binary-rational geometry leaves for the narrow hardware codec.
//! This avoids silently treating round-trip decimals as exact binary values.
use crate::axial_pipe_thread::Error;
use num_bigint::{BigInt, BigUint};

fn parts(value: f64) -> Result<(u64, i32), Error> {
    if !value.is_finite() {
        return Err(Error::Invalid);
    }
    if value == 0. {
        return Ok((0, 0));
    }
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let mut significand = bits & ((1u64 << 52) - 1);
    let mut power = if exponent == 0 {
        -1074
    } else {
        significand |= 1u64 << 52;
        exponent - 1075
    };
    let reduction = significand.trailing_zeros();
    significand >>= reduction;
    power += reduction as i32;
    Ok((significand, power))
}
/// Bound the scalar's formatting visits and output before allocating BigInts.
pub(crate) fn cost(value: f64) -> Result<(usize, usize), Error> {
    let (n, power) = parts(value)?;
    if n == 0 {
        return Ok((1, 2));
    }
    let bits = 64 - n.leading_zeros() as usize;
    let negative = (-power).max(0) as usize;
    // Rational upper approximations to log10(2) and log10(5).
    let digits = if power < 0 {
        (bits * 30_103 + negative * 69_898).div_ceil(100_000) + 1
    } else {
        ((bits + power as usize) * 30_103).div_ceil(100_000) + 1
    };
    Ok((bits + negative + power.max(0) as usize + 1, digits + 12))
}
pub(crate) fn literal(value: f64) -> Result<String, Error> {
    let (n, power) = parts(value)?;
    if n == 0 {
        return Ok("0.".into());
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let n = BigUint::from(n);
    if power >= 0 {
        Ok(format!("{sign}{}.E0", n << (power as usize)))
    } else {
        let exponent = (-power) as u32;
        let integer = n * BigUint::from(5u8).pow(exponent);
        Ok(format!("{sign}{integer}.E-{exponent}"))
    }
}
struct Decimal {
    negative: bool,
    digits: Vec<u8>,
    exponent: i32,
}
fn decimal(text: &str) -> Option<Decimal> {
    if text.is_empty() || text.len() > 2048 {
        return None;
    }
    let negative = text.starts_with('-');
    let text = text.strip_prefix(['+', '-']).unwrap_or(text);
    let (m, e) = text.split_once('E').unwrap_or((text, "0"));
    let mut digits = Vec::with_capacity(m.len());
    let mut fraction = 0i32;
    let mut dotted = false;
    for b in m.bytes() {
        if b == b'.' {
            if dotted {
                return None;
            }
            dotted = true;
        } else if b.is_ascii_digit() {
            digits.push(b);
            if dotted {
                fraction += 1;
            }
        } else {
            return None;
        }
    }
    let mut exponent = e.parse::<i32>().ok()?.checked_sub(fraction)?;
    let first = digits.iter().position(|b| *b != b'0');
    let Some(first) = first else {
        return (!digits.is_empty()).then_some(Decimal {
            negative: false,
            digits: vec![b'0'],
            exponent: 0,
        });
    };
    while digits.last() == Some(&b'0') {
        digits.pop();
        exponent = exponent.checked_add(1)?;
    }
    digits.drain(..first);
    Some(Decimal {
        negative,
        digits,
        exponent,
    })
}
/// Compare actual decimal arithmetic to the exact binary-rational cache value.
/// Equivalent spellings survive; subULP aliases do not. Long unsupported leaves
/// are refused instead of allocating arbitrary-precision unbounded numbers.
pub(crate) fn exact_leaf(token: &str, value: f64) -> Result<bool, Error> {
    let Some(actual) = decimal(token) else {
        return Ok(false);
    };
    let (n, power) = parts(value)?;
    if n == 0 {
        return Ok(actual.digits == b"0");
    }
    if actual.negative != value.is_sign_negative() {
        return Ok(false);
    }
    // For negative binary powers the odd significand times5^k has no decimal
    // trailing zero. Reject short aliases before expensive exact formatting.
    if power < 0 {
        let bits = 64 - n.leading_zeros() as usize;
        let lower = ((bits - 1) * 30_102 + (-power) as usize * 69_897) / 100_000 + 1;
        if actual.digits.len() < lower {
            return Ok(false);
        }
    }
    let expected = decimal(&literal(value)?).ok_or(Error::Invalid)?;
    Ok(actual.digits == expected.digits && actual.exponent == expected.exponent)
}
/// Check the actual endpoint subtraction as a binary-rational identity. An
/// unrepresentable difference is refused, never replaced by a nearby direction.
pub(crate) fn exact_difference(a: f64, b: f64, difference: f64) -> Result<bool, Error> {
    let inputs = [a, b, difference].map(parts);
    let [a_parts, b_parts, d_parts] = inputs;
    let data = [(a, a_parts?), (b, b_parts?), (difference, d_parts?)];
    let power = data
        .iter()
        .map(|(_, (_, p))| *p)
        .min()
        .ok_or(Error::Invalid)?;
    let values = data.map(|(v, (n, p))| {
        let integer = BigInt::from(n) << ((p - power) as usize);
        if v.is_sign_negative() {
            -integer
        } else {
            integer
        }
    });
    Ok(&values[0] - &values[1] == values[2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::BigRational;
    fn rational(value: f64) -> BigRational {
        let bits = value.to_bits();
        let exponent = ((bits >> 52) & 0x7ff) as i32;
        let mantissa = bits & ((1u64 << 52) - 1);
        let (m, e) = if exponent == 0 {
            (mantissa, -1074)
        } else {
            (mantissa | (1u64 << 52), exponent - 1023 - 52)
        };
        let mut n = BigInt::from(m);
        if bits >> 63 != 0 {
            n = -n;
        }
        if e < 0 {
            BigRational::new(n, BigInt::from(1u8) << (-e as usize))
        } else {
            BigRational::from_integer(n << (e as usize))
        }
    }
    #[test]
    fn emitted_decimal_is_the_exact_binary_value_not_only_a_roundtrip() {
        for v in [
            0.,
            -0.,
            0.1,
            -0.1,
            0.866_025_403_784_438_6,
            1.,
            5.,
            f64::from_bits(1),
            f64::MAX,
        ] {
            let text = literal(v).unwrap();
            assert!(text.len() <= cost(v).unwrap().1);
            assert_eq!(text.parse::<f64>().unwrap(), v);
            let sign = if text.starts_with('-') { -1 } else { 1 };
            let t = text.trim_start_matches('-');
            let (mantissa, exponent) = t.split_once('E').unwrap_or((t, "0"));
            let integer = BigInt::parse_bytes(mantissa.trim_end_matches('.').as_bytes(), 10)
                .unwrap()
                * BigInt::from(sign);
            let e: i32 = exponent.parse().unwrap();
            let actual = if e < 0 {
                BigRational::new(integer, BigInt::from(10u8).pow((-e) as u32))
            } else {
                BigRational::from_integer(integer * BigInt::from(10u8).pow(e as u32))
            };
            assert_eq!(actual, rational(v));
            assert!(exact_leaf(&text, v).unwrap());
        }
        let exact = literal(0.1).unwrap();
        let (m, e) = exact.split_once(".E-").unwrap();
        let alias = format!("{m}1.E-{}", e.parse::<usize>().unwrap() + 1);
        assert_eq!(alias.parse::<f64>().unwrap(), 0.1);
        assert!(!exact_leaf(&alias, 0.1).unwrap());
        assert!(!exact_leaf("0.1", 0.1).unwrap());
        for malformed in ["0.Egarbage", "0.E", "0.E+", "0.0E1E2", "+.", "--0"] {
            assert!(!exact_leaf(malformed, 0.).unwrap());
        }
        assert!(exact_leaf("0.0000000000000000E0", 0.).unwrap());
        assert!(exact_leaf("1.0000000000000000E0", 1.).unwrap());
        assert!(exact_leaf("0.5", 0.5).unwrap());
        assert!(exact_leaf("9", 9.).unwrap());
        assert!(!exact_leaf("9.00000000000000000001", 9.).unwrap());
        assert!(exact_difference(9., 6., 3.).unwrap());
        let d = 1e7 - 0.001;
        assert!(!exact_difference(1e7, 0.001, d).unwrap());
    }
}

#[cfg(test)]
mod malformed_controls {
    use super::*;
    #[test]
    fn exact_geometry_leaves_refuse_nonfinite_and_ambiguous_decimal_grammar() {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(cost(v), Err(Error::Invalid)));
            assert!(matches!(literal(v), Err(Error::Invalid)));
        }
        for text in ["", "0..0", "0.5.0", "-0.5", "0.E1E2"] {
            assert!(!exact_leaf(text, 0.5).unwrap());
        }
        assert!(!exact_leaf(&"0".repeat(2049), 0.).unwrap());
        assert!(exact_leaf("+0.5000E0", 0.5).unwrap());
        assert!(exact_leaf("-0.5000E0", -0.5).unwrap());
        assert!(!exact_difference(0.1, 0.2, -0.1_f64.next_up()).unwrap());
    }
}
