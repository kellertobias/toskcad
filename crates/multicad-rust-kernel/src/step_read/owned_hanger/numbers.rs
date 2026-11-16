//! Exact finite-decimal order; binary64 caches never decide endpoint ownership.
use super::*;
use std::cmp::Ordering;
struct Decimal {
    sign: i8,
    order: i64,
    digits: Vec<u8>,
}
fn decimal(value: Scalar<'_>, b: &mut Budget<'_>) -> Result<Decimal, Error> {
    b.charge(value.text.len() + 1)?;
    let text = value.text.strip_prefix('+').unwrap_or(value.text);
    let negative = text.starts_with('-');
    let text = text.strip_prefix('-').unwrap_or(text);
    let (mantissa, exponent) = text.split_once('E').unwrap_or((text, "0"));
    let exponent = exponent.parse::<i32>().map_err(|_| Error::Unsupported)?;
    let fraction = mantissa.split_once('.').map_or(0, |(_, s)| s.len());
    let mut digits = mantissa.bytes().filter(|&v| v != b'.').collect::<Vec<_>>();
    let Some(first) = digits.iter().position(|&v| v != b'0') else {
        return Ok(Decimal {
            sign: 0,
            order: 0,
            digits: Vec::new(),
        });
    };
    let order = i64::from(exponent) - fraction as i64 + (digits.len() - first) as i64;
    digits.drain(..first);
    Ok(Decimal {
        sign: if negative { -1 } else { 1 },
        order,
        digits,
    })
}
pub(super) fn compare(a: Scalar<'_>, c: Scalar<'_>, b: &mut Budget<'_>) -> Result<Ordering, Error> {
    let a = decimal(a, b)?;
    let c = decimal(c, b)?;
    let sign = a.sign.cmp(&c.sign);
    if sign != Ordering::Equal {
        return Ok(sign);
    }
    if a.sign == 0 {
        return Ok(Ordering::Equal);
    }
    let mut order = a.order.cmp(&c.order);
    if order == Ordering::Equal {
        for i in 0..a.digits.len().max(c.digits.len()) {
            b.charge(1)?;
            order = a
                .digits
                .get(i)
                .copied()
                .unwrap_or(b'0')
                .cmp(&c.digits.get(i).copied().unwrap_or(b'0'));
            if order != Ordering::Equal {
                break;
            }
        }
    }
    Ok(if a.sign < 0 { order.reverse() } else { order })
}
pub(super) fn unit(
    value: Scalar<'_>,
    lo: Scalar<'_>,
    hi: Scalar<'_>,
    b: &mut Budget<'_>,
) -> Result<I, Error> {
    let lower = compare(value, lo, b)?;
    let upper = compare(value, hi, b)?;
    if lower == Ordering::Less || upper == Ordering::Greater {
        return Err(Error::Unsupported);
    }
    if lower == Ordering::Equal {
        return Ok(I::point(0.));
    }
    if upper == Ordering::Equal {
        return Ok(I::point(1.));
    }
    let interval = (value.interval() - lo.interval()) / (hi.interval() - lo.interval());
    if !interval.valid() {
        return Err(Error::Unsupported);
    }
    Ok(I::new(interval.lo.max(0.), interval.hi.min(1.)))
}
pub(super) fn unit_coordinate(value: Scalar<'_>, b: &mut Budget<'_>) -> Result<I, Error> {
    unit(value, Scalar::new("0.", 0.), Scalar::new("1.", 1.), b)
}
pub(super) fn rational(
    value: Scalar<'_>,
    b: &mut Budget<'_>,
) -> Result<num_rational::BigRational, Error> {
    use num_bigint::BigInt;
    b.charge(value.text.len() + 1)?;
    let (mantissa, exponent) = value.text.split_once('E').unwrap_or((value.text, "0"));
    let exponent = exponent.parse::<i32>().map_err(|_| Error::Unsupported)?;
    let fraction = mantissa.split_once('.').map_or(0, |(_, s)| s.len());
    let scale = i64::from(exponent) - fraction as i64;
    if value.text.len() > 16384 || scale.unsigned_abs() > 16384 {
        return Err(Error::Limit);
    }
    let bits = value.text.len() * 4 + scale.unsigned_abs() as usize * 4 + 1;
    let words = bits.div_ceil(64);
    b.charge(32 * (1 + words * words) * (1 + (scale.unsigned_abs().max(1)).ilog2() as usize))?;
    let digits = mantissa.bytes().filter(|&v| v != b'.').collect::<Vec<_>>();
    let integer = BigInt::parse_bytes(&digits, 10).ok_or(Error::Invalid)?;
    let power = BigInt::from(10).pow(scale.unsigned_abs() as u32);
    Ok(if scale >= 0 {
        num_rational::BigRational::from_integer(integer * power)
    } else {
        num_rational::BigRational::new(integer, power)
    })
}
pub(super) fn exact_interval(value: Scalar<'_>, b: &mut Budget<'_>) -> Result<I, Error> {
    let actual = rational(value, b)?;
    let cache = num_rational::BigRational::from_float(value.value).ok_or(Error::Invalid)?;
    Ok(match actual.cmp(&cache) {
        Ordering::Equal => I::point(value.value),
        Ordering::Less => I::new(value.value.next_down(), value.value),
        Ordering::Greater => I::new(value.value, value.value.next_up()),
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    fn s(text: &str) -> Scalar<'_> {
        Scalar::new(text, text.parse().unwrap())
    }
    #[test]
    fn actual_underflowing_decimals_refuse_unresolved_chart_and_excessive_exact_scale() {
        let mut b = Budget {
            work: 0,
            remaining: 16_000_000,
            bytes: 0,
            limit: 64 * 1024 * 1024,
            interrupt: None,
        };
        let (lo, value, hi) = (s("1.E-325"), s("2.E-325"), s("3.E-325"));
        assert_eq!([lo.value, value.value, hi.value], [0.; 3]);
        assert_eq!(compare(lo, value, &mut b).unwrap(), Ordering::Less);
        assert_eq!(compare(value, hi, &mut b).unwrap(), Ordering::Less);
        // Distinct exact decimals cannot justify division by their coincident
        // binary caches. The original unresolved chart remains a refusal.
        assert!(matches!(
            unit(value, lo, hi, &mut b),
            Err(Error::Unsupported)
        ));
        assert!(matches!(
            exact_interval(s("1.E-20000"), &mut b),
            Err(Error::Limit)
        ));
        let retry = exact_interval(s("0.25"), &mut b).unwrap();
        assert_eq!([retry.lo, retry.hi], [0.25; 2]);
    }
    #[test]
    fn exact_decimal_order_rejects_hidden_binary64_domain_overshoot() {
        let mut b = Budget {
            work: 0,
            remaining: 16_000_000,
            bytes: 0,
            limit: 64 * 1024 * 1024,
            interrupt: None,
        };
        for (a, c, expected) in [
            ("0.1", "1000.E-4", Ordering::Equal),
            ("-2.", "-1.", Ordering::Less),
            ("1.E2", "99.", Ordering::Greater),
            ("-0.E9", "0.", Ordering::Equal),
            ("1.01", "1.001", Ordering::Greater),
        ] {
            assert_eq!(compare(s(a), s(c), &mut b).unwrap(), expected);
        }
        let inside = s("24999999999999999999.E-20");
        let outside = s("25000000000000000001.E-20");
        assert_eq!(inside.value, outside.value);
        assert!(unit(inside, s("0."), s("0.25"), &mut b).unwrap().hi <= 1.);
        assert!(matches!(
            unit(outside, s("0."), s("0.25"), &mut b),
            Err(Error::Unsupported)
        ));
        assert!(unit_coordinate(s("0."), &mut b).unwrap().lo == 0.);
        assert!(unit_coordinate(s("1."), &mut b).unwrap().hi == 1.);
        for text in ["0.1", "-3.125", "1.E-9", "0.", "1."] {
            let interval = exact_interval(s(text), &mut b).unwrap();
            let q = rational(s(text), &mut b).unwrap();
            assert!(num_rational::BigRational::from_float(interval.lo).unwrap() <= q);
            assert!(num_rational::BigRational::from_float(interval.hi).unwrap() >= q);
        }
    }
}
