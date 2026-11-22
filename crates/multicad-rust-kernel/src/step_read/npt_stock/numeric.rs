//! Original decimal leaves remain the geometry. Binary64 is an outward cache.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct Scalar<'input> {
    value: f64,
    spelling: &'input str,
    bounds: [f64; 2],
}
impl<'input> Scalar<'input> {
    pub fn value(&self) -> f64 {
        self.value
    }
    pub fn spelling(&self) -> &'input str {
        self.spelling
    }
    pub fn bounds(&self) -> [f64; 2] {
        self.bounds
    }
}
pub(super) struct Numbers<'input> {
    tokens: Vec<&'input str>,
    next: usize,
}
impl<'input> Numbers<'input> {
    pub fn new(entity: &Entity<'input>, budget: &mut Budget<'_>) -> Result<Self> {
        let source = entity
            .source
            .expect("neutral Part21 arena retains original records");
        budget.step(source.len() + 1)?;
        let tokens = part21::numeric_spellings(entity).ok_or(Error::Invalid)?;
        Ok(Self { tokens, next: 0 })
    }
    pub fn scalar(&mut self, value: &Value<'_>, budget: &mut Budget<'_>) -> Result<Scalar<'input>> {
        let Value::Number(value) = value else {
            return Err(Error::Invalid);
        };
        let spelling = *self.tokens.get(self.next).ok_or(Error::Invalid)?;
        self.next += 1;
        budget.step(4 + if *value == 0. { spelling.len() } else { 0 })?;
        let zero = *value == 0.
            && spelling
                .split('E')
                .next()
                .expect("numeric mantissa")
                .bytes()
                .filter(|b| !matches!(b, b'+' | b'-' | b'.'))
                .all(|b| b == b'0');
        let integer = value.fract() == 0.
            && value.abs() <= 9_007_199_254_740_992.
            && equal(spelling, &format!("{value:.16E}"), budget)?;
        let bounds = if zero || integer {
            [*value, *value]
        } else {
            [value.next_down(), value.next_up()]
        };
        Ok(Scalar {
            value: *value,
            spelling,
            bounds,
        })
    }
    pub fn integer(&mut self, value: &Value<'_>, budget: &mut Budget<'_>) -> Result<usize> {
        let scalar = self.scalar(value, budget)?;
        let text = scalar.spelling.strip_prefix('+').unwrap_or(scalar.spelling);
        budget.step(text.len() + 1)?;
        if text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
            return Err(Error::Invalid);
        }
        text.parse().map_err(|_| Error::Limit)
    }
    pub fn finish(self) -> Result<()> {
        if self.next == self.tokens.len() {
            Ok(())
        } else {
            Err(Error::Invalid)
        }
    }
}
struct Decimal {
    negative: bool,
    digits: Vec<u8>,
    exponent: i64,
}
fn normalized(text: &str, budget: &mut Budget<'_>) -> Result<Decimal> {
    budget.step(text.len() + 1)?;
    let negative = text.starts_with('-');
    let text = text.strip_prefix(['+', '-']).unwrap_or(text);
    let (mantissa, exponent) = text.split_once('E').unwrap_or((text, "0"));
    let fraction = mantissa.split_once('.').map_or(0, |(_, s)| s.len());
    let mut digits = Vec::new();
    digits
        .try_reserve_exact(mantissa.len())
        .map_err(|_| Error::Limit)?;
    digits.extend(mantissa.bytes().filter(|b| *b != b'.'));
    let leading = digits.iter().take_while(|b| **b == b'0').count();
    digits.drain(..leading);
    if digits.is_empty() {
        return Ok(Decimal {
            negative: false,
            digits,
            exponent: 0,
        });
    }
    let trailing = digits.iter().rev().take_while(|b| **b == b'0').count();
    digits.truncate(digits.len() - trailing);
    let exponent = exponent
        .parse::<i64>()
        .map_err(|_| Error::Limit)?
        .checked_sub(fraction as i64)
        .and_then(|e| e.checked_add(trailing as i64))
        .ok_or(Error::Limit)?;
    Ok(Decimal {
        negative,
        digits,
        exponent,
    })
}
pub(super) fn equal(a: &str, b: &str, budget: &mut Budget<'_>) -> Result<bool> {
    budget.step(a.len() + b.len() + 1)?;
    if a == b {
        return Ok(true);
    }
    let a = normalized(a, budget)?;
    let b = normalized(b, budget)?;
    Ok(a.negative == b.negative && a.exponent == b.exponent && a.digits == b.digits)
}
// Exact decimal addition in the original plane chart. No binary64 cache or
// interval overlap can authorize a changed PCURVE. Extreme exponent spreads
// refuse within the existing arena rather than allocating an unbounded integer.
pub(super) fn translated_equal(
    uv: Scalar<'_>,
    anchor: Scalar<'_>,
    physical: Scalar<'_>,
    budget: &mut Budget<'_>,
) -> Result<bool> {
    if anchor.bounds == [0., 0.] {
        return equal(uv.spelling, physical.spelling, budget);
    }
    let values = [
        normalized(uv.spelling, budget)?,
        normalized(anchor.spelling, budget)?,
        normalized(physical.spelling, budget)?,
    ];
    let Some(lower) = values
        .iter()
        .filter(|v| !v.digits.is_empty())
        .map(|v| v.exponent)
        .min()
    else {
        return Ok(true);
    };
    let upper =
        values
            .iter()
            .filter(|v| !v.digits.is_empty())
            .try_fold(lower, |upper, value| {
                value
                    .exponent
                    .checked_add(value.digits.len() as i64)
                    .map(|end| upper.max(end))
                    .ok_or(Error::Limit)
            })?;
    let width = upper
        .checked_sub(lower)
        .and_then(|n| n.checked_add(1))
        .and_then(|n| usize::try_from(n).ok())
        .filter(|n| *n <= 8192)
        .ok_or(Error::Limit)?;
    let mut positive = budget.allocate::<u8>(width)?;
    let mut negative = budget.allocate::<u8>(width)?;
    budget.step(width.checked_mul(2).ok_or(Error::Limit)?)?;
    positive.resize(width, 0);
    negative.resize(width, 0);
    for (index, value) in values.iter().enumerate() {
        if value.digits.is_empty() {
            continue;
        }
        let offset = usize::try_from(value.exponent.checked_sub(lower).ok_or(Error::Limit)?)
            .map_err(|_| Error::Limit)?;
        let target = if value.negative ^ (index == 2) {
            &mut negative
        } else {
            &mut positive
        };
        // Each source digit, possible propagated carry and destination access
        // is charged before arithmetic. At most width carry positions exist.
        budget.step(
            value
                .digits
                .len()
                .checked_mul(3)
                .and_then(|n| n.checked_add(width))
                .ok_or(Error::Limit)?,
        )?;
        let mut carry = 0;
        for (i, digit) in value.digits.iter().rev().enumerate() {
            let sum = target[offset + i] + (digit - b'0') + carry;
            target[offset + i] = sum % 10;
            carry = sum / 10;
        }
        let mut i = offset + value.digits.len();
        while carry != 0 {
            let sum = target[i] + carry;
            target[i] = sum % 10;
            carry = sum / 10;
            i += 1;
        }
    }
    budget.step(width)?;
    let same = positive == negative;
    budget.release(&positive);
    budget.release(&negative);
    Ok(same)
}
pub(super) fn expected(s: Scalar<'_>, value: f64, budget: &mut Budget<'_>) -> Result<bool> {
    equal(s.spelling, &format!("{value:.16E}"), budget)
}

// The owned geometry copier preserves the same original component spelling,
// binary64 cache and outward bounds. This constructs a borrowing view of that
// immutable copy, never a new numeric or material admission.
pub(super) fn retained_scalar(value: f64, bounds: [f64; 2], spelling: &str) -> Scalar<'_> {
    Scalar {
        value,
        bounds,
        spelling,
    }
}

#[cfg(test)]
mod tests;
