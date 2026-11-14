//! Original numeric leaves remain authoritative; binary64 is an outward cache.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct Scalar<'input> {
    value: f64,
    spelling: &'input str,
}
impl<'input> Scalar<'input> {
    pub fn value(&self) -> f64 {
        self.value
    }
    pub fn spelling(&self) -> &'input str {
        self.spelling
    }
    pub fn bounds(&self) -> [f64; 2] {
        // This proves an exact literal integer, rather than assuming every parsed
        // floating integer had that exact source value (e.g. 11.0000000000000001).
        if self.is_zero()
            || (self.value.fract() == 0.
                && self.value.abs() <= 9_007_199_254_740_992.
                && self.spelling == format!("{:.16E}", self.value))
        {
            [self.value, self.value]
        } else {
            [self.value.next_down(), self.value.next_up()]
        }
    }
    pub(super) fn is_zero(&self) -> bool {
        if self.value != 0. {
            return false;
        }
        self.spelling
            .split('E')
            .next()
            .unwrap()
            .bytes()
            .filter(|b| !matches!(b, b'+' | b'-' | b'.'))
            .all(|b| b == b'0')
    }
    pub(super) fn is_one(&self) -> bool {
        if self.value != 1. {
            return false;
        }
        let text = self.spelling.strip_prefix('+').unwrap_or(self.spelling);
        let (mantissa, exponent) = text.split_once('E').unwrap_or((text, "0"));
        let Ok(exponent) = exponent.parse::<i32>() else {
            return false;
        };
        let fraction = mantissa.split_once('.').map_or(0, |(_, s)| s.len());
        let mut seen = false;
        let mut trailing = 0usize;
        for digit in mantissa.bytes().filter(|&b| b != b'.') {
            if seen {
                if digit != b'0' {
                    return false;
                }
                trailing += 1;
            } else if digit != b'0' {
                if digit != b'1' {
                    return false;
                }
                seen = true;
            }
        }
        seen && i64::from(exponent) - fraction as i64 + trailing as i64 == 0
    }
}
pub(super) struct Numeric<'input> {
    tokens: Vec<&'input str>,
    next: usize,
}
impl<'input> Numeric<'input> {
    pub fn new(entity: &Entity<'input>, w: &mut Ledger<'_>) -> Result<Self> {
        let source = entity.source.ok_or(Error::Invalid)?;
        w.step(source.len())?;
        let tokens = part21::numeric_spellings(entity).ok_or(Error::Invalid)?;
        Ok(Self { tokens, next: 0 })
    }
    pub fn scalar(&mut self, value: &Value<'_>) -> Result<Scalar<'input>> {
        let Value::Number(value) = value else {
            return Err(Error::Invalid);
        };
        if !value.is_finite() {
            return Err(Error::Invalid);
        }
        let spelling = *self.tokens.get(self.next).ok_or(Error::Invalid)?;
        self.next += 1;
        Ok(Scalar {
            value: *value,
            spelling,
        })
    }
    pub fn integer(&mut self, value: &Value<'_>) -> Result<u32> {
        let s = self.scalar(value)?;
        let text = s.spelling.strip_prefix('+').unwrap_or(s.spelling);
        if text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
            return Err(Error::Invalid);
        }
        let text = text.trim_start_matches('0');
        if text.is_empty() {
            return Ok(0);
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
