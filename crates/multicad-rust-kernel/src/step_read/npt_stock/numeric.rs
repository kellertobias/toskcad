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
pub(super) fn expected(s: Scalar<'_>, value: f64, budget: &mut Budget<'_>) -> Result<bool> {
    equal(s.spelling, &format!("{value:.16E}"), budget)
}
