//! Exact finite STEP REAL text; no conversion through a floating-point value.
use super::*;
mod cache;
pub(in crate::step::hanger) use cache::DecimalCache;
use num_bigint::BigInt;
use num_traits::{One, Zero};

pub struct FiniteDecimal(String);
impl FiniteDecimal {
    pub fn as_str(&self) -> &str {
        &self.0
    }
    /// Includes the result descriptor and actual retained String capacity.
    pub fn retained_bytes(&self) -> usize {
        size_of::<Self>() + self.0.capacity()
    }
}

impl EncodedScalar {
    /// Encode this exact scalar as a STEP REAL when its reduced denominator
    /// contains only factors two and five. This is a numeric operation; it
    /// establishes no finite-face, file or material admission.
    pub fn exact_step_real(
        &self,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<FiniteDecimal, Error> {
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(
            options,
            size_of::<Self>() + self.n.len() + self.d.len(),
            interrupt,
        )?;
        let number = self.number(&mut w)?;
        decimal(&number, &mut m, &mut w)
    }
}

/// The caller accounts its borrowed rational. This operation charges scratch
/// and retained output to that SAME ledger. Failure releases its own memory;
/// consumed work is never refunded. Successful output remains charged.
pub(in crate::step::hanger) fn decimal(
    value: &N,
    memory: &mut Memory,
    work: &mut Work<'_>,
) -> Result<FiniteDecimal, Error> {
    let base = memory.used;
    let result = encode(value, memory, work);
    if result.is_err() {
        memory.used = base;
    }
    result
}

fn encode(value: &N, m: &mut Memory, w: &mut Work<'_>) -> Result<FiniteDecimal, Error> {
    w.charge(1)?;
    let numerator_bits = value.0.numer().bits() as usize;
    let denominator_bits = value.0.denom().bits() as usize;
    if numerator_bits.max(denominator_bits) > 16384 {
        return Err(Error::Limit);
    }
    // For d=2^a*5^b, k=max(a,b)<=bits(d). Scaling to denominator
    // 10^k grows the numerator by at most 3*k bits. Reserve conservatively
    // for the cloned divisor, scaling powers/products, bigint arithmetic
    // temporaries, decimal conversion and its temporary text. The 128-fold
    // reserve exceeds these simultaneously live linear-size buffers.
    let bits = numerator_bits + 3 * denominator_bits + 1;
    let scratch = 128 * (256 + bits.div_ceil(8) + denominator_bits.div_ceil(8));
    m.add(scratch)?;
    let words = denominator_bits.div_ceil(64);
    w.charge(1 + words * words)?;
    let mut denominator = value.0.denom().clone();
    let twos = denominator.trailing_zeros().unwrap_or(0) as usize;
    denominator >>= twos;
    let mut fives = 0;
    loop {
        let words = (denominator.bits() as usize).div_ceil(64);
        // Charge division/remainder before allocating their bigint results.
        w.charge(4 * (1 + words * words))?;
        if !(&denominator % 5_u32).is_zero() {
            break;
        }
        denominator /= 5_u32;
        fives += 1;
    }
    if !denominator.is_one() {
        return Err(Error::Unresolved);
    }
    let scale = twos.max(fives);
    let words = bits.div_ceil(64);
    // Bounded exponentiation, multiplication and decimal conversion. Powers
    // have <=16384 exponent and the preflighted <=65537-bit result bound.
    w.charge(16 * (1 + words * words) * (1 + bits.ilog2() as usize))?;
    let mut integer = value.0.numer().clone();
    integer <<= scale - twos;
    integer *= BigInt::from(5_u32).pow((scale - fives) as u32);
    let text = emit(&integer, scale, bits, m, w)?;
    m.used -= scratch;
    Ok(text)
}
// The caller precharges numeric radix conversion and holds its scratch reserve.
fn emit(
    integer: &BigInt,
    scale: usize,
    digits_bound: usize,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<FiniteDecimal, Error> {
    let digits = integer.to_str_radix(10);
    let exponent = scale.to_string();
    // A binary-bit bound also bounds decimal digits. Reserve before writing;
    // count the allocator's actual capacity while all scratch is still live.
    let capacity = digits_bound + 2 + 2 + exponent.len();
    m.add(size_of::<FiniteDecimal>() + capacity)?;
    w.charge(1 + capacity)?;
    // This bounded allocation follows the same preflight as the bigint and
    // decimal scratch allocations; process allocation failure is not a
    // geometric refusal. Count any extra capacity before writing.
    let mut text = String::with_capacity(capacity);
    m.add(text.capacity() - capacity)?;
    text.push_str(&digits);
    text.push('.');
    if scale != 0 {
        text.push_str("E-");
        text.push_str(&exponent);
    }
    w.charge(1)?;
    Ok(FiniteDecimal(text))
}

#[cfg(test)]
mod tests;
