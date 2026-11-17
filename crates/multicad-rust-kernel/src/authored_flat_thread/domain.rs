//! Explicit finite periodic-domain authority, never tolerance equality.
use super::*;
use crate::Rational;
use num_rational::BigRational as Q;
pub(super) fn interval(end: Rational) -> I {
    let (n, d) = end.ratio();
    let bits = n.unsigned_abs();
    if (d as u64).is_power_of_two()
        && (bits == 0 || 64 - bits.leading_zeros() - bits.trailing_zeros() <= 53)
    {
        I::point(end.as_f64())
    } else {
        // Rational::as_f64 performs one correctly rounded exact division.
        let x = end.as_f64();
        I::new(x.next_down(), x.next_up())
    }
}
pub(super) fn validate(spec: Spec, end: Rational) -> Result<I, Error> {
    let (n, d) = end.ratio();
    let q = Q::new(n.into(), d.into());
    if q < Q::new(1.into(), 1000.into()) || q > Q::from_integer(10_000_000.into()) {
        return Err(Error::Invalid);
    }
    // The actual finite phase span, not a rounded quotient, owns the existing
    // turn/cardinality ceiling. Original Spec admission is checked separately.
    let pitch = Q::from_float(spec.pitch_mm).ok_or(Error::Invalid)?;
    if q / pitch + Q::from_integer((2 * i64::from(spec.starts)).into())
        > Q::from_integer(128.into())
    {
        return Err(Error::Limit);
    }
    Ok(interval(end))
}
pub(super) fn contains(end: Rational, z: f64) -> bool {
    if !z.is_finite() || z < 0. {
        return false;
    }
    let (n, d) = end.ratio();
    // Binary64 inputs and the bounded canonical ratio have bounded exact
    // integer representations. No neighbouring float admits an exterior point.
    Q::from_float(z).is_some_and(|p| p <= Q::new(n.into(), d.into()))
}
