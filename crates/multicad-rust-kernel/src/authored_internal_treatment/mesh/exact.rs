//! Bounded exact rational source parameters; never geometric tolerance signs.
use super::*;
use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::{ToPrimitive, Zero};
use std::{cell::Cell, cmp::Ordering, rc::Rc};
const MAX_BITS: u64 = 8192;
struct Value {
    number: BigRational,
    float: Cell<Option<(f64, I)>>,
}
#[derive(Clone)]
pub(super) struct Q(Rc<Value>);
impl Q {
    pub fn float(value: f64, work: &mut Work<'_>) -> Result<Self, Error> {
        work.step()?;
        if !value.is_finite() {
            return Err(Error::Unresolved);
        }
        // Every finite binary64 rational has <=1075 bits per operand.
        let number = BigRational::from_float(value).ok_or(Error::Unresolved)?;
        Ok(Self::make(number, Some((value, I::point(value)))))
    }
    pub fn integer(value: i64, work: &mut Work<'_>) -> Result<Self, Error> {
        work.step()?;
        let f = value as f64;
        Ok(Self::make(
            BigRational::from_integer(value.into()),
            (value.unsigned_abs() <= 1_u64 << 53).then_some((f, I::point(f))),
        ))
    }
    fn make(number: BigRational, float: Option<(f64, I)>) -> Self {
        // Only finite binary64 (<=1075 bits), i64 or preflight-qualified
        // operations reach here. Product/quotient operands have at most the
        // summed bit lengths; addition adds one carry bit. Rational reduction
        // cannot increase either bound, so MAX_BITS is already proved.
        let float = if number.is_zero() {
            Some((0., I::point(0.)))
        } else {
            float
        };
        Self(Rc::new(Value {
            number,
            float: Cell::new(float),
        }))
    }
    fn preflight(
        &self,
        other: &Self,
        add: bool,
        divide: bool,
        work: &mut Work<'_>,
    ) -> Result<(), Error> {
        let [a, b] = [self.0.number.numer().bits(), self.0.number.denom().bits()];
        let [c, d] = [other.0.number.numer().bits(), other.0.number.denom().bits()];
        let upper = if add {
            ((a + d).max(c + b) + 1).max(b + d)
        } else if divide {
            (a + d).max(b + c)
        } else {
            (a + c).max(b + d)
        };
        if upper > MAX_BITS {
            return Err(Error::Limit);
        }
        // Charge limb-product/gcd complexity before library arithmetic. Even
        // extreme admitted float exponents cannot trigger unbounded operands.
        let words = (upper as usize).div_ceil(64);
        work.charge(1 + words * words)
    }
    pub fn add(&self, other: &Self, work: &mut Work<'_>) -> Result<Self, Error> {
        self.preflight(other, true, false, work)?;
        Ok(Self::make(&self.0.number + &other.0.number, None))
    }
    pub fn sub(&self, other: &Self, work: &mut Work<'_>) -> Result<Self, Error> {
        self.preflight(other, true, false, work)?;
        Ok(Self::make(&self.0.number - &other.0.number, None))
    }
    pub fn mul(&self, other: &Self, work: &mut Work<'_>) -> Result<Self, Error> {
        self.preflight(other, false, false, work)?;
        Ok(Self::make(&self.0.number * &other.0.number, None))
    }
    pub fn div(&self, other: &Self, work: &mut Work<'_>) -> Result<Self, Error> {
        work.step()?;
        if other.0.number.is_zero() {
            return Err(Error::Unresolved);
        }
        self.preflight(other, false, true, work)?;
        Ok(Self::make(&self.0.number / &other.0.number, None))
    }
    pub fn cmp(&self, other: &Self, work: &mut Work<'_>) -> Result<Ordering, Error> {
        work.step()?;
        if Rc::ptr_eq(&self.0, &other.0) {
            return Ok(Ordering::Equal);
        }
        let bits = self
            .0
            .number
            .numer()
            .bits()
            .max(self.0.number.denom().bits())
            .max(other.0.number.numer().bits())
            .max(other.0.number.denom().bits());
        work.charge(1 + bits.div_ceil(64) as usize)?;
        if self.0.number.numer() == other.0.number.numer()
            && self.0.number.denom() == other.0.number.denom()
        {
            return Ok(Ordering::Equal);
        }
        self.preflight(other, true, false, work)?;
        Ok(self.0.number.cmp(&other.0.number))
    }
    pub fn sign(&self) -> Ordering {
        match self.0.number.numer().sign() {
            Sign::Minus => Ordering::Less,
            Sign::NoSign => Ordering::Equal,
            Sign::Plus => Ordering::Greater,
        }
    }
    /// Ephemeral exact query value. Callers compact before arena retention.
    pub fn rational(&self, work: &mut Work<'_>) -> Result<BigRational, Error> {
        work.charge(
            1 + (self.0.number.numer().bits() + self.0.number.denom().bits()).div_ceil(64) as usize,
        )?;
        Ok(self.0.number.clone())
    }
    /// A newly allocated compact representation has no historical BigInt
    /// capacity. Headers and magnitude bytes are charged before allocation;
    /// callers account its actual retained boxed length alongside table storage.
    pub fn encode(values: &[&Self], work: &mut Work<'_>) -> Result<Box<[u8]>, Error> {
        let mut length = 0usize;
        for value in values {
            let n = value.0.number.numer().bits().div_ceil(8).max(1) as usize;
            let d = value.0.number.denom().bits().div_ceil(8).max(1) as usize;
            length = length.checked_add(5 + n + d).ok_or(Error::Limit)?;
        }
        work.charge(1 + length.div_ceil(8))?;
        let mut out = Vec::new();
        out.try_reserve_exact(length).map_err(|_| Error::Limit)?;
        for value in values {
            let (sign, n) = value.0.number.numer().to_bytes_le();
            let (_, d) = value.0.number.denom().to_bytes_le();
            out.push(match sign {
                Sign::Minus => 0,
                Sign::NoSign => 1,
                Sign::Plus => 2,
            });
            out.extend_from_slice(&(n.len() as u16).to_le_bytes());
            out.extend_from_slice(&(d.len() as u16).to_le_bytes());
            out.extend_from_slice(&n);
            out.extend_from_slice(&d);
        }
        // num-bigint emits one zero byte for zero magnitudes on some versions;
        // actual retained length, rather than the preflight estimate, is final.
        work.step()?;
        Ok(out.into_boxed_slice())
    }
    /// Internal compact coordinate token, never a geometry/import descriptor.
    pub fn decode(bytes: &[u8], work: &mut Work<'_>) -> Result<Self, Error> {
        work.step()?;
        let header = bytes.get(..5).ok_or(Error::Unresolved)?;
        let n = usize::from(u16::from_le_bytes([header[1], header[2]]));
        let d = usize::from(u16::from_le_bytes([header[3], header[4]]));
        if n > MAX_BITS as usize / 8 || d > MAX_BITS as usize / 8 {
            return Err(Error::Limit);
        }
        if n == 0 || d == 0 || bytes.len() != 5 + n + d {
            return Err(Error::Unresolved);
        }
        let sign = match header[0] {
            0 => Sign::Minus,
            1 => Sign::NoSign,
            2 => Sign::Plus,
            _ => return Err(Error::Unresolved),
        };
        let words = n.max(d).div_ceil(8);
        work.charge(1 + words * words)?;
        let numerator = BigInt::from_bytes_le(sign, &bytes[5..5 + n]);
        let denominator = BigInt::from_bytes_le(Sign::Plus, &bytes[5 + n..]);
        if denominator.is_zero() {
            return Err(Error::Unresolved);
        }
        Ok(Self::make(BigRational::new(numerator, denominator), None))
    }
    /// No rounding assumption about num-rational's nominal conversion is used:
    /// both chosen f64 bounds are checked against the original exact ratio.
    /// Only successful enclosures are cached; cancellation precedes reuse.
    pub fn enclosure(&self, work: &mut Work<'_>) -> Result<(f64, I), Error> {
        work.step()?;
        if let Some(value) = self.0.float.get() {
            return Ok(value);
        }
        let nominal = self
            .0
            .number
            .to_f64()
            .filter(|f| f.is_finite())
            .ok_or(Error::Unresolved)?;
        let out = self.checked_proposal(nominal, work)?;
        work.step()?;
        self.0.float.set(Some(out));
        Ok(out)
    }
    /// Independently validate an external numeric proposal against this exact
    /// source value. No correctness/rounding assumption about its producer.
    pub fn checked_proposal(&self, nominal: f64, work: &mut Work<'_>) -> Result<(f64, I), Error> {
        let center = Self::float(nominal, work)?;
        let order = center.cmp(self, work)?;
        let (lo, hi) = match order {
            Ordering::Equal => (nominal, nominal),
            Ordering::Less => (nominal, nominal.next_up()),
            Ordering::Greater => (nominal.next_down(), nominal),
        };
        let lower = Self::float(lo, work)?;
        let upper = Self::float(hi, work)?;
        if lower.cmp(self, work)? == Ordering::Greater || upper.cmp(self, work)? == Ordering::Less {
            return Err(Error::Unresolved);
        }
        Ok((nominal, I::new(lo, hi)))
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn work() -> Work<'static> {
        Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        }
    }
    #[test]
    fn checked_conversion_handles_subnormal_underflow_exact_identity_and_real_overflow() {
        let w = &mut work();
        let tiny = Q::float(f64::from_bits(1), w).unwrap();
        let two = Q::integer(2, w).unwrap();
        let half = tiny.div(&two, w).unwrap();
        let (nominal, bound) = half.enclosure(w).unwrap();
        assert_eq!(nominal, 0.);
        assert_eq!(bound.lo, 0.);
        assert_eq!(bound.hi, f64::from_bits(1));
        assert_eq!(
            half.mul(&two, w).unwrap().cmp(&tiny, w).unwrap(),
            Ordering::Equal
        );
        let large = Q::float(1e300, w).unwrap();
        let huge = large.mul(&large, w).unwrap();
        assert!(matches!(huge.enclosure(w), Err(Error::Unresolved)));
        let fourth = huge.mul(&huge, w).unwrap();
        let eighth = fourth.mul(&fourth, w).unwrap();
        assert!(matches!(eighth.mul(&eighth, w), Err(Error::Limit)));
        assert!(matches!(
            tiny.div(&Q::integer(0, w).unwrap(), w),
            Err(Error::Unresolved)
        ));
        let zero = Q::float(-0., w).unwrap();
        let plus = Q::integer(0, w).unwrap();
        assert_eq!(
            Q::encode(&[&zero], w).unwrap(),
            Q::encode(&[&plus], w).unwrap()
        );
        let mut cancelled = Work {
            remaining: 0,
            used: 0,
            interrupt: Some(&|| true),
        };
        assert!(matches!(
            tiny.enclosure(&mut cancelled),
            Err(Error::Cancelled)
        ));
        let mut limited = Work {
            remaining: 0,
            used: 0,
            interrupt: None,
        };
        assert!(matches!(tiny.enclosure(&mut limited), Err(Error::Limit)));
        assert_eq!(tiny.enclosure(w).unwrap().1.lo, f64::from_bits(1));
    }
}
#[cfg(test)]
mod malformed_tests {
    use super::*;
    #[test]
    fn external_numeric_inputs_refuse_nonfinite_values_and_bound_integer_rounding() {
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(Q::float(value, w), Err(Error::Unresolved)));
        }
        let integer = Q::integer(i64::MAX, w).unwrap();
        let (_, bounds) = integer.enclosure(w).unwrap();
        assert!(bounds.lo < (i64::MAX as f64) && bounds.hi == i64::MAX as f64);
        assert_eq!(
            integer
                .sub(&Q::integer(i64::MAX, w).unwrap(), w)
                .unwrap()
                .sign(),
            Ordering::Equal
        );
    }
}
#[cfg(test)]
mod proposal_tests {
    use super::*;
    #[test]
    fn exact_source_independently_refuses_wrong_external_numeric_proposals() {
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let third = Q::integer(1, w)
            .unwrap()
            .div(&Q::integer(3, w).unwrap(), w)
            .unwrap();
        let (n, b) = third.checked_proposal(1. / 3., w).unwrap();
        assert_eq!(n, 1. / 3.);
        assert!(b.lo <= n && n <= b.hi);
        for wrong in [0., 1.] {
            assert!(matches!(
                third.checked_proposal(wrong, w),
                Err(Error::Unresolved)
            ));
        }
        assert_eq!(Q::integer(-1, w).unwrap().sign(), Ordering::Less);
    }
}
#[cfg(test)]
mod decode_tests {
    use super::*;
    #[test]
    fn compact_exact_parameter_tokens_decode_with_limits_and_zero_denominator_refusal() {
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        for integer in [-1, 0] {
            let value = Q::integer(integer, w).unwrap();
            let token = Q::encode(&[&value], w).unwrap();
            assert_eq!(
                Q::decode(&token, w).unwrap().cmp(&value, w).unwrap(),
                Ordering::Equal
            );
        }
        let value = Q::float(0.25, w).unwrap();
        let encoded = Q::encode(&[&value], w).unwrap();
        assert_eq!(
            Q::decode(&encoded, w).unwrap().cmp(&value, w).unwrap(),
            Ordering::Equal
        );
        assert!(matches!(Q::decode(&[], w), Err(Error::Unresolved)));
        let mut bad = encoded.to_vec();
        bad[0] = 9;
        assert!(matches!(Q::decode(&bad, w), Err(Error::Unresolved)));
        let mut bad = encoded.to_vec();
        bad[1] = 0;
        bad[2] = 0;
        assert!(matches!(Q::decode(&bad, w), Err(Error::Unresolved)));
        assert!(matches!(
            Q::decode(&encoded[..encoded.len() - 1], w),
            Err(Error::Unresolved)
        ));
        let mut bad = encoded.to_vec();
        let denominator = usize::from(u16::from_le_bytes([bad[3], bad[4]]));
        let end = bad.len();
        bad[end - denominator..].fill(0);
        assert!(matches!(Q::decode(&bad, w), Err(Error::Unresolved)));
        assert!(matches!(Q::decode(&[2, 1, 4, 1, 0], w), Err(Error::Limit)));
    }
}
