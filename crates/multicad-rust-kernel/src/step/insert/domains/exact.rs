//! Ephemeral exact arithmetic inside the preflighted domain scratch arena.
use super::*;
use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::ToPrimitive;
#[derive(Clone)]
pub(crate) struct N(pub BigRational);
impl N {
    pub fn integer(n: i64) -> Self {
        Self(BigRational::from_integer(n.into()))
    }
    pub fn parts(sign: i8, n: &[u8], d: &[u8], w: &mut Work<'_>) -> Result<Self, Error> {
        let words = n.len().max(d.len()).div_ceil(8);
        w.charge(1 + words * words)?;
        let sign = match sign {
            -1 => Sign::Minus,
            0 => Sign::NoSign,
            _ => Sign::Plus,
        };
        Ok(Self(BigRational::new(
            BigInt::from_bytes_le(sign, n),
            BigInt::from_bytes_le(Sign::Plus, d),
        )))
    }
    pub fn view(
        v: crate::authored_internal_treatment::ExactBoundaryRational<'_>,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        Self::parts(
            v.signum(),
            v.numerator_magnitude_le(),
            v.denominator_magnitude_le(),
            w,
        )
    }
    pub fn enclosure(&self, w: &mut Work<'_>) -> Result<I, Error> {
        let bits = self.0.numer().bits().max(self.0.denom().bits());
        let words = (bits as usize).div_ceil(64);
        w.charge(1 + words * words)?;
        let proposal = self
            .0
            .to_f64()
            .filter(|x| x.is_finite())
            .ok_or(Error::Unresolved)?;
        self.enclose_proposal(proposal, w)
    }
    pub fn enclose_proposal(&self, proposal: f64, w: &mut Work<'_>) -> Result<I, Error> {
        let point = Self(BigRational::from_float(proposal).ok_or(Error::Unresolved)?);
        match self.cmp(&point, w)? {
            Ordering::Equal => Ok(I::point(proposal)),
            Ordering::Less => {
                let low = proposal.next_down();
                let exact = BigRational::from_float(low).ok_or(Error::Unresolved)?;
                if self.cmp(&Self(exact), w)? == Ordering::Less {
                    return Err(Error::Unresolved);
                }
                Ok(I::new(low, proposal))
            }
            Ordering::Greater => {
                let high = proposal.next_up();
                let exact = BigRational::from_float(high).ok_or(Error::Unresolved)?;
                if self.cmp(&Self(exact), w)? == Ordering::Greater {
                    return Err(Error::Unresolved);
                }
                Ok(I::new(proposal, high))
            }
        }
    }
    fn preflight(&self, other: &Self, operation: u8, w: &mut Work<'_>) -> Result<(), Error> {
        let [a, b, c, d] = [
            self.0.numer().bits(),
            self.0.denom().bits(),
            other.0.numer().bits(),
            other.0.denom().bits(),
        ];
        let upper = match operation {
            0 => ((a + d).max(c + b) + 1).max(b + d),
            1 => (a + c).max(b + d),
            _ => (a + d).max(b + c),
        };
        // Numerator/denominator cross products, normalization and comparisons
        // cannot exceed this preflight. The16MiB scratch covers all <=12 cell
        // inputs and <=32 live arithmetic values plus bounded bigint temporaries.
        if upper > 16384 {
            return Err(Error::Limit);
        }
        let words = (upper as usize).div_ceil(64);
        w.charge(1 + words * words)
    }
    pub fn cmp(&self, b: &Self, w: &mut Work<'_>) -> Result<std::cmp::Ordering, Error> {
        self.preflight(b, 0, w)?;
        Ok(self.0.cmp(&b.0))
    }
    pub fn add(&self, b: &Self, w: &mut Work<'_>) -> Result<Self, Error> {
        self.preflight(b, 0, w)?;
        Ok(Self(&self.0 + &b.0))
    }
    pub fn sub(&self, b: &Self, w: &mut Work<'_>) -> Result<Self, Error> {
        self.preflight(b, 0, w)?;
        Ok(Self(&self.0 - &b.0))
    }
    pub fn mul(&self, b: &Self, w: &mut Work<'_>) -> Result<Self, Error> {
        self.preflight(b, 1, w)?;
        Ok(Self(&self.0 * &b.0))
    }
    pub fn div(&self, b: &Self, w: &mut Work<'_>) -> Result<Self, Error> {
        self.preflight(b, 2, w)?;
        if b.0 == BigRational::from_integer(0.into()) {
            return Err(Error::Unresolved);
        }
        Ok(Self(&self.0 / &b.0))
    }
}
