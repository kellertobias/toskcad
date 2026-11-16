//! Bounded query-local mathematical constant cache. It never owns topology.
use super::*;
const ENTRIES: usize = 128;
struct Power {
    exponent: usize,
    bytes: Box<[u8]>,
}
pub(in crate::step::hanger) struct DecimalCache {
    powers: Vec<Power>,
    bytes: usize,
}
impl DecimalCache {
    pub(in crate::step::hanger) fn new(m: &mut Memory, w: &mut Work<'_>) -> Result<Self, Error> {
        let base = m.used;
        let result = (|| {
            m.add(size_of::<Self>())?;
            let powers = m.allocate(ENTRIES, w)?;
            Ok(Self {
                powers,
                bytes: m.used - base,
            })
        })();
        if result.is_err() {
            m.used = base;
        }
        result
    }
    pub(in crate::step::hanger) fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub(in crate::step::hanger) fn release(self, m: &mut Memory) {
        let bytes = self.bytes;
        drop(self);
        m.used -= bytes;
    }
    pub(in crate::step::hanger) fn format(
        &mut self,
        value: &N,
        m: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<FiniteDecimal, Error> {
        let base = m.used;
        let result = self.encode(value, m, w);
        // Insertion is the final nonfallible mutation. Any earlier failure
        // leaves the cache unchanged and releases this number's own output.
        if result.is_err() {
            m.used = base;
        }
        result
    }
    fn power(&self, exponent: usize, w: &mut Work<'_>) -> Result<(BigInt, Option<usize>), Error> {
        // At most128 sorted descriptors: eight comparisons, no hash allocation.
        w.charge(9)?;
        match self.powers.binary_search_by_key(&exponent, |p| p.exponent) {
            Ok(i) => {
                w.charge(1 + self.powers[i].bytes.len())?;
                Ok((
                    BigInt::from_bytes_le(Sign::Plus, &self.powers[i].bytes),
                    None,
                ))
            }
            Err(i) => {
                let bits = 3 * exponent + 1;
                let words = bits.div_ceil(64);
                // Retain the audited exponentiation bound, ONCE per new power.
                w.charge(16 * (1 + words * words) * (1 + bits.ilog2() as usize))?;
                Ok((BigInt::from(5_u32).pow(exponent as u32), Some(i)))
            }
        }
    }
    fn encode(
        &mut self,
        value: &N,
        m: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<FiniteDecimal, Error> {
        w.charge(1)?;
        let nb = value.0.numer().bits() as usize;
        let db = value.0.denom().bits() as usize;
        if nb.max(db) > 16384 {
            return Err(Error::Limit);
        }
        let bits = nb + 3 * db + 1;
        let scratch = 128 * (256 + bits.div_ceil(8) + db.div_ceil(8));
        m.add(scratch)?;
        let words = db.div_ceil(64);
        w.charge(1 + words * words)?;
        let mut denominator = value.0.denom().clone();
        let twos = denominator.trailing_zeros().unwrap_or(0) as usize;
        denominator >>= twos;
        let mut fives = 0;
        loop {
            let words = (denominator.bits() as usize).div_ceil(64);
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
        let mut power = None;
        w.charge(1 + (nb + scale - twos).div_ceil(64))?;
        let mut integer = value.0.numer().clone();
        integer <<= scale - twos;
        if scale > fives {
            let (p, index) = self.power(scale - fives, w)?;
            let nw = (integer.bits() as usize).div_ceil(64);
            let pw = (p.bits() as usize).div_ceil(64);
            w.charge(4 * (1 + nw * pw))?;
            integer *= &p;
            power = Some((p, index));
        }
        let output_bits = integer.bits() as usize;
        let words = output_bits.div_ceil(64);
        // Pinned num-bigint0.4.8 to_radix_digits_le uses quadratic limb
        // division (or recursively halved chunks). Exponentiation was charged
        // separately above; a cache hit does not repeat that operation.
        w.charge(16 * (1 + words * words))?;
        // log2(10)>3: ceil(bits/3), plus sign/zero, bounds actual digits.
        let output = emit(&integer, scale, output_bits.div_ceil(3) + 1, m, w)?;
        if let Some((power, Some(index))) = power
            && self.powers.len() < ENTRIES
        {
            let bytes = (power.bits() as usize).div_ceil(8);
            let cost = 1 + 2 * bytes + self.powers.len() - index;
            // Cache retention is optional. It must not reject a successful
            // number solely because this extra constant cannot be retained.
            if w.remaining() >= cost && m.add(bytes).is_ok() {
                w.charge(cost)?;
                let (_, data) = power.to_bytes_le();
                self.powers.insert(
                    index,
                    Power {
                        exponent: scale - fives,
                        bytes: data.into_boxed_slice(),
                    },
                );
                self.bytes += bytes;
            }
        }
        m.used -= scratch;
        Ok(output)
    }
}
#[cfg(test)]
mod tests;
