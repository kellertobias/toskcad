//! Exact physical decimal grid for actual decoded tensor poles. A grid integer
//! always means integer /unit millimetres; no millimetre reinterpretation.
use super::*;
use num_bigint::{BigInt, Sign};
use num_traits::Zero;
pub struct PoleInteger {
    sign: i8,
    bytes: Box<[u8]>,
}
impl PoleInteger {
    pub fn signum(&self) -> i8 {
        self.sign
    }
    pub fn numerator_magnitude_le(&self) -> &[u8] {
        &self.bytes
    }
    fn retain(value: BigInt, b: &mut Budget<'_>) -> Result<Self, Error> {
        b.charge(1 + (value.bits() as usize).div_ceil(8))?;
        let (sign, bytes) = value.to_bytes_le();
        b.retain(bytes.len())?;
        Ok(Self {
            sign: match sign {
                Sign::Minus => -1,
                Sign::NoSign => 0,
                Sign::Plus => 1,
            },
            bytes: bytes.into_boxed_slice(),
        })
    }
}
pub(super) struct Grid {
    pub(super) scale: u32,
    pub(super) five_power: u32,
    pub(super) points: BTreeMap<u32, [PoleInteger; 3]>,
    pub(super) unit: PoleInteger,
}
struct Decimal {
    digits: Vec<u8>,
    exponent: i32,
    negative: bool,
}
fn decimal(value: Scalar<'_>, b: &mut Budget<'_>) -> Result<Decimal, Error> {
    b.charge(value.text.len() * 2 + 1)?;
    let text = value.text.strip_prefix('+').unwrap_or(value.text);
    let negative = text.starts_with('-');
    let text = text.strip_prefix('-').unwrap_or(text);
    let (mantissa, exponent) = text.split_once('E').unwrap_or((text, "0"));
    let fraction = mantissa.split_once('.').map_or(0, |(_, v)| v.len());
    let mut digits = mantissa.bytes().filter(|&v| v != b'.').collect::<Vec<_>>();
    let first = digits.iter().position(|&v| v != b'0');
    let Some(first) = first else {
        return Ok(Decimal {
            digits: vec![b'0'],
            exponent: 0,
            negative: false,
        });
    };
    digits.drain(..first);
    let mut exponent =
        i64::from(exponent.parse::<i32>().map_err(|_| Error::Unsupported)?) - fraction as i64;
    while digits.len() > 1 && digits.last() == Some(&b'0') {
        digits.pop();
        exponent += 1;
    }
    if exponent.unsigned_abs() > 16384 || digits.len() > 16384 {
        return Err(Error::Limit);
    }
    Ok(Decimal {
        digits,
        exponent: exponent as i32,
        negative,
    })
}
struct Powers {
    values: BTreeMap<u32, BigInt>,
    bytes: usize,
}
fn power(exponent: u32, powers: &mut Powers, b: &mut Budget<'_>) -> Result<BigInt, Error> {
    b.charge(16 * (1 + powers.values.len().max(1).ilog2() as usize))?;
    if let Some(value) = powers.values.get(&exponent) {
        b.charge(1 + (value.bits() as usize).div_ceil(8))?;
        return Ok(value.clone());
    }
    let bits = exponent as usize * 4 + 1;
    let words = bits.div_ceil(64);
    b.charge(16 * (1 + words * words) * (1 + exponent.max(1).ilog2() as usize))?;
    let out = BigInt::from(5).pow(exponent);
    let bytes = 512 + 8 * (32 + (out.bits() as usize).div_ceil(8));
    b.retain(bytes)?;
    powers.bytes += bytes;
    powers.values.insert(exponent, out.clone());
    Ok(out)
}
impl Grid {
    pub(super) fn prepare(d: &Decoded<'_>, b: &mut Budget<'_>) -> Result<Self, Error> {
        b.charge(d.tensors.len() + 1)?;
        let count = d
            .tensors
            .iter()
            .try_fold(d.planes.len() + d.vertices.len(), |n, s| {
                n.checked_add(s.points.len()).ok_or(Error::Limit)
            })?;
        let tree_cost = 64 * (1 + (count + 1).ilog2() as usize);
        let mut ids = BTreeSet::new();
        for surface in &d.tensors {
            for &id in &surface.points {
                b.charge(tree_cost)?;
                ids.insert(id);
            }
        }
        for plane in &d.planes {
            b.charge(tree_cost)?;
            ids.insert(plane.origin);
        }
        for &id in d.vertices.values() {
            b.charge(tree_cost)?;
            ids.insert(id);
        }
        // Set nodes and map nodes retain different owners; numeric coincidence is
        // never used to alter the STEP vertex/edge/face registries.
        b.retain(192 + ids.len() * 48)?;
        let mut scale = 0_u32;
        for &id in &ids {
            b.charge(tree_cost)?;
            for coordinate in &d.points[&id].coordinates {
                let value = decimal(*coordinate, b)?;
                scale = scale.max(value.exponent.saturating_neg().max(0) as u32);
            }
        }
        if scale > 16384 {
            return Err(Error::Limit);
        }
        let baseline = b.bytes;
        let grid = match Self::build(d, &ids, scale, scale.min(1), b)? {
            Some(grid) => grid,
            None => {
                b.bytes = baseline;
                Self::build(d, &ids, scale, scale, b)?.ok_or(Error::Invalid)?
            }
        };
        b.bytes -= 192 + ids.len() * 48;
        Ok(grid)
    }
    fn build(
        d: &Decoded<'_>,
        ids: &BTreeSet<u32>,
        scale: u32,
        five_power: u32,
        b: &mut Budget<'_>,
    ) -> Result<Option<Self>, Error> {
        use num_traits::Euclid;
        let mut powers = Powers {
            values: BTreeMap::new(),
            bytes: 0,
        };
        let mut unit = power(five_power, &mut powers, b)?;
        b.charge(1 + (unit.bits() as usize + scale as usize).div_ceil(8))?;
        unit <<= scale as usize;
        let mut points = BTreeMap::new();
        b.retain(1024)?;
        let mut cache = ScalarCache::new(b)?;
        let tree_cost = 64 * (1 + (ids.len() + 1).ilog2() as usize);
        for &id in ids {
            b.charge(tree_cost * 2)?;
            let values = &d.points[&id].coordinates;
            if values.len() != 3 {
                return Err(Error::Invalid);
            }
            let mut coordinates = Vec::with_capacity(3);
            for (axis, scalar) in values.iter().enumerate() {
                let mut lookup = if axis < 2 {
                    Some(cache.lookup(scalar.text, b)?)
                } else {
                    None
                };
                if let Some(value) = lookup.as_mut().and_then(|l| l.value.take()) {
                    coordinates.push(value);
                    continue;
                }
                let value = decimal(*scalar, b)?;
                let two = i64::from(scale) + i64::from(value.exponent);
                let five = i64::from(five_power) + i64::from(value.exponent);
                if !(0..=32768).contains(&two) || five.unsigned_abs() > 32768 {
                    return Err(Error::Limit);
                }
                let words = (value.digits.len() * 4).div_ceil(64);
                let chunks = value.digits.len().div_ceil(9);
                b.charge(8 * value.digits.len() + 16 * (1 + chunks * words))?;
                let mut integer = BigInt::parse_bytes(&value.digits, 10).ok_or(Error::Invalid)?;
                if value.negative {
                    integer = -integer;
                }
                let factor = power(five.unsigned_abs() as u32, &mut powers, b)?;
                let a = (integer.bits() as usize).div_ceil(32);
                let c = (factor.bits() as usize).div_ceil(32);
                b.charge(16 * (1 + a * c))?;
                if five >= 0 {
                    integer *= factor;
                } else {
                    let (quotient, remainder) = integer.div_rem_euclid(&factor);
                    if !remainder.is_zero() {
                        return Ok(None);
                    }
                    integer = quotient;
                }
                b.charge(1 + (integer.bits() as usize + two as usize).div_ceil(8))?;
                integer <<= two as usize;
                let integer = PoleInteger::retain(integer, b)?;
                if let Some(lookup) = lookup {
                    cache.insert(lookup, scalar.text, &integer, b)?;
                }
                coordinates.push(integer);
            }
            let array: [PoleInteger; 3] = coordinates.try_into().map_err(|_| Error::Invalid)?;
            b.retain(224)?;
            points.insert(id, array);
        }
        let unit = PoleInteger::retain(unit, b)?;
        b.bytes -= powers.bytes + cache.bytes;
        drop(powers);
        drop(cache);
        if unit.signum() != 1 {
            return Err(Error::Invalid);
        }
        Ok(Some(Self {
            scale,
            five_power,
            points,
            unit,
        }))
    }
}

/// Exact numeric preparation over actual decoded tensor poles. This grants no
/// occupied-material admission. Integers mean numerator /denominator mm.
#[derive(Clone, Copy)]
pub struct GridOptions {
    pub maximum_preparation_work: usize,
    pub maximum_bytes: usize,
}
impl Default for GridOptions {
    fn default() -> Self {
        Self {
            maximum_preparation_work: 224_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
pub struct ExactPoleGrid<'d, 'i> {
    source: &'d Decoded<'i>,
    grid: Grid,
    work: usize,
    bytes: usize,
}
impl<'d, 'i> ExactPoleGrid<'d, 'i> {
    pub fn prepare(
        source: &'d Decoded<'i>,
        options: GridOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        if options.maximum_preparation_work == 0
            || options.maximum_preparation_work > 256_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let mut b = Budget {
            work: 0,
            remaining: options.maximum_preparation_work,
            bytes: source.bytes,
            limit: options.maximum_bytes,
            interrupt,
        };
        b.charge(1)?;
        b.retain(4 * 1024 * 1024 + size_of::<Self>())?;
        let grid = Grid::prepare(source, &mut b)?;
        b.bytes -= 4 * 1024 * 1024;
        Ok(Self {
            source,
            grid,
            work: b.work,
            bytes: b.bytes,
        })
    }
    pub fn source(&self) -> &Decoded<'i> {
        self.source
    }
    /// Own preparation work; source parsing is a separately earned prerequisite.
    pub fn preparation_work_used(&self) -> usize {
        self.work
    }
    pub fn aggregate_known_work_used(&self) -> usize {
        self.source.work + self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub fn denominator(&self) -> &PoleInteger {
        &self.grid.unit
    }
    /// Exact chosen denominator is2^binary_power *5^five_power.
    pub fn denominator_powers(&self) -> [u32; 2] {
        [self.grid.scale, self.grid.five_power]
    }
    pub fn points(&self) -> &BTreeMap<u32, [PoleInteger; 3]> {
        &self.grid.points
    }
}
const SLOTS: usize = 16384;
struct ScalarEntry<'i> {
    hash: u64,
    text: &'i str,
    value: PoleInteger,
}
struct ScalarCache<'i> {
    slots: Vec<Option<ScalarEntry<'i>>>,
    len: usize,
    bytes: usize,
}
struct Lookup {
    hash: u64,
    slot: Option<usize>,
    value: Option<PoleInteger>,
}
impl PoleInteger {
    fn copy(&self, b: &mut Budget<'_>) -> Result<Self, Error> {
        b.charge(1 + self.bytes.len())?;
        b.retain(self.bytes.len())?;
        Ok(Self {
            sign: self.sign,
            bytes: self.bytes.clone(),
        })
    }
}
impl<'i> ScalarCache<'i> {
    fn new(b: &mut Budget<'_>) -> Result<Self, Error> {
        let bytes = SLOTS * size_of::<Option<ScalarEntry<'_>>>();
        b.retain(bytes)?;
        b.charge(SLOTS)?;
        let mut slots = Vec::with_capacity(SLOTS);
        slots.resize_with(SLOTS, || None);
        let actual = slots.capacity() * size_of::<Option<ScalarEntry<'_>>>();
        b.retain(actual - bytes)?;
        Ok(Self {
            slots,
            len: 0,
            bytes: actual,
        })
    }
    fn lookup(&self, text: &str, b: &mut Budget<'_>) -> Result<Lookup, Error> {
        b.charge(4 * text.len() + 1)?;
        let mut hash = 0xcbf29ce484222325_u64;
        for byte in text.bytes() {
            hash = (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3)
        }
        let mut slot = hash as usize & (SLOTS - 1);
        for _ in 0..SLOTS {
            b.charge(4)?;
            let Some(entry) = &self.slots[slot] else {
                return Ok(Lookup {
                    hash,
                    slot: Some(slot),
                    value: None,
                });
            };
            if entry.hash == hash {
                b.charge(1 + text.len() + entry.text.len())?;
                if entry.text == text {
                    return Ok(Lookup {
                        hash,
                        slot: None,
                        value: Some(entry.value.copy(b)?),
                    });
                }
            }
            slot = (slot + 1) & (SLOTS - 1);
        }
        Ok(Lookup {
            hash,
            slot: None,
            value: None,
        })
    }
    fn insert(
        &mut self,
        lookup: Lookup,
        text: &'i str,
        value: &PoleInteger,
        b: &mut Budget<'_>,
    ) -> Result<(), Error> {
        if self.len >= SLOTS / 2 {
            return Ok(());
        }
        let Some(slot) = lookup.slot else {
            return Ok(());
        };
        if b.bytes
            .checked_add(value.bytes.len())
            .is_none_or(|v| v > b.limit)
            || b.remaining < value.bytes.len() + 1
        {
            return Ok(());
        }
        let copy = value.copy(b)?;
        self.bytes += copy.bytes.len();
        self.slots[slot] = Some(ScalarEntry {
            hash: lookup.hash,
            text,
            value: copy,
        });
        self.len += 1;
        Ok(())
    }
}
