//! Numeric-only literal reuse. Topology owners are never keys in this cache.
use super::*;
use std::{cell::Cell, cmp::Ordering, collections::BTreeMap, rc::Rc};
const ENTRIES: usize = 2048;
#[derive(Eq)]
struct Key {
    sign: i8,
    n: Box<[u8]>,
    d: Box<[u8]>,
    comparisons: Rc<Cell<usize>>,
}
impl PartialEq for Key {
    fn eq(&self, other: &Self) -> bool {
        self.sign == other.sign && self.n == other.n && self.d == other.d
    }
}
fn bytes(a: &[u8], b: &[u8], count: &Cell<usize>) -> Ordering {
    for (a, b) in a.iter().zip(b) {
        count.set(count.get() + 1);
        let order = a.cmp(b);
        if order != Ordering::Equal {
            return order;
        }
    }
    count.set(count.get() + 1);
    a.len().cmp(&b.len())
}
impl Ord for Key {
    fn cmp(&self, b: &Self) -> Ordering {
        self.comparisons.set(self.comparisons.get() + 1);
        self.sign
            .cmp(&b.sign)
            .then_with(|| bytes(&self.n, &b.n, &self.comparisons))
            .then_with(|| bytes(&self.d, &b.d, &self.comparisons))
    }
}
impl PartialOrd for Key {
    fn partial_cmp(&self, b: &Self) -> Option<Ordering> {
        Some(self.cmp(b))
    }
}
pub(super) struct Literals {
    map: BTreeMap<Key, Rc<FiniteDecimal>>,
    comparisons: Rc<Cell<usize>>,
}
impl Literals {
    pub(super) fn new(m: &mut Memory) -> Result<Self, Error> {
        m.add(size_of::<Cell<usize>>() + 2 * size_of::<usize>())?;
        Ok(Self {
            map: BTreeMap::new(),
            comparisons: Rc::new(Cell::new(0)),
        })
    }
    pub(super) fn format(
        &mut self,
        value: &N,
        cache: &mut super::super::decimal::DecimalCache,
        m: &mut Memory,
        w: &mut Work<'_>,
        peak: &mut usize,
    ) -> Result<Rc<FiniteDecimal>, Error> {
        let bits = value.0.numer().bits().max(value.0.denom().bits());
        if bits > 16384 {
            return Err(Error::Limit);
        }
        w.charge(2 + (value.0.numer().bits() + value.0.denom().bits()).div_ceil(8) as usize)?;
        let (sign, n) = value.0.numer().to_bytes_le();
        let (_, d) = value.0.denom().to_bytes_le();
        let key = Key {
            sign: match sign {
                Sign::Minus => -1,
                Sign::NoSign => 0,
                Sign::Plus => 1,
            },
            n: n.into_boxed_slice(),
            d: d.into_boxed_slice(),
            comparisons: self.comparisons.clone(),
        };
        self.comparisons.set(0);
        let hit = self.map.get(&key).cloned();
        w.charge(1 + self.comparisons.get())?;
        if let Some(text) = hit {
            return Ok(text);
        }
        let literal = cache.format(value, m, w)?;
        let nb = value.0.numer().bits() as usize;
        let db = value.0.denom().bits() as usize;
        let scratch = 128 * (256 + (nb + 3 * db + 1).div_ceil(8) + db.div_ceil(8));
        *peak = (*peak).max(m.used + scratch);
        // Rc has two usize counters in addition to the already charged value.
        m.add(2 * size_of::<usize>())?;
        let literal = Rc::new(literal);
        if self.map.len() < ENTRIES {
            // Pinned std BTree B=6 has <=11 key/value slots and12 child slots.
            // This Key+Rc value needs <768 bytes per node. <=one node per entry.
            let retain = 768 + key.n.len() + key.d.len();
            if m.add(retain).is_ok() {
                self.comparisons.set(0);
                self.map.insert(key, literal.clone());
                w.charge(1 + self.comparisons.get())?;
            }
        }
        *peak = (*peak).max(m.used);
        Ok(literal)
    }
}
pub(super) fn release(values: Vec<Rc<FiniteDecimal>>, m: &mut Memory) {
    for value in values {
        if Rc::strong_count(&value) == 1 {
            m.used -= value.retained_bytes() + 2 * size_of::<usize>();
        }
    }
}
// CARTESIAN_POINT is a geometric definition, never a VERTEX_POINT owner.
// Reusing an equal geometric definition preserves every distinct topology ID.
pub(super) struct Points {
    map: BTreeMap<Vec<Key>, Entity>,
    comparisons: Rc<Cell<usize>>,
}
impl Points {
    pub(super) fn new(m: &mut Memory) -> Result<Self, Error> {
        m.add(size_of::<Cell<usize>>() + 2 * size_of::<usize>())?;
        Ok(Self {
            map: BTreeMap::new(),
            comparisons: Rc::new(Cell::new(0)),
        })
    }
    pub(super) fn lookup(&self, values: &[N], w: &mut Work<'_>) -> Result<PointLookup, Error> {
        let mut keys = Vec::with_capacity(values.len());
        for value in values {
            let bits = value.0.numer().bits().max(value.0.denom().bits());
            if bits > 16384 {
                return Err(Error::Limit);
            }
            w.charge(2 + (value.0.numer().bits() + value.0.denom().bits()).div_ceil(8) as usize)?;
            let (sign, n) = value.0.numer().to_bytes_le();
            let (_, d) = value.0.denom().to_bytes_le();
            keys.push(Key {
                sign: match sign {
                    Sign::Minus => -1,
                    Sign::NoSign => 0,
                    Sign::Plus => 1,
                },
                n: n.into_boxed_slice(),
                d: d.into_boxed_slice(),
                comparisons: self.comparisons.clone(),
            });
        }
        self.comparisons.set(0);
        let entity = self.map.get(&keys).copied();
        w.charge(1 + self.comparisons.get())?;
        Ok(PointLookup { keys, entity })
    }
    pub(super) fn insert(
        &mut self,
        lookup: PointLookup,
        entity: Entity,
        m: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<(), Error> {
        if self.map.len() >= 8192 {
            return Ok(());
        }
        // Leaf:16-byte header+11*(Vec24+Entity8)=368; internal adds96=464.
        let retain = 512
            + lookup.keys.capacity() * size_of::<Key>()
            + lookup
                .keys
                .iter()
                .map(|k| k.n.len() + k.d.len())
                .sum::<usize>();
        if m.add(retain).is_ok() {
            self.comparisons.set(0);
            self.map.insert(lookup.keys, entity);
            w.charge(1 + self.comparisons.get())?;
        }
        Ok(())
    }
}
pub(super) struct PointLookup {
    keys: Vec<Key>,
    pub(super) entity: Option<Entity>,
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn numeric_cache_preserves_exact_values_dimensions_and_key_order_without_topology_ownership() {
        let mut m = Memory::new(0, 64 * 1024 * 1024).unwrap();
        let mut w = Work::new(16_000_000, None).unwrap();
        let mut cache = super::super::super::decimal::DecimalCache::new(&mut m, &mut w).unwrap();
        let mut literals = Literals::new(&mut m).unwrap();
        let mut peak = m.used;
        let values = [
            N::integer(-7),
            N::integer(0),
            N::integer(7),
            N(BigRational::new(7.into(), 2.into())),
            N(BigRational::new(14.into(), 4.into())),
        ];
        let text = values
            .iter()
            .map(|q| {
                literals
                    .format(q, &mut cache, &mut m, &mut w, &mut peak)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        assert!(Rc::ptr_eq(&text[3], &text[4]));
        assert!(!Rc::ptr_eq(&text[0], &text[2]));
        for (value, text) in values.iter().zip(&text) {
            let (mantissa, exponent) = text
                .as_str()
                .split_once("E-")
                .unwrap_or((text.as_str(), "0"));
            let integer =
                num_bigint::BigInt::parse_bytes(mantissa.strip_suffix('.').unwrap().as_bytes(), 10)
                    .unwrap();
            let decimal = BigRational::new(
                integer,
                num_bigint::BigInt::from(10).pow(exponent.parse().unwrap()),
            );
            assert_eq!(decimal, value.0);
        }
        for (a, _) in &literals.map {
            for (b, _) in &literals.map {
                assert_eq!(a.partial_cmp(b), Some(a.cmp(b)));
                assert_eq!(a == b, a.cmp(b) == Ordering::Equal);
            }
        }
        let mut points = Points::new(&mut m).unwrap();
        let xy = [N::integer(1), N::integer(2)];
        let first = points.lookup(&xy, &mut w).unwrap();
        assert!(first.entity.is_none());
        points.insert(first, Entity(19), &mut m, &mut w).unwrap();
        assert_eq!(
            points.lookup(&xy, &mut w).unwrap().entity.map(|e| e.0),
            Some(19)
        );
        assert!(
            points
                .lookup(&[N::integer(1), N::integer(2), N::integer(0)], &mut w)
                .unwrap()
                .entity
                .is_none()
        );
        let prior = m.used;
        release(text, &mut m);
        assert_eq!(m.used, prior); // cached Rc values remain owned.
    }
}
