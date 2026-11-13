//! Query-local ownership table. It retains compact exact-coordinate bytes, not
//! arbitrary-capacity BigInt allocations. No coordinate tolerance is a key.
use super::*;
use std::mem::size_of;
const MAX_OWNERS: usize = 131_072;
const BYTE_LIMIT: usize = 64 * 1024 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Constraint {
    Rail(i64),
    Row(u32),
    Station(u8, i64),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Key(pub [Constraint; 2]);
impl Key {
    fn hash(self) -> usize {
        let mut hash = 0xcbf29ce484222325_u64;
        for c in self.0 {
            let (tag, value) = match c {
                Constraint::Rail(v) => (0, v as u64),
                Constraint::Row(v) => (1, u64::from(v)),
                Constraint::Station(v, owner) => (2 + u64::from(v), owner as u64),
            };
            for byte in [tag, value].into_iter().flat_map(u64::to_le_bytes) {
                hash = (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3);
            }
        }
        hash as usize
    }
}
#[derive(Clone, Copy)]
struct Slot {
    key: Key,
    index: usize,
}
struct Owner {
    coordinates: Box<[u8]>,
    vertex: u32,
}
pub(super) struct Owners {
    slots: Vec<Option<Slot>>,
    owners: Vec<Owner>,
    bytes: usize,
    maximum: usize,
}
impl Owners {
    pub fn new(maximum: usize, w: &mut Work<'_>) -> Result<Self, Error> {
        w.step()?;
        let maximum = maximum.min(MAX_OWNERS);
        if maximum == 0 {
            return Err(Error::Limit);
        }
        let count = maximum
            .checked_mul(2)
            .and_then(usize::checked_next_power_of_two)
            .ok_or(Error::Limit)?;
        w.charge(count + maximum)?;
        let mut slots = Vec::new();
        let mut owners = Vec::new();
        slots.try_reserve_exact(count).map_err(|_| Error::Limit)?;
        owners
            .try_reserve_exact(maximum)
            .map_err(|_| Error::Limit)?;
        let allocations = [
            (slots.capacity(), size_of::<Option<Slot>>()),
            (owners.capacity(), size_of::<Owner>()),
        ];
        let bytes = storage::bytes(&allocations, size_of::<Self>(), BYTE_LIMIT)?;
        slots.resize(count, None);
        Ok(Self {
            slots,
            owners,
            bytes,
            maximum,
        })
    }
    /// Coordinates have already undergone the SAME symbolic periodic transfer
    /// as key. A reused key with different exact source data is an error.
    pub fn vertex(
        &mut self,
        key: Key,
        coordinates: Box<[u8]>,
        proposed: u32,
        w: &mut Work<'_>,
    ) -> Result<u32, Error> {
        w.step()?;
        let mask = self.slots.len() - 1;
        let mut at = key.hash() & mask;
        // At most half the fixed table can be occupied. Probing and coordinate
        // comparison are nevertheless charged, including every cache hit.
        loop {
            w.step()?;
            match self.slots[at] {
                Some(slot) if slot.key == key => {
                    let owner = &self.owners[slot.index];
                    w.charge(1 + coordinates.len().max(owner.coordinates.len()).div_ceil(8))?;
                    return if coordinates == owner.coordinates {
                        Ok(owner.vertex)
                    } else {
                        Err(Error::Unresolved)
                    };
                }
                Some(_) => at = (at + 1) & mask,
                None => break,
            }
        }
        if self.owners.len() == self.maximum {
            return Err(Error::Limit);
        }
        let next = self
            .bytes
            .checked_add(coordinates.len())
            .ok_or(Error::Limit)?;
        if next > BYTE_LIMIT {
            return Err(Error::Limit);
        }
        w.step()?;
        let index = self.owners.len();
        self.owners.push(Owner {
            coordinates,
            vertex: proposed,
        });
        self.slots[at] = Some(Slot { key, index });
        self.bytes = next;
        Ok(proposed)
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
    fn exact_ownership_is_bounded_and_failed_reuse_preserves_original() {
        let w = &mut work();
        let mut owners = Owners::new(2, w).unwrap();
        let key = Key([Constraint::Rail(-2), Constraint::Station(0, 0)]);
        assert_eq!(
            owners
                .vertex(key, vec![1, 2, 3].into_boxed_slice(), 7, w)
                .unwrap(),
            7
        );
        assert_eq!(
            owners
                .vertex(key, vec![1, 2, 3].into_boxed_slice(), 9, w)
                .unwrap(),
            7
        );
        assert_eq!(
            owners.vertex(key, vec![1, 2, 4].into_boxed_slice(), 9, w),
            Err(Error::Unresolved)
        );
        assert_eq!(
            owners
                .vertex(key, vec![1, 2, 3].into_boxed_slice(), 9, w)
                .unwrap(),
            7
        );
        let second = Key([Constraint::Rail(3), Constraint::Row(1)]);
        assert_eq!(
            owners
                .vertex(second, vec![4].into_boxed_slice(), 8, w)
                .unwrap(),
            8
        );
        assert_eq!(
            owners.vertex(
                Key([Constraint::Row(2), Constraint::Station(1, 0)]),
                Box::new([]),
                9,
                w
            ),
            Err(Error::Limit)
        );
        let mut cancelled = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: Some(&|| true),
        };
        assert_eq!(
            owners.vertex(key, vec![1, 2, 3].into_boxed_slice(), 9, &mut cancelled),
            Err(Error::Cancelled)
        );
        assert_eq!(owners.owners.len(), 2);
        assert!(owners.bytes <= BYTE_LIMIT);
    }
}
#[cfg(test)]
mod resource_tests {
    use super::*;
    #[test]
    fn retained_exact_bytes_hit_declared_memory_budget_before_owner_count_and_preserve_hits() {
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let mut owners = Owners::new(10_000, w).unwrap();
        let first = Key([Constraint::Rail(0), Constraint::Row(0)]);
        // The generic storage guard receives compact coordinate bytes; its
        // budget must hold even at the maximum bounded rational encoding size.
        let bytes = vec![1; 8192].into_boxed_slice();
        let mut inserted = 0;
        loop {
            let key = Key([Constraint::Rail(inserted), Constraint::Row(0)]);
            match owners.vertex(key, bytes.clone(), inserted as u32, w) {
                Err(Error::Limit) => break,
                result => {
                    assert_eq!(result, Ok(inserted as u32));
                    inserted += 1;
                }
            }
        }
        assert!(inserted > 1 && inserted < 10_000);
        assert!(owners.bytes <= BYTE_LIMIT);
        assert_eq!(owners.vertex(first, bytes, 99, w).unwrap(), 0);
        assert!(matches!(Owners::new(0, w), Err(Error::Limit)));
    }
}
