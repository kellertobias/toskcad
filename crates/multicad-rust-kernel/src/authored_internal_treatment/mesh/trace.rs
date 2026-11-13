//! Exact angular ownership for actual source cap boundaries. Compact immutable
//! bytes share the existing64MiB footprint budget with the request point cache.
use super::{exact::Q, *};
use std::cmp::Ordering;
const BYTES: usize = 64 * 1024 * 1024;
struct Angle {
    token: Box<[u8]>,
    bounds: I,
}
pub(super) struct Trace {
    keys: Vec<Option<Angle>>,
    bytes: usize,
    origin: f64,
}
impl Trace {
    #[cfg(test)]
    pub fn new(maximum: usize, w: &mut Work<'_>) -> Result<Self, Error> {
        Self::at_origin(maximum, 0., w)
    }
    pub fn at_origin(maximum: usize, origin: f64, w: &mut Work<'_>) -> Result<Self, Error> {
        w.step()?;
        if maximum > 131_072 {
            return Err(Error::Limit);
        }
        w.charge(maximum)?;
        let mut keys = Vec::new();
        keys.try_reserve_exact(maximum).map_err(|_| Error::Limit)?;
        let allocations = [(keys.capacity(), std::mem::size_of::<Option<Angle>>())];
        let bytes = storage::bytes(&allocations, std::mem::size_of::<Self>(), BYTES)?;
        keys.resize_with(maximum, || None);
        Ok(Self {
            keys,
            bytes,
            origin,
        })
    }
    pub fn origin(&self) -> f64 {
        self.origin
    }
    pub fn end(&self) -> f64 {
        1. + self.origin
    }
    pub fn memory_bytes(&self) -> usize {
        self.bytes
    }
    pub fn check_memory(&self, owners: usize) -> Result<(), Error> {
        storage::bytes(&[(owners, 1)], self.bytes, BYTES).map(|_| ())
    }
    pub fn remember(
        &mut self,
        id: usize,
        v: &Q,
        owners: usize,
        w: &mut Work<'_>,
    ) -> Result<(), Error> {
        w.step()?;
        if id >= self.keys.len() {
            return Err(Error::Limit);
        }
        let end = if self.origin == 0. {
            Q::integer(1, w)?
        } else {
            Q::float(self.end(), w)?
        };
        let canonical = if v.cmp(&end, w)? == Ordering::Equal {
            if self.origin == 0. {
                Q::integer(0, w)?
            } else {
                Q::float(self.origin, w)?
            }
        } else {
            v.clone()
        };
        let token = Q::encode(&[&canonical], w)?;
        if let Some(old) = &self.keys[id] {
            w.charge(1 + old.token.len().max(token.len()).div_ceil(8))?;
            return if token == old.token {
                Ok(())
            } else {
                Err(Error::Unresolved)
            };
        }
        let next = self.bytes.checked_add(token.len()).ok_or(Error::Limit)?;
        storage::bytes(&[(owners, 1)], next, BYTES)?;
        let bounds = canonical.enclosure(w)?.1;
        w.step()?;
        self.keys[id] = Some(Angle { token, bounds });
        self.bytes = next;
        Ok(())
    }
    pub fn bounds(&self, id: usize, w: &mut Work<'_>) -> Result<I, Error> {
        w.step()?;
        self.keys
            .get(id)
            .and_then(Option::as_ref)
            .map(|a| a.bounds)
            .ok_or(Error::Unresolved)
    }
    pub fn parameter(&self, id: usize, work: &mut Work<'_>) -> Result<Q, Error> {
        work.step()?;
        let value = self
            .keys
            .get(id)
            .and_then(Option::as_ref)
            .ok_or(Error::Unresolved)?;
        Q::decode(&value.token, work)
    }
    pub fn compare(&self, a: usize, b: usize, w: &mut Work<'_>) -> Result<Ordering, Error> {
        w.step()?;
        let [a, b] = [a, b].map(|id| {
            self.keys
                .get(id)
                .and_then(Option::as_ref)
                .ok_or(Error::Unresolved)
        });
        let a = a?;
        let b = b?;
        if a.bounds.hi < b.bounds.lo {
            return Ok(Ordering::Less);
        }
        if b.bounds.hi < a.bounds.lo {
            return Ok(Ordering::Greater);
        }
        w.charge(1 + a.token.len().max(b.token.len()).div_ceil(8))?;
        if a.token == b.token {
            return Ok(Ordering::Equal);
        }
        Q::decode(&a.token, w)?.cmp(&Q::decode(&b.token, w)?, w)
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
    fn canonical_angular_ids_preserve_positive_subulp_separation_and_failed_receipts() {
        let w = &mut work();
        let mut trace = Trace::new(4, w).unwrap();
        let zero = Q::integer(0, w).unwrap();
        let one = Q::integer(1, w).unwrap();
        trace.remember(0, &one, 0, w).unwrap();
        trace.remember(0, &zero, 0, w).unwrap();
        let half = Q::float(0.5, w).unwrap();
        let quarter = Q::float(0.75, w).unwrap();
        let above = half.add(&Q::float(1e-20, w).unwrap(), w).unwrap();
        trace.remember(1, &half, 0, w).unwrap();
        trace.remember(2, &quarter, 0, w).unwrap();
        trace.remember(3, &above, 0, w).unwrap();
        assert_eq!(trace.compare(0, 1, w).unwrap(), Ordering::Less);
        assert_eq!(trace.compare(2, 1, w).unwrap(), Ordering::Greater);
        assert_eq!(trace.compare(1, 1, w).unwrap(), Ordering::Equal);
        assert_eq!(trace.compare(1, 3, w).unwrap(), Ordering::Less);
        assert_eq!(trace.compare(3, 1, w).unwrap(), Ordering::Greater);
        assert_eq!(trace.remember(0, &half, 0, w), Err(Error::Unresolved));
        assert_eq!(trace.compare(0, 0, w).unwrap(), Ordering::Equal);
        assert!(matches!(trace.bounds(4, w), Err(Error::Unresolved)));
        assert_eq!(trace.remember(4, &half, 0, w), Err(Error::Limit));
        assert!(matches!(Trace::new(131_073, w), Err(Error::Limit)));
    }
    #[test]
    fn combined_exact_byte_budget_and_cancel_before_reuse_are_enforced() {
        let w = &mut work();
        let mut trace = Trace::new(2, w).unwrap();
        let zero = Q::integer(0, w).unwrap();
        let initial = trace.memory_bytes();
        assert_eq!(
            trace.remember(0, &zero, BYTES - initial, w),
            Err(Error::Limit)
        );
        assert_eq!(trace.memory_bytes(), initial);
        trace.remember(0, &zero, 0, w).unwrap();
        assert_eq!(trace.check_memory(BYTES), Err(Error::Limit));
        assert!(trace.check_memory(0).is_ok());
        let mut cancelled = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: Some(&|| true),
        };
        assert_eq!(
            trace.remember(0, &zero, 0, &mut cancelled),
            Err(Error::Cancelled)
        );
        assert!(matches!(
            trace.bounds(0, &mut cancelled),
            Err(Error::Cancelled)
        ));
        assert_eq!(trace.compare(0, 0, &mut cancelled), Err(Error::Cancelled));
        assert_eq!(trace.compare(0, 1, w), Err(Error::Unresolved));
        assert_eq!(trace.compare(1, 0, w), Err(Error::Unresolved));
    }
}
