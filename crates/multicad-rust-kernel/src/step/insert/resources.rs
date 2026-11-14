use super::*;
pub(super) struct Work<'a> {
    pub used: usize,
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl<'a> Work<'a> {
    pub fn new(limit: usize, interrupt: Option<&'a dyn Fn() -> bool>) -> Result<Self, Error> {
        if limit == 0 || limit > 16_000_000 {
            return Err(Error::Invalid);
        }
        Ok(Self {
            used: 0,
            remaining: limit,
            interrupt,
        })
    }
    pub fn remaining(&self) -> usize {
        self.remaining
    }
    pub fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
pub(super) struct Memory {
    pub used: usize,
    limit: usize,
}
impl Memory {
    pub fn new(used: usize, limit: usize) -> Result<Self, Error> {
        if used > limit {
            Err(Error::Limit)
        } else {
            Ok(Self { used, limit })
        }
    }
    pub fn add(&mut self, n: usize) -> Result<(), Error> {
        let next = self.used.checked_add(n).ok_or(Error::Limit)?;
        if next > self.limit {
            return Err(Error::Limit);
        }
        self.used = next;
        Ok(())
    }
    pub fn allocate<T>(&mut self, n: usize, w: &mut Work<'_>) -> Result<Vec<T>, Error> {
        w.charge(n.checked_add(1).ok_or(Error::Limit)?)?;
        let request = n.checked_mul(size_of::<T>()).ok_or(Error::Limit)?;
        if self.used.checked_add(request).ok_or(Error::Limit)? > self.limit {
            return Err(Error::Limit);
        }
        let mut out = Vec::new();
        out.try_reserve_exact(n).map_err(|_| Error::Limit)?;
        self.account(out)
    }
    // Admission uses ACTUAL retained capacity, including allocator overcapacity.
    // Keep this final check independently usable for a completed allocation.
    pub(super) fn account<T>(&mut self, out: Vec<T>) -> Result<Vec<T>, Error> {
        let bytes = out
            .capacity()
            .checked_mul(size_of::<T>())
            .ok_or(Error::Limit)?;
        self.add(bytes)?;
        Ok(out)
    }
    pub fn release<T>(&mut self, v: &Vec<T>) {
        self.used -= v.capacity() * size_of::<T>();
    }
}
