use super::*;
use std::mem::size_of;
pub(super) struct Work<'a> {
    pub remaining: usize,
    pub used: usize,
    pub interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    pub fn charge(&mut self, n: usize) -> Result<(), Error> {
        checkpoint(self.interrupt)?;
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
pub(super) struct Memory {
    pub used: usize,
    pub limit: usize,
}
impl Memory {
    pub fn add(&mut self, n: usize) -> Result<(), Error> {
        let next = self.used.checked_add(n).ok_or(Error::Limit)?;
        if next > self.limit {
            return Err(Error::Limit);
        }
        self.used = next;
        Ok(())
    }
    pub fn allocate<T>(&mut self, count: usize, w: &mut Work<'_>) -> Result<Vec<T>, Error> {
        w.charge(count.checked_add(1).ok_or(Error::Limit)?)?;
        let request = count.checked_mul(size_of::<T>()).ok_or(Error::Limit)?;
        if self.used.checked_add(request).ok_or(Error::Limit)? > self.limit {
            return Err(Error::Limit);
        }
        let mut v = Vec::new();
        v.try_reserve_exact(count).map_err(|_| Error::Limit)?;
        self.retain(v)
    }
    /// Account actual allocator capacity before handing an owned buffer on.
    pub fn retain<T>(&mut self, v: Vec<T>) -> Result<Vec<T>, Error> {
        self.add(
            v.capacity()
                .checked_mul(size_of::<T>())
                .ok_or(Error::Limit)?,
        )?;
        Ok(v)
    }
}
pub(super) fn loop_bytes<T>(loops: &Vec<Vec<T>>) -> usize {
    loops.capacity() * size_of::<Vec<T>>()
        + loops
            .iter()
            .map(|v| v.capacity() * size_of::<T>())
            .sum::<usize>()
}
pub(super) fn failure(e: crate::authored_internal_treatment::Error) -> Error {
    use crate::authored_internal_treatment::Error as E;
    match e {
        E::Invalid => Error::Invalid,
        E::Limit => Error::Limit,
        E::Unresolved => Error::Unresolved,
        E::Cancelled => Error::Cancelled,
    }
}
