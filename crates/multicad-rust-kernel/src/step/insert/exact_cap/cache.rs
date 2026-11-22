//! Exact coordinates cached by immutable original POINT leaf identity.
use super::*;
#[derive(Clone, Copy)]
pub(crate) struct LiteralPoint<'a> {
    pub index: usize,
    pub coordinates: [&'a str; 3],
}
struct Entry<'a> {
    source: [&'a str; 3],
    value: [N; 3],
}
pub(crate) struct DecimalMemo<'a> {
    entries: Vec<Option<Entry<'a>>>,
    payload: usize,
}
impl<'a> DecimalMemo<'a> {
    pub(crate) fn new(
        maximum: usize,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        let mut entries = memory.allocate(maximum, w)?;
        entries.resize_with(maximum, || None);
        Ok(Self {
            entries,
            payload: 0,
        })
    }
    pub(super) fn point(
        &mut self,
        point: LiteralPoint<'a>,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<[N; 3], Error> {
        w.charge(4)?;
        let slot = self.entries.get_mut(point.index).ok_or(Error::Invalid)?;
        if let Some(old) = slot {
            // Same arena/index means the SAME immutable original leaves.
            if !old
                .source
                .iter()
                .zip(point.coordinates)
                .all(|(a, b)| std::ptr::eq(*a, b))
            {
                return Err(Error::Invalid);
            }
            let words = old
                .value
                .iter()
                .map(|v| ((v.0.numer().bits() + v.0.denom().bits()) as usize).div_ceil(32))
                .sum::<usize>();
            w.charge(1 + words)?;
            return Ok(old.value.clone());
        }
        // Caller retains the4MiB temporary arena while this query runs.
        let value = [
            decimal(point.coordinates[0], w)?,
            decimal(point.coordinates[1], w)?,
            decimal(point.coordinates[2], w)?,
        ];
        // num-bigint0.4.8 BigDigits::clone copies its normalized slice into
        // Inline0/1 or slice::to_vec. Fresh clones retain this conservative
        //32/64-bit limb bound, not arithmetic temporary capacities.
        let bytes = value
            .iter()
            .flat_map(|v| [v.0.numer().bits(), v.0.denom().bits()])
            .map(|bits| (bits as usize).div_ceil(32) * 4 + 8)
            .sum::<usize>();
        memory.add(bytes)?;
        w.charge(1 + bytes.div_ceil(4))?;
        *slot = Some(Entry {
            source: point.coordinates,
            value: value.clone(),
        });
        self.payload += bytes;
        Ok(value)
    }
    pub(crate) fn release(self, memory: &mut Memory) {
        memory.release(&self.entries);
        memory.used -= self.payload;
    }
}
