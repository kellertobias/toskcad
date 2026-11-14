//! Independent literal-alias registry over a real retained point inventory.
use super::*;
pub(super) struct Registry<'s> {
    source: &'s Boundary,
    slots: Vec<Option<(usize, usize)>>,
}
impl<'s> Registry<'s> {
    pub(super) fn prepare(
        source: &'s Boundary,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        let count = [Side::Outer, Side::Inner]
            .into_iter()
            .map(|side| source.source_boundary(side).points().len())
            .sum();
        let mut slots = memory.allocate(count, w)?;
        slots.resize(count, None);
        Ok(Self { source, slots })
    }
    pub(super) fn bind(
        &mut self,
        vertex: BoundaryVertexKey,
        row: &CoupledRow,
        index: usize,
        endpoint: usize,
        previous: &[CoupledRow],
        w: &mut Work<'_>,
    ) -> Result<(), Error> {
        w.charge(1)?;
        if self.source.point(vertex).is_none() {
            return Err(Error::Invalid);
        }
        let slot = vertex_index(self.source, vertex);
        if let Some((prior, end)) = self.slots[slot] {
            let earlier = if prior == index {
                row
            } else {
                &previous[prior]
            };
            same_point(
                &earlier.coefficients[end * 11],
                &row.coefficients[endpoint * 11],
                w,
            )?;
        } else {
            self.slots[slot] = Some((index, endpoint));
        }
        Ok(())
    }
    pub(super) fn complete(&self, w: &mut Work<'_>) -> Result<(), Error> {
        w.charge(self.slots.len() + 1)?;
        if self.slots.iter().any(Option::is_none) {
            Err(Error::Unresolved)
        } else {
            Ok(())
        }
    }
    pub(super) fn release(&self, memory: &mut Memory) {
        memory.release(&self.slots);
    }
}
