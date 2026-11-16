//! Original component spans survive neutral-arena destruction.
use super::*;
use crate::step_read::literal_member_binding::RetainedRecords;
use std::{mem::size_of, ops::Range};
struct Row {
    id: u32,
    first: usize,
    count: usize,
}
pub(crate) struct Records {
    input: Arc<str>,
    rows: Vec<Row>,
    components: Vec<Range<usize>>,
    pub(crate) member: u32,
}
pub(super) fn storage(rows: usize, components: usize) -> GeometryResult<usize> {
    rows.checked_mul(size_of::<Row>())
        .and_then(|n| n.checked_add(components.checked_mul(size_of::<Range<usize>>())?))
        .ok_or(GeometryError::Limit)
}
impl Records {
    pub(super) fn from_actual(
        input: &str,
        entities: &BTreeMap<u32, Vec<part21::Entity<'_>>>,
        member: u32,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> GeometryResult<Self> {
        let total = entities
            .values()
            .try_fold(0usize, |n, e| n.checked_add(e.len()))
            .ok_or(GeometryError::Limit)?;
        let mut rows = Vec::new();
        let mut components = Vec::new();
        rows.try_reserve_exact(entities.len())
            .map_err(|_| GeometryError::Limit)?;
        components
            .try_reserve_exact(total)
            .map_err(|_| GeometryError::Limit)?;
        for (&id, values) in entities {
            checkpoint(interrupt)?;
            let first = components.len();
            for entity in values {
                checkpoint(interrupt)?;
                let source = entity.source.ok_or(GeometryError::Invalid)?;
                let start = (source.as_ptr() as usize)
                    .checked_sub(input.as_ptr() as usize)
                    .ok_or(GeometryError::Invalid)?;
                let end = start
                    .checked_add(source.len())
                    .ok_or(GeometryError::Limit)?;
                if input.get(start..end) != Some(source) {
                    return Err(GeometryError::Invalid);
                }
                components.push(start..end);
            }
            rows.push(Row {
                id,
                first,
                count: values.len(),
            });
        }
        checkpoint(interrupt)?;
        Ok(Self {
            input: Arc::from(input),
            rows,
            components,
            member,
        })
    }
    pub(super) fn input(&self) -> &str {
        &self.input
    }
    pub(super) fn bytes(&self) -> GeometryResult<usize> {
        storage(self.rows.capacity(), self.components.capacity())?
            .checked_add(self.input.len())
            .and_then(|n| n.checked_add(32 + size_of::<Self>()))
            .ok_or(GeometryError::Limit)
    }
}
impl RetainedRecords for Records {
    fn record_count(&self) -> usize {
        self.rows.len()
    }
    fn record_index(&self, id: u32) -> Option<usize> {
        self.rows.binary_search_by_key(&id, |row| row.id).ok()
    }
    fn component(&self, record: usize, index: usize) -> Option<&str> {
        let row = self.rows.get(record)?;
        if index >= row.count {
            return None;
        }
        self.input
            .get(self.components.get(row.first + index)?.clone())
    }
}
