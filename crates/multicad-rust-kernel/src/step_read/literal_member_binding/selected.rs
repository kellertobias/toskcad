//! Immutable selected closure from the original bounded neutral AST. This
//! owns indices/original component slices only; it confers no geometry/material
//! authority. Original input lifetime is retained; the neutral AST may be dropped.
use super::*;
use crate::step_read::part21::Value;
pub(crate) struct SelectedRecords<'input> {
    ids: Vec<u32>,
    ranges: Vec<std::ops::Range<usize>>,
    components: Vec<&'input str>,
    work: usize,
    retained: usize,
    peak: usize,
}
impl RetainedRecords for SelectedRecords<'_> {
    fn record_count(&self) -> usize {
        self.ids.len()
    }
    fn record_index(&self, id: u32) -> Option<usize> {
        self.ids.binary_search(&id).ok()
    }
    fn component(&self, record: usize, index: usize) -> Option<&str> {
        self.components
            .get(self.ranges.get(record)?.clone().nth(index)?)
            .copied()
    }
}
impl SelectedRecords<'_> {
    pub(crate) fn work_used(&self) -> usize {
        self.work
    }
    pub(crate) fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
    pub(crate) fn query_peak_bytes(&self) -> usize {
        self.peak
    }
}
struct Walk<'a> {
    options: Options,
    interrupt: Option<&'a dyn Fn() -> bool>,
    work: usize,
    ids: Vec<u32>,
    visited: Vec<u8>,
    pending: Vec<usize>,
}
impl Walk<'_> {
    fn charge(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.work = self.work.checked_add(n).ok_or(Error::WorkLimit)?;
        if self.work > self.options.maximum_work {
            return Err(Error::WorkLimit);
        }
        Ok(())
    }
    fn reference(&mut self, id: u32) -> Result<()> {
        // Sorted original IDs, with a maximum of2^20 entries. This uses the
        // existing32-unit lookup convention; no string search finds records.
        self.charge(32)?;
        let index = self.ids.binary_search(&id).map_err(|_| Error::Invalid)?;
        if self.visited[index] == 0 {
            self.charge(1)?;
            self.visited[index] = 1;
            self.pending.push(index);
        }
        Ok(())
    }
    fn values(&mut self, values: &[Value<'_>]) -> Result<()> {
        for value in values {
            self.charge(1)?;
            match value {
                Value::Reference(id) => self.reference(*id)?,
                Value::List(values) => self.values(values)?,
                Value::Typed(entity) => self.values(&entity.args)?,
                _ => {}
            }
        }
        Ok(())
    }
}
/// Root must be an actual represented member of this neutral document. The
/// returned references retain the exact selected original components and IDs;
/// unrelated members and metadata are never copied into its dependency set.
/// Document input/AST storage/work are separate resident scopes. This reports
/// the additional index/reference buffers, with exact reservation preflights.
pub(crate) fn selected_records<'input>(
    document: &CompoundDocument<'input>,
    root: u32,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<SelectedRecords<'input>> {
    if options.maximum_work == 0
        || options.maximum_work > 512_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
        || options.maximum_entities == 0
        || options.maximum_entities > 1_048_576
    {
        return Err(Error::Invalid);
    }
    if interrupt.is_some_and(|f| f()) {
        return Err(Error::Cancelled);
    }
    let count = document.records().len();
    if count > options.maximum_entities {
        return Err(Error::EntityLimit);
    }
    let initial = 4 * count + document.represented_members().len() + 1;
    if initial > options.maximum_work {
        return Err(Error::WorkLimit);
    }
    let total_components = document
        .records()
        .values()
        .try_fold(0usize, |n, r| n.checked_add(r.len()))
        .ok_or(Error::StorageLimit)?;
    let mut bytes = 16 * 1024 + std::mem::size_of::<SelectedRecords<'_>>();
    let mut ids = allocate(count, &mut bytes, options.maximum_bytes)?;
    let mut visited = allocate(count, &mut bytes, options.maximum_bytes)?;
    let pending = allocate(count, &mut bytes, options.maximum_bytes)?;
    // Admit the final reference array while traversal scratch is still live.
    let mut ranges = allocate(count, &mut bytes, options.maximum_bytes)?;
    let mut components = allocate(total_components, &mut bytes, options.maximum_bytes)?;
    visited.resize(count, 0);
    ids.extend(document.records().keys().copied());
    let mut walk = Walk {
        options,
        interrupt,
        work: initial,
        ids,
        visited,
        pending,
    };
    if !document
        .represented_members()
        .iter()
        .any(|m| m.item_entity() == root)
    {
        return Err(Error::Invalid);
    }
    walk.reference(root)?;
    while let Some(index) = walk.pending.pop() {
        walk.charge(32)?;
        for entity in &document.records()[&walk.ids[index]] {
            walk.values(&entity.args)?;
        }
    }
    let mut index = 0;
    walk.ids.retain(|_| {
        let keep = walk.visited[index] != 0;
        index += 1;
        keep
    });
    for index in 0..walk.ids.len() {
        walk.charge(32)?;
        let start = components.len();
        for entity in &document.records()[&walk.ids[index]] {
            walk.charge(1)?;
            components.push(entity.source.expect("original parsed component"));
        }
        ranges.push(start..components.len());
    }
    let retained = std::mem::size_of::<SelectedRecords<'_>>()
        + walk.ids.capacity() * std::mem::size_of::<u32>()
        + ranges.capacity() * std::mem::size_of::<std::ops::Range<usize>>()
        + components.capacity() * std::mem::size_of::<&str>();
    Ok(SelectedRecords {
        ids: walk.ids,
        ranges,
        components,
        work: walk.work,
        retained,
        peak: bytes,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    fn input() -> String {
        let mut writer = crate::step::Writer::new();
        for x in [0., 5.] {
            writer
                .append_sphere(
                    crate::sphere::Sphere::new(1.).unwrap(),
                    [x, 0., 0.],
                    [1., 0., 0.],
                    [0., 1., 0.],
                )
                .unwrap();
        }
        writer.finish().unwrap()
    }
    #[test]
    fn selected_original_components_outlive_arena_and_enforce_real_neutral_frontiers() {
        let input = input();
        let document = CompoundDocument::parse(&input, Default::default(), None).unwrap();
        let root = document.represented_members()[0].item_entity();
        let other = document.represented_members()[1].item_entity();
        let calls = std::cell::Cell::new(0);
        let poll = || {
            calls.set(calls.get() + 1);
            false
        };
        let records = selected_records(&document, root, Default::default(), Some(&poll)).unwrap();
        assert!(records.record_index(root).is_some() && records.record_index(other).is_none());
        let exact = Options {
            maximum_work: records.work_used(),
            maximum_bytes: records.query_peak_bytes(),
            ..Default::default()
        };
        selected_records(&document, root, exact, None).unwrap();
        for o in [
            Options {
                maximum_work: 1,
                ..Default::default()
            },
            Options {
                maximum_work: exact.maximum_work - 1,
                ..Default::default()
            },
        ] {
            assert!(matches!(
                selected_records(&document, root, o, None),
                Err(Error::WorkLimit)
            ));
        }
        for bytes in [1, exact.maximum_bytes - 1] {
            assert!(matches!(
                selected_records(
                    &document,
                    root,
                    Options {
                        maximum_bytes: bytes,
                        ..Default::default()
                    },
                    None
                ),
                Err(Error::StorageLimit)
            ));
        }
        for o in [
            Options {
                maximum_work: 0,
                ..Default::default()
            },
            Options {
                maximum_work: 512_000_001,
                ..Default::default()
            },
            Options {
                maximum_bytes: 0,
                ..Default::default()
            },
            Options {
                maximum_bytes: 64 * 1024 * 1024 + 1,
                ..Default::default()
            },
            Options {
                maximum_entities: 0,
                ..Default::default()
            },
            Options {
                maximum_entities: 1_048_577,
                ..Default::default()
            },
        ] {
            assert!(matches!(
                selected_records(&document, root, o, None),
                Err(Error::Invalid)
            ));
        }
        assert!(matches!(
            selected_records(
                &document,
                root,
                Options {
                    maximum_entities: 1,
                    ..Default::default()
                },
                None
            ),
            Err(Error::EntityLimit)
        ));
        assert!(matches!(
            selected_records(&document, u32::MAX, Default::default(), None),
            Err(Error::Invalid)
        ));
        for boundary in [1, 2, calls.get()] {
            let n = std::cell::Cell::new(0);
            let stop = || {
                n.set(n.get() + 1);
                n.get() >= boundary
            };
            assert!(matches!(
                selected_records(&document, root, Default::default(), Some(&stop)),
                Err(Error::Cancelled)
            ));
        }
        let surface = document
            .records()
            .values()
            .flatten()
            .find(|e| e.name == "SPHERICAL_SURFACE")
            .unwrap();
        let Value::Reference(axis) = surface.args[1] else {
            panic!("real sphere axis");
        };
        let body = surface.source.unwrap();
        let changed = body.replacen(&format!("#{axis}"), "#999999", 1);
        let malformed = input.replacen(body, &changed, 1);
        let neutral = CompoundDocument::parse(&malformed, Default::default(), None).unwrap();
        assert!(matches!(
            selected_records(&neutral, root, Default::default(), None),
            Err(Error::Invalid)
        ));
        drop(document);
        assert!(
            records
                .component(records.record_index(root).unwrap(), 0)
                .unwrap()
                .starts_with("MANIFOLD_SOLID_BREP")
        );
    }
}
