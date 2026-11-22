//! One cumulative adaptive source-integration budget; no partial publication.
use super::*;
use measure::{Density, QUANTITIES, Selection, Values};
use source::Source;
use std::collections::BTreeMap;
pub(super) struct Accuracy {
    pub scales: Values,
    pub epsilon: f64,
    pub selection: Selection,
}
#[derive(Clone, Copy)]
pub(super) struct Limits {
    pub cells: usize,
    pub depth: u8,
}
struct Cell {
    source: Source,
    value: Values,
    depth: u8,
}
pub(super) fn face(
    sources: Vec<Source>,
    density: Density,
    chart: &axial::Chart,
    center: f64,
    accuracy: Accuracy,
    work: &mut Work<'_>,
) -> Result<Values, Error> {
    with_limits(
        sources,
        density,
        chart,
        center,
        accuracy,
        Limits {
            cells: 32768,
            depth: 48,
        },
        work,
    )
}
pub(super) fn with_limits(
    sources: Vec<Source>,
    density: Density,
    chart: &axial::Chart,
    center: f64,
    accuracy: Accuracy,
    limits: Limits,
    work: &mut Work<'_>,
) -> Result<Values, Error> {
    with_acceptance(
        sources, density, chart, center, accuracy, limits, None, work,
    )
}
/// Surface-only integration accepts the real corrected output widths, instead
/// of imposing an independent, uniformly overspecified moment tolerance.
pub(super) fn surface_face(
    sources: Vec<Source>,
    density: Density,
    chart: &axial::Chart,
    center: f64,
    scales: Values,
    acceptance: &dyn Fn(&Values) -> bool,
    work: &mut Work<'_>,
) -> Result<Values, Error> {
    with_acceptance(
        sources,
        density,
        chart,
        center,
        Accuracy {
            scales,
            epsilon: 0.,
            selection: Selection::Surface,
        },
        Limits {
            cells: 32768,
            depth: 48,
        },
        Some(acceptance),
        work,
    )
}
#[allow(clippy::too_many_arguments)]
fn with_acceptance(
    sources: Vec<Source>,
    density: Density,
    chart: &axial::Chart,
    center: f64,
    accuracy: Accuracy,
    limits: Limits,
    acceptance: Option<&dyn Fn(&Values) -> bool>,
    work: &mut Work<'_>,
) -> Result<Values, Error> {
    let Accuracy {
        scales,
        epsilon,
        selection,
    } = accuracy;
    let mut heap: BTreeMap<u64, Vec<Cell>> = BTreeMap::new();
    let mut count = 0;
    let mut total = [I::point(0.); QUANTITIES];
    for source in sources {
        for cell in ready(source, 0, density, chart, center, limits, selection, work)? {
            for i in selection.indices() {
                total[i] = total[i] + cell.value[i];
            }
            insert(&mut heap, cell, &scales);
            count += 1;
            if count > limits.cells {
                return Err(Error::Limit);
            }
        }
    }
    let mut iterations = 0;
    loop {
        iterations += 1;
        // Endpoint subtraction remains outward, but very coarse early
        // cells can leave a cancellation floor. Rebuild from live cells
        // before mistaking that arithmetic history for source uncertainty.
        if iterations % 32 == 0 {
            total = [I::point(0.); QUANTITIES];
            for cell in heap.values().flatten() {
                work.charge(QUANTITIES)?;
                for i in selection.indices() {
                    total[i] = total[i] + cell.value[i];
                }
            }
        }
        work.step()?;
        let accepted = if let Some(acceptance) = acceptance {
            // The fixed 3D acceptance check uses at most ten interval
            // operations for collar corrections, 33 for world transport,
            // and seven for the four directed widths/centroid targets.
            // Reserve 64, including scalar negations and fixed checks,
            // inside the same caller ledger; no epsilon decides acceptance.
            work.charge(64)?;
            acceptance(&total)
        } else {
            total
                .iter()
                .zip(scales)
                .all(|(v, s)| v.width() <= (I::point(epsilon) * s).lo)
        };
        if accepted {
            return Ok(total);
        }
        let (key, mut bucket) = heap.pop_last().ok_or(Error::Unresolved)?;
        // Only nonempty buckets are inserted. Keep the remaining bucket
        // without introducing a geometry-dependent ordering tie-breaker.
        let old = bucket
            .pop()
            .expect("only nonempty score buckets are inserted");
        if !bucket.is_empty() {
            heap.insert(key, bucket);
        }
        count -= 1;
        if old.depth >= limits.depth {
            return Err(Error::Unresolved);
        }
        let along_v = old.source.dimensions() == 2 && old.depth % 2 == 1;
        let children = old.source.subdivide(along_v)?;
        for i in selection.indices() {
            total[i] = I::new(
                (I::point(total[i].lo) - I::point(old.value[i].lo)).lo,
                (I::point(total[i].hi) - I::point(old.value[i].hi)).hi,
            );
        }
        for child in children {
            for cell in ready(
                child,
                old.depth + 1,
                density,
                chart,
                center,
                limits,
                selection,
                work,
            )? {
                for i in selection.indices() {
                    total[i] = total[i] + cell.value[i];
                }
                insert(&mut heap, cell, &scales);
                count += 1;
                if count > limits.cells {
                    return Err(Error::Limit);
                }
            }
        }
    }
}
fn insert(heap: &mut BTreeMap<u64, Vec<Cell>>, cell: Cell, scales: &Values) {
    let score = cell
        .value
        .iter()
        .zip(scales)
        .map(|(v, s)| v.width() / s.lo)
        .fold(0., f64::max);
    heap.entry(score.to_bits()).or_default().push(cell);
}
#[allow(clippy::too_many_arguments)]
fn ready(
    source: Source,
    depth: u8,
    density: Density,
    chart: &axial::Chart,
    center: f64,
    limits: Limits,
    selection: Selection,
    work: &mut Work<'_>,
) -> Result<Vec<Cell>, Error> {
    let mut pending = vec![(source, depth)];
    let mut ready = Vec::new();
    while let Some((source, depth)) = pending.pop() {
        let enclosure = match selection {
            Selection::All => quadrature::enclosure(&source, density, chart, center, work),
            Selection::Volume => {
                quadrature::volume_enclosure(&source, density, chart, center, work)
            }
            Selection::Surface => {
                quadrature::surface_enclosure(&source, density, chart, center, work)
            }
        };
        match enclosure {
            Ok(value) => ready.push(Cell {
                source,
                value,
                depth,
            }),
            Err(Error::Unresolved) => {
                if depth >= limits.depth {
                    return Err(Error::Unresolved);
                }
                let along_v = source.dimensions() == 2 && depth % 2 == 1;
                let [a, b] = source.subdivide(along_v)?;
                pending.push((b, depth + 1));
                pending.push((a, depth + 1));
            }
            Err(error) => return Err(error),
        }
        if pending.len() + ready.len() > limits.cells {
            return Err(Error::Limit);
        }
    }
    Ok(ready)
}
