//! Stable endpoint ordering for already validated indexed edge uses.
use super::*;
fn storage(edges: usize, scratch: usize, counts: usize) -> Result<usize, MeshError> {
    edges
        .checked_add(scratch)
        .and_then(|n| n.checked_mul(std::mem::size_of::<[u32; 3]>()))
        .and_then(|n| {
            counts
                .checked_mul(std::mem::size_of::<usize>())
                .and_then(|c| n.checked_add(c))
        })
        .ok_or(MeshError::BudgetExceeded)
}
/// Endpoints are below `vertices`; original coedge ordinals are increasing.
/// Both invariants are established in certify_impl before this private call.
pub(super) fn sort<E, const COUNTED: bool>(
    edges: &mut Vec<[u32; 3]>,
    vertices: usize,
    budget: &mut Budget<E, COUNTED>,
) -> Result<(), CountedError<E>> {
    budget.storage(storage(edges.capacity(), edges.len(), vertices)?)?;
    // One logical edge-record initialization and one scalar initialization
    // each, before the corresponding bounded allocation/copy.
    budget.work(edges.len() + vertices)?;
    let mut scratch = vec![[0; 3]; edges.len()];
    let mut counts = vec![0usize; vertices];
    budget.storage(storage(
        edges.capacity(),
        scratch.capacity(),
        counts.capacity(),
    )?)?;
    pass(edges, &mut scratch, &mut counts, 1, budget)?;
    budget.work(counts.len())?;
    counts.fill(0);
    pass(&scratch, edges, &mut counts, 0, budget)
}
fn pass<E, const COUNTED: bool>(
    source: &[[u32; 3]],
    target: &mut [[u32; 3]],
    counts: &mut [usize],
    axis: usize,
    budget: &mut Budget<E, COUNTED>,
) -> Result<(), CountedError<E>> {
    for edge in source {
        // One bounded histogram update, using a previously validated key.
        budget.work(1)?;
        counts[edge[axis] as usize] += 1;
    }
    let mut offset = 0;
    for count in counts.iter_mut() {
        // One bounded prefix update. All counts sum to source.len(), so the
        // addition and every later destination remain within slice bounds.
        budget.work(1)?;
        let occurrences = *count;
        *count = offset;
        offset += occurrences;
    }
    for edge in source {
        // Same literal record-copy unit as the old merge sort, plus one
        // bounded offset update. Source order makes each pass stable.
        budget.work(2)?;
        let destination = &mut counts[edge[axis] as usize];
        target[*destination] = *edge;
        *destination += 1;
    }
    Ok(())
}
#[cfg(test)]
mod tests;
