//! Vertex links as cycles of the already validated corner permutation.
use super::*;
/// certify_impl supplies next(opposite(c)) for EVERY corner after validating
/// unique opposite edge uses and winding. This is a permutation preserving the
/// original vertex index. Its cycles equal the old fan-union components.
pub(super) fn validate<E, const COUNTED: bool>(
    successors: &mut Vec<u32>,
    indices: &[u32],
    vertices: usize,
    shell_capacity: usize,
    budget: &mut Budget<E, COUNTED>,
) -> Result<u32, CountedError<E>> {
    let bytes = |links: usize| {
        successors
            .capacity()
            .checked_add(shell_capacity)
            .and_then(|n| n.checked_add(links))
            .and_then(|n| n.checked_mul(std::mem::size_of::<u32>()))
            .ok_or(MeshError::BudgetExceeded)
    };
    budget.storage(bytes(vertices)?)?;
    budget.work(vertices)?;
    let mut links = vec![u32::MAX; vertices];
    budget.storage(bytes(links.capacity())?)?;
    let mut vertex_links = 0;
    for start in 0..successors.len() {
        budget.work(1)?;
        if successors[start] == u32::MAX {
            continue;
        }
        // A second cycle at the SAME vertex is the non-manifold bow tie that
        // the original union/root comparison rejected. Unused vertices remain
        // absent, just as before.
        budget.work(2)?;
        let link = &mut links[indices[start] as usize];
        if *link != u32::MAX {
            return Err(INVALID.into());
        }
        *link = start as u32;
        vertex_links += 1;
        let mut corner = start;
        loop {
            // One bounded link-consumption action: read successor, mark it
            // visited, advance. This is the same scalar-update abstraction as
            // one original path-compression iteration.
            budget.work(1)?;
            let successor = successors[corner];
            successors[corner] = u32::MAX;
            corner = successor as usize;
            if corner == start {
                break;
            }
        }
    }
    Ok(vertex_links)
}
