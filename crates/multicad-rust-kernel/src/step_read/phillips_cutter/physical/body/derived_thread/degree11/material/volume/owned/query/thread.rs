//! Whole(theta,z) ray exclusion on actual coupled degree11 source cells.
//! The occupied constructor already proves a one-sheet positive radial graph
//! covering every angle and strict cap height. For a nonzero query point its
//! unique same-angle/same-z boundary must occur in one of these actual cells.
//! Exclusion tests bound ALL source/query pairs, never an equal-v section.
use super::*;
use crate::tapered_thread::TensorCell11;
pub(super) const NODE_WORK: usize = 8192;
const SOURCE_WORK: usize = 8192;
// Covers144-control copies, both transpose overlaps and complete12-row,
// degree11 de Casteljau loops including interval operations and Vec handling.
const SPLIT_WORK: usize = 262144;
// Same conservative16x16 interval-triple cell envelope as occupied.rs.
const CELL_BYTES: usize = 16384;
type Pending = [Option<(TensorCell11, usize)>; MAX_THREAD_DEPTH + 2];
pub(super) fn scratch_bytes(depth: usize) -> usize {
    // depth+1 DFS cells, eight simultaneous split/transpose/row workspaces.
    // Fixed pending headers/receipt/work state are in METADATA_SCRATCH.
    (depth + 9) * CELL_BYTES
}
fn outside_projection(p: [I; 3], cell: [I; 3]) -> bool {
    if cell[2].hi < p[2].lo || cell[2].lo > p[2].hi {
        return true;
    }
    let cross = p[0] * cell[1] - p[1] * cell[0];
    if cross.lo > 0. || cross.hi < 0. {
        return true;
    }
    let dot = p[0] * cell[0] + p[1] * cell[1];
    // Positive rays require positive dot. Axis points are handled by the
    // guaranteed occupied core; they cannot justify an outside predicate.
    dot.hi <= 0.
}
pub(super) fn classify(
    source: &OwnedDegree11Material,
    p: [I; 3],
    depth: usize,
    work: &mut QueryWork<'_>,
) -> R<C> {
    let [low, high] = source.cap_stations_mm();
    if p[2].hi < low || p[2].lo > high {
        return Ok(C::Outside);
    }
    if p[2].lo <= low || p[2].hi >= high {
        return Ok(C::UnresolvedBoundary);
    }
    let radius = (p[0].square() + p[1].square()).sqrt();
    if radius.hi < source.occupied_core_radius_mm() {
        return Ok(C::Inside);
    }
    if prove(source, p, radius, true, depth, work)? {
        return Ok(C::Inside);
    }
    if radius.lo > 0. && prove(source, p, radius, false, depth, work)? {
        return Ok(C::Outside);
    }
    Ok(C::UnresolvedBoundary)
}
fn prove(
    source: &OwnedDegree11Material,
    p: [I; 3],
    radius: I,
    inside: bool,
    maximum_depth: usize,
    work: &mut QueryWork<'_>,
) -> R<bool> {
    for face in source.thread_faces() {
        let [u, v] = face.source().spans();
        for i in 0..u {
            for j in 0..v {
                work.thread_charge(SOURCE_WORK)?;
                let root = face
                    .source()
                    .cell([i, j])
                    .map_err(Degree11MaterialQueryError::Thread)?;
                let mut pending: Pending = std::array::from_fn(|_| None);
                pending[0] = Some((root, 0));
                let mut count = 1;
                while count > 0 {
                    count -= 1;
                    let (cell, depth) = pending[count].take().expect("tracked bounded DFS entry");
                    work.cell()?;
                    let b = cell
                        .bounds()
                        .map_err(Degree11MaterialQueryError::Thread)?
                        .map(|b| I::new(b[0], b[1]));
                    if outside_projection(p, b) {
                        continue;
                    }
                    let boundary_radius = (b[0].square() + b[1].square()).sqrt();
                    if if inside {
                        boundary_radius.lo > radius.hi
                    } else {
                        boundary_radius.hi < radius.lo
                    } {
                        continue;
                    }
                    // An apparently opposite radial sign can be a false candidate
                    // from coupled(theta,z) hulls. Refine it too; no equal-v or
                    // sampled inverse is allowed to turn it into a witness.
                    if depth == maximum_depth {
                        return Ok(false);
                    }
                    work.thread_charge(SPLIT_WORK)?;
                    let [a, b] = cell
                        .subdivide(depth % 2)
                        .map_err(Degree11MaterialQueryError::Thread)?;
                    pending[count] = Some((b, depth + 1));
                    pending[count + 1] = Some((a, depth + 1));
                    count += 2;
                }
            }
        }
    }
    // Every source cell was excluded from the complete ray/z query or proved
    // on the required radial side. Earned coverage supplies each nonzero ray;
    // axis points are inside the core and outside proof requires radius.lo>0.
    Ok(true)
}

#[cfg(test)]
mod tests;
