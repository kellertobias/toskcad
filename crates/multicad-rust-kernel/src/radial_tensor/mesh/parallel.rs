//! Bounded parallel checks of immutable retained source cells.
use super::*;
mod cells;
mod credits;
const MIN_FACETS: usize = 8192;
const MIN_WORK: usize = 2_097_152;
pub(in crate::radial_tensor) fn refill(balance: &mut usize, needed: usize) -> Result<(), Error> {
    credits::refill(balance, needed)
}
pub(in crate::radial_tensor) fn local_budget<T>(f: impl FnOnce() -> T) -> T {
    credits::local_budget(f)
}
#[allow(clippy::too_many_arguments)]
pub(super) fn check(
    p: &mut Prepared<'_>,
    mesh: &DisplayMesh,
    grid: &Grid,
    ring: &Ring,
    pass: usize,
    changes: &mut Refinement,
    accepted: &mut Vec<[usize; 3]>,
    work: &mut Work<'_>,
) -> Result<bool, Error> {
    // Coverage instrumentation (cargo-llvm-cov sets `cfg(coverage)`) makes
    // every worker update the same profile counters, so this cell check runs
    // several times slower in parallel than serially there; the R 1/8 UI
    // preview went from 9.5 s serial to over 45 s with eight workers. Serial
    // and parallel checks emit bit-identical meshes (paired BSPT test).
    let workers = if cfg!(coverage) {
        1
    } else {
        std::thread::available_parallelism()
            .map_or(1, usize::from)
            .min(8)
    };
    #[cfg(test)]
    let workers = parallel_probe::WORKERS.with(Cell::get).unwrap_or(workers);
    let threshold = MIN_FACETS;
    #[cfg(test)]
    let threshold = parallel_probe::MIN_FACETS
        .with(Cell::get)
        .unwrap_or(threshold);
    if workers < 2
        || !p.tensor_cycle
        || mesh.indices.len() / 3 < threshold
        || work.remaining < MIN_WORK
    {
        return Ok(false);
    }
    cells::check(p, mesh, grid, ring, pass, changes, accepted, workers, work)
}
#[cfg(test)]
use std::cell::Cell;
