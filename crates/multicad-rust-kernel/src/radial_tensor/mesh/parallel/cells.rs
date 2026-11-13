use super::*;
use std::cell::Cell;
use std::{
    rc::Rc,
    sync::{
        Arc, Mutex,
        atomic::{AtomicUsize, Ordering},
        mpsc,
    },
    thread,
    time::Duration,
};
#[cfg(test)]
struct ActiveGuard {
    old: bool,
}
#[cfg(test)]
impl ActiveGuard {
    fn new() -> Self {
        Self {
            old: super::super::parallel_probe::ACTIVE.with(|a| a.replace(true)),
        }
    }
}
#[cfg(test)]
impl Drop for ActiveGuard {
    fn drop(&mut self) {
        super::super::parallel_probe::ACTIVE.with(|a| a.set(self.old));
    }
}
#[derive(Clone, Copy)]
struct Job {
    side: usize,
    u: [usize; 2],
    v: [usize; 2],
    first: usize,
    source_keys: usize,
}
struct Record {
    index: usize,
    at: [usize; 3],
    key: [u64; 18],
    cached: bool,
    verdict: Option<CellCheck>,
    axis: Option<usize>,
}
struct Done {
    job: usize,
    #[cfg(test)]
    worker: usize,
    #[cfg(test)]
    elapsed: Duration,
    result: Result<Vec<Record>, (usize, Error)>,
}
fn cuts(nodes: &[f64], knots: &[f64]) -> Vec<usize> {
    (0..nodes.len())
        .filter(|&i| i == 0 || i + 1 == nodes.len() || knots.contains(&nodes[i]))
        .collect()
}
fn jobs(
    p: &Prepared<'_>,
    grid: &Grid,
    ring: &Ring,
    work: &mut Work<'_>,
) -> Result<Option<Vec<Job>>, Error> {
    let mut jobs = Vec::new();
    let mut owners = BTreeSet::new();
    let mut cells = 0usize;
    for (side, s) in p.sides.iter().enumerate() {
        let Side::Tensor(t) = s else {
            return Err(Error::Unsupported);
        };
        let us = cuts(&grid.u[side], &t.knots[0]);
        let vs = cuts(&grid.v, &t.knots[1]);
        for v in vs.windows(2) {
            for u in us.windows(2) {
                if jobs.len() == p.limits.parallel_jobs {
                    return Ok(None);
                };
                work.step()?;
                let range = [
                    [grid.u[side][u[0]], grid.u[side][u[1]]],
                    [grid.v[v[0]], grid.v[v[1]]],
                ];
                // Position and normal strip keys both contain the actual
                // patch domain and the complete raw U-cell interval.
                let mut source_keys = 0usize;
                for patch in t.patches(range, work)? {
                    let lowers = &grid.u[side][u[0]..u[1]];
                    let uppers = &grid.u[side][u[0] + 1..=u[1]];
                    let first = uppers.partition_point(|&x| x <= patch.domain[0][0]);
                    let last = lowers.partition_point(|&x| x < patch.domain[0][1]);
                    source_keys = source_keys
                        .checked_add(last.checked_sub(first).ok_or(Error::Invalid)?)
                        .ok_or(Error::Limit)?;
                    let identity = (side, patch.domain.map(|r| r.map(f64::to_bits)));
                    if !register_source(&mut owners, identity, p.limits.parallel_sources)? {
                        return Ok(None);
                    }
                }
                cells = cells
                    .checked_add((u[1] - u[0]).checked_mul(v[1] - v[0]).ok_or(Error::Limit)?)
                    .ok_or(Error::Limit)?;
                jobs.push(Job {
                    side,
                    u: [u[0], u[1]],
                    v: [v[0], v[1]],
                    first: v[0] * ring.entries.len() + ring.offsets[side] + u[0],
                    source_keys,
                });
            }
        }
    }
    let expected = (grid.v.len() - 1) * ring.entries.len();
    if cells != expected {
        return Err(Error::Invalid);
    }
    jobs.sort_by_key(|j| j.first);
    Ok(Some(jobs))
}
#[allow(clippy::too_many_arguments)]
fn run_job(
    p: &Prepared<'_>,
    mesh: &DisplayMesh,
    grid: &Grid,
    ring: &Ring,
    pass: usize,
    job: Job,
    current: &Cell<usize>,
    mut records: Vec<Record>,
    cache_limit: usize,
    w: &mut Work<'_>,
) -> Result<Vec<Record>, (usize, Error)> {
    let Side::Tensor(t) = &p.sides[job.side] else {
        return Err((job.first, Error::Unsupported));
    };
    let mut cache = strip::Cache::new(t, cache_limit);
    let result = (|| -> Result<(), Error> {
        for row in job.v[0]..job.v[1] {
            current.set(row * ring.entries.len() + ring.offsets[job.side] + job.u[0]);
            cache.advance(grid.v[row], w)?;
            for u in job.u[0]..job.u[1] {
                let side = job.side;
                let index = row * ring.entries.len() + ring.offsets[side] + u;
                current.set(index);
                w.step()?;
                let range = [
                    [grid.u[side][u], grid.u[side][u + 1]],
                    [grid.v[row], grid.v[row + 1]],
                ];
                let ids = [
                    ring.index(side, u, row),
                    ring.index(side, u + 1, row),
                    ring.index(side, u + 1, row + 1),
                    ring.index(side, u, row + 1),
                ];
                let vertices = ids.map(|i| point(mesh, i));
                let actual = std::array::from_fn::<_, 2, _>(|tri| {
                    std::array::from_fn(|k| {
                        let id = mesh.indices[6 * index + 3 * tri + k] as usize;
                        ids.iter()
                            .position(|&i| i == id)
                            .expect("actual quad corners")
                    })
                });
                let topology = actual
                    .into_iter()
                    .flatten()
                    .enumerate()
                    .fold(0u64, |key, (i, v)| key | ((v as u64) << (2 * i)));
                let mut key = [0; 18];
                key[..6].copy_from_slice(&[
                    side as u64,
                    range[0][0].to_bits(),
                    range[0][1].to_bits(),
                    range[1][0].to_bits(),
                    range[1][1].to_bits(),
                    topology,
                ]);
                for (slot, value) in key[6..].iter_mut().zip(vertices.into_iter().flatten()) {
                    *slot = value.to_bits();
                }
                w.charge(128)?;
                let cached = p.certificates.contains(&key);
                let verdict = if cached {
                    Some(CellCheck::Match)
                } else {
                    match p.cell_check(
                        CandidateCell {
                            side,
                            range,
                            vertices,
                            actual,
                        },
                        None,
                        Some(&mut cache),
                        w,
                    ) {
                        Ok(v) => Some(v),
                        Err(Error::Unresolved) => None,
                        Err(e) => return Err(e),
                    }
                };
                let axis = if cached || verdict == Some(CellCheck::Match) {
                    None
                } else {
                    Some(p.hint(side, range, pass, w)?)
                };
                records.push(Record {
                    index,
                    at: [side, u, row],
                    key,
                    cached,
                    verdict,
                    axis,
                });
            }
        }
        Ok(())
    })();
    result.map(|_| records).map_err(|e| (current.get(), e))
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
    requested: usize,
    work: &mut Work<'_>,
) -> Result<bool, Error> {
    #[cfg(test)]
    let setup_at = std::time::Instant::now();
    let _cache_total = p
        .limits
        .strips
        .checked_mul(p.sides.len())
        .ok_or(Error::Limit)?;

    let Some(jobs) = jobs(p, grid, ring, work)? else {
        return Ok(false);
    };
    let n = jobs.len().min(requested).clamp(1, 8);
    let mut per_owner = vec![0usize; p.sides.len()];
    for job in &jobs {
        per_owner[job.side] += 1;
    }
    // Each live job owns one cache. Even if every worker selects the same
    // side, its two cache tables retain at most the old per-owner entry cap.
    let quotas: Vec<_> = per_owner
        .iter()
        .map(|&count| p.limits.strips / n.min(count).max(1))
        .collect();
    // Full lifetime key coverage is conservative; avoid extra source work
    // when splitting a small cache would force recomputation.
    if p.limits.strips != 0 && jobs.iter().any(|j| j.source_keys > quotas[j.side]) {
        #[cfg(test)]
        eprintln!(
            "PARALLEL_CACHE_FALLBACK strips={} max_keys={} quotas={quotas:?}",
            p.limits.strips,
            jobs.iter().map(|j| j.source_keys).max().unwrap_or(0)
        );
        return Ok(false);
    }
    let initial_bytes = jobs
        .capacity()
        .checked_mul(std::mem::size_of::<Job>())
        .ok_or(Error::Limit)?;
    let Some((buffers, bytes)) = record_buffers(&jobs, initial_bytes, p.limits.parallel_bytes)
    else {
        return Ok(false);
    };
    work.charge(jobs.len() + bytes / 64)?;
    let pool = credits::Pool::new(work.remaining);
    let next = AtomicUsize::new(0);
    let buffers = Mutex::new(buffers);
    let (tx, rx) = mpsc::channel();
    let mut completed = Vec::new();
    let original_interrupt = work.interrupt;
    let stack_bytes = p.limits.parallel_stack;
    let prepared: &Prepared<'_> = p;
    let mut external_cancel = false;
    let mut spawn_failure = false;
    #[cfg(test)]
    let setup_wall = setup_at.elapsed();
    #[cfg(test)]
    let dispatch_at = std::time::Instant::now();
    #[cfg(test)]
    let active = ActiveGuard::new();
    thread::scope(|scope| {
        let mut handles = Vec::new();
        for _worker in 0..n {
            let quotas = &quotas;
            let tx = tx.clone();
            let worker_pool = Arc::clone(&pool);
            let jobs = &jobs;
            let next = &next;
            let buffers = &buffers;
            let spawn = thread::Builder::new()
                .name("retained-certificate".into())
                .stack_size(stack_bytes)
                .spawn_scoped(scope, move || {
                    worker_loop(
                        _worker,
                        prepared,
                        mesh,
                        grid,
                        ring,
                        pass,
                        jobs,
                        next,
                        buffers,
                        tx,
                        worker_pool,
                        quotas,
                    )
                });
            match spawn {
                Ok(h) => handles.push(h),
                Err(_) => {
                    spawn_failure = true;
                    pool.cancel();
                    break;
                }
            }
        }
        drop(tx);
        loop {
            if original_interrupt.is_some_and(|f| f()) {
                external_cancel = true;
                pool.cancel();
            }
            match rx.recv_timeout(Duration::from_millis(1)) {
                Ok(done) => completed.push(done),
                Err(mpsc::RecvTimeoutError::Timeout) => continue,
                Err(mpsc::RecvTimeoutError::Disconnected) => break,
            }
        }
        spawn_failure |= join_workers(handles, &pool);
    });
    #[cfg(test)]
    let dispatch_wall = dispatch_at.elapsed();
    #[cfg(test)]
    drop(active);
    let (remaining, _spent, _refills) = pool.finish()?;
    work.remaining = remaining;
    if original_interrupt.is_some_and(|f| f()) {
        external_cancel = true;
    }
    if external_cancel {
        return Err(Error::Cancelled);
    }
    if spawn_failure {
        return Err(Error::Limit);
    }
    verify_completed(&completed, jobs.len())?;
    #[cfg(test)]
    let mut worker_wall = vec![Duration::ZERO; n];
    #[cfg(test)]
    for done in &completed {
        worker_wall[done.worker] += done.elapsed;
    }
    #[cfg(test)]
    let merge_at = std::time::Instant::now();
    apply_completed(p, grid, ring, changes, accepted, completed, work)?;
    #[cfg(test)]
    let merge_wall = merge_at.elapsed();
    #[cfg(test)]
    eprintln!(
        "PARALLEL_TIMING workers={n} setup={setup_wall:?} dispatch={dispatch_wall:?} coordinator_merge_apply={merge_wall:?} jobs_sum={:?} jobs_per_worker={worker_wall:?} credits={:?}",
        worker_wall.iter().sum::<Duration>(),
        pool.timing()
    );
    #[cfg(test)]
    let expected = (grid.v.len() - 1) * ring.entries.len();
    #[cfg(test)]
    eprintln!(
        "PARALLEL_CELLS jobs={} workers={n} cells={expected} scalar_work={_spent} credit_refills={_refills} record_bytes={bytes}",
        jobs.len()
    );
    Ok(true)
}

fn register_source(
    owners: &mut BTreeSet<(usize, [[u64; 2]; 2])>,
    key: (usize, [[u64; 2]; 2]),
    limit: usize,
) -> Result<bool, Error> {
    if owners.len() >= limit {
        return Ok(false);
    }
    if !owners.insert(key) {
        return Err(Error::Invalid);
    }
    Ok(true)
}
fn account_records(bytes: &mut usize, count: usize, limit: usize) -> Option<()> {
    let total = bytes.checked_add(count.checked_mul(std::mem::size_of::<Record>())?)?;
    if total > limit {
        return None;
    }
    *bytes = total;
    Some(())
}
fn record_buffer(count: usize, bytes: &mut usize, limit: usize) -> Option<Vec<Record>> {
    let mut preview = *bytes;
    account_records(&mut preview, count, limit)?;
    let mut records = Vec::new();
    records.try_reserve_exact(count).ok()?;
    account_records(bytes, records.capacity(), limit)?;
    Some(records)
}
#[cfg(test)]
mod tests;

#[allow(clippy::too_many_arguments)]
fn worker_loop(
    _worker: usize,
    prepared: &Prepared<'_>,
    mesh: &DisplayMesh,
    grid: &Grid,
    ring: &Ring,
    pass: usize,
    jobs: &[Job],
    next: &AtomicUsize,
    buffers: &Mutex<Vec<Option<Vec<Record>>>>,
    tx: mpsc::Sender<Done>,
    pool: Arc<credits::Pool>,
    quotas: &[usize],
) {
    loop {
        let number = next.fetch_add(1, Ordering::AcqRel);
        let Some(&job) = jobs.get(number) else { break };
        if job.first > pool.stop.load(Ordering::Acquire) || pool.cancel.load(Ordering::Acquire) {
            break;
        }
        let records = match buffers.lock() {
            Ok(mut b) => b[number].take().expect("one owner per job"),
            Err(_) => {
                let _ = tx.send(Done {
                    job: number,
                    #[cfg(test)]
                    worker: _worker,
                    #[cfg(test)]
                    elapsed: Duration::ZERO,
                    result: Err((job.first, Error::Unresolved)),
                });
                break;
            }
        };
        let current = Rc::new(Cell::new(job.first));
        #[cfg(test)]
        let job_at = std::time::Instant::now();
        let run = credits::run(Arc::clone(&pool), Rc::clone(&current), |w| {
            run_job(
                prepared,
                mesh,
                grid,
                ring,
                pass,
                job,
                &current,
                records,
                quotas[job.side],
                w,
            )
            .map_err(|(_, e)| e)
        });
        let result = run.map_err(|e| (current.get(), e));
        if let Err((at, e)) = result
            && e != Error::Cancelled
        {
            pool.halt_after(at);
        }
        if tx
            .send(Done {
                job: number,
                #[cfg(test)]
                worker: _worker,
                #[cfg(test)]
                elapsed: job_at.elapsed(),
                result,
            })
            .is_err()
        {
            break;
        }
    }
}
fn join_workers<'s>(handles: Vec<thread::ScopedJoinHandle<'s, ()>>, pool: &credits::Pool) -> bool {
    let mut spawn_failure = false;
    for handle in handles {
        if handle.join().is_err() {
            spawn_failure = true;
            pool.cancel();
        }
    }
    spawn_failure
}
fn verify_completed(completed: &[Done], expected: usize) -> Result<(), Error> {
    let fatal = completed
        .iter()
        .filter_map(|done| done.result.as_ref().err().copied())
        .filter(|(_, e)| *e != Error::Cancelled)
        .min_by_key(|(at, _)| *at);
    if let Some((_, e)) = fatal {
        return Err(e);
    }
    if completed.len() != expected {
        return Err(Error::Unresolved);
    }
    Ok(())
}
fn apply_completed(
    p: &mut Prepared<'_>,
    grid: &Grid,
    ring: &Ring,
    changes: &mut Refinement,
    accepted: &mut Vec<[usize; 3]>,
    mut completed: Vec<Done>,
    work: &mut Work<'_>,
) -> Result<(), Error> {
    completed.sort_by_key(|d| d.job);
    let mut streams = completed
        .into_iter()
        .map(|done| done.result.map(Vec::into_iter).map_err(|(_, e)| e))
        .collect::<Result<Vec<_>, _>>()?;
    let mut heads: Vec<_> = streams.iter_mut().map(Iterator::next).collect();
    let mut heap = std::collections::BinaryHeap::new();
    for (job, head) in heads.iter().enumerate() {
        if let Some(record) = head {
            heap.push(std::cmp::Reverse((record.index, job)));
        }
    }
    let expected = (grid.v.len() - 1) * ring.entries.len();
    let compensate = p.limits.compensate_v && p.tensor_cycle;
    let mut processed = 0;
    while let Some(std::cmp::Reverse((index, job))) = heap.pop() {
        if index != processed {
            return Err(Error::Invalid);
        }
        processed += 1;
        let record = heads[job].take().ok_or(Error::Invalid)?;
        heads[job] = streams[job].next();
        if let Some(next) = &heads[job] {
            heap.push(std::cmp::Reverse((next.index, job)));
        }
        if record.cached {
            if compensate && accepted.len() < p.limits.compensation_cells {
                work.charge(4)?;
                accepted.push(record.at);
            }
            continue;
        }
        if compensate
            && record.verdict == Some(CellCheck::Normal)
            && accepted.len() < p.limits.compensation_cells
        {
            work.charge(4)?;
            accepted.push(record.at);
        }
        if record.verdict == Some(CellCheck::Match) {
            if compensate && accepted.len() < p.limits.compensation_cells {
                work.charge(4)?;
                accepted.push(record.at);
            }
            work.charge(128)?;
            if p.certificates.len() < p.limits.certificates {
                p.certificates.insert(record.key);
            }
        } else {
            changes.mark(p, grid, record.at, record.axis.ok_or(Error::Invalid)?, work)?;
        }
    }
    if processed != expected {
        return Err(Error::Invalid);
    }
    Ok(())
}

fn record_buffers(
    jobs: &[Job],
    mut bytes: usize,
    limit: usize,
) -> Option<(Vec<Option<Vec<Record>>>, usize)> {
    let cells = jobs.iter().try_fold(0usize, |sum, j| {
        sum.checked_add((j.u[1] - j.u[0]).checked_mul(j.v[1] - j.v[0])?)
    })?;
    let mut preview = bytes;
    account_records(&mut preview, cells, limit)?;
    let mut buffers = Vec::with_capacity(jobs.len());
    for job in jobs {
        let count = (job.u[1] - job.u[0]).checked_mul(job.v[1] - job.v[0])?;
        buffers.push(Some(record_buffer(count, &mut bytes, limit)?));
    }
    Some((buffers, bytes))
}
