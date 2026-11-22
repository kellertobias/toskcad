//! Private batched work ledger; only rare refills synchronize worker budgets.
use super::*;
use std::cell::{Cell, RefCell};
use std::{
    rc::Rc,
    sync::{
        Arc, Condvar, Mutex,
        atomic::{AtomicBool, AtomicUsize, Ordering},
    },
};
const BATCH: usize = 262_144;
struct State {
    available: usize,
    leased: usize,
    spent: usize,
    refills: usize,
}
pub(super) struct Pool {
    initial: usize,
    #[cfg(test)]
    sync_ns: std::sync::atomic::AtomicU64,
    #[cfg(test)]
    waits: AtomicUsize,
    state: Mutex<State>,
    wake: Condvar,
    pub cancel: AtomicBool,
    pub stop: AtomicUsize,
}
impl Pool {
    pub fn new(initial: usize) -> Arc<Self> {
        Arc::new(Self {
            initial,
            #[cfg(test)]
            sync_ns: std::sync::atomic::AtomicU64::new(0),
            #[cfg(test)]
            waits: AtomicUsize::new(0),
            state: Mutex::new(State {
                available: initial,
                leased: 0,
                spent: 0,
                refills: 0,
            }),
            wake: Condvar::new(),
            cancel: AtomicBool::new(false),
            stop: AtomicUsize::new(usize::MAX),
        })
    }
    pub fn interrupted(&self, current: usize) -> bool {
        self.cancel.load(Ordering::Acquire) || current > self.stop.load(Ordering::Acquire)
    }
    fn take(&self, needed: usize, current: usize) -> Result<usize, Error> {
        #[cfg(test)]
        let _timing = Timing {
            at: std::time::Instant::now(),
            target: &self.sync_ns,
        };
        let mut s = self.state.lock().map_err(|_| Error::Unresolved)?;
        loop {
            if self.interrupted(current) {
                return Err(Error::Cancelled);
            }
            if s.available >= needed {
                let amount = s.available.min(BATCH.max(needed));
                s.available -= amount;
                s.leased += 1;
                s.refills += 1;
                return Ok(amount);
            }
            if s.leased == 0 {
                return Err(Error::Limit);
            }
            #[cfg(test)]
            self.waits.fetch_add(1, Ordering::Relaxed);
            s = self
                .wake
                .wait_timeout(s, std::time::Duration::from_millis(1))
                .map_err(|_| Error::Unresolved)?
                .0;
        }
    }
    fn release(&self, grant: usize, balance: usize) -> Result<(), Error> {
        #[cfg(test)]
        let _timing = Timing {
            at: std::time::Instant::now(),
            target: &self.sync_ns,
        };
        if grant == 0 {
            return if balance == 0 {
                Ok(())
            } else {
                Err(Error::Invalid)
            };
        };
        let mut s = self.state.lock().map_err(|_| Error::Unresolved)?;
        if balance > grant || s.leased == 0 {
            return Err(Error::Invalid);
        }
        s.leased -= 1;
        s.available = s.available.checked_add(balance).ok_or(Error::Limit)?;
        s.spent = s.spent.checked_add(grant - balance).ok_or(Error::Limit)?;
        drop(s);
        self.wake.notify_all();
        Ok(())
    }
    pub fn halt_after(&self, index: usize) {
        self.stop.fetch_min(index, Ordering::AcqRel);
        self.wake.notify_all();
    }
    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Release);
        self.wake.notify_all();
    }
    #[cfg(test)]
    pub fn timing(&self) -> (std::time::Duration, usize) {
        (
            std::time::Duration::from_nanos(self.sync_ns.load(Ordering::Relaxed)),
            self.waits.load(Ordering::Relaxed),
        )
    }
    pub fn finish(&self) -> Result<(usize, usize, usize), Error> {
        let s = self.state.lock().map_err(|_| Error::Unresolved)?;
        if s.leased != 0 || s.available.checked_add(s.spent) != Some(self.initial) {
            return Err(Error::Invalid);
        }
        Ok((s.available, s.spent, s.refills))
    }
}
#[cfg(test)]
struct Timing<'a> {
    at: std::time::Instant,
    target: &'a std::sync::atomic::AtomicU64,
}
#[cfg(test)]
impl Drop for Timing<'_> {
    fn drop(&mut self) {
        self.target.fetch_add(
            self.at.elapsed().as_nanos().min(u128::from(u64::MAX)) as u64,
            Ordering::Relaxed,
        );
    }
}
struct Local {
    pool: Arc<Pool>,
    grant: Cell<usize>,
    current: Rc<Cell<usize>>,
}
thread_local! {static LOCAL:RefCell<Option<Rc<Local>>>=const{RefCell::new(None)};}
/// An explicitly bounded nested query must neither refill from, nor return
/// its unrelated balance to, an active worker's lease. Restore on panic too.
pub(super) fn local_budget<T>(f: impl FnOnce() -> T) -> T {
    struct Restore(Option<Rc<Local>>);
    impl Drop for Restore {
        fn drop(&mut self) {
            LOCAL.with(|slot| slot.replace(self.0.take()));
        }
    }
    let _restore = Restore(LOCAL.with(|slot| slot.replace(None)));
    f()
}
pub(super) fn refill(balance: &mut usize, needed: usize) -> Result<(), Error> {
    LOCAL.with(|slot| {
        let local = slot.borrow().clone().ok_or(Error::Limit)?;
        let grant = local.grant.replace(0);
        local.pool.release(grant, std::mem::take(balance))?;
        let grant = local.pool.take(needed, local.current.get())?;
        local.grant.set(grant);
        *balance = grant;
        Ok(())
    })
}
pub(super) fn run<T>(
    pool: Arc<Pool>,
    current: Rc<Cell<usize>>,
    f: impl FnOnce(&mut Work<'_>) -> Result<T, Error>,
) -> Result<T, Error> {
    let local = Rc::new(Local {
        pool: Arc::clone(&pool),
        grant: Cell::new(0),
        current: Rc::clone(&current),
    });
    let old = LOCAL.with(|slot| slot.replace(Some(Rc::clone(&local))));
    let interrupt = || pool.interrupted(current.get());
    let mut work = Work {
        remaining: 0,
        interrupt: Some(&interrupt),
    };
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| f(&mut work)))
        .map_err(|_| Error::Unresolved)
        .and_then(|r| r);
    let release = pool.release(local.grant.replace(0), work.remaining);
    LOCAL.with(|slot| {
        slot.replace(old);
    });
    release?;
    result
}

#[test]
fn ledger_refunds_unused_credits_on_success_error_and_panic() {
    let pool = Pool::new(BATCH * 4);
    for failure in 0..3 {
        let result = run(Arc::clone(&pool), Rc::new(Cell::new(0)), |w| {
            w.charge(BATCH + 11)?;
            w.charge(17)?;
            match failure {
                0 => Ok(()),
                1 => Err(Error::Unsupported),
                _ => panic!("injected worker failure"),
            }
        });
        assert_eq!(
            result,
            match failure {
                0 => Ok(()),
                1 => Err(Error::Unsupported),
                _ => Err(Error::Unresolved),
            }
        );
        let (remaining, spent, _) = pool.finish().unwrap();
        assert_eq!(spent, (failure + 1) * (BATCH + 28));
        assert_eq!(remaining + spent, BATCH * 4);
        assert!(LOCAL.with(|s| s.borrow().is_none()));
    }
}

#[test]
fn independent_source_queries_cannot_refill_or_refund_a_workers_lease() {
    let source = crate::radial_tensor::tests::material();
    let bounds = [[0.1, 0.1], [0.2, 0.2], [4., 4.]];
    let query_cost = source
        .classify_box(bounds, Default::default(), None)
        .unwrap()
        .work_used();
    let initial = BATCH * 4;
    let pool = Pool::new(initial);
    run(Arc::clone(&pool), Rc::new(Cell::new(0)), |w| {
        w.charge(7)?;
        let before = w.remaining;
        assert_eq!(
            source
                .classify_box(
                    bounds,
                    crate::radial_tensor::BoxQueryOptions {
                        maximum_work: query_cost - 1
                    },
                    None
                )
                .err(),
            Some(Error::Limit)
        );
        assert_eq!(w.remaining, before);
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| local_budget(|| panic!(
                "independent query callback"
            ))))
            .is_err()
        );
        assert_eq!(
            source
                .classify_box(bounds, Default::default(), None)
                .unwrap()
                .work_used(),
            query_cost
        );
        w.charge(query_cost)?;
        w.charge(BATCH)?;
        Ok(())
    })
    .unwrap();
    let (available, spent, refills) = pool.finish().unwrap();
    assert_eq!(spent, 7 + query_cost + BATCH);
    assert_eq!(available + spent, initial);
    assert_eq!(refills, 2);
    let (_, waits) = pool.timing();
    assert_eq!(
        waits, 0,
        "independent queries must not wait on worker credits"
    );
}

#[test]
fn ledger_waits_for_returnable_credit_and_obeys_terminal_stop() {
    let pool = Pool::new(100);
    let grant = pool.take(70, 0).unwrap();
    assert_eq!(grant, 100);
    assert_eq!(pool.finish(), Err(Error::Invalid));
    let other = Arc::clone(&pool);
    let (tx, rx) = std::sync::mpsc::channel();
    let worker = std::thread::spawn(move || {
        tx.send(()).unwrap();
        let grant = other.take(30, 1).unwrap();
        other.release(grant, grant - 30).unwrap();
    });
    rx.recv().unwrap();
    pool.release(grant, 40).unwrap();
    worker.join().unwrap();
    assert_eq!(pool.finish().unwrap().0, 10);
    assert_eq!(pool.take(11, 0), Err(Error::Limit));
    pool.halt_after(4);
    assert_eq!(pool.take(1, 5), Err(Error::Cancelled));
    let remaining = pool.take(1, 4).unwrap();
    pool.release(remaining, remaining).unwrap();
    pool.cancel();
    assert_eq!(pool.take(1, 0), Err(Error::Cancelled));
    assert_eq!(pool.finish().unwrap().1, 90);
}

#[test]
fn cancelled_and_exhausted_worker_refund_without_publishing_success() {
    let pool = Pool::new(80);
    let result = run(Arc::clone(&pool), Rc::new(Cell::new(0)), |w| {
        w.charge(41)?;
        w.charge(40)
    });
    assert_eq!(result, Err(Error::Limit));
    assert_eq!(pool.finish().unwrap().0, 39);
    let result = run(Arc::clone(&pool), Rc::new(Cell::new(0)), |w| {
        w.charge(7)?;
        pool.cancel();
        w.step()
    });
    assert_eq!(result, Err(Error::Cancelled));
    assert_eq!(pool.finish().unwrap().0, 32);
    assert!(LOCAL.with(|s| s.borrow().is_none()));
    assert_eq!(refill(&mut 0, 1), Err(Error::Limit));
}

#[test]
fn invalid_credit_refunds_are_rejected_without_mutating_the_ledger() {
    let pool = Pool::new(100);
    assert_eq!(pool.release(0, 1), Err(Error::Invalid));
    assert_eq!(pool.release(1, 0), Err(Error::Invalid));
    let grant = pool.take(1, 0).unwrap();
    assert_eq!(pool.release(grant, grant + 1), Err(Error::Invalid));
    pool.release(grant, grant).unwrap();
    assert_eq!(pool.finish().unwrap().0, 100);
}
