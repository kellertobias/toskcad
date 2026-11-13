use super::*;
fn options() -> Tessellation {
    Tessellation {
        deflection: 0.025,
        angle: 0.15,
        triangle_limit: 100_000,
    }
}
fn grid(sides: usize) -> Grid {
    Grid {
        u: vec![vec![0., 0.5, 1.]; sides],
        v: vec![0., 0.5, 1.],
    }
}
#[test]
fn source_owner_and_record_memory_preflight_refuse_before_allocation() {
    let mut owners = BTreeSet::new();
    let key = (2, [[0, 1], [2, 3]]);
    assert_eq!(register_source(&mut owners, key, 0), Ok(false));
    assert_eq!(register_source(&mut owners, key, 2), Ok(true));
    assert_eq!(register_source(&mut owners, key, 2), Err(Error::Invalid));
    let impossible_job = Job {
        side: 0,
        u: [0, usize::MAX],
        v: [0, 2],
        first: 0,
        source_keys: 0,
    };
    assert!(record_buffers(&[impossible_job], 0, usize::MAX).is_none());
    let mut bytes = 0;
    assert!(record_buffer(1, &mut bytes, 0).is_none());
    assert!(record_buffer(usize::MAX, &mut bytes, usize::MAX).is_none());
    let impossible = isize::MAX as usize / std::mem::size_of::<Record>() + 1;
    assert!(record_buffer(impossible, &mut bytes, usize::MAX).is_none());
    let records = record_buffer(3, &mut bytes, 10_000).unwrap();
    assert_eq!(bytes, records.capacity() * std::mem::size_of::<Record>());
    assert!(account_records(&mut bytes, 10000, 0).is_none());
}
#[test]
fn actual_cam_resource_caps_and_mixed_source_keep_serial_fallback_honest() {
    let body = super::super::super::parallel_probe::multi_u_cam();
    for limits in [
        Limits {
            parallel_jobs: 0,
            ..DEFAULT_LIMITS
        },
        Limits {
            parallel_sources: 0,
            ..DEFAULT_LIMITS
        },
        Limits {
            parallel_bytes: 0,
            ..DEFAULT_LIMITS
        },
    ] {
        let mut w = Work {
            remaining: DEFAULT_LIMITS.work,
            interrupt: None,
        };
        let mut p = Prepared::new(&body, options(), limits, &mut w).unwrap();
        let grid = grid(p.sides.len());
        let ring = Ring::new(&grid);
        let mesh = p.emit(&grid, &ring, 0., 100_000, &mut w).unwrap();
        let sides = p.sides.len();
        assert!(
            !check(
                &mut p,
                &mesh,
                &grid,
                &ring,
                0,
                &mut Refinement::new(sides),
                &mut Vec::new(),
                8,
                &mut w
            )
            .unwrap()
        );
    }
    let mixed = super::super::super::tests::cam();
    let mut w = Work {
        remaining: DEFAULT_LIMITS.work,
        interrupt: None,
    };
    let mut p = Prepared::new(&mixed, options(), DEFAULT_LIMITS, &mut w).unwrap();
    let grid = grid(p.sides.len());
    let ring = Ring::new(&grid);
    assert!(matches!(
        jobs(&p, &grid, &ring, &mut w),
        Err(Error::Unsupported)
    ));
    let mesh = p.emit(&grid, &ring, 0., 100_000, &mut w).unwrap();
    let side = p
        .sides
        .iter()
        .position(|s| matches!(s, Side::Cone(_)))
        .unwrap();
    let job = Job {
        side,
        u: [0, 1],
        v: [0, 1],
        first: 0,
        source_keys: 0,
    };
    assert!(matches!(
        run_job(
            &p,
            &mesh,
            &grid,
            &ring,
            0,
            job,
            &Cell::new(0),
            Vec::new(),
            0,
            &mut w
        ),
        Err((0, Error::Unsupported))
    ));
}

fn completed_records(
    p: &Prepared<'_>,
    mesh: &DisplayMesh,
    grid: &Grid,
    ring: &Ring,
    w: &mut Work<'_>,
) -> Vec<Done> {
    jobs(p, grid, ring, w)
        .unwrap()
        .unwrap()
        .into_iter()
        .enumerate()
        .map(|(i, job)| Done {
            job: i,
            worker: 0,
            elapsed: Duration::ZERO,
            result: run_job(
                p,
                mesh,
                grid,
                ring,
                0,
                job,
                &Cell::new(job.first),
                Vec::new(),
                128,
                w,
            ),
        })
        .collect()
}
#[test]
fn actual_cell_degeneracy_and_incomplete_or_reordered_receipts_fail_closed() {
    let body = super::super::super::parallel_probe::multi_u_cam();
    let original = body.input().clone();
    let mut w = Work {
        remaining: DEFAULT_LIMITS.work,
        interrupt: None,
    };
    let mut p = Prepared::new(&body, options(), DEFAULT_LIMITS, &mut w).unwrap();
    let grid = grid(p.sides.len());
    let mut ring = Ring::new(&grid);
    ring.entries.pop();
    assert!(matches!(
        jobs(&p, &grid, &ring, &mut w),
        Err(Error::Invalid)
    ));
    let ring = Ring::new(&grid);
    let mesh = p.emit(&grid, &ring, 0., 100_000, &mut w).unwrap();
    let planned = jobs(&p, &grid, &ring, &mut w).unwrap().unwrap();
    let job = planned[0];
    let mut degenerate = mesh.clone();
    degenerate.indices[6 * job.first + 1] = degenerate.indices[6 * job.first];
    let records = run_job(
        &p,
        &degenerate,
        &grid,
        &ring,
        0,
        job,
        &Cell::new(job.first),
        Vec::new(),
        128,
        &mut w,
    )
    .unwrap();
    assert!(records.iter().any(|r| r.verdict.is_none()));
    let done = completed_records(&p, &mesh, &grid, &ring, &mut w);
    assert_eq!(
        verify_completed(&done, done.len() + 1),
        Err(Error::Unresolved)
    );
    let cancelled = Done {
        job: 0,
        worker: 0,
        elapsed: Duration::ZERO,
        result: Err((0, Error::Cancelled)),
    };
    let sides = p.sides.len();
    assert_eq!(
        apply_completed(
            &mut p,
            &grid,
            &ring,
            &mut Refinement::new(sides),
            &mut Vec::new(),
            vec![cancelled],
            &mut w
        ),
        Err(Error::Cancelled)
    );
    for remove_last in [false, true] {
        let mut done = completed_records(&p, &mesh, &grid, &ring, &mut w);
        let target = if remove_last {
            (grid.v.len() - 1) * ring.entries.len() - 1
        } else {
            0
        };
        for d in &mut done {
            let rows = d.result.as_mut().unwrap();
            if let Some(i) = rows.iter().position(|r| r.index == target) {
                rows.remove(i);
                break;
            }
        }
        assert_eq!(
            apply_completed(
                &mut p,
                &grid,
                &ring,
                &mut Refinement::new(sides),
                &mut Vec::new(),
                done,
                &mut w
            ),
            Err(Error::Invalid)
        );
    }
    let fatal = vec![
        Done {
            job: 1,
            worker: 0,
            elapsed: Duration::ZERO,
            result: Err((4, Error::Unresolved)),
        },
        Done {
            job: 0,
            worker: 0,
            elapsed: Duration::ZERO,
            result: Err((2, Error::Limit)),
        },
    ];
    assert_eq!(verify_completed(&fatal, 2), Err(Error::Limit));
    assert_eq!(body.input(), &original);
}
#[test]
fn poisoned_worker_buffer_and_closed_channel_return_all_unused_work() {
    let body = super::super::super::parallel_probe::multi_u_cam();
    let mut w = Work {
        remaining: DEFAULT_LIMITS.work,
        interrupt: None,
    };
    let mut p = Prepared::new(&body, options(), DEFAULT_LIMITS, &mut w).unwrap();
    let grid = grid(p.sides.len());
    let ring = Ring::new(&grid);
    let mesh = p.emit(&grid, &ring, 0., 100_000, &mut w).unwrap();
    let planned = jobs(&p, &grid, &ring, &mut w).unwrap().unwrap();
    let jobs = &planned[..1];
    let poison = Mutex::new(vec![Some(Vec::new())]);
    let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _guard = poison.lock().unwrap();
        panic!("injected buffer transport failure");
    }));
    let pool = credits::Pool::new(1_000_000);
    let (tx, rx) = mpsc::channel();
    worker_loop(
        0,
        &p,
        &mesh,
        &grid,
        &ring,
        0,
        jobs,
        &AtomicUsize::new(0),
        &poison,
        tx,
        Arc::clone(&pool),
        &vec![128; p.sides.len()],
    );
    assert!(matches!(
        rx.recv().unwrap().result,
        Err((_, Error::Unresolved))
    ));
    assert_eq!(pool.finish().unwrap().0, 1_000_000);
    let buffers = Mutex::new(vec![Some(Vec::new())]);
    let (tx, rx) = mpsc::channel();
    drop(rx);
    worker_loop(
        0,
        &p,
        &mesh,
        &grid,
        &ring,
        0,
        jobs,
        &AtomicUsize::new(0),
        &buffers,
        tx,
        Arc::clone(&pool),
        &vec![128; p.sides.len()],
    );
    assert!(pool.finish().unwrap().0 < 1_000_000);
}
#[test]
fn native_thread_refusal_and_join_failure_cannot_publish_geometry() {
    let pool = credits::Pool::new(100);
    thread::scope(|scope| {
        let handles = vec![scope.spawn(|| panic!("injected worker transport panic"))];
        assert!(join_workers(handles, &pool));
    });
    assert!(pool.cancel.load(Ordering::Acquire));
    assert_eq!(pool.finish().unwrap().0, 100);
    let body = super::super::super::parallel_probe::multi_u_cam();
    let original = body.input().clone();
    let mut w = Work {
        remaining: DEFAULT_LIMITS.work,
        interrupt: None,
    };
    let mut p = Prepared::new(
        &body,
        options(),
        Limits {
            parallel_stack: usize::MAX,
            ..DEFAULT_LIMITS
        },
        &mut w,
    )
    .unwrap();
    let grid = grid(p.sides.len());
    let ring = Ring::new(&grid);
    let mesh = p.emit(&grid, &ring, 0., 100_000, &mut w).unwrap();
    let sides = p.sides.len();
    assert_eq!(
        check(
            &mut p,
            &mesh,
            &grid,
            &ring,
            0,
            &mut Refinement::new(sides),
            &mut Vec::new(),
            8,
            &mut w
        ),
        Err(Error::Limit)
    );
    assert_eq!(body.input(), &original);
}
#[test]
fn coordinator_cancellation_after_worker_completion_still_refuses_publication() {
    let body = super::super::super::parallel_probe::multi_u_cam();
    let original = body.input().clone();
    let saw = Cell::new(false);
    let coordinator = thread::current().id();
    let interrupt = || {
        assert_eq!(thread::current().id(), coordinator);
        let active = super::super::super::parallel_probe::ACTIVE.with(Cell::get);
        if active {
            saw.set(true);
            false
        } else {
            saw.get()
        }
    };
    let mut w = Work {
        remaining: DEFAULT_LIMITS.work,
        interrupt: Some(&interrupt),
    };
    let mut p = Prepared::new(&body, options(), DEFAULT_LIMITS, &mut w).unwrap();
    let grid = grid(p.sides.len());
    let ring = Ring::new(&grid);
    let mesh = p.emit(&grid, &ring, 0., 100_000, &mut w).unwrap();
    let sides = p.sides.len();
    assert_eq!(
        check(
            &mut p,
            &mesh,
            &grid,
            &ring,
            0,
            &mut Refinement::new(sides),
            &mut Vec::new(),
            8,
            &mut w
        ),
        Err(Error::Cancelled)
    );
    assert!(saw.get());
    assert_eq!(body.input(), &original);
}
