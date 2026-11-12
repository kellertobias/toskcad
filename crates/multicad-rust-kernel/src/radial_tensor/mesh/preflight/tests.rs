use super::*;
fn work() -> Work<'static> {
    Work {
        remaining: 100_000_000,
        interrupt: None,
    }
}
fn grid(p: &Prepared<'_>) -> Grid {
    Grid {
        u: p.sides
            .iter()
            .map(|s| {
                if matches!(s, Side::Cone(_)) {
                    vec![0., 0.5, 1.]
                } else {
                    vec![0., 1.]
                }
            })
            .collect(),
        v: vec![0., 0.5, 1.],
    }
}
#[test]
fn actual_stencil_can_only_plan_refinement_and_refuses_uncancellable_or_unbounded_work() {
    let body = super::super::tests::cam();
    let original = body.input().clone();
    let mut w = work();
    let mut p = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let g = grid(&p);
    let ring = Ring::new(&g);
    let mut candidate = p.emit(&g, &ring, 0.001, 50000, &mut w).unwrap();
    let side = p
        .sides
        .iter()
        .position(|s| matches!(s, Side::Tensor(_)))
        .unwrap();
    let cone = p
        .sides
        .iter()
        .position(|s| matches!(s, Side::Cone(_)))
        .unwrap();
    assert_eq!(
        p.probe_normal(&candidate, &g, &ring, [cone, 0, 0], 0, &mut w)
            .unwrap(),
        None
    );
    let cell = ring.offsets[side];
    candidate.indices.swap(6 * cell, 6 * cell + 1);
    assert!(
        p.probe_normal(&candidate, &g, &ring, [side, 0, 0], 0, &mut w)
            .unwrap()
            .is_some()
    );
    let mut changes = Refinement::new(g.u.len());
    assert!(
        p.normal_preflight(&candidate, &g, &ring, 0, &mut changes, &mut w)
            .unwrap()
    );
    assert!(!changes.empty());
    let mut final_changes = Refinement::new(g.u.len());
    p.check(&candidate, &g, &ring, 0, &mut final_changes, &mut w)
        .unwrap();
    assert!(!final_changes.empty());
    for limit in [0, 1] {
        p.limits.probe_cells = limit;
        p.limits.compensation_cells = limit;
        let mut change = Refinement::new(g.u.len());
        let _ = p
            .normal_preflight(&candidate, &g, &ring, 0, &mut change, &mut w)
            .unwrap();
    }
    p.limits = DEFAULT_LIMITS;
    for budget in [0, 100] {
        let mut low = Work {
            remaining: budget,
            interrupt: None,
        };
        assert_eq!(
            p.normal_preflight(
                &candidate,
                &g,
                &ring,
                0,
                &mut Refinement::new(g.u.len()),
                &mut low
            )
            .err(),
            Some(Error::Limit)
        );
    }
    let observed = std::cell::Cell::new(0);
    let observe = || {
        observed.set(observed.get() + 1);
        false
    };
    let mut observed_work = Work {
        remaining: 100_000_000,
        interrupt: Some(&observe),
    };
    p.normal_preflight(
        &candidate,
        &g,
        &ring,
        0,
        &mut Refinement::new(g.u.len()),
        &mut observed_work,
    )
    .unwrap();
    assert!(observed.get() > 2);
    for at in 1..=observed.get() {
        let count = std::cell::Cell::new(0);
        let cancel = || {
            count.set(count.get() + 1);
            count.get() >= at
        };
        let mut cancelled = Work {
            remaining: 100_000_000,
            interrupt: Some(&cancel),
        };
        assert_eq!(
            p.normal_preflight(
                &candidate,
                &g,
                &ring,
                0,
                &mut Refinement::new(g.u.len()),
                &mut cancelled
            )
            .err(),
            Some(Error::Cancelled)
        );
    }
    let first = candidate.indices[6 * cell] as usize;
    for k in 0..3 {
        candidate.vertices[3 * first + k] =
            candidate.vertices[3 * candidate.indices[6 * cell + 1] as usize + k];
    }
    let _ = p
        .probe_normal(&candidate, &g, &ring, [side, 0, 0], 0, &mut w)
        .unwrap();
    candidate.vertices[3 * first] = f64::INFINITY;
    let _ = p
        .probe_normal(&candidate, &g, &ring, [side, 0, 0], 0, &mut w)
        .unwrap();
    assert_eq!(body.input(), &original);
}

#[test]
fn a_collapsed_parameter_candidate_cannot_be_cached_as_certified() {
    let body = super::super::tests::cam();
    let mut w = work();
    let mut p = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let g = Grid {
        u: p.sides
            .iter()
            .map(|s| {
                if matches!(s, Side::Cone(_)) {
                    vec![0., 1e-200, 0.5, 1.]
                } else {
                    vec![0., 1e-200, 1.]
                }
            })
            .collect(),
        v: vec![0., 1e-200, 1.],
    };
    let ring = Ring::new(&g);
    let candidate = p.emit(&g, &ring, 0.001, 50000, &mut w).unwrap();
    let mut changes = Refinement::new(g.u.len());
    p.check(&candidate, &g, &ring, 0, &mut changes, &mut w)
        .unwrap();
    assert!(!changes.empty());
    assert!(
        !p.certificates
            .iter()
            .any(|key| key[2] == 1e-200_f64.to_bits() && key[4] == 1e-200_f64.to_bits())
    );
}
