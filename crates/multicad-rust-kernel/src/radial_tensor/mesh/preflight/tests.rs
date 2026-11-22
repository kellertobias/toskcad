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
        p.probe_normal(
            Some(&candidate),
            &g,
            &ring,
            [cone, 0, 0],
            0,
            &mut BTreeMap::new(),
            &mut w
        )
        .unwrap(),
        None
    );
    let cell = ring.offsets[side];
    candidate.indices.swap(6 * cell, 6 * cell + 1);
    assert!(
        p.probe_normal(
            Some(&candidate),
            &g,
            &ring,
            [side, 0, 0],
            0,
            &mut BTreeMap::new(),
            &mut w
        )
        .unwrap()
        .is_some()
    );
    let mut changes = Refinement::new(g.u.len());
    assert!(
        p.normal_preflight(Some(&candidate), &g, &ring, 0, &mut changes, &mut w)
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
            .normal_preflight(Some(&candidate), &g, &ring, 0, &mut change, &mut w)
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
                Some(&candidate),
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
        Some(&candidate),
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
                Some(&candidate),
                &g,
                &ring,
                0,
                &mut Refinement::new(g.u.len()),
                &mut cancelled
            )
            .err(),
            Some(Error::Cancelled)
        );
        assert_eq!(
            count.get(),
            at,
            "No planning work may follow observed cancellation"
        );
    }
    let first = candidate.indices[6 * cell] as usize;
    for k in 0..3 {
        candidate.vertices[3 * first + k] =
            candidate.vertices[3 * candidate.indices[6 * cell + 1] as usize + k];
    }
    let _ = p
        .probe_normal(
            Some(&candidate),
            &g,
            &ring,
            [side, 0, 0],
            0,
            &mut BTreeMap::new(),
            &mut w,
        )
        .unwrap();
    candidate.vertices[3 * first] = f64::INFINITY;
    let _ = p
        .probe_normal(
            Some(&candidate),
            &g,
            &ring,
            [side, 0, 0],
            0,
            &mut BTreeMap::new(),
            &mut w,
        )
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

#[test]
fn on_demand_tensor_cells_match_emitted_caps_rails_seams_and_diagonals_bitwise() {
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
            .map(|_| vec![0., 0.125, 0.5, 0.875, 1.])
            .collect(),
        v: vec![0., 0.125, 0.5, 0.875, 1.],
    };
    let ring = Ring::new(&g);
    let mesh = p.emit(&g, &ring, 0.001, 50000, &mut w).unwrap();
    let mut points = BTreeMap::new();
    let mut count = 0;
    for side in 0..p.sides.len() {
        if !matches!(p.sides[side], Side::Tensor(_)) {
            continue;
        }
        for row in 0..g.v.len() - 1 {
            for u in 0..g.u[side].len() - 1 {
                let a = p
                    .stencil_cell(None, &g, &ring, [side, u, row], &mut points, &mut w)
                    .unwrap();
                let b = p
                    .stencil_cell(Some(&mesh), &g, &ring, [side, u, row], &mut points, &mut w)
                    .unwrap();
                assert_eq!(a.range, b.range);
                assert_eq!(
                    a.vertices.map(|p| p.map(f64::to_bits)),
                    b.vertices.map(|p| p.map(f64::to_bits))
                );
                assert_eq!(a.actual, b.actual);
                count += 1;
            }
        }
    }
    assert!(count >= 16);
    assert!(!points.is_empty());
}

#[test]
fn width_cohorts_group_rounded_grid_steps_only_for_planning() {
    let a = [0., 1. / 568.];
    let b = [567. / 568., 1.];
    assert_ne!(a[1] - a[0], b[1] - b[0]);
    assert!(same_or_coarser(a, b));
    assert!(same_or_coarser(b, a));
    assert!(same_or_coarser([0., 0.5], a));
    assert!(!same_or_coarser([0., 0.5 / 568.], a));
}

#[test]
fn an_actual_unsampled_flipped_facet_requests_recovery_in_both_parameter_directions() {
    for transverse in [false, true] {
        let body = if transverse {
            super::super::parallel_probe::multi_u_cam()
        } else {
            super::super::tests::cam()
        };
        let original = body.input().clone();
        let mut w = work();
        let mut p = Prepared::new(
            &body,
            Tessellation {
                deflection: 10.,
                angle: 3.,
                triangle_limit: 50000,
            },
            DEFAULT_LIMITS,
            &mut w,
        )
        .unwrap();
        let n = if transverse { 4 } else { 32 };
        let mut v = (0..=16).map(|i| i as f64 / 16.).collect::<Vec<_>>();
        v[4] += 1e-4;
        let g = Grid {
            u: p.sides
                .iter()
                .map(|_| (0..=n).map(|i| i as f64 / n as f64).collect())
                .collect(),
            v,
        };
        let ring = Ring::new(&g);
        let mut mesh = p.emit(&g, &ring, 0.001, 50000, &mut w).unwrap();
        let side = p
            .sides
            .iter()
            .position(|s| matches!(s, Side::Tensor(_)))
            .unwrap();
        let u = 3 * n / 8;
        let row = 3;
        let cell = row * ring.entries.len() + ring.offsets[side] + u;
        // This literal facet is outside the full-certificate stencil's three
        // complete rows and the source-knot rows' boundary columns.
        mesh.indices.swap(6 * cell, 6 * cell + 1);
        let range = [[g.u[side][u], g.u[side][u + 1]], [g.v[row], g.v[row + 1]]];
        assert_eq!(
            p.hint(side, range, 0, &mut w).unwrap(),
            usize::from(!transverse)
        );
        let mut changes = Refinement::new(g.u.len());
        assert!(
            p.normal_preflight(Some(&mesh), &g, &ring, 0, &mut changes, &mut w)
                .unwrap()
        );
        if transverse {
            assert!(!changes.u[side].is_empty());
            assert!(changes.v.is_empty());
            // Exercise cancellation in both the late approximate probe and
            // the source-failure cohort planner, after literal cache invalidation.
            for selected_source_failure in [false, true] {
                if selected_source_failure {
                    let selected = ring.offsets[side];
                    mesh.indices.swap(6 * selected, 6 * selected + 1);
                    let mut changed = Refinement::new(g.u.len());
                    assert!(
                        p.normal_preflight(Some(&mesh), &g, &ring, 0, &mut changed, &mut w)
                            .unwrap()
                    );
                    assert!(!changed.u[side].is_empty());
                }
                let observed = std::cell::Cell::new(0);
                let observe = || {
                    observed.set(observed.get() + 1);
                    false
                };
                p.normal_preflight(
                    Some(&mesh),
                    &g,
                    &ring,
                    0,
                    &mut Refinement::new(g.u.len()),
                    &mut Work {
                        remaining: 100_000_000,
                        interrupt: Some(&observe),
                    },
                )
                .unwrap();
                for at in 1..=observed.get() {
                    let calls = std::cell::Cell::new(0);
                    let cancel = || {
                        calls.set(calls.get() + 1);
                        calls.get() >= at
                    };
                    assert_eq!(
                        p.normal_preflight(
                            Some(&mesh),
                            &g,
                            &ring,
                            0,
                            &mut Refinement::new(g.u.len()),
                            &mut Work {
                                remaining: 100_000_000,
                                interrupt: Some(&cancel)
                            }
                        )
                        .err(),
                        Some(Error::Cancelled)
                    );
                    assert_eq!(calls.get(), at);
                }
            }
        } else {
            assert!(!changes.v.is_empty());
        }
        assert_eq!(body.input(), &original);
    }
}

#[test]
fn bounded_source_row_selection_and_accepted_skew_cells_keep_planning_resource_limits() {
    let body = super::super::parallel_probe::multi_u_cam();
    let mut w = work();
    let mut p = Prepared::new(
        &body,
        Tessellation {
            deflection: 10.,
            angle: 3.,
            triangle_limit: 50000,
        },
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let mut v = (0..=16).map(|i| i as f64 / 16.).collect::<Vec<_>>();
    v[4] += 1e-4;
    let g = Grid {
        u: p.sides.iter().map(|_| vec![0., 0.5, 1.]).collect(),
        v,
    };
    let ring = Ring::new(&g);
    let mesh = p.emit(&g, &ring, 0.001, 50000, &mut w).unwrap();
    p.limits.probe_cells = 4;
    let _ = p
        .normal_preflight(
            Some(&mesh),
            &g,
            &ring,
            0,
            &mut Refinement::new(g.u.len()),
            &mut w,
        )
        .unwrap();
    let g = Grid {
        u: p.sides.iter().map(|_| vec![0., 1.]).collect(),
        v: vec![0., 0.001, 1.],
    };
    let ring = Ring::new(&g);
    let mesh = p.emit(&g, &ring, 0.001, 50000, &mut w).unwrap();
    let cell = p
        .stencil_cell(
            Some(&mesh),
            &g,
            &ring,
            [0, 0, 0],
            &mut BTreeMap::new(),
            &mut w,
        )
        .unwrap();
    assert!(p.cell_check(cell, None, None, &mut w).unwrap() == CellCheck::Match);
    let mut changes = Refinement::new(g.u.len());
    changes.v.insert(0.0005f64.to_bits());
    changes
        .compensate(&p, &g, &mesh, &ring, &[[0, 0, 0]], &mut w)
        .unwrap();
    assert!(!changes.u[0].is_empty());
}

#[test]
fn source_transverse_recovery_leaves_strictly_finer_columns_unchanged() {
    let body = super::super::parallel_probe::multi_u_cam();
    let original = body.input().clone();
    let mut w = work();
    let mut p = Prepared::new(
        &body,
        Tessellation {
            deflection: 10.,
            angle: 3.,
            triangle_limit: 50000,
        },
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let g = Grid {
        u: p.sides
            .iter()
            .map(|_| vec![0., 0.125, 0.25, 0.5, 0.75, 1.])
            .collect(),
        v: vec![0., 0.25, 0.5, 0.75, 1.],
    };
    let ring = Ring::new(&g);
    let mut mesh = p.emit(&g, &ring, 0.001, 50000, &mut w).unwrap();
    let selected = ring.offsets[0] + 2;
    mesh.indices.swap(6 * selected, 6 * selected + 1);
    let mut changes = Refinement::new(g.u.len());
    assert!(
        p.normal_preflight(Some(&mesh), &g, &ring, 0, &mut changes, &mut w)
            .unwrap()
    );
    assert!(!changes.u[0].is_empty());
    assert!(changes.u[0].iter().all(|&bits| f64::from_bits(bits) > 0.25));
    assert!(changes.v.is_empty());
    assert_eq!(body.input(), &original);
}
