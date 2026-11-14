use super::*;
#[test]
fn immutable_source_jet_cache_refuses_partial_construction_and_retries_cleanly() {
    let body = super::super::tests::cam();
    let mut w = Work {
        remaining: 10_000_000,
        interrupt: None,
    };
    let prepared = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let tensors = prepared
        .sides
        .iter()
        .filter_map(|s| {
            if let Side::Tensor(t) = s {
                Some(t)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    let patch = &tensors[0].cells[0];
    let snapshot = patch
        .normal
        .iter()
        .flatten()
        .flatten()
        .map(|v| [v.lo, v.hi])
        .collect::<Vec<_>>();
    assert!(patch.normal_jet.get().is_none());
    assert_eq!(
        patch
            .taylor_normal_box(
                patch.domain,
                &mut Work {
                    remaining: 0,
                    interrupt: None
                }
            )
            .err(),
        Some(Error::Limit)
    );
    assert!(patch.normal_jet.get().is_none());
    assert_eq!(
        patch
            .taylor_normal_box(
                patch.domain,
                &mut Work {
                    remaining: 1_000_000,
                    interrupt: Some(&|| true)
                }
            )
            .err(),
        Some(Error::Cancelled)
    );
    assert!(patch.normal_jet.get().is_none());
    let bound = patch
        .taylor_normal_box(patch.domain, &mut w)
        .unwrap()
        .map(|v| [v.lo, v.hi]);
    assert!(patch.normal_jet.get().unwrap().is_some());
    assert_eq!(
        patch
            .taylor_normal_box(
                patch.domain,
                &mut Work {
                    remaining: 64,
                    interrupt: Some(&|| true)
                }
            )
            .err(),
        Some(Error::Cancelled)
    );
    let warm = patch
        .taylor_normal_box(
            patch.domain,
            &mut Work {
                remaining: 64,
                interrupt: None,
            },
        )
        .unwrap()
        .map(|v| [v.lo, v.hi]);
    assert_eq!(warm, bound);
    assert_eq!(
        patch
            .normal
            .iter()
            .flatten()
            .flatten()
            .map(|v| [v.lo, v.hi])
            .collect::<Vec<_>>(),
        snapshot
    );

    // The same finite source coefficients in an extreme auxiliary parameter
    // chart cannot make an unresolved fast enclosure authoritative. Remember
    // its numerical refusal; callers must retain their original fallback.
    let narrow = Patch {
        index: patch.index,
        controls: patch.controls.clone(),
        normal: patch.normal.clone(),
        domain: [[0., 1e-300], [0., 1.]],
        normal_derivative: patch.normal_derivative,
        normal_hull: patch.normal_hull,
        curvature: patch.curvature,
        normal_jet: std::sync::OnceLock::new(),
    };
    assert_eq!(
        narrow.taylor_normal_box(narrow.domain, &mut w).err(),
        Some(Error::Unresolved)
    );
    assert!(narrow.normal_jet.get().unwrap().is_none());
    assert_eq!(
        narrow
            .taylor_normal_box(
                narrow.domain,
                &mut Work {
                    remaining: 0,
                    interrupt: None
                }
            )
            .err(),
        Some(Error::Unresolved)
    );
}

#[test]
fn normal_stencils_refuse_invalid_queries_and_damaged_derived_coefficients() {
    let body = super::super::tests::cam();
    let original = body.input().clone();
    let mut w = Work {
        remaining: 100_000_000,
        interrupt: None,
    };
    let mut prepared = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let tensor = prepared
        .sides
        .iter_mut()
        .find_map(|s| {
            if let Side::Tensor(t) = s {
                Some(t)
            } else {
                None
            }
        })
        .unwrap();
    for uv in [[0., 0.], [0.2, 0.7], [1., 1.]] {
        assert!(
            tensor
                .normal_hint(uv, &mut w)
                .unwrap()
                .unwrap()
                .iter()
                .all(|v| v.is_finite())
        );
    }
    for uv in [[f64::NAN, 0.], [0., f64::INFINITY], [-0.1, 0.5], [0.5, 1.1]] {
        assert_eq!(tensor.normal_hint(uv, &mut w).unwrap(), None);
    }
    assert_eq!(
        tensor
            .normal_hint(
                [0.5, 0.5],
                &mut Work {
                    remaining: 0,
                    interrupt: None
                }
            )
            .err(),
        Some(Error::Limit)
    );
    assert_eq!(
        tensor
            .normal_hint(
                [0.5, 0.5],
                &mut Work {
                    remaining: 100_000_000,
                    interrupt: Some(&|| true)
                }
            )
            .err(),
        Some(Error::Cancelled)
    );
    let original_normal = tensor.cells[0].normal.clone();
    tensor.cells[0].normal[0][0] = [I::point(f64::MAX); 3];
    tensor.cells[0].normal[17][17] = [I::point(-f64::MAX); 3];
    let uv = tensor.cells[0].domain.map(|d| d[0] + (d[1] - d[0]) * 0.5);
    assert_eq!(tensor.normal_hint(uv, &mut w).unwrap(), None);
    tensor.cells[0].normal = original_normal;
    assert!(tensor.normal_hint(uv, &mut w).unwrap().is_some());
    assert_eq!(body.input(), &original);
}

#[test]
fn uncertain_or_overflowed_stencils_fall_through_to_full_literal_source_checks() {
    let body = super::super::tests::cam();
    let original = body.input().clone();
    let mut w = Work {
        remaining: 100_000_000,
        interrupt: None,
    };
    let mut p = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let g = super::super::Grid {
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
        v: vec![0., 0.4, 1.],
    };
    let ring = Ring::new(&g);
    let mut candidate = p.emit(&g, &ring, 0.001, 50000, &mut w).unwrap();
    // Numerical fault injection affects only the nonauthoritative stencil;
    // the immutable physical boundary and candidate topology stay separate.
    for side in &mut p.sides {
        if let Side::Tensor(t) = side {
            for patch in &mut t.cells {
                patch.normal[0][0] = [I::point(f64::MAX); 3];
                patch.normal[17][17] = [I::point(-f64::MAX); 3];
            }
        }
    }
    for (side, source) in p.sides.iter().enumerate() {
        if matches!(source, Side::Tensor(_)) {
            assert_eq!(
                p.probe_normal(
                    Some(&candidate),
                    &g,
                    &ring,
                    [side, 0, 0],
                    0,
                    &mut BTreeMap::new(),
                    &mut w
                )
                .unwrap(),
                None
            );
        }
    }
    let mut changes = Refinement::new(g.u.len());
    assert!(
        p.normal_preflight(Some(&candidate), &g, &ring, 0, &mut changes, &mut w)
            .unwrap()
    );
    assert!(!changes.empty());
    for side in &mut p.sides {
        if let Side::Tensor(t) = side {
            for patch in &mut t.cells {
                for row in &mut patch.normal {
                    for n in row {
                        *n = [I::point(1e153), I::point(1.3e153), I::point(-1.8e153)];
                    }
                }
            }
        }
    }
    for x in &mut candidate.vertices {
        *x *= 1e74;
    }
    for (side, source) in p.sides.iter().enumerate() {
        if matches!(source, Side::Tensor(_)) {
            assert_eq!(
                p.probe_normal(
                    Some(&candidate),
                    &g,
                    &ring,
                    [side, 0, 0],
                    0,
                    &mut BTreeMap::new(),
                    &mut w
                )
                .unwrap(),
                None
            );
        }
    }
    let mut changes = Refinement::new(g.u.len());
    assert!(
        p.normal_preflight(Some(&candidate), &g, &ring, 0, &mut changes, &mut w)
            .unwrap()
    );
    assert!(!changes.empty());
    assert_eq!(
        p.check(&candidate, &g, &ring, 0, &mut changes, &mut w)
            .err(),
        Some(Error::Unresolved)
    );
    assert!(!changes.empty());
    assert_eq!(body.input(), &original);
}
