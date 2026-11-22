use super::*;
fn work() -> Work<'static> {
    Work {
        remaining: 100_000_000,
        interrupt: None,
    }
}
pub(in crate::radial_tensor::mesh) fn tensor<'a>(p: &'a Prepared<'_>) -> &'a Tensor {
    p.sides
        .iter()
        .find_map(|s| {
            if let Side::Tensor(t) = s {
                Some(t)
            } else {
                None
            }
        })
        .unwrap()
}
fn range(p: &tensor::Patch, v: [f64; 2]) -> [[f64; 2]; 2] {
    let fractions = [[0.11, 0.13], v];
    std::array::from_fn(|k| {
        fractions[k].map(|f| p.domain[k][0] + (p.domain[k][1] - p.domain[k][0]) * f)
    })
}
#[test]
fn actual_source_strips_reuse_u_without_changing_v_domains_or_literal_affines() {
    let body = super::super::tests::cam();
    let mut w = work();
    let prepared = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let owner = tensor(&prepared);
    let patches = owner.patches([[0., 1.]; 2], &mut w).unwrap();
    let patch = patches[0];
    let side = prepared
        .sides
        .iter()
        .position(|s| matches!(s, Side::Tensor(_)))
        .unwrap();
    let mut cache = Cache::new(owner, 2);
    let mut whole_v = range(patch, [0.11, 0.13]);
    whole_v[1] = [0., 1.];
    for range in [
        range(patch, [0.11, 0.13]),
        range(patch, [0.61, 0.63]),
        whole_v,
    ] {
        let uv = [
            [range[0][0], range[1][0]],
            [range[0][1], range[1][0]],
            [range[0][1], range[1][1]],
            [range[0][0], range[1][1]],
        ];
        let points = uv.map(|[u, v]| prepared.point(side, u, v, None, &mut w).unwrap());
        let cell = cache.cell(patch, range, &mut w).unwrap();
        let direct =
            certificate::SourceCell::new(&patch.controls, patch.domain, range, &mut w).unwrap();
        assert_eq!(cache.entries.len(), 1);
        for indices in [[0, 1, 2], [0, 2, 3]] {
            let affine =
                certificate::Affine::new(indices.map(|i| uv[i]), indices.map(|i| points[i]))
                    .unwrap();
            let original = direct.chord_error(&affine, &mut w).unwrap();
            assert_eq!(
                cell.chord_error(&affine, &mut w).unwrap().to_bits(),
                original.to_bits()
            );
            let displaced = certificate::Affine::new(
                indices.map(|i| uv[i]),
                indices.map(|i| {
                    let mut p = points[i];
                    p[0] += 1.;
                    p
                }),
            )
            .unwrap();
            assert!(cell.chord_error(&displaced, &mut w).unwrap() > 0.9);
        }
    }
    // The same parameter domain in another immutable source owner must not
    // obtain this owner's cached coefficient grid, even for equal geometry.
    let other = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let foreign = tensor(&other).patches([[0., 1.]; 2], &mut w).unwrap()[0];
    assert_eq!(
        cache
            .cell(foreign, range(foreign, [0.2, 0.3]), &mut w)
            .err(),
        Some(Error::Invalid)
    );
    assert_eq!(cache.entries.len(), 1);
    let posed = body
        .transformed_with_interrupt(
            [[30., -20., 10.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]],
            1e-7,
            1e-6,
            None,
        )
        .unwrap();
    let posed = Prepared::new(
        &posed,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let posed_patch = tensor(&posed).patches([[0., 1.]; 2], &mut w).unwrap()[0];
    assert_eq!(
        cache
            .cell(posed_patch, range(posed_patch, [0.2, 0.3]), &mut w)
            .err(),
        Some(Error::Invalid)
    );
    assert_eq!(cache.entries.len(), 1);
    cache.advance(patch.domain[1][1], &mut w).unwrap();
    assert_eq!(cache.entries.len(), 1);
    cache.advance(patch.domain[1][1].next_up(), &mut w).unwrap();
    assert!(cache.entries.is_empty());
}
#[test]
fn bounded_cache_caps_and_cold_hot_cancellation_leave_no_failed_entries() {
    let body = super::super::tests::cam();
    let mut w = work();
    let prepared = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let owner = tensor(&prepared);
    let patch = owner.patches([[0., 1.]; 2], &mut w).unwrap()[0];
    let range = range(patch, [0.2, 0.3]);
    for limit in [0, 1] {
        let mut cache = Cache::new(owner, limit);
        cache.cell(patch, range, &mut w).unwrap();
        assert_eq!(cache.entries.len(), limit);
        let mut other = range;
        other[0][0] += 0.001;
        cache.cell(patch, other, &mut w).unwrap();
        assert_eq!(cache.entries.len(), limit);
    }
    let count = std::cell::Cell::new(0);
    let observe = || {
        count.set(count.get() + 1);
        false
    };
    let mut cache = Cache::new(owner, 1);
    cache
        .cell(
            patch,
            range,
            &mut Work {
                remaining: 100_000_000,
                interrupt: Some(&observe),
            },
        )
        .unwrap();
    let cold = count.get();
    assert!(cold > 4);
    for stop in [1, 2, cold / 2, cold] {
        count.set(0);
        let cancel = || {
            count.set(count.get() + 1);
            count.get() >= stop
        };
        let mut cache = Cache::new(owner, 1);
        assert_eq!(
            cache
                .cell(
                    patch,
                    range,
                    &mut Work {
                        remaining: 100_000_000,
                        interrupt: Some(&cancel)
                    }
                )
                .err(),
            Some(Error::Cancelled)
        );
        assert!(cache.entries.is_empty());
    }
    for remaining in [0, 1, 200] {
        let mut cache = Cache::new(owner, 1);
        assert_eq!(
            cache
                .cell(
                    patch,
                    range,
                    &mut Work {
                        remaining,
                        interrupt: None
                    }
                )
                .err(),
            Some(Error::Limit)
        );
        assert!(cache.entries.is_empty());
    }
    let before = w.remaining;
    certificate::SourceCell::new(&patch.controls, patch.domain, range, &mut w).unwrap();
    let original_work = before - w.remaining;
    let before = w.remaining;
    cache.cell(patch, range, &mut w).unwrap();
    let reused_work = before - w.remaining;
    assert!(
        reused_work < original_work,
        "cached work {reused_work} must exclude actual repeated U arithmetic {original_work}"
    );
    // Hot retrieval and eviction both check work/cancellation before reuse or
    // mutation. A valid retained source strip survives an interrupted request.
    let cancel = || true;
    let stop: Option<&dyn Fn() -> bool> = Some(&cancel);
    for (remaining, interrupt, expected) in [
        (0, None, Error::Limit),
        (100_000_000, stop, Error::Cancelled),
    ] {
        assert_eq!(
            cache
                .cell(
                    patch,
                    range,
                    &mut Work {
                        remaining,
                        interrupt
                    }
                )
                .err(),
            Some(expected)
        );
        assert_eq!(cache.entries.len(), 1);
        assert_eq!(
            cache
                .advance(
                    patch.domain[1][1].next_up(),
                    &mut Work {
                        remaining,
                        interrupt
                    }
                )
                .err(),
            Some(expected)
        );
        assert_eq!(cache.entries.len(), 1);
    }
    cache.cell(patch, range, &mut w).unwrap();
    assert_eq!(cache.entries.len(), 1);
}

#[test]
fn actual_normal_strips_retain_each_coefficient_and_each_facet_orientation() {
    let body = super::super::tests::cam();
    let mut w = work();
    let p = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let owner = tensor(&p);
    let patch = owner.patches([[0., 1.]; 2], &mut w).unwrap()[0];
    let mut cache = Cache::new(owner, 1);
    for r in [range(patch, [0.11, 0.13]), range(patch, [0.61, 0.63])] {
        let coefficients = cache.normal_coefficients(patch, r, &mut w).unwrap();
        let direct = certificate::SourceStrip::new(&patch.normal, patch.domain, r[0], &mut w)
            .unwrap()
            .normal_coefficients(r[1], &mut w)
            .unwrap();
        for (a, b) in coefficients.iter().flatten().zip(direct.iter().flatten()) {
            for k in 0..3 {
                assert_eq!(
                    [a[k].lo.to_bits(), a[k].hi.to_bits()],
                    [b[k].lo.to_bits(), b[k].hi.to_bits()]
                );
            }
        }
        let n = patch
            .taylor_normal_box(r, &mut w)
            .unwrap()
            .map(|i| I::point(i.midpoint()));
        for facet in [n, n.map(|x| -x)] {
            assert_eq!(
                certificate::normal_coefficients_match(&coefficients, facet, 0.5, &mut w),
                certificate::normal_matches(
                    patch.normal.clone(),
                    patch.domain,
                    r,
                    facet,
                    0.5,
                    &mut w
                )
            );
        }
        assert!(
            !certificate::normal_coefficients_match(&coefficients, n.map(|x| -x), 0.5, &mut w)
                .unwrap()
        );
    }
    let other = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let foreign = tensor(&other).patches([[0., 1.]; 2], &mut w).unwrap()[0];
    assert_eq!(
        cache
            .normal_coefficients(foreign, range(foreign, [0.2, 0.3]), &mut w)
            .err(),
        Some(Error::Invalid)
    );
    let r = range(patch, [0.2, 0.3]);
    let mut uncached = Cache::new(owner, 0);
    assert!(uncached.normal_coefficients(patch, r, &mut w).is_ok());
    assert!(uncached.normals.is_empty());
    for budget in [0, 100] {
        let mut limited = Work {
            remaining: budget,
            interrupt: None,
        };
        assert_eq!(
            uncached.normal_coefficients(patch, r, &mut limited).err(),
            Some(Error::Limit)
        );
        assert!(uncached.normals.is_empty());
    }
    let mut cancelled = Work {
        remaining: 100_000_000,
        interrupt: Some(&|| true),
    };
    assert_eq!(
        uncached.normal_coefficients(patch, r, &mut cancelled).err(),
        Some(Error::Cancelled)
    );
    assert!(uncached.normals.is_empty());
    assert_eq!(cache.normals.len(), 1);
    cache.advance(patch.domain[1][1].next_up(), &mut w).unwrap();
    assert!(cache.normals.is_empty());
    assert!(cache.normal_coefficients(patch, r, &mut w).is_ok());
}

#[test]
fn continued_position_strips_reset_on_gaps_and_nonmonotone_rows() {
    let body = super::super::tests::cam();
    let mut w = work();
    let p = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let owner = tensor(&p);
    let patch = owner.patches([[0., 1.]; 2], &mut w).unwrap()[0];
    let mut cache = Cache::new(owner, 2).with_continuation(true);
    for v in [[0., 0.125], [0.125, 0.25], [0.5, 0.75], [0., 0.125]] {
        let r = range(patch, v);
        let cell = cache.cell(patch, r, &mut w).unwrap();
        let side = p
            .sides
            .iter()
            .position(|s| matches!(s, Side::Tensor(_)))
            .unwrap();
        let uv = [[r[0][0], r[1][0]], [r[0][1], r[1][0]], [r[0][1], r[1][1]]];
        let vertices = uv.map(|[u, v]| p.point(side, u, v, None, &mut w).unwrap());
        let affine = certificate::Affine::new(uv, vertices).unwrap();
        assert!(cell.chord_error(&affine, &mut w).unwrap().is_finite());
        assert_eq!(cache.entries.len(), 1);
    }
}
