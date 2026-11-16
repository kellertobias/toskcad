use super::*;
use num_rational::BigRational as Q;
#[test]
fn actual_surface_failure_diagnostics_preserve_only_successful_resource_debits() {
    let source = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| (0..12).map(move |j| [i as f64, j as f64, 3.]))
            .collect(),
    )
    .unwrap();
    let invalid = source
        .surface_properties_with_diagnostics(
            SurfaceOptions11 {
                area_width_mm2: 0.,
                ..Default::default()
            },
            None,
        )
        .unwrap_err();
    assert!(!invalid.execution_started);
    assert_eq!(invalid.attempted_stage, None);
    assert_eq!(invalid.work_used, 0);
    assert_eq!(invalid.work_by_stage, [0; 10]);
    let cancelled = source
        .surface_properties_with_diagnostics(Default::default(), Some(&|| true))
        .unwrap_err();
    assert!(cancelled.execution_started);
    assert_eq!(cancelled.error, Error::Cancelled);
    assert_eq!(cancelled.attempted_stage, Some(SurfaceStage11::Setup));
    assert_eq!(cancelled.attempted_work, 1);
    assert_eq!(cancelled.work_used, 0);
    let work = source
        .surface_properties_with_diagnostics(
            SurfaceOptions11 {
                max_work: 1,
                ..Default::default()
            },
            None,
        )
        .unwrap_err();
    assert_eq!(work.error, Error::Limit);
    assert_eq!(work.work_used, 1);
    assert_eq!(work.work_by_stage.iter().sum::<usize>(), 1);
    assert!(work.attempted_work > work.remaining_work);
    let bytes = source
        .surface_properties_with_diagnostics(
            SurfaceOptions11 {
                max_bytes: 1,
                ..Default::default()
            },
            None,
        )
        .unwrap_err();
    assert_eq!(bytes.error, Error::Limit);
    assert_eq!(bytes.attempted_stage, Some(SurfaceStage11::Storage));
    assert!(bytes.attempted_bytes > 1);
    assert_eq!(bytes.accounted_peak_bytes, 0);
    assert_eq!(bytes.work_by_stage.iter().sum::<usize>(), bytes.work_used);
    let result = source
        .surface_properties_with_diagnostics(Default::default(), None)
        .unwrap();
    assert_eq!(result.work_by_stage.iter().sum::<usize>(), result.work_used);
    assert_eq!(
        result.work_by_stage[SurfaceStage11::NormalConstruction as usize],
        NORMAL11_WORK + SQUARED11_WORK
    );
    assert_eq!(
        result.work_by_stage[SurfaceStage11::Publication as usize],
        1
    );
    contains(result.area_mm2, q(121));
    assert_eq!(result.evaluated_cells, 1);
    let gaussian_prefix = result.work_by_stage[..=SurfaceStage11::Taylor as usize]
        .iter()
        .sum();
    let gaussian_failure = source
        .surface_properties_with_diagnostics(
            SurfaceOptions11 {
                max_work: gaussian_prefix,
                ..Default::default()
            },
            None,
        )
        .unwrap_err();
    assert_eq!(gaussian_failure.error, Error::Limit);
    assert!(matches!(
        gaussian_failure.attempted_stage,
        Some(SurfaceStage11::Gauss2 | SurfaceStage11::Gauss4)
    ));
    assert_eq!(gaussian_failure.work_used, gaussian_prefix);
    assert_eq!(gaussian_failure.remaining_work, 0);
    assert!(gaussian_failure.attempted_work > 0);
    assert_eq!(
        gaussian_failure.work_by_stage.iter().sum::<usize>(),
        gaussian_prefix
    );
    assert_eq!(
        source.surface_properties(Default::default(), None).unwrap(),
        result
    );
}
#[test]
fn actual_finite_surface_unattainable_precision_ties_priorities_and_refuses_without_publishing() {
    let source = curved();
    let original = source.clone();
    let smallest = f64::from_bits(1);
    let failure = source
        .surface_properties_with_diagnostics(
            SurfaceOptions11 {
                area_width_mm2: smallest,
                first_moment_width_mm3: smallest,
                centroid_width_mm: smallest,
                max_depth: 0,
                ..Default::default()
            },
            None,
        )
        .unwrap_err();
    assert_eq!(failure.error, Error::Unresolved);
    assert_eq!(failure.evaluated_cells, 1);
    assert_eq!(failure.axis_samples, [1, 1]);
    for r in failure.axis_remainder_max {
        assert!((r[0] / smallest).is_infinite());
    }
    assert_eq!(source, original);
    assert!(source.surface_properties(Default::default(), None).is_ok());
}
fn q(n: i64) -> Q {
    Q::from_integer(n.into())
}
fn contains(bounds: [f64; 2], expected: Q) {
    assert!(
        Q::from_float(bounds[0]).unwrap() <= expected
            && expected <= Q::from_float(bounds[1]).unwrap(),
        "{bounds:?} misses {expected}"
    );
}
fn curved() -> RetainedTensor11 {
    let poles = (0i64..12)
        .flat_map(|i| {
            (0i64..12).map(move |j| {
                [
                    (45 * i - i * (i - 1) * (i - 2) / 6) as f64 / 1024.,
                    45. * j as f64 / 1024.,
                    (9 * i * (i - 1) / 2) as f64 / 1024.,
                ]
            })
        })
        .collect();
    RetainedTensor11::new([1, 1], poles).unwrap()
}
fn choose(n: i64, k: i64) -> Q {
    (0..k).fold(q(1), |x, i| x * q(n - i) / q(i + 1))
}
fn basis(i: i64, t: &Q) -> Q {
    choose(11, i) * t.pow(i as i32) * (q(1) - t).pow((11 - i) as i32)
}
#[test]
fn actual_curved_d11_surface_integrals_enclose_independent_exact_polynomial() {
    let source = curved();
    for u in [Q::new(1.into(), 4.into()), Q::new(3.into(), 4.into())] {
        for v in [Q::new(1.into(), 2.into()), Q::new(7.into(), 8.into())] {
            let mut p = [q(0), q(0), q(0)];
            for i in 0..12 {
                for j in 0..12 {
                    let w = basis(i, &u) * basis(j, &v);
                    for (k, sum) in p.iter_mut().enumerate() {
                        *sum += &w
                            * Q::from_float(source.poles()[12 * i as usize + j as usize][k])
                                .unwrap();
                    }
                }
            }
            assert_eq!(
                p,
                [
                    (q(495) * &u - q(165) * u.pow(3)) / q(1024),
                    q(495) * &v / q(1024),
                    q(495) * u.pow(2) / q(1024)
                ]
            );
        }
    }
    let options = SurfaceOptions11 {
        area_width_mm2: 1e-5,
        first_moment_width_mm3: 1e-5,
        centroid_width_mm: 1e-5,
        ..Default::default()
    };
    let result = source.surface_properties(options, None).unwrap();
    eprintln!("D11_SURFACE_CURVED {result:?}");
    contains(result.area_mm2, q(326700) / q(1 << 20));
    let moments = [Q::new(148240125.into(), 2.into()), q(80858250), q(64686600)];
    for (bound, expected) in result.first_moment_mm3.into_iter().zip(moments) {
        contains(bound, expected / q(1 << 30));
    }
    for (bound, expected) in result.centroid_mm.into_iter().zip([
        Q::new(1815.into(), 8.into()),
        Q::new(495.into(), 2.into()),
        q(198),
    ]) {
        contains(bound, expected / q(1024));
    }
    assert!(result.work_used <= options.max_work);
    assert!(result.accounted_peak_bytes <= options.max_bytes);
    assert!(result.evaluated_cells > 1);
    for tiny in [
        SurfaceOptions11 {
            max_bytes: result.accounted_peak_bytes - 1,
            ..options
        },
        SurfaceOptions11 {
            max_work: result.work_used - 1,
            ..options
        },
    ] {
        assert_eq!(source.surface_properties(tiny, None), Err(Error::Limit));
    }
    let exact = source
        .surface_properties(
            SurfaceOptions11 {
                max_work: result.work_used,
                max_bytes: result.accounted_peak_bytes,
                ..options
            },
            None,
        )
        .unwrap();
    assert_eq!(exact.area_mm2, result.area_mm2);
    assert_eq!(exact.first_moment_mm3, result.first_moment_mm3);
    assert_eq!(exact.centroid_mm, result.centroid_mm);
    let before = source.clone();
    let calls = std::cell::Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    assert_eq!(
        source.surface_properties(options, Some(&count)).unwrap(),
        result
    );
    for stop in [calls.get() / 2, calls.get() - 1, calls.get()] {
        let count = std::cell::Cell::new(0);
        let cancel = || {
            count.set(count.get() + 1);
            count.get() == stop
        };
        assert_eq!(
            source.surface_properties(options, Some(&cancel)),
            Err(Error::Cancelled)
        );
        assert_eq!(count.get(), stop);
        assert_eq!(source, before);
    }
    assert_eq!(source.surface_properties(options, None).unwrap(), result);
}

#[test]
fn actual_two_root_d11_surface_refuses_initial_cell_limit_and_retries_exact_geometry() {
    let source = RetainedTensor11::new(
        [1, 2],
        (0..12)
            .flat_map(|i| (0..23).map(move |j| [i as f64, j as f64, 3.]))
            .collect(),
    )
    .unwrap();
    let before = source.clone();
    assert_eq!(
        source.surface_properties(
            SurfaceOptions11 {
                max_cells: 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    );
    let result = source
        .surface_properties(
            SurfaceOptions11 {
                max_cells: 2,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    contains(result.area_mm2, q(242));
    for (b, v) in result
        .centroid_mm
        .into_iter()
        .zip([Q::new(11.into(), 2.into()), q(11), q(3)])
    {
        contains(b, v);
    }
    for (b, v) in result
        .first_moment_mm3
        .into_iter()
        .zip([q(1331), q(2662), q(726)])
    {
        contains(b, v);
    }
    assert_eq!(result.evaluated_cells, 2);
    assert_eq!(source, before);
}
#[test]
fn actual_planar_d11_surface_resource_cancel_and_retry_are_atomic() {
    let source = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| (0..12).map(move |j| [i as f64, j as f64, 3.]))
            .collect(),
    )
    .unwrap();
    let original = source.clone();
    let options = SurfaceOptions11::default();
    let result = source.surface_properties(options, None).unwrap();
    contains(result.area_mm2, q(121));
    for (b, v) in result.centroid_mm.into_iter().zip([
        Q::new(11.into(), 2.into()),
        Q::new(11.into(), 2.into()),
        q(3),
    ]) {
        contains(b, v);
    }
    for (b, v) in result.first_moment_mm3.into_iter().zip([
        Q::new(1331.into(), 2.into()),
        Q::new(1331.into(), 2.into()),
        q(363),
    ]) {
        contains(b, v);
    }
    for changed in [
        SurfaceOptions11 {
            max_work: result.work_used - 1,
            ..options
        },
        SurfaceOptions11 {
            max_bytes: result.accounted_peak_bytes - 1,
            ..options
        },
    ] {
        assert_eq!(source.surface_properties(changed, None), Err(Error::Limit));
    }
    let exact = source
        .surface_properties(
            SurfaceOptions11 {
                max_work: result.work_used,
                max_bytes: result.accounted_peak_bytes,
                ..options
            },
            None,
        )
        .unwrap();
    assert_eq!(exact.area_mm2, result.area_mm2);
    let calls = std::cell::Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    assert_eq!(
        source.surface_properties(options, Some(&count)).unwrap(),
        result
    );
    for stop in [1, 2, calls.get() / 2, calls.get()] {
        let calls = std::cell::Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert_eq!(
            source.surface_properties(options, Some(&cancel)),
            Err(Error::Cancelled)
        );
        assert_eq!(calls.get(), stop);
        assert_eq!(source, original);
        assert_eq!(source.surface_properties(options, None).unwrap(), result);
    }
}

#[test]
fn actual_curved_d11_surface_proper_pose_preserves_raw_moments_and_requested_widths() {
    let local = curved();
    // Proper orthogonal rotation about Z followed by a literal translation.
    let source = RetainedTensor11::new(
        [1, 1],
        local
            .poles()
            .iter()
            .map(|p| [2. - p[1], -3. + p[0], 5. + p[2]])
            .collect(),
    )
    .unwrap();
    let options = SurfaceOptions11::default();
    let value = source.surface_properties(options, None).unwrap();
    let area = q(326700) / q(1 << 20);
    let local_m =
        [Q::new(148240125.into(), 2.into()), q(80858250), q(64686600)].map(|x| x / q(1 << 30));
    contains(value.area_mm2, area.clone());
    for (b, v) in value.first_moment_mm3.into_iter().zip([
        q(2) * &area - &local_m[1],
        -q(3) * &area + &local_m[0],
        q(5) * &area + &local_m[2],
    ]) {
        contains(b, v);
        assert!(I::new(b[0], b[1]).width() <= options.first_moment_width_mm3);
    }
    for (b, v) in value.centroid_mm.into_iter().zip([
        q(2) - Q::new(495.into(), 2048.into()),
        -q(3) + Q::new(1815.into(), 8192.into()),
        q(5) + Q::new(198.into(), 1024.into()),
    ]) {
        contains(b, v);
        assert!(I::new(b[0], b[1]).width() <= options.centroid_width_mm);
    }
    assert!(I::new(value.area_mm2[0], value.area_mm2[1]).width() <= options.area_width_mm2);
}

#[test]
fn actual_d11_surface_invalid_refinement_frontiers_and_degeneracy_refuse_atomically() {
    let source = curved();
    let original = source.clone();
    let options = SurfaceOptions11::default();
    for bad in [
        SurfaceOptions11 {
            area_width_mm2: f64::NAN,
            ..options
        },
        SurfaceOptions11 {
            area_width_mm2: 0.,
            ..options
        },
        SurfaceOptions11 {
            first_moment_width_mm3: f64::INFINITY,
            ..options
        },
        SurfaceOptions11 {
            centroid_width_mm: -1.,
            ..options
        },
        SurfaceOptions11 {
            max_work: 0,
            ..options
        },
        SurfaceOptions11 {
            max_work: MAX_WORK + 1,
            ..options
        },
        SurfaceOptions11 {
            max_bytes: 0,
            ..options
        },
        SurfaceOptions11 {
            max_bytes: MAX_BYTES + 1,
            ..options
        },
        SurfaceOptions11 {
            max_cells: 0,
            ..options
        },
        SurfaceOptions11 {
            max_cells: MAX_CELLS + 1,
            ..options
        },
        SurfaceOptions11 {
            max_depth: MAX_DEPTH + 1,
            ..options
        },
    ] {
        assert_eq!(source.surface_properties(bad, None), Err(Error::Invalid));
    }
    for tiny in [
        SurfaceOptions11 {
            max_work: 1,
            ..options
        },
        SurfaceOptions11 {
            max_bytes: 1,
            ..options
        },
        SurfaceOptions11 {
            max_cells: 1,
            ..options
        },
    ] {
        assert_eq!(source.surface_properties(tiny, None), Err(Error::Limit));
    }
    assert_eq!(
        source.surface_properties(
            SurfaceOptions11 {
                max_depth: 0,
                ..options
            },
            None
        ),
        Err(Error::Unresolved)
    );
    let degenerate = RetainedTensor11::new([1, 1], vec![[1., 2., 3.]; 144]).unwrap();
    assert_eq!(
        degenerate.surface_properties(
            SurfaceOptions11 {
                max_depth: 0,
                ..options
            },
            None
        ),
        Err(Error::Unresolved)
    );
    assert_eq!(source, original);
    assert!(source.surface_properties(options, None).is_ok());
}
