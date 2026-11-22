//! Independent integer binary-rational oracles, not admitted body fixtures.
use super::*;
fn parts(x: f64) -> (i128, i32) {
    let bits = x.to_bits();
    let exponent = ((bits >> 52) & 2047) as i32;
    let mantissa = (bits & ((1u64 << 52) - 1)) | if exponent == 0 { 0 } else { 1u64 << 52 };
    let sign = if bits >> 63 == 0 { 1 } else { -1 };
    (
        sign * mantissa as i128,
        if exponent == 0 {
            -1074
        } else {
            exponent - 1023 - 52
        },
    )
}
fn integer(x: f64, base: i32) -> i128 {
    let (m, e) = parts(x);
    if m == 0 {
        0
    } else {
        assert!((0..74).contains(&(e - base)));
        m << (e - base)
    }
}
#[test]
fn exact_translation_residuals_match_independent_binary_integer_sums() {
    for a in [
        0.,
        -0.,
        0.1,
        -0.1,
        4.973885390363372,
        -4.398177689636629,
        1.,
        1_f64.next_up(),
        1e6,
    ] {
        for b in [0., 11., -7., 5., 1e6, -1e6, 2_f64.powi(-53)] {
            let (s, r) = two_sum(a, b).unwrap();
            let base = [a, b, s, r]
                .into_iter()
                .map(parts)
                .filter(|(m, _)| *m != 0)
                .map(|(_, e)| e)
                .min()
                .unwrap_or(0);
            // Bounded fixture exponents fit one exact128-bit integer. This
            // oracle does not use the floating TwoSum recurrence.
            assert_eq!(
                integer(a, base) + integer(b, base),
                integer(s, base) + integer(r, base)
            );
        }
    }
    let tiny = f64::from_bits(1);
    assert_eq!(two_sum(1., tiny), Some((1., tiny)));
    assert_eq!(two_sum(-1., tiny), Some((-1., tiny)));
    assert!(two_sum(f64::MAX, f64::MAX).is_none());
    assert!(two_sum(f64::INFINITY, f64::NEG_INFINITY).is_none());
}
#[test]
fn exact_permutation_chain_keeps_coupled_neck_and_bounds_rounding_without_normalization() {
    let frame = PlaneFrame {
        origin: [11., -7., 5.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    };
    let length = 7.865062962962963;
    let p = Placement::new(
        frame,
        PipeHardwareKind::HexHead {
            across_flats_mm: 18.,
            height_mm: 3.,
        },
        length,
        GeometryPolicy::default(),
    )
    .unwrap();
    assert_eq!(p.point_work(), 128);
    assert_eq!(p.sigma, [1., 1.]);
    assert_eq!(p.determinant, I::ONE);
    let (neck, error) = p.point([0., 0., length]).unwrap();
    assert_eq!(neck, frame.origin);
    assert_eq!(error, 0.);
    for point in [
        [4.973885390363372, -0.1, 0.],
        [1., 2., 3.],
        [-4.398177689636629, 0.9407407407407407, 6.5],
    ] {
        let (mapped, error) = p.point(point).unwrap();
        let exact = p.point_bounds(point.map(I::point)).unwrap();
        let mut inverse = p.inverse_point(mapped.map(I::point)).unwrap();
        for k in 0..3 {
            assert!(exact[k].lo <= mapped[k] + error && mapped[k] - error <= exact[k].hi);
            inverse[k] = inverse[k] + I::new(-error, error);
            assert!(inverse[k].lo <= point[k] && point[k] <= inverse[k].hi);
        }
    }
    let ordinary = Placement::new(
        PlaneFrame {
            origin: [10., 20., 30.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        },
        PipeHardwareKind::HexNut {
            across_flats_mm: 18.,
        },
        length,
        GeometryPolicy::default(),
    )
    .unwrap();
    assert_eq!(ordinary.point_work(), 128);
    assert_eq!(ordinary.normal_work(), 160);
    assert!(ordinary.point([f64::NAN, 0., 0.]).is_err());
    assert!(Permutation::new([[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., -1.]]).is_none());
    assert!(Permutation::new([[0.; 3], [1., 0., 0.], [1., 0., 0.], [0., 0., 1.]]).is_none());
}
#[test]
fn actual_generic_affine_fma_residual_bound_contains_exact_binary_rational_dot_products() {
    let columns = [
        [10., 20., 30.],
        [0.8, 0.6, 0.],
        [-0.36, 0.48, 0.8],
        [0.48, -0.64, 0.6],
    ];
    let affine = Affine::new(columns).unwrap();
    for point in [
        [1., 2., 3.],
        [0.1, -4.398177689636629, 5.091262310363371],
        [1e6, -1e6, 1e6],
    ] {
        let (actual, error) = affine.rounded_point(point).unwrap();
        let terms = (0..3)
            .flat_map(|k| {
                std::iter::once((columns[0][k], 1.))
                    .chain((0..3).map(move |j| (columns[j + 1][k], point[j])))
            })
            .collect::<Vec<_>>();
        let base = terms
            .iter()
            .filter_map(|&(a, b)| {
                let (ma, ea) = parts(a);
                let (mb, eb) = parts(b);
                (ma != 0 && mb != 0).then_some(ea + eb)
            })
            .min()
            .unwrap();
        let mut exact_error = 0i128;
        for k in 0..3 {
            let exact = terms[4 * k..4 * k + 4]
                .iter()
                .map(|&(a, b)| {
                    let (ma, ea) = parts(a);
                    let (mb, eb) = parts(b);
                    if ma == 0 || mb == 0 {
                        0
                    } else {
                        assert!((0..21).contains(&(ea + eb - base)));
                        (ma * mb) << (ea + eb - base)
                    }
                })
                .sum::<i128>();
            exact_error += (exact - integer(actual[k], base)).abs();
        }
        assert!(exact_error <= integer(error, base));
    }
}
#[test]
fn actual_raw_facet_normal_fast_bound_and_charged_parallel_fallback_keep_refusals() {
    let pose = Placement::new(
        PlaneFrame::default(),
        PipeHardwareKind::HexNut {
            across_flats_mm: 18.,
        },
        6.,
        GeometryPolicy::default(),
    )
    .unwrap();
    let local = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]];
    let charged = std::cell::Cell::new(0usize);
    pose.require_normal(local, local, 1e-6, |n| {
        charged.set(charged.get() + n);
        Ok(())
    })
    .unwrap();
    assert_eq!(charged.get(), 0);
    // A raw triangle with twice the area has the SAME normal direction. The
    // first sufficient bound is conservative; the full original angular
    // proof must still accept it after explicitly charging its own work.
    let parallel = [[0., 0., 0.], [2., 0., 0.], [0., 1., 0.]];
    pose.require_normal(local, parallel, 1e-6, |n| {
        charged.set(charged.get() + n);
        Ok(())
    })
    .unwrap();
    assert_eq!(charged.get(), 64);
    assert!(matches!(
        pose.require_normal(local, parallel, 1e-6, |_| Err(failure(Cause::Limit))),
        Err(GeometryError {
            kind: GeometryErrorKind::Limit,
            ..
        })
    ));
    assert!(
        pose.require_normal(local, [local[0], local[2], local[1]], 1e-6, |_| Ok(()))
            .is_err()
    );
    assert!(
        pose.require_normal(local, [[0.; 3]; 3], 1e-6, |_| Ok(()))
            .is_err()
    );
}

#[test]
fn cached_original_one_or_two_factor_points_enclose_independent_binary_rational_world_points() {
    fn reduced(x: f64) -> (i128, i32) {
        let (m, e) = parts(x);
        if m == 0 {
            (0, 0)
        } else {
            let shift = m.unsigned_abs().trailing_zeros();
            (m >> shift, e + shift as i32)
        }
    }
    fn scaled(x: f64, base: i32) -> i128 {
        let (m, e) = reduced(x);
        if m == 0 {
            return 0;
        }
        assert!(e >= base, "oracle must retain every exact dyadic bit");
        m.checked_mul(1i128.checked_shl((e - base) as u32).unwrap())
            .unwrap()
    }
    for origin in [
        [10., 20., 30.],
        [1e6, -2e6, 3e6],
        [-838860.8, -629145.6, 0.],
    ] {
        for head in [false, true] {
            let columns = [
                origin,
                [0.8, 0.6, 0.],
                [-0.36, 0.48, 0.8],
                [0.48, -0.64, 0.6],
            ];
            let mut pose = Placement::new(
                PlaneFrame::default(),
                if head {
                    PipeHardwareKind::HexHead {
                        across_flats_mm: 16.,
                        height_mm: 3.,
                    }
                } else {
                    PipeHardwareKind::HexNut {
                        across_flats_mm: 16.,
                    }
                },
                6.5,
                GeometryPolicy::default(),
            )
            .unwrap();
            pose.append(columns).unwrap();
            assert_eq!(pose.point_work(), 128);
            assert!(pose.point([f64::NAN, 0., 0.]).is_err());
            assert!(pose.point([f64::MAX, f64::MAX, f64::MAX]).is_err());
            for point in [
                [1., 2., 3.],
                [-2., 0.5, 6.5],
                [0., 0., 6.5],
                [1048576., 0., 6.5],
            ] {
                let (actual, error) = pose.point(point).unwrap();
                if !head {
                    let legacy = pose.factors[0].rounded_point(point).unwrap().0;
                    eprintln!(
                        "singleOriginalPoint origin={origin:?} local={point:?} literal_changed={} legacy={legacy:?} directed={actual:?}",
                        legacy != actual
                    );
                }
                // Independent integer arithmetic evaluates the ORIGINAL single
                // affine expression, optionally after the exact dyadic head map.
                // No cached interval field or normalized matrix is an oracle.
                let local = if head {
                    [point[0], -point[1], 6.5 - point[2]]
                } else {
                    point
                };
                let terms = (0..3)
                    .flat_map(|k| {
                        std::iter::once((columns[0][k], 1.))
                            .chain((0..3).map(move |j| (columns[j + 1][k], local[j])))
                    })
                    .collect::<Vec<_>>();
                let base = terms
                    .iter()
                    .filter_map(|&(a, b)| {
                        let (ma, ea) = reduced(a);
                        let (mb, eb) = reduced(b);
                        (ma != 0 && mb != 0).then_some(ea + eb)
                    })
                    .chain(
                        actual
                            .into_iter()
                            .chain(std::iter::once(error))
                            .filter_map(|x| {
                                let (m, e) = reduced(x);
                                (m != 0).then_some(e)
                            }),
                    )
                    .min()
                    .unwrap();
                let mut exact_error = 0i128;
                for k in 0..3 {
                    let exact = terms[4 * k..4 * k + 4]
                        .iter()
                        .map(|&(a, b)| {
                            let (ma, ea) = reduced(a);
                            let (mb, eb) = reduced(b);
                            if ma == 0 || mb == 0 {
                                0
                            } else {
                                (ma * mb) << (ea + eb - base)
                            }
                        })
                        .sum::<i128>();
                    exact_error += (exact - scaled(actual[k], base)).abs();
                }
                assert!(exact_error <= scaled(error, base));
            }
        }
    }
}

#[test]
fn cached_three_original_factors_enclose_independent_expanded_world_expression() {
    fn dyadic(x: f64) -> (i128, i32) {
        let (m, e) = parts(x);
        if m == 0 {
            (0, 0)
        } else {
            let shift = m.unsigned_abs().trailing_zeros();
            (m >> shift, e + shift as i32)
        }
    }
    fn scaled(m: i128, e: i32, base: i32) -> i128 {
        if m == 0 {
            return 0;
        }
        assert!(e >= base);
        m.checked_mul(1i128.checked_shl((e - base) as u32).unwrap())
            .unwrap()
    }
    let first = [
        [0.25, -0.5, 0.75],
        [0.8, 0.6, 0.],
        [-0.36, 0.48, 0.8],
        [0.48, -0.64, 0.6],
    ];
    let second = [
        [-0.5, 0.25, 0.125],
        [0.8, 0.6, 0.],
        [-0.6, 0.8, 0.],
        [0., 0., 1.],
    ];
    let mut pose = Placement::new(
        PlaneFrame::default(),
        PipeHardwareKind::HexHead {
            across_flats_mm: 16.,
            height_mm: 3.,
        },
        6.5,
        GeometryPolicy::default(),
    )
    .unwrap();
    pose.append(first).unwrap();
    pose.append(second).unwrap();
    for point in [[1., 2., 3.], [-2., 0.5, 6.5], [0., 0., 6.5]] {
        let (actual, error) = pose.point(point).unwrap();
        let local = [point[0], -point[1], 6.5 - point[2]];
        // Expand the two ORIGINAL arbitrary factors over the exact dyadic
        // head factor. This uses integer polynomial arithmetic, neither the
        // interval composition nor any cached placement field.
        let terms: [Vec<(i128, i32)>; 3] = std::array::from_fn(|k| {
            let mut terms = vec![vec![second[0][k]]];
            for j in 0..3 {
                terms.push(vec![second[j + 1][k], first[0][j]]);
                for l in 0..3 {
                    terms.push(vec![second[j + 1][k], first[l + 1][j], local[l]]);
                }
            }
            terms
                .into_iter()
                .map(|term| {
                    term.into_iter()
                        .map(dyadic)
                        .fold((1i128, 0), |(m, e), (n, f)| {
                            (m.checked_mul(n).unwrap(), e + f)
                        })
                })
                .collect()
        });
        let base = terms
            .iter()
            .flatten()
            .filter_map(|&(m, e)| (m != 0).then_some(e))
            .chain(
                actual
                    .into_iter()
                    .chain(std::iter::once(error))
                    .filter_map(|v| {
                        let (m, e) = dyadic(v);
                        (m != 0).then_some(e)
                    }),
            )
            .min()
            .unwrap();
        let mut exact_error = 0i128;
        for k in 0..3 {
            let exact = terms[k]
                .iter()
                .map(|&(m, e)| scaled(m, e, base))
                .sum::<i128>();
            let (m, e) = dyadic(actual[k]);
            exact_error = exact_error
                .checked_add((exact - scaled(m, e, base)).abs())
                .unwrap();
        }
        let (m, e) = dyadic(error);
        assert!(exact_error <= scaled(m, e, base));
    }
}
