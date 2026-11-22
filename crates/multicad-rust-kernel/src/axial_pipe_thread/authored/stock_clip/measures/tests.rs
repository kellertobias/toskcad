use super::*;
#[test]
fn actual_canonical_faces_have_bounded_source_area_moments() {
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    for face in [
        StockClipFace::LowerFlank,
        StockClipFace::NaturalCrest,
        StockClipFace::StockCylinder,
        StockClipFace::UpperFlank,
        StockClipFace::Root,
        StockClipFace::StartCap,
        StockClipFace::EndCap,
    ] {
        let start = std::time::Instant::now();
        let value = atlas.face_measure(face, Default::default(), None);
        eprintln!(
            "NPT source measure {face:?}: {value:?} {:?}",
            start.elapsed()
        );
        let value = value.unwrap();
        assert!(value.area_mm2()[0] > 0.);
    }
}
#[test]
fn gauss_actual_isolated_nodes_integrate_exact_monomials_and_remainder_constant() {
    let nodes = gauss_nodes();
    let mut positive = nodes[..3]
        .iter()
        .map(|v| v.0.midpoint())
        .collect::<Vec<_>>();
    positive.sort_by(f64::total_cmp);
    assert!(positive[0] > 0. && positive[2] < 1. && positive.windows(2).all(|v| v[0] < v[1]));
    for degree in 0..=12 {
        let q = nodes
            .iter()
            .fold(I::point(0.), |sum, (x, w)| sum + power(*x, degree) * *w);
        let exact = if degree % 2 == 0 {
            I::point(2.) / I::point((degree + 1) as f64)
        } else {
            I::point(0.)
        };
        if degree <= 11 {
            assert!(
                q.lo <= exact.lo && q.hi >= exact.hi,
                "degree{degree}: {q:?} {exact:?}"
            );
        } else {
            let remainder = power(I::point(2.), 13) / I::point(11099088.);
            let actual = exact - q;
            assert!(actual.lo <= remainder.lo && remainder.hi <= actual.hi);
        }
    }
}
#[test]
fn independent_conical_cylindrical_and_cap_area_formulas_match_source_integrals() {
    for quarter in [false, true] {
        for scale in [0.1, 1., 3.] {
            let source = super::super::tests::canonical(quarter, scale);
            let spec = source.spec();
            let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
            let AuthoredSection::NptFlat {
                depth_mm: d,
                flat_width_mm: w,
            } = spec.section
            else {
                unreachable!()
            };
            let r = spec.major_radius_mm;
            let k = spec.radial_taper;
            let l = spec.length_mm;
            let stock = spec.crest_radius_limit_mm.unwrap();
            let zc = (stock - r) / k;
            let delta = l - zc;
            let m = 3_f64.sqrt();
            let slopes = 1. / (m - k) + 1. / (m + k);
            let root_width = spec.pitch_mm - w - d * slopes;
            let omega = std::f64::consts::TAU / spec.pitch_mm;
            let endroot = r - d + k * l;
            let c = stock - endroot;
            let endwidth = w + (d - c) * slopes;
            let examples = [
                (
                    StockClipFace::Root,
                    omega * root_width * l * (r - d + 0.5 * k * l) * (1. + k * k).sqrt(),
                ),
                (
                    StockClipFace::NaturalCrest,
                    omega * w * zc * (r + 0.5 * k * zc) * (1. + k * k).sqrt(),
                ),
                (
                    StockClipFace::StockCylinder,
                    omega * stock * (w * delta + 0.5 * k * delta * delta * slopes),
                ),
                (
                    StockClipFace::StartCap,
                    omega
                        * 0.5
                        * (root_width * (r - d).powi(2)
                            + w * r * r
                            + d * slopes * ((r - d).powi(2) + (r - d) * d + d * d / 3.)),
                ),
                (
                    StockClipFace::EndCap,
                    omega
                        * 0.5
                        * (root_width * endroot * endroot
                            + endwidth * stock * stock
                            + c * slopes * (endroot * endroot + endroot * c + c * c / 3.)),
                ),
            ];
            for face in [StockClipFace::LowerFlank, StockClipFace::UpperFlank] {
                let value = atlas.face_measure(face, Default::default(), None);
                assert!(value.is_ok(), "{quarter}/{scale}/{face:?}: {value:?}");
            }
            for (face, expected) in examples {
                let value = atlas.face_measure(face, Default::default(), None).unwrap();
                let actual = value.area_mm2();
                assert!(
                    actual[0] <= expected && expected <= actual[1],
                    "{quarter}/{scale}/{face:?}: {actual:?} {expected}"
                );
                assert_eq!(value.face(), face);
                assert!(
                    value
                        .first_moments_mm3()
                        .iter()
                        .all(|v| v[0].is_finite() && v[1].is_finite())
                );
                assert!(value.centroid_mm().iter().all(|v| v[1] - v[0] <= 1e-7));
            }
        }
    }
}
#[test]
fn bounded_source_measure_failures_preserve_source_and_successful_retry() {
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let face = StockClipFace::StockCylinder;
    let options = StockClipMeasureOptions::default();
    let good = atlas.face_measure(face, options, None).unwrap();
    for v in [0., -1., f64::NAN, f64::INFINITY] {
        for opts in [
            StockClipMeasureOptions {
                area_width_mm2: v,
                ..options
            },
            StockClipMeasureOptions {
                moment_width_mm3: v,
                ..options
            },
            StockClipMeasureOptions {
                centroid_width_mm: v,
                ..options
            },
        ] {
            assert!(matches!(
                atlas.face_measure(face, opts, None),
                Err(Error::Invalid)
            ));
        }
    }
    for limits in [
        StockClipOptions {
            maximum_work: good.work_used() - 1,
            ..Default::default()
        },
        StockClipOptions {
            maximum_bytes: good.retained_budget_bytes() - 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            atlas.face_measure(face, StockClipMeasureOptions { limits, ..options }, None),
            Err(Error::Limit)
        ));
    }
    let limits = StockClipOptions {
        maximum_work: good.work_used(),
        maximum_bytes: good.retained_budget_bytes(),
    };
    assert_eq!(
        atlas
            .face_measure(face, StockClipMeasureOptions { limits, ..options }, None)
            .unwrap()
            .area_mm2(),
        good.area_mm2()
    );
    assert!(matches!(
        atlas.face_measure(face, options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let polls = std::cell::Cell::new(0usize);
    let observe = || {
        polls.set(polls.get() + 1);
        false
    };
    atlas.face_measure(face, options, Some(&observe)).unwrap();
    let last = polls.get() - 1;
    let count = std::cell::Cell::new(0usize);
    let cancel = || {
        let n = count.get();
        count.set(n + 1);
        n >= last
    };
    assert!(matches!(
        atlas.face_measure(face, options, Some(&cancel)),
        Err(Error::Cancelled)
    ));
    assert!(
        atlas
            .face_measure(
                face,
                StockClipMeasureOptions {
                    centroid_width_mm: 1e-30,
                    ..options
                },
                None
            )
            .is_err()
    );
    assert_eq!(*atlas.source(), source);
    assert_eq!(
        atlas.face_measure(face, options, None).unwrap().area_mm2(),
        good.area_mm2()
    );
}
fn tangent_cross_oracle(atlas: &StockClipPlan, face: StockClipFace, n: usize) -> [f64; 4] {
    let mut sum = [0.; 4];
    for patch in atlas.patches().into_iter().filter(|p| p.face() == face) {
        let mut corners = [[0.; 3]; 4];
        for (i, (u, v)) in [(0., 0.), (1., 0.), (0., 1.), (1., 1.)]
            .into_iter()
            .enumerate()
        {
            let p = atlas
                .point_bounds(patch, [u, u], [v, v], Default::default(), None)
                .unwrap();
            corners[i] = [p.radius_mm(), p.phase_mm(), p.axial_mm()].map(|v| (v[0] + v[1]) * 0.5);
        }
        let omega = std::f64::consts::TAU / atlas.source().spec().pitch_mm;
        for i in 0..=n {
            for j in 0..=n {
                let u = i as f64 / n as f64;
                let v = j as f64 / n as f64;
                let f: [f64; 3] = std::array::from_fn(|a| {
                    (1. - v) * ((1. - u) * corners[0][a] + u * corners[1][a])
                        + v * ((1. - u) * corners[2][a] + u * corners[3][a])
                });
                let du: [f64; 3] = std::array::from_fn(|a| {
                    (1. - v) * (corners[1][a] - corners[0][a]) + v * (corners[3][a] - corners[2][a])
                });
                let dv: [f64; 3] = std::array::from_fn(|a| {
                    (1. - u) * (corners[2][a] - corners[0][a]) + u * (corners[3][a] - corners[1][a])
                });
                let theta = omega * (f[2] - f[1]);
                let (sin, cos) = theta.portable_sin_cos();
                let tangent = |d: [f64; 3]| {
                    [
                        d[0] * cos - f[0] * sin * omega * (d[2] - d[1]),
                        d[0] * sin + f[0] * cos * omega * (d[2] - d[1]),
                        d[2],
                    ]
                };
                let a = tangent(du);
                let b = tangent(dv);
                let cross = [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ];
                let area = cross[0].portable_hypot(cross[1]).portable_hypot(cross[2]);
                let weight = |k| {
                    if k == 0 || k == n {
                        1.
                    } else if k % 2 == 0 {
                        2.
                    } else {
                        4.
                    }
                };
                let w = weight(i) * weight(j) / (9. * (n * n) as f64);
                let p = [1., f[0] * cos, f[0] * sin, f[2]];
                for a in 0..4 {
                    sum[a] += w * area * p[a];
                }
            }
        }
    }
    sum
}
#[test]
fn independent_tangent_cross_quadrature_checks_both_flanks_and_first_moments() {
    for quarter in [false, true] {
        let source = super::super::tests::canonical(quarter, 1.);
        let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
        for face in [StockClipFace::LowerFlank, StockClipFace::UpperFlank] {
            let actual = atlas.face_measure(face, Default::default(), None).unwrap();
            let values = [
                actual.area_mm2(),
                actual.first_moments_mm3()[0],
                actual.first_moments_mm3()[1],
                actual.first_moments_mm3()[2],
            ];
            let coarse = tangent_cross_oracle(&atlas, face, 256);
            let fine = tangent_cross_oracle(&atlas, face, 512);
            for j in 0..4 {
                // This independent floating Simpson oracle is a convergence check,
                // not the authority for the published interval certificate.
                let oracle_error = 2. * (coarse[j] - fine[j]).abs() + 1e-8;
                assert!(
                    fine[j] >= values[j][0] - oracle_error
                        && fine[j] <= values[j][1] + oracle_error,
                    "{face:?}/{j} {fine:?} {values:?}"
                );
                assert!(oracle_error < 2e-4);
            }
        }
    }
}

#[test]
fn immutable_gauss_constants_have_exact_root_and_coupled_weight_authority() {
    use num_rational::BigRational as Q;
    use num_traits::Zero;
    let q = |n: i64| Q::from_integer(n.into());
    let exact_p = |x: &Q| {
        let x2 = x * x;
        (((q(231) * &x2 - q(315)) * &x2 + q(105)) * &x2 - q(5)) / q(16)
    };
    let exact_d = |x: &Q| {
        let x2 = x * x;
        ((q(1386) * &x2 - q(1260)) * &x2 + q(210)) * x / q(16)
    };
    let nodes = gauss_nodes();
    for (i, (x, w)) in nodes[..3].iter().copied().enumerate() {
        assert!(0. < x.lo && x.lo < x.hi && x.hi < 1.);
        if i > 0 {
            assert!(nodes[i - 1].0.hi < x.lo);
        }
        let lo = Q::from_float(x.lo).unwrap();
        let hi = Q::from_float(x.hi).unwrap();
        assert!(exact_p(&lo) * exact_p(&hi) < Q::zero());
        // Positive derivative of D=(1-x²)P6'² proves the endpoint weight
        // ordering on the WHOLE root box, not merely at a fitted root.
        let dp = ((I::point(1386.) * x.square() - I::point(1260.)) * x.square() + I::point(210.))
            * x
            / I::point(16.);
        let dd = ((I::point(6930.) * x.square() - I::point(3780.)) * x.square() + I::point(210.))
            / I::point(16.);
        let derivative =
            -I::point(2.) * x * dp.square() + I::point(2.) * (I::point(1.) - x.square()) * dp * dd;
        assert!(derivative.lo > 0.);
        let exact_weight = |x: &Q| {
            let d = exact_d(x);
            q(2) / ((q(1) - x * x) * &d * &d)
        };
        assert!(Q::from_float(w.lo).unwrap() <= exact_weight(&hi));
        assert!(exact_weight(&lo) <= Q::from_float(w.hi).unwrap());
        assert_eq!(nodes[i + 3].0.lo, -x.hi);
        assert_eq!(nodes[i + 3].0.hi, -x.lo);
        assert_eq!(nodes[i + 3].1.lo, w.lo);
        assert_eq!(nodes[i + 3].1.hi, w.hi);
    }
}
#[test]
fn real_partial_point_work_and_unrepresentable_area_refuse_before_receipt() {
    let source = super::super::tests::canonical(false, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    let point = atlas
        .point_bounds(
            atlas.patches()[0],
            [0., 0.],
            [0., 0.],
            Default::default(),
            None,
        )
        .unwrap();
    let limits = StockClipOptions {
        maximum_work: atlas.work_used() + point.work_used() - 1,
        ..Default::default()
    };
    assert!(matches!(
        atlas.face_measure(
            StockClipFace::LowerFlank,
            StockClipMeasureOptions {
                limits,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        atlas.face_measure(
            StockClipFace::LowerFlank,
            StockClipMeasureOptions {
                area_width_mm2: 1e-30,
                ..Default::default()
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    let mut spec = source.spec();
    let AuthoredSection::NptFlat { depth_mm, .. } = spec.section else {
        unreachable!()
    };
    spec.section = AuthoredSection::NptFlat {
        depth_mm,
        flat_width_mm: f64::from_bits(1),
    };
    let tiny = AuthoredThread::new(spec, 1e-7, 1e-5).unwrap();
    let tiny = tiny.stock_clip_plan(Default::default(), None).unwrap();
    assert!(matches!(
        tiny.face_measure(StockClipFace::NaturalCrest, Default::default(), None),
        Err(Error::Unresolved)
    ));
}
#[test]
fn actual_deep_source_retains_positive_convex_radius_in_derivative_domains() {
    let source = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 1.,
            pitch_mm: 2.,
            length_mm: 2.,
            radial_taper: 1. / 32.,
            section: AuthoredSection::NptFlat {
                depth_mm: 0.99,
                flat_width_mm: 0.02,
            },
            crest_radius_limit_mm: Some(1.02),
        },
        1e-7,
        1e-5,
    )
    .unwrap();
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    // The source is valid, but the default moment width on this small
    // stock face does not establish the tighter default centroid width.
    // Preserve that honest public precision refusal; no Native admission is
    // inferred for this custom source from a looser successful query.
    assert!(matches!(
        atlas.face_measure(StockClipFace::StockCylinder, Default::default(), None),
        Err(Error::Unresolved)
    ));
    let options = StockClipMeasureOptions {
        centroid_width_mm: 1e-5,
        ..Default::default()
    };
    for face in atlas.physical_faces() {
        let result = atlas.face_measure(face, options, None);
        assert!(result.is_ok(), "deep source {face:?}: {result:?}");
        let value = result.unwrap();
        assert!(value.area_mm2()[0] > 0.);
        assert!(
            value
                .centroid_mm()
                .iter()
                .all(|v| v[1] - v[0] <= options.centroid_width_mm)
        );
    }
    assert_eq!(*atlas.source(), source);
}
