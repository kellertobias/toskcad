use super::*;
fn source(rounded: bool, external: bool, limit: Option<f64>) -> AuthoredThread {
    AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 5.,
            pitch_mm: 1.,
            length_mm: 6.,
            radial_taper: if external { 1. / 32. } else { -1. / 32. },
            section: if rounded {
                AuthoredSection::BsptRounded { depth_mm: 0.6 }
            } else {
                AuthoredSection::NptFlat {
                    depth_mm: 0.6,
                    flat_width_mm: 0.2,
                }
            },
            crest_radius_limit_mm: limit,
        },
        1e-7,
        1e-6,
    )
    .unwrap()
}
fn overlap(a: [f64; 2], b: [f64; 2]) {
    assert!(a[0] <= b[1] && b[0] <= a[1], "{a:?} vs {b:?}");
}
fn integral(thread: &AuthoredThread, z: f64) -> I {
    // Independent composite Simpson integration of the ORIGINAL profile jet.
    // For these fixtures r<6,R<1,beta<2: each arc integrand is a finite
    // trig polynomial of frequency<=3beta with fourth derivative<1e6.
    // The line integrands are quadratic. Use that full conservative remainder.
    let n = 2048;
    let mut result = I::point(0.);
    for band in 0..thread.band_count() {
        let mut sum = I::point(0.);
        for i in 0..=n {
            let (q, a) = thread
                .profile(band, Series::<2>::variable(I::point(i as f64 / n as f64)))
                .unwrap();
            let value =
                (a.c[0] + I::point(thread.spec.radial_taper) * I::point(z)).square() * q.c[1];
            let weight = if i == 0 || i == n {
                1.
            } else if i % 2 == 0 {
                2.
            } else {
                4.
            };
            sum = sum + I::point(weight) * value;
        }
        result = result + sum / I::point(3. * n as f64);
    }
    if matches!(thread.spec.section, AuthoredSection::NptFlat { .. }) {
        let gap = I::point(thread.spec.pitch_mm)
            - (thread.endpoint(2, true).unwrap().0 - thread.endpoint(0, false).unwrap().0);
        let root = I::point(thread.spec.major_radius_mm) - I::point(thread.depth())
            + I::point(thread.spec.radial_taper) * I::point(z);
        result = result + gap * root.square();
    }
    let error = I::point(thread.band_count() as f64) * I::point(1e6)
        / I::point(180.)
        / I::point(n as f64).square().square();
    (result + I::new(-error.hi, error.hi)) * pi() / I::point(thread.spec.pitch_mm)
}

#[test]
fn actual_flat_and_rounded_original_laws_give_distinct_headed_and_tapped_material() {
    for rounded in [false, true] {
        for external in [false, true] {
            let thread = source(rounded, external, None);
            let kind = if external {
                Kind::HexHead {
                    across_flats_mm: 16.,
                    height_mm: 3.,
                }
            } else {
                Kind::HexNut {
                    across_flats_mm: 16.,
                }
            };
            let body = Material::new(thread, kind, Default::default(), None).unwrap();
            assert_eq!(body.thread().spec(), thread.spec());
            assert_eq!(body.kind(), kind);
            assert!(body.clipped_thread().is_none());
            let area = [
                integral(&thread, 0.),
                integral(&thread, 3.),
                integral(&thread, 6.),
            ];
            for (expected, actual) in [area[0], area[2]]
                .into_iter()
                .zip(body.original_cap_areas_mm2())
            {
                overlap(bounds(expected), actual);
            }
            // Squared radius is quadratic in z; Simpson in z is exact even
            // though the curved ORIGINAL phase integral is not polynomial.
            let original =
                I::point(6.) * (area[0] + I::point(4.) * area[1] + area[2]) / I::point(6.);
            overlap(bounds(original), body.original_thread_volume_mm3());
            let hex = I::new(body.hex_area_mm2()[0], body.hex_area_mm2()[1]);
            let expected = if external {
                original + hex * I::point(3.)
            } else {
                hex * I::point(6.) - original
            };
            overlap(bounds(expected), body.volume_mm3());
            let planes = body.planes().as_slice();
            assert_eq!(planes.len(), if external { 3 } else { 2 });
            assert!(
                planes
                    .iter()
                    .all(|p| p.area_mm2[0] > 0. && p.height_bounds_mm[0] <= p.height_bounds_mm[1])
            );
            assert_eq!(planes[0].outward_sign, -1);
            assert_eq!(planes.last().unwrap().outward_sign, 1);
            let end = if external { 9. } else { 6. };
            assert_eq!(body.bounds_mm()[2][0], 0.);
            assert!(body.bounds_mm()[2][1] >= end && body.bounds_mm()[2][1] - end < 1e-12);
            assert!(body.preparation_work_used() > 0 && body.preparation_work_used() < 16_000_000);
            assert!(body.retained_budget_bytes() < body.accounted_peak_bytes());
        }
    }
    let flat = Material::new(
        source(false, true, None),
        Kind::HexHead {
            across_flats_mm: 16.,
            height_mm: 3.,
        },
        Default::default(),
        None,
    )
    .unwrap();
    let rounded = Material::new(
        source(true, true, None),
        flat.kind(),
        Default::default(),
        None,
    )
    .unwrap();
    assert!((flat.volume_mm3()[0] - rounded.volume_mm3()[0]).abs() > 0.01);
    assert!(
        flat.thread()
            .derived_rounding_radius_bounds_mm()
            .unwrap()
            .is_none()
    );
    assert!(
        rounded
            .thread()
            .derived_rounding_radius_bounds_mm()
            .unwrap()
            .unwrap()[0]
            > 0.
    );
}

#[test]
fn actual_moving_npt_stock_clip_retains_source_and_matches_original_cap_integrals() {
    let thread = source(false, true, Some(5.1));
    let kind = Kind::HexHead {
        across_flats_mm: 16.,
        height_mm: 3.,
    };
    let body = Material::new(thread, kind, Default::default(), None).unwrap();
    let clip = body.clipped_thread().unwrap();
    assert_eq!(clip.source().spec(), thread.spec());
    for (face, area) in [StockClipFace::StartCap, StockClipFace::EndCap]
        .into_iter()
        .zip(body.original_cap_areas_mm2())
    {
        overlap(
            clip.face_measure(face, Default::default(), None)
                .unwrap()
                .area_mm2(),
            area,
        );
    }
    assert!(body.radial_enclosure_mm()[1] <= 5.1);
    let inactive = Material::new(
        source(false, true, Some(100.)),
        kind,
        Default::default(),
        None,
    )
    .unwrap();
    let none = Material::new(source(false, true, None), kind, Default::default(), None).unwrap();
    assert!(inactive.clipped_thread().is_none());
    assert_eq!(inactive.thread().spec().crest_radius_limit_mm, Some(100.));
    assert_eq!(
        inactive.original_thread_volume_mm3(),
        none.original_thread_volume_mm3()
    );
    assert!(matches!(
        Material::new(
            source(false, true, Some(4.9)),
            kind,
            Default::default(),
            None
        ),
        Err(Error::Unsupported)
    ));
}

#[test]
fn actual_hex_stock_and_original_bore_classify_whole_boxes_without_phase_guessing() {
    let head = Material::new(
        source(false, true, Some(5.1)),
        Kind::HexHead {
            across_flats_mm: 16.,
            height_mm: 3.,
        },
        Default::default(),
        None,
    )
    .unwrap();
    let nut = Material::new(
        source(true, false, None),
        Kind::HexNut {
            across_flats_mm: 16.,
        },
        Default::default(),
        None,
    )
    .unwrap();
    for (body, point, expected) in [
        (&head, [0., 0., 3.], Location::Inside),
        (&head, [7., 0., 7.], Location::Inside),
        (&head, [7., 0., 3.], Location::Outside),
        (&head, [9., 7., 7.], Location::Outside),
        (&head, [7., 0., 6.], Location::Unresolved),
        (&head, [0., 8., 7.], Location::Unresolved),
        (&head, [100., 0., 3.], Location::Outside),
        (&head, [0., 0., 0.], Location::Unresolved),
        (&nut, [0., 0., 3.], Location::Outside),
        (&nut, [7., 0., 3.], Location::Inside),
        (&nut, [9., 7., 3.], Location::Outside),
        (&nut, [4.5, 0., 0.5], Location::Unresolved),
    ] {
        let query = body
            .classify_box(point.map(|v| [v, v]), Default::default(), None)
            .unwrap();
        assert_eq!(query.location(), expected, "{point:?} {:?}", body.kind());
        assert!(query.query_work_used() > 0);
        assert_eq!(
            query.inherited_preparation_work_used(),
            body.preparation_work_used()
        );
        assert!(query.accounted_peak_bytes() <= 64 * 1024 * 1024);
    }
    assert_eq!(
        head.classify_box(
            [[-0.1, 0.1], [-0.1, 0.1], [5.9, 6.1]],
            Default::default(),
            None
        )
        .unwrap()
        .location(),
        Location::Inside
    );
}

#[test]
fn actual_long_clipped_source_preserves_positive_root_margin_at_both_caps() {
    let base = source(false, true, None).spec();
    let spec = AuthoredSpec {
        length_mm: 32.,
        crest_radius_limit_mm: Some(5.4002),
        ..base
    };
    let thread = AuthoredThread::new(spec, 1e-7, 1e-6).unwrap();
    let body = Material::new(
        thread,
        Kind::HexHead {
            across_flats_mm: 16.,
            height_mm: 3.,
        },
        Default::default(),
        None,
    )
    .unwrap();
    let plan = body.clipped_thread().unwrap();
    let end_root = I::point(spec.major_radius_mm) - I::point(thread.depth())
        + I::point(spec.radial_taper) * I::point(spec.length_mm);
    assert!((I::point(plan.stock_radius_limit_mm()) - end_root).lo > thread.minimum);
    assert!((I::point(plan.stock_radius_limit_mm()) - end_root).hi < 0.001);
    for (face, area) in [StockClipFace::StartCap, StockClipFace::EndCap]
        .into_iter()
        .zip(body.original_cap_areas_mm2())
    {
        overlap(
            plan.face_measure(face, Default::default(), None)
                .unwrap()
                .area_mm2(),
            area,
        );
    }
    for limit in [5.4, 5.3998] {
        assert!(matches!(
            AuthoredThread::new(
                AuthoredSpec {
                    crest_radius_limit_mm: Some(limit),
                    ..spec
                },
                1e-7,
                1e-6
            ),
            Err(Error::Unresolved)
        ));
    }
}

#[test]
fn real_source_hardware_refuses_invalid_role_precision_resources_and_cancellation() {
    let thread = source(false, true, None);
    let head = Kind::HexHead {
        across_flats_mm: 16.,
        height_mm: 3.,
    };
    for kind in [
        Kind::HexHead {
            across_flats_mm: f64::NAN,
            height_mm: 3.,
        },
        Kind::HexHead {
            across_flats_mm: 16.,
            height_mm: 0.,
        },
        Kind::HexNut {
            across_flats_mm: 0.,
        },
    ] {
        assert!(matches!(
            Material::new(thread, kind, Default::default(), None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        Material::new(
            thread,
            Kind::HexNut {
                across_flats_mm: 16.
            },
            Default::default(),
            None
        ),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        Material::new(source(false, false, None), head, Default::default(), None),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        Material::new(
            thread,
            Kind::HexHead {
                across_flats_mm: 4.,
                height_mm: 3.
            },
            Default::default(),
            None
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        Material::new(
            thread,
            Kind::HexHead {
                across_flats_mm: 16.,
                height_mm: 1e-12
            },
            Default::default(),
            None
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        Material::new(
            thread,
            head,
            Options {
                volume_width_mm3: 1e-30,
                ..Default::default()
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    for limits in [
        StockClipOptions {
            maximum_work: 1,
            ..Default::default()
        },
        StockClipOptions {
            maximum_work: 4096,
            ..Default::default()
        },
        StockClipOptions {
            maximum_bytes: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            Material::new(
                thread,
                head,
                Options {
                    limits,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        Material::new(thread, head, Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let body = Material::new(thread, head, Default::default(), None).unwrap();
    let query = [[0., 0.], [0., 0.], [3., 3.]];
    assert!(matches!(
        body.classify_box(query, Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        body.classify_box(
            query,
            StockClipOptions {
                maximum_work: 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        body.classify_box(
            query,
            StockClipOptions {
                maximum_bytes: 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        body.classify_box(
            [[f64::NAN, 0.], [0., 0.], [0., 0.]],
            Default::default(),
            None
        ),
        Err(Error::Invalid)
    ));
    assert_eq!(body.thread().spec(), thread.spec());
}
