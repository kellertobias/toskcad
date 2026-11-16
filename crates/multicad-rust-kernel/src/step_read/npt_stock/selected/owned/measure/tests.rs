use super::*;
use crate::step_read::owned_insert::DocumentOptions;
fn source(quarter: bool, mapped: bool) -> OwnedMaterial {
    let text = crate::step_read::npt_stock::tests::fixture(quarter);
    let text = if mapped {
        crate::step_read::npt_stock::selected::tests::mapped(&text, "(-3.,0.,0.)", "1.")
    } else {
        text
    };
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap()
}
#[test]
fn actual_literal_cap_areas_and_centroids_survive_source_destruction_and_mapping() {
    for (quarter, mapped) in [(false, false), (true, false), (false, true)] {
        let source = source(quarter, mapped);
        for (face, owner) in source.faces().iter().enumerate() {
            if !matches!(owner.support(), Support::Plane { .. }) {
                continue;
            }
            let result = source
                .cap_face_measure(face, Default::default(), None)
                .unwrap();
            assert_eq!(result.face_index(), face);
            assert!(Arc::ptr_eq(&result.source().0, &source.0));
            assert!(result.area_mm2()[0] > 0.);
            assert!(result.work_used() < FaceMeasureOptions::default().maximum_work);
            assert!(result.retained_budget_bytes() < FaceMeasureOptions::default().maximum_bytes);
            for (c, b) in result
                .centroid_mm()
                .into_iter()
                .zip(source.world_enclosure_mm())
            {
                assert!(c[0] >= b[0] && c[1] <= b[1]);
            }
            let limited = FaceMeasureOptions {
                maximum_work: result.work_used() - 1,
                ..Default::default()
            };
            assert!(matches!(
                source.cap_face_measure(face, limited, None),
                Err(Error::Limit)
            ));
            assert!(matches!(
                source.cap_face_measure(face, Default::default(), Some(&|| true)),
                Err(Error::Cancelled)
            ));
            let retained = result.source().clone();
            drop(result);
            drop(retained);
        }
    }
}

#[test]
fn literal_green_moments_enclose_an_independent_parabolic_region() {
    let mut work = Work {
        left: 1_000_000,
        used: 0,
        bytes: 64 * 1024,
        interrupt: None,
    };
    let weights = CapWeights::new(&mut work).unwrap();
    // Degree-elevated x=t,y=t², closed by the reversed line y=x.
    // Direct elementary integration gives A=1/6, Mx=1/12, My=1/15.
    let parabola = std::array::from_fn(|i| {
        [
            p(i as f64) / p(11.),
            p((i * (i.saturating_sub(1))) as f64) / p(110.),
            p(0.),
        ]
    });
    let line = std::array::from_fn(|i| [p(i as f64) / p(11.), p(i as f64) / p(11.), p(0.)]);
    let a = weights.curve(parabola, &mut work).unwrap();
    let b = weights.curve(line, &mut work).unwrap();
    for (value, denominator) in a.into_iter().zip(b).map(|(a, b)| a - b).zip([6, 12, 15]) {
        let expected = num_rational::BigRational::new(1.into(), denominator.into());
        assert!(num_rational::BigRational::from_float(value.lo).unwrap() <= expected);
        assert!(num_rational::BigRational::from_float(value.hi).unwrap() >= expected);
        assert!(value.width() < 1e-10);
    }
}

#[test]
fn actual_cap_measure_refuses_bad_settings_and_cancels_before_publication() {
    use std::cell::Cell;
    let source = source(false, false);
    let face = source
        .faces()
        .iter()
        .position(|f| matches!(f.support(), Support::Plane { .. }))
        .unwrap();
    for options in [
        FaceMeasureOptions {
            maximum_work: 0,
            ..Default::default()
        },
        FaceMeasureOptions {
            area_width_mm2: f64::NAN,
            ..Default::default()
        },
        FaceMeasureOptions {
            centroid_width_mm: 0.,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            source.cap_face_measure(face, options, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        source.cap_face_measure(usize::MAX, Default::default(), None),
        Err(Error::Invalid)
    ));
    let lateral = source
        .faces()
        .iter()
        .position(|f| matches!(f.support(), Support::Tensor(_)))
        .unwrap();
    assert!(matches!(
        source.cap_face_measure(lateral, Default::default(), None),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        source.cap_face_measure(
            face,
            FaceMeasureOptions {
                maximum_bytes: 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        source.cap_face_measure(
            face,
            FaceMeasureOptions {
                area_width_mm2: 1e-30,
                ..Default::default()
            },
            None
        ),
        Err(Error::Unsupported)
    ));
    let count = Cell::new(0);
    let original = source
        .cap_face_measure(
            face,
            Default::default(),
            Some(&|| {
                count.set(count.get() + 1);
                false
            }),
        )
        .unwrap();
    for stop in [1, count.get() / 2, count.get() - 1] {
        let n = Cell::new(0);
        assert!(matches!(
            source.cap_face_measure(
                face,
                Default::default(),
                Some(&|| {
                    n.set(n.get() + 1);
                    n.get() == stop
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    let again = source
        .cap_face_measure(face, Default::default(), None)
        .unwrap();
    assert_eq!(original.area_mm2(), again.area_mm2());
    assert_eq!(original.centroid_mm(), again.centroid_mm());
}

#[test]
fn actual_original_literal_surface_face_measure_uses_prepared_normal_owner() {
    use std::cell::Cell;
    let mut baseline = Vec::new();
    for (quarter, mapped) in [(false, false), (true, false), (false, true)] {
        let source = source(quarter, mapped);
        let plan = source.geometry_plan(Default::default(), None).unwrap();
        for (face, _) in source.faces().iter().enumerate() {
            let started = std::time::Instant::now();
            let result = plan.face_measure(face, Default::default(), None).unwrap();
            eprintln!(
                "literal quarter={quarter} mapped={mapped} face={face} area={:?} work={} peak={} elapsed={:?}",
                result.area_mm2(),
                result.work_used(),
                result.accounted_peak_bytes(),
                started.elapsed()
            );
            if !quarter && !mapped {
                baseline.push((result.area, result.centroid));
            }
            if mapped {
                let (area, c) = baseline[face];
                assert!(result.area.lo <= area.hi && area.lo <= result.area.hi);
                let expected = [p(11.) - c[1], p(-5.) + c[0], p(3.) + c[2]];
                for (actual, expected) in result.centroid.into_iter().zip(expected) {
                    assert!(actual.lo <= expected.hi && expected.lo <= actual.hi);
                }
            }
            assert_eq!(result.face_index(), face);
            let area = interval(result.area_mm2());
            for ((anchor, moment), centroid) in result
                .moment_anchor_mm()
                .into_iter()
                .zip(result.first_moments_about_anchor_mm3())
                .zip(result.centroid_mm())
            {
                let recovered = p(anchor) + interval(moment) / area;
                assert!(
                    recovered.lo <= centroid[1] && centroid[0] <= recovered.hi,
                    "public centered moments must reconstruct the reported centroid"
                );
            }

            assert!(Arc::ptr_eq(&result.source().0, &source.0));
            assert!(result.area_mm2()[0] > 0.);
            assert!(
                result.area_mm2()[1] - result.area_mm2()[0]
                    <= FaceMeasureOptions::default().area_width_mm2
            );
            assert!(
                result
                    .centroid_mm()
                    .iter()
                    .all(|v| v[1] - v[0] <= FaceMeasureOptions::default().centroid_width_mm)
            );
            assert!(result.retained_budget_bytes() < result.accounted_peak_bytes());
            assert!(matches!(
                plan.face_measure(
                    face,
                    FaceMeasureOptions {
                        maximum_bytes: result.accounted_peak_bytes() - 1,
                        ..Default::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
            if face == 0 && !quarter && !mapped {
                assert!(matches!(
                    plan.face_measure(
                        face,
                        FaceMeasureOptions {
                            maximum_work: result.work_used() - 1,
                            ..Default::default()
                        },
                        None
                    ),
                    Err(Error::Limit)
                ));
                let calls = Cell::new(0usize);
                assert!(matches!(
                    plan.face_measure(
                        face,
                        Default::default(),
                        Some(&|| {
                            calls.set(calls.get() + 1);
                            calls.get() == 10000
                        })
                    ),
                    Err(Error::Cancelled)
                ));
                assert_eq!(calls.get(), 10000);
            }
        }
        assert!(matches!(
            plan.face_measure(usize::MAX, Default::default(), None),
            Err(Error::Invalid)
        ));
        for face in [
            0,
            source
                .faces()
                .iter()
                .position(|f| matches!(f.support(), Support::Plane { .. }))
                .unwrap(),
        ] {
            for options in [
                FaceMeasureOptions {
                    maximum_bytes: 64 * 1024 * 1024 + 1,
                    ..Default::default()
                },
                FaceMeasureOptions {
                    maximum_work: 0,
                    ..Default::default()
                },
                FaceMeasureOptions {
                    area_width_mm2: f64::NAN,
                    ..Default::default()
                },
                FaceMeasureOptions {
                    centroid_width_mm: 0.,
                    ..Default::default()
                },
            ] {
                assert!(matches!(
                    plan.face_measure(face, options, None),
                    Err(Error::Invalid)
                ));
            }
            assert!(matches!(
                plan.face_measure(face, Default::default(), Some(&|| true)),
                Err(Error::Cancelled)
            ));
        }
    }
}

#[test]
fn actual_far_mapped_surface_refuses_unrepresentable_world_centroid_precision() {
    let original = crate::step_read::npt_stock::tests::fixture(false);
    let mapped =
        crate::step_read::npt_stock::selected::tests::mapped(&original, "(-3.,0.,0.)", "1.");
    let text = mapped.replacen("(11.,-5.,3.)", "(9999990.,-5.,3.)", 1);
    assert_ne!(text, mapped);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let source = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap();
    let plan = source.geometry_plan(Default::default(), None).unwrap();
    assert!(matches!(source.faces()[6].support(), Support::Tensor(_)));
    let ordinary = plan.face_measure(6, Default::default(), None).unwrap();
    assert!(ordinary.centroid_mm()[0][0] > 9_999_980.);
    eprintln!(
        "far centroid ordinary {:?} work {}",
        ordinary.centroid_mm(),
        ordinary.work_used()
    );
    let achievable = plan
        .face_measure(
            6,
            FaceMeasureOptions {
                centroid_width_mm: 1e-8,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    assert!(achievable.centroid_mm().iter().all(|v| v[1] - v[0] <= 1e-8));
    let strict = FaceMeasureOptions {
        centroid_width_mm: 1e-9,
        ..Default::default()
    };
    assert!(matches!(
        plan.face_measure(6, strict, None),
        Err(Error::Unsupported)
    ));
    assert!(plan.face_measure(6, Default::default(), None).is_ok());
}
