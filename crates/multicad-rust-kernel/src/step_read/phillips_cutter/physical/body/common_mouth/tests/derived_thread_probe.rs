use super::*;
#[test]
fn actual_same_body_repaired_thread_source_attempts_existing_occupied_gate() {
    for scale in [1., 1.137] {
        for posed in [false, true] {
            let text =
                super::super::super::super::branches::representative::scaled_text(scale, posed);
            let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
            let outer = body
                .outer_carrier_bounds(
                    OuterCarrierOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 10_000_000,
                    },
                    None,
                )
                .unwrap();
            let (lifts, cubics) = closed_inputs(&body);
            let sewn = lifts
                .sewn_edge_sources(&body, &cubics, Default::default(), None)
                .unwrap();
            let rep = body
                .derived_cavity_representative(&sewn, Default::default(), None)
                .unwrap();
            let closed = rep
                .closed_boundary_certificate(Default::default(), None)
                .unwrap();
            let common = closed
                .common_mouth_head_envelope(&outer, Default::default(), None)
                .unwrap();
            let head = common
                .closed_repaired_head(Default::default(), None)
                .unwrap();
            let source = head
                .derived_thread_carrier(Default::default(), None)
                .unwrap();
            eprintln!(
                "THREADSOURCE scale{scale} pose{posed} work{} bytes{} radius{} caps{:?} bounds{:?}",
                source.work_used(),
                source.retained_budget_bytes(),
                source.entry_radius_upper_mm(),
                source.cap_stations_mm(),
                source.transfers()
            );
            assert!(std::ptr::eq(source.body(), &body));
            assert_eq!(source.chart(), head.chart());
            for (i, b) in source.transfers().iter().enumerate() {
                assert_eq!(b.source_face_index, 2 + i);
                assert_eq!(b.original_v_start_cell, 1);
                assert_eq!(
                    b.original_rail_entities[1],
                    source.transfers()[(i + 1) % 4].original_rail_entities[0]
                );
            }
            let occupied = source.certify_occupied(1e-7, 1e-5, None).unwrap();
            assert!(std::ptr::eq(occupied.body(), &body));
            assert_eq!(occupied.boundary().side_faces().len(), 4);
            eprintln!(
                "OCCUPIED scale{scale} pose{posed} error{}",
                occupied.boundary().boundary_error_mm()
            );
            assert_eq!(body.literal_input(), text);
        }
    }
}
pub(in crate::step_read::phillips_cutter::physical::body) fn with_head(
    text: &str,
    f: &mut dyn FnMut(&ClosedRepairedHead<'_, '_, '_, '_, '_>),
) {
    let body = read_phillips_body(text, 1e-7, 1_000_000).unwrap();
    let outer = body
        .outer_carrier_bounds(
            OuterCarrierOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let (lifts, cubics) = closed_inputs(&body);
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, Default::default(), None)
        .unwrap();
    let rep = body
        .derived_cavity_representative(&sewn, Default::default(), None)
        .unwrap();
    let closed = rep
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    let common = closed
        .common_mouth_head_envelope(&outer, Default::default(), None)
        .unwrap();
    let head = common
        .closed_repaired_head(Default::default(), None)
        .unwrap();
    f(&head);
    assert_eq!(body.literal_input(), text);
}
fn ordinary_text() -> String {
    super::super::super::super::branches::representative::scaled_text(1., false)
}
#[test]
fn repaired_thread_query_limits_precision_and_cancellation_preserve_original_source() {
    let text = ordinary_text();
    with_head(&text, &mut |head| {
        let defaults = DerivedThreadOptions::default();
        let baseline = head.derived_thread_carrier(defaults, None).unwrap();
        let expected = baseline.source().clone();
        let bounds = baseline.transfers().to_vec();
        let work = baseline.work_used();
        for options in [
            DerivedThreadOptions {
                maximum_distance_mm: 0.,
                ..defaults
            },
            DerivedThreadOptions {
                maximum_distance_mm: f64::NAN,
                ..defaults
            },
            DerivedThreadOptions {
                maximum_distance_mm: 1e-6,
                ..defaults
            },
            DerivedThreadOptions {
                maximum_cell_derivative_mm: 0.,
                ..defaults
            },
            DerivedThreadOptions {
                maximum_cell_derivative_mm: f64::NAN,
                ..defaults
            },
            DerivedThreadOptions {
                maximum_cell_derivative_mm: 1e-6,
                ..defaults
            },
            DerivedThreadOptions {
                maximum_controls: 262_145,
                ..defaults
            },
            DerivedThreadOptions {
                maximum_work: 16_000_001,
                ..defaults
            },
            DerivedThreadOptions {
                maximum_bytes: 64 * 1024 * 1024 + 1,
                ..defaults
            },
        ] {
            assert!(matches!(
                head.derived_thread_carrier(options, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        for options in [
            DerivedThreadOptions {
                maximum_distance_mm: 1e-12,
                ..defaults
            },
            DerivedThreadOptions {
                maximum_cell_derivative_mm: 1e-12,
                ..defaults
            },
        ] {
            assert!(matches!(
                head.derived_thread_carrier(options, None),
                Err(PhillipsReadError::Unsupported)
            ));
        }
        assert!(matches!(
            head.derived_thread_carrier(
                DerivedThreadOptions {
                    maximum_work: work - 1,
                    ..defaults
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        assert_eq!(
            head.derived_thread_carrier(
                DerivedThreadOptions {
                    maximum_work: work,
                    ..defaults
                },
                None
            )
            .unwrap()
            .source(),
            &expected
        );
        assert!(matches!(
            head.derived_thread_carrier(
                DerivedThreadOptions {
                    maximum_bytes: 0,
                    ..defaults
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        assert!(matches!(
            head.derived_thread_carrier(
                DerivedThreadOptions {
                    maximum_controls: 1,
                    ..defaults
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        let base = head.retained_budget_bytes()
            + std::mem::size_of::<DerivedThreadCarrier<'_, '_>>()
            + 1024 * 1024;
        for extra in [0, 1, 200, 500, 5000, 10000] {
            assert!(matches!(
                head.derived_thread_carrier(
                    DerivedThreadOptions {
                        maximum_bytes: base + extra,
                        ..defaults
                    },
                    None
                ),
                Err(PhillipsReadError::Limit)
            ));
        }
        let mut lo = 0;
        let mut hi = defaults.maximum_bytes;
        while lo + 1 < hi {
            let mid = lo + (hi - lo) / 2;
            match head.derived_thread_carrier(
                DerivedThreadOptions {
                    maximum_bytes: mid,
                    ..defaults
                },
                None,
            ) {
                Ok(_) => hi = mid,
                Err(PhillipsReadError::Limit) => lo = mid,
                Err(e) => panic!("unexpected {e:?}"),
            }
        }
        assert!(matches!(
            head.derived_thread_carrier(
                DerivedThreadOptions {
                    maximum_bytes: hi - 1,
                    ..defaults
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        let bounded = head
            .derived_thread_carrier(
                DerivedThreadOptions {
                    maximum_bytes: hi,
                    ..defaults
                },
                None,
            )
            .unwrap();
        assert_eq!(bounded.source(), &expected);
        assert!(bounded.retained_budget_bytes() + 1024 * 1024 <= hi);
        let calls = Cell::new(0usize);
        let total = head
            .derived_thread_carrier(
                defaults,
                Some(&|| {
                    calls.set(calls.get() + 1);
                    false
                }),
            )
            .unwrap();
        assert_eq!(total.source(), &expected);
        let total_calls = calls.get();
        let mut checkpoints = (1..=8)
            .chain((0..=32).map(|i| 1 + (total_calls - 1) * i / 32))
            .chain(total_calls - 7..=total_calls)
            .collect::<Vec<_>>();
        checkpoints.sort_unstable();
        checkpoints.dedup();
        for at in checkpoints {
            let calls = Cell::new(0usize);
            let result = head.derived_thread_carrier(
                defaults,
                Some(&|| {
                    calls.set(calls.get() + 1);
                    calls.get() == at
                }),
            );
            assert!(
                matches!(result, Err(PhillipsReadError::Cancelled)),
                "at{at}/{total_calls}"
            );
        }
        let retry = head.derived_thread_carrier(defaults, None).unwrap();
        assert_eq!(retry.source(), &expected);
        assert_eq!(retry.transfers().len(), bounds.len());
        let err = retry.certify_occupied(1e-7, 1e-5, Some(&|| true));
        assert!(matches!(err, Err(crate::tapered_thread::Error::Cancelled)));
        let source = head.derived_thread_carrier(defaults, None).unwrap();
        assert!(matches!(
            source.certify_occupied(0., 1e-5, None),
            Err(crate::tapered_thread::Error::Invalid)
        ));
        let positive = head
            .derived_thread_carrier(defaults, None)
            .unwrap()
            .certify_occupied(1e-7, 1e-5, None)
            .unwrap();
        assert_eq!(positive.chart(), head.chart());
        assert_eq!(
            positive.cap_stations_mm()[0],
            head.virtual_bottom_station()[1]
        );
        assert_eq!(positive.transfers().len(), 4);
        eprintln!("THREAD LIMITS exactwork{work} minbytes{hi} callbacks{total_calls}");
    });
}
#[test]
fn actual_identical_control_with_distinct_entity_keeps_shared_geometric_rail() {
    let text = ordinary_text();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, body.tolerance, true);
    let (surface, _, _) = helical::face(&shell, body.tensor_face_entities()[0]).unwrap();
    let entity = shell.simple(surface).unwrap();
    let rows = list(&entity.args[3]).unwrap();
    assert!(matches!(list(&entity.args[1]), Err(ReadError::Invalid)));
    let ids = refs(&rows[0]).unwrap();
    let id = ids[22];
    let new_id = body.entities.keys().next_back().unwrap() + 1;
    let point = shell.simple(id).unwrap().source.unwrap();
    let old = entity.source.unwrap();
    let new = old.replace(&format!("#{id},"), &format!("#{new_id},"));
    assert_ne!(old, new);
    let changed = text.replacen(old, &new, 1).replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("#{new_id}={point};\nENDSEC;\nEND-ISO-10303-21;"),
        1,
    );
    with_head(&changed, &mut |head| {
        let source = head
            .derived_thread_carrier(Default::default(), None)
            .unwrap();
        let occupied = source.certify_occupied(1e-7, 1e-5, None).unwrap();
        assert_eq!(occupied.boundary().side_faces().len(), 4);
    });
}
#[test]
fn actual_additional_h0_heads_and_multistarts_keep_source_domains() {
    for shape in [
        crate::head_profile::Head::Cylindrical {
            radius_mm: 2.,
            height_mm: 1.6,
        },
        crate::head_profile::Head::Countersunk {
            radius_mm: 3.,
            height_mm: 1.65,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.1,
            neck_radius_mm: 0.1,
        },
    ] {
        let text = super::super::super::closed_head::tests::additional_text(shape);
        with_head(&text, &mut |head| {
            let source = head
                .derived_thread_carrier(Default::default(), None)
                .unwrap();
            assert_eq!(
                source.cap_stations_mm()[0],
                head.virtual_bottom_station()[1]
            );
            source.certify_occupied(1e-7, 1e-5, None).unwrap();
        });
    }
    for (left, starts) in [(false, 2), (true, 2), (false, 8), (true, 8)] {
        let screw = super::super::super::super::tests::fixture(0, left, false, starts);
        let text =
            crate::step::write_phillips_screw(&screw, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
                .unwrap();
        with_head(&text, &mut |head| {
            let source = head
                .derived_thread_carrier(Default::default(), None)
                .unwrap();
            assert_eq!(source.transfers().len(), 4 * starts as usize);
            let occupied = source.certify_occupied(1e-7, 1e-5, None).unwrap();
            assert_eq!(occupied.boundary().side_faces().len(), 4 * starts as usize);
        });
    }
    let nominal = super::super::super::super::tests::fixture(0, false, true, 1);
    let text =
        crate::step::write_phillips_screw(&nominal, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    with_head(&text, &mut |head| {
        assert!(matches!(
            head.derived_thread_carrier(Default::default(), None),
            Err(PhillipsReadError::Unsupported)
        ));
    });
}
