use super::*;
#[test]
fn actual_same_body_repaired_heads_measure_all_raw_supports_without_material_admission() {
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
            assert!(std::ptr::eq(head.body(), &body));
            assert!(std::ptr::eq(head.common_mouth(), &common));
            assert!(std::ptr::eq(
                head.source().profile_proposal(),
                body.retained_meridian_header()
            ));
            assert_eq!(head.chart(), common.chart());
            assert_eq!(
                head.virtual_bottom_station(),
                body.retained_meridian_header().end()
            );
            assert_eq!(
                head.inherited_finite_matching_admission_mm(),
                body.retained_meridian_header().rounding_error_mm()
            );
            let bounds = head.support_transfers().collect::<Vec<_>>();
            assert_eq!(bounds.len(), body.profile_face_entities().len() + 1);
            assert_eq!(bounds[0].meridian_role, None);
            assert_eq!(
                bounds[0].material_face_entity,
                body.virtual_outer_cap().material_face
            );
            assert_eq!(bounds[0].support_entity, body.virtual_outer_cap().surface);
            for ((b, &id), s) in bounds[1..]
                .iter()
                .zip(body.profile_face_entities())
                .zip(body.retained_meridian_header().segments())
            {
                assert_eq!(b.meridian_role, Some(s.role));
                assert_eq!(b.material_face_entity, id);
            }
            for kind in [
                RawHeadSupportKind::Plane,
                RawHeadSupportKind::Cylinder,
                RawHeadSupportKind::Sphere,
                RawHeadSupportKind::Torus,
            ] {
                assert!(bounds.iter().any(|b| b.raw_kind == kind));
            }
            for b in &bounds {
                assert!(b.support_entity > 0 && b.distance_mm > 0. && b.distance_mm <= 1e-7);
            }
            assert!(head.work_used() > common.work_used() + head.source().work_used());
            assert_eq!(
                head.retained_budget_bytes(),
                common.retained_budget_bytes() + std::mem::size_of_val(&head)
            );
            eprintln!(
                "head transfer scale{scale} pose{posed}: {bounds:?}, work{} bytes{}",
                head.work_used(),
                head.retained_budget_bytes()
            );
            assert_eq!(body.literal_input(), text);
        }
    }
}
#[test]
fn actual_counter_and_cylindrical_h0_packets_retain_cones_and_virtual_caps() {
    use crate::head_profile::Head;
    for head_shape in [
        Head::Countersunk {
            radius_mm: 3.,
            height_mm: 1.65,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.1,
            neck_radius_mm: 0.1,
        },
        Head::Cylindrical {
            radius_mm: 2.,
            height_mm: 1.6,
        },
    ] {
        let text = super::super::super::closed_head::tests::additional_text(head_shape);
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
        let bounds = head.support_transfers().collect::<Vec<_>>();
        assert_eq!(bounds.len(), head.source().segment_count() + 1);
        assert!(
            bounds
                .iter()
                .all(|b| b.distance_mm > 0. && b.distance_mm <= 1e-7)
        );
        assert_eq!(
            bounds
                .iter()
                .any(|b| b.raw_kind == RawHeadSupportKind::Cone),
            matches!(head_shape, Head::Countersunk { .. })
        );
        assert!(std::ptr::eq(head.body(), &body));
        eprintln!("additional head{head_shape:?}: {bounds:?}");
    }
}
#[test]
fn same_body_head_support_limits_and_every_cancellation_prefix_publish_no_partial_packet() {
    let text = super::super::super::super::branches::representative::scaled_text(1.137, true);
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
    let good = common
        .closed_repaired_head(Default::default(), None)
        .unwrap();
    for options in [
        ClosedHeadOptions {
            maximum_distance_mm: f64::NAN,
            ..Default::default()
        },
        ClosedHeadOptions {
            maximum_distance_mm: 0.,
            ..Default::default()
        },
        ClosedHeadOptions {
            maximum_distance_mm: 1e-6,
            ..Default::default()
        },
        ClosedHeadOptions {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        ClosedHeadOptions {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            common.closed_repaired_head(options, None),
            Err(PhillipsReadError::Invalid)
        ));
    }
    for options in [
        ClosedHeadOptions {
            maximum_work: 0,
            ..Default::default()
        },
        ClosedHeadOptions {
            maximum_work: common.work_used() + 1,
            ..Default::default()
        },
        ClosedHeadOptions {
            maximum_work: good.work_used() - 1,
            ..Default::default()
        },
        ClosedHeadOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            common.closed_repaired_head(options, None),
            Err(PhillipsReadError::Limit)
        ));
    }
    assert!(
        common
            .closed_repaired_head(
                ClosedHeadOptions {
                    maximum_work: good.work_used(),
                    ..Default::default()
                },
                None
            )
            .is_ok()
    );
    assert!(matches!(
        common.closed_repaired_head(
            ClosedHeadOptions {
                maximum_distance_mm: 1e-30,
                ..Default::default()
            },
            None
        ),
        Err(PhillipsReadError::Unsupported)
    ));
    let source_error = (0..good.source().segment_count())
        .map(|i| {
            let c = good.source().segment_certificate(i).unwrap();
            c.support_transfer_mm + c.finite_trim_transfer_mm
        })
        .fold(0., f64::max);
    let actual = good
        .support_transfers()
        .map(|b| b.distance_mm)
        .fold(0., f64::max);
    assert!(actual > source_error);
    // This precision clears the independently measured local repair and fails
    // on actual whole raw-support transfer, rather than changing geometry.
    assert!(
        common
            .closed_repaired_head(
                ClosedHeadOptions {
                    maximum_distance_mm: (actual + source_error) / 2.,
                    ..Default::default()
                },
                None
            )
            .is_err()
    );
    let (mut lo, mut hi) = (1usize, BYTES);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        match common.closed_repaired_head(
            ClosedHeadOptions {
                maximum_bytes: mid,
                ..Default::default()
            },
            None,
        ) {
            Ok(_) => hi = mid,
            Err(PhillipsReadError::Limit) => lo = mid + 1,
            Err(e) => panic!("unexpected byte frontier error: {e:?}"),
        }
    }
    assert!(lo > good.retained_budget_bytes());
    assert!(
        common
            .closed_repaired_head(
                ClosedHeadOptions {
                    maximum_bytes: lo,
                    ..Default::default()
                },
                None
            )
            .is_ok()
    );
    assert!(matches!(
        common.closed_repaired_head(
            ClosedHeadOptions {
                maximum_bytes: lo - 1,
                ..Default::default()
            },
            None
        ),
        Err(PhillipsReadError::Limit)
    ));
    let calls = Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    common
        .closed_repaired_head(Default::default(), Some(&count))
        .unwrap();
    for phase in 1..=calls.get() {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() == phase
        };
        assert!(
            matches!(
                common.closed_repaired_head(Default::default(), Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ),
            "phase{phase}"
        );
    }
    let retry = common
        .closed_repaired_head(Default::default(), None)
        .unwrap();
    assert_eq!(retry.chart(), good.chart());
    assert_eq!(retry.work_used(), good.work_used());
    assert_eq!(
        retry
            .support_transfers()
            .map(|b| b.distance_mm)
            .collect::<Vec<_>>(),
        good.support_transfers()
            .map(|b| b.distance_mm)
            .collect::<Vec<_>>()
    );
    assert_eq!(body.literal_input(), text);
}
