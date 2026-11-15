use super::*;
#[test]
fn actual_periodic_owner_survives_every_borrowed_source_and_retains_actual_evaluators() {
    let mut result = None;
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let finite = head
            .periodic_thread_finite_curves(&field, Default::default(), None)
            .unwrap();
        let join = head
            .periodic_shared_crest_shoulder(&finite, Default::default(), None)
            .unwrap();
        let material = head
            .periodic_closed_material(&join, Default::default(), None)
            .unwrap();
        let volume = material.volume_bounds(Default::default(), None).unwrap();
        let expected = (
            head.source().station(0),
            head.source().station(head.source().segment_count()),
            field.source().phase_node_expressions_mm(),
            volume.local_volume_mm3(),
            volume.component_volumes_mm3(),
            volume.total_work_used(),
        );
        let peak = volume.owned_peak_budget_bytes();
        let owned = volume.owned_source(Default::default(), None).unwrap();
        assert!(owned.retained_budget_bytes() < peak);
        assert_eq!(owned.head().view().station(0), expected.0);
        assert_eq!(
            owned
                .head()
                .view()
                .station(owned.head().view().segment_count()),
            expected.1
        );
        assert_eq!(owned.thread().phase_node_expressions_mm(), expected.2);
        assert_eq!(owned.local_volume_mm3(), expected.3);
        assert_eq!(owned.component_volumes_mm3(), expected.4);
        assert_eq!(owned.inherited_work_used(), expected.5);
        assert_eq!(
            owned.total_work_used(),
            owned.inherited_work_used() + owned.new_work_used()
        );
        assert_eq!(owned.chart(), head.chart());
        assert_eq!(
            owned.head_bottom_station_mm(),
            join.head_bottom_station_mm()
        );
        assert_eq!(
            owned.cap_expressions_mm(),
            finite.cap_stations_mm().map(|z| [z, 0.])
        );
        assert!(owned.determinant_bounds()[0] > 0.);
        assert!(
            owned
                .containment_clearances_mm()
                .into_iter()
                .all(|v| v > 0.)
        );
        assert!(owned.world_volume_mm3()[0] > 0.);
        assert_eq!(owned.original_curve_owners().unwrap().len(), 17);
        assert_eq!(owned.thread_transfers().unwrap().len(), 4);
        assert!(
            owned
                .original_curve_owners()
                .unwrap()
                .iter()
                .all(|o| o.face_entity > 0
                    && o.edge_entity > 0
                    && o.curve_entity > 0
                    && o.vertex_entities.into_iter().all(|v| v > 0)
                    && o.maximum_distance_mm < 1e-7
                    && o.maximum_cell_derivative_mm < 1e-7
                    && o.endpoint_incidence_reserve_mm < 1e-7)
        );
        for i in 0..owned.face_count() {
            assert!(owned.face_source(i).is_some());
        }
        assert!(owned.face_source(owned.face_count()).is_none());
        assert!(owned.face_source(usize::MAX).is_none());
        result = Some(owned);
    });
    let owned = result.unwrap();
    assert!(owned.original_input().unwrap().contains("ISO-10303-21"));
    let e = owned.local_enclosure_bounds();
    assert!(
        e.into_iter()
            .all(|b| b[0].is_finite() && b[1].is_finite() && b[0] < b[1])
    );
    assert!(owned.cavity().volume_bounds_mm3(1e-5, 100_000).is_ok());
    let bounds = owned
        .thread()
        .bounds(1, [0.25, 0.75], [1., 1.5], 16_000_000, None)
        .unwrap();
    assert!(
        bounds
            .position_mm
            .into_iter()
            .all(|b| b[0].is_finite() && b[1].is_finite())
    );
    eprintln!(
        "OWNED PERIODIC retained{} inherited{} new{} known{} faces{}",
        owned.retained_budget_bytes(),
        owned.inherited_work_used(),
        owned.new_work_used(),
        owned.total_work_used(),
        owned.face_count()
    );
}
#[test]
fn actual_periodic_ownership_obeys_new_total_and_live_peak_budgets() {
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let finite = head
            .periodic_thread_finite_curves(&field, Default::default(), None)
            .unwrap();
        let join = head
            .periodic_shared_crest_shoulder(&finite, Default::default(), None)
            .unwrap();
        let material = head
            .periodic_closed_material(&join, Default::default(), None)
            .unwrap();
        let volume = material.volume_bounds(Default::default(), None).unwrap();
        let o = PeriodicOwnedOptions::default();
        let good = volume.owned_source(o, None).unwrap();
        let peak = volume.owned_peak_budget_bytes();
        let used = good.new_work_used();
        let known = good.total_work_used();
        drop(good);
        for bad in [
            PeriodicOwnedOptions {
                maximum_work: WORK + 1,
                ..o
            },
            PeriodicOwnedOptions {
                maximum_total_work: TOTAL + 1,
                ..o
            },
            PeriodicOwnedOptions {
                maximum_bytes: BYTES + 1,
                ..o
            },
        ] {
            assert!(matches!(
                volume.owned_source(bad, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        for bad in [
            PeriodicOwnedOptions {
                maximum_work: 0,
                ..o
            },
            PeriodicOwnedOptions {
                maximum_work: 1,
                ..o
            },
            PeriodicOwnedOptions {
                maximum_work: used - 1,
                ..o
            },
            PeriodicOwnedOptions {
                maximum_total_work: volume.total_work_used() - 1,
                ..o
            },
            PeriodicOwnedOptions {
                maximum_total_work: volume.total_work_used(),
                ..o
            },
            PeriodicOwnedOptions {
                maximum_total_work: known - 1,
                ..o
            },
            PeriodicOwnedOptions {
                maximum_bytes: peak - 1,
                ..o
            },
        ] {
            assert!(matches!(
                volume.owned_source(bad, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        assert!(
            volume
                .owned_source(
                    PeriodicOwnedOptions {
                        maximum_work: used,
                        maximum_total_work: known,
                        maximum_bytes: peak
                    },
                    None
                )
                .is_ok()
        );
        let calls = Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        volume.owned_source(o, Some(&count)).unwrap();
        let complete = calls.get();
        assert!(complete >= 5);
        for stop in 1..=complete {
            calls.set(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() == stop
            };
            assert!(
                matches!(
                    volume.owned_source(o, Some(&cancel)),
                    Err(PhillipsReadError::Cancelled)
                ),
                "callback{stop}/{complete}"
            );
        }
        assert!(volume.owned_source(o, None).is_ok());
    });
}

#[test]
fn actual_owned_periodic_export_proposal_has_independent_file_derived_material() {
    let original = degree11::tests::fixture::actual_source();
    let source = OwnedPeriodicMaterial::from_authored(&original, Default::default(), None).unwrap();
    let proposal = source.authored_mesh_approximation(WORK, None).unwrap();
    let now = std::time::Instant::now();
    let text = crate::step::write_phillips_screw(
        proposal.proposal(),
        [0.; 3],
        [1., 0., 0.],
        [0., 1., 0.],
        1e-7,
    )
    .unwrap();
    eprintln!(
        "PERIODIC STEP proposal {}bytes {:?}",
        text.len(),
        now.elapsed()
    );
    let mut decoded = None;
    crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head(&text,&mut|head| {
        let field=head.common_chart_periodic_thread_repair(Default::default(),None).unwrap();
        let finite=head.periodic_thread_finite_curves(&field,Default::default(),None).unwrap();
        let join=head.periodic_shared_crest_shoulder(&finite,Default::default(),None).unwrap();
        let material=head.periodic_closed_material(&join,Default::default(),None).unwrap();
        let volume=material.volume_bounds(Default::default(),None).unwrap();
        decoded=Some(volume.owned_source(Default::default(),None).unwrap());
    });
    let decoded = decoded.unwrap();
    let thread_match = source
        .thread()
        .finite_correspondence(
            decoded.thread(),
            source.cap_expressions_mm(),
            decoded.cap_expressions_mm(),
            WORK,
            None,
        )
        .unwrap();
    eprintln!("actual decoded periodic field correspondence {thread_match:?}");
    assert!(
        thread_match.position_mm < 1e-7 && thread_match.du_mm < 1e-7 && thread_match.dv_mm < 1e-7
    );
    let ca = source.cavity().surface_patches().unwrap();
    let cb = decoded.cavity().surface_patches().unwrap();
    let cavity_match = ca
        .source_correspondence(&cb, Default::default(), None)
        .unwrap();
    eprintln!("actual decoded cavity correspondence {:?}", cavity_match);

    let head_match = source
        .head()
        .view()
        .translated_correspondence(
            &decoded.head().view(),
            source.head_z_shift_mm(),
            decoded.head_z_shift_mm(),
            WORK,
            None,
        )
        .unwrap();
    eprintln!("actual decoded head correspondence {:?}", head_match);
    assert!(head_match.position_mm < 1e-7 && head_match.derivative_mm < 1e-7);
    eprintln!(
        "STEP CAVITY specs {:?}/{:?} domain {:?}/{:?}",
        source.cavity().spec(),
        decoded.cavity().spec(),
        source.cavity().domain(),
        decoded.cavity().domain()
    );
    assert_eq!(decoded.original_input(), Some(text.as_str()));
    assert!(
        decoded
            .original_curve_owners()
            .unwrap()
            .iter()
            .all(|c| c.face_entity > 0 && c.edge_entity > 0 && c.curve_entity > 0)
    );
    assert!(decoded.local_volume_mm3()[0] > 0.);
    eprintln!(
        "PERIODIC STEP actual-file derived {:?} knownwork{} bytes{} elapsed{:?}\nsource chart{:?} headShift{} mouth{:?} caps{:?} head {:?}/{:?}\ndecoded chart{:?} shift{} mouth{:?} caps{:?} head {:?}/{:?}\nradii {:?}/{:?}\nphase {:?}/{:?}",
        decoded.world_volume_mm3(),
        decoded.total_work_used(),
        decoded.retained_budget_bytes(),
        now.elapsed(),
        source.chart(),
        source.head_z_shift_mm(),
        source.cavity_mouth_expression_mm(),
        source.cap_expressions_mm(),
        source.head().view().station(0),
        source.head_bottom_station_mm(),
        decoded.chart(),
        decoded.head_z_shift_mm(),
        decoded.cavity_mouth_expression_mm(),
        decoded.cap_expressions_mm(),
        decoded.head().view().station(0),
        decoded.head_bottom_station_mm(),
        source.thread().radius_nodes_mm(),
        decoded.thread().radius_nodes_mm(),
        source.thread().phase_node_expressions_mm(),
        decoded.thread().phase_node_expressions_mm()
    );
}
