use super::*;
fn closed_inputs<'a>(body: &'a Body<'_>) -> (AnalyticEdgeLifts<'a>, CubicSupportBounds<'a>) {
    let b = body.cavity();
    let ends = b
        .endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 1_000_000,
            },
            None,
        )
        .unwrap();
    let analytic = b
        .analytic_support_bounds(
            &ends,
            AnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let lifts = analytic
        .finite_edge_lifts(ends, AnalyticEdgeLiftOptions::default(), None)
        .unwrap();
    let cubics = b
        .cubic_support_bounds(
            CubicSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    (lifts, cubics)
}
#[test]
fn actual_common_chart_carriers_and_complete_head_envelope_are_qualified_without_sewing_claims() {
    for scale in [1., 1.137] {
        for posed in [false, true] {
            let text = super::super::super::branches::representative::scaled_text(scale, posed);
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
            let receipt = closed
                .common_mouth_head_envelope(&outer, Default::default(), None)
                .unwrap();
            assert!(std::ptr::eq(receipt.closed_cavity(), &closed));
            assert_eq!(receipt.chart().columns, closed.actual_affine_columns());
            assert_eq!(
                receipt.chart().outer_top_z,
                body.retained_meridian_header().top()[1]
            );
            assert_eq!(receipt.tensor_transfers().len(), outer.tensors().len());
            for (a, b) in receipt.tensor_transfers().iter().zip(outer.tensors()) {
                assert_eq!(
                    (
                        a.face_entity,
                        a.surface_entity,
                        a.control_rows,
                        a.control_columns
                    ),
                    (
                        b.face_entity,
                        b.surface_entity,
                        b.control_rows,
                        b.control_columns
                    )
                );
                assert!(a.distance_mm > 0. && a.distance_mm < 1e-7);
            }
            assert!(
                receipt.raw_mouth_plane_distance_mm() > 0.
                    && receipt.raw_mouth_plane_distance_mm() < 1e-7
            );
            assert!(
                receipt.whole_head_frame_transfer_mm() > 0.
                    && receipt.whole_head_frame_transfer_mm() < 1e-7
            );
            assert_eq!(
                receipt.inherited_head_reconstruction_reserve_mm(),
                outer.meridian_reconstruction_error_mm()
            );
            let head = receipt.proposed_head_envelope();
            assert!(std::ptr::eq(head.source(), rep.derived_cavity_source()));
            assert!(std::ptr::eq(
                head.profile_proposal(),
                body.retained_meridian_header()
            ));
            assert!(
                head.minimum_radial_clearance_mm() > 1e-7 && head.axial_floor_clearance_mm() > 1e-7
            );
            assert!(head.cells() > 0);
            assert!(head.work_used() > closed.local_source_certificate().work_used());
            assert!(head.retained_budget_bytes() > std::mem::size_of_val(head));
            assert!(
                receipt.work_used() > closed.work_used() + outer.work_used()
                    && receipt.work_used() < WORK
            );
            assert_eq!(
                receipt.retained_budget_bytes(),
                closed.retained_budget_bytes() + std::mem::size_of_val(&receipt)
            );
            assert_eq!(
                receipt
                    .closed_cavity()
                    .representative()
                    .body()
                    .literal_input(),
                text
            );
            eprintln!(
                "common scale{scale} pose{posed} rawPlane={} headFrame={} tensor={} radial={} axial={} cells={} work={}",
                receipt.raw_mouth_plane_distance_mm(),
                receipt.whole_head_frame_transfer_mm(),
                receipt
                    .tensor_transfers()
                    .iter()
                    .map(|t| t.distance_mm)
                    .fold(0., f64::max),
                head.minimum_radial_clearance_mm(),
                head.axial_floor_clearance_mm(),
                head.cells(),
                receipt.work_used()
            );
        }
    }
}
#[test]
fn common_chart_budgets_foreign_receipts_and_all_query_cancellation_phases_refuse_atomically() {
    let text = super::super::super::branches::representative::scaled_text(1., true);
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
    macro_rules! run {
        ($options:expr,$interrupt:expr) => {
            closed.common_mouth_head_envelope(&outer, $options, $interrupt)
        };
    }
    let good = run!(Default::default(), None).unwrap();
    for options in [
        CommonMouthOptions {
            maximum_carrier_distance_mm: f64::NAN,
            ..Default::default()
        },
        CommonMouthOptions {
            maximum_carrier_distance_mm: 0.,
            ..Default::default()
        },
        CommonMouthOptions {
            maximum_carrier_distance_mm: 1e-6,
            ..Default::default()
        },
        CommonMouthOptions {
            minimum_head_clearance_mm: f64::INFINITY,
            ..Default::default()
        },
        CommonMouthOptions {
            minimum_head_clearance_mm: 0.,
            ..Default::default()
        },
        CommonMouthOptions {
            maximum_cells: 0,
            ..Default::default()
        },
        CommonMouthOptions {
            maximum_cells: CELLS + 1,
            ..Default::default()
        },
        CommonMouthOptions {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        CommonMouthOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        CommonMouthOptions {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            run!(options, None),
            Err(PhillipsReadError::Invalid)
        ));
    }
    for options in [
        CommonMouthOptions {
            maximum_work: 0,
            ..Default::default()
        },
        CommonMouthOptions {
            maximum_work: closed.work_used() + outer.work_used(),
            ..Default::default()
        },
        CommonMouthOptions {
            maximum_work: good.work_used() - 1,
            ..Default::default()
        },
        CommonMouthOptions {
            maximum_bytes: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(run!(options, None), Err(PhillipsReadError::Limit)));
    }
    assert!(
        run!(
            CommonMouthOptions {
                maximum_work: good.work_used(),
                ..Default::default()
            },
            None
        )
        .is_ok()
    );
    assert!(
        run!(
            CommonMouthOptions {
                maximum_carrier_distance_mm: good.raw_mouth_plane_distance_mm() / 2.,
                ..Default::default()
            },
            None
        )
        .is_err()
    );
    let initial = good
        .raw_mouth_plane_distance_mm()
        .max(good.whole_head_frame_transfer_mm());
    let tensor = good
        .tensor_transfers()
        .iter()
        .map(|t| t.distance_mm)
        .fold(0., f64::max);
    assert!(tensor > initial);
    assert!(
        run!(
            CommonMouthOptions {
                maximum_carrier_distance_mm: (initial + tensor) / 2.,
                ..Default::default()
            },
            None
        )
        .is_err()
    );
    assert!(
        run!(
            CommonMouthOptions {
                minimum_head_clearance_mm: 1.,
                ..Default::default()
            },
            None
        )
        .is_err()
    );
    // Locate the real lowerable aggregate byte frontier, without reproducing
    // implementation accounting in the test. The exact next-lower budget fails.
    let (mut lo, mut hi) = (1usize, BYTES);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        match run!(
            CommonMouthOptions {
                maximum_bytes: mid,
                ..Default::default()
            },
            None
        ) {
            Ok(_) => hi = mid,
            Err(PhillipsReadError::Limit) => lo = mid + 1,
            Err(e) => panic!("unexpected byte-query error: {e:?}"),
        }
    }
    assert!(lo > good.retained_budget_bytes() + outer.retained_capacity_bytes());
    assert!(
        run!(
            CommonMouthOptions {
                maximum_bytes: lo,
                ..Default::default()
            },
            None
        )
        .is_ok()
    );
    assert!(matches!(
        run!(
            CommonMouthOptions {
                maximum_bytes: lo - 1,
                ..Default::default()
            },
            None
        ),
        Err(PhillipsReadError::Limit)
    ));
    // Distinct parser ownership must never borrow the first arena's earned proof.
    let other = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let foreign = other
        .outer_carrier_bounds(
            OuterCarrierOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    assert!(matches!(
        closed.common_mouth_head_envelope(&foreign, Default::default(), None),
        Err(PhillipsReadError::Invalid)
    ));
    // Find the complete actual path's poll count, then cancel at every publication,
    // factor, original-pole and inherited-head checkpoint. No failed receipt caches.
    let calls = Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    run!(Default::default(), Some(&count)).unwrap();
    let polls = calls.get();
    assert!(polls > 100);
    // Stage-spanning and all final polls; full prefixes are cheap on this tiny packet.
    for phase in (1..=polls)
        .step_by(97)
        .chain(polls.saturating_sub(4)..=polls)
    {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == phase
        };
        assert!(
            matches!(
                run!(Default::default(), Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ),
            "phase {phase}/{polls}"
        );
    }
    let retry = run!(Default::default(), None).unwrap();
    assert_eq!(good.chart(), retry.chart());
    assert_eq!(good.work_used(), retry.work_used());
    assert_eq!(
        good.raw_mouth_plane_distance_mm(),
        retry.raw_mouth_plane_distance_mm()
    );
}
mod closed_head_probe;
mod head_transfer_probe;
