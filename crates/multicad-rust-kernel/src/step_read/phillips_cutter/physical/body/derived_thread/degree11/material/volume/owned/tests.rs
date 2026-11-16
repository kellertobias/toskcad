use super::super::super::super::tests::actual_m3;
use super::*;
use crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head;
fn hash(s: &str) -> u64 {
    s.bytes().fold(14695981039346656037u64, |h, b| {
        (h ^ b as u64).wrapping_mul(1099511628211)
    })
}
#[test]
fn actual_owned_d11_sources_remain_identical_after_all_original_owners_drop() {
    for fitted in [false, true] {
        let mut right_volume: Option<[f64; 2]> = None;
        for left in [false, true] {
            let text = actual_m3(fitted, left);
            let rawhash = hash(&text);
            let rawbytes = text.as_bytes().to_vec();
            let mut detached = None;
            let mut expected = None;
            with_head(&text, &mut |head| {
                let occupied = head
                    .degree11_thread_carrier(Default::default(), None)
                    .unwrap()
                    .certify_occupied(Default::default(), None)
                    .unwrap();
                let material = head
                    .closed_degree11_material(&occupied, Default::default(), None)
                    .unwrap();
                let volume = material.volume_bounds(Default::default(), None).unwrap();
                let owned = volume.owned_snapshot(Default::default(), None).unwrap();
                assert!(
                    owned.accounted_peak_bytes() <= MAX_BYTES && owned.work_used() <= WORK_LIMIT
                );
                assert_eq!(
                    owned.accounted_peak_bytes(),
                    volume.owned_peak_budget_bytes()
                );
                assert!(owned.work_used() > volume.work_used());
                for (a, b) in occupied.carrier().faces().iter().zip(owned.thread_faces()) {
                    assert_eq!(a.source().poles().as_ptr(), b.source().poles().as_ptr());
                    assert_eq!(a.source(), b.source());
                    assert_eq!(a.original_control_entities(), b.original_control_entities());
                    assert_eq!(a.forward(), b.forward());
                    let a = a.transfer();
                    let b = b.transfer();
                    assert_eq!(
                        (
                            a.face_entity,
                            a.surface_entity,
                            a.source_face_index,
                            a.original_v_start_cell,
                            a.original_rail_entities,
                            a.distance_mm,
                            a.derivative_mm,
                            a.original_cells,
                            a.source_cells
                        ),
                        (
                            b.face_entity,
                            b.surface_entity,
                            b.source_face_index,
                            b.original_v_start_cell,
                            b.original_rail_entities,
                            b.distance_mm,
                            b.derivative_mm,
                            b.original_cells,
                            b.source_cells
                        )
                    );
                }
                let view = owned.head().view();
                assert_eq!(view.profile_proposal(), head.source().profile_proposal());
                for i in 0..view.segment_count() {
                    let a = view.segment_certificate(i).unwrap();
                    let b = head.source().segment_certificate(i).unwrap();
                    assert_eq!(
                        (
                            a.role,
                            a.stations,
                            a.progress,
                            a.minimum_progress_mm,
                            a.minimum_radius_mm,
                            a.support_transfer_mm,
                            a.finite_trim_transfer_mm
                        ),
                        (
                            b.role,
                            b.stations,
                            b.progress,
                            b.minimum_progress_mm,
                            b.minimum_radius_mm,
                            b.support_transfer_mm,
                            b.finite_trim_transfer_mm
                        )
                    );
                    let a = view.certified_interval_bounds(i, [0., 1.]);
                    let b = head.source().certified_interval_bounds(i, [0., 1.]);
                    assert_eq!(
                        (a.parameter, a.position_mm, a.derivative_mm),
                        (b.parameter, b.position_mm, b.derivative_mm)
                    );
                }
                let mut enclosure = occupied.local_enclosure_bounds_mm();
                for i in 0..head.source().segment_count() {
                    let p = head
                        .source()
                        .certified_interval_bounds(i, [0., 1.])
                        .position_mm;
                    for axis in enclosure.iter_mut().take(2) {
                        axis[0] = axis[0].min(-p[0][1]);
                        axis[1] = axis[1].max(p[0][1]);
                    }
                    enclosure[2][0] = enclosure[2][0].min(p[1][0]);
                    enclosure[2][1] = enclosure[2][1].max(p[1][1]);
                }
                assert_eq!(owned.local_enclosure_bounds_mm(), enclosure);
                assert_eq!(
                    owned.transfer_collar().local_total_mm,
                    material.transfer_collar().local_total_mm
                );
                assert_eq!(
                    owned.join_radial_margin_mm(),
                    material.join_radial_margin_mm()
                );
                assert_eq!(
                    owned.containment_clearances_mm(),
                    [
                        material.containment().minimum_radial_clearance_mm(),
                        material.containment().axial_floor_clearance_mm()
                    ]
                );
                assert_eq!(owned.chart(), material.chart());
                assert_eq!(owned.determinant_bounds(), material.determinant_bounds());
                assert_eq!(owned.cap_stations_mm(), occupied.cap_stations_mm());
                assert_eq!(
                    owned.occupied_core_radius_mm(),
                    occupied.occupied_core_radius_mm()
                );
                assert_eq!(
                    owned.cap_winding_angle_bounds_rad(),
                    occupied.cap_winding_angle_bounds_rad()
                );
                assert_eq!(owned.cavity(), material.cavity());
                assert_eq!(owned.local_volume_mm3(), volume.local_volume_mm3());
                assert_eq!(owned.world_volume_mm3(), volume.world_volume_mm3());
                assert_eq!(
                    owned.component_volumes_mm3(),
                    volume.component_volumes_mm3()
                );
                assert_eq!(owned.face_count(), material.face_count());
                for i in 0..material.face_count() {
                    assert_eq!(owned.face_source(i), material.face_source(i));
                }
                assert!(owned.face_source(owned.face_count()).is_none());
                assert_eq!(
                    owned.join_inner_uses(),
                    material.join_inner_uses().collect::<Vec<_>>()
                );
                assert!(owned.entry_rim_poles(owned.thread_faces().len()).is_none());
                let points = owned
                    .thread_faces()
                    .iter()
                    .map(|f| f.source().poles().to_vec())
                    .collect::<Vec<_>>();
                let ids = owned
                    .thread_faces()
                    .iter()
                    .map(|f| f.original_control_entities().to_vec())
                    .collect::<Vec<_>>();
                let inner = material.join_inner_uses().collect::<Vec<_>>();
                let rim = (0..owned.thread_faces().len())
                    .map(|i| material.entry_rim_poles(i).unwrap().collect::<Vec<_>>())
                    .collect::<Vec<_>>();
                expected = Some((
                    points,
                    ids,
                    inner,
                    rim,
                    material.chart(),
                    volume.local_volume_mm3(),
                    material.cavity().spec(),
                    material.cavity().domain(),
                ));
                eprintln!(
                    "D11_OWNED fitted={fitted} left={left} retained={} peak={} work={}",
                    owned.retained_budget_bytes(),
                    owned.accounted_peak_bytes(),
                    owned.work_used()
                );
                detached = Some(owned);
            });
            // Body/parser, complete head/cavity receipts, occupied/material/
            // volume views AND the caller's original raw String are all gone.
            drop(text);
            let owned = detached.unwrap();
            let (points, ids, inner, rim, chart, volume, cavity_spec, cavity_domain) =
                expected.unwrap();
            assert_eq!(hash(owned.original_input()), rawhash);
            assert_eq!(owned.original_input().as_bytes(), rawbytes);
            // These evaluations use only the detached exact head and cavity.
            let head = owned.head().view();
            for i in 0..head.segment_count() {
                let bounds = head.certified_interval_bounds(i, [0., 1.]);
                assert!(
                    bounds
                        .position_mm
                        .iter()
                        .chain(bounds.derivative_mm.iter())
                        .flatten()
                        .all(|v| v.is_finite())
                );
            }
            assert!(
                owned
                    .cavity()
                    .envelope()
                    .tip_depth_mm_bounds()
                    .iter()
                    .all(|v| v.is_finite())
            );
            assert_eq!(owned.chart(), chart);
            assert_eq!(owned.local_volume_mm3(), volume);
            assert_eq!(owned.cavity().spec(), cavity_spec);
            assert_eq!(owned.cavity().domain(), cavity_domain);
            assert_eq!(owned.join_inner_uses(), inner);
            for (i, f) in owned.thread_faces().iter().enumerate() {
                assert_eq!(f.source().poles(), points[i]);
                assert_eq!(f.original_control_entities(), ids[i]);
                assert_eq!(
                    owned.entry_rim_poles(i).unwrap().collect::<Vec<_>>(),
                    rim[i]
                );
                assert!(
                    f.source()
                        .cell([0, 0])
                        .unwrap()
                        .first_derivative_bounds()
                        .unwrap()
                        .du
                        .iter()
                        .flatten()
                        .all(|x| x.is_finite())
                );
            }
            assert!(owned.retained_budget_bytes() < owned.accounted_peak_bytes());
            // The actual right/left/posed integrals must overlap; this checks
            // physical hand/pose invariance without choosing a nominal mass.
            if let Some(right) = right_volume {
                assert!(volume[0] <= right[1] && right[0] <= volume[1]);
            } else {
                right_volume = Some(volume);
            }
        }
    }
}
#[test]
fn actual_owned_snapshot_limits_every_cancel_boundary_and_same_receipt_retry() {
    use std::cell::Cell;
    let text = actual_m3(true, false);
    with_head(&text, &mut |head| {
        let occupied = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let material = head
            .closed_degree11_material(&occupied, Default::default(), None)
            .unwrap();
        let volume = material.volume_bounds(Default::default(), None).unwrap();
        let options = Degree11OwnedMaterialOptions::default();
        for bad in [
            Degree11OwnedMaterialOptions {
                maximum_work: WORK_LIMIT + 1,
                ..options
            },
            Degree11OwnedMaterialOptions {
                maximum_bytes: MAX_BYTES + 1,
                ..options
            },
        ] {
            assert!(matches!(
                volume.owned_snapshot(bad, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        let count = Cell::new(0usize);
        let good = volume
            .owned_snapshot(
                options,
                Some(&|| {
                    count.set(count.get() + 1);
                    false
                }),
            )
            .unwrap();
        // Publish real ready geometry against its actual measured capacities,
        // independent of the constructor's deliberately larger reservation.
        let prefix = volume.source().retained_budget_bytes()
            + std::mem::size_of_val(&volume)
            + std::mem::size_of::<OwnedDegree11Material>()
            + SCRATCH;
        let actual_bytes = prefix + good.retained_budget_bytes();
        for exact in [false, true] {
            let ready = volume.owned_snapshot(options, None).unwrap();
            let prior_work = ready.work_used();
            let mut work = Work {
                remaining: options.maximum_work - prior_work,
                interrupt: None,
                cancelled: Cell::new(false),
            };
            let result = publish(
                ready,
                prefix,
                Degree11OwnedMaterialOptions {
                    maximum_bytes: actual_bytes - usize::from(!exact),
                    ..options
                },
                &mut work,
            );
            if exact {
                let value = result.unwrap();
                assert_eq!(value.work_used(), prior_work + 1);
                assert_eq!(value.original_input(), good.original_input());
            } else {
                assert!(matches!(result, Err(PhillipsReadError::Limit)));
            }
        }
        for stop in 1..=count.get() {
            let calls = Cell::new(0usize);
            assert!(
                matches!(
                    volume.owned_snapshot(
                        options,
                        Some(&|| {
                            calls.set(calls.get() + 1);
                            calls.get() >= stop
                        })
                    ),
                    Err(PhillipsReadError::Cancelled)
                ),
                "callback {stop}"
            );
            assert_eq!(calls.get(), stop);
        }
        for bad in [
            Degree11OwnedMaterialOptions {
                maximum_work: volume.work_used() - 1,
                ..options
            },
            Degree11OwnedMaterialOptions {
                maximum_work: good.work_used() - 1,
                ..options
            },
            Degree11OwnedMaterialOptions {
                maximum_bytes: volume.owned_peak_budget_bytes() - 1,
                ..options
            },
        ] {
            assert!(matches!(
                volume.owned_snapshot(bad, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        let retry = volume
            .owned_snapshot(
                Degree11OwnedMaterialOptions {
                    maximum_work: good.work_used(),
                    maximum_bytes: volume.owned_peak_budget_bytes(),
                },
                None,
            )
            .unwrap();
        assert_eq!(retry.work_used(), good.work_used());
        assert_eq!(retry.original_input(), good.original_input());
        assert_eq!(retry.local_volume_mm3(), good.local_volume_mm3());
        for (a, b) in retry.thread_faces().iter().zip(good.thread_faces()) {
            assert_eq!(a.source(), b.source());
            assert_eq!(a.original_control_entities(), b.original_control_entities());
        }
    });
}
#[test]
fn owned_numeric_allocation_preflight_refuses_capacity_and_arithmetic_overflow() {
    assert!(reserve::<u32>(usize::MAX).is_err());
    assert_eq!(reserve::<u32>(4).unwrap().len(), 0);
    assert!(check_peak(usize::MAX, 1, MAX_BYTES).is_err());
    assert!(check_peak(4, 5, 8).is_err());
    assert_eq!(check_peak(4, 5, 9).unwrap(), 9);
}
