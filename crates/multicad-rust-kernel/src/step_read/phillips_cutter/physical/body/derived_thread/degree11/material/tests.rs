use super::super::tests::{actual_m3, m3_source};
use super::*;
use crate::step_read::phillips_cutter::physical::body::common_mouth::tests::{
    closed_inputs, derived_thread_probe::with_head,
};

#[test]
fn actual_flat_and_fitted_m3_same_sources_earn_complete_closed_material_both_hands_and_pose() {
    for fitted in [false, true] {
        for left in [false, true] {
            let text = actual_m3(fitted, left);
            with_head(&text, &mut |head| {
                let occupied = head
                    .degree11_thread_carrier(Default::default(), None)
                    .unwrap()
                    .certify_occupied(Default::default(), None)
                    .unwrap();
                let original: Vec<_> = occupied
                    .carrier()
                    .faces()
                    .iter()
                    .map(|f| f.source().clone())
                    .collect();
                let material = head
                    .closed_degree11_material(&occupied, Default::default(), None)
                    .unwrap();
                eprintln!(
                    "D11_MATERIAL fitted={fitted} left={left} inherited={} new={} total={} peak={} retained={} rimcells={} containmentcells={} margin={} radius={} collar={:?}",
                    occupied.total_work_used(),
                    material.new_work_used(),
                    material.work_used(),
                    material.accounted_peak_bytes(),
                    material.retained_budget_bytes(),
                    material.rim_cells(),
                    material.containment().cells(),
                    material.join_radial_margin_mm(),
                    material.entry_radius_upper_mm(),
                    material.transfer_collar()
                );
                assert!(std::ptr::eq(material.body(), head.body()));
                assert!(std::ptr::eq(material.thread(), &occupied));
                assert!(std::ptr::eq(material.head(), head.source()));
                assert!(std::ptr::eq(
                    material.cavity(),
                    head.common_mouth()
                        .closed_cavity()
                        .local_source_certificate()
                        .source()
                ));
                assert!(std::ptr::eq(
                    material.cavity_certificate().source(),
                    material.cavity()
                ));
                assert_eq!(material.chart(), head.chart());
                assert_eq!(
                    material.determinant_bounds(),
                    head.common_mouth().closed_cavity().determinant_bounds()
                );
                assert_eq!(material.body().literal_input(), text);
                assert!(
                    material.work_used() <= MAX_WORK
                        && material.accounted_peak_bytes() <= MAX_BYTES
                );
                assert!(material.retained_budget_bytes() > occupied.retained_budget_bytes());
                assert!(material.join_radial_margin_mm() > 1e-7);
                let c = material.transfer_collar();
                assert_eq!(
                    c.inherited_sewing_world_mm,
                    head.common_mouth()
                        .closed_cavity()
                        .representative()
                        .inherited_sewing_collar_mm()
                );
                assert!(
                    c.local_total_mm >= c.head_intrinsic_local_mm
                        && c.inverse_component_to_local_l1 >= 3.
                );
                assert!(
                    material.containment().minimum_radial_clearance_mm() > 1e-7 + c.local_total_mm
                );
                assert!(
                    material.containment().axial_floor_clearance_mm() > 1e-7 + c.local_total_mm
                );
                let mut reversed = Vec::new();
                let mut sides = 0;
                let mut headfaces = 0;
                let mut caps = [false; 3];
                for i in 0..material.face_count() {
                    match material.face_source(i).unwrap() {
                        Degree11MaterialFace::HeadMeridian {
                            segment,
                            raw_material_face,
                        } => {
                            assert_eq!(head.body().profile_faces[segment], raw_material_face);
                            headfaces += 1;
                        }
                        Degree11MaterialFace::ThreadSide {
                            index,
                            raw_material_face,
                        } => {
                            assert_eq!(
                                occupied.carrier().faces()[index].transfer().face_entity,
                                raw_material_face
                            );
                            sides += 1;
                        }
                        Degree11MaterialFace::ThreadTip { raw_material_face } => {
                            assert_eq!(
                                raw_material_face,
                                head.body().material_faces[head.body().tensor_face_count + 1]
                            );
                            caps[0] = true;
                        }
                        Degree11MaterialFace::JoinAnnulus { raw_shoulder_face } => {
                            assert_eq!(
                                raw_shoulder_face,
                                head.body().material_faces[head.body().tensor_face_count]
                            );
                            caps[1] = true;
                        }
                        Degree11MaterialFace::TopAnnulus { raw_material_face } => {
                            assert_eq!(raw_material_face, head.body().cap.material_face);
                            caps[2] = true;
                        }
                        Degree11MaterialFace::ReversedCavity {
                            face,
                            raw_material_face,
                        } => {
                            assert_ne!(face, FaceKey::Mouth);
                            assert!(!reversed.contains(&face));
                            assert!(
                                head.common_mouth()
                                    .closed_cavity()
                                    .representative()
                                    .faces()
                                    .iter()
                                    .any(|f| f.source_role == face
                                        && f.raw_material_entity == Some(raw_material_face))
                            );
                            reversed.push(face);
                        }
                    }
                }
                assert_eq!(headfaces, head.source().segment_count());
                assert_eq!(sides, occupied.carrier().faces().len());
                assert_eq!(caps, [true; 3]);
                assert_eq!(reversed.len(), 28);
                assert!(material.face_source(material.face_count()).is_none());
                let uses: Vec<_> = material.join_inner_uses().collect();
                let mut ends = Vec::new();
                for u in uses {
                    let mut p: Vec<_> = material.entry_rim_poles(u.face_index).unwrap().collect();
                    if u.reversed {
                        p.reverse();
                    }
                    ends.push((p[0], *p.last().unwrap()));
                }
                for i in 0..ends.len() {
                    assert_eq!(ends[i].1, ends[(i + 1) % ends.len()].0);
                }
                assert!(material.entry_rim_poles(sides).is_none());
                for (f, old) in occupied.carrier().faces().iter().zip(&original) {
                    assert_eq!(f.source(), old);
                }
            });
        }
    }
}

#[test]
fn strict_original_rim_clearance_refuses_shared_crest_and_refines_real_narrow_annulus() {
    for (shank, success) in [(1.47, false), (1.471, true)] {
        // Independent flat-root fixtures; no canonical6g dimensions are edited.
        let text = m3_source(false, false, shank, 1., false);
        with_head(&text, &mut |head| {
            let occupied = head
                .degree11_thread_carrier(Default::default(), None)
                .unwrap()
                .certify_occupied(Default::default(), None)
                .unwrap();
            let result = head.closed_degree11_material(&occupied, Default::default(), None);
            if success {
                let material = result.unwrap();
                assert!(material.join_radial_margin_mm() > 1e-7);
                assert!(
                    material.rim_cells()
                        > occupied
                            .carrier()
                            .faces()
                            .iter()
                            .map(|f| f.source().spans()[0])
                            .sum::<usize>()
                );
                assert!(matches!(
                    head.closed_degree11_material(
                        &occupied,
                        Degree11MaterialOptions {
                            maximum_rim_depth: 0,
                            ..Default::default()
                        },
                        None
                    ),
                    Err(PhillipsReadError::Unsupported)
                ));
            } else {
                assert!(matches!(result, Err(PhillipsReadError::Unsupported)));
            }
        });
    }
}

#[test]
fn actual_foreign_body_and_independent_binding_metadata_are_rejected() {
    let text = actual_m3(true, false);
    with_head(&text, &mut |head| {
        let occupied = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        with_head(&text, &mut |foreign| {
            assert!(matches!(
                foreign.closed_degree11_material(&occupied, Default::default(), None),
                Err(PhillipsReadError::Unsupported)
            ))
        });
        let rebuilt = head
            .common_mouth()
            .closed_repaired_head(Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(rebuilt.body(), head.body()));
        assert!(std::ptr::eq(
            rebuilt.source().profile_proposal(),
            head.source().profile_proposal()
        ));
        assert!(!std::ptr::eq(rebuilt.source(), head.source()));
        assert!(matches!(
            rebuilt.closed_degree11_material(&occupied, Default::default(), None),
            Err(PhillipsReadError::Unsupported)
        ));
        assert!(matches!(
            rebuilt.degree11_thread_contacts(occupied.carrier(), Default::default(), None),
            Err(PhillipsReadError::Invalid)
        ));
        let body = head.body();
        let profile = head.source().profile_proposal();
        let chart = head.chart();
        let z = head.virtual_bottom_station()[1];
        assert!(compatible(head, body, profile, chart, z, &occupied));
        // Scalar provenance validator inputs, never fabricated occupied types.
        let mut wrong = chart;
        wrong.columns[0][0] += 1.;
        assert!(!compatible(head, body, profile, wrong, z, &occupied));
        assert!(!compatible(head, body, profile, chart, z + 1., &occupied));
        let other_profile = *profile;
        assert!(!compatible(head, body, &other_profile, chart, z, &occupied));
    });
}

#[test]
fn actual_core_entering_counter_cavity_has_no_head_only_material_receipt() {
    let text = m3_source(false, false, 1.5, 0., true);
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
    let expected = sewn
        .edges()
        .iter()
        .map(|e| {
            let r = e.raw_owner_distances_mm();
            (I::point(e.position_transfer_mm()) + I::point(r[0].max(r[1]))).hi
        })
        .fold(0f64, f64::max);
    assert_eq!(rep.inherited_sewing_collar_mm(), expected);
    assert!(expected > 0. && expected < 1e-7);
    let closed = rep
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    assert!(matches!(
        closed.common_mouth_head_envelope(&outer, Default::default(), None),
        Err(PhillipsReadError::Unsupported)
    ));
    assert_eq!(body.literal_input(), text);
}

#[test]
fn actual_inverse_chart_bound_is_directed_and_not_a_rigid_projection() {
    let chart = CommonMouthChart {
        columns: [[0.; 3], [1., 0., 0.], [0.1, 1., 0.], [0., 0., 1.]],
        outer_top_z: 0.,
    };
    let bound = inverse_bound(chart).unwrap();
    // Exact inverse rows(1,-.1,0),(0,1,0),(0,0,1): entrywise sum3.1.
    assert!((3.1..3.10000000001).contains(&bound));
    let mut bad = chart;
    bad.columns[2] = bad.columns[1];
    assert!(inverse_bound(bad).is_err());
}
#[test]
fn transfer_collar_composition_checks_real_numeric_bound_domain_and_overflow() {
    let b = [0.1, 0.2, 0.3, 0.4, 0.5];
    let sum = compose_collar(b, 3.).unwrap();
    assert!((3.5..3.50000000001).contains(&sum));
    for index in 0..5 {
        for invalid in [-1., f64::NAN, f64::INFINITY] {
            let mut wrong = b;
            wrong[index] = invalid;
            assert!(compose_collar(wrong, 3.).is_err());
        }
    }
    for invalid in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(compose_collar(b, invalid).is_err());
    }
    assert!(compose_collar([f64::MAX, 0., 0., 0., 0.], 3.).is_err());
}

#[test]
fn closed_d11_material_resource_cancellation_and_clean_retry_contract() {
    let text = actual_m3(true, false);
    with_head(&text, &mut |head| {
        let occupied = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let ticks = Cell::new(0usize);
        let observe = || {
            ticks.set(ticks.get() + 1);
            false
        };
        let reference = head
            .closed_degree11_material(&occupied, Default::default(), Some(&observe))
            .unwrap();
        let calls = ticks.get();
        let o = Degree11MaterialOptions::default();
        for limit in [
            Degree11MaterialOptions {
                maximum_work: reference.work_used() - 1,
                ..o
            },
            Degree11MaterialOptions {
                maximum_work: occupied.total_work_used() - 1,
                ..o
            },
            Degree11MaterialOptions {
                maximum_bytes: reference.accounted_peak_bytes() - 1,
                ..o
            },
            Degree11MaterialOptions {
                maximum_rim_cells: reference.rim_cells() - 1,
                ..o
            },
            Degree11MaterialOptions {
                maximum_containment_cells: reference.containment().cells() - 1,
                ..o
            },
        ] {
            assert!(matches!(
                head.closed_degree11_material(&occupied, limit, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        // This fixed fixture has a small metadata/containment frontier. Sweep
        // every checkpoint, including the nested cavity provider, so error
        // attribution and cancellation remain observable across the handoff.
        assert!(calls < 1024);
        for at in 1..=calls {
            let n = Cell::new(0);
            let stop = || {
                n.set(n.get() + 1);
                n.get() == at
            };
            assert!(matches!(
                head.closed_degree11_material(&occupied, o, Some(&stop)),
                Err(PhillipsReadError::Cancelled)
            ));
        }
        assert!(matches!(
            head.closed_degree11_material(
                &occupied,
                Degree11MaterialOptions {
                    minimum_clearance_mm: f64::NAN,
                    ..o
                },
                None
            ),
            Err(PhillipsReadError::Invalid)
        ));
        let retry = head
            .closed_degree11_material(
                &occupied,
                Degree11MaterialOptions {
                    maximum_work: reference.work_used(),
                    maximum_bytes: reference.accounted_peak_bytes(),
                    ..o
                },
                None,
            )
            .unwrap();
        assert_eq!(
            retry.join_radial_margin_mm(),
            reference.join_radial_margin_mm()
        );
        assert_eq!(retry.work_used(), reference.work_used());
        assert_eq!(
            retry.entry_radius_upper_mm(),
            reference.entry_radius_upper_mm()
        );
        assert_eq!(retry.body().literal_input(), text);
    });
}

#[test]
fn requested_cavity_wall_margin_cannot_be_borrowed_from_a_wider_thread_join() {
    // An independent wide-shank flat-root source, not an altered6g recipe.
    let text = m3_source(false, false, 2., 1., false);
    with_head(&text, &mut |head| {
        let occupied = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let material = head
            .closed_degree11_material(&occupied, Default::default(), None)
            .unwrap();
        let entry = material.cavity().envelope().entry_radius_mm_bounds();
        let upper_mouth_gap =
            (I::point(material.head().station(0).unwrap()[0]) - I::point(entry[0])).hi;
        assert!(material.join_radial_margin_mm() > upper_mouth_gap + 0.01);
        let request = (upper_mouth_gap + material.join_radial_margin_mm()) / 2.;
        assert!(matches!(
            head.closed_degree11_material(
                &occupied,
                Degree11MaterialOptions {
                    minimum_clearance_mm: request,
                    ..Default::default()
                },
                None
            ),
            Err(PhillipsReadError::Unsupported)
        ));
    });
}
