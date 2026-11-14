use super::*;
#[test]
fn actual_metric_outer_tensors_and_finite_trims_match_without_a_recess_recipe() {
    for index in 0..6 {
        for left in [false, true] {
            let source = super::super::super::tests::fixture(index, left, false, 1);
            let (o, u, v) = if left {
                ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
            } else {
                ([0.; 3], [1., 0., 0.], [0., 1., 0.])
            };
            let text = crate::step::write_phillips_screw(&source, o, u, v, 1e-7).unwrap();
            let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
            let receipt = body
                .outer_carrier_bounds(
                    OuterCarrierOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 1_000_000,
                    },
                    None,
                )
                .unwrap_or_else(|e| panic!("{index} left{left}: {e:?}"));
            assert!(std::ptr::eq(receipt.body(), &body));
            assert_eq!(receipt.tensors().len(), 4);
            assert!(receipt.work_used() > 0);
            assert!(receipt.meridian_reconstruction_error_mm() > 0.);
            assert!(receipt.retained_capacity_bytes() > std::mem::size_of_val(&receipt));
            for tensor in receipt.tensors() {
                assert!(tensor.control_rows >= 12 && tensor.control_columns >= 12);
                assert!(tensor.tip_plane_distance_mm > 0. && tensor.tip_plane_distance_mm < 1e-7);
                if let Some(error) = tensor.entry_circle_distance_mm {
                    assert!(error > 0. && error < 1e-7);
                }
                assert!(tensor.control_deviation_mm > 0. && tensor.control_deviation_mm < 1e-7);
                eprintln!(
                    "H{index} left{left} tensor{} residual{}",
                    tensor.face_entity, tensor.control_deviation_mm
                );
            }
        }
    }
}
fn fixture() -> String {
    let source = super::super::super::tests::fixture(1, false, false, 1);
    crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap()
}
fn options(error: f64) -> OuterCarrierOptions {
    OuterCarrierOptions {
        maximum_distance_mm: error,
        work_limit: 1_000_000,
    }
}
#[test]
fn actual_interior_pole_change_is_retained_and_rejected_or_measured_never_replaced() {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, 1e-7, true);
    let (surface, _, _) = helical::face(&shell, body.tensor_face_entities()[0]).unwrap();
    let Value::List(rows) = &shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS").unwrap()[3] else {
        panic!("rows")
    };
    let point = refs(&rows[1]).unwrap()[1];
    let p = shell.point(point).unwrap();
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{point}=")))
        .unwrap();
    let changed = text.replacen(
        old,
        &format!(
            "#{point}=CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}));",
            p[0] + 0.001,
            p[1],
            p[2]
        ),
        1,
    );
    let changed = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    assert_eq!(changed.inferred_thread_spec(), body.inferred_thread_spec());
    assert!(changed.outer_carrier_bounds(options(1e-7), None).is_err());
    let receipt = changed.outer_carrier_bounds(options(0.01), None).unwrap();
    let actual = (p[0] + 0.001) - p[0];
    assert!(receipt.tensors()[0].control_deviation_mm >= actual);
    assert!(receipt.tensors()[0].control_deviation_mm < 0.001000001);
    assert_eq!(
        receipt.body().literal_tensor_control_points(0).unwrap()[1][1],
        [p[0] + 0.001, p[1], p[2]]
    );
    // A tighter numerical request refuses rather than assigning its value as a
    // measured error or silently snapping retained controls to the header law.
    assert!(body.outer_carrier_bounds(options(1e-20), None).is_err());
}
#[test]
fn malformed_later_nets_and_work_cancellation_refuse_before_a_carrier_receipt() {
    let source = super::super::super::tests::fixture(1, false, true, 2);
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    body.outer_carrier_bounds(options(1e-7), None).unwrap();
    let shell = Shell::new(&body.entities, 1e-7, true);
    let (surface, _, _) = helical::face(&shell, body.tensor_face_entities()[4]).unwrap();
    for value in ["0.", "(0.)"] {
        let changed = super::super::tests::argument(&text, surface, 3, value);
        let parsed = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
        assert!(parsed.outer_carrier_bounds(options(1e-7), None).is_err());
    }
    for error in [0., -1., f64::INFINITY, f64::NAN] {
        assert!(body.outer_carrier_bounds(options(error), None).is_err());
    }
    assert!(
        body.outer_carrier_bounds(
            OuterCarrierOptions {
                work_limit: 0,
                ..options(1e-7)
            },
            None
        )
        .is_err()
    );
    assert!(
        body.outer_carrier_bounds(
            OuterCarrierOptions {
                work_limit: 10,
                ..options(1e-7)
            },
            None
        )
        .is_err()
    );
    assert!(matches!(
        body.outer_carrier_bounds(options(1e-7), Some(&|| true)),
        Err(PhillipsReadError::Cancelled)
    ));
    let calls = Cell::new(0);
    let interrupt = || {
        calls.set(calls.get() + 1);
        calls.get() > 30
    };
    assert!(matches!(
        body.outer_carrier_bounds(options(1e-7), Some(&interrupt)),
        Err(PhillipsReadError::Cancelled)
    ));
    for visible in [f64::NAN, -1., body.thread.spec().length_mm] {
        assert!(prepare_retained_reference(&body.thread, visible, 1e-7).is_err());
    }
    assert!(prepare_retained_reference(&body.thread, 0., 0.).is_err());
    assert!(norm([I::point(f64::INFINITY); 3]).is_err());
}
#[test]
fn true_rounded_root_phase_and_rolled_fully_threaded_datum_are_matched_as_file_carriers() {
    for left in [false, true] {
        let base = super::super::super::tests::fixture(1, left, false, 1);
        let spec = crate::helical_thread::Spec {
            major_radius_mm: 1.4635,
            pitch_mm: 0.5,
            length_mm: 9.,
            depth_mm: 0.29908756729740676,
            crest_width_mm: 0.07131009168581062,
            root_radius_mm: Some(0.07216878364870323),
            ..base.thread().spec()
        };
        let thread = crate::helical_thread::Thread::new(spec).unwrap();
        let source = crate::phillips_drive::PhillipsScrew::new(
            thread,
            *base.profile(),
            base.spec(),
            base.cavity().clone(),
            100_000,
            &|| false,
        )
        .unwrap();
        let text = crate::step::write_phillips_screw(
            &source,
            [105., -20., 60.],
            [0.8, 0.6, 0.],
            [0.36, -0.48, -0.8],
            1e-7,
        )
        .unwrap();
        let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
        let receipt = body.outer_carrier_bounds(options(1e-7), None).unwrap();
        assert!(receipt.body().inferred_visible_interval_mm()[0] > 0.);
        assert!(
            receipt
                .body()
                .inferred_thread_spec()
                .root_radius_mm
                .is_some()
        );
        assert_eq!(receipt.body().inferred_thread_spec().left_hand, left);
        for t in receipt.tensors() {
            assert!(t.control_deviation_mm < 1e-10);
        }
    }
}
fn replace_direction(text: &str, id: u32, p: Point) -> String {
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{id}=")))
        .unwrap();
    text.replacen(
        old,
        &format!(
            "#{id}=DIRECTION('',({:.16E},{:.16E},{:.16E}));",
            p[0], p[1], p[2]
        ),
        1,
    )
}
#[test]
fn structurally_retained_tilted_tip_and_head_supports_do_not_pass_geometric_matching() {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, 1e-7, true);
    for face_id in [
        body.material_faces[body.tensor_face_count + 1],
        body.profile_faces[0],
    ] {
        let (surface, _, _) = helical::face(&shell, face_id).unwrap();
        let axis = reference(&shell.simple(surface).unwrap().args[1]).unwrap();
        let normal = reference(&shell.args(axis, "AXIS2_PLACEMENT_3D").unwrap()[2]).unwrap();
        let p = shell.direction(normal).unwrap();
        let changed = replace_direction(&text, normal, [p[0] + 1e-8, p[1], p[2]]);
        let retained = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
        assert!(retained.outer_carrier_bounds(options(1e-7), None).is_err());
    }
    // The shared primitive is also used by the legacy physical reader, where
    // raw malformed cap schema and wrong finite wire are rejected directly.
    let id = body.material_faces[body.tensor_face_count + 1];
    let (surface, _, bounds) = helical::face(&shell, id).unwrap();
    let basis = prepare_retained_reference(&body.thread, body.visible[0], 1e-7).unwrap();
    let proof = helical::Proof {
        error: Some(1e-7),
        extent: 2.,
    };
    let mut proof_shell = Shell::new(&body.entities, 1e-7, true);
    let prefix = helical::verify_prefix(
        &mut proof_shell,
        body.tensor_face_entities(),
        &basis,
        body.frame,
        helical::PrefixContract {
            proof,
            forward: false,
            legacy_shank_entry: false,
            all_entry_crests: true,
            exit_crests: false,
        },
        |_, _| Ok(()),
    )
    .unwrap();
    let raw = text
        .lines()
        .find(|l| l.starts_with(&format!("#{surface}=")))
        .unwrap();
    let changed = text.replacen(raw, &format!("#{surface}=PLANE('');"), 1);
    let parsed = crate::step_read::part21::parse(&changed).unwrap();
    assert!(
        helical::verify_plain_tip(
            proof,
            &Shell::new(&parsed, 1e-7, true),
            id,
            body.frame,
            (2., false),
            &prefix.rings[1]
        )
        .is_err()
    );
    let mut wrong = prefix.rings[1].clone();
    wrong.swap(0, 1);
    assert!(helical::verify_plain_tip(proof, &shell, id, body.frame, (2., false), &wrong).is_err());
    assert_eq!(bounds.len(), 1);
}
#[test]
fn a_request_below_measured_actual_support_residual_refuses_even_when_the_net_reference_matches() {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let receipt = body.outer_carrier_bounds(options(1e-7), None).unwrap();
    let actual = receipt
        .tensors()
        .iter()
        .map(|t| {
            t.tip_plane_distance_mm
                .max(t.entry_circle_distance_mm.unwrap_or(0.))
        })
        .fold(0_f64, f64::max);
    let source = receipt.meridian_reconstruction_error_mm();
    eprintln!("actual cap support residual {actual}, head reserve {source}");
    assert!(actual > source);
    let request = source.midpoint(actual);
    assert!(
        receipt
            .tensors()
            .iter()
            .all(|t| t.control_deviation_mm < request)
    );
    assert!(body.outer_carrier_bounds(options(request), None).is_err());
}
#[test]
fn approximate_headers_do_not_erase_displaced_join_or_tilted_finite_circle_supports() {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, 1e-7, true);
    let (surface, _, _) =
        helical::face(&shell, body.material_faces[body.tensor_face_count]).unwrap();
    let axis = reference(&shell.simple(surface).unwrap().args[1]).unwrap();
    let origin = reference(&shell.args(axis, "AXIS2_PLACEMENT_3D").unwrap()[1]).unwrap();
    let p = shell.point(origin).unwrap();
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{origin}=")))
        .unwrap();
    let changed = text.replacen(
        old,
        &format!(
            "#{origin}=CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}));",
            p[0] + 1e-10,
            p[1],
            p[2]
        ),
        1,
    );
    let retained = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    assert!(retained.outer_carrier_bounds(options(1e-7), None).is_err());
    let (_, _, bounds) = helical::face(&shell, body.cap.material_face).unwrap();
    let ring = helical::ring(&shell, bounds[0], true).unwrap();
    let curve = reference(&shell.args(ring[0].0, "EDGE_CURVE").unwrap()[3]).unwrap();
    let axis = reference(&shell.args(curve, "CIRCLE").unwrap()[1]).unwrap();
    let normal = reference(&shell.args(axis, "AXIS2_PLACEMENT_3D").unwrap()[2]).unwrap();
    let p = shell.direction(normal).unwrap();
    let changed = replace_direction(&text, normal, [p[0] + 1e-11, p[1], p[2]]);
    let retained = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    assert!(retained.outer_carrier_bounds(options(1e-7), None).is_err());
}
#[test]
fn original_spline_integer_metadata_cannot_be_replaced_by_nearby_reals_with_the_same_float() {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, 1e-7, true);
    let (surface, forward, bounds) = helical::face(&shell, body.tensor_face_entities()[0]).unwrap();
    let wire = helical::orient(&helical::ring(&shell, bounds[0], true).unwrap(), forward);
    let curve = reference(&shell.args(wire[1].0, "EDGE_CURVE").unwrap()[3]).unwrap();
    assert_eq!("11.0000000000000001".parse::<f64>().unwrap(), 11.);
    for (entity, index) in [(surface, 1), (surface, 2), (curve, 1)] {
        let changed = super::super::tests::argument(&text, entity, index, "11.0000000000000001");
        let retained = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
        assert!(retained.outer_carrier_bounds(options(1e-7), None).is_err());
    }
    assert_eq!("12.0000000000000001".parse::<f64>().unwrap(), 12.);
    for (entity, index) in [(surface, 8), (surface, 9), (curve, 6)] {
        let Value::List(values) = &shell.simple(entity).unwrap().args[index] else {
            panic!("multiplicities")
        };
        let fields: Vec<_> = values
            .iter()
            .enumerate()
            .map(|(i, value)| {
                if i == 0 {
                    "12.0000000000000001".to_owned()
                } else {
                    let Value::Number(v) = value else {
                        panic!("integer")
                    };
                    format!("{}", *v as u64)
                }
            })
            .collect();
        let changed =
            super::super::tests::argument(&text, entity, index, &format!("({})", fields.join(",")));
        let retained = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
        assert!(retained.outer_carrier_bounds(options(1e-7), None).is_err());
    }
    let source = super::super::super::tests::fixture(1, false, true, 2);
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let (surface, _, _) = helical::face(
        &Shell::new(&body.entities, 1e-7, true),
        body.tensor_face_entities()[4],
    )
    .unwrap();
    let changed = super::super::tests::argument(&text, surface, 8, "0.");
    let retained = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    assert!(retained.outer_carrier_bounds(options(1e-7), None).is_err());
}
