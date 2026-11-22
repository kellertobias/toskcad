//! Real separate clamp/wood inputs; source-only staging never claims a validated connection.
use multicad_core::{
    commands::{Actor, Channel},
    model::*,
};
use multicad_model::Compiler;
#[path = "support/clamped_insert.rs"]
mod stack_fixture;
use stack_fixture::request;

#[test]
fn actual_two_stock_sources_expose_independent_clock_and_clearance_contracts() {
    let (session, r) = request();
    let pending = prepare_screw_connection(&session, 0, Actor::local(Channel::Mcp), &r).unwrap();
    let roots: Vec<_> = r
        .receivers
        .iter()
        .map(|row| row.source.output.body)
        .chain(r.hardware.iter().map(|row| row.instance))
        .collect();
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_sources(pending.prepared.document(), 1, &roots)
    .unwrap();
    assert!(model.snapshot().screw_connection_checks.is_empty());
    let hardware = |role| {
        let row = r.hardware.iter().find(|row| row.role == role).unwrap();
        model
            .leaf_solid(&OccurrencePath {
                root: row.instance,
                members: vec![],
            })
            .unwrap()
    };
    let insert = hardware(ScrewConnectionHardwareRole::Insert);
    let screw = hardware(ScrewConnectionHardwareRole::Screw);
    let data = insert.rust_wood_insert_data().unwrap();
    let (actual, screw_frame) = screw.rust_machine_screw_geometry().unwrap();
    eprintln!(
        "INSERT frame={:?}; SCREW frame={:?}, source={:?}",
        data.frame,
        screw_frame,
        actual.thread().spec()
    );
    assert_eq!(data.source.full_length_mm(), 4.);
    assert_eq!(actual.thread().spec().length_mm, 8.);
    let wood = model.leaf_solid(&r.receivers[1].source.occurrence).unwrap();
    let policy = wood.geometry_policy();
    let req = multicad_kernel::rust_kernel::wood_thread_contact::Requirements {
        minimum_wall_mm: 1.,
        minimum_floor_mm: 1.,
        minimum_engagement_mm: 2.,
        tip_clearance_mm: 1.,
        linear_tolerance_mm: policy.linear_tolerance_mm,
        minimum_feature_mm: policy.minimum_feature_mm,
    };
    let contact = wood.rust_wood_insert_contact(insert, req, None);
    let mating = screw.rust_wood_insert_mating(insert, 2., None);
    let clearance = wood.rust_wood_insert_screw_clearance(screw, 1., None);
    eprintln!(
        "actual contact={:?}; mating={:?}; screw-floor={:?}",
        contact, mating, clearance
    );
    let frame =
        |f: multicad_kernel::PlaneFrame| multicad_kernel::rust_kernel::wood_thread_contact::Frame {
            origin_mm: f.origin,
            axes: [f.u, f.v, f.normal()],
        };
    let new = multicad_kernel::rust_kernel::wood_thread_contact::certify_clocked_insert_mating(
        &data.source,
        frame(data.frame),
        actual,
        frame(screw_frame),
        multicad_kernel::rust_kernel::wood_thread_contact::ClockedMatingOptions {
            requirements: multicad_kernel::rust_kernel::wood_thread_contact::MatingRequirements {
                minimum_engagement_mm: 2.,
                minimum_clearance_mm: 4. * policy.linear_tolerance_mm,
            },
            maximum_cells: 4096,
            linear_tolerance_mm: policy.linear_tolerance_mm,
        },
        &|| false,
    );
    match &new {
        Ok(v) => eprintln!(
            "actual CLOCKED cells={}, radial={}, axial={}, angle={}, mating={:?}",
            v.cells_used, v.radial_map_error_mm, v.axial_map_error_mm, v.angle_error_rad, v.mating
        ),
        Err(e) => eprintln!("actual CLOCKED refusal={e:?}"),
    };
    assert!(new.is_ok());
    let general_floor = wood.rust_clocked_wood_insert_screw_clearance(screw, 1., None);
    eprintln!("actual CLOCKED floor={general_floor:?}");
    assert!(general_floor.is_ok());
    let clamp = model.leaf_solid(&r.receivers[0].source.occurrence).unwrap();
    let baseline = clamp
        .rust_wood_insert_clamp_contact(screw, insert, 1., None)
        .unwrap();
    eprintln!(
        "MANAGED CLAMP land={}, head-gap={:?}, insert-gap={:?}",
        baseline.bearing_land_lower_mm,
        baseline.bearing_face_gap_mm,
        baseline.insert_to_clamp_gap_mm
    );
    let actual_inset = wood
        .rust_inset_wood_insert_contact(insert, req, policy.linear_tolerance_mm, None)
        .unwrap();
    eprintln!(
        "MANAGED WOOD inset={:?}, engagement={:?}, gap={:?}",
        actual_inset.inset_mm,
        actual_inset.contact.engagement_mm,
        actual_inset.contact.tip_clearance_mm
    );
    assert!(
        actual_inset.inset_mm[0] > 0. && actual_inset.inset_mm[1] <= policy.linear_tolerance_mm
    );
    let original_front = insert
        .transformed(multicad_kernel::RigidTransform {
            translation_mm: [0., 0., -data.frame.origin[2]],
            ..Default::default()
        })
        .unwrap();
    assert!(
        clamp
            .rust_wood_insert_clamp_contact(screw, &original_front, 1., None)
            .is_err()
    );
    assert!(contact.is_err());
    assert!(mating.is_err());
    assert!(clearance.is_err());
}

#[test]
fn managed_actual_stack_receipts_survive_both_hands_history_and_reopen() {
    use multicad_core::project_file::{self, ProjectFile};
    for left in [false, true] {
        let (mut session, mut r) = request();
        for row in &mut r.hardware {
            let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
                unreachable!()
            };
            match &mut part.parameters {
                MechanicalPartParameters::MachineScrew { thread, .. } => {
                    thread.hand = if left {
                        ThreadHand::Left
                    } else {
                        ThreadHand::Right
                    }
                }
                MechanicalPartParameters::CustomConnection { parameters } => {
                    let CustomConnectionHardwareParameters::WoodInsert {
                        outer_thread,
                        inner_thread,
                        ..
                    } = parameters.as_mut()
                    else {
                        unreachable!()
                    };
                    for thread in [outer_thread, inner_thread] {
                        thread.hand = if left {
                            ThreadHand::Left
                        } else {
                            ThreadHand::Right
                        };
                    }
                }
                _ => unreachable!(),
            }
        }
        let before = session.document().clone();
        let prepared =
            prepare_screw_connection(&session, 0, Actor::local(Channel::Mcp), &r).unwrap();
        let compiler = Compiler::default();
        let model = compiler.compile(prepared.prepared.document(), 1).unwrap();
        let check = &model.snapshot().screw_connection_checks[0];
        let gap = check.actual_head_bearing_gap_mm.unwrap();
        let inset = check.actual_insert_inset_mm.unwrap();
        assert!(gap[0] >= 0. && gap[0] <= gap[1] && gap[1] <= 1e-7);
        assert!(inset[0] > 0. && inset[0] <= inset[1] && inset[1] <= 1e-7);
        assert!(
            check
                .scope
                .contains("stock-stock disjointness is not certified")
        );
        assert!(check.native_component_collision_checked);
        let wire = serde_json::to_value(check).unwrap();
        assert_eq!(wire["actual_head_bearing_gap_mm"], serde_json::json!(gap));
        assert_eq!(wire["actual_insert_inset_mm"], serde_json::json!(inset));
        let mut older = wire.clone();
        older
            .as_object_mut()
            .unwrap()
            .remove("actual_head_bearing_gap_mm");
        older
            .as_object_mut()
            .unwrap()
            .remove("actual_insert_inset_mm");
        let compatible: ScrewConnectionCheck = serde_json::from_value(older).unwrap();
        assert_eq!(compatible.actual_head_bearing_gap_mm, None);
        assert_eq!(compatible.actual_insert_inset_mm, None);
        let compatible = serde_json::to_value(compatible).unwrap();
        assert!(compatible.get("actual_head_bearing_gap_mm").is_none());
        assert!(compatible.get("actual_insert_inset_mm").is_none());
        let mut frames = vec![];
        for row in &r.hardware {
            let solid = model
                .leaf_solid(&OccurrencePath {
                    root: row.instance,
                    members: vec![],
                })
                .unwrap();
            if row.role == ScrewConnectionHardwareRole::Screw {
                let (source, frame) = solid.rust_machine_screw_geometry().unwrap();
                assert_eq!(source.thread().spec().length_mm, 8.);
                assert_eq!(source.thread().spec().left_hand, left);
                frames.push(frame);
            } else {
                let data = solid.rust_wood_insert_data().unwrap();
                assert_eq!(data.source.full_length_mm(), 4.);
                frames.push(data.frame);
            }
            let mesh = solid
                .manufacturing_mesh_with_receipt(Default::default())
                .unwrap();
            assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
        }
        for receiver in &r.receivers {
            let solid = model.leaf_solid(&receiver.source.occurrence).unwrap();
            let mesh = solid
                .manufacturing_mesh_with_receipt(Default::default())
                .unwrap();
            assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
            let directory = tempfile::tempdir().unwrap();
            let path = directory.path().join("actual-receiver.step");
            solid.write_step(&path).unwrap();
            let text = std::fs::read_to_string(path).unwrap();
            let decoded = multicad_kernel::Solid::import_rust_step(
                &text,
                multicad_core::document::EntityId::new(),
                Default::default(),
            )
            .unwrap();
            assert_eq!(decoded.solid_count().unwrap(), 1);
            decoded
                .manufacturing_mesh_with_receipt(Default::default())
                .unwrap();
        }
        session.commit(prepared.prepared).unwrap();
        let created = session.document().clone();
        session.undo(1, Actor::local(Channel::Mcp)).unwrap();
        assert_eq!(*session.document(), before);
        session.redo(2, Actor::local(Channel::Mcp)).unwrap();
        assert_eq!(*session.document(), created);
        r.receivers[1].source.extent = ConnectionReceiverExtent::Blind {
            depth: length(6.25),
            floor: length(1.),
        };
        let edit =
            prepare_screw_connection(&session, session.revision(), Actor::local(Channel::Mcp), &r)
                .unwrap();
        compiler
            .compile(edit.prepared.document(), session.revision() + 1)
            .unwrap();
        session.commit(edit.prepared).unwrap();
        let mut bytes = vec![];
        project_file::write(
            &ProjectFile::new(session.document().clone()).unwrap(),
            &mut bytes,
            Default::default(),
        )
        .unwrap();
        let reopened = project_file::read(
            &mut bytes.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap()
        .project
        .document;
        assert_eq!(reopened, *session.document());
        let restored = compiler.compile(&reopened, 0).unwrap();
        assert_eq!(
            restored.snapshot().screw_connection_checks[0].actual_insert_inset_mm,
            Some(inset)
        );
        for (row, expected) in r.hardware.iter().zip(frames) {
            let solid = restored
                .leaf_solid(&OccurrencePath {
                    root: row.instance,
                    members: vec![],
                })
                .unwrap();
            let actual = if row.role == ScrewConnectionHardwareRole::Screw {
                solid.rust_machine_screw_geometry().unwrap().1
            } else {
                solid.rust_wood_insert_data().unwrap().frame
            };
            assert_eq!(actual, expected);
        }
    }
}

#[test]
fn unrelated_actual_material_and_unmatched_receiver_inputs_are_never_waived() {
    for (origin, success) in [([100.; 3], true), ([9., 9., -2.], false)] {
        let (mut session, r) = request();
        let extra = PrimitiveRequest::new(
            "Independent real stock",
            PrimitiveParameters::Box {
                size: std::array::from_fn(|_| length(2.)),
                frame: PrimitiveFrame {
                    origin: origin.map(length),
                    ..Default::default()
                },
            },
        );
        let q = prepare_primitive(&session, 0, Actor::local(Channel::Mcp), &extra).unwrap();
        session.commit(q.prepared).unwrap();
        let before = session.document().clone();
        let q =
            prepare_screw_connection(&session, session.revision(), Actor::local(Channel::Mcp), &r)
                .unwrap();
        let result = Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile(q.prepared.document(), session.revision() + 1);
        assert_eq!(
            result.is_ok(),
            success,
            "{}",
            result.err().map(|e| format!("{e:?}")).unwrap_or_default()
        );
        assert_eq!(*session.document(), before);
    }
    for case in 0..4 {
        let (session, mut r) = request();
        let before = session.document().clone();
        match case {
            0 => {
                r.receivers[1].source.extent = ConnectionReceiverExtent::Blind {
                    depth: length(4.5),
                    floor: length(1.),
                }
            }
            1 => r.receivers[1].source.station = length(4.),
            2 => r.receivers.swap(0, 1),
            3 => {
                r.receivers[0].source.seat = Some(ConnectionReceiverSeat::Counterbore {
                    diameter: length(5.),
                    depth: length(1.),
                })
            }
            _ => unreachable!(),
        }
        if let Ok(q) = prepare_screw_connection(&session, 0, Actor::local(Channel::Mcp), &r) {
            assert!(
                Compiler {
                    omit_display: true,
                    ..Default::default()
                }
                .compile(q.prepared.document(), 1)
                .is_err()
            );
        }
        assert_eq!(*session.document(), before);
    }
}

#[test]
fn actual_shared_proper_stock_pose_retains_both_source_laws_and_measured_gaps() {
    for left in [false, true] {
        let (mut session, mut r) = request();
        for (index, row) in r.receivers.iter().enumerate() {
            let feature = row.source.input.producer;
            let mut parameters = read_primitive(session.document(), feature)
                .unwrap()
                .parameters;
            let PrimitiveParameters::Box { frame, .. } = &mut parameters else {
                unreachable!()
            };
            *frame = PrimitiveFrame {
                origin: [
                    length(10.),
                    length(20.),
                    length(if index == 0 { 30. } else { 10. }),
                ],
                u: [0.8, 0.6, 0.],
                v: [-0.6, 0.8, 0.],
            };
            let edit = prepare_primitive_edit(
                &session,
                session.revision(),
                Actor::local(Channel::Mcp),
                feature,
                &parameters,
            )
            .unwrap();
            session.commit(edit.prepared).unwrap();
        }
        // Explicit authored radial clearance for actual rounded common-pose
        // eccentricity; neither query nor numerical seater changes this size.
        r.receivers[1].source.cut = ConnectionReceiverCut::Round {
            diameter: length(8.1),
        };
        for row in &mut r.hardware {
            let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
                unreachable!()
            };
            match &mut part.parameters {
                MechanicalPartParameters::MachineScrew { thread, .. } => {
                    thread.hand = if left {
                        ThreadHand::Left
                    } else {
                        ThreadHand::Right
                    }
                }
                MechanicalPartParameters::CustomConnection { parameters } => {
                    let CustomConnectionHardwareParameters::WoodInsert {
                        outer_thread,
                        inner_thread,
                        ..
                    } = parameters.as_mut()
                    else {
                        unreachable!()
                    };
                    for t in [outer_thread, inner_thread] {
                        t.hand = if left {
                            ThreadHand::Left
                        } else {
                            ThreadHand::Right
                        };
                    }
                }
                _ => unreachable!(),
            }
        }
        let prepared =
            prepare_screw_connection(&session, session.revision(), Actor::local(Channel::Mcp), &r)
                .unwrap();
        let model = Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile(prepared.prepared.document(), session.revision() + 1)
        .unwrap();
        let check = &model.snapshot().screw_connection_checks[0];
        assert!(check.actual_head_bearing_gap_mm.unwrap()[0] >= 0.);
        assert!(check.actual_insert_inset_mm.unwrap()[0] > 0.);
        assert!(check.actual_insert_inset_mm.unwrap()[1] <= 1e-7);
        assert_eq!(check.actual_engagement_mm.unwrap().round(), 4.);
    }
}
