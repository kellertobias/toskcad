//! Normal canonical connection compilation through a physical blind counterbore.
#[path = "../../multicad-mcp/tests/support/screw_connection.rs"]
mod fixture;
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityId},
    model::*,
};
use multicad_model::{Compiler, Model};
fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}
fn request(
    left: bool,
    washer_height: Option<f64>,
    plain: f64,
) -> (Document, ScrewConnectionRequest) {
    let (document, mut request, _) = fixture::tapped_fixture();
    request.parameters.minimum_engagement = length(1.8);
    request.receivers[0].source.seat = Some(ConnectionReceiverSeat::Counterbore {
        diameter: length(7.),
        depth: length(1.),
    });
    let hand = if left {
        ThreadHand::Left
    } else {
        ThreadHand::Right
    };
    let ConnectionReceiverCut::Tapped { thread } = &mut request.receivers[0].source.cut else {
        unreachable!()
    };
    thread.hand = hand;
    let ScrewConnectionHardwareSource::Generated { part } = &mut request.hardware[0].source else {
        unreachable!()
    };
    let MechanicalPartParameters::MachineScrew {
        thread,
        length_input,
        ..
    } = &mut part.parameters
    else {
        unreachable!()
    };
    thread.hand = hand;
    let underhead = 2. + washer_height.unwrap_or(0.);
    *length_input = Some(ScrewLengthInput::Underhead {
        value: length(underhead),
    });
    thread.length = length(underhead - plain);
    if let Some(height) = washer_height {
        request.hardware.push(ScrewConnectionHardwareRequest {
            role_id: EntityId::new(),
            role: ScrewConnectionHardwareRole::HeadWasher,
            name: "Seat washer".into(),
            instance: EntityId::new(),
            driver: EntityId::new(),
            source: ScrewConnectionHardwareSource::Generated {
                part: Box::new(MechanicalPartRequest::new(
                    "Seat washer",
                    MechanicalPartParameters::Washer {
                        bore_diameter: length(3.4),
                        outside_diameter: length(6.),
                        thickness: length(height),
                        frame: Default::default(),
                    },
                )),
            },
        });
    }
    (document, request)
}
fn apply(session: &mut Session, request: &ScrewConnectionRequest) -> Model {
    let before = session.document().clone();
    let candidate =
        prepare_screw_connection(session, session.revision(), actor(), request).unwrap();
    let model = compiler()
        .compile(candidate.prepared.document(), session.revision() + 1)
        .unwrap();
    assert_eq!(session.document(), &before);
    session.commit(candidate.prepared).unwrap();
    model
}
fn assert_actual(model: &Model, request: &ScrewConnectionRequest, plain: f64) {
    let path = OccurrencePath {
        root: request.hardware[0].instance,
        members: vec![],
    };
    let screw = model.leaf_solid(&path).unwrap();
    let receiver = model
        .leaf_solid(&request.receivers[0].source.occurrence)
        .unwrap();
    let actual = screw.rust_machine_screw_geometry().unwrap().0;
    assert!((actual.spec().shank_length_mm - plain).abs() < 1e-12);
    let proof = screw.rust_screw_blind_thread(receiver).unwrap();
    let seat = proof.counterbore.unwrap();
    let Some(ConnectionReceiverSeat::Counterbore {
        depth: multicad_core::constraint_sketch::Dimension::Literal { value, .. },
        ..
    }) = &request.receivers[0].source.seat
    else {
        panic!("literal test seat")
    };
    assert_eq!((seat.radius_mm, seat.depth_mm), (3.5, *value));
    assert_eq!(proof.threads[1].0.spec().length_mm, 3.5);
    assert!(proof.entry_is_support_plane);
    let [receipt] = model.snapshot().screw_connection_checks.as_slice() else {
        panic!("one physical connection")
    };
    assert!((receipt.actual_engagement_mm.unwrap() - 2.).abs() < 1e-7);
    assert_eq!(receipt.hardware_instances.len(), request.hardware.len());
    assert!(receipt.native_component_collision_checked);
    assert!(model.snapshot().references.all_resolved());
}
#[test]
fn canonical_seated_blind_metric_hands_pose_and_bearing_stacks() {
    for left in [false, true] {
        for (posed, height, plain) in [
            (false, None, 0.),
            (false, Some(0.5), 0.),
            (true, Some(0.625), 0.25),
        ] {
            let (document, request) = request(left, height, plain);
            let mut session = Session::new(document).unwrap();
            if posed {
                let stock = request.receivers[0].source.input.producer;
                let mut primitive = read_primitive(session.document(), stock).unwrap();
                let PrimitiveParameters::Box { frame, .. } = &mut primitive.parameters else {
                    unreachable!()
                };
                *frame = PrimitiveFrame {
                    origin: [length(100.), length(-30.), length(50.)],
                    u: [0.8, 0.6, 0.],
                    v: [-0.36, 0.48, 0.8],
                };
                let candidate = prepare_primitive_edit(
                    &session,
                    session.revision(),
                    actor(),
                    stock,
                    &primitive.parameters,
                )
                .unwrap();
                session.commit(candidate.prepared).unwrap();
            }
            let model = apply(&mut session, &request);
            assert_actual(&model, &request, plain);
            for leaf in model.snapshot().leaf_occurrences.values() {
                let receipt = model
                    .leaf_solid(&leaf.path)
                    .unwrap()
                    .manufacturing_mesh_with_receipt(Default::default())
                    .unwrap();
                assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
            }
        }
    }
}
#[test]
fn seated_blind_edits_history_and_reopen_keep_original_tool_length() {
    use multicad_core::project_file::{self, ProjectFile};
    let (document, mut request) = request(false, Some(0.5), 0.);
    let mut session = Session::new(document).unwrap();
    assert_actual(&apply(&mut session, &request), &request, 0.);
    let before = session.document().clone();
    let old = compiler()
        .compile(session.document(), session.revision())
        .unwrap();
    let topology = &old
        .leaf_snapshot(&request.receivers[0].source.occurrence)
        .unwrap()
        .topology;
    let floor = topology
        .elements
        .iter()
        .find(|face| {
            face.origins
                .iter()
                .any(|o| o.role == "connection.receiver/thread.floor")
        })
        .unwrap();
    let intent = multicad_core::references::GeometryIntent::capture(topology, floor.key).unwrap();
    request.parameters.minimum_engagement = length(1.9);
    request.receivers[0].source.seat = Some(ConnectionReceiverSeat::Counterbore {
        diameter: length(7.),
        depth: length(0.875),
    });
    let edited = apply(&mut session, &request);
    assert_actual(&edited, &request, 0.);
    assert!(matches!(
        multicad_core::references::resolve(
            &intent,
            &edited
                .leaf_snapshot(&request.receivers[0].source.occurrence)
                .unwrap()
                .topology,
            session.revision()
        ),
        multicad_core::references::Resolution::Resolved { .. }
    ));
    let saved = session.document().clone();
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &before);
    session.redo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &saved);
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(saved).unwrap(),
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
    .project;
    assert_eq!(
        read_screw_connection_request(&reopened.document, request.connection).unwrap(),
        request
    );
    assert_actual(
        &compiler()
            .compile(&reopened.document, session.revision())
            .unwrap(),
        &request,
        0.,
    );
}
#[test]
fn seated_blind_invalid_bearing_engagement_and_tip_are_atomic() {
    for case in 0..8 {
        let (document, mut request) = request(false, Some(0.5), 0.);
        let expected = match case {
            0 => {
                request.parameters.minimum_engagement = length(2.2);
                "insufficient engagement"
            }
            1 => {
                request.parameters.blind_tip_clearance = length(0.75);
                "tip violates"
            }
            2 => {
                request.hardware.remove(1);
                request.receivers[0].source.seat = Some(ConnectionReceiverSeat::Countersink {
                    diameter: length(7.),
                    depth: length(1.),
                });
                "conical"
            }
            3 | 4 => {
                let ScrewConnectionHardwareSource::Generated { part } =
                    &mut request.hardware[1].source
                else {
                    unreachable!()
                };
                let MechanicalPartParameters::Washer {
                    bore_diameter,
                    outside_diameter,
                    ..
                } = &mut part.parameters
                else {
                    unreachable!()
                };
                if case == 3 {
                    *bore_diameter = length(2.9);
                    "bore"
                } else {
                    *outside_diameter = length(7.);
                    "counterbore"
                }
            }
            6 => {
                let ScrewConnectionHardwareSource::Generated { part } =
                    &mut request.hardware[1].source
                else {
                    unreachable!()
                };
                let MechanicalPartParameters::Washer { bore_diameter, .. } = &mut part.parameters
                else {
                    unreachable!()
                };
                *bore_diameter = length(5.8);
                "annular land"
            }
            7 => {
                request.receivers[0].source.seat = Some(ConnectionReceiverSeat::Counterbore {
                    diameter: length(5.5),
                    depth: length(1.),
                });
                "head envelope"
            }
            _ => {
                let ScrewConnectionHardwareSource::Generated { part } =
                    &mut request.hardware[0].source
                else {
                    unreachable!()
                };
                let MechanicalPartParameters::MachineScrew { thread, .. } = &mut part.parameters
                else {
                    unreachable!()
                };
                thread.length = length(1.75);
                "plain shank"
            }
        };
        let session = Session::new(document).unwrap();
        let before = session.document().clone();
        let cause = match prepare_screw_connection(&session, 0, actor(), &request) {
            Ok(candidate) => compiler()
                .compile(candidate.prepared.document(), 1)
                .err()
                .expect("invalid geometry must refuse")
                .to_string(),
            Err(cause) => cause.to_string(),
        };
        assert!(cause.contains(expected), "case {case}: {cause}");
        assert_eq!(session.document(), &before);
    }
}

#[test]
fn seated_blind_displaced_and_tilted_actual_washers_refuse_without_commit() {
    use multicad_core::{
        commands::Request,
        document::{EntityData, Mutation},
    };
    let (document, request) = request(false, Some(0.5), 0.);
    let mut session = Session::new(document).unwrap();
    apply(&mut session, &request);
    let before = session.document().clone();
    for (index, displacement, tilted) in [
        (0, 0.499, false),
        (0, 0.501, false),
        (1, 0.499, false),
        (1, 0.5, true),
        (0, 0.5, true),
    ] {
        let row = &request.hardware[index];
        let mut driver = read_screw_hardware_placement(&before, row.driver).unwrap();
        driver.connection = None;
        driver.station = length(displacement);
        if tilted {
            let angle: f64 = 0.01;
            driver.local_mount.rotation = [
                [1., 0., 0.],
                [0., angle.cos(), -angle.sin()],
                [0., angle.sin(), angle.cos()],
            ];
        }
        let mut entity = before.get(row.driver).unwrap().clone();
        entity.data = EntityData::Extension {
            definition: driver.definition().unwrap(),
            inputs: driver.inputs().unwrap(),
        };
        let candidate = session
            .prepare(Request {
                expected_revision: session.revision(),
                actor: actor(),
                label: "Move physical bearing off its seat".into(),
                edits: vec![Mutation::Replace(entity)],
            })
            .unwrap();
        let failure = compiler()
            .compile(candidate.document(), session.revision() + 1)
            .err()
            .expect("nonseated actual material must refuse")
            .to_string();
        assert!(
            failure.contains("bearing") || failure.contains("station"),
            "{failure}"
        );
        assert_eq!(session.document(), &before);
    }
    let failure = compiler()
        .compile_with_backend(
            &before,
            session.revision(),
            multicad_model::GeometryBackend::Rust,
            Some(&|| true),
        )
        .err()
        .expect("cancel must interrupt");
    assert!(multicad_model::is_superseded(&failure));
    assert_eq!(session.document(), &before);
}

#[test]
fn unseated_blind_keeps_physical_flat_washer_support_and_original_thread_frames() {
    let (document, mut request) = request(false, Some(0.5), 0.);
    request.receivers[0].source.seat = None;
    let mut session = Session::new(document).unwrap();
    let model = apply(&mut session, &request);
    let path = OccurrencePath {
        root: request.hardware[0].instance,
        members: vec![],
    };
    let proof = model
        .leaf_solid(&path)
        .unwrap()
        .rust_screw_blind_thread(
            model
                .leaf_solid(&request.receivers[0].source.occurrence)
                .unwrap(),
        )
        .unwrap();
    assert!(proof.counterbore.is_none());
    assert_eq!(proof.threads[1].0.spec().length_mm, 3.5);
    assert!(
        (model.snapshot().screw_connection_checks[0]
            .actual_engagement_mm
            .unwrap()
            - 2.)
            .abs()
            < 1e-7
    );
}

#[test]
fn seated_blind_step_export_refuses_before_publishing_partial_manufacturing_files() {
    use multicad_core::exchange::{Recipe, Scope, ScopeRoot, SetDefinition, prepare_set};
    use multicad_model::export::Exporter;
    let (document, request) = request(false, Some(0.5), 0.);
    let mut session = Session::new(document).unwrap();
    let model = apply(&mut session, &request);
    assert!(
        !model
            .leaf_solid(&request.receivers[0].source.occurrence)
            .unwrap()
            .supports_step_export()
    );
    let set = EntityId::new();
    let candidate = prepare_set(
        &session,
        session.revision(),
        actor(),
        set,
        "Blind seated assembly STEP",
        &SetDefinition {
            scope: Scope {
                root: ScopeRoot::Project,
                tags: None,
            },
            recipes: vec![Recipe::step("STEP")],
            ..Default::default()
        },
    )
    .unwrap();
    session.commit(candidate).unwrap();
    let before = session.document().clone();
    let exporter = Exporter::default();
    let preview = exporter
        .preview(&before, session.revision(), set, &actor(), None)
        .unwrap();
    assert!(!preview.can_run(), "{preview:?}");
    let directory = tempfile::tempdir().unwrap();
    assert!(
        exporter
            .run(
                &before,
                session.revision(),
                &preview.key,
                directory.path(),
                &actor(),
                None
            )
            .is_err()
    );
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
    assert_eq!(session.document(), &before);
}

#[test]
fn external_hex_actual_head_envelope_and_flat_land_fit_the_blind_seat() {
    for left in [false, true] {
        let (document, mut request) = request(left, None, 0.);
        let ScrewConnectionHardwareSource::Generated { part } = &mut request.hardware[0].source
        else {
            unreachable!()
        };
        let MechanicalPartParameters::MachineScrew { head, drive, .. } = &mut part.parameters
        else {
            unreachable!()
        };
        *head = ScrewHead::ExternalHex {
            across_flats: length(5.5),
            height: length(2.),
        };
        *drive = ScrewDrive::ExternalHex;
        let mut session = Session::new(document).unwrap();
        let model = apply(&mut session, &request);
        assert_actual(&model, &request, 0.);
    }
}
