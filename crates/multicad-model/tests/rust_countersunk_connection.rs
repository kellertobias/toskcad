//! Complete canonical generated hardware, receiver, occurrence and bearing path.
#[path = "../../multicad-mcp/tests/support/screw_connection.rs"]
mod fixture;
use multicad_core::{
    commands::{Actor, Channel, Session},
    model::*,
};
use multicad_model::Compiler;
fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn countersunk(left: bool) -> (multicad_core::document::Document, ScrewConnectionRequest) {
    let (document, mut request, _) = fixture::fixture();
    request.receivers[0].source.seat = Some(ConnectionReceiverSeat::Countersink {
        diameter: length(7.),
        depth: length(1.8),
    });
    for row in &mut request.hardware {
        let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
            unreachable!()
        };
        match &mut part.parameters {
            MechanicalPartParameters::MachineScrew {
                head,
                drive,
                length_input,
                underhead_length,
                thread,
                ..
            } => {
                *head = ScrewHead::Countersunk {
                    profile: Box::new(CountersunkScrewHeadParameters {
                        diameter: length(6.),
                        height: length(1.65),
                        included_angle: angle(90.),
                        edge_chamfer: length(0.1),
                        neck_radius: length(0.1),
                    }),
                };
                *drive = ScrewDrive::HexSocket {
                    across_flats: length(2.5),
                    depth: length(0.8),
                };
                *underhead_length = None;
                *length_input = Some(ScrewLengthInput::Underhead { value: length(8.) });
                thread.length = length(7.5);
                thread.hand = if left {
                    ThreadHand::Left
                } else {
                    ThreadHand::Right
                };
            }
            MechanicalPartParameters::HexNut { thread, .. } => {
                thread.hand = if left {
                    ThreadHand::Left
                } else {
                    ThreadHand::Right
                };
            }
            _ => unreachable!(),
        }
    }
    (document, request)
}
#[test]
fn canonical_countersunk_metric_through_nut_compiles_actual_head_and_receiver() {
    for left in [false, true] {
        let (document, request) = countersunk(left);
        let session = Session::new(document).unwrap();
        let candidate = prepare_screw_connection(&session, 0, actor(), &request).unwrap();
        let model = Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile(candidate.prepared.document(), 1)
        .unwrap();
        let [check] = model.snapshot().screw_connection_checks.as_slice() else {
            panic!("actual connection proof");
        };
        assert!(check.native_component_collision_checked);
        assert!((check.actual_engagement_mm.unwrap() - 2.4).abs() < 1e-7);
        assert!(check.scope.contains("source-certified conical bearing"));
    }
}

fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}
fn valid(left: bool, posed: bool) -> (Session, ScrewConnectionRequest, multicad_model::Model) {
    let (document, request) = countersunk(left);
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
    let original = session.document().clone();
    let candidate =
        prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
    let model = compiler()
        .compile(candidate.prepared.document(), session.revision() + 1)
        .unwrap();
    assert_eq!(session.document(), &original);
    session.commit(candidate.prepared).unwrap();
    (session, request, model)
}
fn path(instance: multicad_core::document::EntityId) -> OccurrencePath {
    OccurrencePath {
        root: instance,
        members: vec![],
    }
}
fn assert_actual_bearing(model: &multicad_model::Model, request: &ScrewConnectionRequest) {
    use multicad_kernel::rust_kernel::conical_screw_seat as seat;
    let screw = model
        .leaf_solid(&path(request.hardware[0].instance))
        .unwrap();
    let receiver = model
        .leaf_solid(&request.receivers[0].source.occurrence)
        .unwrap();
    let (head, frame) = screw.rust_profiled_screw_geometry().unwrap();
    let cavities = receiver.rust_countersink_geometry().unwrap();
    let [cavity] = cavities.as_slice() else {
        panic!("actual single receiving cone");
    };
    assert_eq!(cavity.origin.source, request.receivers[0].feature);
    assert_eq!(head.profile().end()[0], 1.5);
    assert!(head.thread().spec().major_radius_mm < 1.5);
    assert_eq!(head.spec().shank_length_mm, 0.5);
    assert!(!head.has_head_thread_overlap());
    assert!(head.thread().spec().root_radius_mm.is_some());
    let policy = multicad_kernel::GeometryPolicy::default();
    let proof = seat::certify_between_frames(
        head.profile(),
        seat::Seat {
            bore_radius_mm: cavity.bore_radius_mm,
            mouth_radius_mm: cavity.mouth_radius_mm,
            depth_mm: cavity.depth_mm,
        },
        seat::Options {
            tolerance_mm: policy.linear_tolerance_mm,
            minimum_land_mm: 5.
                * policy
                    .minimum_feature_mm
                    .max(16. * policy.linear_tolerance_mm),
        },
        seat::Placement {
            origin_mm: frame.origin,
            axes: [frame.u, frame.v, frame.normal()],
            axial_offset_mm: -head.spec().shank_length_mm,
        },
        seat::Frame {
            origin_mm: cavity.frame.origin,
            axes: [cavity.frame.u, cavity.frame.v, cavity.frame.normal()],
            offset_mm: cavity.offset_mm,
            source_error_mm: cavity.source_error_mm,
        },
    )
    .unwrap();
    assert!(
        proof.maximum_physical_gap_mm > 0.
            && proof.maximum_physical_gap_mm <= policy.linear_tolerance_mm
    );
    assert!(proof.minimum_physical_land_mm >= 5. * policy.minimum_feature_mm);
    let [check] = model.snapshot().screw_connection_checks.as_slice() else {
        panic!("complete connection proof");
    };
    assert!(check.native_component_collision_checked);
    assert!((check.actual_engagement_mm.unwrap() - 2.4).abs() < policy.linear_tolerance_mm);
    assert!(model.snapshot().references.all_resolved());
    assert_eq!(model.snapshot().leaf_occurrences.len(), 3);
}

#[test]
fn canonical_counter_history_reopen_and_physical_outputs_survive_both_hands_and_pose() {
    use multicad_core::{
        document::EntityId,
        exchange::{Recipe, Scope, ScopeRoot, SetDefinition, prepare_set},
        project_file::{self, ProjectFile},
        references::{AnalyticGeometry, GeometryIntent, Resolution, TopologyKind, resolve},
    };
    use multicad_model::export::Exporter;
    for left in [false, true] {
        for posed in [false, true] {
            let (mut session, mut request, original) = valid(left, posed);
            assert_actual_bearing(&original, &request);
            let head_path = path(request.hardware[0].instance);
            let cavity_path = request.receivers[0].source.occurrence.clone();
            let old_head = &original.leaf_snapshot(&head_path).unwrap().topology;
            let cap = old_head
                .elements
                .iter()
                .find(|e| {
                    e.key.kind == TopologyKind::Face
                        && e.origins.iter().any(|o| o.role == "screw.head/top.land")
                })
                .unwrap();
            let cap_intent = GeometryIntent::capture(old_head, cap.key).unwrap();
            let old_cavity = &original.leaf_snapshot(&cavity_path).unwrap().topology;
            let mouth=old_cavity.elements.iter().find(|e|e.key.kind==TopologyKind::Edge && matches!(e.geometry,Some(AnalyticGeometry::CircularArc{radius_mm,..})if (radius_mm-3.5).abs()<1e-7)).unwrap();
            let mouth_intent = GeometryIntent::capture(old_cavity, mouth.key).unwrap();
            let before = session.document().clone();
            let ScrewConnectionHardwareSource::Generated { part } = &mut request.hardware[0].source
            else {
                unreachable!()
            };
            let MechanicalPartParameters::MachineScrew {
                length_input,
                thread,
                ..
            } = &mut part.parameters
            else {
                unreachable!()
            };
            *length_input = Some(ScrewLengthInput::Underhead { value: length(9.) });
            thread.length = length(8.5);
            let candidate =
                prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
            let edited = compiler()
                .compile(candidate.prepared.document(), session.revision() + 1)
                .unwrap();
            assert_actual_bearing(&edited, &request);
            assert!(matches!(
                resolve(
                    &cap_intent,
                    &edited.leaf_snapshot(&head_path).unwrap().topology,
                    session.revision() + 1
                ),
                Resolution::Resolved { .. }
            ));
            assert!(matches!(
                resolve(
                    &mouth_intent,
                    &edited.leaf_snapshot(&cavity_path).unwrap().topology,
                    session.revision() + 1
                ),
                Resolution::Resolved { .. }
            ));
            assert_eq!(session.document(), &before);
            session.commit(candidate.prepared).unwrap();
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
            let reopened_model = compiler()
                .compile(&reopened.document, session.revision())
                .unwrap();
            assert_actual_bearing(&reopened_model, &request);
            for leaf in reopened_model.snapshot().leaf_occurrences.values() {
                let solid = reopened_model.leaf_solid(&leaf.path).unwrap();
                let receipt = solid
                    .manufacturing_mesh_with_receipt(Default::default())
                    .unwrap();
                assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
            }
            let set = EntityId::new();
            let candidate = prepare_set(
                &session,
                session.revision(),
                actor(),
                set,
                "Countersunk assembly STEP",
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
            let exporter = Exporter::default();
            let preview = exporter
                .preview(session.document(), session.revision(), set, &actor(), None)
                .unwrap();
            assert!(preview.can_run(), "{preview:?}");
            assert_eq!(preview.artifacts.len(), 3);
            let directory = tempfile::tempdir().unwrap();
            let result = exporter
                .run(
                    session.document(),
                    session.revision(),
                    &preview.key,
                    directory.path(),
                    &actor(),
                    None,
                )
                .unwrap();
            assert_eq!(result.files.len(), 3);
            for file in result.files {
                let filename =
                    std::path::Path::new(&result.directory).join(file.artifact.output.filename);
                let actual =
                    multicad_kernel::Solid::read_step(&filename, Default::default()).unwrap();
                assert_eq!(actual.solid_count().unwrap(), 1);
                assert!((actual.volume_mm3().unwrap() - file.artifact.volume_mm3).abs() < 1e-4);
            }
        }
    }
}

#[test]
fn canonical_counter_refuses_flat_washer_insufficient_engagement_and_floating_hardware() {
    use multicad_core::{
        commands::Request,
        document::{EntityData, EntityId, Mutation},
    };
    for left in [false, true] {
        let (document, mut request) = countersunk(left);
        let session = Session::new(document.clone()).unwrap();
        request.hardware.insert(
            1,
            ScrewConnectionHardwareRequest {
                role_id: EntityId::new(),
                role: ScrewConnectionHardwareRole::HeadWasher,
                name: "Incompatible flat head washer".into(),
                source: ScrewConnectionHardwareSource::Generated {
                    part: Box::new(MechanicalPartRequest::new(
                        "Flat washer",
                        MechanicalPartParameters::Washer {
                            bore_diameter: length(3.4),
                            outside_diameter: length(7.),
                            thickness: length(0.2),
                            frame: Default::default(),
                        },
                    )),
                },
                instance: EntityId::new(),
                driver: EntityId::new(),
            },
        );
        let error = prepare_screw_connection(&session, 0, actor(), &request)
            .err()
            .unwrap();
        assert!(
            error.to_string().to_lowercase().contains("washer"),
            "{error}"
        );
        assert_eq!(session.document(), &document);
        request.hardware.remove(1);
        request.parameters.minimum_engagement = length(4.);
        let candidate = prepare_screw_connection(&session, 0, actor(), &request).unwrap();
        let error = compiler()
            .compile(candidate.prepared.document(), 1)
            .err()
            .unwrap();
        assert!(
            error.to_string().to_lowercase().contains("engagement"),
            "{error}"
        );
        assert_eq!(session.document(), &document);

        let (session, request, _) = valid(left, true);
        let saved = session.document().clone();
        // A positive/negative station displacement changes the actual bearing,
        // not just the connection recipe. Neither candidate may be committed.
        for (station, reason) in [(1.999, "Unresolved"), (2.001, "Interference")] {
            let mut driver =
                read_screw_hardware_placement(&saved, request.hardware[0].driver).unwrap();
            driver.connection = None;
            driver.station = length(station);
            driver.roll = angle(0.);
            let mut entity = saved.get(request.hardware[0].driver).unwrap().clone();
            entity.data = EntityData::Extension {
                definition: driver.definition().unwrap(),
                inputs: driver.inputs().unwrap(),
            };
            let candidate = session
                .prepare(Request {
                    expected_revision: session.revision(),
                    actor: actor(),
                    label: "Displace screw from conical bearing".into(),
                    edits: vec![Mutation::Replace(entity)],
                })
                .unwrap();
            let error = compiler()
                .compile(candidate.document(), session.revision() + 1)
                .err()
                .unwrap();
            let message = error.to_string();
            assert!(
                message.contains("Actual finite conical bearing cannot be certified"),
                "{message}"
            );
            assert!(message.contains(reason), "{message}");
            assert_eq!(session.document(), &saved);
        }
    }
}
