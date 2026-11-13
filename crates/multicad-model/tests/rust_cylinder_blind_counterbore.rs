use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    document::{Document, Entity, EntityData, EntityId, Mutation},
    model::*,
    references::AnalyticGeometry,
};
use multicad_model::Compiler;
fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}
fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}
fn fixture(end: bool) -> (Session, ConnectionReceiverRequest) {
    let mut s = Session::new(Document::new("Native receiving cuts").unwrap()).unwrap();
    let stock = PrimitiveRequest::new(
        "Stock",
        PrimitiveParameters::Cylinder {
            radius: length(10.),
            height: length(10.),
            frame: PrimitiveFrame {
                origin: [length(10.), length(10.), length(0.)],
                ..Default::default()
            },
        },
    );
    s.commit(prepare_primitive(&s, 0, actor(), &stock).unwrap().prepared)
        .unwrap();
    let place = ComponentPlacementRequest::new(
        stock.body,
        NumericTransform {
            translation: [length(30.), length(-20.), length(40.)],
            axis: [1., 1., 1.],
            angle: angle(37.),
            ..Default::default()
        },
    );
    s.commit(
        prepare_place_component(&s, s.revision(), actor(), &place)
            .unwrap()
            .prepared,
    )
    .unwrap();
    let model = compiler().compile(s.document(), s.revision()).unwrap();
    let snap = &model.snapshot().shapes[&stock.feature].topology;
    let face = snap
        .elements
        .iter()
        .find(|e| {
            matches!(e.geometry,Some(AnalyticGeometry::Plane{frame})if if end {frame.normal()[2]>0.99}else{frame.normal()[2]< -0.99})
                && !e.frames.is_empty()
        })
        .unwrap();
    let capture = capture_profile_face(snap, face.key, None).unwrap();
    let frame = face.frames[0].frame;
    let point = [10., 10., if end { 10. } else { 0. }];
    let d: [f64; 3] = std::array::from_fn(|i| point[i] - frame.origin_mm[i]);
    let dot = |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let path = OccurrencePath {
        root: place.instance,
        members: vec![],
    };
    let mut axis = ScrewConnectionAxisDefinition::new(
        BodyInput {
            body: stock.body,
            producer: stock.feature,
        },
        path.clone(),
        capture.intent.id,
    );
    axis.reverse = true;
    axis.offset = [length(dot(d, frame.u)), length(dot(d, frame.v)), length(0.)];
    let entity = Entity::new(
        "Connection axis",
        EntityData::Extension {
            definition: axis.definition().unwrap(),
            inputs: axis.inputs(&capture).unwrap(),
        },
    );
    let aid = entity.id;
    s.apply(Request {
        expected_revision: s.revision(),
        actor: actor(),
        label: "Capture axis".into(),
        edits: vec![Mutation::Insert(entity)],
    })
    .unwrap();
    (
        s,
        ConnectionReceiverRequest {
            feature: EntityId::new(),
            name: "Receiving cut".into(),
            source: ConnectionReceiverDefinition {
                entry_reference: None,
                output: BodyOutput { body: stock.body },
                input: BodyInput {
                    body: stock.body,
                    producer: stock.feature,
                },
                occurrence: path,
                axis: aid,
                station: length(0.),
                cut: ConnectionReceiverCut::Round {
                    diameter: length(4.),
                },
                extent: ConnectionReceiverExtent::Through {},
                minimum_wall: length(1.),
                seat: None,
            },
        },
    )
}

#[test]
fn canonical_cylinder_stock_blind_counterbore_has_a_real_shoulder_and_floor() {
    let (session, mut request) = fixture(true);
    request.source.extent = ConnectionReceiverExtent::Blind {
        depth: length(4.),
        floor: length(1.),
    };
    request.source.seat = Some(ConnectionReceiverSeat::Counterbore {
        diameter: length(8.),
        depth: length(2.),
    });
    let before = session.document().clone();
    let candidate =
        prepare_connection_receiver(&session, session.revision(), actor(), &request).unwrap();
    let model = compiler().compile_with_backend(
        candidate.prepared.document(),
        session.revision() + 1,
        multicad_model::GeometryBackend::Rust,
        None,
    );
    assert_eq!(session.document(), &before);
    let model =
        model.unwrap_or_else(|cause| panic!("normal cylinder-stock blind counterbore: {cause}"));
    assert_eq!(
        model.solid(request.feature).unwrap().solid_count().unwrap(),
        1
    );
}

#[test]
fn both_caps_plain_and_seated_sources_keep_geometry_display_and_verified_step() {
    for end in [false, true] {
        for seat in [false, true] {
            let (mut session, mut request) = fixture(end);
            request.source.extent = ConnectionReceiverExtent::Blind {
                depth: length(4.),
                floor: length(1.),
            };
            if seat {
                request.source.seat = Some(ConnectionReceiverSeat::Counterbore {
                    diameter: length(8.),
                    depth: length(2.),
                });
            }
            let candidate =
                prepare_connection_receiver(&session, session.revision(), actor(), &request)
                    .unwrap();
            let model = Compiler::default()
                .compile_with_backend(
                    candidate.prepared.document(),
                    session.revision() + 1,
                    multicad_model::GeometryBackend::Rust,
                    None,
                )
                .unwrap();
            let solid = model.solid(request.feature).unwrap();
            let (source, _) = solid.rust_cylindrical_blind_pocket_geometry().unwrap();
            assert_eq!(source.end(), end);
            assert_eq!(source.spec().seat.is_some(), seat);
            assert!(
                !model.snapshot().shapes[&request.feature]
                    .mesh
                    .indices
                    .is_empty()
            );
            solid
                .manufacturing_mesh_with_receipt(Default::default())
                .unwrap();
            session.commit(candidate.prepared).unwrap();
            use multicad_core::exchange::{Recipe, Scope, ScopeRoot, SetDefinition, prepare_set};
            use multicad_model::export::Exporter;
            let set = EntityId::new();
            let prepared = prepare_set(
                &session,
                session.revision(),
                actor(),
                set,
                "Blind pocket STEP",
                &SetDefinition {
                    scope: Scope {
                        root: ScopeRoot::Explicit {
                            entities: vec![request.source.occurrence.root],
                        },
                        tags: None,
                    },
                    recipes: vec![Recipe::step("STEP")],
                    ..Default::default()
                },
            )
            .unwrap();
            session.commit(prepared).unwrap();
            let before = session.document().clone();
            let exporter = Exporter {
                compiler: compiler(),
            };
            let preview = exporter
                .preview(&before, session.revision(), set, &actor(), None)
                .unwrap();
            assert!(preview.can_run(), "{preview:?}");
            let directory = tempfile::tempdir().unwrap();
            let receipt = exporter
                .run(
                    &before,
                    session.revision(),
                    &preview.key,
                    directory.path(),
                    &actor(),
                    None,
                )
                .unwrap();
            assert_eq!(receipt.files.len(), 1);
            let file = std::path::Path::new(&receipt.directory)
                .join(&receipt.files[0].artifact.output.filename);
            let text = std::fs::read_to_string(file).unwrap();
            let (read, _) =
                multicad_kernel::rust_kernel::step_read::read_import_solid(&text, 1e-7).unwrap();
            assert_eq!(read.face_count, if seat { 7 } else { 5 });
            let volume = std::f64::consts::PI * if seat { 960. } else { 984. };
            assert!((read.volume_mm3 - volume).abs() < 1e-8);
            assert_eq!(session.document(), &before);
        }
    }
}
#[test]
fn pocket_edits_are_immutable_and_history_and_archive_rebuild_current_source() {
    use multicad_core::project_file::{self, ProjectFile};
    let (mut session, mut request) = fixture(true);
    request.source.extent = ConnectionReceiverExtent::Blind {
        depth: length(4.),
        floor: length(1.),
    };
    request.source.seat = Some(ConnectionReceiverSeat::Counterbore {
        diameter: length(8.),
        depth: length(2.),
    });
    let candidate =
        prepare_connection_receiver(&session, session.revision(), actor(), &request).unwrap();
    let old = compiler()
        .compile(candidate.prepared.document(), session.revision() + 1)
        .unwrap();
    session.commit(candidate.prepared).unwrap();
    let before = session.document().clone();
    request.source.extent = ConnectionReceiverExtent::Blind {
        depth: length(5.),
        floor: length(1.),
    };
    let edit = prepare_connection_receiver_edit(
        &session,
        session.revision(),
        actor(),
        request.feature,
        &request.source,
    )
    .unwrap();
    let edited = compiler()
        .compile(edit.prepared.document(), session.revision() + 1)
        .unwrap();
    assert_eq!(
        old.solid(request.feature)
            .unwrap()
            .rust_cylindrical_blind_pocket_geometry()
            .unwrap()
            .0
            .spec()
            .depth_mm,
        4.
    );
    assert_eq!(
        edited
            .solid(request.feature)
            .unwrap()
            .rust_cylindrical_blind_pocket_geometry()
            .unwrap()
            .0
            .spec()
            .depth_mm,
        5.
    );
    session.commit(edit.prepared).unwrap();
    let after = session.document().clone();
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &before);
    session.redo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &after);
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(after.clone()).unwrap(),
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
    assert_eq!(reopened.document, after);
    assert_eq!(
        compiler()
            .compile(&reopened.document, 0)
            .unwrap()
            .solid(request.feature)
            .unwrap()
            .rust_cylindrical_blind_pocket_geometry()
            .unwrap()
            .0
            .spec()
            .depth_mm,
        5.
    );
}
#[test]
fn invalid_material_conical_seats_and_cancellation_leave_the_session_unchanged() {
    for variant in 0..4 {
        let (session, mut request) = fixture(true);
        request.source.extent = ConnectionReceiverExtent::Blind {
            depth: length(if variant == 0 { 9.5 } else { 4. }),
            floor: length(1.),
        };
        request.source.minimum_wall = length(if variant == 1 { 7. } else { 1. });
        request.source.seat = Some(if variant == 2 {
            ConnectionReceiverSeat::Countersink {
                diameter: length(8.),
                depth: length(2.),
            }
        } else {
            ConnectionReceiverSeat::Counterbore {
                diameter: length(8.),
                depth: length(2.),
            }
        });
        let before = session.document().clone();
        let candidate =
            prepare_connection_receiver(&session, session.revision(), actor(), &request).unwrap();
        let cancel = || variant == 3;
        assert!(
            compiler()
                .compile_with_backend(
                    candidate.prepared.document(),
                    session.revision() + 1,
                    multicad_model::GeometryBackend::Rust,
                    Some(&cancel)
                )
                .is_err()
        );
        assert_eq!(session.document(), &before);
    }
}
