use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    document::{Document, Entity, EntityData, EntityId, Mutation},
    model::*,
    project_file::{self, ProjectFile},
    references::{AnalyticGeometry, TopologyKind},
};
use multicad_kernel::{GeometryPolicy, Solid};
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
fn fixture() -> (Session, ConnectionReceiverRequest) {
    let mut s = Session::new(Document::new("Native receiving cuts").unwrap()).unwrap();
    let stock = PrimitiveRequest::new(
        "Stock",
        PrimitiveParameters::Box {
            size: [length(20.), length(20.), length(10.)],
            frame: Default::default(),
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
            matches!(e.geometry,Some(AnalyticGeometry::Plane{frame})if frame.normal()[2]>0.99)
                && !e.frames.is_empty()
        })
        .unwrap();
    let capture = capture_profile_face(snap, face.key, None).unwrap();
    let frame = face.frames[0].frame;
    let point = [10., 10., 10.];
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
fn add(s: &mut Session, r: &ConnectionReceiverRequest) {
    let q = prepare_connection_receiver(s, s.revision(), actor(), r).unwrap();
    s.commit(q.prepared).unwrap();
}
#[test]
fn real_cylinders_blind_floors_and_seats_keep_local_geometry_history_and_step() {
    let dir = tempfile::tempdir().unwrap();
    for (i, seat, blind, removed) in [
        (0, None, false, 40. * std::f64::consts::PI),
        (
            1,
            Some(ConnectionReceiverSeat::Counterbore {
                diameter: length(8.),
                depth: length(2.),
            }),
            false,
            64. * std::f64::consts::PI,
        ),
        (
            2,
            Some(ConnectionReceiverSeat::Countersink {
                diameter: length(8.),
                depth: length(2.),
            }),
            false,
            152. * std::f64::consts::PI / 3.,
        ),
        (3, None, true, 16. * std::f64::consts::PI),
    ] {
        let (mut s, mut r) = fixture();
        r.source.seat = seat;
        if blind {
            r.source.extent = ConnectionReceiverExtent::Blind {
                depth: length(4.),
                floor: length(1.),
            };
        }
        add(&mut s, &r);
        let original = s.document().clone();
        let model = compiler().compile(s.document(), s.revision()).unwrap();
        let shape = model.solid(r.feature).unwrap();
        let before = model.solid(r.source.input.producer).unwrap();
        assert!((before.volume_mm3().unwrap() - 4000.).abs() < 1e-6);
        assert!(
            (shape.volume_mm3().unwrap() - (4000. - removed)).abs() < 1e-5,
            "case{i}:{}",
            shape.volume_mm3().unwrap()
        );
        assert_eq!(shape.solid_count().unwrap(), 1);
        let bounds = shape.bounds().unwrap();
        for k in 0..3 {
            assert!((bounds.min[k]).abs() < 1e-6);
            assert!((bounds.max[k] - [20., 20., 10.][k]).abs() < 1e-6);
        }
        let topology = shape.topology_snapshot(r.feature, s.revision()).unwrap();
        assert!(
            topology
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Face)
                .all(|e| !e.origins.is_empty())
        );
        assert!(topology.elements.iter().any(|e| {
            e.origins
                .iter()
                .any(|o| o.source == r.feature && o.role.contains("connection.receiver/"))
        }));
        let path = dir.path().join(format!("receiving-{i}.step"));
        shape.write_step(&path).unwrap();
        let back = Solid::read_step(&path, GeometryPolicy::default()).unwrap();
        let area = shape
            .faces()
            .unwrap()
            .iter()
            .map(|f| f.area_mm2)
            .sum::<f64>();
        assert!(
            (back.volume_mm3().unwrap() - shape.volume_mm3().unwrap()).abs()
                <= area * 1e-7 + 1e-7 * shape.volume_mm3().unwrap()
        );
        assert_eq!(*s.document(), original);
        let mut bytes = vec![];
        project_file::write(
            &ProjectFile::new(s.document().clone()).unwrap(),
            &mut bytes,
            Default::default(),
        )
        .unwrap();
        let reopen = project_file::read(
            &mut std::io::Cursor::new(bytes),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap();
        let build = compiler().compile(&reopen.project.document, 0).unwrap();
        assert!(
            (build.solid(r.feature).unwrap().volume_mm3().unwrap() - shape.volume_mm3().unwrap())
                .abs()
                < 1e-6
        );
    }
}
#[test]
fn off_part_inaccessible_entry_breakout_floor_and_tapped_extent_fail_without_live_commit() {
    for case in 0..5 {
        let (s, mut r) = fixture();
        match case {
            0 => r.source.station = length(1.),
            1 => r.source.station = length(-1.),
            2 => r.source.minimum_wall = length(20.),
            3 => {
                r.source.extent = ConnectionReceiverExtent::Blind {
                    depth: length(9.5),
                    floor: length(1.),
                }
            }
            _ => {
                let thread = ThreadParameters {
                    major_diameter: length(3.),
                    pitch: length(0.5),
                    length: length(3.),
                    fit: ThreadFit::Metric6HInternal,
                    ..Default::default()
                };
                r.source.cut = ConnectionReceiverCut::Tapped {
                    thread: Box::new(thread),
                };
            }
        }
        let before = s.document().clone();
        let candidate = prepare_connection_receiver(&s, s.revision(), actor(), &r).unwrap();
        assert!(
            compiler()
                .compile(candidate.prepared.document(), s.revision() + 1)
                .is_err(),
            "case{case}"
        );
        assert_eq!(*s.document(), before);
    }
}
#[test]
fn real_tapped_receiver_uses_same_axis_and_retains_source_length() {
    let (mut s, mut r) = fixture();
    let thread = ThreadParameters {
        major_diameter: length(3.),
        pitch: length(0.5),
        length: length(3.),
        fit: ThreadFit::Metric6HInternal,
        minimum_wall: length(1.),
        ..Default::default()
    };
    r.source.cut = ConnectionReceiverCut::Tapped {
        thread: Box::new(thread),
    };
    r.source.extent = ConnectionReceiverExtent::Blind {
        depth: length(3.),
        floor: length(1.),
    };
    add(&mut s, &r);
    let model = compiler().compile(s.document(), s.revision()).unwrap();
    let shape = model.solid(r.feature).unwrap();
    assert_eq!(shape.solid_count().unwrap(), 1);
    assert!(shape.volume_mm3().unwrap() < 4000.);
    assert!(
        model
            .solid(r.source.input.producer)
            .unwrap()
            .volume_mm3()
            .unwrap()
            > shape.volume_mm3().unwrap()
    );
    assert_eq!(
        read_connection_receiver(s.document(), r.feature).unwrap(),
        r.source
    );
}

#[test]
fn plain_blind_receiving_depth_edits_keep_actual_floor_references_history_and_reopen() {
    use multicad_core::references::{GeometryIntent, Resolution, resolve};
    let (mut session, mut request) = fixture();
    request.source.extent = ConnectionReceiverExtent::Blind {
        depth: length(4.),
        floor: length(1.),
    };
    let q = prepare_connection_receiver(&session, session.revision(), actor(), &request).unwrap();
    let original = compiler()
        .compile(q.prepared.document(), session.revision() + 1)
        .unwrap();
    session.commit(q.prepared).unwrap();
    let old = &original.snapshot().shapes[&request.feature].topology;
    let floor = old
        .elements
        .iter()
        .find(|e| {
            e.key.kind == TopologyKind::Face
                && e.origins
                    .iter()
                    .any(|o| o.role.starts_with("pocket.floor/") && o.source == request.feature)
        })
        .unwrap();
    let intent = GeometryIntent::capture(old, floor.key).unwrap();
    let before = session.document().clone();
    request.source.extent = ConnectionReceiverExtent::Blind {
        depth: length(5.),
        floor: length(1.),
    };
    let q = prepare_connection_receiver_edit(
        &session,
        session.revision(),
        actor(),
        request.feature,
        &request.source,
    )
    .unwrap();
    let edited = compiler()
        .compile(q.prepared.document(), session.revision() + 1)
        .unwrap();
    assert!(matches!(
        resolve(
            &intent,
            &edited.snapshot().shapes[&request.feature].topology,
            session.revision() + 1
        ),
        Resolution::Resolved { .. }
    ));
    assert!(
        (edited.solid(request.feature).unwrap().volume_mm3().unwrap()
            - (4000. - 20. * std::f64::consts::PI))
            .abs()
            < 1e-7
    );
    assert_eq!(session.document(), &before);
    session.commit(q.prepared).unwrap();
    let saved = session.document().clone();
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &before);
    session.redo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &saved);
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(saved.clone()).unwrap(),
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
    assert_eq!(reopened.document, saved);
    let model = compiler()
        .compile(&reopened.document, session.revision())
        .unwrap();
    assert!(
        (model.solid(request.feature).unwrap().volume_mm3().unwrap()
            - edited.solid(request.feature).unwrap().volume_mm3().unwrap())
        .abs()
            < 1e-9
    );
    let receipt = model
        .solid(request.feature)
        .unwrap()
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
}
#[test]
fn plain_blind_receivers_refuse_uncertified_counterbore_pose_cones_wrong_entry_and_material_failures_atomically()
 {
    for case in 0..6 {
        let (session, mut request) = fixture();
        request.source.extent = ConnectionReceiverExtent::Blind {
            depth: length(4.),
            floor: length(1.),
        };
        match case {
            0 => {
                request.source.seat = Some(ConnectionReceiverSeat::Counterbore {
                    diameter: length(8.),
                    depth: length(2.),
                })
            }
            1 => {
                request.source.seat = Some(ConnectionReceiverSeat::Countersink {
                    diameter: length(8.),
                    depth: length(2.),
                })
            }
            2 => request.source.station = length(0.1),
            3 => request.source.station = length(-0.1),
            4 => request.source.minimum_wall = length(20.),
            _ => {
                request.source.extent = ConnectionReceiverExtent::Blind {
                    depth: length(9.5),
                    floor: length(1.),
                }
            }
        };
        let before = session.document().clone();
        let q =
            prepare_connection_receiver(&session, session.revision(), actor(), &request).unwrap();
        let failure = compiler()
            .compile(q.prepared.document(), session.revision() + 1)
            .err()
            .expect("blind source guard");
        assert_eq!(failure.entity, Some(request.feature));
        if case == 0 {
            assert!(
                failure.message.contains("exact cardinal entry chart"),
                "{failure}"
            );
        } else if case == 1 {
            assert!(failure.message.contains("conical seats"), "{failure}");
        }
        assert_eq!(session.document(), &before);
    }
}

#[test]
fn plain_blind_receiver_cancellation_does_not_publish_a_candidate() {
    let (session, mut request) = fixture();
    request.source.extent = ConnectionReceiverExtent::Blind {
        depth: length(4.),
        floor: length(1.),
    };
    let before = session.document().clone();
    let q = prepare_connection_receiver(&session, session.revision(), actor(), &request).unwrap();
    let failure = compiler()
        .compile_with_backend(
            q.prepared.document(),
            session.revision() + 1,
            multicad_model::GeometryBackend::Rust,
            Some(&|| true),
        )
        .err()
        .unwrap();
    assert!(multicad_model::is_superseded(&failure));
    assert_eq!(session.document(), &before);
}
