use multicad_core::{
    commands::{Actor, Channel, Session},
    model::*,
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, GeometryBackend};
#[path = "support/custom_screw_connections.rs"]
mod support;
fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}
fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}
fn fixture(posed: bool, left: bool, starts: u32) -> (Session, ScrewConnectionRequest) {
    let (doc, mut request) = support::fixture(ScrewConnectionMode::DirectWood);
    let mut session = Session::new(doc).unwrap();
    if posed {
        let feature = request.receivers[0].source.input.producer;
        let mut parameters = read_primitive(session.document(), feature)
            .unwrap()
            .parameters;
        let PrimitiveParameters::Box { frame, .. } = &mut parameters else {
            unreachable!()
        };
        *frame = PrimitiveFrame {
            origin: [length(10.), length(20.), length(30.)],
            u: [0.8, 0.6, 0.],
            v: [-0.6, 0.8, 0.],
        };
        session
            .commit(
                prepare_primitive_edit(&session, 0, actor(), feature, &parameters)
                    .unwrap()
                    .prepared,
            )
            .unwrap();
        // Explicit radial clearance accommodates the actual rounded relative
        // placement; the certificate never repairs a shared-radius pose.
        request.receivers[0].source.cut = ConnectionReceiverCut::Round {
            diameter: length(4.2),
        };
    }
    let ScrewConnectionHardwareSource::Generated { part } = &mut request.hardware[0].source else {
        unreachable!()
    };
    let MechanicalPartParameters::CustomConnection { parameters } = &mut part.parameters else {
        unreachable!()
    };
    let CustomConnectionHardwareParameters::WoodScrew { wood_thread, .. } = parameters.as_mut()
    else {
        unreachable!()
    };
    wood_thread.hand = if left {
        ThreadHand::Left
    } else {
        ThreadHand::Right
    };
    wood_thread.starts = starts;
    (Session::new(session.document().clone()).unwrap(), request)
}
#[test]
fn actual_direct_wood_contact_survives_both_hands_starts_pose_edit_history_and_reopen() {
    for (posed, left, starts) in [(false, false, 1), (true, true, 1)] {
        let (mut session, mut request) = fixture(posed, left, starts);
        let before = session.document().clone();
        let prepared = prepare_screw_connection(&session, 0, actor(), &request).unwrap();
        let candidate = compiler()
            .compile_with_backend(prepared.prepared.document(), 1, GeometryBackend::Rust, None)
            .unwrap();
        let check = &candidate.snapshot().screw_connection_checks[0];
        assert_eq!(check.mode, ScrewConnectionMode::DirectWood);
        assert!(check.native_component_collision_checked);
        assert!((check.actual_wood_engagement_mm.unwrap() - 12.).abs() < 1e-7);
        let path = OccurrencePath {
            root: request.hardware[0].instance,
            members: vec![],
        };
        let screw = candidate.leaf_solid(&path).unwrap();
        let (source, frame) = screw.rust_machine_screw_geometry().unwrap();
        assert_eq!(source.thread().spec().starts, starts);
        assert_eq!(source.thread().spec().left_hand, left);
        assert_eq!(source.thread().spec().length_mm, 12.);
        for solid in [
            screw,
            candidate
                .leaf_solid(&request.receivers[0].source.occurrence)
                .unwrap(),
        ] {
            let mesh = solid
                .manufacturing_mesh_with_receipt(Default::default())
                .unwrap();
            assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("actual-leaf.step");
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
        session.undo(1, actor()).unwrap();
        assert_eq!(*session.document(), before);
        session.redo(2, actor()).unwrap();
        assert_eq!(*session.document(), created);
        let ConnectionReceiverExtent::Blind { depth, .. } = &mut request.receivers[0].source.extent
        else {
            unreachable!()
        };
        *depth = length(15.75);
        let edit =
            prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
        compiler()
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
        let model = compiler().compile(&reopened, 0).unwrap();
        assert_eq!(
            model
                .leaf_solid(&path)
                .unwrap()
                .rust_machine_screw_geometry()
                .unwrap()
                .1,
            frame
        );
        assert_eq!(
            read_screw_connection_request(&reopened, request.connection).unwrap(),
            request
        );
        assert_eq!(*session.document(), reopened);
    }
}

#[test]
fn foreign_unknown_source_cannot_be_declared_clear_by_an_old_bounding_box() {
    let (mut session, request) = fixture(false, false, 1);
    let unknown = PrimitiveRequest::new(
        "Unqualified curved counterpart",
        PrimitiveParameters::Cylinder {
            radius: length(1.),
            height: length(2.),
            frame: PrimitiveFrame {
                origin: std::array::from_fn(|_| length(100.)),
                ..Default::default()
            },
        },
    );
    session
        .commit(
            prepare_primitive(&session, 0, actor(), &unknown)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let before = session.document().clone();
    let candidate =
        prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
    let failure = compiler()
        .compile(candidate.prepared.document(), session.revision() + 1)
        .err()
        .unwrap();
    assert!(failure.message.contains("source separation"), "{failure:?}");
    assert_eq!(*session.document(), before);
}

#[test]
fn a_far_actual_rectangular_counterpart_has_an_owned_separation_proof() {
    let (mut session, request) = fixture(false, false, 1);
    let known = PrimitiveRequest::new(
        "Known rectangular counterpart",
        PrimitiveParameters::Box {
            size: std::array::from_fn(|_| length(2.)),
            frame: PrimitiveFrame {
                origin: std::array::from_fn(|_| length(100.)),
                ..Default::default()
            },
        },
    );
    session
        .commit(
            prepare_primitive(&session, 0, actor(), &known)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let before = session.document().clone();
    let candidate =
        prepare_screw_connection(&session, session.revision(), actor(), &request).unwrap();
    let model = compiler()
        .compile(candidate.prepared.document(), session.revision() + 1)
        .unwrap();
    assert!(model.snapshot().screw_connection_checks[0].native_component_collision_checked);
    assert_eq!(*session.document(), before);
}
