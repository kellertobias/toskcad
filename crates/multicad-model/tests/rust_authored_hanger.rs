use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::*,
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, GeometryBackend};
fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}
fn create() -> (Session, MechanicalPartRequest) {
    let mut session = Session::new(Document::new("Actual two-ended hanger").unwrap()).unwrap();
    let request = MechanicalPartRequest::new(
        "Wood and metric hanger",
        default_custom_connection_hanger_bolt(),
    );
    session
        .commit(
            prepare_mechanical_part(&session, 0, actor(), &request)
                .unwrap()
                .prepared,
        )
        .unwrap();
    (session, request)
}
#[test]
fn canonical_hanger_compiles_actual_visible_geometry_and_manufacturing_without_occt() {
    let (session, request) = create();
    let before = session.document().clone();
    let started = std::time::Instant::now();
    let model = Compiler::default()
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let shape = model.solid(request.feature).unwrap();
    let data = shape.rust_hanger_data().unwrap();
    assert_eq!(data.source.wood().full_length_mm(), 12.);
    assert_eq!(data.source.metric().full_length_mm(), 12.);
    assert_eq!(data.source.spec().plain_length_mm, 2.);
    assert!(data.source.metric().spec().major_radius_mm < 1.5);
    assert_eq!(shape.faces().unwrap().len(), 13);
    assert_eq!(
        shape
            .topology_snapshot(request.feature, session.revision())
            .unwrap()
            .elements
            .len(),
        57
    );
    let manufacturing = shape
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert!(manufacturing.edge_pairs > 0 && manufacturing.vertex_links > 0);
    assert!(shape.volume_mm3().unwrap() > 0.);
    assert_eq!(*session.document(), before);
    eprintln!(
        "canonical hanger visible compile and actual MFG {:?}",
        started.elapsed()
    );
}
#[test]
fn hanger_hand_lengths_pose_history_and_saved_source_reopen_preserve_real_components() {
    let (mut session, mut request) = create();
    let compiler = Compiler::default();
    let original = compiler
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap()
        .solid(request.feature)
        .unwrap()
        .rust_hanger_data()
        .unwrap();
    let MechanicalPartParameters::CustomConnection { parameters } = &mut request.parameters else {
        panic!()
    };
    let CustomConnectionHardwareParameters::HangerBolt {
        wood_thread,
        metric_thread,
        plain_length,
        frame,
        ..
    } = parameters.as_mut()
    else {
        panic!()
    };
    wood_thread.length = length(12.3);
    metric_thread.length = length(12.7);
    wood_thread.hand = ThreadHand::Left;
    metric_thread.hand = ThreadHand::Left;
    *plain_length = length(2.3);
    *frame = PrimitiveFrame {
        origin: [length(11.), length(-7.), length(5.)],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    session
        .commit(
            prepare_mechanical_part_edit(
                &session,
                session.revision(),
                actor(),
                request.feature,
                &request.parameters,
            )
            .unwrap()
            .prepared,
        )
        .unwrap();
    let model = compiler
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let shape = model.solid(request.feature).unwrap();
    let edited = shape.rust_hanger_data().unwrap();
    assert_ne!(edited, original);
    assert!(edited.source.metric().spec().left_hand && edited.source.wood().spec().left_hand);
    assert_eq!(edited.source.metric().full_length_mm(), 12.7);
    assert_eq!(edited.source.wood().full_length_mm(), 12.3);
    assert_eq!(edited.source.spec().plain_length_mm, 2.3);
    shape
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(
        compiler
            .compile_with_backend(
                session.document(),
                session.revision(),
                GeometryBackend::Rust,
                None
            )
            .unwrap()
            .solid(request.feature)
            .unwrap()
            .rust_hanger_data(),
        Some(original)
    );
    session.redo(session.revision(), actor()).unwrap();
    assert_eq!(
        compiler
            .compile_with_backend(
                session.document(),
                session.revision(),
                GeometryBackend::Rust,
                None
            )
            .unwrap()
            .solid(request.feature)
            .unwrap()
            .rust_hanger_data(),
        Some(edited)
    );
    let before = session.document().clone();
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(before.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    assert_eq!(
        compiler
            .compile_with_backend(&opened.project.document, 0, GeometryBackend::Rust, None)
            .unwrap()
            .solid(request.feature)
            .unwrap()
            .rust_hanger_data(),
        Some(edited)
    );
    assert_eq!(*session.document(), before);
}
#[test]
fn tiny_positive_hanger_shoulder_refuses_before_current_document_or_body_changes() {
    let (session, mut request) = create();
    let compiler = Compiler::default();
    let current = compiler
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let before = session.document().clone();
    let actual = current
        .solid(request.feature)
        .unwrap()
        .rust_hanger_data()
        .unwrap();
    let MechanicalPartParameters::CustomConnection { parameters } = &mut request.parameters else {
        panic!()
    };
    let CustomConnectionHardwareParameters::HangerBolt { plain_diameter, .. } = parameters.as_mut()
    else {
        panic!()
    };
    *plain_diameter = length(
        2. * (actual.source.metric().maximum_radius_mm()
            + multicad_kernel::GeometryPolicy::default().minimum_feature_mm / 10.),
    );
    let pending = prepare_mechanical_part_edit(
        &session,
        session.revision(),
        actor(),
        request.feature,
        &request.parameters,
    )
    .unwrap();
    assert!(
        compiler
            .compile_with_backend(
                pending.prepared.document(),
                session.revision() + 1,
                GeometryBackend::Rust,
                None
            )
            .is_err()
    );
    assert_eq!(*session.document(), before);
    assert_eq!(
        current.solid(request.feature).unwrap().rust_hanger_data(),
        Some(actual)
    );
}
