use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::*,
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, GeometryBackend};
#[test]
fn canonical_default_insert_compiles_visible_source_geometry_and_reopens_without_native_fallback() {
    let mut session = Session::new(Document::new("Actual default wood insert").unwrap()).unwrap();
    let request =
        MechanicalPartRequest::new("Female insert", default_custom_connection_wood_insert());
    session
        .commit(
            prepare_mechanical_part(&session, 0, Actor::local(Channel::Mcp), &request)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let before = session.document().clone();
    let compiler = Compiler::default();
    let timer = std::time::Instant::now();
    let model = compiler
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let shape = model.solid(request.feature).unwrap();
    let data = shape.rust_wood_insert_data().unwrap();
    assert_eq!(data.source.outer().spec().major_radius_mm, 5.);
    assert_eq!(data.source.inner().spec().major_radius_mm, 1.5);
    assert_eq!(data.source.full_length_mm(), 4.);
    assert_eq!(shape.solid_count().unwrap(), 1);
    assert!(shape.volume_mm3().unwrap() > 0.);
    assert_eq!(shape.faces().unwrap().len(), 10);
    let receipt = shape
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(before.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let reopened = project_file::read(
        &mut std::io::Cursor::new(bytes),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let fresh = compiler
        .compile_with_backend(&reopened.project.document, 0, GeometryBackend::Rust, None)
        .unwrap();
    assert_eq!(
        fresh
            .solid(request.feature)
            .unwrap()
            .rust_wood_insert_data(),
        Some(data)
    );
    assert_eq!(*session.document(), before);
    eprintln!(
        "actual canonical insert visible compile+MFG+reopen {:?}",
        timer.elapsed()
    );
}

#[test]
fn actual_insert_edit_history_and_saved_pose_keep_both_original_thread_sources() {
    let actor = || Actor::local(Channel::Mcp);
    let mut s = Session::new(Document::new("Insert history").unwrap()).unwrap();
    let mut r =
        MechanicalPartRequest::new("Actual insert", default_custom_connection_wood_insert());
    s.commit(
        prepare_mechanical_part(&s, 0, actor(), &r)
            .unwrap()
            .prepared,
    )
    .unwrap();
    let compiler = Compiler::default();
    let original = compiler
        .compile_with_backend(s.document(), s.revision(), GeometryBackend::Rust, None)
        .unwrap()
        .solid(r.feature)
        .unwrap()
        .rust_wood_insert_data()
        .unwrap();
    let MechanicalPartParameters::CustomConnection { parameters } = &mut r.parameters else {
        panic!()
    };
    let CustomConnectionHardwareParameters::WoodInsert {
        outer_thread,
        inner_thread,
        frame,
    } = parameters.as_mut()
    else {
        panic!()
    };
    outer_thread.length = length(4.3);
    inner_thread.length = length(4.3);
    outer_thread.hand = ThreadHand::Left;
    inner_thread.hand = ThreadHand::Left;
    *frame = PrimitiveFrame {
        origin: [length(11.), length(-7.), length(5.)],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    s.commit(
        prepare_mechanical_part_edit(&s, s.revision(), actor(), r.feature, &r.parameters)
            .unwrap()
            .prepared,
    )
    .unwrap();
    let model = compiler
        .compile_with_backend(s.document(), s.revision(), GeometryBackend::Rust, None)
        .unwrap();
    let shape = model.solid(r.feature).unwrap();
    let changed = shape.rust_wood_insert_data().unwrap();
    assert_ne!(changed, original);
    assert!(changed.source.outer().spec().left_hand && changed.source.inner().spec().left_hand);
    assert_eq!(changed.source.full_length_mm(), 4.3);
    let mesh = shape
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert!(mesh.edge_pairs > 0 && mesh.vertex_links > 0);
    s.undo(s.revision(), actor()).unwrap();
    assert_eq!(
        compiler
            .compile_with_backend(s.document(), s.revision(), GeometryBackend::Rust, None)
            .unwrap()
            .solid(r.feature)
            .unwrap()
            .rust_wood_insert_data(),
        Some(original)
    );
    s.redo(s.revision(), actor()).unwrap();
    assert_eq!(
        compiler
            .compile_with_backend(s.document(), s.revision(), GeometryBackend::Rust, None)
            .unwrap()
            .solid(r.feature)
            .unwrap()
            .rust_wood_insert_data(),
        Some(changed)
    );
    let before = s.document().clone();
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(before.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let opened = project_file::read(
        &mut std::io::Cursor::new(bytes),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    assert_eq!(
        compiler
            .compile_with_backend(&opened.project.document, 0, GeometryBackend::Rust, None)
            .unwrap()
            .solid(r.feature)
            .unwrap()
            .rust_wood_insert_data(),
        Some(changed)
    );
    assert_eq!(*s.document(), before);
}
