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
