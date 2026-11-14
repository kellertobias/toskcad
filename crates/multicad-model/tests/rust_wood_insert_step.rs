use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::*,
};
use multicad_model::{
    Compiler, GeometryBackend,
    export::{StepVerificationMethod, verify_rust_step_shape_receipt},
};

#[test]
fn canonical_insert_project_export_has_independently_decoded_file_mass() {
    let mut session = Session::new(Document::new("Actual insert export").unwrap()).unwrap();
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
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(
        session.document(),
        session.revision(),
        GeometryBackend::Rust,
        None,
    )
    .unwrap();
    let shape = &model.snapshot().shapes[&request.feature];
    let solid = model.solid(request.feature).unwrap();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("actual-insert.step");
    solid.write_step(&path).unwrap();
    let receipt = verify_rust_step_shape_receipt(&path, solid, shape, 1e-7)
        .expect("normal project export must verify the actual retained insert file");
    assert_eq!(
        receipt.method,
        StepVerificationMethod::DecodedBoundaryEnclosure
    );
    assert!(receipt.decoded_volume_mm3.is_none());
    assert_eq!(receipt.retained_members.len(), 1);
    assert_eq!(receipt.source_correspondence.len(), 1);
    assert!(receipt.source_correspondence[0].maximum_position_error_mm <= 1e-7);
    assert!(receipt.source_correspondence[0].row_derivative_error_mm <= 1e-6);
    let member = &receipt.retained_members[0];
    assert!(member.decoded_volume_bounds_mm3[0] > 0.);
    assert!(
        member.source_volume_bounds_mm3[0] <= shape.volume_mm3
            && shape.volume_mm3 <= member.source_volume_bounds_mm3[1]
    );
    let mut stale = shape.as_ref().clone();
    stale.volume_mm3 *= 1.01;
    assert!(verify_rust_step_shape_receipt(&path, solid, &stale, 1e-7).is_err());
    assert_eq!(session.document(), &before);
}
#[test]
fn canonical_insert_runs_the_complete_project_export_job() {
    use multicad_core::{
        document::EntityId,
        exchange::{Scope, ScopeRoot, SetDefinition, prepare_set},
    };
    use multicad_model::export::Exporter;
    let actor = Actor::local(Channel::Mcp);
    let mut s = Session::new(Document::new("Project insert export job").unwrap()).unwrap();
    let request =
        MechanicalPartRequest::new("Female insert", default_custom_connection_wood_insert());
    s.commit(
        prepare_mechanical_part(&s, 0, actor.clone(), &request)
            .unwrap()
            .prepared,
    )
    .unwrap();
    let set = EntityId::new();
    s.commit(
        prepare_set(
            &s,
            s.revision(),
            actor.clone(),
            set,
            "Actual insert STEP",
            &SetDefinition {
                scope: Scope {
                    root: ScopeRoot::Explicit {
                        entities: vec![request.body],
                    },
                    tags: None,
                },
                ..Default::default()
            },
        )
        .unwrap(),
    )
    .unwrap();
    let before = s.document().clone();
    let exporter = Exporter::default();
    let preview = exporter
        .preview(s.document(), s.revision(), set, &actor, None)
        .unwrap();
    assert!(preview.can_run(), "{preview:?}");
    let directory = tempfile::tempdir().unwrap();
    let result = exporter
        .run(
            s.document(),
            s.revision(),
            &preview.key,
            directory.path(),
            &actor,
            None,
        )
        .unwrap();
    assert_eq!(result.files.len(), 1);
    let path =
        std::path::Path::new(&result.directory).join(&result.files[0].artifact.output.filename);
    assert!(
        std::fs::read_to_string(path)
            .unwrap()
            .contains("MANIFOLD_SOLID_BREP")
    );
    assert_eq!(s.document(), &before);
}
