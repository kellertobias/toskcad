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

#[test]
fn saved_nonbinary_placed_insert_project_runs_normal_export_for_both_hands() {
    use multicad_core::{
        document::EntityId,
        exchange::{Scope, ScopeRoot, SetDefinition, prepare_set},
        project_file::{self, ProjectFile},
    };
    use multicad_model::export::Exporter;
    let actor = Actor::local(Channel::Mcp);
    let mut session = Session::new(Document::new("Placed insert project export").unwrap()).unwrap();
    let mut bodies = Vec::new();
    for (left, origin) in [(false, [11., -7., 5.]), (true, [-23., 17., 31.])] {
        let mut parameters = default_custom_connection_wood_insert();
        let MechanicalPartParameters::CustomConnection { parameters: custom } = &mut parameters
        else {
            panic!()
        };
        let CustomConnectionHardwareParameters::WoodInsert {
            outer_thread,
            inner_thread,
            frame,
        } = custom.as_mut()
        else {
            panic!()
        };
        outer_thread.length = length(4.3);
        inner_thread.length = length(4.3);
        outer_thread.hand = if left {
            ThreadHand::Left
        } else {
            ThreadHand::Right
        };
        inner_thread.hand = outer_thread.hand;
        *frame = PrimitiveFrame {
            origin: origin.map(length),
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        };
        let request = MechanicalPartRequest::new(
            if left { "Left insert" } else { "Right insert" },
            parameters,
        );
        session
            .commit(
                prepare_mechanical_part(&session, session.revision(), actor.clone(), &request)
                    .unwrap()
                    .prepared,
            )
            .unwrap();
        bodies.push(request.body);
    }
    let set = EntityId::new();
    session
        .commit(
            prepare_set(
                &session,
                session.revision(),
                actor.clone(),
                set,
                "Placed insert exports",
                &SetDefinition {
                    scope: Scope {
                        root: ScopeRoot::Explicit { entities: bodies },
                        tags: None,
                    },
                    ..Default::default()
                },
            )
            .unwrap(),
        )
        .unwrap();
    let before = session.document().clone();
    let mut bytes = Vec::new();
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
    assert_eq!(opened.project.document, before);
    let exporter = Exporter::default();
    let preview = exporter
        .preview(&opened.project.document, 0, set, &actor, None)
        .unwrap();
    assert!(preview.can_run(), "{preview:?}");
    let directory = tempfile::tempdir().unwrap();
    let result = exporter
        .run(
            &opened.project.document,
            0,
            &preview.key,
            directory.path(),
            &actor,
            None,
        )
        .unwrap();
    assert_eq!(result.files.len(), 2);
    for artifact in &result.files {
        let path = std::path::Path::new(&result.directory).join(&artifact.artifact.output.filename);
        let text = std::fs::read_to_string(path).unwrap();
        let file = multicad_kernel::rust_kernel::step_read::owned_insert::File::parse(
            &text,
            Default::default(),
            None,
        )
        .unwrap();
        let [member] = file.represented_members() else {
            panic!("exactly one actual exported occurrence")
        };
        assert!(member.mapping_link().is_some());
    }
    assert_eq!(opened.project.document, before);
    assert_eq!(session.document(), &before);
}
