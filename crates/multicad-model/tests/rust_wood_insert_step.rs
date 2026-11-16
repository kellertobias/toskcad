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

#[test]
fn composite_entry_boundary_receipts_never_invent_decoded_mass() {
    use multicad_core::document::EntityId;
    use multicad_kernel::{GeometryPolicy, HelicalThread, Mesh, PlaneFrame, Solid};
    use multicad_model::ShapeSnapshot;
    use std::sync::Arc;
    let policy = GeometryPolicy::default();
    let root = 1. / (4. * 3_f64.sqrt());
    let spec = HelicalThread {
        major_radius_mm: 3.,
        pitch_mm: 1.,
        length_mm: 6.6,
        depth_mm: 0.875 * 3_f64.sqrt() / 2. - root,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand: true,
    };
    let entry = Solid::rust_helical_thread_entry_chamfer(
        spec,
        Some(root),
        0.3,
        None,
        EntityId::new(),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let sphere = Solid::rust_sphere_named(
        1.,
        PlaneFrame {
            origin: [20., 0., 0.],
            ..Default::default()
        },
        EntityId::new(),
        policy,
    )
    .unwrap();
    // This is an actual Native source snapshot consumed through the public
    // Model receipt API. Empty display arrays match omit_display compilation.
    let owner = EntityId::new();
    let solid = Solid::compound_named(&[&entry, &sphere], owner, policy).unwrap();
    assert!(!solid.has_wood_insert_step_source());
    assert!(solid.has_partitioned_step_source());
    let shape = ShapeSnapshot {
        producer: owner,
        bounds: solid.bounds().unwrap(),
        volume_mm3: solid.volume_mm3().unwrap(),
        solid_count: solid.solid_count().unwrap(),
        mesh: Arc::new(Mesh {
            vertices: vec![],
            indices: vec![],
            face_indices: vec![],
        }),
        topology: solid.topology_snapshot(owner, 0).unwrap(),
        edges: Arc::new(vec![]),
        vertices: Arc::new(vec![]),
    };
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("entry-composite.step");
    solid.write_step(&path).unwrap();
    let receipt = verify_rust_step_shape_receipt(&path, &solid, &shape, 1e-7).unwrap();
    assert_eq!(receipt.method, StepVerificationMethod::MixedMemberBoundary);
    assert!(receipt.decoded_volume_mm3.is_none());
    assert!(receipt.retained_members.is_empty());
    assert_eq!(receipt.decoded_scalar_members.len(), 1);
    assert_eq!(receipt.decoded_scalar_members[0].0, 1);
    assert_eq!(receipt.source_known_members.len(), 1);
    assert_eq!(receipt.source_known_members[0].member_index, 0);
    assert!(
        receipt.source_known_members[0].source_volume_bounds_mm3[0] <= entry.volume_mm3().unwrap()
    );
    assert!(
        receipt.source_known_members[0].source_volume_bounds_mm3[1] >= entry.volume_mm3().unwrap()
    );
    assert!(multicad_model::export::verify_rust_step_shape(&path, &solid, &shape, 1e-7).is_err());
    let encoded = serde_json::to_value(&receipt).unwrap();
    assert!(
        encoded["source_known_members"][0]
            .get("decoded_volume_mm3")
            .is_none()
    );
    let reopened: multicad_model::export::StepVerificationReceipt =
        serde_json::from_value(encoded).unwrap();
    assert_eq!(reopened, receipt);
}

#[test]
fn adjacent_binary_length_insert_projects_use_the_normal_export_and_reopen_contract() {
    use multicad_core::project_file::{self, ProjectFile};
    for (extent, left, posed) in [
        (4_f64.next_up(), false, false),
        (4_f64.next_up(), true, true),
        (4_f64.next_down(), false, true),
        (4.5, true, false),
    ] {
        let mut session =
            Session::new(Document::new("Exact finite insert export").unwrap()).unwrap();
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
        outer_thread.length = length(extent);
        inner_thread.length = length(extent);
        outer_thread.hand = if left {
            ThreadHand::Left
        } else {
            ThreadHand::Right
        };
        inner_thread.hand = outer_thread.hand;
        if posed {
            *frame = PrimitiveFrame {
                origin: [40., -18., 7.].map(length),
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            };
        }
        let request = MechanicalPartRequest::new("Actual finite insert", parameters);
        session
            .commit(
                prepare_mechanical_part(
                    &session,
                    session.revision(),
                    Actor::local(Channel::Mcp),
                    &request,
                )
                .unwrap()
                .prepared,
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
        let reopened = project_file::read(
            &mut bytes.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap()
        .project;
        assert_eq!(reopened.document, before);
        let model = Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile_with_backend(&reopened.document, 0, GeometryBackend::Rust, None)
        .unwrap();
        let solid = model.solid(request.feature).unwrap();
        let data = solid.rust_wood_insert_data().unwrap();
        assert_eq!(
            data.source.outer().full_length_mm().to_bits(),
            extent.to_bits()
        );
        assert_eq!(
            data.source.inner().full_length_mm().to_bits(),
            extent.to_bits()
        );
        assert_eq!(data.source.outer().spec().left_hand, left);
        assert_eq!(data.source.inner().spec().left_hand, left);
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("exact-insert.step");
        solid.write_step(&path).unwrap();
        let shape = &model.snapshot().shapes[&request.feature];
        let proof = verify_rust_step_shape_receipt(&path, solid, shape, 1e-7).unwrap();
        assert_eq!(
            proof.method,
            StepVerificationMethod::DecodedBoundaryEnclosure
        );
        assert_eq!(proof.retained_members.len(), 1);
        assert_eq!(proof.source_correspondence.len(), 1);
        assert!(proof.source_correspondence[0].maximum_position_error_mm <= 1e-7);
        assert!(proof.source_correspondence[0].row_derivative_error_mm <= 1e-6);
        let member = &proof.retained_members[0];
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
}
