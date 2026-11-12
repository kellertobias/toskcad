use multicad_core::{
    commands::{Actor, Channel, Session},
    constraint_sketch::Dimension,
    document::Document,
    model::{
        self, CustomThreadProfile, ThreadFamily, ThreadHand, ThreadParameters, ThreadRequest,
        angle, length,
    },
    project_file::{self, ProjectFile},
    references::{GeometryIntent, Resolution, TopologyKind, resolve},
};
use multicad_kernel::MeshOptions;
use multicad_model::{
    Compiler,
    export::{StepVerificationMethod, verify_rust_step_shape_receipt},
};

#[test]
fn custom_taper_normal_compiler_preserves_law_references_history_and_export() {
    let actor = Actor::local(Channel::Ui);
    let mut session = Session::new(Document::new("Owned custom taper").unwrap()).unwrap();
    let parameters = ThreadParameters {
        family: ThreadFamily::Custom,
        custom: Some(CustomThreadProfile {
            depth: length(0.5),
            crest_width: length(0.25),
            lower_flank: angle(12.),
            upper_flank: angle(18.),
            radial_taper: Dimension::Literal {
                value: 0.03125,
                unit: "ratio".into(),
            },
        }),
        length: length(1.25),
        starts: 2,
        ..Default::default()
    };
    let request = ThreadRequest::new("Two-start cone-spine", parameters.clone());
    let compiler = Compiler {
        omit_display: true,
        ..Default::default()
    };
    let prepared = model::prepare_thread(&session, 0, actor.clone(), &request).unwrap();
    let preview = compiler.compile(prepared.prepared.document(), 1).unwrap();
    let solid = preview.solid(request.feature).unwrap();
    let (source, _) = solid.rust_retained_helical_thread_source().unwrap();
    assert_eq!(source.profile.starts, 2);
    assert_eq!(source.profile.radial_taper, 0.03125);
    let shape = &preview.snapshot().shapes[&request.feature];
    let intents = [
        (TopologyKind::Face, "cap.end"),
        (TopologyKind::Edge, "thread/start/0/crest/leading_rail"),
    ]
    .map(|(kind, role)| {
        GeometryIntent::capture(
            &shape.topology,
            shape
                .topology
                .elements
                .iter()
                .find(|e| e.key.kind == kind && e.origins.iter().any(|o| o.role == role))
                .unwrap()
                .key,
        )
        .unwrap()
    });
    session.commit(prepared.prepared).unwrap();
    let original = session.document().clone();
    let edited = ThreadParameters {
        length: length(1.5),
        hand: ThreadHand::Left,
        starts: 1,
        ..parameters
    };
    let prepared =
        model::prepare_thread_edit(&session, 1, actor.clone(), request.feature, &edited).unwrap();
    let model = compiler.compile(prepared.prepared.document(), 2).unwrap();
    let shape = &model.snapshot().shapes[&request.feature];
    for intent in &intents {
        assert!(matches!(
            resolve(intent, &shape.topology, 2),
            Resolution::Resolved { .. }
        ));
    }
    let solid = model.solid(request.feature).unwrap();
    let data = solid.rust_retained_helical_thread_source().unwrap().0;
    assert!(data.profile.left_hand);
    assert_eq!(data.profile.starts, 1);
    assert_eq!(data.profile.length_mm, 1.5);
    let receipt = solid
        .manufacturing_mesh_with_receipt(MeshOptions {
            linear_deflection_mm: 0.05,
            angular_deflection_rad: 0.4,
            maximum_triangles: 200000,
        })
        .unwrap();
    assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
    // A display/manufacturing budget refusal never changes canonical history.
    assert!(
        solid
            .manufacturing_mesh(MeshOptions {
                maximum_triangles: 1,
                ..MeshOptions {
                    linear_deflection_mm: 0.05,
                    angular_deflection_rad: 0.4,
                    maximum_triangles: 200000
                }
            })
            .is_err()
    );
    assert_eq!(session.document(), &original);
    let directory = tempfile::tempdir().unwrap();
    let file = directory.path().join("custom.step");
    solid.write_step(&file).unwrap();
    let receipt = verify_rust_step_shape_receipt(&file, solid, shape, 1e-7).unwrap();
    assert_eq!(
        receipt.method,
        StepVerificationMethod::DecodedBoundaryEnclosure
    );
    assert!(receipt.decoded_volume_mm3.is_none());
    assert_eq!(receipt.retained_members.len(), 1);
    session.commit(prepared.prepared).unwrap();
    let saved = session.document().clone();
    session.undo(2, actor.clone()).unwrap();
    assert_eq!(session.document(), &original);
    session.redo(3, actor.clone()).unwrap();
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
    assert_eq!(
        model::read_thread(&reopened.document, request.feature)
            .unwrap()
            .parameters,
        edited
    );
    let rebuilt = compiler.compile(&reopened.document, 0).unwrap();
    assert_eq!(
        rebuilt.snapshot().shapes[&request.feature].volume_mm3,
        shape.volume_mm3
    );
    let viewport = Compiler::default().compile(&saved, 4).unwrap();
    let visible = &viewport.snapshot().shapes[&request.feature];
    assert!(!visible.mesh.indices.is_empty());
    assert_eq!(visible.edges.len(), 14);
    assert_eq!(visible.vertices.len(), 10);
    let canceled = compiler
        .compile_with_backend(
            &saved,
            4,
            multicad_model::GeometryBackend::Rust,
            Some(&|| true),
        )
        .err()
        .unwrap();
    assert!(canceled.message.contains("superseded"));
    assert_eq!(session.document(), &saved);
    let unsupported = ThreadParameters {
        end_treatments: Some(multicad_core::model::ThreadEndTreatments {
            start: multicad_core::model::ThreadStartTreatment::Chamfer45 {
                radial_depth: length(0.1),
            },
            end: multicad_core::model::ThreadEndTreatment::FlatCut,
        }),
        ..edited
    };
    assert!(
        model::prepare_thread_edit(&session, 4, actor.clone(), request.feature, &unsupported)
            .is_err()
    );
    let mut valid = model::read_thread(&saved, request.feature)
        .unwrap()
        .parameters;
    valid.length = length(1.6);
    let candidate =
        model::prepare_thread_edit(&session, 4, actor, request.feature, &valid).unwrap();
    let limited = Compiler {
        policy: multicad_kernel::GeometryPolicy {
            maximum_coordinate_mm: 1.,
            ..Default::default()
        },
        omit_display: true,
        ..Default::default()
    };
    assert!(limited.compile(candidate.prepared.document(), 5).is_err());
    let count = std::cell::Cell::new(0);
    let cancel = || {
        let n = count.get();
        count.set(n + 1);
        n >= 80
    };
    let failure = compiler
        .compile_with_backend(
            candidate.prepared.document(),
            5,
            multicad_model::GeometryBackend::Rust,
            Some(&cancel),
        )
        .err()
        .unwrap();
    assert!(failure.message.contains("superseded"));
    assert_eq!(session.document(), &saved);
}
