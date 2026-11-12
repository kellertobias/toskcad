use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{self, *},
};
use multicad_model::Compiler;
#[test]
fn metric_receiving_terminal_relief_uses_the_actual_neck_source() {
    let actor = Actor::local(Channel::Ui);
    let mut s = Session::new(Document::new("M3 receiving terminal relief").unwrap()).unwrap();
    let stock = PrimitiveRequest::new(
        "Blank",
        PrimitiveParameters::Cylinder {
            radius: length(1.5),
            height: length(3.3),
            frame: Default::default(),
        },
    );
    s.commit(
        model::prepare_primitive(&s, 0, actor.clone(), &stock)
            .unwrap()
            .prepared,
    )
    .unwrap();
    let request = ThreadRequest::receiving(
        "Relieved M3",
        ThreadMode::CutExternal,
        BodyInput {
            body: stock.body,
            producer: stock.feature,
        },
        ThreadParameters {
            major_diameter: length(3.),
            pitch: length(0.5),
            length: length(3.3),
            fit: ThreadFit::Metric6gExternal,
            end_treatments: Some(ThreadEndTreatments {
                start: ThreadStartTreatment::FlatCut,
                end: ThreadEndTreatment::CylindricalRelief {
                    length: length(0.3),
                    radial_depth: length(0.05),
                },
            }),
            ..Default::default()
        },
    )
    .unwrap();
    let q = model::prepare_thread(&s, s.revision(), actor, &request).unwrap();
    let result = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile(q.prepared.document(), s.revision() + 1)
    .unwrap();
    assert!(
        result
            .solid(request.feature)
            .unwrap()
            .rust_thread_end_relief_data()
            .is_some()
    );
}

#[test]
fn canonical_terminal_receiver_edits_reopen_export_and_refusals_are_atomic() {
    use multicad_core::{
        document::EntityId,
        exchange::{Scope, ScopeRoot, SetDefinition, prepare_set},
        project_file::{self, ProjectFile},
        references::{GeometryIntent, Resolution, TopologyKind, resolve},
    };
    use multicad_model::export::{
        Exporter, StepVerificationMethod, verify_rust_step_shape_receipt,
    };
    let actor = Actor::local(Channel::Ui);
    let compiler = Compiler {
        omit_display: true,
        ..Default::default()
    };
    for left in [false, true] {
        let mut s = Session::new(Document::new("Receiving terminal history").unwrap()).unwrap();
        let stock = PrimitiveRequest::new(
            "M3 blank",
            PrimitiveParameters::Cylinder {
                radius: length(1.5),
                height: length(3.3),
                frame: PrimitiveFrame {
                    origin: [length(100.), length(-30.), length(50.)],
                    u: [0.8, 0.6, 0.],
                    v: [-0.36, 0.48, 0.8],
                },
            },
        );
        s.commit(
            model::prepare_primitive(&s, 0, actor.clone(), &stock)
                .unwrap()
                .prepared,
        )
        .unwrap();
        let stock_model = compiler.compile(s.document(), s.revision()).unwrap();
        let old = &stock_model.snapshot().shapes[&stock.feature].topology;
        let intents = old
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face && e.key.ordinal < 2)
            .map(|e| {
                let mut i = GeometryIntent::capture(old, e.key).unwrap();
                i.boundary = None;
                i
            })
            .collect::<Vec<_>>();
        let parameters = ThreadParameters {
            major_diameter: length(3.),
            pitch: length(0.5),
            length: length(3.3),
            fit: ThreadFit::Metric6gExternal,
            hand: if left {
                ThreadHand::Left
            } else {
                ThreadHand::Right
            },
            frame: PrimitiveFrame {
                origin: [length(100.), length(-30.), length(50.)],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
            end_treatments: Some(ThreadEndTreatments {
                start: ThreadStartTreatment::FlatCut,
                end: ThreadEndTreatment::CylindricalRelief {
                    length: length(0.2),
                    radial_depth: length(0.05),
                },
            }),
            ..Default::default()
        };
        let request = ThreadRequest::receiving(
            "M3 terminal",
            ThreadMode::CutExternal,
            BodyInput {
                body: stock.body,
                producer: stock.feature,
            },
            parameters.clone(),
        )
        .unwrap();
        let q = model::prepare_thread(&s, s.revision(), actor.clone(), &request).unwrap();
        let first = compiler
            .compile(q.prepared.document(), s.revision() + 1)
            .unwrap();
        assert!(
            first
                .solid(request.feature)
                .unwrap()
                .rust_thread_end_relief_data()
                .is_some()
        );
        s.commit(q.prepared).unwrap();
        let before = s.document().clone();
        let edited = ThreadParameters {
            end_treatments: Some(ThreadEndTreatments {
                start: ThreadStartTreatment::FlatCut,
                end: ThreadEndTreatment::CylindricalRelief {
                    length: length(0.3),
                    radial_depth: length(0.07),
                },
            }),
            ..parameters
        };
        let q =
            model::prepare_thread_edit(&s, s.revision(), actor.clone(), request.feature, &edited)
                .unwrap();
        let model = compiler
            .compile(q.prepared.document(), s.revision() + 1)
            .unwrap();
        let solid = model.solid(request.feature).unwrap();
        for intent in &intents {
            assert!(matches!(
                resolve(
                    intent,
                    &model.snapshot().shapes[&request.feature].topology,
                    s.revision() + 1
                ),
                Resolution::Resolved { .. }
            ));
        }
        assert!(
            solid.volume_mm3().unwrap()
                < first.solid(request.feature).unwrap().volume_mm3().unwrap()
        );
        assert_eq!(s.document(), &before);
        s.commit(q.prepared).unwrap();
        let saved = s.document().clone();
        s.undo(s.revision(), actor.clone()).unwrap();
        assert_eq!(s.document(), &before);
        s.redo(s.revision(), actor.clone()).unwrap();
        assert_eq!(s.document(), &saved);
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
        let reopened_model = compiler.compile(&reopened.document, s.revision()).unwrap();
        assert_eq!(
            solid.rust_thread_end_relief_data(),
            reopened_model
                .solid(request.feature)
                .unwrap()
                .rust_thread_end_relief_data()
        );
        for invalid in [
            ThreadParameters {
                length: length(3.),
                ..edited.clone()
            },
            ThreadParameters {
                length: length(3.6),
                ..edited.clone()
            },
        ] {
            let q = model::prepare_thread_edit(
                &s,
                s.revision(),
                actor.clone(),
                request.feature,
                &invalid,
            )
            .unwrap();
            assert!(
                compiler
                    .compile(q.prepared.document(), s.revision() + 1)
                    .is_err()
            );
            assert_eq!(s.document(), &saved);
        }
        let set = EntityId::new();
        s.commit(
            prepare_set(
                &s,
                s.revision(),
                actor.clone(),
                set,
                "Receiving STEP",
                &SetDefinition {
                    scope: Scope {
                        root: ScopeRoot::Explicit {
                            entities: vec![stock.body],
                        },
                        tags: None,
                    },
                    ..Default::default()
                },
            )
            .unwrap(),
        )
        .unwrap();
        let exporter = Exporter::default();
        let preview = exporter
            .preview(s.document(), s.revision(), set, &actor, None)
            .unwrap();
        assert!(preview.can_run(), "{preview:?}");
        let dir = tempfile::tempdir().unwrap();
        let result = exporter
            .run(
                s.document(),
                s.revision(),
                &preview.key,
                dir.path(),
                &actor,
                None,
            )
            .unwrap();
        assert_eq!(result.files.len(), 1);
        let file =
            std::path::Path::new(&result.directory).join(&result.files[0].artifact.output.filename);
        let receipt = verify_rust_step_shape_receipt(
            &file,
            solid,
            &model.snapshot().shapes[&request.feature],
            1e-7,
        )
        .unwrap();
        assert_eq!(receipt.method, StepVerificationMethod::DecodedBoundary);
        assert!((receipt.decoded_volume_mm3.unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-5);
        let imported = multicad_kernel::Solid::read_step(&file, Default::default()).unwrap();
        assert!(
            imported
                .rust_thread_end_relief_data()
                .unwrap()
                .0
                .intent()
                .is_none()
        );
        assert!((imported.volume_mm3().unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-5);
        assert!(
            compiler
                .compile_with_backend(
                    s.document(),
                    s.revision(),
                    multicad_model::GeometryBackend::Rust,
                    Some(&|| true)
                )
                .is_err()
        );
    }
}
