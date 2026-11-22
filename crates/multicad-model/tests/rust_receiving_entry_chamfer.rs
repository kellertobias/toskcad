use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{self, *},
};
use multicad_model::Compiler;
#[test]
fn metric_receiving_entry_uses_the_actual_supported_chamfer_source() {
    let actor = Actor::local(Channel::Ui);
    let mut s = Session::new(Document::new("M6 receiving entry").unwrap()).unwrap();
    let stock = PrimitiveRequest::new(
        "Shaft",
        PrimitiveParameters::Cylinder {
            radius: length(3.),
            height: length(6.6),
            frame: Default::default(),
        },
    );
    let q = model::prepare_primitive(&s, 0, actor.clone(), &stock).unwrap();
    s.commit(q.prepared).unwrap();
    let request = ThreadRequest::receiving(
        "M6 chamfered shaft",
        ThreadMode::CutExternal,
        BodyInput {
            body: stock.body,
            producer: stock.feature,
        },
        ThreadParameters {
            major_diameter: length(6.),
            pitch: length(1.),
            length: length(6.6),
            fit: ThreadFit::Metric6gExternal,
            end_treatments: Some(ThreadEndTreatments {
                start: ThreadStartTreatment::Chamfer45 {
                    radial_depth: length(0.3),
                },
                end: ThreadEndTreatment::FlatCut,
            }),
            ..Default::default()
        },
    )
    .unwrap();
    let q = model::prepare_thread(&s, s.revision(), actor, &request).unwrap();
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile(q.prepared.document(), s.revision() + 1)
    .unwrap();
    assert!(
        model
            .solid(request.feature)
            .unwrap()
            .rust_thread_entry_chamfer_data()
            .is_some()
    );
}

#[test]
fn canonical_receiving_entry_edits_reopen_export_and_refusals_are_atomic() {
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
        let mut s = Session::new(Document::new("Receiving entry history").unwrap()).unwrap();
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
                start: ThreadStartTreatment::Chamfer45 {
                    radial_depth: length(0.15),
                },
                end: ThreadEndTreatment::FlatCut,
            }),
            ..Default::default()
        };
        let request = ThreadRequest::receiving(
            "M3 entry",
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
                .rust_thread_entry_chamfer_data()
                .is_some()
        );
        s.commit(q.prepared).unwrap();
        let before = s.document().clone();
        let edited = ThreadParameters {
            end_treatments: Some(ThreadEndTreatments {
                start: ThreadStartTreatment::Chamfer45 {
                    radial_depth: length(0.2),
                },
                end: ThreadEndTreatment::CylindricalRelief {
                    length: length(0.3),
                    radial_depth: length(0.05),
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
            solid.rust_thread_entry_chamfer_data(),
            reopened_model
                .solid(request.feature)
                .unwrap()
                .rust_thread_entry_chamfer_data()
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
        assert_eq!(receipt.method, StepVerificationMethod::SourceKnownBoundary);
        assert_eq!(receipt.decoded_volume_mm3, None);
    }
}

#[test]
fn unsupported_receiving_families_and_internal_entry_do_not_change_the_document() {
    let actor = Actor::local(Channel::Ui);
    let mut s = Session::new(Document::new("Unsupported receiving entry").unwrap()).unwrap();
    let stock = PrimitiveRequest::new(
        "Blank",
        PrimitiveParameters::Cylinder {
            radius: length(5.),
            height: length(6.6),
            frame: Default::default(),
        },
    );
    s.commit(
        model::prepare_primitive(&s, 0, actor.clone(), &stock)
            .unwrap()
            .prepared,
    )
    .unwrap();
    let saved = s.document().clone();
    for (mode, family, fit) in [
        (
            ThreadMode::CutInternal,
            ThreadFamily::IsoMetric,
            ThreadFit::Metric6HInternal,
        ),
        (
            ThreadMode::CutExternal,
            ThreadFamily::Unc,
            ThreadFit::Unified2A,
        ),
    ] {
        let parameters = ThreadParameters {
            major_diameter: length(if family == ThreadFamily::Unc {
                6.35
            } else {
                6.
            }),
            pitch: length(if family == ThreadFamily::Unc {
                1.27
            } else {
                1.
            }),
            length: length(6.6),
            family,
            fit,
            end_treatments: Some(ThreadEndTreatments {
                start: ThreadStartTreatment::Chamfer45 {
                    radial_depth: length(0.2),
                },
                end: ThreadEndTreatment::FlatCut,
            }),
            ..Default::default()
        };
        let request = ThreadRequest::receiving(
            "Unsupported entry",
            mode,
            BodyInput {
                body: stock.body,
                producer: stock.feature,
            },
            parameters,
        )
        .unwrap();
        let candidate = model::prepare_thread(&s, s.revision(), actor.clone(), &request);
        if family == ThreadFamily::Unc {
            let error = candidate.err().unwrap();
            assert!(
                error
                    .to_string()
                    .contains("parallel single-start metric thread"),
                "{error}"
            );
            assert_eq!(s.document(), &saved);
            continue;
        }
        let q = candidate.unwrap();
        let error = Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile(q.prepared.document(), s.revision() + 1)
        .err()
        .unwrap();
        // A 0.2 mm chamfer on the M6x1 6H minor radius still starts inside
        // the thread cutter, so the receiving entry is refused for that reason.
        assert!(
            error.to_string().contains(
                "The receiving entry chamfer circle (radius 2.7175000000000002 mm) lies inside the thread cutter (maximum radius 3 mm)"
            ),
            "{error}"
        );
        assert_eq!(s.document(), &saved);
    }
}
