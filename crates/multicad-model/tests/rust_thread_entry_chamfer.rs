use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{
        self, ThreadEndTreatment, ThreadEndTreatments, ThreadFit, ThreadParameters, ThreadRequest,
        ThreadStartTreatment, length,
    },
    project_file::{self, ProjectFile},
};
use multicad_model::{
    Compiler,
    export::{StepVerificationMethod, verify_rust_step_shape, verify_rust_step_shape_receipt},
};
#[test]
fn canonical_entry_and_terminal_edits_keep_literal_source_and_honest_export_mass() {
    let actor = Actor::local(Channel::Ui);
    let mut session = Session::new(Document::new("Metric conical entry").unwrap()).unwrap();
    let ends = |c, w| ThreadEndTreatments {
        start: ThreadStartTreatment::Chamfer45 {
            radial_depth: length(c),
        },
        end: if w > 0. {
            ThreadEndTreatment::CylindricalRelief {
                length: length(w),
                radial_depth: length(0.1),
            }
        } else {
            ThreadEndTreatment::FlatCut
        },
    };
    let parameters = ThreadParameters {
        major_diameter: length(6.),
        pitch: length(1.),
        length: length(6.6),
        fit: ThreadFit::Metric6gExternal,
        end_treatments: Some(ends(0.3, 0.)),
        ..Default::default()
    };
    let request = ThreadRequest::new("M6 conical entry", parameters.clone());
    let compiler = Compiler {
        omit_display: true,
        ..Default::default()
    };
    let prepared = model::prepare_thread(&session, 0, actor.clone(), &request).unwrap();
    let preview = compiler.compile(prepared.prepared.document(), 1).unwrap();
    assert!(
        preview
            .solid(request.feature)
            .unwrap()
            .rust_thread_entry_chamfer_data()
            .is_some()
    );
    session.commit(prepared.prepared).unwrap();
    let original = session.document().clone();
    let edited = ThreadParameters {
        end_treatments: Some(ends(0.4, 0.6)),
        ..parameters
    };
    let prepared =
        model::prepare_thread_edit(&session, 1, actor.clone(), request.feature, &edited).unwrap();
    let model = compiler.compile(prepared.prepared.document(), 2).unwrap();
    let solid = model.solid(request.feature).unwrap();
    let data = solid.rust_thread_entry_chamfer_data().unwrap();
    assert_eq!(data.source.source_thread().spec().length_mm, 6.6);
    assert_eq!(data.source.intent().unwrap().radial_depth_mm, 0.4);
    assert_eq!(
        data.source
            .terminal_relief()
            .unwrap()
            .intent()
            .unwrap()
            .width_mm,
        0.6
    );
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("entry.step");
    solid.write_step(&file).unwrap();
    let shape = &model.snapshot().shapes[&request.feature];
    let receipt = verify_rust_step_shape_receipt(&file, solid, shape, 1e-7).unwrap();
    assert_eq!(receipt.method, StepVerificationMethod::SourceKnownBoundary);
    assert_eq!(receipt.decoded_volume_mm3, None);
    let range = receipt.source_volume_bounds_mm3.unwrap();
    assert!(range[0] <= shape.volume_mm3 && shape.volume_mm3 <= range[1]);
    assert!(verify_rust_step_shape(&file, solid, shape, 1e-7).is_err());
    for tolerance in [f64::NAN, f64::INFINITY, 0., -1., 1e-30] {
        assert!(verify_rust_step_shape_receipt(&file, solid, shape, tolerance).is_err());
    }
    let mut stale = shape.as_ref().clone();
    stale.bounds.min[0] += 1.;
    assert!(verify_rust_step_shape_receipt(&file, solid, &stale, 1e-7).is_err());
    let text = std::fs::read_to_string(&file).unwrap();
    // The retained imported boundary reads the export back as the literal
    // boundary it states: no thread semantics are recovered, and its exact
    // volume, bounds and certified mesh agree with the thread body's
    // certified mass enclosure and exact bounds.
    let imported =
        multicad_kernel::Solid::import_rust_step(&text, request.feature, Default::default())
            .unwrap();
    assert!(imported.rust_thread_entry_chamfer_data().is_none());
    assert!(imported.rust_rounded_thread_data().is_none());
    assert!(matches!(
        imported.verify_imported_brep_step(&text),
        Some(Ok(_))
    ));
    // The retained volume is integrated with its own certified error: its
    // interval meets the thread body's mass enclosure.
    let retained =
        multicad_kernel::rust_kernel::step_read::brep::read_solids(&text, 1e-7, None).unwrap();
    let error: f64 = retained.iter().map(|s| s.volume_error_mm3()).sum();
    let volume = imported.volume_mm3().unwrap();
    assert!(error <= 1e-6 * volume, "{error}");
    assert!(
        volume + error >= range[0] && volume - error <= range[1],
        "{volume} {error} {range:?}"
    );
    let bounds = imported.bounds().unwrap();
    for axis in 0..3 {
        assert!((bounds.min[axis] - shape.bounds.min[axis]).abs() <= 1e-6);
        assert!((bounds.max[axis] - shape.bounds.max[axis]).abs() <= 1e-6);
    }
    let options = multicad_kernel::MeshOptions::default();
    let mesh = imported.manufacturing_mesh(options).unwrap();
    let corner = |index: u32| -> [f64; 3] {
        std::array::from_fn(|axis| mesh.vertices[3 * index as usize + axis])
    };
    let signed: f64 = mesh
        .indices
        .chunks(3)
        .map(|t| {
            let [a, b, c] = [t[0], t[1], t[2]].map(corner);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.
        })
        .sum();
    let area: f64 = imported.faces().unwrap().iter().map(|f| f.area_mm2).sum();
    let slack = options.linear_deflection_mm * area;
    assert!(
        range[0] - slack <= signed && signed <= range[1] + slack,
        "{signed} {range:?}"
    );
    for point in mesh.vertices.chunks(3) {
        for (axis, &x) in point.iter().enumerate() {
            assert!(x >= shape.bounds.min[axis] - 1e-6);
            assert!(x <= shape.bounds.max[axis] + 1e-6);
        }
    }
    std::fs::write(
        &file,
        text.replace("CONICAL_SURFACE", "CYLINDRICAL_SURFACE"),
    )
    .unwrap();
    assert!(verify_rust_step_shape_receipt(&file, solid, shape, 1e-7).is_err());
    session.commit(prepared.prepared).unwrap();
    let saved = session.document().clone();
    session.undo(2, actor.clone()).unwrap();
    assert_eq!(session.document(), &original);
    session.redo(3, actor).unwrap();
    assert_eq!(session.document(), &saved);
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(saved).unwrap(),
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
    let reopened = compiler.compile(&reopened.document, 0).unwrap();
    assert_eq!(
        reopened
            .solid(request.feature)
            .unwrap()
            .rust_thread_entry_chamfer_data()
            .unwrap()
            .source,
        data.source
    );
    // The ordinary background exporter must consume the discriminated proof,
    // rather than the decoded-only compatibility API above.
    {
        use multicad_core::document::EntityId;
        use multicad_core::exchange::{ExportResult, Scope, ScopeRoot, SetDefinition, prepare_set};
        use multicad_model::export::{ExportWorker, Exporter};
        let set = EntityId::new();
        let prepared = prepare_set(
            &session,
            session.revision(),
            Actor::local(Channel::Ui),
            set,
            "Entry STEP",
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
        .unwrap();
        session.commit(prepared).unwrap();
        let worker = ExportWorker::new(Exporter::default()).unwrap();
        let wait = |id| {
            let until = std::time::Instant::now() + std::time::Duration::from_secs(90);
            loop {
                if let Some(reply) = worker.poll() {
                    assert_eq!(reply.id, id);
                    break reply.result;
                }
                assert!(
                    std::time::Instant::now() < until,
                    "Entry export worker did not reply"
                );
                std::thread::sleep(std::time::Duration::from_millis(2));
            }
        };
        let preview = worker
            .preview(
                session.document(),
                session.revision(),
                set,
                Actor::local(Channel::Ui),
            )
            .unwrap();
        let ExportResult::Preview(preview) = wait(preview) else {
            panic!("Expected export preview")
        };
        let preview = preview.unwrap();
        assert!(preview.can_run(), "{preview:?}");
        let run = worker
            .run(
                session.document(),
                session.revision(),
                &preview.key,
                dir.path(),
                Actor::local(Channel::Ui),
            )
            .unwrap();
        let ExportResult::Run(run) = wait(run) else {
            panic!("Expected source-verified publication")
        };
        let run = run.unwrap();
        assert_eq!(run.files.len(), 1);
        assert_eq!(run.files[0].artifact.producer, request.feature);
        let actual =
            std::path::Path::new(&run.directory).join(&run.files[0].artifact.output.filename);
        assert!(verify_rust_step_shape_receipt(&actual, solid, shape, 1e-7).is_ok());
    }
}

#[test]
fn enabling_editing_and_disabling_entry_preserves_original_finite_far_rims_and_band_intents() {
    use multicad_core::{
        commands::Request,
        document::{Entity, EntityData, Link, Mutation},
        model::{PrimitiveFrame, ThreadHand},
        references::{
            AnalyticGeometry, Attachment, GeometryIntent, Resolution, TopologyKind, resolve,
        },
    };
    for (hand, starts, fit) in [
        (ThreadHand::Right, 1, ThreadFit::Metric6gExternal),
        (ThreadHand::Left, 1, ThreadFit::Metric6gExternal),
        (ThreadHand::Right, 1, ThreadFit::NominalBasic),
        (ThreadHand::Left, 1, ThreadFit::NominalBasic),
    ] {
        let mut session =
            Session::new(Document::new("Original thread references").unwrap()).unwrap();
        let compiler = Compiler {
            omit_display: true,
            ..Default::default()
        };
        let actor = || Actor::local(Channel::Ui);
        let parameters = ThreadParameters {
            major_diameter: length(6.),
            pitch: length(1.),
            length: length(6.6),
            fit,
            starts,
            hand,
            frame: PrimitiveFrame {
                origin: [length(17.), length(-3.), length(21.)],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
            ..Default::default()
        };
        let request = ThreadRequest::new("Original M6", parameters.clone());
        let prepared = model::prepare_thread(&session, 0, actor(), &request).unwrap();
        session.commit(prepared.prepared).unwrap();
        let original = compiler.compile(session.document(), 1).unwrap();
        let original_snapshot = &original.snapshot().shapes[&request.feature].topology;
        let intents = original_snapshot
            .elements
            .iter()
            .filter(|e| {
                e.origins.iter().any(|o| {
                    o.role == "cap.end"
                        || o.role.starts_with("cap.end/rim/")
                        || o.role.starts_with("cap.end/corner/")
                        || (e.key.kind == TopologyKind::Face && o.role.starts_with("thread/start/"))
                        || o.role.ends_with("/root/leading_rail")
                        || o.role.ends_with("/lower/leading_rail")
                })
            })
            .map(|e| {
                (
                    GeometryIntent::capture(original_snapshot, e.key).unwrap(),
                    e.geometry,
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(intents.len(), 14 * starts as usize + 1);
        let entrance = original_snapshot
            .elements
            .iter()
            .find(|e| {
                e.key.kind == TopologyKind::Edge
                    && e.origins
                        .iter()
                        .any(|o| o.role.starts_with("cap.start/rim/") && o.role.ends_with("/crest"))
            })
            .unwrap();
        let entrance = GeometryIntent::capture(original_snapshot, entrance.key).unwrap();
        let cap = original_snapshot
            .elements
            .iter()
            .find(|e| {
                e.key.kind == TopologyKind::Face && e.origins.iter().any(|o| o.role == "cap.end")
            })
            .unwrap();
        let mut attachment = GeometryIntent::capture(original_snapshot, cap.key).unwrap();
        attachment.attachment = Some(Attachment {
            anchor: cap.frames[0].origin.clone(),
            offset_mm: [0.2, -0.3, 0.4],
            angle_rad: 0.37,
        });
        session
            .apply(Request {
                expected_revision: 1,
                actor: actor(),
                label: "Attach original far cap".into(),
                edits: vec![Mutation::Insert(Entity::new(
                    "Far cap datum",
                    EntityData::Plane {
                        origin_mm: [0.; 3],
                        normal: [0., 0., 1.],
                        anchor: Some(Link {
                            target: request.feature,
                            role: "face_attachment".into(),
                            geometry: Some(attachment.clone()),
                        }),
                    },
                ))],
            })
            .unwrap();
        let check = |built: &multicad_model::Model, revision| {
            let snapshot = &built.snapshot().shapes[&request.feature].topology;
            for (intent, old) in &intents {
                let Resolution::Resolved { geometry, .. } = resolve(intent, snapshot, revision)
                else {
                    panic!("Lost unchanged finite reference {:?}", intent.origins)
                };
                if let (
                    Some(AnalyticGeometry::Point { position_mm: a }),
                    Some(AnalyticGeometry::Point { position_mm: b }),
                ) = (old, geometry)
                {
                    assert!((0..3).map(|i| (a[i] - b[i]).abs()).sum::<f64>() < 1e-7);
                }
            }
            assert!(matches!(
                resolve(&attachment, snapshot, revision),
                Resolution::Resolved { frame: Some(_), .. }
            ));
        };
        for c in [0.3, 0.4, 0.] {
            let mut edited = parameters.clone();
            if c > 0. {
                edited.end_treatments = Some(ThreadEndTreatments {
                    start: ThreadStartTreatment::Chamfer45 {
                        radial_depth: length(c),
                    },
                    end: ThreadEndTreatment::FlatCut,
                });
            }
            let revision = session.revision();
            let prepared =
                model::prepare_thread_edit(&session, revision, actor(), request.feature, &edited)
                    .unwrap();
            let built = compiler
                .compile(prepared.prepared.document(), revision + 1)
                .unwrap();
            check(&built, revision + 1);
            if c > 0. {
                assert!(matches!(
                    resolve(
                        &entrance,
                        &built.snapshot().shapes[&request.feature].topology,
                        revision + 1
                    ),
                    Resolution::Missing
                ));
            }
            session.commit(prepared.prepared).unwrap();
        }
        let revision = session.revision();
        session.undo(revision, actor()).unwrap();
        let built = compiler
            .compile(session.document(), session.revision())
            .unwrap();
        check(&built, session.revision());
        session.redo(session.revision(), actor()).unwrap();
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
        .project;
        let built = compiler.compile(&reopened.document, 0).unwrap();
        check(&built, 0);
    }
}
