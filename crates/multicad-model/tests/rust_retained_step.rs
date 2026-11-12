use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityId},
    exchange::{ExportReply, ExportResult, Scope, ScopeRoot, SetDefinition, prepare_set},
    model::imports::prepare_step_import,
    project_file::{self, ProjectFile},
};
use multicad_kernel::{
    Solid,
    rust_kernel::{self, helical_thread::Spec, tapered_thread::Thread},
};
use multicad_model::{
    Compiler, GeometryBackend,
    export::{
        ExportWorker, Exporter, StepVerificationMethod, verify_rust_step_shape,
        verify_rust_step_shape_receipt,
    },
    imports::read_step_file_with_backend,
};
use std::{
    cell::Cell,
    path::Path,
    time::{Duration, Instant},
};

fn reply(worker: &ExportWorker, id: EntityId) -> ExportReply {
    let deadline = Instant::now() + Duration::from_secs(90);
    loop {
        if let Some(reply) = worker.poll() {
            assert_eq!(reply.id, id);
            return reply;
        }
        assert!(
            Instant::now() < deadline,
            "Retained export worker did not reply"
        );
        std::thread::sleep(Duration::from_millis(2));
    }
}
fn source() -> rust_kernel::radial_tensor::RadialBoundary {
    Thread::new(
        Spec {
            major_radius_mm: 4.,
            pitch_mm: 1.25,
            length_mm: 1.25,
            depth_mm: 0.5,
            crest_width_mm: 0.25,
            lower_flank_deg: 12.,
            upper_flank_deg: 18.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        },
        0.03125,
        1e-7,
        1e-6,
    )
    .unwrap()
    .retained_boundary_with_interrupt([[0.; 3], [1., 0., 0.], [0., 1., 0.]], 1e-7, 1e-6, None)
    .unwrap()
    .boundary
}
fn exchange(source: &rust_kernel::radial_tensor::RadialBoundary, mixed: bool) -> String {
    let mut writer = rust_kernel::step::Writer::new();
    writer.append_radial_boundary(source, 1e-7).unwrap();
    if mixed {
        writer
            .append_box(
                [2, 3, 4].map(|v| rust_kernel::Rational::new(v, 1).unwrap()),
                [20., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
    }
    writer.finish().unwrap()
}

#[test]
fn retained_step_worker_publishes_reopens_and_keeps_interval_mass_and_atomic_rollback() {
    for mixed in [false, true] {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("retained.step");
        let boundary = source();
        let original = exchange(&boundary, mixed);
        std::fs::write(&input, &original).unwrap();
        let file = read_step_file_with_backend(&input, GeometryBackend::Rust).unwrap();
        let mut session = Session::new(Document::new("Retained finite source").unwrap()).unwrap();
        let actor = Actor::local(Channel::Ui);
        let prepared = prepare_step_import(&session, 0, actor.clone(), &file).unwrap();
        let part = prepared.model.outputs[0].body;
        session.commit(prepared.model.prepared).unwrap();
        let set = EntityId::new();
        let prepared = prepare_set(
            &session,
            session.revision(),
            actor.clone(),
            set,
            "Retained STEP",
            &SetDefinition {
                scope: Scope {
                    root: ScopeRoot::Explicit {
                        entities: vec![part],
                    },
                    tags: None,
                },
                ..Default::default()
            },
        )
        .unwrap();
        session.commit(prepared).unwrap();
        std::fs::remove_file(input).unwrap();
        let mut saved = Vec::new();
        project_file::write(
            &ProjectFile::new(session.document().clone()).unwrap(),
            &mut saved,
            Default::default(),
        )
        .unwrap();
        let reopened = project_file::read(
            &mut saved.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap()
        .project;
        assert_eq!(&reopened.document, session.document());
        let compiler = Compiler {
            omit_display: true,
            ..Default::default()
        };
        let model = compiler
            .compile_with_backend(
                &reopened.document,
                session.revision(),
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        let snapshot = model.snapshot();
        let producer = *snapshot
            .bodies
            .get(&part)
            .or_else(|| snapshot.body_sources.get(&part))
            .unwrap();
        let solid = model.solid(producer).unwrap();
        let shape = &snapshot.shapes[&producer];
        assert!(solid.has_retained_radial_source());
        let worker =
            ExportWorker::with_geometry_backend(Exporter::default(), GeometryBackend::Rust)
                .unwrap();
        let id = worker
            .preview(session.document(), session.revision(), set, actor.clone())
            .unwrap();
        let ExportResult::Preview(preview) = reply(&worker, id).result else {
            panic!()
        };
        let preview = preview.unwrap();
        assert!(preview.can_run(), "{preview:?}");
        let output = tempfile::tempdir().unwrap();
        let id = worker
            .run(
                session.document(),
                session.revision(),
                &preview.key,
                output.path(),
                actor.clone(),
            )
            .unwrap();
        let ExportResult::Run(run) = reply(&worker, id).result else {
            panic!()
        };
        let run = run.unwrap();
        assert_eq!(run.files.len(), 1);
        let path = Path::new(&run.directory).join(&run.files[0].artifact.output.filename);
        let published = std::fs::read_to_string(&path).unwrap();
        let proof = verify_rust_step_shape_receipt(&path, solid, shape, 1e-7).unwrap();
        assert_eq!(
            proof.method,
            StepVerificationMethod::DecodedBoundaryEnclosure
        );
        assert_eq!(proof.decoded_volume_mm3, None);
        assert_eq!(proof.retained_members.len(), 1);
        let member = &proof.retained_members[0];
        assert!(member.decoded_volume_bounds_mm3[0] > 0.);
        assert!(member.decoded_volume_bounds_mm3[0] <= member.source_volume_bounds_mm3[1]);
        assert!(member.source_volume_bounds_mm3[0] <= member.decoded_volume_bounds_mm3[1]);
        assert_eq!(proof.decoded_scalar_members.len(), usize::from(mixed));
        assert!(verify_rust_step_shape(&path, solid, shape, 1e-7).is_err());
        let reread =
            Solid::import_rust_step(&published, EntityId::new(), Default::default()).unwrap();
        assert!(reread.has_retained_radial_source());
        assert_eq!(reread.solid_count().unwrap(), if mixed { 2 } else { 1 });
        let decoded = read_step_file_with_backend(&path, GeometryBackend::Rust).unwrap();
        assert_eq!(decoded.solids.len(), 1);
        let mut stale = shape.as_ref().clone();
        stale.bounds.min[0] += 1.;
        assert!(verify_rust_step_shape_receipt(&path, solid, &stale, 1e-7).is_err());
        let mut stale = shape.as_ref().clone();
        stale.volume_mm3 = f64::NAN;
        assert!(verify_rust_step_shape_receipt(&path, solid, &stale, 1e-7).is_err());
        let moved = boundary
            .transformed_with_interrupt(
                [[1., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
                1e-7,
                1e-6,
                None,
            )
            .unwrap();
        let wrong = exchange(&moved, mixed);
        let failed_root = tempfile::tempdir().unwrap();
        let poisoned = Cell::new(false);
        let poison = || {
            if !poisoned.get() {
                for entry in std::fs::read_dir(failed_root.path()).unwrap().flatten() {
                    if entry
                        .file_name()
                        .to_string_lossy()
                        .starts_with(".multicad-export-")
                    {
                        let staged = entry.path().join(&preview.artifacts[0].output.filename);
                        if staged.is_file() {
                            std::fs::write(staged, &wrong).unwrap();
                            poisoned.set(true);
                        }
                    }
                }
            }
            false
        };
        let before = session.document().clone();
        let failure = Exporter::default()
            .run(
                session.document(),
                session.revision(),
                &preview.key,
                failed_root.path(),
                &actor,
                Some(&poison),
            )
            .unwrap_err();
        assert_eq!(failure.kind, multicad_core::exchange::ErrorKind::Invalid);
        assert!(poisoned.get(), "Actual staged STEP was not reached");
        assert!(std::fs::read_dir(failed_root.path()).unwrap().all(|entry| {
            entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .ends_with(".lock")
        }));
        assert_eq!(std::fs::read_to_string(&path).unwrap(), published);
        assert_eq!(session.document(), &before);
    }
}
