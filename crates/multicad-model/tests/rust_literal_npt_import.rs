use multicad_core::references::TopologyKind;
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::imports::prepare_step_import,
    project_file::{self, ProjectFile},
};
use multicad_kernel::MeshOptions;
use multicad_model::{Compiler, GeometryBackend, imports::read_step_file_with_backend};
use multicad_precision::PortableMath as _;
fn stock(
    quarter: bool,
) -> multicad_kernel::rust_kernel::axial_pipe_thread::authored::StockClipPlan {
    use multicad_kernel::rust_kernel::axial_pipe_thread::authored::{
        AuthoredSection, AuthoredSpec, AuthoredThread,
    };
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let pitch = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let s = 1. / (m - k) + 1. / (m + k);
    let depth = pitch / s - (0.033 + max) * pitch;
    let spec = AuthoredSpec {
        major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + depth) / 2.,
        pitch_mm: pitch,
        length_mm: l2 * 25.4,
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: depth,
            flat_width_mm: (pitch - depth * s) / 2.,
        },
        crest_radius_limit_mm: Some(od * 25.4 / 2.),
    };
    AuthoredThread::new(spec, 1e-7, 1e-6)
        .unwrap()
        .stock_clip_plan(Default::default(), None)
        .unwrap()
}
fn text(quarter: bool) -> String {
    use multicad_kernel::rust_kernel::step::npt_stock as step;
    let source = stock(quarter);
    let curves = step::CurvePlan::new(&source, Default::default(), None).unwrap();
    let surfaces = step::SurfacePlan::new(curves, Default::default(), None).unwrap();
    step::write_source_candidate(&surfaces, 16_000_000, 64 * 1024 * 1024, None)
        .unwrap()
        .into_string()
}

#[test]
fn actual_literal_npt_import_compiles_displays_and_reopens_without_external_source() {
    for quarter in [false, true] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pipe.step");
        let original = text(quarter);
        std::fs::write(&path, &original).unwrap();
        let imported = read_step_file_with_backend(&path, GeometryBackend::Rust).unwrap();
        assert_eq!(imported.solids.len(), 1);
        assert_eq!(imported.solids[0].1, original);
        let expected = imported.solids[0].2.volume_mm3;
        let mut session = Session::new(Document::new("Literal pipe import").unwrap()).unwrap();
        let prepared =
            prepare_step_import(&session, 0, Actor::local(Channel::Ui), &imported).unwrap();
        let part = prepared.model.outputs[0].body;
        session.commit(prepared.model.prepared).unwrap();
        std::fs::remove_file(&path).unwrap();
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
        let compiler = Compiler::default();
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
        assert!(solid.rust_literal_npt_source().is_some());
        assert!(solid.rust_clipped_npt_source().is_none());
        assert!((solid.volume_mm3().unwrap() - expected).abs() < 1e-8);
        let display = solid
            .topology_with_display(producer, session.revision(), Some(0.05))
            .unwrap();
        for (kind, count) in [
            (TopologyKind::Face, 10),
            (TopologyKind::Edge, 20),
            (TopologyKind::Vertex, 12),
        ] {
            assert_eq!(
                display
                    .topology
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == kind)
                    .count(),
                count
            );
        }
        assert_eq!(display.edges.len(), 20);
        assert_eq!(display.vertices.len(), 12);
        let before = solid.rust_literal_npt_runtime_evidence().unwrap();
        let manufactured = solid
            .manufacturing_mesh_with_receipt(MeshOptions {
                linear_deflection_mm: 0.1,
                angular_deflection_rad: 0.3,
                maximum_triangles: 100_000,
            })
            .unwrap();
        assert!(manufactured.edge_pairs > 0 && manufactured.vertex_links > 0);
        assert!(manufactured.mesh.face_indices.iter().all(|f| *f < 10));
        let after = solid.rust_literal_npt_runtime_evidence().unwrap();
        assert_eq!(after.mesh_builds, before.mesh_builds + 1);
        assert_eq!(after.resident_budget_bytes, before.resident_budget_bytes);
        let cached = solid.mesh(Default::default()).unwrap();
        assert_eq!(cached.vertices, snapshot.shapes[&producer].mesh.vertices);
        let hot = solid.rust_literal_npt_runtime_evidence().unwrap();
        assert_eq!(hot.mesh_builds, after.mesh_builds);
        assert!(hot.last_mesh_actual_work.unwrap() < after.last_mesh_actual_work.unwrap());
        assert!(
            compiler
                .compile_with_backend(
                    &reopened.document,
                    session.revision(),
                    GeometryBackend::Rust,
                    Some(&|| true)
                )
                .is_err()
        );
        assert!(solid.rust_literal_npt_contains([0., 0., 2.]).unwrap());
        assert!(!solid.rust_literal_npt_contains([100., 0., 2.]).unwrap());
        let stable = session.document().clone();
        let revision = session.revision();
        std::fs::write(&path, original.replace("CLOSED_SHELL", "OPEN_SHELL")).unwrap();
        assert!(read_step_file_with_backend(&path, GeometryBackend::Rust).is_err());
        assert_eq!(session.document(), &stable);
        assert_eq!(session.revision(), revision);

        use multicad_core::model::{
            self, BodyInput, NumericTransform, TransformRequest, angle, length,
        };
        let request = TransformRequest::new(
            BodyInput {
                body: part,
                producer,
            },
            NumericTransform {
                translation: [length(11.), length(-5.), length(3.)],
                axis: [0., 1., 0.],
                angle: angle(37.),
                ..Default::default()
            },
        );
        let prepared =
            model::prepare_transform(&session, revision, Actor::local(Channel::Ui), &request)
                .unwrap();
        let posed = compiler
            .compile_with_backend(
                prepared.prepared.document(),
                revision + 1,
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        assert_eq!(session.document(), &stable);
        let moved = posed.solid(request.feature).unwrap();
        assert!(moved.rust_literal_npt_occurrence().is_some());
        assert_eq!(
            moved
                .rust_literal_npt_source()
                .unwrap()
                .world_enclosure_mm(),
            solid
                .rust_literal_npt_source()
                .unwrap()
                .world_enclosure_mm()
        );
        assert!((moved.volume_mm3().unwrap() - expected).abs() < 1e-8);
        let a = 37f64.to_radians();
        assert!(
            moved
                .rust_literal_npt_contains([
                    11. + 2. * a.portable_sin(),
                    -5.,
                    3. + 2. * a.portable_cos()
                ])
                .unwrap()
        );
        assert!(!moved.rust_literal_npt_contains([111., -5., 3.]).unwrap());
        let proof = moved
            .manufacturing_mesh_with_receipt(MeshOptions {
                linear_deflection_mm: 0.1,
                angular_deflection_rad: 0.3,
                maximum_triangles: 100_000,
            })
            .unwrap();
        assert!(proof.edge_pairs > 0 && proof.vertex_links > 0);
        session.commit(prepared.prepared).unwrap();
        let mut moved_saved = Vec::new();
        project_file::write(
            &ProjectFile::new(session.document().clone()).unwrap(),
            &mut moved_saved,
            Default::default(),
        )
        .unwrap();
        let moved_reopened = project_file::read(
            &mut moved_saved.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap()
        .project;
        assert_eq!(&moved_reopened.document, session.document());
        let rebuilt = compiler
            .compile_with_backend(
                &moved_reopened.document,
                session.revision(),
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        assert!(
            (rebuilt
                .solid(request.feature)
                .unwrap()
                .volume_mm3()
                .unwrap()
                - expected)
                .abs()
                < 1e-8
        );
    }
}

#[test]
fn actual_literal_export_worker_keeps_decoded_mass_intervals_and_atomic_publication() {
    use multicad_core::{
        document::EntityId,
        exchange::{ExportReply, ExportResult, Scope, ScopeRoot, SetDefinition, prepare_set},
    };
    use multicad_kernel::{RigidTransform, Solid};
    use multicad_model::export::{
        ExportWorker, Exporter, StepVerificationMethod, verify_rust_step_shape_receipt,
    };
    use std::{
        cell::Cell,
        path::Path,
        time::{Duration, Instant},
    };
    fn reply(worker: &ExportWorker, id: EntityId) -> ExportReply {
        let deadline = Instant::now() + Duration::from_secs(90);
        loop {
            if let Some(r) = worker.poll() {
                assert_eq!(r.id, id);
                return r;
            }
            assert!(Instant::now() < deadline, "Literal export worker timed out");
            std::thread::sleep(Duration::from_millis(2));
        }
    }
    for mixed in [false, true] {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("actual.step");
        let original = text(false);
        let source =
            Solid::import_rust_step(&original, EntityId::new(), Default::default()).unwrap();
        let pose = RigidTransform {
            translation_mm: [11., -5., 3.],
            ..Default::default()
        };
        let moved = source.transformed(pose).unwrap();
        let mut writer = multicad_kernel::rust_kernel::step::Writer::new();
        let candidate = multicad_kernel::rust_kernel::step::literal_npt::verify_placed_replay(
            source.rust_literal_npt_source().unwrap(),
            moved.rust_literal_npt_occurrence().unwrap(),
            Default::default(),
            None,
        )
        .unwrap();
        writer
            .append_verified_literal_npt(&candidate, Default::default(), None)
            .unwrap();
        if mixed {
            writer
                .append_box(
                    [2, 3, 4].map(|n| multicad_kernel::rust_kernel::Rational::new(n, 1).unwrap()),
                    [-20., 0., 0.],
                    [1., 0., 0.],
                    [0., 1., 0.],
                )
                .unwrap();
            let candidate = multicad_kernel::rust_kernel::step::literal_npt::verify_exact_replay(
                source.rust_literal_npt_source().unwrap(),
                Default::default(),
                None,
            )
            .unwrap();
            writer
                .append_verified_literal_npt(&candidate, Default::default(), None)
                .unwrap();
        }
        std::fs::write(&input, writer.finish().unwrap()).unwrap();
        let imported = read_step_file_with_backend(&input, GeometryBackend::Rust).unwrap();
        let mut session = Session::new(Document::new("Literal export workflow").unwrap()).unwrap();
        let actor = Actor::local(Channel::Ui);
        let prepared = prepare_step_import(&session, 0, actor.clone(), &imported).unwrap();
        let part = prepared.model.outputs[0].body;
        session.commit(prepared.model.prepared).unwrap();
        std::fs::remove_file(&input).unwrap();
        use multicad_core::model::{
            self, BodyInput, NumericTransform, TransformRequest, angle, length,
        };
        let original_model = Compiler {
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
        let producer = original_model.snapshot().bodies[&part];
        let transform = TransformRequest::new(
            BodyInput {
                body: part,
                producer,
            },
            NumericTransform {
                translation: [length(3.), length(-2.), length(5.)],
                axis: [1., 2., 3.],
                angle: angle(19.),
                ..Default::default()
            },
        );
        let prepared =
            model::prepare_transform(&session, session.revision(), actor.clone(), &transform)
                .unwrap();
        session.commit(prepared.prepared).unwrap();

        let set = EntityId::new();
        let prepared = prepare_set(
            &session,
            session.revision(),
            actor.clone(),
            set,
            "Actual STEP",
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
        session = Session::new(reopened.document).unwrap();
        let model = Compiler::default()
            .compile_with_backend(
                session.document(),
                session.revision(),
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        let snapshot = model.snapshot();
        let producer = snapshot.bodies[&part];
        let solid = model.solid(producer).unwrap();
        let shape = &snapshot.shapes[&producer];
        let direct = directory.path().join("direct.step");
        solid.write_step(&direct).unwrap();
        let actual = std::fs::read_to_string(&direct).unwrap();
        if !mixed {
            solid
                .verify_literal_npt_step(&actual, 1e-7)
                .unwrap_or_else(|e| panic!("Fresh direct mixed={mixed}: {e:?}"));
        }
        let kinds = multicad_kernel::rust_kernel::step_read::read_source_compound(&actual, 1e-7)
            .unwrap()
            .iter()
            .map(|m| match m {
                multicad_kernel::rust_kernel::step_read::SourceImport::ActualMappedHead(_) => {
                    "actual_mapped_head"
                }
                multicad_kernel::rust_kernel::step_read::SourceImport::LiteralNpt(_) => "literal",
                multicad_kernel::rust_kernel::step_read::SourceImport::Analytic(_) => "analytic",
                multicad_kernel::rust_kernel::step_read::SourceImport::Radial(_) => "radial",
                multicad_kernel::rust_kernel::step_read::SourceImport::PipeSleeve(_) => {
                    "pipe_sleeve"
                }
                multicad_kernel::rust_kernel::step_read::SourceImport::PeriodicPhillips(_) => {
                    "periodic_phillips"
                }
            })
            .collect::<Vec<_>>();
        solid
            .verify_literal_npt_step_compound(&actual, 1e-7)
            .unwrap()
            .unwrap_or_else(|e| panic!("Fresh compound mixed={mixed} kinds={kinds:?}: {e:?}"));
        let worker =
            ExportWorker::with_geometry_backend(Exporter::default(), GeometryBackend::Rust)
                .unwrap();
        let id = worker
            .preview(session.document(), session.revision(), set, actor.clone())
            .unwrap();
        let ExportResult::Preview(preview) = reply(&worker, id).result else {
            panic!("Expected preview");
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
            panic!("Expected run");
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
        assert_eq!(proof.retained_members.len(), if mixed { 2 } else { 1 });
        assert_eq!(proof.decoded_scalar_members.len(), usize::from(mixed));
        for member in &proof.retained_members {
            assert!(member.decoded_volume_bounds_mm3[0] > 0.);
            assert!(member.decoded_volume_bounds_mm3[1] >= member.decoded_volume_bounds_mm3[0]);
            assert!(member.source_volume_bounds_mm3[0] > 0.);
        }
        let reopened =
            Solid::import_rust_step(&published, EntityId::new(), Default::default()).unwrap();
        assert_eq!(reopened.solid_count().unwrap(), if mixed { 3 } else { 1 });
        assert!((reopened.volume_mm3().unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-5);
        let mut stale = shape.as_ref().clone();
        stale.bounds.min[0] += 1.;
        assert!(verify_rust_step_shape_receipt(&path, solid, &stale, 1e-7).is_err());
        let failed_root = tempfile::tempdir().unwrap();
        let poisoned = Cell::new(false);
        let mut wrong = multicad_kernel::rust_kernel::step::Writer::new();
        wrong
            .append_box(
                [7, 8, 9].map(|n| multicad_kernel::rust_kernel::Rational::new(n, 1).unwrap()),
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
        let wrong = wrong.finish().unwrap();
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
        let error = Exporter::default()
            .run(
                session.document(),
                session.revision(),
                &preview.key,
                failed_root.path(),
                &actor,
                Some(&poison),
            )
            .unwrap_err();
        assert!(poisoned.get());
        assert_eq!(error.kind, multicad_core::exchange::ErrorKind::Invalid);
        assert_eq!(session.document(), &before);
        assert_eq!(std::fs::read_to_string(&path).unwrap(), published);
        assert!(
            std::fs::read_dir(failed_root.path()).unwrap().all(|e| e
                .unwrap()
                .file_name()
                .to_string_lossy()
                .ends_with(".lock"))
        );
    }
}
