use super::*;
use multicad_core::model::{BodyInput, NumericTransform, TransformRequest};
use multicad_core::{
    document::EntityId,
    exchange::{ExportReply, ExportResult, Scope, ScopeRoot, SetDefinition, prepare_set},
};
use multicad_model::export::{ExportWorker, Exporter};
use multicad_model::export::{StepVerificationMethod, verify_rust_step_shape_receipt};
fn reply(worker: &ExportWorker, id: EntityId) -> ExportReply {
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(120);
    loop {
        if let Some(result) = worker.poll() {
            assert_eq!(result.id, id);
            return result;
        }
        assert!(
            std::time::Instant::now() < deadline,
            "Periodic export worker timed out"
        );
        std::thread::sleep(std::time::Duration::from_millis(2));
    }
}

#[test]
fn actual_nominal_phillips_file_receipts_survive_transform_and_project_reopen() {
    for (kind, left) in [(0, false), (0, true), (5, false), (5, true)] {
        let mut session =
            Session::new(Document::new("Phillips export lifecycle").unwrap()).unwrap();
        let request = MechanicalPartRequest::new(
            "Nominal physical Phillips",
            recipe(kind, ThreadFit::NominalBasic, left, 3.5),
        );
        let prepared = model::prepare_mechanical_part(&session, 0, actor(), &request).unwrap();
        session.commit(prepared.prepared).unwrap();
        if left {
            let transformed = model::prepare_transform(
                &session,
                session.revision(),
                actor(),
                &TransformRequest::new(
                    BodyInput {
                        body: request.body,
                        producer: request.feature,
                    },
                    NumericTransform {
                        translation: [length(11.), length(-7.), length(5.)],
                        axis: [1., 2., 3.],
                        angle: angle(19.),
                        ..Default::default()
                    },
                ),
            )
            .unwrap();
            session.commit(transformed.prepared).unwrap();
        }
        let set = EntityId::new();
        session
            .commit(
                prepare_set(
                    &session,
                    session.revision(),
                    actor(),
                    set,
                    "Physical Phillips",
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
        let session = Session::new(reopened.document).unwrap();
        let built = compiler()
            .compile_with_backend(
                session.document(),
                session.revision(),
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        let snapshot = built.snapshot();
        let producer = snapshot.bodies[&request.body];
        let source = built.solid(producer).unwrap();
        assert!(source.rust_periodic_repaired_phillips_source().is_some());
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("actual.step");
        assert!(source.supports_step_export());
        source.write_step(&path).unwrap();
        let original = std::fs::read_to_string(&path).unwrap();
        let imported = multicad_kernel::Solid::import_rust_step(
            &original,
            multicad_core::document::EntityId::new(),
            Default::default(),
        )
        .unwrap();
        assert!(imported.rust_periodic_repaired_phillips_source().is_some());
        let file =
            multicad_model::imports::read_step_file_with_backend(&path, GeometryBackend::Rust)
                .unwrap();
        let mut imported_session =
            Session::new(Document::new("Cold periodic import").unwrap()).unwrap();
        let prepared = multicad_core::model::imports::prepare_step_import(
            &imported_session,
            0,
            actor(),
            &file,
        )
        .unwrap();
        let imported_body = prepared.model.outputs[0].body;
        imported_session.commit(prepared.model.prepared).unwrap();
        let mut saved_import = Vec::new();
        project_file::write(
            &ProjectFile::new(imported_session.document().clone()).unwrap(),
            &mut saved_import,
            Default::default(),
        )
        .unwrap();
        std::fs::remove_file(&path).unwrap();
        let cold = project_file::read(
            &mut saved_import.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap()
        .project;
        let cold_session = Session::new(cold.document).unwrap();
        let cold_model = Compiler::default()
            .compile_with_backend(
                cold_session.document(),
                cold_session.revision(),
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        let cold_snapshot = cold_model.snapshot();
        let cold_solid = cold_model
            .solid(cold_snapshot.bodies[&imported_body])
            .unwrap();
        assert!(
            cold_solid
                .rust_periodic_repaired_phillips_source()
                .is_some()
        );
        assert_eq!(cold_solid.solid_count().unwrap(), 1);
        assert!(
            !cold_snapshot.shapes[&cold_snapshot.bodies[&imported_body]]
                .mesh
                .indices
                .is_empty()
        );
        assert!((cold_solid.volume_mm3().unwrap() - source.volume_mm3().unwrap()).abs() < 1e-5);
        std::fs::write(&path, &original).unwrap();
        let receipt =
            verify_rust_step_shape_receipt(&path, source, &snapshot.shapes[&producer], 1e-7)
                .unwrap();
        assert_eq!(
            receipt.method,
            StepVerificationMethod::DecodedPeriodicRepair
        );
        assert_eq!(receipt.decoded_volume_mm3, None);
        assert!(receipt.healed_members.is_empty() && receipt.retained_members.is_empty());
        assert_eq!(receipt.periodic_members.len(), 1);
        let member = &receipt.periodic_members[0];
        assert_eq!(member.member_index, 0);
        assert!(member.decoded_volume_bounds_mm3[0] > 0.);
        assert!(member.maximum_world_position_error_mm <= member.verified_tolerance_mm);
        assert!(member.maximum_world_chart_derivative_error_mm <= member.verified_tolerance_mm);
        assert!(member.reader_reserved_work > 0 && member.correspondence_cavity_cells > 0);
        let worker =
            ExportWorker::with_geometry_backend(Exporter::default(), GeometryBackend::Rust)
                .unwrap();
        let id = worker
            .preview(session.document(), session.revision(), set, actor())
            .unwrap();
        let ExportResult::Preview(preview) = reply(&worker, id).result else {
            panic!("Expected export preview")
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
                actor(),
            )
            .unwrap();
        let ExportResult::Run(result) = reply(&worker, id).result else {
            panic!("Expected export run")
        };
        let result = result.unwrap();
        let exported =
            std::path::Path::new(&result.directory).join(&result.files[0].artifact.output.filename);
        let verified =
            verify_rust_step_shape_receipt(&exported, source, &snapshot.shapes[&producer], 1e-7)
                .unwrap();
        assert_eq!(
            verified.method,
            StepVerificationMethod::DecodedPeriodicRepair
        );
        let mut stale = snapshot.shapes[&producer].as_ref().clone();
        stale.bounds.min[0] += 1.;
        assert!(verify_rust_step_shape_receipt(&path, source, &stale, 1e-7).is_err());
        let mut wrong = multicad_kernel::rust_kernel::step::Writer::new();
        wrong
            .append_box(
                [2, 3, 4].map(|n| multicad_kernel::rust_kernel::Rational::new(n, 1).unwrap()),
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
        std::fs::write(&path, wrong.finish().unwrap()).unwrap();
        assert!(
            verify_rust_step_shape_receipt(&path, source, &snapshot.shapes[&producer], 1e-7)
                .is_err()
        );
        std::fs::write(&path, original).unwrap();
        assert!(
            verify_rust_step_shape_receipt(&path, source, &snapshot.shapes[&producer], 1e-7)
                .is_ok()
        );
    }
}

#[test]
fn actual_mixed_periodic_and_ordinary_export_keeps_explicit_member_mass_provenance() {
    let mut session = Session::new(Document::new("Mixed periodic sources").unwrap()).unwrap();
    let mut requests = Vec::new();
    for left in [false, true] {
        let request = MechanicalPartRequest::new(
            "Physical M2",
            recipe(0, ThreadFit::NominalBasic, left, 3.5),
        );
        let candidate =
            model::prepare_mechanical_part(&session, session.revision(), actor(), &request)
                .unwrap();
        session.commit(candidate.prepared).unwrap();
        requests.push(request);
    }
    let model = compiler()
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let right = model.solid(requests[0].feature).unwrap();
    let left = model
        .solid(requests[1].feature)
        .unwrap()
        .transformed(multicad_kernel::RigidTransform {
            translation_mm: [20., -5., 8.],
            ..Default::default()
        })
        .unwrap();
    let mut ordinary = multicad_kernel::rust_kernel::step::Writer::new();
    ordinary
        .append_box(
            [2, 3, 4].map(|n| multicad_kernel::rust_kernel::Rational::new(n, 1).unwrap()),
            [-20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    let ordinary = multicad_kernel::Solid::import_rust_step(
        &ordinary.finish().unwrap(),
        EntityId::new(),
        Default::default(),
    )
    .unwrap();
    let compound = multicad_kernel::Solid::compound_named(
        &[right, &ordinary, &left],
        EntityId::new(),
        Default::default(),
    )
    .unwrap();
    assert!(compound.supports_step_export());
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("physical-mixed.step");
    compound.write_step(&path).unwrap();
    let actual = std::fs::read_to_string(&path).unwrap();
    let proof = compound
        .verify_partitioned_step_source(&actual, 1e-7, None)
        .unwrap()
        .unwrap();
    assert_eq!(proof.members.solid_count, 3);
    assert!(proof.document_ownership_work > 0 && proof.topology_entity_count > 0);
    assert_eq!(proof.periodic_members.len(), 2);
    assert_eq!(proof.members.decoded_scalar_members, vec![(1, 24.)]);
    assert!(proof.members.retained_members.is_empty() && proof.healed_members.is_empty());
    let mut orphan = actual.clone();
    orphan.insert_str(
        orphan.rfind("ENDSEC;").unwrap(),
        "#999999=VERTEX_POINT('',#1);\n",
    );
    assert_eq!(
        compound
            .verify_partitioned_step_source(&orphan, 1e-7, None)
            .unwrap()
            .unwrap_err()
            .kind,
        multicad_kernel::GeometryErrorKind::InvalidResult
    );
    assert_eq!(
        multicad_kernel::Solid::import_rust_step(&orphan, EntityId::new(), Default::default())
            .err()
            .unwrap()
            .kind,
        multicad_kernel::GeometryErrorKind::InvalidResult
    );
    assert_eq!(
        compound
            .verify_partitioned_step_source(&actual, 1e-7, Some(&|| true))
            .unwrap()
            .unwrap_err()
            .kind,
        multicad_kernel::GeometryErrorKind::Cancelled
    );

    for ((index, source), (actual_index, member)) in [(0, right), (2, &left)]
        .into_iter()
        .zip(&proof.periodic_members)
    {
        assert_eq!(index, *actual_index);
        assert_eq!(member.represented_member_count, 3);
        assert_eq!(member.member_index, index as usize);
        assert!(member.selection_work > 0 && member.decoded_volume_bounds_mm3[0] > 0.);
        source
            .verify_periodic_phillips_step_member(&actual, index as usize, 1e-7, None)
            .unwrap()
            .unwrap();
        assert!(
            source
                .verify_periodic_phillips_step_member(&actual, 1, 1e-7, None)
                .unwrap()
                .is_err()
        );
    }
    assert!(
        right
            .verify_periodic_phillips_step_source(&actual, 1e-7, None)
            .unwrap()
            .is_err()
    );
    // Import the complete actual file as one canonical Part. The project must
    // own its physical source after the external STEP file disappears.
    let file =
        multicad_model::imports::read_step_file_with_backend(&path, GeometryBackend::Rust).unwrap();
    let mut imported = Session::new(Document::new("Cold mixed Phillips").unwrap()).unwrap();
    let prepared =
        multicad_core::model::imports::prepare_step_import(&imported, 0, actor(), &file).unwrap();
    assert_eq!(prepared.model.outputs.len(), 1);
    let body = prepared.model.outputs[0].body;
    imported.commit(prepared.model.prepared).unwrap();
    let set = EntityId::new();
    imported
        .commit(
            prepare_set(
                &imported,
                imported.revision(),
                actor(),
                set,
                "Reexport physical members",
                &SetDefinition {
                    scope: Scope {
                        root: ScopeRoot::Explicit {
                            entities: vec![body],
                        },
                        tags: None,
                    },
                    ..Default::default()
                },
            )
            .unwrap(),
        )
        .unwrap();
    let mut saved = Vec::new();
    project_file::write(
        &ProjectFile::new(imported.document().clone()).unwrap(),
        &mut saved,
        Default::default(),
    )
    .unwrap();
    std::fs::remove_file(&path).unwrap();
    let reopened = project_file::read(
        &mut saved.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    let imported = Session::new(reopened.document).unwrap();
    let built = Compiler::default()
        .compile_with_backend(
            imported.document(),
            imported.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let snapshot = built.snapshot();
    let producer = snapshot.bodies[&body];
    let solid = built.solid(producer).unwrap();
    assert_eq!(solid.solid_count().unwrap(), 3);
    assert!(!snapshot.shapes[&producer].mesh.indices.is_empty());
    assert!((solid.volume_mm3().unwrap() - compound.volume_mm3().unwrap()).abs() < 1e-5);
    let worker =
        ExportWorker::with_geometry_backend(Exporter::default(), GeometryBackend::Rust).unwrap();
    let id = worker
        .preview(imported.document(), imported.revision(), set, actor())
        .unwrap();
    let ExportResult::Preview(preview) = reply(&worker, id).result else {
        panic!("Expected mixed export preview")
    };
    let preview = preview.unwrap();
    assert!(preview.can_run(), "{preview:?}");
    let output = tempfile::tempdir().unwrap();
    let id = worker
        .run(
            imported.document(),
            imported.revision(),
            &preview.key,
            output.path(),
            actor(),
        )
        .unwrap();
    let ExportResult::Run(result) = reply(&worker, id).result else {
        panic!("Expected mixed export run")
    };
    let result = result.unwrap();
    assert_eq!(result.files.len(), 1);
    let exported =
        std::path::Path::new(&result.directory).join(&result.files[0].artifact.output.filename);
    let verified =
        verify_rust_step_shape_receipt(&exported, solid, &snapshot.shapes[&producer], 1e-7)
            .unwrap();
    assert_eq!(verified.method, StepVerificationMethod::MixedMemberBoundary);
    assert!(verified.document_ownership.as_ref().unwrap().work_used > 0);
    assert!(verified.decoded_volume_mm3.is_none());
    assert_eq!(verified.periodic_members.len(), 2);
    assert_eq!(verified.decoded_scalar_members.len(), 1);
}
