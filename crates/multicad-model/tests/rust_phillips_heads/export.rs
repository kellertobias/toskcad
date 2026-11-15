use super::*;
use multicad_core::model::{BodyInput, NumericTransform, TransformRequest};
use multicad_model::export::{StepVerificationMethod, verify_rust_step_shape_receipt};

#[test]
fn actual_nominal_phillips_file_receipts_survive_transform_and_project_reopen() {
    for left in [false, true] {
        let mut session =
            Session::new(Document::new("Phillips export lifecycle").unwrap()).unwrap();
        let request = MechanicalPartRequest::new(
            "M2 Phillips",
            recipe(0, ThreadFit::NominalBasic, left, 3.5),
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
        source
            .write_periodic_phillips_step_with_verification(&path, 1e-7, None)
            .unwrap()
            .unwrap();
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
