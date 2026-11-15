use multicad_core::references::TopologyKind;
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::imports::prepare_step_import,
    project_file::{self, ProjectFile},
};
use multicad_kernel::MeshOptions;
use multicad_model::{Compiler, GeometryBackend, imports::read_step_file_with_backend};
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
                .rust_literal_npt_contains([11. + 2. * a.sin(), -5., 3. + 2. * a.cos()])
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
