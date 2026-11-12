use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::imports::{prepare_step_import, prepare_step_update},
    project_file::{self, ProjectFile},
    references::{GeometryIntent, Resolution, TopologyKey, TopologyKind, resolve},
};
use multicad_kernel::rust_kernel::{
    Circle, Rational,
    chamfered_bores::{ChamferedBores, Countersink},
    polygon::PolygonPrism,
    step::Writer,
};
use multicad_model::{Compiler, GeometryBackend, Model, imports::read_step_file_with_backend};
fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn compile(s: &Session) -> Model {
    Compiler::default()
        .compile_with_backend(s.document(), s.revision(), GeometryBackend::Rust, None)
        .unwrap()
}
fn text(depth: f64) -> String {
    let q = |n| Rational::new(n, 1).unwrap();
    let p = PolygonPrism::with_circular_holes(
        vec![
            vec![[0, 0], [20, 0], [20, 16], [0, 16]]
                .into_iter()
                .map(|v| v.map(q))
                .collect(),
        ],
        vec![Circle::new([q(10), q(8)], q(1)).unwrap()],
        q(6),
    )
    .unwrap();
    let b = ChamferedBores::countersunk(
        &p,
        &[Countersink {
            hole: 0,
            end: false,
            radius_mm: 3.,
            depth_mm: depth,
        }],
        1e-7,
        1e-4,
    )
    .unwrap();
    let mut w = Writer::new();
    w.append_treated_polygon_prism(&p, &b, [0.; 3], [1., 0., 0.], [0., 1., 0.], false)
        .unwrap();
    w.finish().unwrap()
}
fn reopened(s: &Session) -> Session {
    let mut bytes = Vec::new();
    project_file::write(
        &ProjectFile::new(s.document().clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    Session::new(
        project_file::read(
            &mut bytes.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap()
        .project
        .document,
    )
    .unwrap()
}
#[test]
fn conical_import_replacement_history_and_offline_reopen_preserve_unchanged_faces_across_45_degrees()
 {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("seated-plate.step");
    std::fs::write(&path, text(1.25)).unwrap();
    let file = read_step_file_with_backend(&path, GeometryBackend::Rust).unwrap();
    let mut session = Session::new(Document::new("Conical seat import").unwrap()).unwrap();
    let import = prepare_step_import(&session, 0, actor(), &file).unwrap();
    let source = import.source;
    let part = import.model.outputs[0].clone();
    session.commit(import.model.prepared).unwrap();
    let first = compile(&session);
    let topo = &first.snapshot().shapes[&part.feature].topology;
    let fixed: Vec<_> = (0..6)
        .map(|ordinal| {
            GeometryIntent::capture(
                topo,
                TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal,
                },
            )
            .unwrap()
        })
        .collect();
    let changing: Vec<_> = (6..8)
        .map(|ordinal| {
            GeometryIntent::capture(
                topo,
                TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal,
                },
            )
            .unwrap()
        })
        .collect();
    for depth in [2., 2.5] {
        std::fs::write(&path, text(depth)).unwrap();
        let file = read_step_file_with_backend(&path, GeometryBackend::Rust).unwrap();
        let before = session.document().clone();
        let revision = session.revision();
        let (plan, update) =
            prepare_step_update(&session, revision, actor(), source, &file).unwrap();
        assert_eq!(plan.matched.len(), 1);
        assert!(plan.missing.is_empty());
        assert!(plan.added.is_empty());
        let preview = Compiler::default()
            .compile_with_backend(
                update.model.prepared.document(),
                revision + 1,
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        assert_eq!(session.document(), &before);
        session.commit(update.model.prepared).unwrap();
        let after = session.document().clone();
        let model = compile(&session);
        let shape = &model.snapshot().shapes[&part.feature];
        let expected = (320. - std::f64::consts::PI) * 6. - std::f64::consts::PI * depth * 10. / 3.;
        assert!((model.solid(part.body).unwrap().volume_mm3().unwrap() - expected).abs() < 1e-7);
        assert!((preview.solid(part.body).unwrap().volume_mm3().unwrap() - expected).abs() < 1e-7);
        for intent in &fixed {
            assert!(matches!(
                resolve(intent, &shape.topology, session.revision()),
                Resolution::Resolved { .. }
            ));
        }
        for intent in &changing {
            assert!(matches!(
                resolve(intent, &shape.topology, session.revision()),
                Resolution::Missing
            ));
        }
        session.undo(session.revision(), actor()).unwrap();
        assert_eq!(session.document(), &before);
        compile(&session);
        session.redo(session.revision(), actor()).unwrap();
        assert_eq!(session.document(), &after);
        compile(&session);
    }
    std::fs::remove_file(path).unwrap();
    let reopened = reopened(&session);
    let final_model = compile(&reopened);
    let before = compile(&session);
    assert_eq!(
        before.solid(part.body).unwrap().volume_mm3().unwrap(),
        final_model.solid(part.body).unwrap().volume_mm3().unwrap()
    );
    let out = dir.path().join("self-contained.step");
    final_model
        .solid(part.body)
        .unwrap()
        .write_step(&out)
        .unwrap();
    assert!(
        multicad_kernel::Solid::read_step(&out, Default::default())
            .unwrap()
            .volume_mm3()
            .unwrap()
            > 0.
    );
}
