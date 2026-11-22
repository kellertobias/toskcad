use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{self, PrimitiveFrame, ThreadParameters, ThreadRepresentation, ThreadRequest, length},
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, GeometryBackend};
fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn compile(s: &Session) -> multicad_model::Model {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(s.document(), s.revision(), GeometryBackend::Rust, None)
    .unwrap()
}

#[test]
fn cosmetic_shaft_is_display_only_until_modeled_and_survives_history_and_reopen() {
    let mut s = Session::new(Document::new("Cosmetic thread").unwrap()).unwrap();
    let parameters = ThreadParameters {
        representation: ThreadRepresentation::CosmeticEnvelope,
        length: length(3.7),
        frame: PrimitiveFrame {
            origin: [length(10.), length(20.), length(30.)],
            ..Default::default()
        },
        ..Default::default()
    };
    let request = ThreadRequest::new("Cosmetic M8", parameters.clone());
    s.commit(
        model::prepare_thread(&s, 0, actor(), &request)
            .unwrap()
            .prepared,
    )
    .unwrap();
    let cosmetic = s.document().clone();
    let model = compile(&s);
    let solid = model.solid(request.feature).unwrap();
    assert!(solid.uses_rust_geometry());
    assert!((solid.volume_mm3().unwrap() - std::f64::consts::PI * 16. * 3.7).abs() < 1e-9);
    assert_eq!(solid.bounds().unwrap().min, [6., 16., 30.]);
    assert_eq!(
        model::read_thread(s.document(), request.feature)
            .unwrap()
            .parameters,
        parameters
    );
    assert!(
        model::manufacturing_limitation_in_history(s.document(), request.feature)
            .unwrap()
            .is_some()
    );
    let modeled = ThreadParameters {
        representation: ThreadRepresentation::Modeled,
        ..parameters.clone()
    };
    s.commit(
        model::prepare_thread_edit(&s, s.revision(), actor(), request.feature, &modeled)
            .unwrap()
            .prepared,
    )
    .unwrap();
    let finished = s.document().clone();
    assert!(
        compile(&s)
            .solid(request.feature)
            .unwrap()
            .volume_mm3()
            .unwrap()
            < solid.volume_mm3().unwrap()
    );
    assert!(
        model::manufacturing_limitation_in_history(s.document(), request.feature)
            .unwrap()
            .is_none()
    );
    s.undo(s.revision(), actor()).unwrap();
    assert_eq!(s.document(), &cosmetic);
    compile(&s);
    s.redo(s.revision(), actor()).unwrap();
    assert_eq!(s.document(), &finished);
    compile(&s);
    s.undo(s.revision(), actor()).unwrap();
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(s.document().clone()).unwrap(),
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
    let reopened = Session::new(reopened.document).unwrap();
    assert_eq!(
        model::read_thread(reopened.document(), request.feature)
            .unwrap()
            .parameters,
        parameters
    );
    assert!(
        model::manufacturing_limitation_in_history(reopened.document(), request.feature)
            .unwrap()
            .is_some()
    );
    assert!(
        (compile(&reopened)
            .solid(request.feature)
            .unwrap()
            .volume_mm3()
            .unwrap()
            - solid.volume_mm3().unwrap())
        .abs()
            < 1e-9
    );
}

#[test]
fn print_permission_checks_selected_producers_without_blocking_unrelated_cosmetic_parts() {
    use multicad_core::{
        exchange::{AidPolicy, PrintRequest, PrintSource},
        model::{PrimitiveParameters, PrimitiveRequest},
    };
    use multicad_model::export::Exporter;
    let mut session = Session::new(Document::new("Selected print permission").unwrap()).unwrap();
    let cosmetic = ThreadRequest::new(
        "Display-only M8",
        ThreadParameters {
            representation: ThreadRepresentation::CosmeticEnvelope,
            ..Default::default()
        },
    );
    session
        .commit(
            model::prepare_thread(&session, session.revision(), actor(), &cosmetic)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let printable = PrimitiveRequest::new(
        "Printable box",
        PrimitiveParameters::Box {
            size: [length(2.), length(3.), length(4.)],
            frame: PrimitiveFrame {
                origin: [length(20.), length(0.), length(0.)],
                ..Default::default()
            },
        },
    );
    session
        .commit(
            model::prepare_primitive(&session, session.revision(), actor(), &printable)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let unchanged = session.document().clone();
    let request = |parts| PrintRequest {
        source: PrintSource::Parts { parts },
        aids: AidPolicy::IncludeAll,
        quality: Default::default(),
        stl_per_part: true,
    };
    let exporter = Exporter::default();
    let denied = exporter
        .print_preview(
            session.document(),
            session.revision(),
            &request(vec![cosmetic.body]),
            &actor(),
            None,
        )
        .unwrap_err();
    assert!(denied.message.contains("cosmetic") && denied.message.contains("Modeled"));
    let allowed = exporter
        .print_preview(
            session.document(),
            session.revision(),
            &request(vec![printable.body]),
            &actor(),
            None,
        )
        .unwrap();
    assert!(allowed.can_write(), "{allowed:?}");
    assert_eq!(allowed.objects.len(), 1);
    assert_eq!(allowed.objects[0].part, printable.body);
    assert!((allowed.objects[0].volume_mm3 - 24.).abs() < 1e-9);
    assert_eq!(session.document(), &unchanged);
}
