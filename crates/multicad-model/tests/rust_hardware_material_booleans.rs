use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::*,
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, GeometryBackend};

fn actor() -> Actor {
    Actor::local(Channel::CommandLine)
}
fn compile(session: &Session) -> multicad_model::Model {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(
        session.document(),
        session.revision(),
        GeometryBackend::Rust,
        None,
    )
    .unwrap()
}

#[test]
fn actual_insert_and_hanger_material_intersections_rebuild_through_history_and_saved_source() {
    for parameters in [
        default_custom_connection_wood_insert(),
        default_custom_connection_hanger_bolt(),
    ] {
        let mut session =
            Session::new(Document::new("Hardware material Boolean").unwrap()).unwrap();
        let hardware = MechanicalPartRequest::new("Actual hardware", parameters);
        session
            .commit(
                prepare_mechanical_part(&session, session.revision(), actor(), &hardware)
                    .unwrap()
                    .prepared,
            )
            .unwrap();
        let model = compile(&session);
        let solid = model.solid(hardware.feature).unwrap();
        let (frame, center) = if let Some(data) = solid.rust_wood_insert_data() {
            let r = (data.source.inner().maximum_radius_mm()
                + data.source.outer().minimum_radius_mm())
                / 2.;
            (data.frame, [r, 0., data.source.full_length_mm() / 2.])
        } else {
            let data = solid.rust_hanger_data().unwrap();
            (
                data.frame,
                [
                    data.source.spec().plain_radius_mm / 2.,
                    0.,
                    -data.source.spec().plain_length_mm / 2.,
                ],
            )
        };
        let n = frame.normal();
        let world: [f64; 3] = std::array::from_fn(|i| {
            frame.origin[i] + frame.u[i] * center[0] + frame.v[i] * center[1] + n[i] * center[2]
        });
        let cube = PrimitiveRequest::new(
            "Material probe",
            PrimitiveParameters::Box {
                size: std::array::from_fn(|_| length(0.02)),
                frame: PrimitiveFrame {
                    origin: world.map(|v| length(v - 0.01)),
                    ..Default::default()
                },
            },
        );
        session
            .commit(
                prepare_primitive(&session, session.revision(), actor(), &cube)
                    .unwrap()
                    .prepared,
            )
            .unwrap();
        let before = session.document().clone();
        let intersection = BooleanRequest::new(
            Operation::Intersect,
            BodyInput {
                body: hardware.body,
                producer: hardware.feature,
            },
            vec![BodyInput {
                body: cube.body,
                producer: cube.feature,
            }],
        );
        let prepared =
            prepare_boolean(&session, session.revision(), actor(), &intersection).unwrap();
        // Validate the complete candidate before committing it, like Preview.
        let candidate = Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile_with_backend(
            prepared.prepared.document(),
            session.revision() + 1,
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
        assert!(
            (candidate
                .solid(intersection.feature)
                .unwrap()
                .volume_mm3()
                .unwrap()
                - 0.02_f64.powi(3))
            .abs()
                < 1e-15
        );
        assert_eq!(session.document(), &before);
        session.commit(prepared.prepared).unwrap();
        let applied = session.document().clone();
        let model = compile(&session);
        assert_eq!(
            model.snapshot().bodies[&hardware.body],
            intersection.feature
        );
        assert!(model.snapshot().consumed_bodies.contains(&cube.body));
        session.undo(session.revision(), actor()).unwrap();
        assert_eq!(session.document(), &before);
        compile(&session);
        session.redo(session.revision(), actor()).unwrap();
        assert_eq!(session.document(), &applied);
        compile(&session);
        let mut bytes = vec![];
        project_file::write(
            &ProjectFile::new(applied.clone()).unwrap(),
            &mut bytes,
            Default::default(),
        )
        .unwrap();
        let reopened = project_file::read(
            &mut std::io::Cursor::new(bytes),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap();
        assert_eq!(reopened.project.document, applied);
        let fresh = Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile_with_backend(&reopened.project.document, 0, GeometryBackend::Rust, None)
        .unwrap();
        assert!(
            (fresh
                .solid(intersection.feature)
                .unwrap()
                .volume_mm3()
                .unwrap()
                - 0.02_f64.powi(3))
            .abs()
                < 1e-15
        );
    }
}
