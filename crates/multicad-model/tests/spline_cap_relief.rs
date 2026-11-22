use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::{Curve, CurveKind, Dimension, SplineBasis, source::Source},
    document::{Document, Entity, EntityData, EntityId, Mutation},
    model::{
        self, BodyInput, BooleanRequest, Direction, Extent, ExtrudeParameters, ExtrudeRequest,
        Operation, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, Termination, length,
    },
    parameters::BoundExpression,
    project_file::{self, ProjectFile},
    references::{GeometryIntent, Resolution, TopologyKind, resolve},
    sketch::Datum,
};
use multicad_kernel::{GeometryPolicy, NativeRegionQueries, Solid};
use multicad_model::{Compiler, Model};
fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}
fn compile(session: &Session) -> Model {
    Compiler::default()
        .compile(session.document(), session.revision())
        .unwrap()
}
fn fixture(operation: Operation) -> (Session, EntityId, EntityId, f64) {
    let (sketch, mut edits) = Source::create(Datum::Xy, "Counter columns").unwrap();
    edits.push(Mutation::Insert(Entity::new(
        "depth",
        EntityData::Parameter {
            expression: "1".into(),
            unit: "mm".into(),
            exposed: true,
            inputs: vec![],
        },
    )));
    let document = Document::new("Editable spline counters")
        .unwrap()
        .edited(&edits)
        .unwrap();
    let mut source = Source::read(&document, sketch).unwrap().0;
    source.add_rectangle([2., 2.], [9., 7.]).unwrap();
    for shift in [0., 3.] {
        let points = [
            [3. + shift, 4.],
            [4. + shift, 5.],
            [5. + shift, 4.],
            [4. + shift, 3.],
        ]
        .map(|p| source.add_point(p));
        let controls = [
            [3. + shift, 5.],
            [5. + shift, 5.],
            [5. + shift, 3.],
            [3. + shift, 3.],
        ]
        .map(|p| source.add_point(p));
        for i in 0..4 {
            source.sketch.curves.push(Curve {
                id: EntityId::new(),
                construction: false,
                kind: CurveKind::BSpline {
                    poles: vec![points[i], controls[i], points[(i + 1) % 4]],
                    basis: SplineBasis::open_uniform(2, 3).unwrap(),
                },
            });
        }
    }
    let mut session = Session::new(
        document
            .edited(&source.replace(&document, sketch).unwrap())
            .unwrap(),
    )
    .unwrap();
    let request = PrimitiveRequest::new(
        "Plate",
        PrimitiveParameters::Box {
            size: [20., 20., 5.].map(length),
            frame: PrimitiveFrame::default(),
        },
    );
    let prepared =
        model::prepare_primitive(&session, session.revision(), actor(), &request).unwrap();
    session.commit(prepared.prepared).unwrap();
    let stock = BodyInput {
        body: request.body,
        producer: request.feature,
    };
    let sketch_eval = multicad_solver::evaluate_document_sketch(
        session.document(),
        session.revision(),
        sketch,
        None,
    )
    .unwrap();
    let regions = sketch_eval
        .regions_with_queries(&NativeRegionQueries::analytic(GeometryPolicy::default()).unwrap())
        .unwrap();
    let region = regions.regions.iter().find(|r| r.holes.len() == 2).unwrap();
    let area = region.area_mm2;
    let depth = Dimension::Expression {
        expression: BoundExpression::bind(session.document(), "=$depth", "mm", None).unwrap(),
    };
    let extrusion = ExtrudeRequest::new(
        sketch,
        vec![region.intent()],
        ExtrudeParameters {
            direction: Direction {
                reverse: operation == Operation::Join,
                ..Default::default()
            },
            extent: Extent::OneSided {
                end: Termination::Distance { distance: depth },
            },
            ..Default::default()
        },
    );
    let prepared =
        model::prepare_extrude(&session, session.revision(), actor(), &extrusion).unwrap();
    session.commit(prepared.prepared).unwrap();
    let output = &extrusion.outputs[0];
    let boolean = BooleanRequest::new(
        operation,
        stock,
        vec![BodyInput {
            body: output.body,
            producer: output.feature,
        }],
    );
    let prepared = model::prepare_boolean(&session, session.revision(), actor(), &boolean).unwrap();
    session.commit(prepared.prepared).unwrap();
    (session, request.body, boolean.feature, area)
}
#[test]
fn saved_parameter_edits_rebuild_counter_columns_keep_face_intents_and_undo_redo() {
    let directory = tempfile::tempdir().unwrap();
    for operation in [Operation::Join, Operation::Cut] {
        let (mut session, body, feature, area) = fixture(operation);
        let sign = if operation == Operation::Join {
            1.
        } else {
            -1.
        };
        let original = compile(&session);
        let solid = original.solid(body).unwrap();
        assert!((solid.volume_mm3().unwrap() - (2000. + sign * area)).abs() < 1e-6);
        assert!(solid.rust_spline_cap_relief_data().is_some());
        let before = solid
            .topology_snapshot(feature, session.revision())
            .unwrap();
        let intents = before
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face)
            .map(|f| GeometryIntent::capture(&before, f.key).unwrap())
            .collect::<Vec<_>>();
        let initial = session.document().clone();
        let mut parameter = session
            .document()
            .entities()
            .find(|e| e.name == "depth")
            .unwrap()
            .clone();
        let EntityData::Parameter { expression, .. } = &mut parameter.data else {
            panic!()
        };
        *expression = "1.5".into();
        session
            .apply(Request {
                expected_revision: session.revision(),
                actor: actor(),
                label: "Edit retained column depth".into(),
                edits: vec![Mutation::Replace(parameter)],
            })
            .unwrap();
        let changed = session.document().clone();
        let rebuilt = compile(&session);
        let solid = rebuilt.solid(body).unwrap();
        assert!((solid.volume_mm3().unwrap() - (2000. + sign * 1.5 * area)).abs() < 1e-6);
        let topology = solid
            .topology_snapshot(feature, session.revision())
            .unwrap();
        for intent in &intents {
            assert!(matches!(
                resolve(intent, &topology, session.revision()),
                Resolution::Resolved { .. }
            ));
        }
        solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        let path = directory.path().join("saved-counter-columns.step");
        solid.write_step(&path).unwrap();
        Solid::read_step(&path, Default::default()).unwrap();
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
        .unwrap();
        let reopened = Session::new(reopened.project.document).unwrap();
        assert!(
            (compile(&reopened)
                .solid(body)
                .unwrap()
                .volume_mm3()
                .unwrap()
                - (2000. + sign * 1.5 * area))
                .abs()
                < 1e-6
        );
        session.undo(session.revision(), actor()).unwrap();
        assert_eq!(session.document(), &initial);
        assert!(
            (compile(&session).solid(body).unwrap().volume_mm3().unwrap() - (2000. + sign * area))
                .abs()
                < 1e-6
        );
        session.redo(session.revision(), actor()).unwrap();
        assert_eq!(session.document(), &changed);
        assert!(
            (compile(&session).solid(body).unwrap().volume_mm3().unwrap()
                - (2000. + sign * 1.5 * area))
                .abs()
                < 1e-6
        );
    }
}
