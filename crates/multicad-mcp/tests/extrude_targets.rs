//! CAD-405: MCP extrude without target_part discovers its Join/Cut bodies,
//! including the body whose face the sketch lies on. Native evaluation of the
//! final candidate decides every reported volume.
use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::{Curve, CurveKind, source::Source},
    document::{Document, EntityData, EntityId},
    model::{self, BodyInput, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, length},
    references::TopologyKind,
    sketch::{self, CreationPlane, PlaneOptions},
};
use multicad_mcp::{
    Headless, Workspace,
    tools::{Call, dispatch},
};
use serde_json::{Value, json};
use std::sync::{Arc, atomic::AtomicBool};

fn call() -> Call {
    Call {
        actor: Actor {
            id: "mcp:extrude-targets".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(false)),
    }
}
fn mm(v: f64) -> Value {
    json!({ "value": v, "unit": "mm" })
}
fn near(actual: &Value, expected: f64) {
    let actual = actual
        .as_f64()
        .unwrap_or_else(|| panic!("number: {actual}"));
    assert!((actual - expected).abs() < 1e-4, "{actual} != {expected}");
}
fn box_part(s: &mut Session, name: &str, origin: [f64; 3], size: [f64; 3]) -> BodyInput {
    let request = PrimitiveRequest::new(
        name,
        PrimitiveParameters::Box {
            size: size.map(length),
            frame: PrimitiveFrame {
                origin: origin.map(length),
                ..Default::default()
            },
        },
    );
    let prepared = model::prepare_primitive(s, s.revision(), call().actor, &request).unwrap();
    s.commit(prepared.prepared).unwrap();
    BodyInput {
        body: request.body,
        producer: request.feature,
    }
}

/// A 20 × 20 × 10 box with `draw` applied to a sketch on its top face.
fn face_sketch(draw: impl FnOnce(&mut Source, sketch::Frame)) -> (Session, BodyInput, EntityId) {
    let mut session = Session::new(Document::new("CAD-405 MCP").unwrap()).unwrap();
    let stock = box_part(&mut session, "Box", [0.; 3], [20., 20., 10.]);
    let native = multicad_model::Compiler::default()
        .compile(session.document(), session.revision())
        .unwrap();
    let top = native
        .solid(stock.producer)
        .unwrap()
        .faces()
        .unwrap()
        .into_iter()
        .find(|f| {
            f.plane
                .is_some_and(|(p, n)| (p[2] - 10.).abs() < 1e-7 && n[2] > 0.9)
        })
        .unwrap();
    let snapshot = &native.snapshot().shapes[&stock.producer].topology;
    let element = snapshot
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == top.index)
        .unwrap();
    let created = sketch::prepare_new(
        &session,
        session.revision(),
        call().actor,
        "Top face",
        CreationPlane::Face {
            snapshot,
            key: element.key,
            anchor: Some(&element.frames[0].origin),
        },
        PlaneOptions::default(),
    )
    .unwrap();
    let (sketch, frame) = (created.sketch, created.frame);
    session.commit(created.prepared).unwrap();
    let EntityData::Sketch { definition, .. } = &session.document().get(sketch).unwrap().data
    else {
        panic!()
    };
    let mut source = Source::from_definition(definition).unwrap();
    draw(&mut source, frame);
    session
        .apply(Request {
            expected_revision: session.revision(),
            actor: call().actor,
            label: "Draw on the face".into(),
            edits: source.replace(session.document(), sketch).unwrap(),
        })
        .unwrap();
    (session, stock, sketch)
}
fn host(session: &Session) -> Headless {
    Headless::new(
        Workspace::new(session.document().clone()).unwrap(),
        None,
        false,
    )
}
fn region(host: &Headless, sketch: EntityId, area: impl Fn(f64) -> bool) -> usize {
    let inspected = dispatch(host, &call(), "inspect_sketch", json!({ "sketch": sketch })).unwrap();
    inspected["regions"]
        .as_array()
        .unwrap()
        .iter()
        .position(|r| area(r["area_mm2"].as_f64().unwrap()))
        .expect("region")
}
fn extrude(host: &Headless, args: Value) -> Result<Value, String> {
    dispatch(host, &call(), "extrude", args).map_err(|e| format!("{e:?}"))
}

#[test]
fn join_and_cut_without_target_part_use_the_box_that_owns_the_sketch_face() {
    let (session, stock, sketch) = face_sketch(|source, frame| {
        source
            .add_rectangle(
                frame.coordinates([5., 5., 10.]),
                frame.coordinates([15., 15., 10.]),
            )
            .unwrap();
    });
    for (operation, distance, expected) in [("join", 5., 4500.), ("cut", -5., 3500.)] {
        let host = host(&session);
        let index = region(&host, sketch, |a| (a - 100.).abs() < 1e-6);
        let result = extrude(
            &host,
            json!({ "expected_revision": 0, "sketch": sketch,
                "regions": [index], "distance": mm(distance), "operation": operation }),
        )
        .unwrap();
        assert_eq!(result["created"][0]["part"], json!(stock.body));
        near(&result["bodies"][0]["volume_mm3"], expected);
        assert_eq!(result["bodies"][0]["solid_count"], 1);
    }
}

#[test]
fn corner_circle_join_no_longer_reports_no_affected_bodies() {
    let (session, stock, sketch) = face_sketch(|source, frame| {
        let center = source.add_point(frame.coordinates([18., 18., 10.]));
        source.sketch.curves.push(Curve {
            id: EntityId::new(),
            construction: false,
            kind: CurveKind::Circle {
                center,
                radius_mm: 6.,
            },
        });
    });
    let host = host(&session);
    let index = region(&host, sketch, |a| a < 100.);
    let before = host.document();
    let result = extrude(
        &host,
        json!({ "expected_revision": 0, "sketch": sketch,
            "regions": [index], "distance": mm(5.), "operation": "join" }),
    );
    match result {
        // The box is the discovered target; the joined body gains the prism.
        Ok(result) => {
            eprintln!("curved Join applied: {}", result["bodies"]);
            assert_eq!(result["created"][0]["part"], json!(stock.body));
            near(
                &result["bodies"][0]["volume_mm3"],
                4000. + 55.822171121641276 * 5.,
            );
        }
        // Target discovery succeeded; the exact curved-prism union itself is
        // CAD-406's kernel work. Nothing may change in the meantime.
        Err(error) => {
            eprintln!("curved Join target found; Boolean refused: {error}");
            assert!(!error.contains("No affected bodies"), "{error}");
            assert!(
                error.contains("no exact provider") || error.contains("certified positive convex"),
                "{error}"
            );
            assert_eq!(host.document(), before);
        }
    }
}

#[test]
fn several_touching_bodies_are_listed_and_intersect_still_needs_target_part() {
    let (mut session, _, sketch) = face_sketch(|source, frame| {
        source
            .add_rectangle(
                frame.coordinates([5., 5., 10.]),
                frame.coordinates([15., 15., 10.]),
            )
            .unwrap();
    });
    // The lid rests exactly on the 5 mm extrusion's end cap.
    let lid = box_part(&mut session, "Lid", [0., 0., 15.], [20., 20., 4.]);
    let host = host(&session);
    let index = region(&host, sketch, |a| (a - 100.).abs() < 1e-6);
    let before = host.document();
    let args = |operation: &str| {
        json!({ "expected_revision": 0, "sketch": sketch,
            "regions": [index], "distance": mm(5.), "operation": operation })
    };
    let error = extrude(&host, args("join")).unwrap_err();
    assert!(error.contains("Several bodies touch or overlap"), "{error}");
    assert!(error.contains("Box (") && error.contains(&format!("Lid ({})", lid.body)));
    let error = extrude(&host, args("intersect")).unwrap_err();
    assert!(error.contains("intersect requires target_part"), "{error}");
    // Contact alone has nothing to cut.
    let error = extrude(&host, args("cut")).unwrap_err();
    assert!(error.contains("No affected bodies"), "{error}");
    assert_eq!(host.document(), before);
    // An explicit target still chooses between them.
    let mut chosen = args("join");
    chosen["target_part"] = json!(lid.body);
    let joined = extrude(&host, chosen).unwrap();
    assert_eq!(joined["created"][0]["part"], json!(lid.body));
    near(&joined["bodies"][0]["volume_mm3"], 1600. + 500.);
}
