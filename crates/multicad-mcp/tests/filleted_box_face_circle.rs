//! CAD-424: a circle sketched on the flat top face of a box whose twelve
//! edges are all filleted extrudes as Join and Cut through MCP. The results
//! are exact Rust solids whose volumes follow the closed-form rounded box.
use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::source::Source,
    document::{Document, EntityData},
    model::{
        self, BodyInput, EdgeTreatmentParameters, EdgeTreatmentProfile, EdgeTreatmentRequest,
        PrimitiveParameters, PrimitiveRequest, length,
    },
    references::{AnalyticGeometry, GeometryIntent, TopologyKind},
    sketch::{self, CreationPlane, PlaneOptions},
};
use multicad_mcp::{
    Headless, Workspace,
    tools::{Call, dispatch},
};
use multicad_model::{Compiler, GeometryBackend};
use serde_json::{Value, json};
use std::f64::consts::PI;
use std::sync::{Arc, atomic::AtomicBool};

fn call() -> Call {
    Call {
        actor: Actor {
            id: "mcp:filleted-box-face-circle".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(false)),
    }
}

fn mm(v: f64) -> Value {
    json!({ "value": v, "unit": "mm" })
}

fn number(value: &Value) -> f64 {
    value.as_f64().unwrap_or_else(|| panic!("number: {value}"))
}

const SIZE: [f64; 3] = [40., 30., 20.];
const RADIUS: f64 = 3.;

fn rounded_volume() -> f64 {
    let i = SIZE.map(|s| s - 2. * RADIUS);
    i[0] * i[1] * i[2]
        + 2. * RADIUS * (i[0] * i[1] + i[1] * i[2] + i[0] * i[2])
        + PI * RADIUS * RADIUS * (i[0] + i[1] + i[2])
        + 4. / 3. * PI * RADIUS.powi(3)
}

/// A box with all edges filleted and a sketch on its flat top face holding
/// one circle. Returns the host, the body and the sketch.
fn rounded_box_with_face_circle(center: [f64; 2], radius: f64) -> (Headless, String, String) {
    let mut session = Session::new(Document::new("Rounded box circle").unwrap()).unwrap();
    let stock = PrimitiveRequest::new(
        "Stock",
        PrimitiveParameters::Box {
            size: SIZE.map(length),
            frame: Default::default(),
        },
    );
    let prepared = model::prepare_primitive(&session, 0, call().actor, &stock).unwrap();
    session.commit(prepared.prepared).unwrap();
    let compile = |session: &Session| {
        Compiler::default()
            .compile_with_backend(
                session.document(),
                session.revision(),
                GeometryBackend::Rust,
                None,
            )
            .unwrap()
    };
    let native = compile(&session);
    let topology = &native.snapshot().shapes[&stock.feature].topology;
    let edges: Vec<_> = topology
        .elements
        .iter()
        .filter(|e| {
            e.key.kind == TopologyKind::Edge
                && matches!(e.geometry, Some(AnalyticGeometry::Line { .. }))
        })
        .map(|e| GeometryIntent::capture(topology, e.key).unwrap())
        .collect();
    let fillet = EdgeTreatmentRequest::new(
        BodyInput {
            body: stock.body,
            producer: stock.feature,
        },
        EdgeTreatmentParameters {
            profile: EdgeTreatmentProfile::Fillet,
            size: length(RADIUS),
        },
        edges,
    );
    let prepared =
        model::prepare_edge_treatment(&session, session.revision(), call().actor, &fillet).unwrap();
    session.commit(prepared.prepared).unwrap();
    let native = compile(&session);
    let top = native
        .solid(fillet.feature)
        .unwrap()
        .faces()
        .unwrap()
        .into_iter()
        .find(|f| {
            f.plane
                .is_some_and(|(p, n)| (p[2] - SIZE[2]).abs() < 1e-7 && n[2] > 0.9)
        })
        .unwrap();
    let snapshot = &native.snapshot().shapes[&fillet.feature].topology;
    let element = snapshot
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == top.index)
        .unwrap();
    let created = sketch::prepare_new(
        &session,
        session.revision(),
        call().actor,
        "Face circle",
        CreationPlane::Face {
            snapshot,
            key: element.key,
            anchor: Some(&element.frames[0].origin),
        },
        PlaneOptions::default(),
    )
    .unwrap();
    let sketch = created.sketch;
    let frame = created.frame;
    session.commit(created.prepared).unwrap();
    let EntityData::Sketch { definition, .. } = &session.document().get(sketch).unwrap().data
    else {
        panic!()
    };
    let mut source = Source::from_definition(definition).unwrap();
    let c = source.add_point(frame.coordinates([center[0], center[1], SIZE[2]]));
    source.add_circle(c, radius);
    session
        .apply(Request {
            expected_revision: session.revision(),
            actor: call().actor,
            label: "Draw circle".into(),
            edits: source.replace(session.document(), sketch).unwrap(),
        })
        .unwrap();
    let host = Headless::new(
        Workspace::new(session.document().clone()).unwrap(),
        None,
        false,
    );
    (host, stock.body.to_string(), sketch.to_string())
}

#[test]
fn mcp_extrude_joins_and_cuts_a_circle_on_a_fully_filleted_box() {
    let area = PI * 36.;
    for (operation, distance, expected) in [
        ("join", 5., rounded_volume() + area * 5.),
        ("cut", -5., rounded_volume() - area * 5.),
        ("cut", -SIZE[2], rounded_volume() - area * SIZE[2]),
    ] {
        let (host, body, sketch) = rounded_box_with_face_circle([20., 15.], 6.);
        let inspected = dispatch(
            &host,
            &call(),
            "inspect_sketch",
            json!({ "sketch": sketch }),
        )
        .unwrap();
        let regions = inspected["regions"].as_array().unwrap();
        let index = regions
            .iter()
            .position(|r| (number(&r["area_mm2"]) - area).abs() < 1e-6)
            .unwrap_or_else(|| panic!("{regions:?}"));
        let result = dispatch(
            &host,
            &call(),
            "extrude",
            json!({ "expected_revision": 0, "sketch": sketch, "regions": [index],
                "distance": mm(distance), "operation": operation, "target_part": body }),
        )
        .unwrap_or_else(|e| panic!("{operation} {distance}: {}", e.message));
        let volume = number(&result["bodies"][0]["volume_mm3"]);
        assert!(
            (volume - expected).abs() < 1e-6 * expected,
            "{operation} {distance}: {volume} != {expected}"
        );
    }
}
