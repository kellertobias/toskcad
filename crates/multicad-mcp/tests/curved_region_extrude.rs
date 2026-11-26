//! CAD-406 / CAD-407: a circle sketched on a box face and clipped by that
//! face extrudes as Join and Cut into the box through MCP. The results are
//! exact Rust solids whose volumes follow from the clipped region's area.
use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::source::Source,
    document::{Document, EntityData},
    model::{self, PrimitiveParameters, PrimitiveRequest, length},
    references::TopologyKind,
    sketch::{self, CreationPlane, PlaneOptions},
};
use multicad_mcp::{
    Headless, Workspace,
    tools::{Call, dispatch},
};
use multicad_precision::PortableMath as _;
use serde_json::{Value, json};
use std::f64::consts::PI;
use std::sync::{Arc, atomic::AtomicBool};

fn call() -> Call {
    Call {
        actor: Actor {
            id: "mcp:curved-region-extrude".into(),
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

/// A box with a sketch on its top face holding one circle of `radius` about
/// the world point `center` (on the face plane). Returns the session, the
/// box primitive and the sketch.
fn box_with_face_circle(
    center: [f64; 2],
    radius: f64,
) -> (Headless, PrimitiveRequest, String, u64) {
    let mut session = Session::new(Document::new("Curved region").unwrap()).unwrap();
    let stock = PrimitiveRequest::new(
        "Stock",
        PrimitiveParameters::Box {
            size: SIZE.map(length),
            frame: Default::default(),
        },
    );
    let prepared = model::prepare_primitive(&session, 0, call().actor, &stock).unwrap();
    session.commit(prepared.prepared).unwrap();
    let native = multicad_model::Compiler::default()
        .compile(session.document(), session.revision())
        .unwrap();
    let top = native
        .solid(stock.feature)
        .unwrap()
        .faces()
        .unwrap()
        .into_iter()
        .find(|f| {
            f.plane
                .is_some_and(|(p, n)| (p[2] - SIZE[2]).abs() < 1e-7 && n[2] > 0.9)
        })
        .unwrap();
    let snapshot = &native.snapshot().shapes[&stock.feature].topology;
    let element = snapshot
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == top.index)
        .unwrap();
    let created = sketch::prepare_new(
        &session,
        session.revision(),
        call().actor,
        "Corner",
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
    (host, stock, sketch.to_string(), 0)
}

/// Area of a disk of radius `r` beyond a chord at distance `d` from its centre.
fn segment_area(r: f64, d: f64) -> f64 {
    r * r * (d / r).portable_acos() - d * (r * r - d * d).sqrt()
}

/// Area of the part of a disk with x > a and y > b relative to its centre.
fn corner_area(r: f64, a: f64, b: f64) -> f64 {
    let f = |x: f64| 0.5 * (x * (r * r - x * x).sqrt() + r * r * (x / r).portable_asin()) - b * x;
    f((r * r - b * b).sqrt()) - f(a)
}

#[test]
fn mcp_extrude_joins_and_cuts_a_face_clipped_circle_with_the_box() {
    let block = SIZE[0] * SIZE[1] * SIZE[2];
    for (center, radius, area) in [
        ([40., 30.], 10., PI * 25.),
        (
            [37., 27.],
            8.,
            PI * 64. - 2. * segment_area(8., 3.) + corner_area(8., 3., 3.),
        ),
        ([20., 27.], 8., PI * 64. - segment_area(8., 3.)),
        ([20., 15.], 6., PI * 36.),
    ] {
        let (host, _, sketch, _) = box_with_face_circle(center, radius);
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
            .unwrap_or_else(|| panic!("{center:?}: {regions:?}"));
        for (operation, distance, expected) in [
            ("join", 5., block + area * 5.),
            ("cut", -5., block - area * 5.),
            ("cut", -SIZE[2], block - area * SIZE[2]),
        ] {
            let (host, stock, sketch, revision) = box_with_face_circle(center, radius);
            let result = dispatch(
                &host,
                &call(),
                "extrude",
                json!({ "expected_revision": revision, "sketch": sketch, "regions": [index],
                    "distance": mm(distance), "operation": operation,
                    "target_part": stock.body }),
            )
            .unwrap_or_else(|e| panic!("{center:?} {operation} {distance}: {}", e.message));
            let volume = number(&result["bodies"][0]["volume_mm3"]);
            assert!(
                (volume - expected).abs() < 1e-6 * expected,
                "{center:?} {operation} {distance}: {volume} != {expected}"
            );
        }
    }
}
