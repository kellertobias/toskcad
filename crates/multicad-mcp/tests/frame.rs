//! CAD-439: external control builds frames through the same frame edits and
//! transaction as the UI, and reads the same measured cut list with mass.
use multicad_core::{
    commands::{Actor, Channel},
    document::Document,
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
            id: "mcp:frames".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(false)),
    }
}
fn near(v: &Value, expected: f64) {
    let v = v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"));
    assert!(
        (v - expected).abs() <= 1e-6 * expected.abs().max(1.),
        "{v} != {expected}"
    );
}
fn mm(value: f64) -> Value {
    json!({ "value": value, "unit": "mm" })
}
fn member<'a>(list: &'a Value, name: &str) -> &'a Value {
    list["members"]
        .as_array()
        .unwrap()
        .iter()
        .find(|m| m["name"] == name)
        .unwrap_or_else(|| panic!("no member {name}"))
}

#[test]
fn mcp_draws_a_mitred_rectangle_adds_legs_and_reads_the_cut_list_with_mass() {
    // A steel Material with a density, so the cut list reports mass.
    let mut session =
        multicad_core::commands::Session::new(Document::new("MCP frames").unwrap()).unwrap();
    let steel = multicad_core::document::Entity::new(
        "S235 steel",
        multicad_core::document::EntityData::Material {
            density_kg_m3: Some(7850.),
            color_rgba: [120, 124, 130, 255],
            section_hatch: None,
        },
    );
    let material = steel.id.to_string();
    let prepared = session
        .prepare(multicad_core::commands::Request {
            expected_revision: 0,
            actor: Actor::local(Channel::Ui),
            label: "Add steel".into(),
            edits: vec![multicad_core::document::Mutation::Insert(steel)],
        })
        .unwrap();
    session.commit(prepared).unwrap();
    let host = Headless::with_geometry_backend(
        Workspace::new(session.document().clone()).unwrap(),
        None,
        false,
        multicad_model::GeometryBackend::Rust,
    );
    let run = |name: &str, arguments: Value| dispatch(&host, &call(), name, arguments);
    run(
        "set_parameter",
        json!({ "expected_revision": 0, "name": "width", "expression": "1000", "unit": "mm" }),
    )
    .unwrap();
    let shs = run(
        "create_profile_stock",
        json!({ "expected_revision": 1, "name": "SHS 40×40×2", "material": material,
        "shape": "rectangular", "width_mm": 40, "height_mm": 40, "wall_mm": 2 }),
    )
    .unwrap();
    let revision = shs["revision"].clone();
    let shs = shs["stock"].clone();

    let frame = run(
        "create_frame",
        json!({ "expected_revision": revision, "name": "Table", "stock": shs, "side": "inner",
        "start": { "u": mm(0.), "v": mm(0.) },
        "segments": [
            { "length": { "expression": "=$width", "unit": "mm" }, "angle_deg": 0 },
            { "length": mm(600.), "angle_deg": 90 },
            { "length": { "expression": "=$width", "unit": "mm" }, "angle_deg": 180 },
        ],
        "close": true, "name_prefix": "Rail" }),
    )
    .unwrap();
    assert_eq!(frame["members"], 4);
    let list = run("weldment_cut_list", json!({})).unwrap();
    let area = 1600. - 36. * 36.;
    for (name, length) in [
        ("Rail 1", 1000.),
        ("Rail 2", 600.),
        ("Rail 3", 1000.),
        ("Rail 4", 600.),
    ] {
        let r = member(&list, name);
        near(&r["length_mm"], length);
        assert_eq!(r["start"]["treatment"], "miter", "{name}");
        near(&r["start"]["angle_deg"], 45.);
        assert_eq!(r["start"]["long_point"]["face"], "-x");
        assert_eq!(r["path"]["frame"], frame["frame"]);
        near(&r["volume_mm3"], area * (length - 40.));
        near(&r["mass_kg"], area * (length - 40.) * 1e-9 * 7850.);
    }

    // Four legs from the corners, 700 mm down along the normal.
    let corners: Vec<Value> = frame["frame_members"]
        .as_array()
        .unwrap()
        .iter()
        .map(|m| m["start"].clone())
        .collect();
    let operations: Vec<Value> = corners
        .iter()
        .map(|c| {
            json!({ "op": "add_chain", "from_node": c,
            "segments": [{ "length": mm(0.), "angle_deg": 0, "rise": mm(-700.) }] })
        })
        .collect();
    let edited = run(
        "edit_frame",
        json!({ "expected_revision": frame["revision"], "frame": frame["frame"], "stock": shs,
        "operations": operations, "name_prefix": "Leg" }),
    )
    .unwrap();
    assert_eq!(edited["members"], 8);
    let list = run("weldment_cut_list", json!({})).unwrap();
    for i in 5..=8 {
        let leg = member(&list, &format!("Leg {i}"));
        near(&leg["length_mm"], 720.);
        near(&leg["volume_mm3"], area * 720.);
    }
    near(&member(&list, "Rail 1")["length_mm"], 920.);

    // An override that cannot apply is refused with the reason; nothing
    // changes.
    let before = host.document();
    let e = run(
        "edit_frame",
        json!({ "expected_revision": edited["revision"], "frame": frame["frame"],
        "operations": [{ "op": "set_joint", "node": corners[0], "kind": "miter",
            "through": frame["frame_members"][0]["member"] }] }),
    )
    .unwrap_err();
    assert!(
        e.message.contains("A mitre has no through member"),
        "{}",
        e.message
    );
    assert_eq!(host.document(), before);

    // A leg corner built as rails mitred over the leg: the same transaction
    // as the Frame tool's corner dropdown.
    let operations: Vec<Value> = corners
        .iter()
        .map(|c| json!({ "op": "set_joint", "node": c, "corner": "mitre_on_leg" }))
        .collect();
    run(
        "edit_frame",
        json!({ "expected_revision": edited["revision"], "frame": frame["frame"],
        "operations": operations }),
    )
    .unwrap();
    let list = run("weldment_cut_list", json!({})).unwrap();
    near(&member(&list, "Rail 1")["length_mm"], 1000.);
    assert_eq!(member(&list, "Rail 1")["start"]["treatment"], "miter");
    for i in 5..=8 {
        near(&member(&list, &format!("Leg {i}"))["length_mm"], 680.);
    }
}
