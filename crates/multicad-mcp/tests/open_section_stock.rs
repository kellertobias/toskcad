//! CAD-434 (Q11): external control creates parametric channel, T and I/H
//! stock with a flange and a web thickness, builds members of them with the
//! same transactions as the UI, and is refused in plain words, without
//! changing the document, for an invalid web or a joint frames cannot build.
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
            id: "mcp:open-sections".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(false)),
    }
}
fn near(v: &Value, expected: f64) {
    let v = v.as_f64().unwrap();
    assert!(
        (v - expected).abs() <= 1e-7 * expected.abs().max(1.),
        "{v} != {expected}"
    );
}

#[test]
fn mcp_creates_channel_and_beam_stock_builds_members_and_refuses_what_frames_cannot() {
    let host = Headless::with_geometry_backend(
        Workspace::new(Document::new("Open sections").unwrap()).unwrap(),
        None,
        false,
        multicad_model::GeometryBackend::Rust,
    );
    let run = |name: &str, arguments: Value| dispatch(&host, &call(), name, arguments);
    let channel = run(
        "create_profile_stock",
        json!({ "expected_revision": 0, "name": "UPE 100", "material_name": "S235",
        "shape": "channel", "width_mm": 55, "height_mm": 100, "wall_mm": 7.5, "web_mm": 4.5 }),
    )
    .unwrap();
    let material = channel["material"].clone();
    let channel = channel["stock"].clone();

    // A web as wide as the section is refused with the reason; nothing is
    // created.
    let before = host.document();
    let e = run(
        "create_profile_stock",
        json!({ "expected_revision": 1, "name": "Bad T", "material": material,
        "shape": "tee", "width_mm": 50, "height_mm": 50, "wall_mm": 6, "web_mm": 50 }),
    )
    .unwrap_err();
    assert!(
        e.message
            .contains("The web must be thinner than the width (50 mm)"),
        "{}",
        e.message
    );
    // A web on a tube is refused too.
    let e = run(
        "create_profile_stock",
        json!({ "expected_revision": 1, "name": "Webbed tube", "material": material,
        "shape": "rectangular", "width_mm": 40, "height_mm": 40, "wall_mm": 2, "web_mm": 2 }),
    )
    .unwrap_err();
    assert!(
        e.message
            .contains("Only channel, T and I/H sections have a web thickness"),
        "{}",
        e.message
    );
    assert_eq!(host.document(), before);

    let beam = run(
        "create_profile_stock",
        json!({ "expected_revision": 1, "name": "IPE 200", "material": material,
        "shape": "beam", "width_mm": 100, "height_mm": 200, "wall_mm": 8.5, "web_mm": 5.6 }),
    )
    .unwrap()["stock"]
        .clone();
    let shs = run(
        "create_profile_stock",
        json!({ "expected_revision": 2, "name": "SHS 40×40×2", "material": material,
        "shape": "rectangular", "width_mm": 40, "height_mm": 40, "wall_mm": 2 }),
    )
    .unwrap()["stock"]
        .clone();
    let rail = run(
        "create_members",
        json!({ "expected_revision": 3, "stock": beam,
        "segments": [{ "name": "Beam rail", "start_mm": [0, 0, 0], "end_mm": [800, 0, 0] }] }),
    )
    .unwrap()["parts"][0]["part"]
        .clone();
    run(
        "create_members",
        json!({ "expected_revision": 4, "stock": channel,
        "segments": [{ "name": "Channel", "start_mm": [0, 500, 0], "end_mm": [600, 500, 0] }] }),
    )
    .unwrap();
    let post = run(
        "create_members",
        json!({ "expected_revision": 5, "stock": shs,
        "segments": [{ "name": "Post", "start_mm": [400, 0, 0], "end_mm": [400, 0, 600] }] }),
    )
    .unwrap()["parts"][0]["part"]
        .clone();

    let list = run("weldment_cut_list", json!({})).unwrap();
    let members = list["members"].as_array().unwrap();
    let channel = members.iter().find(|m| m["name"] == "Channel").unwrap();
    assert_eq!(channel["shape"], "channel");
    assert_eq!(channel["issues"], json!([]));
    near(&channel["length_mm"], 600.);
    near(&channel["volume_mm3"], (55. * 100. - 50.5 * 85.) * 600.);
    let rail_record = members.iter().find(|m| m["name"] == "Beam rail").unwrap();
    assert_eq!(rail_record["shape"], "beam");
    near(
        &rail_record["volume_mm3"],
        (2. * 100. * 8.5 + 5.6 * (200. - 17.)) * 800.,
    );

    // A post butted onto the beam running through is refused; nothing changes.
    let before = host.document();
    let e = run(
        "connect_members",
        json!({ "expected_revision": 6, "name": "Post to beam", "kind": "butt",
        "members": [post, rail] }),
    )
    .unwrap_err();
    assert!(
        e.message
            .contains("Butt joints with an I/H beam running through are not supported yet"),
        "{}",
        e.message
    );
    assert_eq!(host.document(), before);
}
