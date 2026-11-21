//! CAD-365: an MCP panel without an explicit grain choice defaults to Any,
//! explicit choices are kept, and an invalid grain is refused atomically.
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
            id: "mcp:grain".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(false)),
    }
}
fn mm(v: f64) -> Value {
    json!({ "value": v, "unit": "mm" })
}

#[test]
fn mcp_panels_default_to_any_grain_keep_explicit_grain_and_refuse_invalid_values() {
    let host = Headless::new(
        Workspace::new(Document::new("Grain").unwrap()).unwrap(),
        None,
        false,
    );
    let run = |name: &str, args: Value| dispatch(&host, &call(), name, args);
    let stock = run(
        "create_panel_stock",
        json!({ "expected_revision": 0, "name": "Birch 18", "material_name": "Birch",
        "thickness_mm": 18, "sheet_length_mm": 2500, "sheet_width_mm": 1250 }),
    )
    .unwrap()["stock"]
        .clone();
    let panel = |revision: u64, name: &str, grain: Option<&str>| {
        let mut args = json!({ "expected_revision": revision, "name": name, "stock": stock,
            "plane": "xy", "origin_mm": [0, 0, revision as f64 * 50.],
            "length": mm(600.), "width": mm(400.) });
        if let Some(grain) = grain {
            args["grain"] = json!(grain);
        }
        run("create_panel", args)
    };
    panel(1, "Default", None).unwrap();
    panel(2, "Length", Some("u")).unwrap();
    panel(3, "Width", Some("v")).unwrap();
    let before = host.document();
    assert!(panel(4, "Bad", Some("diagonal")).is_err());
    assert_eq!(host.document(), before, "refused atomically");
    let grains: Vec<_> = multicad_core::model::panels(&before)
        .map(|e| {
            let part = multicad_core::model::read_panel(&before, e.id).unwrap();
            (
                before.get(part.output.body).unwrap().name.clone(),
                part.grain,
            )
        })
        .collect();
    use multicad_core::model::GrainAxis;
    for (name, grain) in [
        ("Default", GrainAxis::Free),
        ("Length", GrainAxis::U),
        ("Width", GrainAxis::V),
    ] {
        assert!(
            grains.contains(&(name.to_owned(), grain)),
            "{name}: {grains:?}"
        );
    }
}
