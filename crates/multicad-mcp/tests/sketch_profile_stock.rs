//! CAD-385 over MCP: sketch-based profile stock is drawn from section
//! elements or given as its canonical constraint sketch, drives real member
//! solids, edits rebuild members in one undoable transaction, and invalid
//! sections are refused without changing the document.
use multicad_core::{
    commands::{Actor, Channel},
    document::Document,
    materials::{LibraryDefinition, ProfileShape, StockKind},
};
use multicad_mcp::{
    Headless, Workspace,
    tools::{Call, dispatch},
};
use serde_json::{Value, json};
use std::f64::consts::PI;
use std::sync::{Arc, atomic::AtomicBool};

fn call() -> Call {
    Call {
        actor: Actor {
            id: "mcp:sketch-profile".into(),
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
fn host() -> Headless {
    Headless::with_geometry_backend(
        Workspace::new(Document::new("Sketch profile MCP").unwrap()).unwrap(),
        None,
        false,
        multicad_model::GeometryBackend::Rust,
    )
}
fn stored_sketch(host: &Headless, stock: &Value) -> Value {
    let document = host.document();
    let id = serde_json::from_value(stock.clone()).unwrap();
    let entity = document.get(id).unwrap();
    let Some(LibraryDefinition::Stock(stock)) = multicad_core::materials::read(entity).unwrap()
    else {
        panic!("stock")
    };
    let StockKind::Profile {
        shape: ProfileShape::Sketch,
        section: Some(section),
        width_mm,
        height_mm,
        ..
    } = stock.specification
    else {
        panic!("sketch stock")
    };
    json!({ "sketch": section.sketch, "size": [width_mm, height_mm] })
}

#[test]
fn drawn_tube_section_creates_members_and_edits_rebuild_with_history() {
    let host = host();
    let run = |name: &str, arguments: Value| dispatch(&host, &call(), name, arguments);
    let stock = run(
        "create_profile_stock",
        json!({
            "expected_revision":0,"name":"Custom tube","material_name":"S235","shape":"sketch",
            "section":[
                {"kind":"rectangle","corner_mm":[0,0],"opposite_mm":[40,30]},
                {"kind":"circle","center_mm":[20,15],"radius_mm":8},
            ],
        }),
    )
    .unwrap();
    assert_eq!(
        stored_sketch(&host, &stock["stock"])["size"],
        json!([40., 30.])
    );
    run(
        "create_members",
        json!({
            "expected_revision":1,"stock":stock["stock"],
            "segments":[{"name":"Rail","start_mm":[0,0,0],"end_mm":[500,0,0]}],
        }),
    )
    .unwrap();
    let rail = |run: &dyn Fn(&str, Value) -> Result<Value, _>| {
        run("weldment_cut_list", json!({})).unwrap()["members"][0].clone()
    };
    let first = rail(&run);
    near(&first["volume_mm3"], (1200. - PI * 64.) * 500.);
    assert_eq!(first["shape"], "sketch");
    // Replace the section: the member keeps its Part and is rebuilt.
    run(
        "edit_profile_stock_section",
        json!({
            "expected_revision":2,"stock":stock["stock"],
            "section":[
                {"kind":"rectangle","corner_mm":[0,0],"opposite_mm":[40,30]},
                {"kind":"circle","center_mm":[20,15],"radius_mm":10},
            ],
        }),
    )
    .unwrap();
    let edited = rail(&run);
    assert_eq!(edited["part"], first["part"]);
    near(&edited["volume_mm3"], (1200. - PI * 100.) * 500.);
    run("undo", json!({"expected_revision":3})).unwrap();
    near(&rail(&run)["volume_mm3"], (1200. - PI * 64.) * 500.);
    run("redo", json!({"expected_revision":4})).unwrap();
    near(&rail(&run)["volume_mm3"], (1200. - PI * 100.) * 500.);
}

#[test]
fn lines_and_arcs_join_and_the_canonical_sketch_edits_dimensions() {
    let host = host();
    let run = |name: &str, arguments: Value| dispatch(&host, &call(), name, arguments);
    // A 40 × 20 slot from separately given pieces sharing endpoints.
    let stock = run(
        "create_profile_stock",
        json!({
            "expected_revision":0,"name":"Slot","material_name":"Aluminium","shape":"sketch",
            "section":[
                {"kind":"line","start_mm":[10,0],"end_mm":[30,0]},
                {"kind":"arc","center_mm":[30,10],"start_mm":[30,0],"end_mm":[30,20]},
                {"kind":"line","start_mm":[30,20],"end_mm":[10,20]},
                {"kind":"arc","center_mm":[10,10],"start_mm":[10,20],"end_mm":[10,0]},
            ],
        }),
    )
    .unwrap();
    assert_eq!(
        stored_sketch(&host, &stock["stock"])["size"],
        json!([40., 20.])
    );
    run(
        "create_members",
        json!({
            "expected_revision":1,"stock":stock["stock"],
            "segments":[{"name":"Slot rail","start_mm":[0,0,0],"end_mm":[300,0,0]}],
        }),
    )
    .unwrap();
    let list = run("weldment_cut_list", json!({})).unwrap();
    near(&list["members"][0]["volume_mm3"], (400. + PI * 100.) * 300.);
    // A rectangle carries its horizontal/vertical constraints: adding a
    // 60 mm length dimension to its bottom edge in the canonical sketch
    // re-solves to a 60 × 30 section.
    let rect = run(
        "create_profile_stock",
        json!({
            "expected_revision":2,"name":"Flat","material":stock["material"],"shape":"sketch",
            "section":[{"kind":"rectangle","corner_mm":[0,0],"opposite_mm":[40,30]}],
        }),
    )
    .unwrap();
    let stored = stored_sketch(&host, &rect["stock"]);
    let mut sketch = stored["sketch"].clone();
    let y = |id: &Value| {
        sketch["points"]
            .as_array()
            .unwrap()
            .iter()
            .find(|p| p["id"] == *id)
            .unwrap()["position_mm"][1]
            .as_f64()
            .unwrap()
    };
    let bottom = sketch["curves"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| y(&c["kind"]["start"]) == 0. && y(&c["kind"]["end"]) == 0.)
        .unwrap()["id"]
        .clone();
    sketch["constraints"].as_array_mut().unwrap().push(json!({
        "id": uuid(),
        "reference": false,
        "kind": {"kind":"length","line":bottom,"dimension":{"value":60.,"unit":"mm"}},
    }));
    run(
        "edit_profile_stock_section",
        json!({"expected_revision":3,"stock":rect["stock"],"section_sketch":sketch}),
    )
    .unwrap();
    let size = stored_sketch(&host, &rect["stock"])["size"].clone();
    near(&size[0], 60.);
    near(&size[1], 30.);
}

fn uuid() -> String {
    multicad_core::document::EntityId::new().to_string()
}

#[test]
fn invalid_sections_and_sizes_are_refused_without_changes() {
    let host = host();
    let run = |name: &str, arguments: Value| dispatch(&host, &call(), name, arguments);
    let before = host.document();
    let open = run(
        "create_profile_stock",
        json!({
            "expected_revision":0,"name":"Open","material_name":"S235","shape":"sketch",
            "section":[
                {"kind":"line","start_mm":[0,0],"end_mm":[10,0]},
                {"kind":"line","start_mm":[10,0],"end_mm":[10,10]},
            ],
        }),
    )
    .unwrap_err()
    .to_string();
    assert!(open.contains("open"), "{open}");
    let crossing = run(
        "create_profile_stock",
        json!({
            "expected_revision":0,"name":"Bow tie","material_name":"S235","shape":"sketch",
            "section":[{"kind":"polygon","vertices_mm":[[0,0],[10,10],[10,0],[0,10]]}],
        }),
    )
    .unwrap_err()
    .to_string();
    assert!(crossing.contains("crosses"), "{crossing}");
    let sized = run(
        "create_profile_stock",
        json!({
            "expected_revision":0,"name":"Sized","material_name":"S235","shape":"sketch",
            "width_mm":41,
            "section":[{"kind":"rectangle","corner_mm":[0,0],"opposite_mm":[40,30]}],
        }),
    )
    .unwrap_err()
    .to_string();
    assert!(sized.contains("derived from the section"), "{sized}");
    let walled = run(
        "create_profile_stock",
        json!({
            "expected_revision":0,"name":"Walled","material_name":"S235","shape":"sketch",
            "wall_mm":2,
            "section":[{"kind":"rectangle","corner_mm":[0,0],"opposite_mm":[40,30]}],
        }),
    );
    assert!(walled.is_err());
    let missing = run(
        "create_profile_stock",
        json!({"expected_revision":0,"name":"Empty","material_name":"S235","shape":"sketch"}),
    );
    assert!(missing.is_err());
    assert_eq!(host.document(), before);
    // A plain rectangular stock cannot take a section, and its section is
    // not editable.
    let plain = run(
        "create_profile_stock",
        json!({
            "expected_revision":0,"name":"SHS","material_name":"S235","shape":"rectangular",
            "width_mm":40,"height_mm":40,"wall_mm":2,
        }),
    )
    .unwrap();
    let error = run(
        "edit_profile_stock_section",
        json!({
            "expected_revision":1,"stock":plain["stock"],
            "section":[{"kind":"rectangle","corner_mm":[0,0],"opposite_mm":[40,30]}],
        }),
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("Only sketch-based"), "{error}");
}
