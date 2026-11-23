//! CAD-364 over MCP: read and edit a generic symbol's fixed plan with the same
//! typed update transaction the Inspector's plan editor applies.
use multicad_core::{
    commands::{Actor, Channel, Session},
    constraint_sketch::source::Source,
    document::{Document, EntityId},
    model::{self, PrimitiveParameters, PrimitiveRequest, length},
    project_file::{self, ProjectFile},
    sketch::Datum,
    symbols::{self, PlanAlignment, SymbolAuthorRequest, SymbolKind, SymbolPlan},
};
use multicad_mcp::{
    Headless, Storage, Workspace,
    tools::{Call, dispatch},
};
use serde_json::{Value, json};
use std::{
    path::Path,
    sync::{Arc, atomic::AtomicBool},
};
fn actor() -> Actor {
    Actor {
        id: "mcp:symbol-plans".into(),
        channel: Channel::Mcp,
    }
}
fn call() -> Call {
    Call {
        actor: actor(),
        cancel: Arc::new(AtomicBool::new(false)),
    }
}
fn run(h: &Headless, name: &str, args: Value) -> Result<Value, String> {
    dispatch(h, &call(), name, args).map_err(|e| format!("{e:?}"))
}
struct Store;
impl Storage for Store {
    fn open(&mut self, path: &Path) -> Result<ProjectFile, String> {
        project_file::read(
            &mut std::fs::File::open(path).map_err(|e| e.to_string())?,
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .map(|r| r.project)
        .map_err(|e| e.to_string())
    }
    fn save(&mut self, path: &Path, project: &ProjectFile) -> Result<(), String> {
        project_file::write(
            project,
            &mut std::fs::File::create(path).map_err(|e| e.to_string())?,
            Default::default(),
        )
        .map_err(|e| e.to_string())
    }
}
/// A library file with one 600 × 400 × 750 table whose plan is a separate
/// 600 × 400 rectangle sketch.
fn library(dir: &Path) -> EntityId {
    let (sketch, edits) = Source::create(Datum::Xy, "Table plan symbol").unwrap();
    let document = Document::new("Table library")
        .unwrap()
        .edited(&edits)
        .unwrap();
    let mut source = Source::read(&document, sketch).unwrap().0;
    source.add_rectangle([0., 0.], [600., 400.]).unwrap();
    let document = document
        .edited(&source.replace(&document, sketch).unwrap())
        .unwrap();
    let mut s = Session::new(document).unwrap();
    let table = PrimitiveRequest::new(
        "Table",
        PrimitiveParameters::Box {
            size: [length(600.), length(400.), length(750.)],
            frame: Default::default(),
        },
    );
    s.commit(
        model::prepare_primitive(&s, s.revision(), actor(), &table)
            .unwrap()
            .prepared,
    )
    .unwrap();
    let author = SymbolAuthorRequest::new(
        "Table",
        table.body,
        SymbolKind::Static,
        vec![],
        Some(SymbolPlan::new(sketch)),
    );
    let prepared = symbols::prepare_author(&s, s.revision(), actor(), &author).unwrap();
    s.commit(prepared.prepared).unwrap();
    Store
        .save(
            &dir.join("tables.tcd"),
            &ProjectFile::new(s.document().clone()).unwrap(),
        )
        .unwrap();
    sketch
}

#[test]
fn mcp_reads_and_edits_symbol_plans_as_one_undoable_revision_with_typed_refusals() {
    let dir = tempfile::tempdir().unwrap();
    let template_sketch = library(dir.path());
    let h = Headless::new(
        Workspace::project(dir.path().join("consumer.tcd"), Box::new(Store)).unwrap(),
        Some(dir.path().into()),
        false,
    );
    let placed = run(
        &h,
        "place_library_symbol",
        json!({"expected_revision":0,"file_name":"tables.tcd","name":"Table 1","overrides":{},
            "placement":{"translation":[{"value":1000.,"unit":"mm"},{"value":0.,"unit":"mm"},{"value":0.,"unit":"mm"}]}}),
    )
    .unwrap();
    let instance = placed["instance"].as_str().unwrap().to_string();
    let revision = placed["revision"].as_u64().unwrap();
    let read = run(&h, "symbol_plan", json!({"symbol": instance})).unwrap();
    assert_eq!(read["plan"]["sketch"], json!(template_sketch));
    assert_eq!(read["plan"]["rotation_deg"], json!(0.));
    assert_eq!(read["template_sketches"].as_array().unwrap().len(), 1);
    assert_eq!(read["instances"].as_array().unwrap().len(), 1);
    let definition = EntityId::parse(read["definition"].as_str().unwrap()).unwrap();

    // Refusals leave the document unchanged.
    let before = h.document();
    for (plan, reason) in [
        (
            json!({"sketch": EntityId::new()}),
            "captured in this symbol",
        ),
        (
            json!({"sketch": template_sketch, "offset_mm":[1e9, 0.]}),
            "Plan alignment",
        ),
    ] {
        let error = run(
            &h,
            "edit_symbol_plan",
            json!({"expected_revision":revision,"symbol":instance,"plan":plan}),
        )
        .unwrap_err();
        assert!(error.contains(reason), "{error}");
        assert_eq!(h.document(), before);
    }
    let stale = run(
        &h,
        "edit_symbol_plan",
        json!({"expected_revision":revision + 7,"symbol":instance,"plan":{"sketch":template_sketch,"rotation_deg":90.}}),
    )
    .unwrap_err();
    assert!(stale.contains("stale_revision"), "{stale}");

    let edited = run(
        &h,
        "edit_symbol_plan",
        json!({"expected_revision":revision,"symbol":instance,
            "plan":{"sketch":template_sketch,"offset_mm":[0., -500.],"rotation_deg":90.}}),
    )
    .unwrap();
    assert_eq!(edited["native_geometry_checked"], json!(true));
    assert_ne!(edited["old_pin"], edited["new_pin"]);
    let stored = symbols::read_definition(&h.document(), definition).unwrap();
    assert_eq!(stored.revision, "2");
    assert_eq!(
        stored.plan.unwrap().alignment,
        PlanAlignment {
            offset_mm: [0., -500.],
            rotation_deg: 90.
        }
    );
    let after = h.document();
    let reread = run(&h, "symbol_plan", json!({"symbol": definition})).unwrap();
    assert_eq!(reread["plan"]["offset_mm"], json!([0., -500.]));
    // One attributable, undoable transaction.
    let revision = edited["revision"].as_u64().unwrap();
    run(&h, "undo", json!({"expected_revision": revision})).unwrap();
    assert_eq!(h.document(), before);
    run(&h, "redo", json!({"expected_revision": revision + 1})).unwrap();
    assert_eq!(h.document(), after);
    // Removing the plan is explicit too.
    let removed = run(
        &h,
        "edit_symbol_plan",
        json!({"expected_revision": revision + 2, "symbol": instance, "plan": null}),
    )
    .unwrap();
    assert!(removed["revision"].as_u64().is_some());
    assert!(
        symbols::read_definition(&h.document(), definition)
            .unwrap()
            .plan
            .is_none()
    );
}
