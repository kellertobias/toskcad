//! Canonical pipe creation uses the real Rust kernel transaction path.
use multicad_core::{
    commands::{Actor, Channel},
    document::Document,
    model,
    project_file::{self, ProjectFile},
};
use multicad_mcp::{
    Headless, Storage, Workspace,
    tools::{Call, dispatch},
};
use multicad_model::{Compiler, GeometryBackend};
use serde_json::{Value, json};
use std::sync::{Arc, atomic::AtomicBool};
fn call(cancel: bool) -> Call {
    Call {
        actor: Actor {
            id: "mcp:rust-pipe".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(cancel)),
    }
}
fn mm(v: f64) -> Value {
    json!({"value":v,"unit":"mm"})
}
fn run(h: &Headless, name: &str, args: Value) -> Value {
    dispatch(h, &call(false), name, args).unwrap()
}
struct Store;
impl Storage for Store {
    fn open(&mut self, path: &std::path::Path) -> Result<ProjectFile, String> {
        project_file::read(
            &mut std::fs::File::open(path).map_err(|e| e.to_string())?,
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .map(|r| r.project)
    }
    fn save(&mut self, path: &std::path::Path, p: &ProjectFile) -> Result<(), String> {
        project_file::write(
            p,
            &mut std::fs::File::create(path).map_err(|e| e.to_string())?,
            Default::default(),
        )
        .map_err(|e| e.to_string())
    }
}
fn intent() -> model::PipeThreadIntent {
    model::PipeThreadIntent {
        form: model::PipeThreadForm::R,
        size: model::PipeThreadSize::OneEighth,
        gauge: model::PipeGaugeConvention::Iso7Nominal1994,
    }
}
fn args() -> Value {
    let i = intent();
    json!({"expected_revision":0,"name":"Rust R 1/8","family":"bspt","fit":"bspt_nominal","pipe":i,"major_diameter":mm(i.base_major_mm()),"pitch":mm(i.preset().pitch_mm),"length":mm(i.default_length_mm()),"representation":"modeled"})
}
#[test]
fn rust_bspt_mcp_creation_is_physical_atomic_and_reopens() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("pipe.tcd");
    Store
        .save(
            &path,
            &ProjectFile::new(Document::new("Rust pipe").unwrap()).unwrap(),
        )
        .unwrap();
    let h = Headless::with_geometry_backend(
        Workspace::project(path.clone(), Box::new(Store)).unwrap(),
        Some(dir.path().into()),
        false,
        GeometryBackend::Rust,
    );
    let original = h.document();
    assert_eq!(
        dispatch(&h, &call(true), "set_thread", args())
            .unwrap_err()
            .code,
        "cancelled"
    );
    assert_eq!(h.document(), original);
    let created = run(&h, "set_thread", args());
    let feature = serde_json::from_value(created["thread"].clone()).unwrap();
    let inspected = run(&h, "inspect_thread", json!({"thread":feature}));
    assert_eq!(inspected["source"]["parameters"]["pipe"], json!(intent()));
    let model = Compiler {
        omit_display: false,
        ..Default::default()
    }
    .compile_with_backend(&h.document(), 1, GeometryBackend::Rust, None)
    .unwrap();
    let solid = model.solid(feature).unwrap();
    assert!(solid.rust_authored_axial_pipe_thread_source().is_some());
    assert!((solid.volume_mm3().unwrap() - 423.7556150148513).abs() < 1e-3);
    let exported = run(
        &h,
        "export_step",
        json!({"part":created["part"],"file_name":"pipe.step"}),
    );
    assert!(dir.path().join("pipe.step").is_file());
    assert!(!dir.path().join(".pipe.step.partial.step").exists());
    assert!(exported["readback_volume_mm3"].is_null());
    let decoded =
        &exported["step_verification"]["retained_members"][0]["decoded_volume_bounds_mm3"];
    assert!(decoded[0].as_f64().unwrap() > 400.);
    assert!(decoded[1].as_f64().unwrap() >= decoded[0].as_f64().unwrap());
    let initial = h.document();
    let bad = json!({"expected_revision":1,"thread":feature,"pipe":{"form":"rc","size":"1/8","gauge":"iso7_nominal1994"}});
    assert!(dispatch(&h, &call(false), "set_thread", bad).is_err());
    assert_eq!(h.document(), initial);
    run(
        &h,
        "set_thread",
        json!({"expected_revision":1,"thread":feature,"length":mm(6.7)}),
    );
    let edited = h.document();
    run(&h, "undo", json!({"expected_revision":2}));
    assert_eq!(h.document(), initial);
    run(&h, "redo", json!({"expected_revision":3}));
    assert_eq!(h.document(), edited);
    run(&h, "save_project", json!({"expected_revision":4}));
    drop(h);
    let reopened = Headless::with_geometry_backend(
        Workspace::project(path, Box::new(Store)).unwrap(),
        Some(dir.path().into()),
        false,
        GeometryBackend::Rust,
    );
    assert_eq!(reopened.document(), edited);
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(&reopened.document(), 0, GeometryBackend::Rust, None)
    .unwrap();
    assert!(
        model
            .solid(feature)
            .unwrap()
            .rust_authored_axial_pipe_thread_source()
            .is_some()
    );
}
#[test]
fn active_nominal_npt_stock_clips_keep_actual_material_and_closed_manufacturing() {
    let dir = tempfile::tempdir().unwrap();
    let h = Headless::with_geometry_backend(
        Workspace::new(Document::new("NPT clip gap").unwrap()).unwrap(),
        Some(dir.path().into()),
        false,
        GeometryBackend::Rust,
    );
    let i = model::NptThreadIntent {
        form: model::NptThreadForm::External,
        size: model::NptThreadSize::OneEighth,
        gauge: model::NptGaugeConvention::Nbs1959PlainFace,
        truncation: model::NptTruncationTarget::MidpointConeParallel,
    };
    let before = h.document();
    let arguments = json!({"expected_revision":0,"name":"Nominal NPT","family":"npt","fit":"npt_historical_midpoint","npt":i,"major_diameter":mm(i.base_major_mm()),"pitch":mm(i.preset().pitch_mm),"length":mm(i.default_length_mm())});
    assert!(dispatch(&h, &call(true), "set_thread", arguments.clone()).is_err());
    assert_eq!(h.document(), before);
    let created = run(&h, "set_thread", arguments);
    let feature = serde_json::from_value(created["thread"].clone()).unwrap();
    let compiled = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(&h.document(), 1, GeometryBackend::Rust, None)
    .unwrap();
    let solid = compiled.solid(feature).unwrap();
    let clipped = solid
        .rust_clipped_npt_source()
        .expect("Stock clipping must remain physical");
    assert_eq!(
        clipped.source().spec().major_radius_mm,
        i.base_major_mm() / 2.
    );
    assert_eq!(clipped.source().spec().length_mm, i.default_length_mm());
    let receipt = solid
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
    assert_eq!(receipt.removed_exact_zero_triangles, 0);
    let inspected = run(&h, "inspect_thread", json!({"thread":feature}));
    assert_eq!(inspected["source"]["parameters"]["npt"], json!(i));
    let before_export = h.document();
    let exported = run(
        &h,
        "export_step",
        json!({"part":created["part"],"file_name":"nominal-npt.step"}),
    );
    assert!(dir.path().join("nominal-npt.step").is_file());
    assert!(!dir.path().join(".nominal-npt.step.partial.step").exists());
    assert!(exported["readback_volume_mm3"].is_null());
    let decoded =
        &exported["step_verification"]["retained_members"][0]["decoded_volume_bounds_mm3"];
    assert!(decoded[0].as_f64().unwrap() > 470.);
    assert!(decoded[1].as_f64().unwrap() < 471.);
    assert_eq!(
        exported["step_verification"]["source_correspondence"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(h.document(), before_export);
}
