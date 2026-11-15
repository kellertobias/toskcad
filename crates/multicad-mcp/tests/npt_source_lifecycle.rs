//! Source/native lifecycle only; the unchanged rounded_threads NPT regression
//! additionally requires independent STEP export/readback and remains separate.
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
use serde_json::{Value, json};
use std::sync::{Arc, atomic::AtomicBool};
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
    fn save(&mut self, path: &std::path::Path, project: &ProjectFile) -> Result<(), String> {
        project_file::write(
            project,
            &mut std::fs::File::create(path).map_err(|e| e.to_string())?,
            Default::default(),
        )
        .map_err(|e| e.to_string())
    }
}
fn call(cancel: bool) -> Call {
    Call {
        actor: Actor {
            id: "mcp:clipped-npt-source".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(cancel)),
    }
}
fn run(host: &Headless, tool: &str, args: Value) -> Value {
    dispatch(host, &call(false), tool, args).unwrap()
}
fn mm(v: f64) -> Value {
    json!({"value":v,"unit":"mm"})
}
#[test]
fn canonical_clipped_npt_source_native_edit_cancel_history_and_reopen() {
    for size in [
        model::NptThreadSize::OneEighth,
        model::NptThreadSize::OneQuarter,
    ] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("clipped-source.tcd");
        Store
            .save(
                &path,
                &ProjectFile::new(Document::new("Clipped NPT source").unwrap()).unwrap(),
            )
            .unwrap();
        let host = Headless::new(
            Workspace::project(path.clone(), Box::new(Store)).unwrap(),
            Some(dir.path().into()),
            false,
        );
        let intent = model::NptThreadIntent {
            form: model::NptThreadForm::External,
            size,
            gauge: model::NptGaugeConvention::Nbs1959PlainFace,
            truncation: model::NptTruncationTarget::MidpointConeParallel,
        };
        let args = json!({"expected_revision":0,"name":"Clipped NPT source","family":"npt","fit":"npt_historical_midpoint","npt":intent,"major_diameter":mm(intent.base_major_mm()),"pitch":mm(intent.preset().pitch_mm),"length":mm(intent.default_length_mm())});
        let before = host.document();
        assert_eq!(
            dispatch(&host, &call(true), "set_thread", args.clone())
                .unwrap_err()
                .code,
            "cancelled"
        );
        assert_eq!(host.document(), before);
        let made = run(&host, "set_thread", args);
        let feature = made["thread"].clone();
        let original = host.document();
        assert!(
            dispatch(
                &host,
                &call(false),
                "set_thread",
                json!({"expected_revision":1,"thread":feature,"axial_allowance":mm(0.01)})
            )
            .is_err()
        );
        assert_eq!(host.document(), original);
        run(
            &host,
            "set_thread",
            json!({"expected_revision":1,"thread":feature,"radial_allowance":mm(0.001)}),
        );
        let edited = host.document();
        let info = run(&host, "inspect_thread", json!({"thread":feature}));
        assert_eq!(info["source"]["parameters"]["npt"], json!(intent));
        run(&host, "undo", json!({"expected_revision":2}));
        assert_eq!(host.document(), original);
        run(&host, "redo", json!({"expected_revision":3}));
        assert_eq!(host.document(), edited);
        run(&host, "save_project", json!({"expected_revision":4}));
        drop(host);
        let reopened = Headless::new(
            Workspace::project(path, Box::new(Store)).unwrap(),
            Some(dir.path().into()),
            false,
        );
        assert_eq!(reopened.document(), edited);
        let info = run(&reopened, "inspect_thread", json!({"thread":feature}));
        assert_eq!(info["source"]["parameters"]["npt"], json!(intent));
    }
}
