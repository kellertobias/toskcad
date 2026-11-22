//! Shared normal MCP transactions compile the same actual combined material.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{self, PrimitiveParameters, PrimitiveRequest, length},
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
        actor: Actor::local(Channel::Mcp),
        cancel: Arc::new(AtomicBool::new(cancel)),
    }
}
fn mm(v: f64) -> Value {
    json!({"value":v,"unit":"mm"})
}
fn run(h: &Headless, name: &str, args: Value) -> Value {
    dispatch(h, &call(false), name, args.clone()).unwrap_or_else(|e| panic!("{name} {args}: {e:?}"))
}
struct Store;
impl Storage for Store {
    fn open(&mut self, path: &std::path::Path) -> Result<ProjectFile, String> {
        project_file::read(
            &mut std::fs::File::open(path).map_err(|e| e.to_string())?,
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .map(|p| p.project)
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
fn ends(c: f64, w: f64) -> Value {
    json!({"start":{"chamfer45":{"radial_depth":mm(c)}},"end":{"cylindrical_relief":{"length":mm(w),"radial_depth":mm(0.125)}}})
}
#[test]
fn normal_mcp_combined_internal_treatments_preview_edit_reopen_and_refuse_step_atomically() {
    for (fit, hand) in ["nominal_basic", "metric_6h_internal"]
        .into_iter()
        .flat_map(|fit| ["right", "left"].map(|hand| (fit, hand)))
    {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("combined.tcd");
        let mut source =
            Session::new(Document::new("MCP actual treated receiver").unwrap()).unwrap();
        let stock = PrimitiveRequest::new(
            "Cylinder stock",
            PrimitiveParameters::Cylinder {
                radius: length(2.5),
                height: length(4.),
                frame: Default::default(),
            },
        );
        let p = model::prepare_primitive(&source, 0, call(false).actor, &stock).unwrap();
        source.commit(p.prepared).unwrap();
        Store
            .save(&path, &ProjectFile::new(source.document().clone()).unwrap())
            .unwrap();
        let open = || {
            Headless::with_geometry_backend(
                Workspace::project(path.clone(), Box::new(Store)).unwrap(),
                Some(dir.path().into()),
                false,
                GeometryBackend::Rust,
            )
        };
        let host = open();
        let original = host.document();
        let catalog = run(&host, "thread_profiles", json!({}));
        assert_eq!(catalog["end_treatments"]["standalone_thread_only"], false);
        let scope = &catalog["active_geometry_support"]["combined_internal_treatments"];
        assert_eq!(scope["starts"], json!([1]));
        assert_eq!(
            scope["fits"],
            json!(["nominal_basic", "metric_6h_internal"])
        );
        assert_eq!(scope["step"], false);
        assert_eq!(scope["manufacturing_mesh"], true);
        let args = json!({"expected_revision":0,"name":"M3 combined internal","mode":"cut_internal","target_part":stock.body,"family":"iso_metric","fit":fit,"major_diameter":mm(3.),"pitch":mm(0.5),"length":mm(4.),"minimum_wall":mm(0.5),"hand":hand,"end_treatments":ends(0.5,0.5)});
        assert_eq!(
            dispatch(&host, &call(true), "set_thread", args.clone())
                .unwrap_err()
                .code,
            "cancelled"
        );
        assert_eq!(host.document(), original);
        for (key, value) in [
            ("length", mm(3.)),
            ("origin", json!([mm(0.1), mm(0.), mm(0.)])),
            ("minimum_wall", mm(1.)),
            ("starts", json!(2)),
        ] {
            let mut invalid = args.clone();
            invalid[key] = value;
            assert!(dispatch(&host, &call(false), "set_thread", invalid).is_err());
            assert_eq!(host.document(), original);
        }
        let created = run(&host, "set_thread", args);
        let feature = serde_json::from_value(created["thread"].clone()).unwrap();
        let initial = host.document();
        let edit = json!({"expected_revision":1,"thread":created["thread"],"end_treatments":ends(0.625,0.625)});
        assert_eq!(
            dispatch(&host, &call(true), "set_thread", edit.clone())
                .unwrap_err()
                .code,
            "cancelled"
        );
        let mut stale = edit;
        stale["expected_revision"] = json!(0);
        assert_eq!(
            dispatch(&host, &call(false), "set_thread", stale)
                .unwrap_err()
                .code,
            "stale_revision"
        );
        assert_eq!(host.document(), initial);

        let preview = run(&host, "evaluate_model", json!({}));
        assert_eq!(preview["parts"][0]["feature"], created["thread"]);
        assert_eq!(preview["parts"][0]["solid_count"], 1);
        let compiled = Compiler::default()
            .compile_with_backend(&initial, 1, GeometryBackend::Rust, None)
            .unwrap();
        let solid = compiled.solid(feature).unwrap();
        let data = solid.rust_internal_thread_treatment_data().unwrap();
        assert_eq!(
            data.material.cavity().source_thread().spec().left_hand,
            hand == "left"
        );
        assert_eq!(data.material.cavity().full_length_mm(), 4.);
        assert!(!compiled.snapshot().shapes[&feature].mesh.indices.is_empty());
        let receipt = solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert_eq!(receipt.edge_pairs as usize * 2, receipt.mesh.indices.len());
        let dest = dir.path().join("unsupported.step");
        std::fs::write(&dest, b"preserve destination").unwrap();
        assert!(
            dispatch(
                &host,
                &call(false),
                "export_step",
                json!({"part":created["part"],"file_name":"unsupported.step"})
            )
            .is_err()
        );
        assert_eq!(std::fs::read(dest).unwrap(), b"preserve destination");
        assert!(
            std::fs::read_dir(dir.path()).unwrap().all(|e| !e
                .unwrap()
                .file_name()
                .to_string_lossy()
                .contains("partial"))
        );
        assert_eq!(host.document(), initial);
        run(
            &host,
            "set_thread",
            json!({"expected_revision":1,"thread":created["thread"],"end_treatments":ends(0.625,0.625)}),
        );
        let edited = host.document();
        run(&host, "undo", json!({"expected_revision":2}));
        assert_eq!(host.document(), initial);
        run(&host, "redo", json!({"expected_revision":3}));
        assert_eq!(host.document(), edited);
        run(&host, "save_project", json!({"expected_revision":4}));
        let reopened = open();
        assert_eq!(reopened.document(), edited);
        let rebuilt = Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile_with_backend(&reopened.document(), 0, GeometryBackend::Rust, None)
        .unwrap();
        let data = rebuilt
            .solid(feature)
            .unwrap()
            .rust_internal_thread_treatment_data()
            .unwrap();
        assert_eq!(data.material.cavity().intent().entry_depth_mm, Some(0.625));
        assert_eq!(data.material.cavity().terminal_station_mm(), Some(3.375));
    }
}
