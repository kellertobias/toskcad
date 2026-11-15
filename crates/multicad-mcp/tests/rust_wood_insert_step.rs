//! Real MCP export verification for canonical placed insert sources.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::*,
};
use multicad_mcp::{
    Headless, Storage, Workspace,
    tools::{Call, dispatch},
};
use multicad_model::GeometryBackend;
use serde_json::json;
use std::sync::{Arc, atomic::AtomicBool};
struct Store;
impl Storage for Store {
    fn open(
        &mut self,
        path: &std::path::Path,
    ) -> Result<multicad_core::project_file::ProjectFile, String> {
        multicad_core::project_file::read(
            &mut std::fs::File::open(path).map_err(|e| e.to_string())?,
            Default::default(),
            &multicad_core::project_file::supported_capabilities(),
        )
        .map(|r| r.project)
    }
    fn save(
        &mut self,
        path: &std::path::Path,
        project: &multicad_core::project_file::ProjectFile,
    ) -> Result<(), String> {
        multicad_core::project_file::write(
            project,
            &mut std::fs::File::create(path).map_err(|e| e.to_string())?,
            Default::default(),
        )
        .map(|_| ())
    }
}
#[test]
fn canonical_insert_recipe_preserves_real_threads_placement_history_and_saved_source() {
    use multicad_core::project_file::ProjectFile;
    for hand in ["right", "left"] {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("insert.tcd");
        Store
            .save(
                &path,
                &ProjectFile::new(Document::new("Authored insert").unwrap()).unwrap(),
            )
            .unwrap();
        let host = Headless::with_geometry_backend(
            Workspace::project(path.clone(), Box::new(Store)).unwrap(),
            Some(directory.path().into()),
            false,
            GeometryBackend::Rust,
        );
        let call = |cancel| Call {
            actor: Actor::local(Channel::Mcp),
            cancel: Arc::new(AtomicBool::new(cancel)),
        };
        let run = |name: &str, args| dispatch(&host, &call(false), name, args).unwrap();
        let quantity = |v| json!({"value":v,"unit":"mm"});
        let mut recipe = json!({"kind":"wood_insert","outer_thread":{"hand":hand,"length":quantity(4.3)},"inner_thread":{"hand":hand,"length":quantity(4.3)}});
        if hand == "left" {
            recipe["origin"] = json!([quantity(11.), quantity(-7.), quantity(5.)]);
            recipe["axis_u"] = json!([0.8, 0.6, 0.]);
            recipe["axis_v"] = json!([-0.36, 0.48, 0.8]);
        }
        let initial = host.document();
        let args = json!({"expected_revision":0,"name":"Physical insert","recipe":recipe});
        assert_eq!(
            dispatch(&host, &call(true), "set_mechanical_part", args.clone())
                .unwrap_err()
                .code,
            "cancelled"
        );
        assert_eq!(host.document(), initial);
        let created = run("set_mechanical_part", args);
        let feature = serde_json::from_value(created["feature"].clone()).unwrap();
        let build = |doc: &Document, revision| {
            multicad_model::Compiler {
                omit_display: true,
                ..Default::default()
            }
            .compile_with_backend(doc, revision, GeometryBackend::Rust, None)
            .unwrap()
        };
        let original = host.document();
        let first = build(&original, 1);
        let first_solid = first.solid(feature).unwrap();
        let before = first_solid.volume_mm3().unwrap();
        assert!(first_solid.rust_wood_insert_data().is_some());
        let inspected = run("inspect_mechanical_part", json!({"feature":feature}));
        assert_eq!(inspected["resolved"]["connection_kind"], "wood_insert");
        let needs_context = dispatch(&host, &call(false), "set_mechanical_part", json!({"expected_revision":1,"feature":feature,"recipe":{"kind":"wood_insert","outer_thread":{"length":quantity(5.1)},"inner_thread":{"length":quantity(5.1)}}})).unwrap_err();
        assert!(needs_context.message.contains("normal engagement"));
        assert_eq!(host.document(), original);
        let edited = run(
            "set_mechanical_part",
            json!({"expected_revision":1,"feature":feature,"recipe":{"kind":"wood_insert","outer_thread":{"length":quantity(5.1)},"inner_thread":{"length":quantity(5.1),"metric_fit_engagement":{"kind":"intended_length","length":quantity(3.)}}}}),
        );
        assert_eq!(edited["part"], created["part"]);
        let edited_doc = host.document();
        let after = build(&edited_doc, 2);
        let after_solid = after.solid(feature).unwrap();
        assert!((after_solid.volume_mm3().unwrap() - before * 5.1 / 4.3).abs() < 1e-6);
        let mesh = after_solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert!(mesh.edge_pairs > 0 && mesh.vertex_links > 0);
        run("undo", json!({"expected_revision":2}));
        assert_eq!(host.document(), original);
        run("redo", json!({"expected_revision":3}));
        assert_eq!(host.document(), edited_doc);
        for args in [
            json!({"expected_revision":3,"feature":feature,"recipe":{"kind":"wood_insert"}}),
            json!({"expected_revision":4,"feature":feature,"recipe":{"kind":"hanger_bolt"}}),
            json!({"expected_revision":4,"feature":feature,"recipe":{"kind":"wood_insert","inner_thread":{"major_diameter":quantity(30.)}}}),
            json!({"expected_revision":4,"feature":feature,"recipe":{"kind":"wood_insert","inner_thread":{"length":quantity(9.)}}}),
            json!({"expected_revision":4,"feature":feature,"recipe":{"kind":"wood_insert","outer_thread":{"origin":[quantity(1.),quantity(0.),quantity(0.)]}}}),
            json!({"expected_revision":4,"feature":feature,"recipe":{"kind":"wood_insert","inner_thread":{"plane":"xz"}}}),
        ] {
            assert!(dispatch(&host, &call(false), "set_mechanical_part", args).is_err());
            assert_eq!(host.document(), edited_doc);
        }
        run("save_project", json!({"expected_revision":4}));
        let bytes = std::fs::read(&path).unwrap();
        let reopened = Headless::with_geometry_backend(
            Workspace::project(path.clone(), Box::new(Store)).unwrap(),
            Some(directory.path().into()),
            false,
            GeometryBackend::Rust,
        );
        assert_eq!(reopened.document(), edited_doc);
        let receipt = dispatch(
            &reopened,
            &call(false),
            "export_step",
            json!({"part":created["part"],"file_name":"authored.step"}),
        )
        .unwrap();
        assert_eq!(
            receipt["step_verification"]["method"],
            "decoded_boundary_enclosure"
        );
        assert!(receipt["readback_volume_mm3"].is_null());
        assert_eq!(std::fs::read(&path).unwrap(), bytes);
    }
}
#[test]
fn actual_mcp_export_keeps_mapped_insert_mass_provenance_and_document_unchanged() {
    let directory = tempfile::tempdir().unwrap();
    for left in [false, true] {
        let mut s = Session::new(Document::new("MCP placed insert").unwrap()).unwrap();
        let actor = Actor::local(Channel::Mcp);
        let mut parameters = default_custom_connection_wood_insert();
        let MechanicalPartParameters::CustomConnection { parameters: custom } = &mut parameters
        else {
            panic!()
        };
        let CustomConnectionHardwareParameters::WoodInsert {
            outer_thread,
            inner_thread,
            frame,
        } = custom.as_mut()
        else {
            panic!()
        };
        outer_thread.length = length(4.3);
        inner_thread.length = length(4.3);
        outer_thread.hand = if left {
            ThreadHand::Left
        } else {
            ThreadHand::Right
        };
        inner_thread.hand = outer_thread.hand;
        *frame = PrimitiveFrame {
            origin: [length(11.), length(-7.), length(5.)],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        };
        let request = MechanicalPartRequest::new("Actual export insert", parameters);
        s.commit(
            prepare_mechanical_part(&s, 0, actor.clone(), &request)
                .unwrap()
                .prepared,
        )
        .unwrap();
        let before = s.document().clone();
        let host = Headless::with_geometry_backend(
            Workspace::new(before.clone()).unwrap(),
            Some(directory.path().into()),
            false,
            GeometryBackend::Rust,
        );
        let call = |cancel| Call {
            actor: actor.clone(),
            cancel: Arc::new(AtomicBool::new(cancel)),
        };
        let name = if left { "left.step" } else { "right.step" };
        let args = json!({"part":request.body,"file_name":name});
        let cancelled = dispatch(&host, &call(true), "export_step", args.clone()).unwrap_err();
        assert_eq!(cancelled.code, "cancelled");
        assert!(!directory.path().join(name).exists());
        let result = dispatch(&host, &call(false), "export_step", args).unwrap();
        assert!(result["readback_volume_mm3"].is_null());
        let proof = &result["step_verification"];
        assert_eq!(proof["method"], "decoded_boundary_enclosure");
        assert_eq!(proof["retained_members"].as_array().unwrap().len(), 1);
        assert_eq!(proof["source_correspondence"].as_array().unwrap().len(), 1);
        let bounds = &proof["retained_members"][0]["decoded_volume_bounds_mm3"];
        assert!(
            bounds[0].as_f64().unwrap() > 0.
                && bounds[1].as_f64().unwrap() >= bounds[0].as_f64().unwrap()
        );
        assert!(
            proof["source_correspondence"][0]["maximum_position_error_mm"]
                .as_f64()
                .unwrap()
                <= 1e-7
        );
        let bytes = std::fs::read(directory.path().join(name)).unwrap();
        assert_eq!(result["bytes"].as_u64().unwrap(), bytes.len() as u64);
        assert!(String::from_utf8(bytes).unwrap().contains("MAPPED_ITEM"));
        assert!(
            !directory
                .path()
                .join(format!(".{name}.partial.step"))
                .exists()
        );
        assert_eq!(host.document(), before);
    }
}
