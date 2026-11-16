//! Normal MCP exports preserve explicit support-healing provenance.
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
fn call(cancel: bool) -> Call {
    Call {
        actor: Actor::local(Channel::Mcp),
        cancel: Arc::new(AtomicBool::new(cancel)),
    }
}
fn run(host: &Headless, name: &str, arguments: serde_json::Value) -> serde_json::Value {
    dispatch(host, &call(false), name, arguments).unwrap()
}
#[test]
fn hanger_recipe_creation_partial_edit_history_and_saved_reopen_use_the_actual_source() {
    use multicad_core::project_file::ProjectFile;
    for hand in ["right", "left"] {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("hanger.tcd");
        Store
            .save(
                &path,
                &ProjectFile::new(Document::new("MCP authored hanger").unwrap()).unwrap(),
            )
            .unwrap();
        let host = Headless::with_geometry_backend(
            Workspace::project(path.clone(), Box::new(Store)).unwrap(),
            Some(directory.path().into()),
            false,
            GeometryBackend::Rust,
        );
        let profiles = run(&host, "mechanical_part_profiles", json!({}));
        assert!(
            profiles["supported_geometry"]
                .as_array()
                .unwrap()
                .iter()
                .any(|v| v == "hanger_bolt")
        );
        let initial = host.document();
        let args = json!({"expected_revision":0,"name":"Authored hanger","recipe":{"kind":"hanger_bolt","wood_thread":{"hand":hand},"metric_thread":{"hand":hand}}});
        assert_eq!(
            dispatch(&host, &call(true), "set_mechanical_part", args.clone())
                .unwrap_err()
                .code,
            "cancelled"
        );
        assert_eq!(host.document(), initial);
        let created = run(&host, "set_mechanical_part", args);
        let feature = serde_json::from_value(created["feature"].clone()).unwrap();
        let build = |doc: &Document, revision| {
            multicad_model::Compiler {
                omit_display: true,
                ..Default::default()
            }
            .compile_with_backend(doc, revision, GeometryBackend::Rust, None)
            .unwrap()
        };
        let first_doc = host.document();
        let first = build(&first_doc, 1);
        let first_solid = first.solid(feature).unwrap();
        let before = first_solid.volume_mm3().unwrap();
        let data = first_solid.rust_hanger_data().unwrap();
        assert_eq!(data.source.spec().plain_length_mm, 2.);
        assert_eq!(data.source.wood().spec().left_hand, hand == "left");
        assert_eq!(data.source.metric().spec().left_hand, hand == "left");
        let edited = run(
            &host,
            "set_mechanical_part",
            json!({"expected_revision":1,"feature":feature,"recipe":{"kind":"hanger_bolt","plain_length":{"value":3.3,"unit":"mm"}}}),
        );
        assert_eq!(edited["part"], created["part"]);
        let edited_doc = host.document();
        let after = build(&edited_doc, 2);
        let after_solid = after.solid(feature).unwrap();
        let data = after_solid.rust_hanger_data().unwrap();
        assert_eq!(data.source.spec().plain_length_mm, 3.3);
        assert_eq!(data.source.wood().spec().left_hand, hand == "left");
        assert!(
            (after_solid.volume_mm3().unwrap()
                - before
                - std::f64::consts::PI * 1.5f64.powi(2) * 1.3)
                .abs()
                < 1e-6
        );
        run(&host, "undo", json!({"expected_revision":2}));
        assert_eq!(host.document(), first_doc);
        run(&host, "redo", json!({"expected_revision":3}));
        assert_eq!(host.document(), edited_doc);
        let before_invalid = host.document();
        for args in [
            json!({"expected_revision":3,"feature":feature,"recipe":{"kind":"hanger_bolt"}}),
            json!({"expected_revision":4,"feature":feature,"recipe":{"kind":"washer"}}),
            json!({"expected_revision":4,"feature":feature,"recipe":{"kind":"hanger_bolt","plain_length":{"value":-1.,"unit":"mm"}}}),
            json!({"expected_revision":4,"feature":feature,"recipe":{"kind":"hanger_bolt","plain_length":{"value":2.,"unit":"deg"}}}),
        ] {
            assert!(dispatch(&host, &call(false), "set_mechanical_part", args).is_err());
            assert_eq!(host.document(), before_invalid);
        }
        for thread in ["wood_thread", "metric_thread"] {
            for bad in [
                json!({"origin":[{"value":1.,"unit":"mm"},{"value":0.,"unit":"mm"},{"value":0.,"unit":"mm"}]}),
                json!({"plane":"yz"}),
                json!({"axis_u":[1.,0.,0.]}),
                json!({"axis_v":[0.,1.,0.]}),
            ] {
                let mut recipe = json!({"kind":"hanger_bolt"});
                recipe[thread] = bad;
                let error = dispatch(
                    &host,
                    &call(false),
                    "set_mechanical_part",
                    json!({"expected_revision":4,"feature":feature,"recipe":recipe}),
                )
                .unwrap_err();
                assert!(error.message.contains("nested thread frames"));
                assert_eq!(host.document(), before_invalid);
            }
        }
        run(&host, "save_project", json!({"expected_revision":4}));
        let bytes = std::fs::read(&path).unwrap();
        let reopened = Headless::with_geometry_backend(
            Workspace::project(path.clone(), Box::new(Store)).unwrap(),
            Some(directory.path().into()),
            false,
            GeometryBackend::Rust,
        );
        assert_eq!(reopened.document(), edited_doc);
        let proof = run(
            &reopened,
            "export_step",
            json!({"part":created["part"],"file_name":"authored.step"}),
        );
        assert_eq!(
            proof["step_verification"]["method"],
            "decoded_support_healing"
        );
        assert_eq!(std::fs::read(&path).unwrap(), bytes);
    }
}
#[test]
fn actual_mcp_hanger_export_reports_healed_mass_and_preserves_source_and_existing_output() {
    let directory = tempfile::tempdir().unwrap();
    for left in [false, true] {
        let mut session = Session::new(Document::new("Actual hanger export").unwrap()).unwrap();
        let actor = Actor::local(Channel::Mcp);
        let mut parameters = default_custom_connection_hanger_bolt();
        let MechanicalPartParameters::CustomConnection { parameters: custom } = &mut parameters
        else {
            panic!()
        };
        let CustomConnectionHardwareParameters::HangerBolt {
            wood_thread,
            metric_thread,
            ..
        } = custom.as_mut()
        else {
            panic!()
        };
        let hand = if left {
            ThreadHand::Left
        } else {
            ThreadHand::Right
        };
        wood_thread.hand = hand;
        metric_thread.hand = hand;
        let request = MechanicalPartRequest::new("Physical hanger", parameters);
        session
            .commit(
                prepare_mechanical_part(&session, 0, actor.clone(), &request)
                    .unwrap()
                    .prepared,
            )
            .unwrap();
        let before = session.document().clone();
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
        let path = directory.path().join(name);
        std::fs::write(&path, "previous output").unwrap();
        let args = json!({"part":request.body,"file_name":name});
        assert_eq!(
            dispatch(&host, &call(true), "export_step", args.clone())
                .unwrap_err()
                .code,
            "cancelled"
        );
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "previous output");
        let result = dispatch(&host, &call(false), "export_step", args).unwrap();
        assert!(result["readback_volume_mm3"].is_null());
        let proof = &result["step_verification"];
        assert_eq!(proof["method"], "decoded_support_healing");
        assert!(proof["retained_members"].is_null());
        assert!(proof["source_known_members"].is_null());
        let healed = proof["healed_members"].as_array().unwrap();
        assert_eq!(healed.len(), 1);
        assert!(healed[0]["healed_volume_bounds_mm3"][0].as_f64().unwrap() > 0.);
        assert!(
            healed[0]["source_to_healed_position_error_mm"]
                .as_f64()
                .unwrap()
                <= 1e-7
        );
        assert!(healed[0]["matched_entities"].as_u64().unwrap() > 100);
        assert_eq!(
            result["bytes"].as_u64().unwrap(),
            std::fs::metadata(&path).unwrap().len()
        );
        assert_eq!(host.document(), before);
        assert!(
            !directory
                .path()
                .join(format!(".{name}.partial.step"))
                .exists()
        );
    }
}
