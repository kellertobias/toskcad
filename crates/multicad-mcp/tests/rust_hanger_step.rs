//! Normal MCP exports preserve explicit support-healing provenance.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::*,
};
use multicad_mcp::{
    Headless, Workspace,
    tools::{Call, dispatch},
};
use multicad_model::GeometryBackend;
use serde_json::json;
use std::sync::{Arc, atomic::AtomicBool};
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
