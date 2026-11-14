//! Real MCP export verification for canonical placed insert sources.
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
