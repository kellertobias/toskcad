//! External control of Rust double-helical gears and pairs.
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

fn call(cancelled: bool) -> Call {
    Call {
        actor: Actor {
            id: "mcp:double-helical".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(cancelled)),
    }
}

fn scalar(value: f64, unit: &str) -> Value {
    json!({"value": value, "unit": unit})
}

fn run(h: &Headless, name: &str, args: Value) -> Value {
    dispatch(h, &call(false), name, args).unwrap_or_else(|e| panic!("{name}: {e:?}"))
}

fn member(gap: f64, hand: &str) -> Value {
    json!({"family": "double_helical", "helix_angle": scalar(15., "deg"), "hand": hand,
        "face_width": scalar(10., "mm"), "center_gap": scalar(gap, "mm"),
        "bore_diameter": scalar(8., "mm")})
}

#[test]
fn rust_pair_evaluates_exports_edits_and_refuses_invalid_gaps() {
    let directory = tempfile::tempdir().unwrap();
    let h = Headless::with_geometry_backend(
        Workspace::new(Document::new("Rust double helical").unwrap()).unwrap(),
        Some(directory.path().to_path_buf()),
        false,
        multicad_model::GeometryBackend::Rust,
    );
    let created = run(
        &h,
        "set_gear_pair",
        json!({"expected_revision": 0, "name": "Double helical pair",
            "first": member(3., "right"), "second": member(3., "left")}),
    );
    let pair = created["pair"].clone();
    let evaluated = run(&h, "evaluate_model", json!({}));
    let check = &evaluated["gear_pairs"][0];
    assert_eq!(check["current_pose_checked"], true, "{evaluated}");
    assert_eq!(check["cycle_drive_angles_rad"].as_array().unwrap().len(), 9);
    let gear: multicad_core::document::EntityId =
        serde_json::from_value(created["members"][0]["gear"].clone()).unwrap();
    let gap = |h: &Headless| {
        multicad_core::model::read_gear(&h.document(), gear)
            .unwrap()
            .parameters
            .center_gap
    };
    assert_eq!(gap(&h), multicad_core::model::length(3.));

    let exported = run(
        &h,
        "export_step",
        json!({"part": created["members"][0]["part"], "file_name": "pinion.step"}),
    );
    let volume = exported["readback_volume_mm3"].as_f64().unwrap();
    assert!(volume > 0.);
    let text = std::fs::read_to_string(directory.path().join("pinion.step")).unwrap();
    // The core and both shoulder rows are real STEP faces.
    assert!(text.matches("=PLANE(").count() > 40);

    // Widening both gaps removes tooth material and keeps the pair proved.
    let edited = run(
        &h,
        "set_gear_pair",
        json!({"expected_revision": 1, "pair": pair,
            "first": {"center_gap": scalar(4., "mm")}, "second": {"center_gap": scalar(4., "mm")}}),
    );
    assert_eq!(edited["members"], created["members"]);
    let widened = run(&h, "evaluate_model", json!({}));
    assert_eq!(widened["gear_pairs"][0]["current_pose_checked"], true);
    let after = run(
        &h,
        "export_step",
        json!({"part": created["members"][0]["part"], "file_name": "widened.step"}),
    );
    assert!(after["readback_volume_mm3"].as_f64().unwrap() < volume);

    let before = h.document();
    for args in [
        json!({"expected_revision": 2, "pair": pair, "first": {"center_gap": scalar(10., "mm")},
            "second": {"center_gap": scalar(10., "mm")}}),
        json!({"expected_revision": 2, "pair": pair, "first": {"center_gap": scalar(-1., "mm")}}),
        json!({"expected_revision": 1, "pair": pair, "first": {"center_gap": scalar(2., "mm")}}),
    ] {
        assert!(dispatch(&h, &call(false), "set_gear_pair", args).is_err());
        assert_eq!(h.document(), before);
    }
    assert!(
        dispatch(
            &h,
            &call(true),
            "set_gear_pair",
            json!({"expected_revision": 2, "pair": pair, "first": {"center_gap": scalar(2., "mm")}})
        )
        .is_err()
    );
    assert_eq!(h.document(), before);
    run(&h, "undo", json!({"expected_revision": 2}));
    assert_eq!(gap(&h), multicad_core::model::length(3.));
    run(&h, "redo", json!({"expected_revision": 3}));
    assert_eq!(h.document(), before);
}
