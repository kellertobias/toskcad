//! Real two-stock managed numerical authoring via the shared MCP transaction.
#[path = "../../multicad-model/tests/support/clamped_insert.rs"]
mod fixture;
use multicad_core::{
    commands::{Actor, Channel},
    model::*,
    project_file::{self, ProjectFile},
};
use multicad_mcp::{
    Headless, Workspace,
    tools::{Call, dispatch},
};
use serde_json::json;
use std::sync::{Arc, atomic::AtomicBool};
fn call(cancel: bool) -> Call {
    Call {
        actor: Actor::local(Channel::Mcp),
        cancel: Arc::new(AtomicBool::new(cancel)),
    }
}
#[test]
fn real_clamped_insert_receipts_edit_cancel_and_reopen_both_hands() {
    for left in [false, true] {
        let (session, mut r) = fixture::request();
        let original = session.document().clone();
        for row in &mut r.hardware {
            let ScrewConnectionHardwareSource::Generated { part } = &mut row.source else {
                unreachable!()
            };
            match &mut part.parameters {
                MechanicalPartParameters::MachineScrew { thread, .. } => {
                    thread.hand = if left {
                        ThreadHand::Left
                    } else {
                        ThreadHand::Right
                    }
                }
                MechanicalPartParameters::CustomConnection { parameters } => {
                    let CustomConnectionHardwareParameters::WoodInsert {
                        outer_thread,
                        inner_thread,
                        ..
                    } = parameters.as_mut()
                    else {
                        unreachable!()
                    };
                    for t in [outer_thread, inner_thread] {
                        t.hand = if left {
                            ThreadHand::Left
                        } else {
                            ThreadHand::Right
                        };
                    }
                }
                _ => unreachable!(),
            }
        }
        let host = Headless::new(Workspace::new(original.clone()).unwrap(), None, false);
        assert_eq!(
            dispatch(
                &host,
                &call(true),
                "set_screw_connection",
                json!({"expected_revision":0,"request":r})
            )
            .unwrap_err()
            .code,
            "cancelled"
        );
        assert_eq!(host.document(), original);
        let created = dispatch(
            &host,
            &call(false),
            "set_screw_connection",
            json!({"expected_revision":0,"request":r}),
        )
        .unwrap();
        assert_eq!(created["revision"], 1);
        assert_eq!(created["native_geometry_checked"], true);
        for name in ["actual_head_bearing_gap_mm", "actual_insert_inset_mm"] {
            let bounds = created["checks"][0][name].as_array().unwrap();
            let lo = bounds[0].as_f64().unwrap();
            let hi = bounds[1].as_f64().unwrap();
            assert!(lo >= 0. && lo <= hi && hi <= 1e-7);
        }
        let good = host.document();
        dispatch(&host, &call(false), "undo", json!({"expected_revision":1})).unwrap();
        assert_eq!(host.document(), original);
        dispatch(&host, &call(false), "redo", json!({"expected_revision":2})).unwrap();
        assert_eq!(host.document(), good);
        let mut bad = r.clone();
        bad.parameters.blind_tip_clearance = length(10.);
        assert!(
            dispatch(
                &host,
                &call(false),
                "set_screw_connection",
                json!({"expected_revision":3,"request":bad})
            )
            .is_err()
        );
        assert_eq!(host.document(), good);
        r.receivers[1].source.extent = ConnectionReceiverExtent::Blind {
            depth: length(6.25),
            floor: length(1.),
        };
        let edited = dispatch(
            &host,
            &call(false),
            "set_screw_connection",
            json!({"expected_revision":3,"request":r}),
        )
        .unwrap();
        assert_eq!(edited["revision"], 4);
        let inspected = dispatch(
            &host,
            &call(false),
            "inspect_screw_connection",
            json!({"connection":r.connection}),
        )
        .unwrap();
        assert_eq!(inspected["request"], serde_json::to_value(&r).unwrap());
        let mut bytes = vec![];
        project_file::write(
            &ProjectFile::new(host.document()).unwrap(),
            &mut bytes,
            Default::default(),
        )
        .unwrap();
        let restored = project_file::read(
            &mut bytes.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap()
        .project
        .document;
        let model = multicad_model::Compiler::default()
            .compile(&restored, 0)
            .unwrap();
        assert!(
            model.snapshot().screw_connection_checks[0]
                .actual_insert_inset_mm
                .unwrap()[0]
                > 0.
        );
        let removed = dispatch(
            &host,
            &call(false),
            "remove_screw_connection",
            json!({"expected_revision":4,"connection":r.connection}),
        )
        .unwrap();
        assert_eq!(removed["revision"], 5);
        assert_eq!(host.document(), original);
    }
}
