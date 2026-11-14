//! Native custom connection lifecycle through the real shared MCP dispatcher.
#[path = "../../multicad-model/tests/support/custom_screw_connections.rs"]
mod fixture;
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityId, Placement},
    model::*,
    project_file::{self, ProjectFile},
    symbols::{self, SymbolAuthorRequest, SymbolKind, SymbolPlaceRequest},
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
fn lifecycle(mode: ScrewConnectionMode) {
    lifecycle_with_pins(mode, false)
}
fn lifecycle_with_pins(mode: ScrewConnectionMode, pinned: bool) {
    let (original, mut request) = fixture::fixture(mode);
    if pinned {
        for row in &mut request.hardware {
            let ScrewConnectionHardwareSource::Generated { part } = &row.source else {
                panic!()
            };
            let part =
                MechanicalPartRequest::new("Pinned current material", part.parameters.clone());
            let mut author = Session::new(Document::new("Library material").unwrap()).unwrap();
            let p = prepare_mechanical_part(&author, 0, Actor::local(Channel::Mcp), &part).unwrap();
            author.commit(p.prepared).unwrap();
            let definition = SymbolAuthorRequest::new(
                "Physical pin",
                part.body,
                SymbolKind::Static,
                vec![],
                None,
            );
            let p = symbols::prepare_author(
                &author,
                author.revision(),
                Actor::local(Channel::Mcp),
                &definition,
            )
            .unwrap();
            author.commit(p.prepared).unwrap();
            let definition =
                symbols::read_definition(author.document(), definition.definition).unwrap();
            let placement = SymbolPlaceRequest::from_definition(
                EntityId::new(),
                &definition,
                "Hardware",
                Default::default(),
                Placement::default(),
            )
            .unwrap();
            row.instance = placement.instance;
            row.source = ScrewConnectionHardwareSource::Library {
                definition: Box::new(definition),
                placement: Box::new(placement),
            };
        }
    }
    let stock_features = request
        .receivers
        .iter()
        .map(|r| r.source.input.producer)
        .collect::<Vec<_>>();
    let host = Headless::new(Workspace::new(original.clone()).unwrap(), None, false);
    let cancelled = dispatch(
        &host,
        &call(true),
        "set_screw_connection",
        json!({"expected_revision":0,"request":request}),
    )
    .unwrap_err();
    assert_eq!(cancelled.code, "cancelled");
    assert_eq!(host.document(), original);
    let created = dispatch(
        &host,
        &call(false),
        "set_screw_connection",
        json!({"expected_revision":0,"request":request}),
    )
    .unwrap();
    assert_eq!(created["revision"], 1);
    assert_eq!(created["native_geometry_checked"], true);
    assert!(
        created["checks"][0]["actual_wood_engagement_mm"]
            .as_f64()
            .unwrap()
            >= 4. - 1e-7
    );
    assert_eq!(
        created["checks"][0]["native_component_collision_checked"],
        true
    );
    let inspected = dispatch(
        &host,
        &call(false),
        "inspect_screw_connection",
        json!({"connection":request.connection}),
    )
    .unwrap();
    assert_eq!(
        inspected["request"],
        serde_json::to_value(&request).unwrap()
    );
    let good = host.document();
    assert_eq!(
        dispatch(
            &host,
            &call(true),
            "set_screw_connection",
            json!({"expected_revision":1,"request":request})
        )
        .unwrap_err()
        .code,
        "cancelled"
    );
    assert_eq!(host.document(), good);
    let mut invalid = request.clone();
    invalid.parameters.minimum_engagement = length(100.);
    assert!(
        dispatch(
            &host,
            &call(false),
            "set_screw_connection",
            json!({"expected_revision":1,"request":invalid})
        )
        .is_err()
    );
    assert_eq!(host.document(), good);
    let mut incompatible = request.clone();
    incompatible.mode = ScrewConnectionMode::ThroughNut;
    assert!(
        dispatch(
            &host,
            &call(false),
            "set_screw_connection",
            json!({"expected_revision":1,"request":incompatible})
        )
        .is_err()
    );
    assert_eq!(host.document(), good);
    assert_eq!(
        dispatch(
            &host,
            &call(false),
            "set_screw_connection",
            json!({"expected_revision":0,"request":request})
        )
        .unwrap_err()
        .code,
        "stale_revision"
    );
    assert_eq!(host.document(), good);
    let pilot = request.receivers.last_mut().unwrap();
    let ConnectionReceiverExtent::Blind { depth, .. } = &mut pilot.source.extent else {
        panic!()
    };
    *depth = length(if mode == ScrewConnectionMode::WoodInsert {
        5.75
    } else {
        15.75
    });
    let edited = dispatch(
        &host,
        &call(false),
        "set_screw_connection",
        json!({"expected_revision":1,"request":request}),
    )
    .unwrap();
    assert_eq!(edited["revision"], 2);
    assert_eq!(
        read_screw_connection_request(&host.document(), request.connection).unwrap(),
        request
    );
    for id in stock_features {
        assert_eq!(host.document().get(id), original.get(id));
    }
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(host.document()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let reopened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    assert_eq!(
        read_screw_connection_request(&reopened, request.connection).unwrap(),
        request
    );
    let native = multicad_model::Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile(&reopened, 0)
    .unwrap();
    assert_eq!(native.snapshot().screw_connection_checks.len(), 1);
    assert!(
        native.snapshot().screw_connection_checks[0]
            .actual_wood_engagement_mm
            .unwrap()
            >= 4. - 1e-7
    );
    let removed = dispatch(
        &host,
        &call(false),
        "remove_screw_connection",
        json!({"expected_revision":2,"connection":request.connection}),
    )
    .unwrap();
    assert_eq!(removed["revision"], 3);
    if pinned {
        // Library definitions are reusable assets; removing a connection
        // removes its instances/drivers/cuts, retaining the imported pins.
        let mut expected = original
            .entities()
            .map(|e| (e.id, e.clone()))
            .collect::<std::collections::BTreeMap<_, _>>();
        for row in &request.hardware {
            let ScrewConnectionHardwareSource::Library { placement, .. } = &row.source else {
                panic!()
            };
            let entity = good.get(placement.definition).unwrap();
            expected.insert(entity.id, entity.clone());
        }
        assert_eq!(
            host.document()
                .entities()
                .map(|e| (e.id, e.clone()))
                .collect::<std::collections::BTreeMap<_, _>>(),
            expected
        );
    } else {
        assert_eq!(host.document(), original);
    }
}
#[test]
fn direct_wood_full_mcp_lifecycle() {
    lifecycle(ScrewConnectionMode::DirectWood)
}
#[test]
fn female_insert_full_mcp_lifecycle() {
    lifecycle(ScrewConnectionMode::WoodInsert)
}
#[test]
fn distinct_hanger_full_mcp_lifecycle() {
    lifecycle(ScrewConnectionMode::HangerBolt)
}
#[test]
fn pinned_hanger_nut_and_washer_full_mcp_lifecycle() {
    lifecycle_with_pins(ScrewConnectionMode::HangerBolt, true)
}

#[test]
fn hanger_step_refusal_is_atomic_for_document_and_destination() {
    let (original, request) = fixture::fixture(ScrewConnectionMode::HangerBolt);
    let dir = tempfile::tempdir().unwrap();
    let host = Headless::new(
        Workspace::new(original).unwrap(),
        Some(dir.path().into()),
        false,
    );
    dispatch(
        &host,
        &call(false),
        "set_screw_connection",
        json!({"expected_revision":0,"request":request}),
    )
    .unwrap();
    let row = request
        .hardware
        .iter()
        .find(|h| h.role == ScrewConnectionHardwareRole::HangerBolt)
        .unwrap();
    let ScrewConnectionHardwareSource::Generated { part } = &row.source else {
        panic!()
    };
    let before = host.document();
    for existing in [false, true] {
        let file_name = if existing {
            "existing.step"
        } else {
            "new.step"
        };
        let path = dir.path().join(file_name);
        if existing {
            std::fs::write(&path, b"Existing user data").unwrap();
        }
        let result = dispatch(
            &host,
            &call(false),
            "export_step",
            json!({"part":part.body,"file_name":file_name}),
        );
        assert!(
            result.is_err(),
            "unsupported Hanger exchange must not publish"
        );
        assert_eq!(host.document(), before);
        if existing {
            assert_eq!(std::fs::read(&path).unwrap(), b"Existing user data");
        } else {
            assert!(!path.exists());
        }
    }
}
