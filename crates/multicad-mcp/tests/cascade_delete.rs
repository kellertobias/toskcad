//! CAD-426: MCP delete refuses to cascade unless asked, lists the dependents,
//! and with `cascade: true` removes them in one Undo step.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityData, EntityId},
    model::*,
    references::{GeometryIntent, Origin, TopologyKind},
};
use multicad_mcp::{
    Headless, Workspace,
    tools::{Call, dispatch},
};
use serde_json::json;
use std::{
    collections::BTreeSet,
    sync::{Arc, atomic::AtomicBool},
};

fn call() -> Call {
    Call {
        actor: Actor {
            id: "mcp:cascade-delete".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(false)),
    }
}

/// A box with a fillet and a chamfer after it, plus an unrelated box.
fn document() -> (Document, BodyInput, EntityId, EntityId, BodyInput) {
    let actor = Actor::local(Channel::Ui);
    let mut session = Session::new(Document::new("MCP cascade").unwrap()).unwrap();
    let boxed = |session: &mut Session, name: &str| {
        let request = PrimitiveRequest::new(
            name,
            PrimitiveParameters::Box {
                size: [length(30.), length(20.), length(10.)],
                frame: Default::default(),
            },
        );
        let prepared =
            prepare_primitive(session, session.revision(), actor.clone(), &request).unwrap();
        session.commit(prepared.prepared).unwrap();
        BodyInput {
            body: request.body,
            producer: request.feature,
        }
    };
    let other = boxed(&mut session, "Spacer");
    let part = boxed(&mut session, "Bracket");
    let mut producer = part.producer;
    let mut treatments = vec![];
    for profile in [EdgeTreatmentProfile::Fillet, EdgeTreatmentProfile::Chamfer] {
        let request = EdgeTreatmentRequest::new(
            BodyInput {
                body: part.body,
                producer,
            },
            EdgeTreatmentParameters {
                profile,
                size: length(1.),
            },
            vec![GeometryIntent {
                id: EntityId::new(),
                kind: TopologyKind::Edge,
                origins: BTreeSet::from([Origin {
                    source: producer,
                    role: "selected-source-edge".into(),
                }]),
                boundary: None,
                attachment: None,
            }],
        );
        let prepared =
            prepare_edge_treatment(&session, session.revision(), actor.clone(), &request).unwrap();
        session.commit(prepared.prepared).unwrap();
        producer = request.feature;
        treatments.push(request.feature);
    }
    (
        session.document().clone(),
        part,
        treatments[0],
        treatments[1],
        other,
    )
}

#[test]
fn mcp_delete_refuses_without_cascade_and_cascades_in_one_undo_step_when_asked() {
    let (document, part, fillet, chamfer, other) = document();
    let host = Headless::new(Workspace::new(document.clone()).unwrap(), None, false);
    let run = |name: &str, arguments| dispatch(&host, &call(), name, arguments);

    let refused = run(
        "delete_entity",
        json!({ "expected_revision": 0, "id": part.body.to_string() }),
    )
    .unwrap_err();
    assert_eq!(refused.code, "has_dependents", "{}", refused.message);
    assert!(refused.message.contains("cascade"), "{}", refused.message);
    for id in [fillet, chamfer] {
        assert!(refused.involved.contains(&id), "{:?}", refused.involved);
    }
    // An explicit false is the same refusal; neither changes anything.
    let explicit = run(
        "delete_entity",
        json!({ "expected_revision": 0, "id": part.body.to_string(), "cascade": false }),
    )
    .unwrap_err();
    assert_eq!(explicit.code, "has_dependents");
    assert_eq!(host.document(), document);

    let deleted = run(
        "delete_entity",
        json!({ "expected_revision": 0, "id": part.body.to_string(), "cascade": true }),
    )
    .unwrap();
    assert_eq!(deleted["revision"], 1, "one transaction");
    let listed: Vec<_> = deleted["dependents"]
        .as_array()
        .unwrap()
        .iter()
        .map(|d| d["id"].as_str().unwrap().to_owned())
        .collect();
    assert_eq!(listed, [fillet.to_string(), chamfer.to_string()]);
    let after = host.document();
    for id in [part.body, part.producer, fillet, chamfer] {
        assert!(after.get(id).is_none());
    }
    for id in [other.body, other.producer] {
        assert_eq!(after.get(id), document.get(id));
    }

    run("undo", json!({ "expected_revision": 1 })).unwrap();
    assert_eq!(host.document(), document, "one Undo step restores all");

    // A leaf feature has no dependents and needs no option.
    let leaf = run(
        "delete_entity",
        json!({ "expected_revision": 2, "id": chamfer.to_string() }),
    )
    .unwrap();
    assert_eq!(leaf["dependents"], json!([]));
    let Some(EntityData::Part { features, .. }) =
        host.document().get(part.body).map(|e| e.data.clone())
    else {
        panic!("the part stays");
    };
    assert_eq!(features, [part.producer, fillet]);

    // Unknown and root IDs are refused without changing the document.
    let before = host.document();
    assert!(
        run(
            "delete_entity",
            json!({ "expected_revision": 3, "id": EntityId::new().to_string(), "cascade": true }),
        )
        .is_err()
    );
    assert!(
        run(
            "delete_entity",
            json!({ "expected_revision": 3, "id": before.id().to_string(), "cascade": true }),
        )
        .is_err()
    );
    assert_eq!(host.document(), before);
}
