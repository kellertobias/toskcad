//! CAD-426: Delete removes an entry together with its dependents in one
//! canonical transaction, names them first, and never cascades silently.
use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    document::{Document, Entity, EntityData, EntityId, Link, Mutation, Placement},
    hierarchy::{self, Action},
    model::*,
    project_file::{self, ProjectFile},
    references::{GeometryIntent, Origin, TopologyKind},
};
use std::collections::BTreeSet;

fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}
fn edge(source: EntityId) -> GeometryIntent {
    GeometryIntent {
        id: EntityId::new(),
        kind: TopologyKind::Edge,
        origins: BTreeSet::from([Origin {
            source,
            role: "selected-source-edge".into(),
        }]),
        boundary: None,
        attachment: None,
    }
}
fn named(session: &mut Session, id: EntityId, name: &str) {
    let mut entity = session.document().get(id).unwrap().clone();
    entity.name = name.into();
    apply(session, vec![Mutation::Replace(entity)]);
}
fn apply(session: &mut Session, edits: Vec<Mutation>) {
    session
        .apply(Request {
            expected_revision: session.revision(),
            actor: actor(),
            label: "Fixture".into(),
            edits,
        })
        .unwrap();
}
fn boxed(session: &mut Session, name: &str) -> BodyInput {
    let request = PrimitiveRequest::new(
        name,
        PrimitiveParameters::Box {
            size: [length(30.), length(20.), length(10.)],
            frame: Default::default(),
        },
    );
    let prepared = prepare_primitive(session, session.revision(), actor(), &request).unwrap();
    session.commit(prepared.prepared).unwrap();
    named(session, request.body, name);
    BodyInput {
        body: request.body,
        producer: request.feature,
    }
}
fn treat(session: &mut Session, input: &BodyInput, name: &str) -> BodyInput {
    let request = EdgeTreatmentRequest::new(
        input.clone(),
        EdgeTreatmentParameters {
            profile: EdgeTreatmentProfile::Fillet,
            size: length(1.),
        },
        vec![edge(input.producer)],
    );
    let prepared = prepare_edge_treatment(session, session.revision(), actor(), &request).unwrap();
    session.commit(prepared.prepared).unwrap();
    named(session, request.feature, name);
    BodyInput {
        body: input.body,
        producer: request.feature,
    }
}
fn delete(session: &mut Session, id: EntityId, cascade: bool) -> Result<(), String> {
    let action = if cascade {
        Action::DeleteWithDependents { id }
    } else {
        Action::Delete { id }
    };
    let prepared = hierarchy::prepare(session, session.revision(), actor(), action)
        .map_err(|error| error.to_string())?;
    session
        .commit(prepared.prepared)
        .map(|_| ())
        .map_err(|error| error.to_string())
}
fn names(session: &Session, id: EntityId) -> Vec<String> {
    hierarchy::deletion_plan(session, id)
        .unwrap()
        .dependents
        .into_iter()
        .map(|d| d.name)
        .collect()
}
fn features(document: &Document, part: EntityId) -> Vec<EntityId> {
    match document.get(part).map(|e| &e.data) {
        Some(EntityData::Part { features, .. }) => features.clone(),
        other => panic!("not a part: {other:?}"),
    }
}
fn reopened(document: &Document) -> Document {
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(document.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document
}

#[test]
fn deleting_a_part_removes_its_whole_history_in_one_undo_step_and_keeps_unrelated_parts() {
    let mut session = Session::new(Document::new("Part cascade").unwrap()).unwrap();
    let other = boxed(&mut session, "Spacer");
    let other_fillet = treat(&mut session, &other, "Fillet 2");
    let part = boxed(&mut session, "Extrude 1");
    let fillet = treat(&mut session, &part, "Fillet 1");
    let chamfer = treat(&mut session, &fillet, "Chamfer 1");
    let before = session.document().clone();
    let revision = session.revision();

    let plan = hierarchy::deletion_plan(&session, part.body).unwrap();
    assert_eq!(
        plan.question(),
        "Delete Extrude 1 and 2 dependent features (Fillet 1, Chamfer 1)?"
    );
    // Without the cascade option nothing is removed and the refusal names them.
    let refusal = delete(&mut session, part.body, false).unwrap_err();
    assert!(
        refusal.contains("“Fillet 1”, “Chamfer 1” depend on this entry"),
        "{refusal}"
    );
    let reason = hierarchy::deletion_available(&session, part.body).unwrap_err();
    assert!(reason.contains("Fillet 1") && reason.contains("Chamfer 1"));
    assert_eq!(session.document(), &before);
    assert_eq!(session.revision(), revision);

    delete(&mut session, part.body, true).unwrap();
    assert_eq!(session.revision(), revision + 1, "one transaction");
    assert_eq!(
        session.activity().last().unwrap().label,
        "Delete model entry and dependents"
    );
    for id in [part.body, part.producer, fillet.producer, chamfer.producer] {
        assert!(session.document().get(id).is_none());
    }
    for id in [other.body, other.producer, other_fillet.producer] {
        assert_eq!(session.document().get(id), before.get(id), "unrelated");
    }
    assert!(session.document().reference_diagnostics().is_empty());
    let after = session.document().clone();
    assert_eq!(&reopened(&after), &after);

    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &before, "Undo restores exactly");
    session.redo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &after);
}

#[test]
fn deleting_a_feature_removes_nested_dependents_and_keeps_the_earlier_history() {
    let mut session = Session::new(Document::new("Feature cascade").unwrap()).unwrap();
    let part = boxed(&mut session, "Block");
    let fillet = treat(&mut session, &part, "Fillet 1");
    let chamfer = treat(&mut session, &fillet, "Chamfer 1");
    let last = treat(&mut session, &chamfer, "Fillet 2");
    let before = session.document().clone();

    assert_eq!(names(&session, fillet.producer), ["Chamfer 1", "Fillet 2"]);
    assert_eq!(
        hierarchy::deletion_plan(&session, fillet.producer)
            .unwrap()
            .question(),
        "Delete Fillet 1 and 2 dependent features (Chamfer 1, Fillet 2)?"
    );
    assert!(delete(&mut session, fillet.producer, false).is_err());
    assert_eq!(session.document(), &before);

    // The last feature has no dependents and deletes without the option.
    hierarchy::deletion_available(&session, last.producer).unwrap();
    delete(&mut session, last.producer, false).unwrap();
    assert_eq!(
        features(session.document(), part.body),
        [part.producer, fillet.producer, chamfer.producer]
    );
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &before);

    delete(&mut session, fillet.producer, true).unwrap();
    assert_eq!(features(session.document(), part.body), [part.producer]);
    for id in [fillet.producer, chamfer.producer, last.producer] {
        assert!(session.document().get(id).is_none());
    }
    assert_eq!(
        session.document().get(part.producer),
        before.get(part.producer)
    );
    assert!(session.document().reference_diagnostics().is_empty());
    let after = session.document().clone();
    assert_eq!(&reopened(&after), &after);
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &before);
}

#[test]
fn deleting_the_creating_feature_removes_its_part_and_deleting_a_tool_removes_the_boolean() {
    let mut session = Session::new(Document::new("Consumers").unwrap()).unwrap();
    let target = boxed(&mut session, "Base");
    let tool = boxed(&mut session, "Boss");
    let boolean = BooleanRequest::new(Operation::Join, target.clone(), vec![tool.clone()]);
    let prepared = prepare_boolean(&session, session.revision(), actor(), &boolean).unwrap();
    session.commit(prepared.prepared).unwrap();
    named(&mut session, boolean.feature, "Join 1");
    let joined = BodyInput {
        body: target.body,
        producer: boolean.feature,
    };
    let fillet = treat(&mut session, &joined, "Fillet 1");
    let copy = Entity::new(
        "Linked base",
        EntityData::Instance {
            component: target.body,
            placement: Placement::default(),
            placement_source: None,
        },
    );
    // A drawing and a face-anchored plane only point at geometry: they stay
    // and report the lost reference.
    let plane = Entity::new(
        "Face plane",
        EntityData::Plane {
            origin_mm: [0.; 3],
            normal: [0., 0., 1.],
            anchor: Some(Link {
                target: tool.producer,
                role: "face".into(),
                geometry: Some(GeometryIntent {
                    id: EntityId::new(),
                    kind: TopologyKind::Face,
                    origins: BTreeSet::from([Origin {
                        source: tool.producer,
                        role: "top".into(),
                    }]),
                    boundary: None,
                    attachment: None,
                }),
            }),
        },
    );
    apply(
        &mut session,
        vec![
            Mutation::Insert(copy.clone()),
            Mutation::Insert(plane.clone()),
        ],
    );
    let before = session.document().clone();

    // Deleting a Boolean tool removes the Boolean and what consumed it.
    assert_eq!(names(&session, tool.body), ["Join 1", "Fillet 1"]);
    delete(&mut session, tool.body, true).unwrap();
    assert_eq!(features(session.document(), target.body), [target.producer]);
    assert_eq!(session.document().get(copy.id), Some(&copy));
    assert_eq!(session.document().get(plane.id), Some(&plane));
    let lost: Vec<_> = session
        .document()
        .reference_diagnostics()
        .into_iter()
        .map(|d| (d.owner, d.target))
        .collect();
    assert_eq!(lost, [(plane.id, tool.producer)]);
    assert_eq!(&reopened(session.document()), session.document());
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &before);

    // The creating feature takes its part, the part's later features and the
    // instances of that part with it.
    let plan = hierarchy::deletion_plan(&session, target.producer).unwrap();
    let listed: BTreeSet<_> = plan.dependents.iter().map(|d| d.id).collect();
    assert_eq!(
        listed,
        BTreeSet::from([target.body, copy.id, boolean.feature, fillet.producer])
    );
    delete(&mut session, target.producer, true).unwrap();
    for id in [
        target.body,
        target.producer,
        boolean.feature,
        fillet.producer,
        copy.id,
    ] {
        assert!(session.document().get(id).is_none());
    }
    assert_eq!(session.document().get(tool.body), before.get(tool.body));
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &before);
}

#[test]
fn deleting_a_container_removes_its_contents_and_stale_or_root_deletes_are_refused() {
    let mut session = Session::new(Document::new("Containers").unwrap()).unwrap();
    let outside = boxed(&mut session, "Outside");
    let inside = boxed(&mut session, "Inside");
    let fillet = treat(&mut session, &inside, "Fillet 1");
    let group = Entity::new("Panels", EntityData::Group {});
    apply(&mut session, vec![Mutation::Insert(group.clone())]);
    let moved = hierarchy::prepare(
        &session,
        session.revision(),
        actor(),
        Action::Move {
            entries: vec![inside.body],
            parent: Some(group.id),
            before: None,
        },
    )
    .unwrap();
    session.commit(moved.prepared).unwrap();
    let before = session.document().clone();
    assert!(names(&session, group.id).contains(&"Inside".to_string()));
    assert!(delete(&mut session, group.id, false).is_err());
    delete(&mut session, group.id, true).unwrap();
    for id in [group.id, inside.body, inside.producer, fillet.producer] {
        assert!(session.document().get(id).is_none());
    }
    assert_eq!(
        session.document().get(outside.body),
        before.get(outside.body)
    );
    session.undo(session.revision(), actor()).unwrap();
    assert_eq!(session.document(), &before);

    let root = session.document().id();
    assert!(hierarchy::deletion_plan(&session, root).is_err());
    let stale = hierarchy::prepare(
        &session,
        session.revision() + 1,
        actor(),
        Action::DeleteWithDependents { id: inside.body },
    );
    assert!(stale.is_err());
    assert_eq!(session.document(), &before);
}
