use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    document::{Definition, Document, Entity, EntityData, Mutation},
    model::{self, PrimitiveParameters, PrimitiveRequest, length},
};
use multicad_model::{Compiler, EvaluationCache, GeometryBackend};

#[test]
fn full_compile_validates_current_dependency_closure_instead_of_retained_history() {
    let mut session = Session::new(Document::new("Retained source history").unwrap()).unwrap();
    let actor = Actor::local(Channel::Mcp);
    let request = PrimitiveRequest::new(
        "Current box",
        PrimitiveParameters::Box {
            size: [2.; 3].map(length),
            frame: Default::default(),
        },
    );
    let draft =
        model::prepare_primitive(&session, session.revision(), actor.clone(), &request).unwrap();
    session.commit(draft.prepared).unwrap();
    let compiler = Compiler::default();
    let previous = compiler
        .compile(session.document(), session.revision())
        .unwrap();
    let legacy = Entity::new(
        "Preserved unavailable producer",
        EntityData::Feature {
            definition: Definition {
                type_id: "legacy.unavailable".into(),
                schema_version: 1,
                data: serde_json::json!({}),
            },
            inputs: vec![],
            suppressed: false,
        },
    );
    let mut body = session.document().get(request.body).unwrap().clone();
    let EntityData::Part { features, .. } = &mut body.data else {
        panic!()
    };
    features.insert(0, legacy.id);
    session
        .apply(Request {
            expected_revision: session.revision(),
            actor: actor.clone(),
            label: "Preserve history".into(),
            edits: vec![
                Mutation::Insert(legacy.clone()),
                Mutation::Replace(body.clone()),
            ],
        })
        .unwrap();
    let document = session.document().clone();
    let scoped = compiler
        .compile_sources_with_backend(
            &document,
            session.revision(),
            &[request.feature],
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    assert_eq!(scoped.snapshot().shapes[&request.feature].volume_mm3, 8.);
    let full = compiler.compile(&document, session.revision()).unwrap();
    assert_eq!(full.snapshot().bodies[&request.body], request.feature);
    assert_eq!(full.snapshot().shapes[&request.feature].volume_mm3, 8.);
    assert!(!full.snapshot().shapes.contains_key(&legacy.id));
    assert_eq!(previous.snapshot().shapes[&request.feature].volume_mm3, 8.);
    assert_eq!(session.document(), &document);
    // The same unknown producer is required when the current feature names it.
    let mut active = document.get(request.feature).unwrap().clone();
    let required = model::TransformDefinition {
        input: model::BodyInput {
            body: request.body,
            producer: legacy.id,
        },
        output: model::BodyOutput { body: request.body },
        parameters: Default::default(),
    };
    let EntityData::Feature {
        definition, inputs, ..
    } = &mut active.data
    else {
        panic!()
    };
    *definition = required.definition().unwrap();
    *inputs = required.inputs().unwrap();
    session
        .apply(Request {
            expected_revision: session.revision(),
            actor: actor.clone(),
            label: "Require preserved producer".into(),
            edits: vec![Mutation::Replace(active)],
        })
        .unwrap();
    let invalid = session.document().clone();
    let mut cache = EvaluationCache::default();
    compiler
        .compile_cached_with_backend(
            &mut cache,
            &document,
            session.revision() - 1,
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let failure = compiler
        .compile_cached_with_backend(
            &mut cache,
            &invalid,
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .err()
        .unwrap();
    assert_eq!(failure.entity, Some(legacy.id));
    assert!(
        failure
            .message
            .contains("no provider for legacy.unavailable")
    );
    assert_eq!(session.document(), &invalid);
    let retry = compiler
        .compile_cached_with_backend(
            &mut cache,
            &document,
            session.revision() - 1,
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    assert_eq!(retry.snapshot().shapes[&request.feature].volume_mm3, 8.);
    // An unavailable current head is likewise rejected.
    let EntityData::Part { features, .. } = &mut body.data else {
        panic!()
    };
    features.swap(0, 1);
    session
        .apply(Request {
            expected_revision: session.revision(),
            actor,
            label: "Unavailable current head".into(),
            edits: vec![
                Mutation::Replace(body),
                Mutation::Replace(document.get(request.feature).unwrap().clone()),
            ],
        })
        .unwrap();
    let failure = compiler
        .compile(session.document(), session.revision())
        .err()
        .unwrap();
    assert_eq!(failure.entity, Some(legacy.id));
}
