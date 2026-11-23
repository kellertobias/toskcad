//! CAD-396: Assembly and Instance compositions reuse an earlier result only
//! for the same member solids and placements, and the reused model equals a
//! cache-free compile of the same document and revision.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{EntityData, EntityId},
    model::{self, InstanceMoveRequest, NumericTransform, length},
    project_file,
};
use multicad_model::{Compiler, EvaluationCache, GeometryBackend, Model};
use std::{cell::Cell, collections::BTreeSet};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn toolcase() -> Session {
    let document = project_file::read(
        &mut include_bytes!("../../../examples/toolcase.tcd").as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    Session::new(document).unwrap()
}
fn cached(cache: &mut EvaluationCache, session: &Session, revision: u64) -> Model {
    Compiler::default()
        .compile_cached_with_backend(
            cache,
            session.document(),
            revision,
            GeometryBackend::Rust,
            None,
        )
        .unwrap()
}
fn fresh(session: &Session, revision: u64) -> Model {
    cached(&mut EvaluationCache::new(0, 0), session, revision)
}
fn compositions(document: &multicad_core::document::Document) -> BTreeSet<EntityId> {
    document
        .entities()
        .filter(|e| {
            matches!(&e.data, EntityData::Instance { .. })
                || matches!(&e.data, EntityData::Assembly { instances } if !instances.is_empty())
        })
        .map(|e| e.id)
        .collect()
}
fn compare(a: &Model, b: &Model, revision: u64) {
    let (a, b) = (a.snapshot(), b.snapshot());
    assert_eq!(a.stamp, b.stamp);
    assert_eq!(a.stamp.revision, revision);
    assert_eq!(a.references, b.references);
    assert_eq!(a.bodies, b.bodies);
    assert_eq!(a.body_sources, b.body_sources);
    assert_eq!(
        a.shapes.keys().collect::<Vec<_>>(),
        b.shapes.keys().collect::<Vec<_>>()
    );
    for (id, x) in &a.shapes {
        let y = &b.shapes[id];
        assert_eq!(x.topology, y.topology, "{id}");
        assert_eq!(x.topology.revision, revision);
        assert_eq!(x.solid_count, y.solid_count);
        assert_eq!(x.volume_mm3.to_bits(), y.volume_mm3.to_bits(), "{id}");
        assert_eq!(x.bounds, y.bounds);
        assert_eq!(x.mesh.vertices, y.mesh.vertices, "{id}");
        assert_eq!(x.mesh.indices, y.mesh.indices, "{id}");
        assert_eq!(x.mesh.face_indices, y.mesh.face_indices, "{id}");
        assert_eq!(x.edges.len(), y.edges.len(), "{id}");
        for (e, f) in x.edges.iter().zip(y.edges.iter()) {
            assert_eq!(e.points_mm, f.points_mm, "{id}");
        }
    }
    assert_eq!(
        a.occurrences.keys().collect::<Vec<_>>(),
        b.occurrences.keys().collect::<Vec<_>>()
    );
    for (id, x) in &a.occurrences {
        let y = &b.occurrences[id];
        assert_eq!(
            (x.instance, x.component, x.producer, x.pose, x.bounds),
            (y.instance, y.component, y.producer, y.pose, y.bounds)
        );
        assert_eq!(x.topology, y.topology, "{id}");
        assert_eq!(x.topology.revision, revision);
    }
    assert_eq!(
        a.leaf_occurrences.keys().collect::<Vec<_>>(),
        b.leaf_occurrences.keys().collect::<Vec<_>>()
    );
    for (path, x) in &a.leaf_occurrences {
        let y = &b.leaf_occurrences[path];
        assert_eq!(
            (&x.path, x.component, x.producer, x.pose, x.bounds),
            (&y.path, y.component, y.producer, y.pose, y.bounds)
        );
        assert_eq!(x.topology, y.topology);
        assert_eq!(x.topology.revision, revision);
    }
}

#[test]
fn unchanged_toolcase_compositions_are_reused_and_equal_a_fresh_compile() {
    let mut session = toolcase();
    let nodes = compositions(session.document());
    assert!(!nodes.is_empty());
    let mut cache = EvaluationCache::default();
    cached(&mut cache, &session, 1);
    let stats = cache.last_compile();
    assert!(stats.composition_reused.is_empty());
    assert!(!stats.composition_evaluated.is_empty());
    assert!(stats.composition_evaluated.is_subset(&nodes));
    assert!(cache.composition_entries() > 0);

    // A new revision misses the full-model replay; every feature and every
    // composition is reused and the model equals a cache-free compile.
    let reused = cached(&mut cache, &session, 2);
    let stats = cache.last_compile().clone();
    assert!(stats.evaluated.is_empty(), "{:?}", stats.evaluated);
    assert!(!stats.reused_full_model);
    assert!(stats.composition_evaluated.is_empty());
    assert!(!stats.composition_reused.is_empty());
    compare(&reused, &fresh(&session, 2), 2);

    // Moving one placed occurrence re-places it and recomposes only what
    // contains it; everything else is still reused.
    // The first reused occurrence that accepts a free placement edit.
    let (moved, prepared) = session
        .document()
        .entities()
        .filter(|e| {
            matches!(&e.data, EntityData::Instance { .. })
                && stats.composition_reused.contains(&e.id)
        })
        .find_map(|e| {
            let movement = InstanceMoveRequest::new(
                e.id,
                NumericTransform {
                    translation: [length(7.), length(-3.), length(2.)],
                    ..Default::default()
                },
            );
            model::prepare_instance_moves(&session, session.revision(), actor(), &[movement])
                .ok()
                .map(|prepared| (e.id, prepared.prepared))
        })
        .unwrap();
    session.commit(prepared).unwrap();
    let revision = session.revision();
    let after = cached(&mut cache, &session, revision);
    let stats = cache.last_compile().clone();
    assert!(stats.composition_evaluated.contains(&moved), "{moved}");
    assert!(!stats.composition_reused.contains(&moved));
    compare(&after, &fresh(&session, revision), revision);
}

#[test]
fn superseded_and_zero_capacity_compiles_publish_no_compositions() {
    let session = toolcase();
    let mut cache = EvaluationCache::default();
    let count = Cell::new(0usize);
    Compiler::default()
        .compile_cached_with_backend(
            &mut EvaluationCache::default(),
            session.document(),
            1,
            GeometryBackend::Rust,
            Some(&|| {
                count.set(count.get() + 1);
                false
            }),
        )
        .unwrap();
    // Cancel at the last checkpoint, after every composition was staged.
    let observed = Cell::new(0usize);
    assert!(
        Compiler::default()
            .compile_cached_with_backend(
                &mut cache,
                session.document(),
                1,
                GeometryBackend::Rust,
                Some(&|| {
                    observed.set(observed.get() + 1);
                    observed.get() >= count.get()
                }),
            )
            .is_err()
    );
    assert_eq!(cache.composition_entries(), 0);
    cached(&mut cache, &session, 2);
    assert!(cache.last_compile().composition_reused.is_empty());

    let mut disabled = EvaluationCache::new(0, 0);
    for revision in 1..=2 {
        cached(&mut disabled, &session, revision);
        assert_eq!(disabled.composition_entries(), 0);
        assert!(disabled.last_compile().composition_reused.is_empty());
    }
}
