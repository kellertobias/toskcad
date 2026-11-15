//! Full native connection validation must stay live around cached exact None
//! predicates. These controls do not assert an operator-latency target.
#[path = "../../multicad-mcp/tests/support/screw_connection.rs"]
mod fixture;
use multicad_core::{
    commands::{Actor, Channel, Session},
    model::*,
};
use multicad_model::{Compiler, EvaluationCache, GeometryBackend};
use std::time::Instant;
fn actor() -> Actor {
    Actor::local(Channel::CommandLine)
}
fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}

#[test]
fn warm_metric_connections_reuse_none_but_recheck_engagement_and_new_foreign_geometry() {
    check_warm_connections(GeometryBackend::OpenCascade);
}
#[test]
fn rust_warm_connections_reuse_sources_but_recheck_engagement_and_foreign_geometry() {
    check_warm_connections(GeometryBackend::Rust);
}
fn check_warm_connections(backend: GeometryBackend) {
    for (document, request, _) in [fixture::fixture(), fixture::tapped_fixture()] {
        let mut session = Session::new(document).unwrap();
        session
            .commit(
                prepare_screw_connection(&session, 0, actor(), &request)
                    .unwrap()
                    .prepared,
            )
            .unwrap();
        let compiler = compiler();
        let compile = |cache: &mut EvaluationCache,
                       document: &multicad_core::document::Document,
                       revision| {
            compiler.compile_cached_with_backend(cache, document, revision, backend, None)
        };
        let mut cache = EvaluationCache::default();
        let cold = Instant::now();
        let first = compile(&mut cache, session.document(), 0).unwrap();
        let cold = cold.elapsed();
        let initially_evaluated = cache.last_compile().evaluated.clone();
        let initially_uncacheable = cache.last_compile().uncacheable.clone();
        let first_stats = cache.native_intersection_proof_stats();
        if backend == GeometryBackend::OpenCascade {
            assert!(
                first_stats.published > 0 && first_stats.native_queries > 0,
                "The OCCT connection must publish actual native None queries: {first_stats:?}"
            );
        } else {
            assert!(!cache.last_compile().evaluated.is_empty());
            assert!(!first.snapshot().screw_connection_checks.is_empty());
            assert!(
                first
                    .snapshot()
                    .screw_connection_checks
                    .iter()
                    .all(|c| c.native_component_collision_checked)
            );
            let unchanged = compile(&mut cache, session.document(), 0).unwrap();
            assert!(cache.last_compile().reused_full_model);
            assert!(cache.last_compile().evaluated.is_empty());
            assert_eq!(
                unchanged.snapshot().screw_connection_checks,
                first.snapshot().screw_connection_checks
            );
        }
        let warm = Instant::now();
        let second = compile(&mut cache, session.document(), 1).unwrap();
        let warm = warm.elapsed();
        let warm_stats = cache.native_intersection_proof_stats();
        if backend == GeometryBackend::OpenCascade {
            assert!(
                warm_stats.hits > first_stats.hits,
                "Unchanged exact geometry should reuse completed None evidence: {warm_stats:?}"
            );
        } else {
            let warm = cache.last_compile();
            assert!(
                warm.reused_full_model
                    || initially_evaluated
                        .difference(&initially_uncacheable)
                        .all(|feature| warm.reused.contains(feature)),
                "Unchanged cacheable Rust geometry must be reused: {warm:?}"
            );
            assert!(
                warm.evaluated.is_subset(&initially_uncacheable),
                "Only providers requiring live checks may run again: {warm:?}"
            );
        }
        assert_eq!(
            second.snapshot().screw_connection_checks,
            first.snapshot().screw_connection_checks
        );
        assert_eq!(second.snapshot().bodies, first.snapshot().bodies);
        assert_eq!(
            second
                .snapshot()
                .leaf_occurrences
                .keys()
                .collect::<Vec<_>>(),
            first.snapshot().leaf_occurrences.keys().collect::<Vec<_>>()
        );
        eprintln!(
            "{:?} {backend:?} cache cold={cold:?} warm={warm:?}, {warm_stats:?}",
            request.mode
        );
        let mut invalid =
            read_screw_connection_request(session.document(), request.connection).unwrap();
        invalid.parameters.minimum_engagement = length(4.);
        let invalid =
            prepare_screw_connection(&session, session.revision(), actor(), &invalid).unwrap();
        // Reuse the last successful revision deliberately: the complete source
        // hash must invalidate readiness even if a caller repeats the stamp.
        let error = compile(&mut cache, invalid.prepared.document(), 1)
            .err()
            .expect("A cached None must not skip current engagement validation");
        assert!(
            error.message.to_lowercase().contains("engagement"),
            "{error}"
        );
        assert_eq!(
            cache.native_intersection_proof_stats().published,
            warm_stats.published
        );
        if request.mode != ScrewConnectionMode::ThroughNut {
            continue;
        }
        let collider = PrimitiveRequest::new(
            "Independent foreign component",
            PrimitiveParameters::Box {
                size: std::array::from_fn(|_| length(0.5)),
                frame: PrimitiveFrame {
                    origin: std::array::from_fn(|_| length(100.)),
                    ..Default::default()
                },
            },
        );
        session
            .commit(
                prepare_primitive(&session, session.revision(), actor(), &collider)
                    .unwrap()
                    .prepared,
            )
            .unwrap();
        compile(&mut cache, session.document(), 3).unwrap();
        let before = cache.native_intersection_proof_stats();
        let before_move = session.document().clone();
        let moved = PrimitiveParameters::Box {
            size: std::array::from_fn(|_| length(0.5)),
            frame: PrimitiveFrame {
                origin: [length(9.75), length(9.75), length(4.5)],
                ..Default::default()
            },
        };
        session
            .commit(
                prepare_primitive_edit(
                    &session,
                    session.revision(),
                    actor(),
                    collider.feature,
                    &moved,
                )
                .unwrap()
                .prepared,
            )
            .unwrap();
        let source = session.document().clone();
        let error = compile(&mut cache, session.document(), 3).err().expect(
            "A newly intersecting foreign component must invalidate full connection readiness",
        );
        let after = cache.native_intersection_proof_stats();
        if backend == GeometryBackend::OpenCascade {
            assert!(
                error.message.contains("unintended component interference"),
                "{error}"
            );
            assert!(
                after.native_queries > before.native_queries,
                "The changed actual foreign solid must reach a new native query"
            );
        } else {
            assert!(
                error
                    .message
                    .contains("clearance to another nearby physical component"),
                "{error}"
            );
        }
        assert_eq!(
            after.published, before.published,
            "A later failed query cannot publish earlier staged results"
        );
        assert_eq!(session.document(), &source);
        compile(&mut cache, &before_move, 5)
            .expect("A failed candidate must not poison prior valid connection proofs");
    }
}
