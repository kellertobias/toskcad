//! Timing probe for a fully cached full compile of the shipped Toolcase
//! (CAD-396). Each compile uses a new revision, so the full-model replay
//! misses and every feature is reused from the evaluation cache.
//!
//! cargo test -p multicad-model --release --test cached_compile_probe -- --ignored --nocapture
//! Optional: CAD396_PROJECT=path.tcd, CAD396_COMPILES=n.
use multicad_core::project_file;
use multicad_model::{Compiler, EvaluationCache, GeometryBackend};

#[test]
#[ignore = "timing probe; run explicitly in a release build"]
fn toolcase_fully_cached_full_compile() {
    let path = std::env::var("CAD396_PROJECT").unwrap_or_else(|_| {
        concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/toolcase.tcd").into()
    });
    let compiles: u64 = std::env::var("CAD396_COMPILES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(16);
    let document = project_file::read(
        &mut std::fs::read(&path).unwrap().as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    let compiler = Compiler::default();
    let mut cache = EvaluationCache::default();
    let mut times = vec![];
    let mut previous = None;
    for revision in 0..=compiles {
        let started = std::time::Instant::now();
        let model = compiler
            .compile_cached_with_backend(
                &mut cache,
                &document,
                revision,
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        let compile = started.elapsed().as_secs_f64() * 1000.;
        let stats = cache.last_compile().clone();
        let dropped = std::time::Instant::now();
        drop(previous.replace(model));
        let drop_ms = dropped.elapsed().as_secs_f64() * 1000.;
        let t = previous.as_ref().unwrap().timings();
        eprintln!(
            "revision {revision}: compile {compile:.1} ms, drop previous {drop_ms:.1} ms, evaluated {}, reused {}, compositions evaluated {} reused {}, full {}, assembly+occurrence {:.1}, publication {:.1}, preparation {:.1}",
            stats.evaluated.len(),
            stats.reused.len(),
            stats.composition_evaluated.len(),
            stats.composition_reused.len(),
            stats.reused_full_model,
            t.assembly_and_occurrence_ms,
            t.shape_publication_ms,
            t.source_preparation_ms,
        );
        eprintln!(
            "  compositions retained {} ({:.1} MB estimated), feature cache {:.1} MB",
            cache.composition_entries(),
            cache.composition_bytes() as f64 / 1e6,
            cache.estimated_bytes() as f64 / 1e6
        );
        eprintln!("  timings {}", serde_json::to_string(t).unwrap());
        if revision > 0 {
            assert!(stats.evaluated.is_empty(), "{:?}", stats.evaluated);
            times.push(compile);
        }
    }
    times.sort_by(f64::total_cmp);
    eprintln!(
        "cached compiles {}: p50 {:.1} ms, max {:.1} ms",
        times.len(),
        times[times.len() / 2],
        times[times.len() - 1]
    );
}
