use super::*;
const WORK: usize = 512_000_000;
const BYTES: usize = 64 * 1024 * 1024;
fn quality() -> Tessellation {
    Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    }
}
fn plan(quarter: bool) -> StockClipCellPlan {
    super::super::super::tests::canonical(quarter, 1.)
        .stock_clip_plan(Default::default(), None)
        .unwrap()
        .cell_plan(128, 128, WORK, BYTES, None)
        .unwrap()
}
fn check(quarter: bool) {
    let source = plan(quarter);
    let plain = source.mesh_candidate(quality(), WORK, BYTES, None).unwrap();
    let traced = source
        .mesh_candidate_with_cap_trace(quality(), WORK, BYTES, None)
        .unwrap();
    assert!(std::ptr::eq(traced.source(), &source));
    assert_eq!(traced.candidate().mesh(), plain.mesh());
    assert_eq!(traced.cap_edges().len(), source.cap_edge_count());
    assert_eq!(
        traced.work_used(),
        plain.work_used() + 513 * source.cap_edge_count()
    );
    assert!(traced.accounted_peak_bytes() > plain.retained_budget_bytes());
    assert!(traced.accounted_peak_bytes() <= BYTES);
    for station in [StockClipStation::Start, StockClipStation::End] {
        let edges = traced
            .cap_edges()
            .iter()
            .filter(|e| e.station == station)
            .collect::<Vec<_>>();
        let mut outgoing = std::collections::BTreeMap::new();
        for edge in &edges {
            assert!(
                outgoing
                    .insert(edge.vertices[0], edge.vertices[1])
                    .is_none()
            );
            assert!(edge.chord_error_mm <= quality().deflection);
            for j in 0..2 {
                assert_eq!(
                    Some(edge.owners[j]),
                    source.vertex_owner(edge.vertices[j] as usize)
                );
                let p = &plain.mesh().vertices[3 * edge.vertices[j] as usize..][..3];
                let distance = (0..3)
                    .map(|k| (I::point(p[k]) - b(edge.source_positions_mm[j][k])).square())
                    .fold(I::point(0.), |a, b| a + b)
                    .sqrt()
                    .hi;
                assert!(distance <= edge.chord_error_mm);
                let z = if station == StockClipStation::Start {
                    0.
                } else {
                    source.source().source().spec().length_mm
                };
                assert_eq!(p[2], z);
                assert!(edge.parameters[j][1][0] <= z && z <= edge.parameters[j][1][1]);
            }
        }
        let start = edges[0].vertices[0];
        let mut current = start;
        let mut visited = std::collections::BTreeSet::new();
        loop {
            assert!(visited.insert(current));
            current = outgoing[&current];
            if current == start {
                break;
            }
        }
        assert_eq!(visited.len(), edges.len());
    }
}
#[test]
fn actual_eighth_cap_trace_keeps_exact_owners_and_closed_literal_rings() {
    check(false);
}
#[test]
fn actual_quarter_cap_trace_keeps_exact_owners_and_closed_literal_rings() {
    check(true);
}
#[test]
fn actual_cap_trace_has_atomic_work_storage_cancel_and_input_frontiers() {
    let source = plan(false);
    let plain = source.mesh_candidate(quality(), WORK, BYTES, None).unwrap();
    let callbacks = std::cell::Cell::new(0usize);
    let count = || {
        callbacks.set(callbacks.get() + 1);
        false
    };
    source
        .mesh_candidate(quality(), WORK, BYTES, Some(&count))
        .unwrap();
    let candidate_callbacks = callbacks.get();
    callbacks.set(0);
    let cancel_in_trace = || {
        callbacks.set(callbacks.get() + 1);
        callbacks.get() > candidate_callbacks + 1
    };
    assert!(matches!(
        source.mesh_candidate_with_cap_trace(quality(), WORK, BYTES, Some(&cancel_in_trace)),
        Err(Error::Cancelled)
    ));
    assert_eq!(callbacks.get(), candidate_callbacks + 2);
    assert!(matches!(
        source.mesh_candidate_with_cap_trace(quality(), plain.work_used(), BYTES, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        source.mesh_candidate_with_cap_trace(quality(), WORK, 1, None),
        Err(Error::Limit)
    ));
    let exact_metadata = size_of::<StockClipTracedMesh<'_>>()
        + source.cap_edge_count() * size_of::<StockClipCapTraceEdge>();
    assert!(matches!(
        source.mesh_candidate_with_cap_trace(quality(), WORK, exact_metadata, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        source.mesh_candidate_with_cap_trace(quality(), WORK, BYTES, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    for (work, bytes) in [(0, BYTES), (WORK + 1, BYTES), (WORK, 0), (WORK, BYTES + 1)] {
        assert!(matches!(
            source.mesh_candidate_with_cap_trace(quality(), work, bytes, None),
            Err(Error::Invalid)
        ));
    }
    let invalid = Tessellation {
        deflection: f64::NAN,
        ..quality()
    };
    assert!(matches!(
        source.mesh_candidate_with_cap_trace(invalid, WORK, BYTES, None),
        Err(Error::Invalid)
    ));
    assert_eq!(
        source
            .mesh_candidate_with_cap_trace(quality(), WORK, BYTES, None)
            .unwrap()
            .candidate()
            .mesh(),
        plain.mesh()
    );
}
