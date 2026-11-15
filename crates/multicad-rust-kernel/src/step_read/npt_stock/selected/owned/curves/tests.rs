use super::*;
use crate::step_read::owned_insert::DocumentOptions;
fn source(quarter: bool, mapped: bool) -> OwnedMaterial {
    let text = super::super::super::super::tests::fixture(quarter);
    let text = if mapped {
        super::super::super::tests::mapped(&text, "(-3.,0.,0.)", "1.")
    } else {
        text
    };
    let doc = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    Selected::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap()
}
fn independent_point(source: &OwnedMaterial, edge: usize, parameter: f64) -> [I; 3] {
    let curve = &source.curves()[source.edges()[edge].curve_index()];
    let span = ((parameter * curve.spans() as f64) as usize).min(curve.spans() - 1);
    let t = I::point(parameter * curve.spans() as f64 - span as f64);
    let mut controls = std::array::from_fn::<_, 12, _>(|i| {
        source
            .point(curve.point_indices()[span * 11 + i])
            .unwrap()
            .coordinates()
            .map(|s| interval(s.bounds()))
    });
    for remaining in (1..12).rev() {
        for i in 0..remaining {
            controls[i] = std::array::from_fn(|k| {
                (I::point(1.) - t) * controls[i][k] + t * controls[i + 1][k]
            });
        }
    }
    source
        .map_box_mm(controls[0].map(bounds))
        .unwrap()
        .map(interval)
}
#[test]
fn actual_original_edges_keep_owners_knots_world_mapping_and_whole_curve_error() {
    for (quarter, mapped) in [(false, false), (true, false), (false, true)] {
        let source = source(quarter, mapped);
        let vertices = (0..source.vertices().len())
            .map(|v| source.sample_vertex(v, Default::default(), None).unwrap())
            .collect::<Vec<_>>();
        for (index, sample) in vertices.iter().enumerate() {
            assert_eq!(sample.vertex_index(), index);
            assert_eq!(sample.vertex_entity(), source.vertices()[index].entity());
            assert_eq!(
                sample.point_entity(),
                source
                    .point(source.vertices()[index].point_index())
                    .unwrap()
                    .entity()
            );
            assert!(std::sync::Arc::ptr_eq(&source.0, &sample.source().0));
            assert!(
                sample.work_used() > 0
                    && sample.accounted_peak_bytes() >= source.retained_budget_bytes()
            );
            for (p, b) in sample
                .position_mm()
                .into_iter()
                .zip(sample.position_bounds_mm())
            {
                assert!(p >= b[0] && p <= b[1]);
            }
        }
        for (edge, owner) in source.edges().iter().enumerate() {
            let sampled = source
                .sample_edge(edge, 0.05, POINT_LIMIT, Default::default(), None)
                .unwrap();
            assert!(std::sync::Arc::ptr_eq(&source.0, &sampled.source().0));
            assert_eq!(sampled.edge_index(), edge);
            assert_eq!(sampled.edge_entity(), owner.entity());
            assert_eq!(
                sampled.curve_entity(),
                source.curves()[owner.curve_index()].entity()
            );
            assert_eq!(
                sampled.vertex_entities(),
                owner
                    .vertex_indices()
                    .map(|v| source.vertices()[v].entity())
            );
            assert_eq!(
                sampled.points_mm().first(),
                Some(&vertices[owner.vertices[0]].position_mm())
            );
            assert_eq!(
                sampled.points_mm().last(),
                Some(&vertices[owner.vertices[1]].position_mm())
            );
            assert_eq!(sampled.parameters().len(), sampled.points_mm().len());
            assert!(sampled.parameters().windows(2).all(|p| p[0] < p[1]));
            let spans = source.curves()[owner.curve_index()].spans();
            for knot in 0..=spans {
                assert!(sampled.parameters().contains(&(knot as f64 / spans as f64)));
            }
            for index in
                (0..sampled.points_mm().len() - 1).step_by((sampled.points_mm().len() / 8).max(1))
            {
                let t = sampled.parameters()[index].midpoint(sampled.parameters()[index + 1]);
                let actual = independent_point(&source, edge, t);
                let line = std::array::from_fn(|k| {
                    sampled.points_mm()[index][k].midpoint(sampled.points_mm()[index + 1][k])
                });
                assert!(
                    l1(subtract(actual, line.map(I::point))) <= sampled.position_error_mm() + 1e-10
                );
            }
            let len = sampled.length_bounds_mm();
            assert!(len[0] > 0. && len[1] >= len[0] && len[1].is_finite());
            assert!(sampled.position_error_mm() <= 0.05);
            assert_eq!(sampled.work_used() + sampled.remaining_work(), 16_000_000);
            assert!(sampled.accounted_peak_bytes() >= source.retained_budget_bytes());
        }
    }
}
#[test]
fn actual_edge_queries_refuse_precision_and_resource_frontiers_without_source_changes() {
    let source = source(false, false);
    let edge = source
        .edges()
        .iter()
        .position(|e| source.curves()[e.curve].spans() > 1)
        .unwrap();
    assert_eq!(
        source
            .sample_vertex(source.vertices().len(), Default::default(), None)
            .err(),
        Some(Error::Invalid)
    );
    assert_eq!(
        source
            .sample_edge(
                source.edges().len(),
                0.1,
                POINT_LIMIT,
                Default::default(),
                None
            )
            .err(),
        Some(Error::Invalid)
    );
    for deflection in [f64::NAN, f64::INFINITY, 0., -1.] {
        assert_eq!(
            source
                .sample_edge(edge, deflection, POINT_LIMIT, Default::default(), None)
                .err(),
            Some(Error::Invalid)
        );
    }
    for limit in [0, 1, POINT_LIMIT + 1] {
        assert_eq!(
            source
                .sample_edge(edge, 0.1, limit, Default::default(), None)
                .err(),
            Some(Error::Invalid)
        );
    }
    for options in [
        Options {
            maximum_work: 0,
            ..Default::default()
        },
        Options {
            maximum_work: 16_000_001,
            ..Default::default()
        },
        Options {
            maximum_bytes: 0,
            ..Default::default()
        },
        Options {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
    ] {
        assert_eq!(
            source.sample_vertex(0, options, None).err(),
            Some(Error::Invalid)
        );
        assert_eq!(
            source
                .sample_edge(edge, 0.1, POINT_LIMIT, options, None)
                .err(),
            Some(Error::Invalid)
        );
    }
    for options in [
        Options {
            maximum_work: 1,
            ..Default::default()
        },
        Options {
            maximum_bytes: 1,
            ..Default::default()
        },
        Options {
            maximum_bytes: source.retained_budget_bytes() + SCRATCH + 1,
            ..Default::default()
        },
    ] {
        assert_eq!(
            source
                .sample_edge(edge, 0.1, POINT_LIMIT, options, None)
                .err(),
            Some(Error::Limit)
        );
    }
    assert_eq!(
        source
            .sample_vertex(
                0,
                Options {
                    maximum_work: 1,
                    ..Default::default()
                },
                None
            )
            .err(),
        Some(Error::Limit)
    );
    assert_eq!(
        source
            .sample_vertex(0, Default::default(), Some(&|| true))
            .err(),
        Some(Error::Cancelled)
    );
    assert_eq!(
        source
            .sample_edge(edge, 0.1, POINT_LIMIT, Default::default(), Some(&|| true))
            .err(),
        Some(Error::Cancelled)
    );
    let calls = std::cell::Cell::new(0usize);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() > 20
    };
    assert_eq!(
        source
            .sample_edge(edge, 0.1, POINT_LIMIT, Default::default(), Some(&cancel))
            .err(),
        Some(Error::Cancelled)
    );
    assert_eq!(
        source
            .sample_edge(edge, 0.05, 2, Default::default(), None)
            .err(),
        Some(Error::Limit)
    );
    let limit = source.curves()[source.edges()[edge].curve].spans() + 1;
    assert_eq!(
        source
            .sample_edge(edge, 0.05, limit, Default::default(), None)
            .err(),
        Some(Error::Limit)
    );
    assert_eq!(
        source
            .sample_edge(edge, 1e-30, POINT_LIMIT, Default::default(), None)
            .err(),
        Some(Error::Unsupported)
    );
    let coarse = source
        .sample_edge(edge, 0.1, POINT_LIMIT, Default::default(), None)
        .unwrap();
    let fine = source
        .sample_edge(edge, 0.01, POINT_LIMIT, Default::default(), None)
        .unwrap();
    assert!(fine.points_mm().len() >= coarse.points_mm().len());
    assert!(
        fine.length_bounds_mm()[1] - fine.length_bounds_mm()[0]
            <= coarse.length_bounds_mm()[1] - coarse.length_bounds_mm()[0] + 1e-8
    );
    assert!(std::sync::Arc::ptr_eq(&source.0, &fine.source().0));
}
