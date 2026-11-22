use super::*;
fn solid() -> Solid {
    Solid::rust_phillips_source_for_test(
        super::super::tests::fixture(0, 0.2, true),
        EntityId::from_uuid_bits(71),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        &|| false,
    )
    .unwrap()
}
#[test]
fn actual_source_curves_keep_topology_indices_vertices_and_literal_pose() {
    let solid = solid();
    let (body, frame) = solid.helical_body().unwrap();
    let (source, _) = solid.rust_phillips_screw_geometry().unwrap();
    let mut remaining = MAX_DISPLAY_EDGE_POINTS;
    let edges = super::edges(
        body,
        source,
        frame,
        0.05,
        GeometryPolicy::default(),
        &mut remaining,
    )
    .unwrap();
    assert_eq!(edges.len(), 76);
    assert!(remaining < MAX_DISPLAY_EDGE_POINTS);
    let topology = body.topology(frame, body.feature, 0).unwrap();
    let cutter = source.cavity().source().boundary_topology().unwrap();
    for (edge, source_edge) in edges.iter().zip(cutter.edges()) {
        let origin =
            super::super::topology::curve_origin(body, &cutter, source_edge.source_identity());
        assert!(
            topology
                .elements
                .iter()
                .any(|e| e.key.kind == TopologyKind::Edge
                    && e.key.ordinal == edge.index
                    && e.origins.contains(&origin))
        );
        assert!(edge.points_mm.len() >= 2);
        for (point, t) in [
            (edge.points_mm[0], 0.),
            (*edge.points_mm.last().unwrap(), 1.),
        ] {
            let expected = super::super::topology::placed_point(
                frame,
                cutter
                    .curve_point_mm_bounds(
                        cutter.edge_index(source_edge.source_identity()).unwrap(),
                        t,
                    )
                    .unwrap(),
                source.head_top_z(),
                GeometryPolicy::default(),
            )
            .unwrap();
            assert_eq!(point, expected);
        }
    }
    let shifted = PlaneFrame {
        origin: [32., -64., 128.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    };
    let mut budget = MAX_DISPLAY_EDGE_POINTS;
    let posed = super::edges(
        body,
        source,
        shifted,
        0.05,
        GeometryPolicy::default(),
        &mut budget,
    )
    .unwrap();
    assert_eq!(posed.len(), edges.len());
    for (a, b) in edges.iter().zip(posed) {
        assert_eq!(a.index, b.index);
        assert_eq!(a.points_mm.len(), b.points_mm.len());
    }
}
#[test]
fn source_quality_and_shared_point_budgets_refuse_without_publishing_partial_edges() {
    let solid = solid();
    let (body, frame) = solid.helical_body().unwrap();
    let (source, _) = solid.rust_phillips_screw_geometry().unwrap();
    for deflection in [0., f64::NAN, 1e-300] {
        let mut remaining = MAX_DISPLAY_EDGE_POINTS;
        assert!(
            super::edges(
                body,
                source,
                frame,
                deflection,
                GeometryPolicy::default(),
                &mut remaining
            )
            .is_err()
        );
    }
    for budget in [0, 1, 2, 12, 64, 128] {
        let mut remaining = budget;
        assert!(
            super::edges(
                body,
                source,
                frame,
                0.005,
                GeometryPolicy::default(),
                &mut remaining
            )
            .is_err(),
            "budget{budget}"
        );
    }
    let mut remaining = MAX_DISPLAY_EDGE_POINTS;
    let far = PlaneFrame {
        origin: [1_000_000.1; 3],
        ..frame
    };
    assert!(
        super::edges(
            body,
            source,
            far,
            1e-10,
            GeometryPolicy::default(),
            &mut remaining
        )
        .is_err()
    );
    let mut remaining = MAX_DISPLAY_EDGE_POINTS;
    assert!(
        super::edges(
            body,
            source,
            PlaneFrame {
                u: [0.; 3],
                ..frame
            },
            0.1,
            GeometryPolicy::default(),
            &mut remaining
        )
        .is_err()
    );
    let mut remaining = MAX_DISPLAY_EDGE_POINTS;
    assert_eq!(
        super::edges(
            body,
            source,
            frame,
            0.1,
            GeometryPolicy::default(),
            &mut remaining
        )
        .unwrap()
        .len(),
        76
    );
}

#[test]
fn conservative_actual_curve_bounds_and_parameter_resolution_cannot_publish_unproved_chords() {
    let solid = solid();
    let (body, frame) = solid.helical_body().unwrap();
    let (source, _) = solid.rust_phillips_screw_geometry().unwrap();
    let cutter = source.cavity().source().boundary_topology().unwrap();
    let (edge, curve) = cutter
        .edges()
        .iter()
        .enumerate()
        .find_map(|(i, _)| match cutter.curve_analytic_support(i).unwrap() {
            Some(curve @ AnalyticCurve::CircleArc { .. }) => Some((i, curve)),
            _ => None,
        })
        .unwrap();
    let query = Query {
        topology: &cutter,
        source,
        frame,
        policy: GeometryPolicy::default(),
        deflection: 0.01,
        operator: operator_bound(frame),
    };
    let bound = analytic_span(curve).second.max(1.);
    // The original full-domain curvature bound remains a valid conservative
    // bound after restricting to either tiny normalized source interval.
    for span in [
        Span {
            domain: [0., 1.],
            second: bound,
            source_error: 0.01,
        },
        Span {
            domain: [f64::from_bits(1), f64::from_bits(2)],
            second: bound,
            source_error: 0.,
        },
        Span {
            domain: [0.5, 0.5_f64.next_up()],
            second: bound,
            source_error: 0.,
        },
    ] {
        let mut points = Vec::new();
        let mut remaining = MAX_DISPLAY_EDGE_POINTS;
        assert!(
            query
                .append(edge, span, &mut points, &mut remaining)
                .is_err()
        );
    }
    let far = Query {
        frame: PlaneFrame {
            origin: [1_000_000.1; 3],
            ..frame
        },
        policy: GeometryPolicy {
            linear_tolerance_mm: 1e-12,
            ..GeometryPolicy::default()
        },
        ..query
    };
    assert!(far.point(edge, 0.).is_err());
    let mut remaining = MAX_DISPLAY_EDGE_POINTS;
    assert!(
        super::edges(
            body,
            source,
            frame,
            f64::from_bits(1),
            GeometryPolicy::default(),
            &mut remaining
        )
        .is_err()
    );
}
#[test]
fn actual_trim_exhaustion_after_valid_prefix_does_not_return_partial_drive_edges() {
    let solid = solid();
    let (body, frame) = solid.helical_body().unwrap();
    let (source, _) = solid.rust_phillips_screw_geometry().unwrap();
    let cutter = source.cavity().source().boundary_topology().unwrap();
    let mut remaining = MAX_DISPLAY_EDGE_POINTS;
    let full = super::edges(
        body,
        source,
        frame,
        0.01,
        GeometryPolicy::default(),
        &mut remaining,
    )
    .unwrap();
    let index = cutter
        .edges()
        .iter()
        .enumerate()
        .find(|(i, e)| {
            matches!(e.source_identity(), CurveKey::Trim { .. }) && full[*i].points_mm.len() > 3
        })
        .unwrap()
        .0;
    let mut budget = full[..index]
        .iter()
        .map(|e| e.points_mm.len())
        .sum::<usize>()
        + full[index].points_mm.len()
        - 1;
    assert!(
        super::edges(
            body,
            source,
            frame,
            0.01,
            GeometryPolicy::default(),
            &mut budget
        )
        .is_err()
    );
    let mut remaining = MAX_DISPLAY_EDGE_POINTS;
    assert_eq!(
        super::edges(
            body,
            source,
            frame,
            0.01,
            GeometryPolicy::default(),
            &mut remaining
        )
        .unwrap()
        .len(),
        76
    );
}
