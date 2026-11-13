use super::*;
#[test]
fn actual_metric_vertices_couple_to_literal_curves_without_any_requested_residual_expansion() {
    for index in 0..6 {
        for left in [false, true] {
            let source = super::super::tests::fixture(index, left, false, 1);
            let (o, u, v) = if left {
                ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
            } else {
                ([0.; 3], [1., 0., 0.], [0., 1., 0.])
            };
            let text = crate::step::write_phillips_screw(&source, o, u, v, 1e-7).unwrap();
            let boundary = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
            let receipt = boundary
                .endpoint_coupling_bounds(
                    EndpointCouplingOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 10_000,
                    },
                    None,
                )
                .unwrap_or_else(|e| panic!("index{index} left{left}: {e:?}"));
            assert!(std::ptr::eq(receipt.boundary(), &boundary));
            assert_eq!(receipt.edge_distances_mm().len(), 76);
            assert_eq!(receipt.work_used(), 608);
            assert!(receipt.retained_capacity_bytes() < 4096);
            eprintln!(
                "index{index} left{left} numericEndpointBound{}",
                receipt
                    .edge_distances_mm()
                    .iter()
                    .flatten()
                    .copied()
                    .fold(0_f64, f64::max)
            );
        }
    }
}
#[test]
fn a_large_requested_distance_cannot_hide_nonzero_literal_endpoint_residuals() {
    let source = super::super::tests::fixture(0, false, false, 1);
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let boundary = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let edge = boundary
        .edges()
        .iter()
        .find(|e| matches!(e.curve, Curve::Circle { .. }))
        .unwrap();
    let vertex = &boundary.vertices()[edge.vertices[0]];
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("#{}=", vertex.point_entity)))
        .unwrap();
    let p = vertex.point;
    let replacement = format!(
        "#{}=CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}));",
        vertex.point_entity,
        p[0] + 0.0001,
        p[1],
        p[2]
    );
    let changed =
        read_phillips_screw_boundary(&text.replacen(line, &replacement, 1), 1e-7, 1_000_000)
            .unwrap();
    assert!(
        changed
            .endpoint_coupling_bounds(
                EndpointCouplingOptions {
                    maximum_distance_mm: 100.,
                    work_limit: 10_000
                },
                None
            )
            .is_err()
    );
    let options = |distance, work| EndpointCouplingOptions {
        maximum_distance_mm: distance,
        work_limit: work,
    };
    for distance in [0., f64::NAN, 1e-30] {
        assert!(
            boundary
                .endpoint_coupling_bounds(options(distance, 10_000), None)
                .is_err()
        );
    }
    assert!(
        boundary
            .endpoint_coupling_bounds(options(1e-7, 0), None)
            .is_err()
    );
    assert_eq!(
        boundary
            .endpoint_coupling_bounds(options(1e-7, 10_000), Some(&|| true))
            .unwrap_err(),
        PhillipsReadError::Cancelled
    );
}
#[test]
fn ill_conditioned_primitive_couplings_refuse_instead_of_claiming_tiny_reserves() {
    let axis = Axis {
        origin: [0.; 3],
        normal: [0., 0., 1.],
        reference: [1., 0., 0.],
    };
    let line = Curve::Line {
        origin: [0.; 3],
        direction: [f64::from_bits(1), 0., 0.],
        magnitude: 1.,
        vector_entity: 0,
        direction_entity: 0,
    };
    assert!(analytic_endpoint(&line, [0.; 3]).is_err());
    let ellipse = Curve::Ellipse {
        axis,
        semiaxes_mm: [f64::from_bits(1); 2],
    };
    assert!(analytic_endpoint(&ellipse, [0.; 3]).is_err());
    let ellipse = Curve::Ellipse {
        axis,
        semiaxes_mm: [1e308; 2],
    };
    assert!(analytic_endpoint(&ellipse, [1e308, 0., 0.]).is_err());
    let mut bad = axis;
    bad.normal = [0.; 3];
    assert!(
        analytic_endpoint(
            &Curve::Circle {
                axis: bad,
                radius_mm: 1.
            },
            [1., 0., 0.]
        )
        .is_err()
    );
}

#[test]
fn an_arbitrary_cubic_point_query_cannot_borrow_the_arena_alias_certificate() {
    let curve = Curve::CubicIntersection {
        spline_entity: 0,
        associated_surface_entities: [0; 2],
        pole_entities: vec![],
        spans: vec![],
    };
    assert!(analytic_endpoint(&curve, [123.; 3]).is_err());
}
