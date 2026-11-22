use super::*;
#[test]
fn actual_metric_analytic_curves_bind_whole_supports_without_an_authored_recipe() {
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
            let endpoints = boundary
                .endpoint_coupling_bounds(
                    EndpointCouplingOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 10_000,
                    },
                    None,
                )
                .unwrap();
            let now = std::time::Instant::now();
            let receipt = boundary
                .analytic_support_bounds(
                    &endpoints,
                    AnalyticSupportOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 1_000_000,
                    },
                    None,
                )
                .unwrap_or_else(|e| panic!("index{index} left{left}: {e:?}"));
            assert!(std::ptr::eq(receipt.boundary(), &boundary));
            assert_eq!(receipt.edges().len(), 44);
            assert!(receipt.retained_capacity_bytes() < 4096);
            eprintln!(
                "index{index} left{left} analyticbounds{} work{} time{:?}",
                receipt
                    .edges()
                    .iter()
                    .flat_map(|e| e.owner_distances_mm)
                    .fold(0_f64, f64::max),
                receipt.work_used(),
                now.elapsed()
            );
        }
    }
}
#[test]
fn far_literal_frames_keep_the_same_independent_analytic_contract() {
    for index in [0, 5] {
        let source = super::super::tests::fixture(index, true, false, 1);
        let text = crate::step::write_phillips_screw(
            &source,
            [1e4, -2e4, 3e4],
            [0.8, 0.6, 0.],
            [0.36, -0.48, -0.8],
            1e-7,
        )
        .unwrap();
        let b = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
        let e = b
            .endpoint_coupling_bounds(
                EndpointCouplingOptions {
                    maximum_distance_mm: 1e-7,
                    work_limit: 10_000,
                },
                None,
            )
            .unwrap();
        let a = b
            .analytic_support_bounds(
                &e,
                AnalyticSupportOptions {
                    maximum_distance_mm: 1e-7,
                    work_limit: 1_000_000,
                },
                None,
            )
            .unwrap();
        eprintln!(
            "farindex{index} analyticbound{}",
            a.edges()
                .iter()
                .flat_map(|e| e.owner_distances_mm)
                .fold(0_f64, f64::max)
        );
    }
}
#[test]
fn malformed_owner_supports_cannot_borrow_a_large_requested_coupling_error() {
    let source = super::super::tests::fixture(0, false, false, 1);
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let b = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let e = b
        .endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000,
            },
            None,
        )
        .unwrap();
    let options = |distance, work| AnalyticSupportOptions {
        maximum_distance_mm: distance,
        work_limit: work,
    };
    let cloned = b.clone();
    assert!(
        cloned
            .analytic_support_bounds(&e, options(1e-7, 1_000_000), None)
            .is_err()
    );
    for distance in [0., f64::NAN, 1e-30] {
        assert!(
            b.analytic_support_bounds(&e, options(distance, 1_000_000), None)
                .is_err()
        );
    }
    assert!(
        b.analytic_support_bounds(&e, options(1e-7, 0), None)
            .is_err()
    );
    assert_eq!(
        b.analytic_support_bounds(&e, options(1e-7, 1_000_000), Some(&|| true))
            .unwrap_err(),
        PhillipsReadError::Cancelled
    );
    let face = b
        .faces()
        .iter()
        .find(|f| matches!(f.support, Surface::Cylinder { .. }))
        .unwrap();
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{}=", face.surface_entity)))
        .unwrap();
    let prefix = old.rsplit_once(',').unwrap().0;
    let Surface::Cylinder { radius_mm, .. } = face.support else {
        panic!("cylinder selected")
    };
    let changed = text.replacen(old, &format!("{prefix},{:.16E});", radius_mm + 0.001), 1);
    let changed = read_phillips_screw_boundary(&changed, 1e-7, 1_000_000).unwrap();
    let e = changed
        .endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000,
            },
            None,
        )
        .unwrap();
    assert!(
        changed
            .analytic_support_bounds(&e, options(100., 1_000_000), None)
            .is_err()
    );
}
#[test]
fn exact_rational_primitive_values_and_cone_ruling_bounds_have_independent_oracles() {
    let axis = Axis {
        origin: [0.; 3],
        normal: [0., 0., 1.],
        reference: [1., 0., 0.],
    };
    let mut work = Work {
        remaining: 1_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let cylinder = Chart::new(&Surface::Cylinder {
        axis,
        radius_mm: 2.,
    })
    .unwrap();
    assert!(arc(axis, [2., 2.], &cylinder, &mut work).unwrap() < 1e-10);
    assert!(arc(axis, [2.01, 2.01], &cylinder, &mut work).is_err());
    assert!(!hull(&[I::point(0.), I::new(f64::NAN, f64::NAN)]).valid());
    let cone = Chart::new(&Surface::Cone {
        axis,
        radius_mm: 1.,
        semi_angle_rad: std::f64::consts::FRAC_PI_4,
    })
    .unwrap();
    // Source points are a true straight generator, including its apex.
    let p = [[0., 0., -1.], [2., 0., 1.]].map(|p| p.map(I::point));
    assert!(cone.cone_line_distance(p).unwrap() < 1e-10);
    let shifted = p.map(|p| [p[0] + I::point(0.01), p[1], p[2]]);
    let bound = cone.cone_line_distance(shifted).unwrap();
    assert!((0.01..0.02).contains(&bound));
    let plane = Chart::new(&Surface::Plane { axis }).unwrap();
    assert!(plane.cone_line_distance(p).is_err());
    assert!(cone.cone_line_distance([[I::point(0.); 3]; 2]).is_err());
    assert!(
        cone.cone_line_distance([[I::point(f64::MAX); 3]; 2])
            .is_err()
    );
    let point = [I::point(1.), I::point(0.), I::point(f64::MAX)];
    assert!(cone.cone_line_distance([point; 2]).is_err());
    let tiny = Chart::new(&Surface::Cone {
        axis,
        radius_mm: 1.,
        semi_angle_rad: 1e-18,
    })
    .unwrap();
    assert!(tiny.cone_line_distance(p).is_err());
    let large = Chart::new(&Surface::Cone {
        axis,
        radius_mm: 1e307,
        semi_angle_rad: 1e-14,
    })
    .unwrap();
    assert!(
        large
            .cone_line_distance([[I::point(1.), I::point(0.), I::point(0.)]; 2])
            .is_err()
    );
    assert!(
        cone.uniform_distance(I::point(0.), [I::point(0.), I::point(0.), I::point(-2.)])
            .is_err()
    );
    let zero = Chart::new(&Surface::Cylinder {
        axis,
        radius_mm: 0.,
    })
    .unwrap();
    assert!(
        zero.uniform_distance(I::point(0.), [I::point(0.); 3])
            .is_err()
    );
    let zero_line = line(
        [0.; 3],
        [1., 0., 0.],
        0.,
        [[0.; 3], [1., 0., 0.]],
        &plane,
        &mut work,
        false,
    );
    assert!(zero_line.is_err());
}
