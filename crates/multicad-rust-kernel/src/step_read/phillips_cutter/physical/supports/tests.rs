use super::*;
#[test]
fn actual_metric_cubic_carriers_have_measured_bounds_against_both_literal_surfaces() {
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
            let now = std::time::Instant::now();
            let bounds = boundary
                .cubic_support_bounds(
                    CubicSupportOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 10_000_000,
                    },
                    None,
                )
                .unwrap();
            eprintln!(
                "index{index} left{left} measuredC0 work{} time{:?} maximum{}",
                bounds.work_used(),
                now.elapsed(),
                bounds
                    .edges()
                    .iter()
                    .flat_map(|e| e.owner_distances_mm)
                    .fold(0_f64, f64::max)
            );
            assert!(std::ptr::eq(bounds.boundary(), &boundary));
            assert_eq!(bounds.edges().len(), 32);
            assert!(
                bounds.retained_capacity_bytes()
                    >= std::mem::size_of_val(&bounds)
                        + 32 * std::mem::size_of::<EdgeSupportBound>()
            );
            assert!(bounds.retained_capacity_bytes() < 4096);
            assert!(bounds.edges().iter().all(|e| {
                e.owner_distances_mm
                    .into_iter()
                    .all(|v| v.is_finite() && v > 0. && v <= 1e-7)
            }));
        }
    }
}
fn axis() -> Axis {
    Axis {
        origin: [100., -20., 60.],
        normal: [0., 0., 2.],
        reference: [3., 0., 0.],
    }
}
#[test]
fn literal_implicit_values_gradients_and_roundtrip_charts_enclose_independent_primitive_oracles() {
    let a = axis();
    for (surface, local, value, gradient) in [
        (Surface::Plane { axis: a }, [2., 1., 0.], 0., [0., 0., 1.]),
        (
            Surface::Cylinder {
                axis: a,
                radius_mm: 2.,
            },
            [2., 0., 0.],
            0.,
            [4., 0., 0.],
        ),
        (
            Surface::Cone {
                axis: a,
                radius_mm: 2.,
                semi_angle_rad: std::f64::consts::FRAC_PI_4,
            },
            [3., 0., 1.],
            0.,
            [6., 0., -6.],
        ),
        (
            Surface::Torus {
                axis: a,
                major_radius_mm: 3.,
                minor_radius_mm: 1.,
            },
            [4., 0., 0.],
            0.,
            [96., 0., 0.],
        ),
    ] {
        let world = std::array::from_fn(|i| a.origin[i] + local[i]);
        let receipt = surface.implicit_enclosure(world.map(|v| [v, v])).unwrap();
        assert!(receipt.value[0] <= value && value <= receipt.value[1]);
        for i in 0..3 {
            assert!(
                receipt.world_gradient[i][0] <= gradient[i]
                    && gradient[i] <= receipt.world_gradient[i][1]
            );
            assert!(receipt.local_point[i][0] <= local[i] && local[i] <= receipt.local_point[i][1]);
            assert!(
                receipt.world_reconstruction[i][0] <= world[i]
                    && world[i] <= receipt.world_reconstruction[i][1]
            );
        }
    }
    let plane = Surface::Plane { axis: a };
    assert!(
        plane
            .implicit_enclosure([[1., 0.], [0., 0.], [0., 0.]])
            .is_err()
    );
    let mut huge = a;
    huge.origin = [-f64::MAX; 3];
    assert!(
        Surface::Plane { axis: huge }
            .implicit_enclosure([[f64::MAX, f64::MAX]; 3])
            .is_err()
    );
    let cone = Surface::Cone {
        axis: a,
        radius_mm: 1.,
        semi_angle_rad: std::f64::consts::FRAC_PI_2,
    };
    assert!(cone.implicit_enclosure([[0., 0.]; 3]).is_err());
    let torus = Surface::Torus {
        axis: a,
        major_radius_mm: 1_f64.next_up(),
        minor_radius_mm: 1.,
    };
    assert!(torus.implicit_enclosure([[0., 0.]; 3]).is_err());
}
#[test]
fn finite_polynomial_support_distance_bounds_cover_known_radial_offsets_and_apex_refusal() {
    let a = Axis {
        origin: [0.; 3],
        normal: [0., 0., 1.],
        reference: [1., 0., 0.],
    };
    for (surface, p, distance) in [
        (Surface::Plane { axis: a }, [0., 0., 0.25], 0.25),
        (
            Surface::Cylinder {
                axis: a,
                radius_mm: 2.,
            },
            [3., 0., 0.],
            1.,
        ),
        (
            Surface::Cone {
                axis: a,
                radius_mm: 2.,
                semi_angle_rad: std::f64::consts::FRAC_PI_4,
            },
            [3.25, 0., 1.],
            0.25,
        ),
        (
            Surface::Torus {
                axis: a,
                major_radius_mm: 3.,
                minor_radius_mm: 1.,
            },
            [4.25, 0., 0.],
            0.25,
        ),
    ] {
        let chart = Chart::new(&surface).unwrap();
        let result = chart.distance([chart.local(p); 4]).unwrap();
        assert!(distance <= result && result < distance + 1e-8);
    }
    let cone = Chart::new(&Surface::Cone {
        axis: a,
        radius_mm: 0.,
        semi_angle_rad: 0.5,
    })
    .unwrap();
    assert!(cone.distance([cone.local([0.; 3]); 4]).is_err());
    let cylinder = Chart::new(&Surface::Cylinder {
        axis: a,
        radius_mm: 0.,
    })
    .unwrap();
    assert!(cylinder.distance([cylinder.local([0.; 3]); 4]).is_err());
    assert!(finite_bound(f64::INFINITY).is_err());
}
#[test]
fn measured_receipts_enforce_precision_resource_and_cancellation_without_minting_uncertainty() {
    let source = super::super::tests::fixture(1, false, false, 1);
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let boundary = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let option = |distance, work| CubicSupportOptions {
        maximum_distance_mm: distance,
        work_limit: work,
    };
    for distance in [0., f64::NAN] {
        assert!(
            boundary
                .cubic_support_bounds(option(distance, 10_000_000), None)
                .is_err()
        );
    }
    assert!(
        boundary
            .cubic_support_bounds(option(1e-7, 0), None)
            .is_err()
    );
    assert_eq!(
        boundary
            .cubic_support_bounds(option(1e-7, 10_000_000), Some(&|| true))
            .unwrap_err(),
        PhillipsReadError::Cancelled
    );
    assert!(
        boundary
            .cubic_support_bounds(option(1e-30, 20_000), None)
            .is_err()
    );
    let edge = boundary
        .edges()
        .iter()
        .find(|e| matches!(e.curve, Curve::CubicIntersection { .. }))
        .unwrap();
    let Curve::CubicIntersection {
        pole_entities,
        spans,
        ..
    } = &edge.curve
    else {
        panic!("cubic selected")
    };
    let pole = pole_entities[1];
    let p = spans[0].poles_mm[1];
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("#{pole}=")))
        .unwrap();
    let changed = text.replacen(
        line,
        &format!(
            "#{pole}=CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}));",
            p[0] + 0.01,
            p[1],
            p[2]
        ),
        1,
    );
    let changed = read_phillips_screw_boundary(&changed, 1e-7, 1_000_000).unwrap();
    assert!(
        changed
            .cubic_support_bounds(option(1e-7, 2000), None)
            .is_err()
    );
}
