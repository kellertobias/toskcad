use super::*;

fn drive(scale: f64) -> Drive {
    Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8 * scale,
            minor_diameter_mm: 2.05 * scale,
            outer_lobe_radius_mm: 0.28 * scale,
        },
        Cone {
            tip_depth_mm: (1.1 + 1.4 / 3_f64.sqrt()) * scale,
            depth_per_radius: 1. / 3_f64.sqrt(),
        },
    )
    .unwrap()
}
fn point(span: CubicIntersectionSpan, t: f64) -> [f64; 3] {
    let a = 1. - t;
    let weights = [a * a * a, 3. * a * a * t, 3. * a * t * t, t * t * t];
    std::array::from_fn(|i| {
        span.poles_mm
            .iter()
            .zip(weights)
            .map(|(p, w)| p[i] * w)
            .sum()
    })
}
#[test]
fn retained_cylinder_cone_carriers_bind_whole_curves_and_keep_literal_shared_contacts() {
    for scale in [0.125, 1., 8.] {
        let source = drive(scale);
        for index in 0..12 {
            let limit = 1e-7 * scale;
            let spans = source.intersection_carrier(index, limit, 10_000).unwrap();
            let arc = source.arcs[index];
            assert_eq!(spans[0].parameters[0], 0.);
            assert_eq!(spans.last().unwrap().parameters[1], 1.);
            let contact = |[x, y]: [f64; 2]| {
                [
                    x,
                    y,
                    source.cone.tip_depth_mm - source.cone.depth_per_radius * x.hypot(y),
                ]
            };
            assert_eq!(spans[0].poles_mm[0], contact(arc.start_mm));
            assert_eq!(spans.last().unwrap().poles_mm[3], contact(arc.end_mm));
            for pair in spans.windows(2) {
                assert_eq!(pair[0].parameters[1], pair[1].parameters[0]);
                assert_eq!(pair[0].poles_mm[3], pair[1].poles_mm[0]);
            }
            for span in spans {
                let proof = source.certify_intersection_span(index, &span).unwrap();
                assert!(proof.total_error_bound_mm() <= limit);
                assert_eq!(
                    proof.source_error_mm,
                    source.spatial_source_error_bound_mm()
                );
                for i in 0..=100 {
                    let t = f64::from(i) / 100.;
                    let u = span.parameters[0] + (span.parameters[1] - span.parameters[0]) * t;
                    let theta = arc.start_angle_rad + arc.sweep_rad * u;
                    let x = arc.center_mm[0] + arc.radius_mm * theta.cos();
                    let y = arc.center_mm[1] + arc.radius_mm * theta.sin();
                    let exact = [
                        x,
                        y,
                        source.cone.tip_depth_mm - source.cone.depth_per_radius * x.hypot(y),
                    ];
                    let p = point(span, t);
                    let distance = p
                        .into_iter()
                        .zip(exact)
                        .map(|(a, b)| (a - b).powi(2))
                        .sum::<f64>()
                        .sqrt();
                    assert!(distance <= proof.total_error_bound_mm());
                    let cylinder = (p[0] - arc.center_mm[0]).hypot(p[1] - arc.center_mm[1]);
                    assert!((cylinder - arc.radius_mm).abs() <= proof.total_error_bound_mm());
                    let cone =
                        source.cone.tip_depth_mm - source.cone.depth_per_radius * p[0].hypot(p[1]);
                    assert!(
                        (p[2] - cone).abs()
                            <= (1. + source.cone.depth_per_radius) * proof.total_error_bound_mm()
                    );
                }
            }
        }
    }
}

#[test]
fn literal_carrier_tampering_and_work_failures_remain_visible_without_source_reconstruction() {
    let source = drive(1.);
    let snapshot = source.clone();
    let original = source.intersection_carrier(1, 1e-7, 10_000).unwrap()[0];
    for axis in 0..3 {
        let mut changed = original;
        changed.poles_mm[1][axis] += 0.01;
        assert!(
            source
                .certify_intersection_span(1, &changed)
                .unwrap()
                .total_error_bound_mm()
                > 0.005
        );
    }
    for parameters in [[f64::NAN, 1.], [-0.1, 1.], [0., 1.1], [0.5, 0.5], [1., 0.]] {
        assert_eq!(
            source.certify_intersection_span(
                1,
                &CubicIntersectionSpan {
                    parameters,
                    ..original
                }
            ),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        source.certify_intersection_span(12, &original),
        Err(Error::Invalid)
    );
    let mut nonfinite = original;
    nonfinite.poles_mm[0][0] = f64::NAN;
    assert_eq!(
        source.certify_intersection_span(1, &nonfinite),
        Err(Error::Invalid)
    );
    let mut huge = original;
    huge.poles_mm[0][0] = 1e10;
    assert_eq!(
        source.certify_intersection_span(1, &huge),
        Err(Error::Limit)
    );
    for width in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            source.intersection_carrier(0, width, 10_000),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        source.intersection_carrier(12, 1e-7, 10_000),
        Err(Error::Invalid)
    );
    for limit in [0, 1, 1_000_001] {
        assert_eq!(
            source.intersection_carrier(0, 1e-7, limit),
            Err(Error::Limit)
        );
    }
    assert_eq!(
        source.intersection_carrier(0, source.spatial_source_error_bound_mm(), 10_000),
        Err(Error::Unresolved)
    );
    assert_eq!(
        source.intersection_carrier_with_cancel(0, 1e-7, 10_000, &|| true),
        Err(Error::Cancelled)
    );
    assert_eq!(source, snapshot);
    assert!(source.intersection_carrier(0, 1e-7, 10_000).is_ok());
}

#[test]
fn ill_conditioned_retained_supports_refuse_a_carrier_instead_of_hiding_chart_uncertainty() {
    let points: [[f64; 2]; 12] = std::array::from_fn(|i| {
        let t = i as f64 * TAU / 12.;
        [t.cos(), t.sin()]
    });
    let arcs = std::array::from_fn(|i| {
        let t = (i as f64 + 0.5) * TAU / 12.;
        let c = if i.is_multiple_of(2) {
            [0.; 2]
        } else {
            [1e6 * t.cos(), 1e6 * t.sin()]
        };
        let from = points[i];
        let to = points[(i + 1) % 12];
        let a = (from[1] - c[1]).atan2(from[0] - c[0]);
        let b = (to[1] - c[1]).atan2(to[0] - c[0]);
        Arc {
            center_mm: c,
            radius_mm: (from[0] - c[0]).hypot(from[1] - c[1]),
            start_mm: from,
            end_mm: to,
            start_angle_rad: a,
            sweep_rad: if i.is_multiple_of(2) {
                (b - a).rem_euclid(TAU)
            } else {
                -(a - b).rem_euclid(TAU)
            },
        }
    });
    let source = Drive::from_retained(
        arcs,
        Cone {
            tip_depth_mm: 10.,
            depth_per_radius: 1.,
        },
        TransportError::ZERO,
    )
    .unwrap();
    assert_eq!(
        source.intersection_carrier(1, 8. * source.spatial_source_error_bound_mm(), 10_000),
        Err(Error::Unresolved)
    );
    assert!(source.intersection_carrier(0, 1e-5, 10_000).is_ok());
}
