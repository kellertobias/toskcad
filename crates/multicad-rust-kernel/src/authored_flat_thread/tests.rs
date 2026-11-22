use super::*;
use std::cell::Cell;
fn spec(left: bool, starts: u32) -> Spec {
    let pitch_mid = (2.675 + 2.775) / 2.;
    let minor_mid = (2.459 + 2.599) / 2.;
    Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.3,
        depth_mm: (3. - minor_mid) / 2.,
        crest_width_mm: 0.25 - (3. - pitch_mid) / 3_f64.sqrt(),
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts,
        left_hand: left,
        root_radius_mm: None,
    }
}
fn options() -> SourceOptions {
    SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    }
}
fn norm(p: [f64; 3]) -> f64 {
    p[0].portable_hypot(p[1]).portable_hypot(p[2])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[test]
fn canonical_shared_sources_retain_intent_and_exact_common_endpoint_expressions() {
    for left in [false, true] {
        for starts in [1, 2, 3, 8] {
            let request = spec(left, starts);
            let source = AuthoredThread::new(request, options()).unwrap();
            assert_eq!(source.spec(), request);
            let (legacy, receipt) = source.legacy_approximation();
            assert_eq!(legacy.spec(), request);
            assert_eq!(receipt, source.transfer_bounds());
            assert!(
                receipt.maximum_position_error_mm > 0. && receipt.maximum_position_error_mm < 1e-7
            );
            assert!(
                receipt.maximum_normal_angle_error_rad > 0.
                    && receipt.maximum_normal_angle_error_rad < 1e-7
            );
            assert_eq!(
                source.minimum_radius_mm(),
                request.major_radius_mm - request.depth_mm
            );
            assert_eq!(source.maximum_radius_mm(), request.major_radius_mm);
            assert_eq!(source.full_length_mm(), request.length_mm);
            let bands = source.bands();
            for i in 0..4 {
                assert_eq!(bands[i].phase_mm[1], bands[i + 1].phase_mm[0]);
                assert_eq!(bands[i].radius_mm[1], bands[i + 1].radius_mm[0]);
                for z in [0., 0.37, request.length_mm] {
                    assert_eq!(
                        source.point_in_band(i, 0, 1., z).unwrap().point_mm,
                        source.point_in_band(i + 1, 0, 0., z).unwrap().point_mm
                    );
                }
            }
        }
    }
}
#[test]
fn actual_stored_affine_discrepancy_and_whole_band_angle_are_inside_derived_receipts() {
    for left in [false, true] {
        for starts in [1, 3] {
            let source = AuthoredThread::new(spec(left, starts), options()).unwrap();
            let (legacy, receipt) = source.legacy_approximation();
            let spec = source.spec();
            for (band, b) in source.bands().into_iter().enumerate() {
                let slope = (b.radius_mm[1] - b.radius_mm[0]) / (b.phase_mm[1] - b.phase_mm[0]);
                let old_slope = legacy.slope_at_phase(b.phase_mm[0], band);
                assert!((slope - old_slope).abs() <= receipt.bands[band].slope_error);
                for u in [0., 0.13, 0.5, 0.91, 1.] {
                    for z in [0., 0.47, spec.length_mm] {
                        let q = (b.phase_mm[1] - b.phase_mm[0]).mul_add(u, b.phase_mm[0]);
                        let r = (b.radius_mm[1] - b.radius_mm[0]).mul_add(u, b.radius_mm[0]);
                        let anchor = if band == 1 {
                            legacy.profile_breaks()[1]
                        } else {
                            legacy.profile_breaks()[4]
                        };
                        let old_r = if band == 0 || band == 4 {
                            spec.major_radius_mm
                        } else if band == 2 {
                            legacy.root_radius_mm()
                        } else {
                            old_slope.mul_add(q - anchor, spec.major_radius_mm)
                        };
                        let old_theta = (z - q) / legacy.helical_rate();
                        let sample = source.point_in_band(band, 0, u, z).unwrap();
                        let old_point = [
                            old_r * old_theta.portable_cos(),
                            old_r * old_theta.portable_sin(),
                            z,
                        ];
                        assert!(
                            norm(std::array::from_fn(|i| old_point[i] - sample.point_mm[i]))
                                <= receipt.bands[band].position_error_mm + sample.error_mm
                        );
                        let theta = (if left { -1. } else { 1. }) * std::f64::consts::TAU * (z - q)
                            / (spec.pitch_mm * f64::from(starts));
                        let rate =
                            (if left { -1. } else { 1. }) * spec.pitch_mm * f64::from(starts)
                                / std::f64::consts::TAU;
                        let normal = |t: f64, beta: f64, d: f64| {
                            [
                                t.portable_cos() - beta * t.portable_sin(),
                                t.portable_sin() + beta * t.portable_cos(),
                                -d,
                            ]
                        };
                        let n = normal(theta, rate * slope / r, slope);
                        let o = normal(
                            old_theta,
                            legacy.helical_rate() * old_slope / old_r,
                            old_slope,
                        );
                        let angle = norm(cross(n, o))
                            .portable_atan2(n.into_iter().zip(o).map(|(x, y)| x * y).sum());
                        assert!(
                            angle <= receipt.bands[band].normal_angle_error_rad,
                            "{angle} {receipt:?}"
                        );
                    }
                }
            }
        }
    }
}
#[test]
fn constructor_budget_cancellation_and_precision_refusals_are_atomic() {
    let original = spec(false, 1);
    for invalid in [
        SourceOptions {
            position_limit_mm: 0.,
            ..options()
        },
        SourceOptions {
            normal_angle_limit_rad: f64::NAN,
            ..options()
        },
    ] {
        assert_eq!(AuthoredThread::new(original, invalid), Err(Error::Invalid));
    }
    assert_eq!(
        AuthoredThread::new(
            Spec {
                root_radius_mm: Some(0.05),
                ..original
            },
            options()
        ),
        Err(Error::Unsupported)
    );
    assert_eq!(
        AuthoredThread::new(
            Spec {
                pitch_mm: -1.,
                ..original
            },
            options()
        ),
        Err(Error::Invalid)
    );
    assert_eq!(
        AuthoredThread::new(
            Spec {
                starts: 9,
                ..original
            },
            options()
        ),
        Err(Error::Limit)
    );
    assert_eq!(
        AuthoredThread::new(
            original,
            SourceOptions {
                position_limit_mm: 1e-300,
                ..options()
            }
        ),
        Err(Error::Unresolved)
    );
    let tiny = Spec {
        major_radius_mm: 1.,
        depth_mm: 1_f64.next_down(),
        pitch_mm: 4.,
        length_mm: 20.,
        crest_width_mm: 1.,
        ..original
    };
    assert_eq!(AuthoredThread::new(tiny, options()), Err(Error::Unresolved));
    let calls = Cell::new(0);
    let stop = || {
        calls.set(calls.get() + 1);
        calls.get() >= 4
    };
    assert_eq!(
        AuthoredThread::new_with_interrupt(original, options(), Some(&stop)),
        Err(Error::Cancelled)
    );
    assert_eq!(calls.get(), 4);
    let source = AuthoredThread::new(original, options()).unwrap();
    for (band, start, u, z) in [
        (5, 0, 0., 0.),
        (0, 1, 0., 0.),
        (0, 0, -1., 0.),
        (0, 0, 0., 4.),
    ] {
        assert_eq!(source.point_in_band(band, start, u, z), Err(Error::Invalid));
    }
}
#[test]
fn exact_turn_chart_and_large_scale_preserve_common_source_joins_with_finite_evaluation_error() {
    for scale in [0.1, 1., 1000.] {
        let base = spec(true, 3);
        let request = Spec {
            major_radius_mm: base.major_radius_mm * scale,
            pitch_mm: base.pitch_mm * scale,
            length_mm: base.length_mm * scale,
            depth_mm: base.depth_mm * scale,
            crest_width_mm: base.crest_width_mm * scale,
            ..base
        };
        let source = AuthoredThread::new(
            request,
            SourceOptions {
                position_limit_mm: 1e-7 * scale,
                ..options()
            },
        )
        .unwrap();
        for z in [0., request.length_mm] {
            let a = source.point_in_band(0, 0, 0., z).unwrap();
            let b = source.point_in_band(4, request.starts - 1, 1., z).unwrap();
            assert!(
                norm(std::array::from_fn(|i| a.point_mm[i] - b.point_mm[i]))
                    <= a.error_mm + b.error_mm
            );
        }
    }
}

#[test]
fn positive_but_ill_conditioned_legacy_root_refuses_unprovable_normal_transport() {
    let source = Spec {
        major_radius_mm: 1.,
        depth_mm: 1. - 1e-8,
        pitch_mm: 4.,
        length_mm: 20.,
        crest_width_mm: 1.,
        ..spec(false, 1)
    };
    // The legacy source is valid and its end radii remain strictly positive;
    // division by the very small radius prevents the requested normal proof.
    assert!(Thread::new(source).is_ok());
    assert_eq!(
        AuthoredThread::new(source, options()),
        Err(Error::Unresolved)
    );
}
