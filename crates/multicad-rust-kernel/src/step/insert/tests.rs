use super::*;
use multicad_precision::PortableMath as _;
#[test]
fn invalid_and_cancelled_requests_fail_before_publication() {
    assert!(matches!(Work::new(0, None), Err(Error::Invalid)));
    assert!(matches!(Work::new(16_000_001, None), Err(Error::Invalid)));
    let stop = || true;
    let mut w = Work::new(10, Some(&stop)).unwrap();
    assert_eq!(w.charge(1), Err(Error::Cancelled));
    let mut w = Work::new(1, None).unwrap();
    assert_eq!(w.charge(2), Err(Error::Limit));
    assert!(matches!(Memory::new(2, 1), Err(Error::Limit)));
}
pub(super) fn graph(left: bool, starts: u32) -> Boundary {
    use crate::{
        authored_flat_thread::{AuthoredThread, SourceOptions},
        helical_thread::Spec,
        wood_insert::{BoundaryOptions, Insert},
    };
    let spec = Spec {
        major_radius_mm: 5.,
        pitch_mm: 1.9,
        length_mm: 4.7,
        depth_mm: 1.,
        crest_width_mm: 0.4,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts,
        left_hand: left,
        root_radius_mm: None,
    };
    let options = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let outer = AuthoredThread::new(spec, options).unwrap();
    let inner = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            depth_mm: 0.3,
            crest_width_mm: 0.1,
            starts: 1,
            ..spec
        },
        options,
    )
    .unwrap();
    Insert::new(outer, inner, 0.5, None)
        .unwrap()
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap()
}
#[test]
fn actual_paired_curve_carriers_share_endpoints_and_bound_tangents() {
    for (left, starts) in [(false, 1), (true, 1), (false, 2)] {
        let body = graph(left, starts);
        let plan = CurvePlan::prepare(&body, CurveOptions::default(), None).unwrap();
        assert!(std::ptr::eq(plan.boundary(), &body));
        assert_eq!(
            plan.curves().len(),
            body.source_boundary(Side::Outer).source_curve_count()
                + body.source_boundary(Side::Inner).source_curve_count()
        );
        for curve in plan.curves() {
            assert_eq!(
                curve.poles()[0],
                vertex(plan.vertices(), &body, curve.endpoints()[0])
            );
            assert_eq!(
                *curve.poles().last().unwrap(),
                vertex(plan.vertices(), &body, curve.endpoints()[1])
            );
            assert!(curve.position_error_mm() <= 1e-7);
            assert!(curve.derivative_error_mm() <= 1e-6);
            assert_eq!(curve.knots()[0], 0.);
            assert_eq!(*curve.knots().last().unwrap(), 1.);
            assert!(matches!(curve.degree(), 1 | 11));
        }
        assert!(plan.work_used() <= 16_000_000);
        assert!(plan.retained_bytes() <= 64 * 1024 * 1024);
    }
}
fn bezier(mut values: Vec<[f64; 3]>, t: f64) -> [f64; 3] {
    for n in (1..values.len()).rev() {
        for i in 0..n {
            let next = values[i + 1];
            for (value, other) in values[i].iter_mut().zip(next) {
                *value = *value * (1. - t) + other * t;
            }
        }
    }
    values[0]
}
#[test]
fn independent_sampled_carriers_and_edited_interior_poles_have_measured_transfer() {
    let body = graph(true, 1);
    let plan = CurvePlan::prepare(&body, CurveOptions::default(), None).unwrap();
    let mut edited = false;
    for (index, curve) in plan.curves().iter().enumerate() {
        let source = body.curve(curve.source_key()).unwrap();
        let [a, d] = source
            .affine_rzv_bounds()
            .map(|v| v.map(|b| (b.lower + b.upper) * 0.5));
        let sign = if body
            .source_boundary(curve.source_key().owner())
            .source()
            .spec()
            .left_hand
        {
            -1.
        } else {
            1.
        };
        for (span, pair) in curve.knots().windows(2).enumerate() {
            let degree = curve.degree();
            let controls = &curve.poles()[span * degree..=span * degree + degree];
            for i in 0..=16 {
                let u = i as f64 / 16.;
                let t = pair[0] + u * (pair[1] - pair[0]);
                let r = a[0] + d[0] * t;
                let theta = sign * std::f64::consts::TAU * (a[2] + d[2] * t);
                let source = [
                    r * theta.portable_cos(),
                    r * theta.portable_sin(),
                    a[1] + d[1] * t,
                ];
                let value = bezier(controls.to_vec(), u);
                let error = (0..3).map(|k| (value[k] - source[k]).abs()).sum::<f64>();
                assert!(
                    error <= curve.position_error_mm() + 1e-12,
                    "{error} > {}",
                    curve.position_error_mm()
                );
                let tangent = [
                    d[0] * theta.portable_cos()
                        - sign * std::f64::consts::TAU * d[2] * r * theta.portable_sin(),
                    d[0] * theta.portable_sin()
                        + sign * std::f64::consts::TAU * d[2] * r * theta.portable_cos(),
                    d[1],
                ];
                let derivative = bezier(
                    controls
                        .windows(2)
                        .map(|p| {
                            std::array::from_fn(|k| {
                                degree as f64 * (p[1][k] - p[0][k]) / (pair[1] - pair[0])
                            })
                        })
                        .collect(),
                    u,
                );
                let error = (0..3)
                    .map(|k| (derivative[k] - tangent[k]).abs())
                    .sum::<f64>();
                assert!(error <= curve.derivative_error_mm() + 1e-11);
            }
        }
        let receipt = plan
            .certify_controls(index, curve.poles(), 16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        assert!(receipt.position_error_mm <= 1e-7);
        assert!(receipt.derivative_error_mm <= 1e-6);
        assert!(receipt.work_used > 0);
        if !edited && curve.degree() == 11 {
            let mut poles = curve.poles().to_vec();
            poles[3][0] += 0.001;
            let bad = plan
                .certify_controls(index, &poles, 16_000_000, 64 * 1024 * 1024, None)
                .unwrap();
            assert!(bad.position_error_mm > 1e-4);
            assert!(bad.derivative_error_mm > 1e-3);
            edited = true;
        }
    }
    assert!(edited);
}
#[test]
fn bounded_curve_queries_refuse_invalid_precision_resources_and_cancellation() {
    let body = graph(false, 1);
    for options in [
        CurveOptions {
            position_error_mm: 0.,
            ..CurveOptions::default()
        },
        CurveOptions {
            derivative_error_mm: f64::NAN,
            ..CurveOptions::default()
        },
        CurveOptions {
            maximum_spans_per_curve: 0,
            ..CurveOptions::default()
        },
        CurveOptions {
            maximum_spans_per_curve: 1025,
            ..CurveOptions::default()
        },
        CurveOptions {
            bytes: 0,
            ..CurveOptions::default()
        },
        CurveOptions {
            bytes: 64 * 1024 * 1024 + 1,
            ..CurveOptions::default()
        },
    ] {
        assert!(matches!(
            CurvePlan::prepare(&body, options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        CurveOptions {
            work: 1,
            ..CurveOptions::default()
        },
        CurveOptions {
            bytes: 1,
            ..CurveOptions::default()
        },
    ] {
        assert!(matches!(
            CurvePlan::prepare(&body, options, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        CurvePlan::prepare(
            &body,
            CurveOptions {
                position_error_mm: 1e-100,
                maximum_spans_per_curve: 1,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        CurvePlan::prepare(&body, CurveOptions::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let plan = CurvePlan::prepare(&body, CurveOptions::default(), None).unwrap();
    assert!(matches!(
        plan.certify_controls(usize::MAX, &[], 1, 64 * 1024 * 1024, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        plan.certify_controls(0, &[], 100, 64 * 1024 * 1024, None),
        Err(Error::Invalid)
    ));
    let mut poles = plan.curves()[0].poles().to_vec();
    poles[0][0] += 1.;
    assert!(matches!(
        plan.certify_controls(0, &poles, 100, 64 * 1024 * 1024, None),
        Err(Error::Unresolved)
    ));
    poles[0][0] = f64::INFINITY;
    assert!(matches!(
        plan.certify_controls(0, &poles, 100, 64 * 1024 * 1024, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        plan.certify_controls(
            0,
            plan.curves()[0].poles(),
            100,
            64 * 1024 * 1024,
            Some(&|| true)
        ),
        Err(Error::Cancelled)
    ));
}
pub(super) fn metric_graph() -> Boundary {
    use crate::{
        authored_flat_thread::{AuthoredThread, SourceOptions},
        helical_thread::Spec,
        wood_insert::{BoundaryOptions, Insert},
    };
    let source = |radius, pitch, depth, crest| {
        AuthoredThread::new(
            Spec {
                major_radius_mm: radius,
                pitch_mm: pitch,
                length_mm: 4.,
                depth_mm: depth,
                crest_width_mm: crest,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                starts: 1,
                left_hand: false,
                root_radius_mm: None,
            },
            SourceOptions {
                position_limit_mm: 1e-8,
                normal_angle_limit_rad: 1e-6,
            },
        )
        .unwrap()
    };
    Insert::new(
        source(5., 2., 1., 0.4),
        source(
            1.5,
            0.5,
            (3. - (2.459 + 2.599) / 2.) / 2.,
            0.25 - (3. - (2.675 + 2.775) / 2.) / 3_f64.sqrt(),
        ),
        0.5,
        None,
    )
    .unwrap()
    .boundary_with_interrupt(BoundaryOptions::default(), None)
    .unwrap()
}
#[test]
fn real_metric_fixture_has_bounded_curve_plan_and_planar_cap_controls() {
    let body = metric_graph();
    let plan = CurvePlan::prepare(&body, CurveOptions::default(), None).unwrap();
    for curve in plan.curves() {
        if let Some(cap) = body
            .curve(curve.source_key())
            .unwrap()
            .faces()
            .into_iter()
            .find(|&f| f < 2)
        {
            let z = if cap == 0 { 0. } else { 4. };
            assert!(curve.poles().iter().all(|p| p[2] == z));
        }
    }
    eprintln!(
        "metric carrier plan: {} curves {} work {} bytes",
        plan.curves().len(),
        plan.work_used(),
        plan.retained_bytes()
    );
}
#[test]
fn supplied_large_literal_controls_fail_closed_and_budget_growth_is_bounded() {
    let body = graph(false, 1);
    let plan = CurvePlan::prepare(&body, CurveOptions::default(), None).unwrap();
    let (index, curve) = plan
        .curves()
        .iter()
        .enumerate()
        .find(|(_, c)| c.degree() == 11)
        .unwrap();
    let mut poles = curve.poles().to_vec();
    poles[3][0] = f64::MAX;
    assert!(matches!(
        plan.certify_controls(index, &poles, 16_000_000, 64 * 1024 * 1024, None),
        Err(Error::Unresolved)
    ));
    let mut memory = Memory::new(0, 1).unwrap();
    assert_eq!(memory.add(2), Err(Error::Limit));
    let mut work = Work::new(16_000_000, None).unwrap();
    assert!(matches!(
        memory.allocate::<[f64; 3]>(1, &mut work),
        Err(Error::Limit)
    ));
    assert!(matches!(
        CurvePlan::prepare(
            &body,
            CurveOptions {
                bytes: body.retained_bytes() + 1,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
}
#[test]
fn actual_affine_source_precision_and_actual_allocated_capacity_remain_checked() {
    let body = graph(false, 1);
    let plan = CurvePlan::prepare(&body, CurveOptions::default(), None).unwrap();
    let curve = plan.curves().iter().find(|c| c.degree() == 1).unwrap();
    let source = body.curve(curve.source_key()).unwrap();
    let mut work = Work::new(16_000_000, None).unwrap();
    let basis = curves::derivative_basis(&mut work).unwrap();
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let ends = [curve.poles()[0], *curve.poles().last().unwrap()];
    assert!(matches!(
        curves::prepare(
            curve.source_key(),
            curve.endpoints(),
            source,
            false,
            ends,
            None,
            CurveOptions {
                position_error_mm: 1e-100,
                ..CurveOptions::default()
            },
            &basis,
            &mut memory,
            &mut work
        ),
        Err(Error::Unresolved)
    ));
    let mut memory = Memory::new(0, 1).unwrap();
    let allocated = vec![0_f64; 10];
    assert!(matches!(memory.account(allocated), Err(Error::Limit)));
}

#[test]
fn control_rechecks_charge_the_original_arena_external_poles_and_workspace() {
    let body = graph(false, 1);
    let plan = CurvePlan::prepare(&body, CurveOptions::default(), None).unwrap();
    let poles = plan.curves()[0].poles();
    for bytes in [0, 64 * 1024 * 1024 + 1] {
        assert!(matches!(
            plan.certify_controls(0, poles, 16_000_000, bytes, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        plan.certify_controls(0, poles, 16_000_000, plan.retained_bytes(), None),
        Err(Error::Limit)
    ));
    let separate = poles.to_vec();
    assert!(matches!(
        plan.certify_controls(
            0,
            &separate,
            16_000_000,
            plan.retained_bytes() + WORKSPACE_BYTES,
            None
        ),
        Err(Error::Limit)
    ));
}
#[test]
fn frozen_original_curve_control_nets_survive_shared_surface_extraction() {
    for (left, starts) in [(false, 1), (true, 1), (false, 2)] {
        let body = graph(left, starts);
        let plan = CurvePlan::prepare(&body, CurveOptions::default(), None).unwrap();
        let mut hash = 0xcbf29ce484222325u64;
        for curve in plan.curves() {
            for number in curve.knots().iter().chain(curve.poles().iter().flatten()) {
                for byte in number.to_bits().to_le_bytes() {
                    hash ^= u64::from(byte);
                    hash = hash.wrapping_mul(0x100000001b3);
                }
            }
        }
        // Re-recorded under portable libm (bc7d6efc). Curve, knot and pole
        // counts are identical to the original macOS-libm capture.
        let expected = match (left, starts) {
            (false, 1) => 0x18c8c29a431750bd,
            (true, 1) => 0x759d8621cfa42cde,
            _ => 0x6acf3276ccc4c32c,
        };
        assert_eq!(hash, expected);
    }
}
