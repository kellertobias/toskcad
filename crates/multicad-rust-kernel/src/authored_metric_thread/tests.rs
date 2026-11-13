use super::*;
use crate::helical_thread::{Spec as LegacySpec, Thread};
use std::cell::Cell;
fn source(left: bool, nb: bool, large: bool) -> (Spec, Thread) {
    let (p, major, pd) = if large {
        (1., (5.794 + 5.974) / 4., (5.212 + 5.324) / 2.)
    } else {
        (0.5, (2.874 + 2.980) / 4., (2.580 + 2.655) / 2.)
    };
    let crest = p / 2. - (2. * major - pd) / 3_f64.sqrt();
    let r = p / (4. * 3_f64.sqrt());
    let spec = Spec {
        major_radius_mm: major,
        pitch_mm: p,
        length_mm: if nb { 4.3 } else { 4. },
        crest_width_mm: crest,
        starts: 1,
        left_hand: left,
    };
    let old = Thread::new(LegacySpec {
        major_radius_mm: major,
        pitch_mm: p,
        length_mm: spec.length_mm,
        depth_mm: (p - crest) / (2. / 3_f64.sqrt()) - r,
        crest_width_mm: crest,
        root_radius_mm: Some(r),
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
    })
    .unwrap();
    (spec, old)
}
fn options() -> TransferOptions {
    TransferOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    }
}
fn encloses(i: [f64; 2], x: f64) {
    assert!(i[0] <= x && x <= i[1], "{i:?} vs{x}");
}
#[test]
fn real_metric_engineering_dimensions_own_continuous_tangent_profile_and_bounded_legacy_transfer() {
    for large in [false, true] {
        for left in [false, true] {
            for nb in [false, true] {
                let (spec, legacy) = source(left, nb, large);
                let a = AuthoredMetricThread::new(spec).unwrap();
                assert_eq!(a.spec(), spec);
                assert_eq!(a.full_length_mm(), spec.length_mm);
                assert_eq!(a.maximum_radius_mm(), spec.major_radius_mm);
                let knots = [
                    Knot::Begin,
                    Knot::CrestEnd,
                    Knot::UpperRoot,
                    Knot::LowerRoot,
                    Knot::CrestStart,
                    Knot::End,
                ];
                for band in 0..5 {
                    assert_eq!(
                        a.profile_bounds(band, [0., 0.]).unwrap().radius_mm,
                        a.knot(knots[band]).radius_mm
                    );
                    assert_eq!(
                        a.profile_bounds(band, [1., 1.]).unwrap().radius_mm,
                        a.knot(knots[band + 1]).radius_mm
                    );
                    let all = a.profile_bounds(band, [0., 1.]).unwrap();
                    assert!(all.radius_mm[0] > 0.);
                }
                let upper = a.profile_bounds(2, [0., 0.]).unwrap();
                let lower = a.profile_bounds(2, [1., 1.]).unwrap();
                encloses(upper.slope, -3_f64.sqrt());
                encloses(lower.slope, 3_f64.sqrt());
                let root = a.profile_bounds(2, [0.5, 0.5]).unwrap();
                assert!(root.slope[0] <= 0. && root.slope[1] >= 0.);
                encloses(
                    root.radius_mm,
                    spec.major_radius_mm
                        - (spec.pitch_mm - spec.crest_width_mm) * 3_f64.sqrt() / 2.
                        + interval(a.root_radius_mm()).midpoint(),
                );
                let t = a.certify_legacy(&legacy, options(), None).unwrap();
                assert!(t.maximum_position_error_mm > 0. && t.maximum_position_error_mm < 1e-7);
                assert!(t.maximum_regular_normal_angle_error_rad < 1e-7);
                assert!(t.material_volume_difference_bound_mm3 > 0.);
                assert!(t.root_radius_error_mm > 0. && t.root_radius_error_mm < 1e-14);
                encloses(a.knot(Knot::UpperRoot).phase_mm, spec.pitch_mm * 0.375);
                encloses(a.knot(Knot::LowerRoot).phase_mm, spec.pitch_mm * 0.625);
                for band in 0..5 {
                    for j in 0..101 {
                        let u = j as f64 / 100.;
                        let new = a.profile_bounds(band, [u, u]).unwrap();
                        let q = new.phase_mm[0].midpoint(new.phase_mm[1]);
                        let r = new.radius_mm[0].midpoint(new.radius_mm[1]);
                        let radius = legacy.radius_at_phase(q);
                        assert!((radius - r).abs() < t.maximum_position_error_mm);
                    }
                }
                assert_eq!(a.bounds_mm()[2], [0., spec.length_mm]);
            }
        }
    }
}
#[test]
fn whole_boxes_follow_the_new_continuous_law_and_retain_boundary_uncertainty() {
    use Classification::*;
    for left in [false, true] {
        let (spec, _) = source(left, true, false);
        let a = AuthoredMetricThread::new(spec).unwrap();
        assert_eq!(
            a.classify_box([[0., 0.], [0., 0.], [1., 1.1]], None),
            Ok(Inside)
        );
        assert_eq!(
            a.classify_box([[0., 0.], [0., 0.], [0., 0.]], None),
            Ok(UnresolvedBoundary)
        );
        assert_eq!(
            a.classify_box([[4., 4.], [0., 0.], [1., 1.]], None),
            Ok(Outside)
        );
        assert_eq!(
            a.classify_box([[0., 0.], [0., 0.], [-2., -1.]], None),
            Ok(Outside)
        );
        for [x, y] in [[1_f64, 0.2], [-1., 0.2], [-1., -0.2], [0., 1.], [0., -1.]] {
            let theta = y.atan2(x);
            let norm = x.hypot(y);
            let z = 1.13;
            let rate = (if left { -1. } else { 1. }) * spec.pitch_mm / std::f64::consts::TAU;
            let q = (z - rate * theta).rem_euclid(spec.pitch_mm);
            let ks = a.knots();
            let band = (0..5).find(|&i| q <= ks[i + 1][0].midpoint()).unwrap();
            let u = (q - ks[band][0].midpoint())
                / (ks[band + 1][0].midpoint() - ks[band][0].midpoint());
            let r = a.profile_bounds(band, [u, u]).unwrap().radius_mm;
            let r = r[0].midpoint(r[1]);
            for (delta, expected) in [(-0.005, Inside), (0.005, Outside)] {
                let xx = (r + delta) * x / norm;
                let yy = (r + delta) * y / norm;
                assert_eq!(
                    a.classify_box([[xx, xx], [yy, yy], [z, z]], None),
                    Ok(expected)
                );
            }
        }
    }
}
#[test]
fn source_and_transfer_refuse_invalid_precision_ownership_resource_and_cancellation() {
    let (spec, legacy) = source(false, false, false);
    let a = AuthoredMetricThread::new(spec).unwrap();
    for x in [0., f64::NAN, f64::INFINITY] {
        assert_eq!(
            AuthoredMetricThread::new(Spec {
                crest_width_mm: x,
                ..spec
            }),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        AuthoredMetricThread::new(Spec { starts: 9, ..spec }),
        Err(Error::Limit)
    );
    assert_eq!(
        AuthoredMetricThread::new(Spec {
            crest_width_mm: spec.pitch_mm * 0.8,
            ..spec
        }),
        Err(Error::Unresolved)
    );
    assert_eq!(a.profile_bounds(5, [0., 1.]), Err(Error::Invalid));
    assert_eq!(a.profile_bounds(0, [-1., 1.]), Err(Error::Invalid));
    assert_eq!(
        a.classify_box([[f64::NAN, 0.]; 3], None),
        Err(Error::Invalid)
    );
    assert_eq!(
        a.certify_legacy(
            &legacy,
            TransferOptions {
                position_limit_mm: 0.,
                ..options()
            },
            None
        ),
        Err(Error::Invalid)
    );
    let other = source(true, false, false).1;
    assert_eq!(
        a.certify_legacy(&other, options(), None),
        Err(Error::Invalid)
    );
    assert_eq!(
        a.certify_legacy(
            &legacy,
            TransferOptions {
                position_limit_mm: 1e-30,
                ..options()
            },
            None
        ),
        Err(Error::Unresolved)
    );
    let count = Cell::new(0);
    let poll = || {
        count.set(count.get() + 1);
        false
    };
    a.certify_legacy(&legacy, options(), Some(&poll)).unwrap();
    for at in 1..=count.get() {
        let c = Cell::new(0);
        let cancel = || {
            c.set(c.get() + 1);
            c.get() == at
        };
        assert_eq!(
            a.certify_legacy(&legacy, options(), Some(&cancel)),
            Err(Error::Cancelled)
        );
    }
    assert_eq!(
        a.classify_box([[0., 0.], [0., 0.], [1., 1.]], Some(&|| true)),
        Err(Error::Cancelled)
    );
    assert_eq!(a.spec(), spec);
}

#[test]
fn exact_shared_algebraic_tangencies_are_not_fitted_from_rounded_supports() {
    use num_rational::BigRational as Q;
    let add = |a: [Q; 2], b: [Q; 2]| [&a[0] + &b[0], &a[1] + &b[1]];
    let mul = |a: [Q; 2], b: [Q; 2]| {
        [
            &a[0] * &b[0] + Q::from_integer(3.into()) * &a[1] * &b[1],
            &a[0] * &b[1] + &a[1] * &b[0],
        ]
    };
    let q = |x| Q::from_float(x).unwrap();
    let zero = q(0.);
    let half = q(0.5);
    for large in [false, true] {
        let (s, _) = source(false, false, large);
        let root = [zero.clone(), q(s.pitch_mm) / q(12.)];
        let center = add(
            [
                q(s.major_radius_mm),
                -(q(s.pitch_mm) - q(s.crest_width_mm)) * &half,
            ],
            mul([q(2.), zero.clone()], root.clone()),
        );
        let tangent = add(
            center.clone(),
            mul([-half.clone(), zero.clone()], root.clone()),
        );
        for sign in [-1_f64, 1.] {
            let qt = [
                q(s.pitch_mm) * &half + q(sign) * q(s.pitch_mm) / q(8.),
                zero.clone(),
            ];
            let base = if sign < 0. {
                [q(s.crest_width_mm) * &half, zero.clone()]
            } else {
                [q(s.pitch_mm) - q(s.crest_width_mm) * &half, zero.clone()]
            };
            let difference = add(qt.clone(), [-base[0].clone(), -base[1].clone()]);
            let line = add(
                [q(s.major_radius_mm), zero.clone()],
                mul([zero.clone(), q(sign)], difference),
            );
            assert_eq!(line, tangent);
            let x = [qt[0].clone() - q(s.pitch_mm) * &half, zero.clone()];
            let y = add(tangent.clone(), [-center[0].clone(), -center[1].clone()]);
            assert_eq!(
                add(mul(x.clone(), x), mul(y.clone(), y)),
                mul(root.clone(), root.clone())
            );
        }
    }
}
#[test]
fn actual_regular_normal_correspondence_is_bounded_without_claiming_old_closure_identity() {
    for left in [false, true] {
        let (s, old) = source(left, true, false);
        let a = AuthoredMetricThread::new(s).unwrap();
        let t = a.certify_legacy(&old, options(), None).unwrap();
        let ks = a.knots();
        let ob = old.profile_breaks();
        for band in 0..5 {
            for u in [0.1, 0.4, 0.9] {
                let n = a.profile_bounds(band, [u, u]).unwrap();
                let qn = n.phase_mm[0].midpoint(n.phase_mm[1]);
                let rn = n.radius_mm[0].midpoint(n.radius_mm[1]);
                let mn = n.slope[0].midpoint(n.slope[1]);
                let qo = ob[band] + u * (ob[band + 1] - ob[band]);
                let ro = old.radius_at_phase(qo);
                let mo = old.slope_at_phase(qo, band);
                let z = 1.23;
                let normal = |q: f64, r: f64, m: f64, rate: f64| {
                    let theta = (z - q) / rate;
                    let (s, c) = theta.sin_cos();
                    let beta = rate * m / r;
                    let v = [c - beta * s, s + beta * c, -m];
                    let l = v.into_iter().map(|x| x * x).sum::<f64>().sqrt();
                    v.map(|x| x / l)
                };
                let nn = normal(qn, rn, mn, a.rate().midpoint());
                let no = normal(qo, ro, mo, old.helical_rate());
                let d = nn
                    .into_iter()
                    .zip(no)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                assert!(d <= t.bands[band].regular_normal_angle_error_rad + 1e-14);
                assert!(ks[band][0].lo <= qn && qn <= ks[band + 1][0].hi);
            }
        }
    }
}

#[test]
fn finite_but_unresolvable_box_arithmetic_does_not_publish_material_membership() {
    let (s, _) = source(false, false, false);
    let a = AuthoredMetricThread::new(s).unwrap();
    assert_eq!(
        a.classify_box([[f64::MIN, f64::MAX], [f64::MIN, f64::MAX], [1., 1.]], None),
        Ok(Classification::UnresolvedBoundary)
    );
    assert_eq!(
        a.classify_box([[0., 0.], [0., 0.], [f64::MIN, f64::MAX]], None),
        Err(Error::Unresolved)
    );
}

#[test]
fn actual_material_radial_transfer_includes_rate_rounding_and_old_one_sided_rails() {
    for large in [false, true] {
        for left in [false, true] {
            let (s, old) = source(left, true, large);
            let a = AuthoredMetricThread::new(s).unwrap();
            let proof = a.certify_legacy(&old, options(), None).unwrap();
            for i in 0..257 {
                let theta = i as f64 * std::f64::consts::TAU / 256.;
                for z in [0., 0.23, 1.73, s.length_mm] {
                    let radius = a
                        .phase_radius(I::point(z) - a.rate() * I::point(theta))
                        .unwrap();
                    let old_radius = old.radius(theta, z);
                    assert!(
                        old_radius >= radius.lo - proof.material_radial_error_mm
                            && old_radius <= radius.hi + proof.material_radial_error_mm
                    );
                }
            }
        }
    }
}

#[test]
fn extreme_aspect_analytic_support_transfer_refuses_unproved_accuracy_and_bounds_finite_phase_span()
{
    let p = 0.001;
    let crest = 0.000125;
    let r = p / (4. * 3_f64.sqrt());
    let spec = Spec {
        major_radius_mm: 1e7,
        pitch_mm: p,
        length_mm: 0.12,
        crest_width_mm: crest,
        starts: 1,
        left_hand: false,
    };
    let a = AuthoredMetricThread::new(spec).unwrap();
    let old = Thread::new(LegacySpec {
        major_radius_mm: spec.major_radius_mm,
        pitch_mm: p,
        length_mm: spec.length_mm,
        crest_width_mm: crest,
        root_radius_mm: Some(r),
        depth_mm: (p - crest) / (2. / 3_f64.sqrt()) - r,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
    })
    .unwrap();
    assert_eq!(
        a.certify_legacy(&old, options(), None),
        Err(Error::Unresolved)
    );
    let proof = a
        .certify_legacy(
            &old,
            TransferOptions {
                position_limit_mm: 0.01,
                normal_angle_limit_rad: 1e-5,
            },
            None,
        )
        .unwrap();
    assert!(proof.maximum_position_error_mm > options().position_limit_mm);
    assert!(proof.maximum_position_error_mm < 0.01);
    assert!(proof.maximum_regular_normal_angle_error_rad < 1e-5);
    // Stress the finite-phase chart at ~750rad. These scalar diagnostics are
    // not a generic bound on either floating evaluator; the interval proof is.
    let breaks = old.profile_breaks();
    for band in 0..5 {
        let u = 0.37;
        let n = a.profile_bounds(band, [u, u]).unwrap();
        let qn = n.phase_mm[0].midpoint(n.phase_mm[1]);
        let rn = n.radius_mm[0].midpoint(n.radius_mm[1]);
        let qo = breaks[band] + u * (breaks[band + 1] - breaks[band]);
        let ro = old.radius_at_phase(qo);
        let tn = (spec.length_mm - qn) / a.rate().midpoint();
        let to = (spec.length_mm - qo) / old.helical_rate();
        let distance = (rn * tn.cos() - ro * to.cos()).hypot(rn * tn.sin() - ro * to.sin());
        assert!(distance < proof.bands[band].position_error_mm);
    }
}
#[test]
fn explicit_multistart_source_uses_the_same_profile_and_actual_retained_rate_correspondence() {
    let (s, old) = source(true, true, false);
    for starts in [2, 8] {
        let a = AuthoredMetricThread::new(Spec { starts, ..s }).unwrap();
        let old = Thread::new(LegacySpec {
            starts,
            ..old.spec()
        })
        .unwrap();
        let bounds = a.certify_legacy(&old, options(), None).unwrap();
        assert!(bounds.maximum_position_error_mm < 1e-7);
        assert_eq!(
            a.classify_box([[0., 0.], [0., 0.], [1., 1.]], None),
            Ok(Classification::Inside)
        );
    }
}

#[test]
fn detached_shell_error_boundary_preserves_precision_resource_and_cancellation_types() {
    use crate::authored_internal_treatment::Error as MeshError;
    for (input, expected) in [
        (MeshError::Invalid, Error::Invalid),
        (MeshError::Unresolved, Error::Unresolved),
        (MeshError::Limit, Error::Limit),
        (MeshError::Cancelled, Error::Cancelled),
    ] {
        assert_eq!(shell_error(input), expected);
    }
}
