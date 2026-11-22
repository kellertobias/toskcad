use super::*;
use std::cell::Cell;
fn npt(k: f64, quarter: bool, clipped: bool) -> Spec {
    let (tpi, e1, l1, l2, od, truncation) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.0605)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.0645)
    };
    let p = 25.4 / tpi;
    let m = 3_f64.sqrt();
    let h = p / (1. / (m - k) + 1. / (m + k));
    let depth = h - 2. * truncation * p;
    let width = (p - depth * (1. / (m - k) + 1. / (m + k))) / 2.;
    let base = e1 * 25.4 - if k > 0. { l1 * 25.4 / 16. } else { 0. };
    Spec {
        major_radius_mm: (base + depth) / 2.,
        pitch_mm: p,
        length_mm: if k > 0. {
            l2 * 25.4
        } else {
            l1 * 25.4 + 4. * p
        },
        radial_taper: k,
        section: Section::NptFlat {
            depth_mm: depth,
            flat_width_mm: width,
        },
        crest_radius_limit_mm: clipped.then_some(od * 25.4 / 2.),
    }
}
fn bspt(k: f64, quarter: bool) -> Spec {
    let (tpi, major, gauge, length) = if quarter {
        (19., 13.157, 6., 9.7)
    } else {
        (28., 9.728, 4., 6.5)
    };
    let p = 25.4 / tpi;
    let a = 27.5_f64.to_radians();
    let m = 1. / a.portable_tan();
    let depth = p / (3. * a.portable_tan());
    let h = p * (m * m - k * k) / (2. * m);
    let rounding = (h - depth) / (2. * (1. / a.portable_sin() - (1. + k * k).sqrt()));
    Spec {
        major_radius_mm: major / 2. - k * if k > 0. { gauge } else { p / 2. },
        pitch_mm: p,
        length_mm: if k > 0. {
            length
        } else {
            if quarter { 11. } else { 7.4 }
        },
        radial_taper: k,
        section: Section::BsptRounded {
            depth_mm: depth,
            rounding_radius_mm: rounding,
        },
        crest_radius_limit_mm: None,
    }
}
fn source(spec: Spec) -> Thread {
    Thread::new(spec, 1e-7, 1e-6).unwrap()
}
fn contains(interval: [f64; 2], x: f64) {
    assert!(
        interval[0] <= x && x <= interval[1],
        "{interval:?} excludes {x}"
    );
}

#[test]
fn all_core_pipe_rows_retain_real_axial_slopes_circles_and_closed_junctions() {
    for quarter in [false, true] {
        for k in [-1. / 32., 1. / 32.] {
            for spec in [npt(k, quarter, false), bspt(k, quarter)] {
                let t = source(spec);
                assert_eq!(t.spec(), spec);
                assert_eq!(t.tolerance_mm(), 1e-7);
                assert_eq!(t.minimum_feature_mm(), 1e-6);
                assert!(t.coefficient_error_mm() > 0.);
                for j in t.junctions() {
                    assert!(
                        j.radial_bounds_mm[1] - j.radial_bounds_mm[0]
                            < 4. * t.coefficient_error_mm()
                    );
                }
                let z = spec.length_mm * 0.7;
                for (band, b) in t.bands().iter().enumerate() {
                    let q = (b.phase_mm[0] + b.phase_mm[1]) / 2.;
                    let bound = t.radial_bounds([q; 2], [z; 2]).unwrap();
                    let expected = match b.curve {
                        Curve::Line {
                            phase_mm,
                            radius_mm,
                            slope,
                        } => radius_mm + slope * (q - phase_mm),
                        Curve::Circle {
                            phase_mm,
                            radius_mm,
                            circle_radius_mm: r,
                            sign,
                        } => radius_mm + sign * (r * r - (q - phase_mm).powi(2)).sqrt(),
                    } + k * (z - q);
                    contains(bound, expected);
                    assert_eq!(
                        t.locate_radial(q, z, expected - 0.001).unwrap(),
                        Location::Inside
                    );
                    assert_eq!(
                        t.locate_radial(q, z, expected + 0.001).unwrap(),
                        Location::Outside
                    );
                    let point = t.surface_point_bounds(band, [q; 2], [z; 2]).unwrap();
                    let theta = std::f64::consts::TAU * (z - q) / spec.pitch_mm;
                    contains(point[0], expected * theta.portable_cos());
                    contains(point[1], expected * theta.portable_sin());
                    contains(point[2], z);
                    // At fixed theta, q and z advance together. Actual axial
                    // sections keep unsheared lines and literal circular radii.
                    let dq = (b.phase_mm[1] - b.phase_mm[0]) * 0.1;
                    for sign in [-1., 1.] {
                        let q1 = q + sign * dq;
                        let z1 = z + sign * dq;
                        let rr = t.radial_bounds([q1; 2], [z1; 2]).unwrap();
                        match b.curve {
                            Curve::Line { slope, .. } => contains(rr, expected + slope * sign * dq),
                            Curve::Circle {
                                phase_mm,
                                radius_mm,
                                circle_radius_mm: r,
                                sign: radial_sign,
                            } => {
                                let center_r = radius_mm + k * (z - q);
                                let radius = center_r
                                    + radial_sign * (r * r - (q1 - phase_mm).powi(2)).sqrt();
                                contains(rr, radius);
                            }
                        }
                    }
                }
                assert_eq!(t.locate_radial(0., -1., 0.).unwrap(), Location::Outside);
                assert_eq!(t.locate_radial(0., 0., 0.).unwrap(), Location::Boundary);
                assert_eq!(
                    t.locate_radial(0., spec.length_mm, 0.).unwrap(),
                    Location::Boundary
                );
            }
        }
    }
}

#[test]
fn npt_stock_ceiling_keeps_imperfect_crests_and_true_positive_material() {
    for quarter in [false, true] {
        let spec = npt(1. / 32., quarter, true);
        let t = source(spec);
        let limit = spec.crest_radius_limit_mm.unwrap();
        let z = (limit - spec.major_radius_mm) / spec.radial_taper;
        assert!(0. < z && z < spec.length_mm);
        for zz in [z - spec.pitch_mm / 2., (z + spec.length_mm) / 2.] {
            let expected = (spec.major_radius_mm + spec.radial_taper * zz).min(limit);
            contains(t.radial_bounds([0.; 2], [zz; 2]).unwrap(), expected);
            assert_eq!(
                t.locate_radial(0., zz, expected - 1e-5).unwrap(),
                Location::Inside
            );
            assert_eq!(
                t.locate_radial(0., zz, expected + 1e-5).unwrap(),
                Location::Outside
            );
        }
        let full = source(Spec {
            crest_radius_limit_mm: None,
            ..spec
        });
        let clipped = t.volume_bounds_mm3(1e-6, None).unwrap();
        let original = full.volume_bounds_mm3(1e-6, None).unwrap();
        assert!(clipped[1] < original[0]);
    }
}

#[test]
fn mass_and_support_bound_independent_section_integration_and_actual_points() {
    for quarter in [false, true] {
        for k in [-1. / 32., 1. / 32.] {
            for spec in [
                npt(k, quarter, false),
                npt(k, quarter, true),
                bspt(k, quarter),
            ] {
                let t = source(spec);
                let volume = t.volume_bounds_mm3(1e-6, None).unwrap();
                assert!(volume[1] - volume[0] <= 1e-6);
                // Composite midpoint cubature is an independent approximation,
                // used with its convergence difference rather than as authority.
                let approximate = |n: usize| {
                    let mut sum = 0.;
                    for b in t.bands() {
                        for i in 0..n {
                            let q = b.phase_mm[0]
                                + (b.phase_mm[1] - b.phase_mm[0]) * (i as f64 + 0.5) / n as f64;
                            for j in 0..n {
                                let z = spec.length_mm * (j as f64 + 0.5) / n as f64;
                                let r = b.curve.radius(I::point(q))
                                    + I::point(k) * (I::point(z) - I::point(q));
                                sum += ceiling(r, spec.crest_radius_limit_mm).midpoint().powi(2)
                                    * (b.phase_mm[1] - b.phase_mm[0])
                                    * spec.length_mm
                                    / (n * n) as f64;
                            }
                        }
                    }
                    PI / spec.pitch_mm * sum
                };
                let coarse = approximate(32);
                let fine = approximate(64);
                assert!((volume[0] - fine).abs() <= 2. * (fine - coarse).abs() + 1e-5);
                contains(
                    t.support_bounds([0., 0., 1.], 1e-6, None).unwrap(),
                    spec.length_mm,
                );
                contains(t.support_bounds([0., 0., -1.], 1e-6, None).unwrap(), 0.);
                assert_eq!(t.support_bounds([0.; 3], 1e-6, None).unwrap(), [0.; 2]);
                let bound = t.support_bounds([1., 0., 0.3], 0.01, None).unwrap();
                assert!(bound[1] - bound[0] <= 0.01);
                for b in t.bands() {
                    for q in b.phase_mm {
                        for z in [0., spec.length_mm] {
                            let r = ceiling(
                                b.curve.radius(I::point(q))
                                    + I::point(k) * (I::point(z) - I::point(q)),
                                spec.crest_radius_limit_mm,
                            )
                            .midpoint();
                            let theta = std::f64::consts::TAU * (z - q) / spec.pitch_mm;
                            assert!(r * theta.portable_cos() + 0.3 * z <= bound[1]);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn invalid_tight_unsupported_and_cancelled_requests_preserve_retries() {
    let s = bspt(1. / 32., false);
    let t = source(s);
    let before = t;
    for accuracy in [0., f64::NAN, f64::INFINITY] {
        assert_eq!(t.volume_bounds_mm3(accuracy, None), Err(Error::Invalid));
        assert_eq!(
            t.support_bounds([1., 0., 0.], accuracy, None),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        t.support_bounds([f64::NAN, 0., 0.], 1e-3, None),
        Err(Error::Invalid)
    );
    assert_eq!(
        t.volume_bounds_mm3(1e-6, Some(&|| true)),
        Err(Error::Cancelled)
    );
    assert_eq!(
        t.support_bounds([1.; 3], 1e-3, Some(&|| true)),
        Err(Error::Cancelled)
    );
    let calls = Cell::new(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() > 3
    };
    assert_eq!(
        t.volume_bounds_mm3(1e-12, Some(&cancel)),
        Err(Error::Cancelled)
    );
    assert!(t.volume_bounds_mm3(1e-18, None).is_err());
    assert!(t.support_bounds([1.; 3], 1e-14, None).is_err());
    for bad in [
        Spec { pitch_mm: 0., ..s },
        Spec {
            radial_taper: 0.,
            ..s
        },
        Spec {
            length_mm: 1e8,
            ..s
        },
        Spec {
            length_mm: 0.1,
            ..s
        },
        Spec {
            crest_radius_limit_mm: Some(1e8),
            ..s
        },
        Spec {
            crest_radius_limit_mm: Some(f64::NAN),
            ..s
        },
        Spec {
            crest_radius_limit_mm: Some(10.),
            ..s
        },
        Spec {
            section: Section::BsptRounded {
                depth_mm: 0.1,
                rounding_radius_mm: 0.1,
            },
            ..s
        },
    ] {
        assert!(Thread::new(bad, 1e-7, 1e-6).is_err());
    }
    assert!(Thread::new(s, 1e-20, 1e-6).is_err());
    assert!(Thread::new(s, 1e-7, 10.).is_err());
    for (q, z) in [
        ([-1., 0.], [0., 1.]),
        ([0., 2.], [0., 1.]),
        ([0., 0.], [0., 10.]),
        ([1., 0.], [0., 1.]),
        ([f64::NAN; 2], [0., 1.]),
    ] {
        assert!(t.radial_bounds(q, z).is_err());
    }
    for (q, z, r) in [
        (0., 0., -1.),
        (0., 0., f64::NAN),
        (f64::NAN, 0., 1.),
        (0., f64::NAN, 1.),
    ] {
        assert!(t.locate_radial(q, z, r).is_err());
    }
    assert!(t.surface_point_bounds(5, [0.; 2], [0.; 2]).is_err());
    assert!(t.surface_point_bounds(0, [-1.; 2], [0.; 2]).is_err());
    assert_eq!(
        t.support_bounds([f64::MAX; 3], 1e-3, None),
        Err(Error::Unresolved)
    );
    assert_eq!(
        t.support_bounds([0., 0., f64::MAX], 1e-3, None),
        Err(Error::Unresolved)
    );
    assert_eq!(
        t.support_bounds([0., 0., 1.], 1e-30, None),
        Err(Error::Unresolved)
    );
    assert_eq!(t, before);
    assert!(t.volume_bounds_mm3(1e-6, None).is_ok());
}

#[test]
fn independent_gauge_cones_and_axial_flat_mass_match_core_conventions() {
    for quarter in [false, true] {
        for k in [-1. / 32., 1. / 32.] {
            let spec = bspt(k, quarter);
            let t = source(spec);
            let Section::BsptRounded {
                depth_mm: depth,
                rounding_radius_mm: r,
            } = spec.section
            else {
                panic!()
            };
            let m = 1. / 27.5_f64.to_radians().portable_tan();
            let normal = (1. + k * k).sqrt();
            let h = spec.pitch_mm * (m * m - k * k) / (2. * m);
            let z = if k > 0. {
                if quarter { 6. } else { 4. }
            } else {
                spec.pitch_mm / 2.
            };
            // These are the cone/circle tangencies from the independent axial
            // construction, not the crest/root circle centers.
            let crest = (-k * r / normal).rem_euclid(spec.pitch_mm);
            let root = h / (m + k) + k * r / normal;
            contains(
                t.radial_bounds([crest; 2], [z; 2]).unwrap(),
                spec.major_radius_mm + k * z,
            );
            contains(
                t.radial_bounds([root; 2], [z; 2]).unwrap(),
                spec.major_radius_mm - depth + k * z,
            );
            for clipped in [false, true] {
                let s = npt(k, quarter, clipped);
                let t = source(s);
                let Section::NptFlat {
                    depth_mm: depth,
                    flat_width_mm: width,
                } = s.section
                else {
                    panic!()
                };
                let slopes = 1. / (3_f64.sqrt() - k) + 1. / (3_f64.sqrt() + k);
                let root_width = s.pitch_mm - width - depth * slopes;
                let area = |z: f64| {
                    let major = s.major_radius_mm + k * z;
                    let root = major - depth;
                    let ceiling = s.crest_radius_limit_mm.map_or(major, |r| r.min(major));
                    PI / s.pitch_mm
                        * (width * ceiling.powi(2)
                            + root_width * root.powi(2)
                            + slopes
                                * ((ceiling.powi(3) - root.powi(3)) / 3.
                                    + (major - ceiling) * ceiling.powi(2)))
                };
                let mut stations = vec![0., s.length_mm];
                if let Some(limit) = s.crest_radius_limit_mm {
                    let crossing = (limit - s.major_radius_mm) / k;
                    if crossing > 0. && crossing < s.length_mm {
                        stations.push(crossing);
                    }
                }
                stations.sort_by(f64::total_cmp);
                let expected: f64 = stations
                    .windows(2)
                    .map(|ab| {
                        let center = (ab[0] + ab[1]) / 2.;
                        let half = (ab[1] - ab[0]) / 2.;
                        half * (area(center - half / 3_f64.sqrt())
                            + area(center + half / 3_f64.sqrt()))
                    })
                    .sum();
                contains(t.volume_bounds_mm3(1e-6, None).unwrap(), expected);
            }
        }
    }
}

#[test]
fn exact_ceiling_contact_and_uncertain_radial_queries_are_distinct() {
    let spec = npt(1. / 32., false, true);
    let t = source(spec);
    let limit = spec.crest_radius_limit_mm.unwrap();
    assert_eq!(
        t.locate_radial(0., spec.length_mm * 0.9, limit).unwrap(),
        Location::Boundary
    );
    let t = source(bspt(-1. / 32., false));
    let r = t.radial_bounds([0.; 2], [1.; 2]).unwrap();
    assert_eq!(
        t.locate_radial(0., 1., (r[0] + r[1]) / 2.).unwrap(),
        Location::Unresolved
    );
    assert_eq!(
        t.radial_bounds([0.; 2], [1.; 2]).unwrap(),
        t.radial_bounds([t.spec().pitch_mm; 2], [1.; 2]).unwrap()
    );
    assert!(t.surface_point_bounds(0, [0.; 2], [-1.; 2]).is_err());
    for spec in [
        Spec {
            crest_radius_limit_mm: Some(1.),
            ..npt(1. / 32., false, false)
        },
        Spec {
            section: Section::NptFlat {
                depth_mm: 0.5,
                flat_width_mm: 2.,
            },
            ..npt(1. / 32., false, false)
        },
    ] {
        assert_eq!(Thread::new(spec, 1e-7, 1e-6), Err(Error::Unresolved));
    }
}

#[test]
fn literal_root_circle_changes_keep_tiny_positive_material_and_junction_closures() {
    let s = bspt(1. / 32., false);
    let original = source(s);
    let Section::BsptRounded {
        depth_mm,
        rounding_radius_mm,
    } = s.section
    else {
        panic!()
    };
    let delta = original.coefficient_error_mm() / 2.;
    let changed = source(Spec {
        section: Section::BsptRounded {
            depth_mm: depth_mm - delta,
            rounding_radius_mm,
        },
        ..s
    });
    let root = changed.bands()[2];
    let q = (root.phase_mm[0] + root.phase_mm[1]) / 2.;
    let z = 2.;
    let a = original.radial_bounds([q; 2], [z; 2]).unwrap();
    let b = changed.radial_bounds([q; 2], [z; 2]).unwrap();
    assert!(a[1] < b[0]);
    let radius = (a[1] + b[0]) / 2.;
    assert_eq!(
        original.locate_radial(q, z, radius).unwrap(),
        Location::Outside
    );
    assert_eq!(
        changed.locate_radial(q, z, radius).unwrap(),
        Location::Inside
    );
    assert_eq!(
        changed.spec().section,
        Section::BsptRounded {
            depth_mm: depth_mm - delta,
            rounding_radius_mm
        }
    );
    let closure = changed.junctions()[1];
    assert!(closure.radial_bounds_mm[1] - closure.radial_bounds_mm[0] > delta / 2.);
    let Curve::Circle {
        circle_radius_mm, ..
    } = root.curve
    else {
        panic!("retained root circle")
    };
    assert_eq!(circle_radius_mm, rounding_radius_mm);
}

#[test]
fn whole_sections_all_clipped_crests_and_intrinsic_junction_limit_remain_honest() {
    let s = npt(1. / 32., false, false);
    let t = source(Spec {
        crest_radius_limit_mm: Some(s.major_radius_mm - 0.05),
        ..s
    });
    let clipped = t.volume_bounds_mm3(1e-6, None).unwrap();
    assert!(clipped[1] < source(s).volume_bounds_mm3(1e-6, None).unwrap()[0]);
    let all = t
        .radial_bounds([0., s.pitch_mm], [0., s.length_mm])
        .unwrap();
    for b in t.bands() {
        for q in b.phase_mm {
            for z in [0., s.length_mm] {
                let r = t.radial_bounds([q; 2], [z; 2]).unwrap();
                assert!(all[0] <= r[0] && r[1] <= all[1]);
            }
        }
    }
    let base = bspt(1. / 32., false);
    let original = source(base);
    let Section::BsptRounded {
        depth_mm: _,
        rounding_radius_mm: r,
    } = base.section
    else {
        panic!()
    };
    let angle = 27.5_f64.to_radians();
    let m = 1. / angle.portable_tan();
    let k = base.radial_taper;
    let expected = base.pitch_mm * (m * m - k * k) / (2. * m)
        - 2. * r * (1. / angle.portable_sin() - (1. + k * k).sqrt());
    let mut refused = false;
    for scale in [0.99, 0.995, 0.999, 0.9999, 1.] {
        let spec = Spec {
            section: Section::BsptRounded {
                depth_mm: expected - original.coefficient_error_mm() * scale,
                rounding_radius_mm: r,
            },
            ..base
        };
        refused |= Thread::new(spec, 1e-7, 1e-6) == Err(Error::Unresolved);
    }
    assert!(
        refused,
        "Actual junction enclosure must fit its declared intrinsic budget"
    );
    assert_eq!(
        original.volume_bounds_mm3(f64::MAX, None),
        Err(Error::Unresolved)
    );
}

#[test]
fn clipping_roots_at_band_and_end_planes_keep_the_entire_uncertain_strip() {
    for k in [-1. / 32., 1. / 32.] {
        let s = npt(k, false, false);
        let limit = s.major_radius_mm + (k * s.length_mm).min(0.);
        let t = source(Spec {
            crest_radius_limit_mm: Some(limit),
            ..s
        });
        let bounds = t.volume_bounds_mm3(1e-6, None).unwrap();
        assert!(bounds[0] > 0. && bounds[1] - bounds[0] <= 1e-6);
        assert!(bounds[1] < source(s).volume_bounds_mm3(1e-6, None).unwrap()[0]);
        for z in [0., s.length_mm] {
            let r = t.radial_bounds([0.; 2], [z; 2]).unwrap();
            assert!(r[1] <= limit);
        }
    }
}

#[test]
fn cancellation_through_every_clipped_piece_leaves_the_literal_source_reusable() {
    let t = source(npt(1. / 32., false, true));
    let before = t;
    let original = t.volume_bounds_mm3(1e-6, None).unwrap();
    let mut completed = false;
    for stop in 0..64 {
        let calls = Cell::new(0);
        let cancel = || {
            let n = calls.get();
            calls.set(n + 1);
            n >= stop
        };
        match t.volume_bounds_mm3(1e-6, Some(&cancel)) {
            Err(error) => assert_eq!(error, Error::Cancelled),
            Ok(result) => {
                assert_eq!(result, original);
                completed = true;
            }
        }
        assert_eq!(t, before);
    }
    assert!(completed);
    assert_eq!(t.volume_bounds_mm3(1e-6, None).unwrap(), original);
}
