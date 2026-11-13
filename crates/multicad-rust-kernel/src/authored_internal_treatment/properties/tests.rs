use super::*;
use crate::{
    Circle, Rational, authored_flat_thread::SourceOptions, cylinder::AxialCylinder,
    helical_thread::Spec,
};
fn thread(left: bool, starts: u32) -> AuthoredThread {
    AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 4.125,
            depth_mm: 0.3,
            crest_width_mm: 0.1,
            lower_flank_deg: 18.,
            upper_flank_deg: 12.,
            starts,
            left_hand: left,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap()
}
fn body(left: bool, starts: u32, entry: bool, tail: bool) -> Cavity {
    Cavity::new(
        thread(left, starts),
        Intent {
            entry_depth_mm: entry.then_some(0.7),
            terminal: tail.then_some(Terminal {
                width_mm: 0.35,
                radial_depth_mm: 0.15,
            }),
        },
    )
    .unwrap()
}
fn material(body: Cavity) -> Material {
    let q = |n, d| Rational::new(n, d).unwrap();
    let stock = AxialCylinder::new(
        7,
        Circle::new([q(3, 7), q(-5, 11)], q(5, 2)).unwrap(),
        q(7, 13),
        q(33, 8),
    )
    .unwrap();
    Material::new(stock, body, 0.2).unwrap()
}
fn mid(v: Enclosure) -> f64 {
    (v.lower + v.upper) * 0.5
}
fn near(v: Enclosure, expected: f64) {
    assert!((mid(v) - expected).abs() < 2e-8, "{v:?} vs{expected}");
    assert!(v.lower <= v.upper);
}
#[test]
fn source_roles_ordinals_and_periodic_ownership_are_explicit() {
    for left in [false, true] {
        for starts in [1, 2, 3, 8] {
            for flags in [(true, false), (false, true), (true, true)] {
                let body = body(left, starts, flags.0, flags.1);
                let material = material(body);
                for i in 0..material.face_count() {
                    let key = material.face_key(i).unwrap();
                    assert_eq!(material.face_index(key), Some(i));
                }
                assert_eq!(material.face_key(material.face_count()), None);
                assert_eq!(body.face_index(FaceKey::MaterialOuterWall), None);
                for p in -10..10 {
                    for band in 0..5 {
                        assert_eq!(
                            body.thread_face_key(p, band).unwrap(),
                            body.thread_face_key(p + i64::from(starts), band).unwrap()
                        );
                    }
                }
                assert_eq!(body.thread_face_key(0, 5), Err(Error::Invalid));
            }
        }
    }
}
#[test]
fn every_positive_entry_local_face_proves_area_moments_and_centroid() {
    for left in [false, true] {
        for starts in [1, 2] {
            let body = body(left, starts, true, true);
            let original = body;
            for index in 0..body.face_count() {
                let key = body.face_key(index).unwrap();
                let p = body
                    .local_face_properties(key, FacePropertyOptions::default())
                    .unwrap();
                assert_eq!(p.face, key);
                assert!(p.area_mm2.lower > 0.);
                for c in p.centroid_mm {
                    assert!(
                        c.lower >= -body.maximum_radius_mm() - 1e-6
                            && c.upper <= body.full_length_mm() + 1e-6
                    );
                }
            }
            let material = material(body);
            for index in 0..material.face_count() {
                let p = material
                    .local_face_properties(
                        material.face_key(index).unwrap(),
                        FacePropertyOptions::default(),
                    )
                    .unwrap();
                assert!(p.area_mm2.lower > 0.);
            }
            assert_eq!(body, original);
        }
    }
}
#[test]
fn circles_cylinders_and_flat_root_have_independent_closed_forms() {
    let body = body(false, 1, true, true);
    let opts = FacePropertyOptions::default();
    let material = material(body);
    let pi = std::f64::consts::PI;
    let a = body.entry_radius_mm().unwrap();
    let cap = body.local_face_properties(FaceKey::StartCap, opts).unwrap();
    near(cap.area_mm2, pi * a * a);
    near(cap.centroid_mm[0], 0.);
    near(cap.centroid_mm[2], 0.);
    let outer = material
        .local_face_properties(FaceKey::MaterialOuterWall, opts)
        .unwrap();
    near(outer.area_mm2, 2. * pi * 2.5 * 4.125);
    near(outer.centroid_mm[2], 4.125 / 2.);
    let annulus = material
        .local_face_properties(FaceKey::StartCap, opts)
        .unwrap();
    near(annulus.area_mm2, pi * (2.5 - a) * (2.5 + a));
    let tail = body
        .local_face_properties(FaceKey::TerminalCylinder, opts)
        .unwrap();
    let z = body.terminal_station_mm().unwrap();
    let r = body.terminal_radius_mm().unwrap();
    near(tail.area_mm2, 2. * pi * r * (4.125 - z));
    near(tail.centroid_mm[2], (4.125 + z) / 2.);
    let root = body
        .local_face_properties(
            FaceKey::Thread {
                start: 0,
                band: ThreadBandRole::Root,
            },
            opts,
        )
        .unwrap();
    let band = body.thread.bands()[2];
    near(
        root.area_mm2,
        (z - 0.7) * 2. * pi / 0.5 * band.radius_mm[0] * (band.phase_mm[1] - band.phase_mm[0]),
    );
    near(root.centroid_mm[2], (z + 0.7) / 2.);
}
#[test]
fn profile_caps_keep_their_non_axis_centroids_and_actual_end_rotation() {
    let start = body(false, 1, false, true)
        .local_face_properties(FaceKey::StartCap, FacePropertyOptions::default())
        .unwrap();
    let endbody = body(false, 1, true, false);
    let end = endbody
        .local_face_properties(FaceKey::EndCap, FacePropertyOptions::default())
        .unwrap();
    assert!(mid(start.centroid_mm[0]).hypot(mid(start.centroid_mm[1])) > 1e-4);
    near(end.area_mm2, mid(start.area_mm2));
    near(end.first_moment_mm3[0], -mid(start.first_moment_mm3[1]));
    near(end.first_moment_mm3[1], mid(start.first_moment_mm3[0]));
    near(end.centroid_mm[2], 4.125);
    let cap = material(endbody)
        .local_face_properties(FaceKey::EndCap, FacePropertyOptions::default())
        .unwrap();
    near(cap.first_moment_mm3[0], -mid(end.first_moment_mm3[0]));
    near(cap.first_moment_mm3[1], -mid(end.first_moment_mm3[1]));
}

#[test]
fn flank_and_clipped_cone_areas_match_independent_antiderivatives() {
    for starts in [1, 2] {
        let body = body(false, starts, true, true);
        let bands = body.thread.bands();
        let a = body.entry.unwrap().radius;
        let k = (a - body.thread.minimum_radius_mm()) / body.entry.unwrap().depth;
        let s = 0.5 * f64::from(starts) / std::f64::consts::TAU;
        let f = body.tail.unwrap().station;
        for (index, role) in [(1, ThreadBandRole::Upper), (3, ThreadBandRole::Lower)] {
            let b = bands[index];
            let m = (b.radius_mm[1] - b.radius_mm[0]) / (b.phase_mm[1] - b.phase_mm[0]);
            let d = (m * s).abs() / (1. + m * m).sqrt();
            let p0 = |r: f64| 0.5 * (r * (r * r + d * d).sqrt() + d * d * (r / d).asinh());
            let p1 = |r: f64| (r * r + d * d).powf(1.5) / 3.;
            let [r0, r1] = b.radius_mm;
            let expected = (1. + m * m).sqrt() / s / m
                * ((f - a / k) * (p0(r1) - p0(r0)) + (p1(r1) - p1(r0)) / k);
            near(
                body.local_face_properties(
                    FaceKey::Thread {
                        start: 0,
                        band: role,
                    },
                    FacePropertyOptions::default(),
                )
                .unwrap()
                .area_mm2,
                expected,
            );
        }
        let radial = bands
            .into_iter()
            .map(|b| {
                let [r0, r1] = b.radius_mm;
                (b.phase_mm[1] - b.phase_mm[0]) * (a * a - (r0 * r0 + r0 * r1 + r1 * r1) / 3.)
            })
            .sum::<f64>()
            * f64::from(starts);
        near(
            body.local_face_properties(FaceKey::EntryCone, FacePropertyOptions::default())
                .unwrap()
                .area_mm2,
            (1. + k * k).sqrt() / (2. * s * k) * radial,
        );
    }
}
#[test]
fn handedness_and_each_actual_start_transform_first_moments() {
    let options = FacePropertyOptions::default();
    for starts in [1, 2, 3] {
        let right = body(false, starts, true, true);
        let left = body(true, starts, true, true);
        for role in [
            ThreadBandRole::Lower,
            ThreadBandRole::Crest,
            ThreadBandRole::Upper,
            ThreadBandRole::Root,
        ] {
            for start in 0..starts {
                let r = right
                    .local_face_properties(FaceKey::Thread { start, band: role }, options)
                    .unwrap();
                let l = left
                    .local_face_properties(
                        FaceKey::Thread {
                            start: (starts - start) % starts,
                            band: role,
                        },
                        options,
                    )
                    .unwrap();
                near(l.area_mm2, mid(r.area_mm2));
                near(l.first_moment_mm3[0], mid(r.first_moment_mm3[0]));
                near(l.first_moment_mm3[1], -mid(r.first_moment_mm3[1]));
                near(l.first_moment_mm3[2], mid(r.first_moment_mm3[2]));
                let next = right
                    .local_face_properties(
                        FaceKey::Thread {
                            start: (start + 1) % starts,
                            band: role,
                        },
                        options,
                    )
                    .unwrap();
                let a = std::f64::consts::TAU / f64::from(starts);
                let (s, c) = a.sin_cos();
                near(
                    next.first_moment_mm3[0],
                    c * mid(r.first_moment_mm3[0]) - s * mid(r.first_moment_mm3[1]),
                );
                near(
                    next.first_moment_mm3[1],
                    s * mid(r.first_moment_mm3[0]) + c * mid(r.first_moment_mm3[1]),
                );
            }
        }
    }
}
#[test]
fn invalid_precision_resource_absent_faces_and_cancel_are_atomic() {
    let body = body(false, 1, true, true);
    let original = body;
    let material = material(body);
    let face = FaceKey::Thread {
        start: 0,
        band: ThreadBandRole::Upper,
    };
    let options = FacePropertyOptions::default();
    for options in [
        FacePropertyOptions {
            area_tolerance_mm2: f64::NAN,
            ..options
        },
        FacePropertyOptions {
            first_moment_tolerance_mm3: 0.,
            ..options
        },
        FacePropertyOptions {
            centroid_tolerance_mm: -1.,
            ..options
        },
        FacePropertyOptions {
            maximum_cells: 0,
            ..options
        },
        FacePropertyOptions {
            maximum_cells: 65537,
            ..options
        },
    ] {
        assert_eq!(
            body.local_face_properties(face, options),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        body.local_face_properties(
            FaceKey::Thread {
                start: 1,
                band: ThreadBandRole::Root
            },
            options
        ),
        Err(Error::Invalid)
    );
    assert_eq!(
        body.local_face_properties(FaceKey::MaterialOuterWall, options),
        Err(Error::Invalid)
    );
    assert_eq!(
        body.local_face_properties(
            face,
            FacePropertyOptions {
                maximum_cells: 1,
                ..options
            }
        ),
        Err(Error::Limit)
    );
    assert_eq!(
        body.local_face_properties(
            FaceKey::StartCap,
            FacePropertyOptions {
                area_tolerance_mm2: 1e-30,
                ..options
            }
        ),
        Err(Error::Unresolved)
    );
    let small = Cavity::new(
        thread(false, 1),
        Intent {
            entry_depth_mm: Some(0.01),
            terminal: None,
        },
    )
    .unwrap();
    assert_eq!(
        small.local_face_properties(face, options),
        Err(Error::Unresolved)
    );
    let count = std::cell::Cell::new(0);
    let probe = || {
        count.set(count.get() + 1);
        false
    };
    let normal = body
        .local_face_properties_with_interrupt(face, options, Some(&probe))
        .unwrap();
    let total = count.get();
    for at in [1, total / 2, total] {
        count.set(0);
        let interrupt = || {
            count.set(count.get() + 1);
            count.get() >= at
        };
        assert_eq!(
            body.local_face_properties_with_interrupt(face, options, Some(&interrupt)),
            Err(Error::Cancelled)
        );
        assert_eq!(count.get(), at);
    }
    assert_eq!(body.local_face_properties(face, options).unwrap(), normal);
    assert_eq!(
        material.local_face_properties_with_interrupt(
            FaceKey::MaterialOuterWall,
            options,
            Some(&|| true)
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(body, original);
}
#[test]
fn interval_series_enclose_actual_linear_angle_and_square_root_derivatives() {
    use super::series::{ORDER, S, sin_cos};
    let u = S::variable(I::new(0., 0.5));
    let p = u.add(S::constant(I::point(1.))).square().sqrt().unwrap();
    assert!(p.0[0].lo <= 1. && p.0[0].hi >= 1.5);
    assert!(p.0[1].lo <= 1. && p.0[1].hi >= 1.);
    assert!(S::constant(I::point(-1.)).sqrt().is_err());
    assert!(S::constant(I::point(f64::NAN)).checked().is_err());
    let (s, c) = sin_cos(I::point(0.2), I::point(-2.), I::new(0., 0.5));
    for point in [0., 0.25, 0.5] {
        let angle: f64 = 0.2 - 2. * point;
        let mut power = 1.;
        for n in 0..=ORDER {
            if n > 0 {
                power *= -2. / n as f64;
            }
            let (ss, cc) = angle.sin_cos();
            let expected_s = power * [ss, cc, -ss, -cc][n % 4];
            let expected_c = power * [cc, -ss, -cc, ss][n % 4];
            assert!(s.0[n].lo <= expected_s && expected_s <= s.0[n].hi);
            assert!(c.0[n].lo <= expected_c && expected_c <= c.0[n].hi);
        }
    }
}

#[test]
fn analytic_queries_check_final_cancellation_without_masking_failed_precision() {
    let body = body(false, 1, true, true);
    let material = material(body);
    let count = std::cell::Cell::new(0);
    let cancel = || {
        count.set(count.get() + 1);
        count.get() >= 2
    };
    assert_eq!(
        body.local_face_properties_with_interrupt(
            FaceKey::StartCap,
            FacePropertyOptions::default(),
            Some(&cancel)
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(count.get(), 2);
    count.set(0);
    assert_eq!(
        material.local_face_properties_with_interrupt(
            FaceKey::MaterialOuterWall,
            FacePropertyOptions::default(),
            Some(&cancel)
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(count.get(), 2);
    count.set(0);
    assert_eq!(
        body.local_face_properties_with_interrupt(
            FaceKey::StartCap,
            FacePropertyOptions {
                area_tolerance_mm2: 1e-30,
                ..FacePropertyOptions::default()
            },
            Some(&cancel)
        ),
        Err(Error::Unresolved)
    );
    assert_eq!(count.get(), 1);
}

#[test]
fn dimensional_scaling_and_exact_stock_origin_leave_local_measure_authority_clear() {
    let body = body(false, 2, true, true);
    let mut spec = body.thread.spec();
    spec.major_radius_mm *= 2.;
    spec.pitch_mm *= 2.;
    spec.length_mm *= 2.;
    spec.depth_mm *= 2.;
    spec.crest_width_mm *= 2.;
    let scaled = Cavity::new(
        AuthoredThread::new(
            spec,
            SourceOptions {
                position_limit_mm: 2e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap(),
        Intent {
            entry_depth_mm: Some(1.4),
            terminal: Some(Terminal {
                width_mm: 0.7,
                radial_depth_mm: 0.3,
            }),
        },
    )
    .unwrap();
    let overlap = |a: Enclosure, b: Enclosure, k: f64| {
        assert!(a.lower <= b.upper * k && b.lower * k <= a.upper)
    };
    for index in 0..body.face_count() {
        let key = body.face_key(index).unwrap();
        let a = body
            .local_face_properties(key, FacePropertyOptions::default())
            .unwrap();
        let b = scaled
            .local_face_properties(key, FacePropertyOptions::default())
            .unwrap();
        overlap(b.area_mm2, a.area_mm2, 4.);
        for axis in 0..3 {
            overlap(b.first_moment_mm3[axis], a.first_moment_mm3[axis], 8.);
            overlap(b.centroid_mm[axis], a.centroid_mm[axis], 2.);
        }
    }
    let shifted = material(body);
    let q = |n, d| Rational::new(n, d).unwrap();
    let stock = AxialCylinder::new(
        8,
        Circle::new([q(0, 1); 2], q(5, 2)).unwrap(),
        q(0, 1),
        q(33, 8),
    )
    .unwrap();
    let centered = Material::new(stock, body, 0.2).unwrap();
    assert_ne!(shifted.stock_chart_origin(), centered.stock_chart_origin());
    for index in 0..shifted.face_count() {
        let key = shifted.face_key(index).unwrap();
        assert_eq!(
            shifted
                .local_face_properties(key, FacePropertyOptions::default())
                .unwrap(),
            centered
                .local_face_properties(key, FacePropertyOptions::default())
                .unwrap()
        );
    }
}
#[test]
fn nonpositive_or_nonfinite_receipts_and_unresolvable_adaptive_accuracy_refuse() {
    let opts = FacePropertyOptions::default();
    let face = FaceKey::StartCap;
    assert_eq!(
        finish(face, [I::point(0.); 4], opts),
        Err(Error::Unresolved)
    );
    assert_eq!(
        finish(
            face,
            [
                I::point(1.),
                I::point(f64::INFINITY),
                I::point(0.),
                I::point(0.)
            ],
            opts
        ),
        Err(Error::Unresolved)
    );
    assert_eq!(
        finish(
            face,
            [I::new(1., 2.), I::point(0.), I::point(0.), I::point(0.)],
            opts
        ),
        Err(Error::Unresolved)
    );
    let body = body(false, 1, true, true);
    assert_eq!(
        body.local_face_properties(
            FaceKey::Thread {
                start: 0,
                band: ThreadBandRole::Upper
            },
            FacePropertyOptions {
                area_tolerance_mm2: 1e-30,
                maximum_cells: 128,
                ..opts
            }
        ),
        Err(Error::Unresolved)
    );
}

#[test]
fn absent_source_supports_and_subnormal_accuracy_never_create_geometry() {
    let only_tail = body(false, 1, false, true);
    let only_entry = body(false, 1, true, false);
    let options = FacePropertyOptions::default();
    let mut work = integrate::Work::new(options, None).unwrap();
    for (body, key) in [
        (&only_tail, FaceKey::EntryCone),
        (&only_entry, FaceKey::TerminalShoulder),
        (&only_entry, FaceKey::TerminalCylinder),
        (&only_entry, FaceKey::MaterialOuterWall),
    ] {
        assert!(matches!(
            source::evaluate(source::Domain::treated(body), None, key, options, &mut work),
            Err(Error::Invalid)
        ));
    }
    let key = FaceKey::Thread {
        start: 0,
        band: ThreadBandRole::Upper,
    };
    for tolerance in [f64::from_bits(1), 1e-322, 1e-320] {
        assert_eq!(
            only_entry.local_face_properties(
                key,
                FacePropertyOptions {
                    centroid_tolerance_mm: tolerance,
                    ..options
                }
            ),
            Err(Error::Unresolved)
        );
    }
}
