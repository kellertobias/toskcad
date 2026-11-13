use super::*;
use std::cell::Cell;
fn fixture(index: usize) -> (Spec, Domain) {
    let (m, r, b, g, inner, beta, head_r, head_h) = match index {
        0 => (
            1.9,
            0.3,
            0.61,
            0.81,
            Inner::Rounded { radius_mm: 0.25 },
            7.,
            2.,
            1.6,
        ),
        1 => (
            3.,
            0.5,
            0.955,
            1.27,
            Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 0.435,
            },
            7.,
            2.8,
            2.4,
        ),
        2 => (
            3.2,
            0.5,
            0.955,
            1.27,
            Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 0.435,
            },
            7.,
            2.75,
            1.65,
        ),
        3 => (
            4.4,
            0.6,
            1.455,
            2.29,
            Inner::ShallowV {
                angle_deg: 140.,
                e_width_mm: 0.815,
            },
            5.75,
            4.,
            3.1,
        ),
        4 => (
            6.9,
            0.8,
            2.395,
            3.81,
            Inner::ShallowV {
                angle_deg: 146.,
                e_width_mm: 2.005,
            },
            5.75,
            6.,
            4.6,
        ),
        _ => (
            9.,
            1.,
            3.465,
            5.08,
            Inner::ShallowV {
                angle_deg: 153.,
                e_width_mm: 2.415,
            },
            7.,
            8.,
            6.,
        ),
    };
    let spec = Spec {
        mouth_diameter_mm: m,
        mouth_radius_mm: r,
        b_width_mm: b,
        g_diameter_mm: g,
        tilt_angle_deg: beta,
        outer_flank_angle_deg: 26.5,
        floor_angle_deg: 28.,
        outer_corner_angle_deg: 92.,
        inner,
    };
    let tip = (m - g) / (2. * 26.5_f64.to_radians().tan()) + g / 2. * 28_f64.to_radians().tan();
    (
        spec,
        Domain {
            corner_start_depth_mm: -2. * m,
            corner_end_depth_mm: tip + m,
            corner_reach_mm: 4. * (head_r + m + (tip - head_h).abs() + 1.),
        },
    )
}
fn enclosed(bounds: Bounds, value: f64) {
    assert!(
        bounds[0] <= value && bounds[1] >= value,
        "{bounds:?} must enclose {value}"
    );
}
fn rotate_point([x, y, z]: [f64; 3], quarter: usize) -> [f64; 3] {
    match quarter {
        0 => [x, y, z],
        1 => [-y, x, z],
        2 => [-x, -y, z],
        _ => [y, -x, z],
    }
}
fn constants(s: Spec) -> [f64; 5] {
    let a = s.outer_flank_angle_deg.to_radians();
    let c = s.floor_angle_deg.to_radians();
    let flank = (s.mouth_diameter_mm - s.g_diameter_mm) / (2. * a.tan());
    let tip = flank + s.g_diameter_mm / 2. * c.tan();
    let tb = flank + (s.g_diameter_mm - s.b_width_mm) / 2. * c.tan();
    let entry = s.mouth_diameter_mm / 2. + s.mouth_radius_mm * (1. - a.sin()) / a.cos();
    [entry, s.mouth_radius_mm * (1. - a.sin()), flank, tip, tb]
}
fn cc_world(s: Spec, d: Domain, [section_x, y, parameter]: [f64; 3]) -> [f64; 3] {
    let beta = s.tilt_angle_deg.to_radians();
    let tb = constants(s)[4];
    let origin = s.b_width_mm / 2. + (tb - d.corner_start_depth_mm) * beta.tan();
    [
        origin + section_x * beta.cos() - parameter * beta.tan(),
        y,
        d.corner_start_depth_mm + section_x * beta.sin() + parameter,
    ]
}
fn oracle(s: Spec, domain: Domain, p: [f64; 3]) -> bool {
    let [entry, join, flank, tip, tb] = constants(s);
    if p[2] < 0. || p[2] > tip {
        return false;
    }
    let radius = if p[2] <= join {
        entry
            - (s.mouth_radius_mm.powi(2) - (p[2] - s.mouth_radius_mm).powi(2))
                .max(0.)
                .sqrt()
    } else if p[2] <= flank {
        s.mouth_diameter_mm / 2. - p[2] * s.outer_flank_angle_deg.to_radians().tan()
    } else {
        (tip - p[2]) / s.floor_angle_deg.to_radians().tan()
    };
    if p[0].hypot(p[1]) > radius {
        return false;
    }
    let beta = s.tilt_angle_deg.to_radians();
    let theta = (s.outer_corner_angle_deg / 2.).to_radians();
    let (js, jy) = match s.inner {
        Inner::ShallowV {
            angle_deg,
            e_width_mm,
        } => (
            e_width_mm / (2. * (angle_deg / 2.).to_radians().tan()),
            e_width_mm / 2.,
        ),
        Inner::Rounded { radius_mm } => {
            (radius_mm - radius_mm * theta.sin(), radius_mm * theta.cos())
        }
    };
    let far = js + domain.corner_reach_mm * theta.cos();
    for quarter in 0..4 {
        let angle = quarter as f64 * PI / 2.;
        let (sin, cos) = angle.sin_cos();
        let x = p[0] * cos + p[1] * sin;
        let y = -p[0] * sin + p[1] * cos;
        let ox = s.b_width_mm / 2. + (tb - domain.corner_start_depth_mm) * beta.tan();
        let delta = [x - ox, y, p[2] - domain.corner_start_depth_mm];
        let section_x = delta[0] * beta.cos() + delta[2] * beta.sin();
        let parameter = delta[2] - section_x * beta.sin();
        if parameter < 0.
            || parameter > domain.corner_end_depth_mm - domain.corner_start_depth_mm
            || section_x > far
        {
            continue;
        }
        let front = if y.abs() <= jy {
            match s.inner {
                Inner::Rounded { radius_mm } => radius_mm - (radius_mm * radius_mm - y * y).sqrt(),
                Inner::ShallowV { angle_deg, .. } => y.abs() / (angle_deg / 2.).to_radians().tan(),
            }
        } else {
            js + (y.abs() - jy) / theta.tan()
        };
        if section_x >= front {
            return false;
        }
    }
    true
}

#[test]
fn metric_sources_retain_real_meridians_sections_caps_and_orthogonal_corner_owners() {
    for index in 0..6 {
        let (spec, domain) = fixture(index);
        let cutter = Cutter::new(spec, domain).unwrap();
        assert_eq!(cutter.spec(), spec);
        assert_eq!(cutter.domain(), domain);
        assert_eq!(cutter, cutter.clone());
        let e = cutter.envelope();
        let values = constants(spec);
        for (bounds, value) in [
            e.entry_radius_mm_bounds(),
            e.roundover_join_depth_mm_bounds(),
            e.flank_join_depth_mm_bounds(),
            e.tip_depth_mm_bounds(),
        ]
        .into_iter()
        .zip(values)
        {
            enclosed(bounds, value);
        }
        let segments = e.segments();
        assert_eq!(
            segments.map(|s| s.role),
            [
                EnvelopeRole::EntryRoundover,
                EnvelopeRole::OuterFlank,
                EnvelopeRole::Floor
            ]
        );
        assert_eq!(segments[0].end_mm_bounds, segments[1].start_mm_bounds);
        assert_eq!(segments[1].end_mm_bounds, segments[2].start_mm_bounds);
        let Curve::CircularArc {
            radius_mm,
            sweep_rad_bounds,
            ..
        } = segments[0].curve
        else {
            panic!("actual roundover circle");
        };
        assert_eq!(radius_mm, spec.mouth_radius_mm);
        enclosed(
            sweep_rad_bounds,
            spec.outer_flank_angle_deg.to_radians() - PI / 2.,
        );
        assert!(sweep_rad_bounds[1] < 0.);
        assert_eq!(segments[1].curve, Curve::Line);
        assert_eq!(segments[2].curve, Curve::Line);
        for (quarter, corner) in cutter.corners().iter().enumerate() {
            assert_eq!(corner.quarter_turns(), quarter as u8);
            assert_eq!(corner.domain(), domain);
            let segments = corner.section_segments().collect::<Vec<_>>();
            assert_eq!(segments.len(), if index == 0 { 4 } else { 5 });
            assert!(segments.iter().any(|s| s.role == SectionRole::Closure));
            for i in 0..segments.len() {
                assert_eq!(
                    segments[i].end_mm_bounds,
                    segments[(i + 1) % segments.len()].start_mm_bounds
                );
            }
            assert_eq!(
                segments[0].role,
                if index == 0 {
                    SectionRole::InnerRound
                } else {
                    SectionRole::InnerUpper
                }
            );
            let beta = spec.tilt_angle_deg.to_radians();
            let axis = rotate_point([-beta.sin(), 0., beta.cos()], quarter);
            for (bounds, value) in corner.extrusion_axis_bounds().into_iter().zip(axis) {
                enclosed(bounds, value);
            }
            let local = [
                0.2,
                0.13,
                (domain.corner_end_depth_mm - domain.corner_start_depth_mm) * 0.4,
            ];
            let world = rotate_point(cc_world(spec, domain, local), quarter);
            for (bounds, value) in corner
                .section_coordinates_mm_bounds(world)
                .unwrap()
                .into_iter()
                .zip(local)
            {
                enclosed(bounds, value);
            }
        }
        let bbox = cutter.bounding_box_mm();
        assert!(bbox[0][0] <= -values[0] && bbox[1][0] >= values[0] && bbox[1][2] >= values[3]);
        let mut counts = [0; 3];
        for k in 1..15 {
            let z = values[3] * k as f64 / 15.;
            for i in -9..=9 {
                for j in -9..=9 {
                    let p = [values[0] * i as f64 / 8.7, values[0] * j as f64 / 8.7, z];
                    let expected = oracle(spec, domain, p);
                    let actual = cutter.classify_point_mm(p).unwrap();
                    match actual {
                        Classification::Inside => {
                            counts[0] += 1;
                            assert!(expected, "{index} {p:?}");
                        }
                        Classification::Outside => {
                            counts[1] += 1;
                            assert!(!expected, "{index} {p:?}");
                        }
                        Classification::UnresolvedBoundary => counts[2] += 1,
                    }
                    for q in 1..4 {
                        assert_eq!(
                            cutter.classify_point_mm(rotate_point(p, q)).unwrap(),
                            actual
                        );
                    }
                }
            }
        }
        assert!(counts[0] > 50 && counts[1] > 50 && counts[2] < 100);
    }
}

#[test]
fn finite_corner_closures_and_caps_change_the_actual_cutter_and_nonconvex_sections_survive() {
    let (mut spec, mut domain) = fixture(1);
    spec.outer_corner_angle_deg = 170.;
    domain.corner_reach_mm = 20.;
    let long = Cutter::new(spec, domain).unwrap();
    let p = [1.1, 0., 0.1];
    assert_eq!(long.classify_point_mm(p).unwrap(), Classification::Outside);
    let short = Cutter::new(
        spec,
        Domain {
            corner_reach_mm: 1.,
            ..domain
        },
    )
    .unwrap();
    assert_ne!(short, long);
    assert_eq!(short.classify_point_mm(p).unwrap(), Classification::Inside);
    let beta = spec.tilt_angle_deg.to_radians();
    let tb = constants(spec)[4];
    let s = 0.12;
    let nonconvex = [
        spec.b_width_mm / 2. + (tb - 0.1) * beta.tan() + s / beta.cos(),
        0.4,
        0.1,
    ];
    assert!(s < 0.4 / (138_f64 / 2.).to_radians().tan());
    assert_eq!(
        long.classify_point_mm(nonconvex).unwrap(),
        Classification::Outside
    );
    assert!(!oracle(spec, domain, nonconvex));
    let cc_s = beta.cos() * (p[0] - spec.b_width_mm / 2.) + beta.sin() * (p[2] - tb);
    let cap_depth = p[2] - cc_s * beta.sin();
    let after = Cutter::new(
        spec,
        Domain {
            corner_end_depth_mm: cap_depth - 0.01,
            ..domain
        },
    )
    .unwrap();
    let before = Cutter::new(
        spec,
        Domain {
            corner_start_depth_mm: cap_depth + 0.01,
            ..domain
        },
    )
    .unwrap();
    assert_eq!(after.classify_point_mm(p).unwrap(), Classification::Inside);
    assert_eq!(before.classify_point_mm(p).unwrap(), Classification::Inside);
    assert_eq!(
        Cutter::new(
            spec,
            Domain {
                corner_end_depth_mm: cap_depth + 0.01,
                ..domain
            }
        )
        .unwrap()
        .classify_point_mm(p)
        .unwrap(),
        Classification::Outside
    );
    let on_cap = Cutter::new(
        spec,
        Domain {
            corner_end_depth_mm: cap_depth,
            ..domain
        },
    )
    .unwrap();
    assert_eq!(
        on_cap.classify_point_mm(p).unwrap(),
        Classification::UnresolvedBoundary
    );
    assert_eq!(
        long.classify_point_mm([0., 0., 0.]).unwrap(),
        Classification::UnresolvedBoundary
    );
    assert_eq!(
        long.classify_point_mm([1e300, 1e300, 1e300]).unwrap(),
        Classification::Outside
    );
    assert_eq!(
        long.classify_point_mm([0., 0., -0.01]).unwrap(),
        Classification::Outside
    );
}

#[test]
fn rounded_h0_uses_the_true_normal_section_and_large_tilts_can_trim_before_the_envelope_tip() {
    let (spec, domain) = fixture(0);
    let cutter = Cutter::new(spec, domain).unwrap();
    let Inner::Rounded { radius_mm: r } = spec.inner else {
        panic!("H0 circle");
    };
    let beta = spec.tilt_angle_deg.to_radians();
    let tb = constants(spec)[4];
    let y = 0.15;
    let depth = 0.4;
    let rb = spec.b_width_mm / 2. + (tb - depth) * beta.tan();
    let true_front = rb + (r - (r * r - y * y).sqrt()) / beta.cos();
    let projected_wrong = rb + r / beta.cos() - r * beta.cos() * (1. - y * y / (r * r)).sqrt();
    assert!(projected_wrong - true_front > 0.002);
    let p = [(true_front + projected_wrong) / 2., y, depth];
    assert_eq!(
        cutter.classify_point_mm(p).unwrap(),
        Classification::Outside
    );
    assert!(!oracle(spec, domain, p));
    assert_eq!(
        cutter
            .classify_point_mm([true_front - 0.001, y, depth])
            .unwrap(),
        Classification::Inside
    );
    assert_eq!(
        cutter.classify_point_mm([true_front, y, depth]).unwrap(),
        Classification::UnresolvedBoundary
    );
    let (mut steep, domain) = fixture(1);
    steep.tilt_angle_deg = 80.;
    let steep = Cutter::new(steep, domain).unwrap();
    let tip = constants(steep.spec())[3];
    assert_eq!(
        steep.classify_point_mm([0., 0., tip - 0.01]).unwrap(),
        Classification::Outside
    );
    assert_eq!(
        steep.classify_point_mm([0., 0., 0.1]).unwrap(),
        Classification::Inside
    );
}

#[test]
fn interval_classification_preserves_real_meridian_and_section_joins() {
    for index in [0, 1] {
        let (spec, domain) = fixture(index);
        let cutter = Cutter::new(spec, domain).unwrap();
        let [entry, join, flank, _, tb] = constants(spec);
        assert_eq!(
            cutter.classify_point_mm([0., 0., join]).unwrap(),
            Classification::Inside
        );
        assert_eq!(
            cutter.classify_point_mm([0., 0., flank]).unwrap(),
            Classification::Inside
        );
        let mut unresolved = 0;
        for depth in [join, flank] {
            let nominal = if depth == join {
                entry
                    - (spec.mouth_radius_mm * spec.mouth_radius_mm
                        - (depth - spec.mouth_radius_mm).powi(2))
                    .sqrt()
            } else {
                spec.g_diameter_mm / 2.
            };
            for offset in -128_i32..=128 {
                let mut radius = nominal;
                for _ in 0..offset.unsigned_abs() {
                    radius = if offset < 0 {
                        radius.next_down()
                    } else {
                        radius.next_up()
                    };
                }
                let p = [radius / 2_f64.sqrt(), radius / 2_f64.sqrt(), depth];
                match cutter.classify_point_mm(p).unwrap() {
                    Classification::Inside => assert!(oracle(spec, domain, p)),
                    Classification::Outside => assert!(!oracle(spec, domain, p)),
                    Classification::UnresolvedBoundary => unresolved += 1,
                }
            }
        }
        let beta = spec.tilt_angle_deg.to_radians();
        let theta = (spec.outer_corner_angle_deg / 2.).to_radians();
        let (js, jy) = match spec.inner {
            Inner::Rounded { radius_mm: r } => (r - r * theta.sin(), r * theta.cos()),
            Inner::ShallowV {
                angle_deg,
                e_width_mm: e,
            } => (e / (2. * (angle_deg / 2.).to_radians().tan()), e / 2.),
        };
        for offset in -64_i32..=64 {
            let mut x = spec.b_width_mm / 2. + (tb - 0.1) * beta.tan() + js / beta.cos();
            for _ in 0..offset.unsigned_abs() {
                x = if offset < 0 {
                    x.next_down()
                } else {
                    x.next_up()
                };
            }
            if cutter.classify_point_mm([x, jy, 0.1]).unwrap() == Classification::UnresolvedBoundary
            {
                unresolved += 1;
            }
        }
        assert!(unresolved > 0);
    }
}

#[test]
fn malformed_domains_source_dimensions_and_unresolved_constructions_are_explicit() {
    let (spec, domain) = fixture(1);
    let mut bad = Vec::new();
    for value in [0., -1., f64::NAN, f64::INFINITY] {
        bad.push(Spec {
            mouth_diameter_mm: value,
            ..spec
        });
        bad.push(Spec {
            mouth_radius_mm: value,
            ..spec
        });
        bad.push(Spec {
            b_width_mm: value,
            ..spec
        });
        bad.push(Spec {
            g_diameter_mm: value,
            ..spec
        });
        bad.push(Spec {
            tilt_angle_deg: value,
            ..spec
        });
        bad.push(Spec {
            outer_flank_angle_deg: value,
            ..spec
        });
        bad.push(Spec {
            floor_angle_deg: value,
            ..spec
        });
        bad.push(Spec {
            outer_corner_angle_deg: value,
            ..spec
        });
        bad.push(Spec {
            inner: Inner::Rounded { radius_mm: value },
            ..spec
        });
        bad.push(Spec {
            inner: Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: value,
            },
            ..spec
        });
    }
    bad.extend([
        Spec {
            tilt_angle_deg: 90.,
            ..spec
        },
        Spec {
            outer_flank_angle_deg: 90.,
            ..spec
        },
        Spec {
            floor_angle_deg: 90.,
            ..spec
        },
        Spec {
            outer_corner_angle_deg: 180.,
            ..spec
        },
        Spec {
            b_width_mm: spec.g_diameter_mm,
            ..spec
        },
        Spec {
            g_diameter_mm: spec.mouth_diameter_mm,
            ..spec
        },
        Spec {
            mouth_radius_mm: spec.mouth_diameter_mm / 2.,
            ..spec
        },
        Spec {
            inner: Inner::ShallowV {
                angle_deg: 90.,
                e_width_mm: 0.435,
            },
            ..spec
        },
        Spec {
            inner: Inner::ShallowV {
                angle_deg: f64::NAN,
                e_width_mm: 0.435,
            },
            ..spec
        },
        Spec {
            inner: Inner::ShallowV {
                angle_deg: 180.,
                e_width_mm: 0.435,
            },
            ..spec
        },
        Spec {
            inner: Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: spec.b_width_mm,
            },
            ..spec
        },
        Spec {
            g_diameter_mm: 3. - 0.1,
            mouth_diameter_mm: 3.,
            mouth_radius_mm: 1.4,
            ..spec
        },
        Spec {
            inner: Inner::Rounded { radius_mm: 1. },
            ..spec
        },
    ]);
    for bad in bad {
        assert!(
            matches!(Cutter::new(bad, domain), Err(Error::Invalid)),
            "{bad:?}"
        );
    }
    for bad in [
        Domain {
            corner_start_depth_mm: f64::NAN,
            ..domain
        },
        Domain {
            corner_end_depth_mm: f64::INFINITY,
            ..domain
        },
        Domain {
            corner_end_depth_mm: domain.corner_start_depth_mm,
            ..domain
        },
        Domain {
            corner_reach_mm: 0.,
            ..domain
        },
    ] {
        assert!(matches!(Cutter::new(spec, bad), Err(Error::Invalid)));
    }
    for bad in [
        Domain {
            corner_reach_mm: 1e8,
            ..domain
        },
        Domain {
            corner_start_depth_mm: -1e8,
            ..domain
        },
        Domain {
            corner_end_depth_mm: 1e8,
            ..domain
        },
    ] {
        assert!(matches!(Cutter::new(spec, bad), Err(Error::Limit)));
    }
    assert!(matches!(
        Cutter::new(
            Spec {
                mouth_diameter_mm: 1e8,
                ..spec
            },
            domain
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        Cutter::new(
            Spec {
                inner: Inner::Rounded { radius_mm: 1e8 },
                ..spec
            },
            domain
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        Cutter::new(
            Spec {
                outer_flank_angle_deg: 0.000001,
                ..spec
            },
            domain
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        Cutter::new(
            Spec {
                outer_flank_angle_deg: f64::MIN_POSITIVE,
                ..spec
            },
            domain
        ),
        Err(Error::Unresolved)
    ));
    let join = spec.mouth_radius_mm * (1. - spec.outer_flank_angle_deg.to_radians().sin());
    let near = Spec {
        g_diameter_mm: spec.mouth_diameter_mm
            - 2. * spec.outer_flank_angle_deg.to_radians().tan() * join,
        ..spec
    };
    assert!(matches!(Cutter::new(near, domain), Err(Error::Unresolved)));
    let radius = spec.b_width_mm / (2. * (spec.outer_corner_angle_deg / 2.).to_radians().cos());
    assert!(matches!(
        Cutter::new(
            Spec {
                inner: Inner::Rounded { radius_mm: radius },
                ..spec
            },
            domain
        ),
        Err(Error::Unresolved)
    ));
    let cutter = Cutter::new(spec, domain).unwrap();
    assert_eq!(
        cutter.classify_point_mm([f64::NAN, 0., 0.]),
        Err(Error::Invalid)
    );
    assert_eq!(
        cutter.corners()[0].section_coordinates_mm_bounds([f64::INFINITY, 0., 0.]),
        Err(Error::Invalid)
    );
    assert_eq!(
        cutter.corners()[0].section_coordinates_mm_bounds([f64::MAX; 3]),
        Err(Error::Unresolved)
    );
}

#[test]
fn cancellation_at_every_source_checkpoint_preserves_the_prior_immutable_source_and_retry() {
    for index in [0, 1] {
        let (spec, domain) = fixture(index);
        let source = Cutter::new(spec, domain).unwrap();
        let calls = Cell::new(0);
        let tested = Cutter::new_with_cancel(spec, domain, &|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
        assert_eq!(tested, source);
        assert!(calls.get() >= 8);
        for checkpoint in 1..=calls.get() {
            let visits = Cell::new(0);
            assert!(matches!(
                Cutter::new_with_cancel(spec, domain, &|| {
                    visits.set(visits.get() + 1);
                    visits.get() >= checkpoint
                }),
                Err(Error::Cancelled)
            ));
            assert_eq!(source.spec(), spec);
            assert_eq!(source.domain(), domain);
            assert_eq!(source, Cutter::new(spec, domain).unwrap());
        }
    }
}

#[test]
fn tiny_positive_authored_features_are_retained_as_source_expressions_without_claiming_a_mesh() {
    let (spec, domain) = fixture(1);
    let tiny = f64::from_bits(1);
    let tiny_angle = Cutter::new(
        Spec {
            outer_corner_angle_deg: tiny,
            ..spec
        },
        domain,
    )
    .unwrap();
    assert_eq!(
        tiny_angle.spec().outer_corner_angle_deg.to_bits(),
        tiny.to_bits()
    );
    assert_eq!(tiny_angle.corners()[0].section_segments().count(), 5);
    let tiny_round = Cutter::new(
        Spec {
            mouth_radius_mm: tiny,
            ..spec
        },
        domain,
    )
    .unwrap();
    let Curve::CircularArc { radius_mm, .. } = tiny_round.envelope().segments()[0].curve else {
        panic!("positive roundover is retained");
    };
    assert_eq!(radius_mm.to_bits(), tiny.to_bits());
    assert_eq!(tiny_round.envelope().segments().len(), 3);
    assert_eq!(
        tiny_round.classify_point_mm([0., 0., 0.1]).unwrap(),
        Classification::Inside
    );
    let below_resolution = Spec {
        mouth_diameter_mm: 5. * tiny,
        mouth_radius_mm: 2. * tiny,
        b_width_mm: 2. * tiny,
        g_diameter_mm: 4. * tiny,
        inner: Inner::ShallowV {
            angle_deg: 138.,
            e_width_mm: tiny,
        },
        ..spec
    };
    assert!(matches!(
        Cutter::new(below_resolution, domain),
        Err(Error::Unresolved)
    ));
}
