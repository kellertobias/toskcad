use super::*;

#[test]
fn retained_transport_precision_remains_explicit_and_consumable() {
    let p = Profile::new(counter(0.1, 0.1), 1.5, 1.).unwrap();
    let segments: Vec<_> = p.segments().collect();
    let literal = Profile::from_retained_segments(&segments).unwrap();
    assert_eq!(literal.transport_error(), TransportError::ZERO);
    assert_eq!(literal.transport_error_mm(), 0.);
    let error = TransportError {
        numeric_conversion_mm: 1e-11,
        arithmetic_mm: 1e-13,
    };
    let transported = Profile::from_transported_segments(&segments, error).unwrap();
    assert_eq!(transported.transport_error(), error);
    assert!(transported.transport_error_mm() >= 1.01e-11);
    assert!(transported.rounding_error_mm() >= literal.rounding_error_mm() + error.total_mm());
    assert_eq!(transported.volume_mm3(), literal.volume_mm3());
    assert_eq!(transported.segments().collect::<Vec<_>>(), segments);
    for bad in [f64::NAN, f64::INFINITY, -1., 1e300] {
        for error in [
            TransportError {
                numeric_conversion_mm: bad,
                arithmetic_mm: 0.,
            },
            TransportError {
                numeric_conversion_mm: 0.,
                arithmetic_mm: bad,
            },
        ] {
            assert!(Profile::from_transported_segments(&segments, error).is_err());
        }
    }
    assert!(
        Profile::from_transported_segments(
            &segments,
            TransportError {
                numeric_conversion_mm: f64::MAX,
                arithmetic_mm: f64::MAX
            }
        )
        .is_err()
    );
}

fn counter(neck: f64, chamfer: f64) -> Head {
    Head::Countersunk {
        radius_mm: 3.,
        height_mm: 1.65,
        included_angle_deg: 90.,
        edge_chamfer_mm: chamfer,
        neck_radius_mm: neck,
    }
}

fn crown(neck: f64) -> Head {
    Head::Crowned {
        radius_mm: 3.,
        height_mm: 2.5,
        crown_radius_mm: 4.,
        top_land_radius_mm: 0.8,
        neck_radius_mm: neck,
    }
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} != {b}");
}

#[test]
fn metric_counter_retains_every_source_segment_and_tangent_neck() {
    let p = Profile::new(counter(0.1, 0.1), 1.5, 1.).unwrap();
    assert_eq!(p.top(), [2.9, -1.65]);
    assert_eq!(p.end(), [1.5, 1.]);
    let s = p.segments().collect::<Vec<_>>();
    assert_eq!(
        s.iter().map(|s| s.role).collect::<Vec<_>>(),
        [
            Role::Chamfer,
            Role::Rim,
            Role::Cone,
            Role::Neck,
            Role::Shank
        ]
    );
    for pair in s.windows(2) {
        assert_eq!(pair[0].end, pair[1].start);
    }
    for s in &s {
        for i in 0..2 {
            close(s.point(0.)[i], s.start[i]);
            close(s.point(1.)[i], s.end[i]);
        }
        assert!(s.area_mm2() > 0.);
    }
    let neck = s[3];
    close(neck.end[1], 0.1 * (PI / 8.).tan());
    close(
        (neck.start[0] - 1.6).hypot(neck.start[1] - neck.end[1]),
        0.1,
    );
    close(p.radius_at(-1.6).unwrap(), 2.95);
    close(p.radius_at(-1.52).unwrap(), 3.);
    close(p.radius_at(-0.5).unwrap(), 2.);
    close(p.radius_at(0.5).unwrap(), 1.5);
    close(p.radius_at(neck.point(0.5)[1]).unwrap(), neck.point(0.5)[0]);
    close(p.support([1., 0., 0.]), 3.);
    close(p.support([0., 0., -1.]), 1.65);
    close(p.support([0., 0., 1.]), 1.);
    assert!(p.admits_socket(2. / 3_f64.sqrt(), 0.6, 0.1));
    assert!(!p.admits_socket(2.7, 1., 0.1));
    assert!(p.radius_at(-2.).is_none());
    assert!(p.radius_at(2.).is_none());
    assert!(p.radius_at(f64::NAN).is_none());
}

#[test]
fn straight_counter_matches_closed_form_and_zero_segments_are_omitted() {
    let p = Profile::new(counter(0., 0.), 1.5, 0.).unwrap();
    assert_eq!(p.segments().count(), 2);
    let rim_height = 1.65 - 1.5;
    close(
        p.volume_mm3(),
        PI * 9. * rim_height + PI * 1.5 * (9. + 4.5 + 2.25) / 3.,
    );
    close(p.radius_at(0.).unwrap(), 1.5);
    let p = Profile::new(counter(0., 0.), 1.5, 1.).unwrap();
    close(
        p.volume_mm3(),
        PI * 9. * rim_height + PI * 1.5 * (9. + 4.5 + 2.25) / 3. + PI * 2.25,
    );
}

#[test]
fn crowned_source_retains_crown_bearing_and_real_neck() {
    for n in [0., 0.1] {
        let p = Profile::new(crown(n), 1.5, 1.).unwrap();
        assert_eq!(p.top(), [0.8, -2.5]);
        close(p.radius_at(-2.5).unwrap(), 0.8);
        close(p.radius_at(0.).unwrap(), 3.);
        let s = p.segments().collect::<Vec<_>>();
        assert_eq!(s[0].role, Role::Crown);
        let sample = s[0].point(0.4);
        close(p.radius_at(sample[1]).unwrap(), sample[0]);
        assert!(s[0].area_mm2() > 0.);
        let d = [1., 0., -1.];
        let support = p.support(d);
        for i in 0..=1000 {
            let q = s[0].point(i as f64 / 1000.);
            assert!(support + 1e-12 >= q[0] - q[1]);
        }
        // A separate midpoint disk integration checks both circular-arc
        // primitives, with an error decreasing quadratically with step size.
        let steps = 100_000;
        let sum: f64 = p
            .segments()
            .map(|s| {
                let dz = (s.end[1] - s.start[1]) / steps as f64;
                (0..steps)
                    .map(|i| PI * s.radius_at(s.start[1] + (i as f64 + 0.5) * dz).powi(2) * dz)
                    .sum::<f64>()
            })
            .sum();
        assert!(
            (p.volume_mm3() - sum).abs() < 1e-7,
            "{} != {sum}",
            p.volume_mm3()
        );
        assert!(p.admits_socket(0.4, 0.5, 0.1));
        assert!(!p.admits_socket(0.75, 0.5, 0.1));
    }
}

#[test]
fn support_arc_interior_and_endpoint_and_areas_are_analytic() {
    let arc = Segment {
        start: [0., -2.],
        end: [2., 0.],
        curve: Curve::Arc {
            center: [0., 0.],
            radius_mm: 2.,
            sweep: FRAC_PI_2,
        },
        role: Role::Crown,
    };
    close(arc.support(1., -1.), 8_f64.sqrt());
    close(arc.support(1., 1.), 2.);
    close(arc.area_mm2(), 8. * PI);
    close(arc.volume_mm3(), 16. * PI / 3.);
    let neck = Segment {
        start: [3., 0.],
        end: [2., 1.],
        curve: Curve::Arc {
            center: [3., 1.],
            radius_mm: 1.,
            sweep: -FRAC_PI_2,
        },
        role: Role::Neck,
    };
    close(neck.support(1., 1.), 3.);
    close(neck.support(0., -1.), 0.);
    close(neck.area_mm2(), TAU * (3. * FRAC_PI_2 - 1.));
    close(neck.radius_at(0.5), 3. - (0.75_f64).sqrt());
    let line = Segment {
        start: [3., 0.],
        end: [2., 0.],
        curve: Curve::Line,
        role: Role::Bearing,
    };
    close(line.radius_at(0.), 3.);
    close(line.area_mm2(), 5. * PI);
    close(line.volume_mm3(), 0.);
    assert_eq!(line.point(0.5), [2.5, 0.]);
}

#[test]
fn invalid_source_never_builds_an_envelope() {
    let bad = |h, r, l| assert_eq!(Profile::new(h, r, l), Err(ThreadError::Invalid));
    for r in [f64::NAN, f64::INFINITY, 0., -1., 3.] {
        bad(counter(0., 0.), r, 1.);
    }
    for l in [f64::NAN, -1.] {
        bad(counter(0., 0.), 1.5, l);
    }
    for h in [f64::NAN, 0., -1.] {
        bad(
            Head::Countersunk {
                radius_mm: 3.,
                height_mm: h,
                included_angle_deg: 90.,
                edge_chamfer_mm: 0.,
                neck_radius_mm: 0.,
            },
            1.5,
            1.,
        );
    }
    for n in [f64::NAN, -0.1, 10.] {
        bad(counter(n, 0.), 1.5, 1.);
    }
    for ch in [f64::NAN, -0.1, 0.2, 3.] {
        bad(counter(0., ch), 1.5, 1.);
    }
    for angle in [f64::NAN, 0., -90., 180.] {
        bad(
            Head::Countersunk {
                radius_mm: 3.,
                height_mm: 1.65,
                included_angle_deg: angle,
                edge_chamfer_mm: 0.,
                neck_radius_mm: 0.,
            },
            1.5,
            1.,
        );
    }
    bad(counter(0.1, 0.), 1.5, 0.);
    for (c, l, h, n, plain) in [
        (f64::NAN, 0.8, 2.5, 0.1, 1.),
        (2., 0.8, 2.5, 0.1, 1.),
        (4., f64::NAN, 2.5, 0.1, 1.),
        (4., -1., 2.5, 0.1, 1.),
        (4., 3., 2.5, 0.1, 1.),
        (4., 0.8, 2.5, 2., 2.),
        (4., 0.8, 2.5, 0.1, 0.),
        (4., 0.8, 0.1, 0.1, 1.),
    ] {
        bad(
            Head::Crowned {
                radius_mm: 3.,
                height_mm: h,
                crown_radius_mm: c,
                top_land_radius_mm: l,
                neck_radius_mm: n,
            },
            1.5,
            plain,
        );
    }
    assert_eq!(
        Profile::new(counter(0., 0.), 1.5, 1e8),
        Err(ThreadError::Limit)
    );
    assert_eq!(
        Profile::new(
            Head::Crowned {
                radius_mm: 3.,
                height_mm: 2.5,
                crown_radius_mm: 1e8,
                top_land_radius_mm: 0.8,
                neck_radius_mm: 0.
            },
            1.5,
            1.
        ),
        Err(ThreadError::Limit)
    );
    let p = Profile::new(counter(0.1, 0.1), 1.5, 1.).unwrap();
    for (r, d, w) in [
        (f64::NAN, 0.6, 0.1),
        (0., 0.6, 0.1),
        (1., 0., 0.1),
        (1., 0.6, 0.),
        (1., 3., 0.1),
    ] {
        assert!(!p.admits_socket(r, d, w));
    }
}

#[test]
fn named_surface_moments_follow_the_actual_meridian() {
    let mut names = std::collections::BTreeSet::new();
    for head in [counter(0.1, 0.1), crown(0.1)] {
        let profile = Profile::new(head, 1.5, 1.).unwrap();
        for segment in profile.segments() {
            names.insert(segment.role.source_role());
            let [x, y, z] = segment.center_mm();
            assert_eq!([x, y], [0., 0.]);
            assert!(z >= segment.start[1] - 1e-12 && z <= segment.end[1] + 1e-12);
            let mut moment = 0.;
            let mut area = 0.;
            for i in 0..10_000 {
                let a = segment.point(i as f64 / 10_000.);
                let b = segment.point((i + 1) as f64 / 10_000.);
                let [r, z] = segment.point((i as f64 + 0.5) / 10_000.);
                let weight = r * (b[0] - a[0]).hypot(b[1] - a[1]);
                area += weight;
                moment += z * weight;
            }
            assert!((moment / area - z).abs() < 1e-8);
        }
    }
    assert_eq!(
        names,
        [
            "top.chamfer",
            "rim",
            "cone",
            "crown",
            "skirt",
            "bearing",
            "neck.blend",
            "shank"
        ]
        .into_iter()
        .collect()
    );
}

#[test]
fn erased_positive_features_and_marginal_socket_walls_are_refused() {
    for head in [
        counter(1e-30, 0.),
        counter(0., 1e-30),
        Head::Crowned {
            radius_mm: 3.,
            height_mm: 2.5,
            crown_radius_mm: 4.,
            top_land_radius_mm: 1e-30,
            neck_radius_mm: 0.,
        },
        Head::Countersunk {
            radius_mm: 1.5 + 1e-15,
            height_mm: 1.65,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.,
            neck_radius_mm: 0.,
        },
        Head::Countersunk {
            radius_mm: 3.,
            height_mm: 1.5 + 1e-15,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.,
            neck_radius_mm: 0.,
        },
        Head::Crowned {
            radius_mm: 3.,
            height_mm: 2.5,
            crown_radius_mm: 4.,
            top_land_radius_mm: 3. - 1e-13,
            neck_radius_mm: 0.,
        },
    ] {
        assert_eq!(Profile::new(head, 1.5, 1.), Err(ThreadError::Invalid));
    }
    let neck_z = 0.1 * (PI / 8.).tan();
    assert_eq!(
        Profile::new(counter(0.1, 0.), 1.5, neck_z + 1e-15),
        Err(ThreadError::Invalid)
    );
    assert_eq!(
        Profile::new(crown(0.1), 1.5, 0.1 + 1e-15),
        Err(ThreadError::Invalid)
    );
    let p = Profile::new(counter(0.1, 0.1), 1.5, 1.).unwrap();
    assert!(p.rounding_error_mm() > 0. && p.rounding_error_mm() < 1e-11);
    let floor = p.top()[1] + 0.6;
    let wall = 0.1;
    let marginal = p.radius_at(floor).unwrap() - wall;
    assert!(!p.admits_socket(marginal, 0.6, wall));
    assert!(p.admits_socket(marginal - 4. * p.rounding_error_mm(), 0.6, wall));
}

#[test]
fn representable_small_land_keeps_its_authored_circle_endpoint() {
    let land = 1e-12;
    let p = Profile::new(
        Head::Crowned {
            radius_mm: 3.,
            height_mm: 2.5,
            crown_radius_mm: 4.,
            top_land_radius_mm: land,
            neck_radius_mm: 0.,
        },
        1.5,
        1.,
    )
    .unwrap();
    assert!(land > p.rounding_error_mm());
    assert_eq!(p.radius_at(p.top()[1]), Some(land));
    let c = p.segments().next().unwrap();
    assert_eq!(c.point(0.), [land, -2.5]);
    assert_eq!(c.radius_at(c.end[1]), c.end[0]);
    assert!(c.area_mm2().is_finite() && c.volume_mm3() > 0.);
}

#[test]
fn physical_meridians_retain_actual_patches_without_an_authored_recipe() {
    for head in [counter(0.1, 0.1), crown(0.1)] {
        let original = Profile::new(head, 1.5, 1.).unwrap();
        let segments: Vec<_> = original.segments().collect();
        let retained = Profile::from_retained_segments(&segments).unwrap();
        assert_eq!(retained.segments().collect::<Vec<_>>(), segments);
        assert_eq!(retained.top(), original.top());
        assert_eq!(retained.end(), original.end());
        close(retained.volume_mm3(), original.volume_mm3());
        assert!(retained.rounding_error_mm() > 0.);
    }
    let stepped = [
        Segment {
            start: [5., -3.],
            end: [5., -1.],
            curve: Curve::Line,
            role: Role::Rim,
        },
        Segment {
            start: [5., -1.],
            end: [3., -1.],
            curve: Curve::Line,
            role: Role::Bearing,
        },
        Segment {
            start: [3., -1.],
            end: [1.5, 0.5],
            curve: Curve::Line,
            role: Role::Cone,
        },
        Segment {
            start: [1.5, 0.5],
            end: [1.5, 1.],
            curve: Curve::Line,
            role: Role::Shank,
        },
    ];
    let p = Profile::from_retained_segments(&stepped).unwrap();
    assert_eq!(p.radius_at(-1.), Some(5.));
    assert!(p.support([1., 0., 0.]) > 3.);
    // The upper step protrudes beyond the cone extended above its finite rim.
    assert!(p.radius_at(-1.1).unwrap() > 3.1);
}

#[test]
fn physical_graphs_refuse_gaps_axis_pinches_false_roles_and_invalid_arc_domains() {
    let source = Profile::new(counter(0.1, 0.1), 1.5, 1.).unwrap();
    let segments: Vec<_> = source.segments().collect();
    assert_eq!(
        Profile::from_retained_segments(&[]),
        Err(ThreadError::Invalid)
    );
    assert_eq!(
        Profile::from_retained_segments(&[segments[0]; 6]),
        Err(ThreadError::Invalid)
    );
    for (index, field, value) in [(0, 0, f64::NAN), (0, 1, f64::INFINITY), (2, 0, -1.)] {
        let mut bad = segments.clone();
        bad[index].start[field] = value;
        assert_eq!(
            Profile::from_retained_segments(&bad),
            Err(ThreadError::Invalid)
        );
    }
    let mut bad = segments.clone();
    bad[1].start[1] += 0.01;
    assert_eq!(
        Profile::from_retained_segments(&bad),
        Err(ThreadError::Invalid)
    );
    let mut bad = segments.clone();
    bad[1].role = bad[0].role;
    assert_eq!(
        Profile::from_retained_segments(&bad),
        Err(ThreadError::Invalid)
    );
    for role in [Role::Crown, Role::Neck] {
        let mut bad = segments.clone();
        bad[0].role = role;
        assert_eq!(
            Profile::from_retained_segments(&bad),
            Err(ThreadError::Invalid)
        );
    }
    for sweep in [f64::NAN, 0., FRAC_PI_2, -PI] {
        let mut bad = segments.clone();
        let Curve::Arc { center, .. } = bad[3].curve else {
            unreachable!()
        };
        bad[3].curve = Curve::Arc {
            center,
            radius_mm: 0.1,
            sweep,
        };
        assert_eq!(
            Profile::from_retained_segments(&bad),
            Err(ThreadError::Invalid)
        );
    }
    for center in [[-1., 0.], [f64::NAN, 0.], [1.5, 0.], segments[3].start] {
        let mut bad = segments.clone();
        bad[3].curve = Curve::Arc {
            center,
            radius_mm: 0.1,
            sweep: -PI / 4.,
        };
        assert_eq!(
            Profile::from_retained_segments(&bad),
            Err(ThreadError::Invalid)
        );
    }
    let mut bad = segments.clone();
    bad[3].role = Role::Bearing;
    assert_eq!(
        Profile::from_retained_segments(&bad),
        Err(ThreadError::Invalid)
    );
    let mut bad = segments.clone();
    bad[3].end[0] += 1e-5;
    bad[4].start = bad[3].end;
    assert_eq!(
        Profile::from_retained_segments(&bad),
        Err(ThreadError::Invalid)
    );
    for (radius, top, end) in [
        (1e8, -1., 0.),
        (1., 1., 2.),
        (1., -2., -1.),
        (0., -1., 0.),
        (1e-200, -1e-200, 0.),
    ] {
        let line = Segment {
            start: [radius, top],
            end: [radius, end],
            curve: Curve::Line,
            role: Role::Rim,
        };
        let result = Profile::from_retained_segments(&[line]);
        assert_eq!(
            result,
            Err(if radius > 1e7 {
                ThreadError::Limit
            } else {
                ThreadError::Invalid
            })
        );
    }
    let pinched = [
        Segment {
            start: [2., -1.],
            end: [0., 0.],
            curve: Curve::Line,
            role: Role::Cone,
        },
        Segment {
            start: [0., 0.],
            end: [1., 1.],
            curve: Curve::Line,
            role: Role::Chamfer,
        },
    ];
    assert_eq!(
        Profile::from_retained_segments(&pinched),
        Err(ThreadError::Invalid)
    );
}

#[test]
fn literal_support_radius_is_preserved_and_never_fitted_to_trim_vertices() {
    let p = Profile::new(counter(0.1, 0.1), 1.5, 1.).unwrap();
    let parts: Vec<_> = p.segments().collect();
    let Curve::Arc {
        center,
        radius_mm,
        sweep,
    } = parts[3].curve
    else {
        unreachable!()
    };
    assert_eq!(radius_mm, 0.1);
    for radius_mm in [f64::NAN, 0., -0.1, 1e-30, 0.100001, 1e8] {
        let mut bad = parts.clone();
        bad[3].curve = Curve::Arc {
            center,
            radius_mm,
            sweep,
        };
        assert!(Profile::from_retained_segments(&bad).is_err());
    }
    let arc = Segment {
        start: [0., -1.],
        end: [1., 0.],
        curve: Curve::Arc {
            center: [0.; 2],
            radius_mm: 1. + 1e-14,
            sweep: FRAC_PI_2,
        },
        role: Role::Crown,
    };
    let actual = Profile::from_retained_segments(&[arc]).unwrap();
    let q = arc.point(0.5);
    assert!((q[0].hypot(q[1]) - (1. + 1e-14)).abs() < 1e-15);
    assert!(actual.support([1., 0., 0.]) > 1.);
    assert!(actual.support([0., 0., -1.]) > 1.);
}

mod cylindrical;

#[test]
fn expanded_union_neck_respects_the_complete_axial_source_envelope() {
    assert_eq!(
        Profile::new_for_thread_union(
            Head::Crowned {
                radius_mm: 8e6,
                height_mm: 5e6,
                crown_radius_mm: 9e6,
                top_land_radius_mm: 4e6,
                neck_radius_mm: 6e6
            },
            1e6,
            0.
        ),
        Err(ThreadError::Limit)
    );
}
