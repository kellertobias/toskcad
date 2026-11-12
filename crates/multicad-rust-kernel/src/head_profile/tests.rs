use super::*;

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
