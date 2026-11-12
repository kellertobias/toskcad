use super::*;

// Literal twelve-arc source for existing PanTorx A2.8/B2.05/re0.28 fixture.
fn fixture() -> [Arc; 12] {
    type SourceArc = ([f64; 2], f64, [f64; 2], [f64; 2], bool);
    let source: [SourceArc; 12] = [
        (
            [1.1199999999999999, 0.0],
            0.28,
            [1.198258138673755, -0.26884133560767637],
            [1.198258138673755, 0.26884133560767637],
            true,
        ),
        (
            [1.3462562439180192, 0.777261404824283],
            0.5295228096485661,
            [1.198258138673755, 0.26884133560767637],
            [0.8319524955604634, 0.9033013205790903],
            false,
        ),
        (
            [0.56, 0.9699484522385711],
            0.28,
            [0.8319524955604634, 0.9033013205790903],
            [0.366305643113292, 1.172142656186767],
            true,
        ),
        (
            [1.554312234475219e-16, 1.554522809648566],
            0.5295228096485661,
            [0.366305643113292, 1.172142656186767],
            [-0.3663056431132916, 1.1721426561867667],
            false,
        ),
        (
            [-0.5599999999999997, 0.9699484522385712],
            0.28,
            [-0.3663056431132916, 1.1721426561867667],
            [-0.831952495560463, 0.9033013205790905],
            true,
        ),
        (
            [-1.3462562439180192, 0.7772614048242832],
            0.5295228096485661,
            [-0.831952495560463, 0.9033013205790905],
            [-1.198258138673755, 0.2688413356076765],
            false,
        ),
        (
            [-1.1199999999999999, 1.3716044150450355e-16],
            0.28,
            [-1.198258138673755, 0.2688413356076765],
            [-1.198258138673755, -0.2688413356076762],
            true,
        ),
        (
            [-1.3462562439180195, -0.7772614048242829],
            0.5295228096485661,
            [-1.198258138673755, -0.2688413356076762],
            [-0.8319524955604636, -0.90330132057909],
            false,
        ),
        (
            [-0.5600000000000004, -0.9699484522385707],
            0.28,
            [-0.8319524955604636, -0.90330132057909],
            [-0.36630564311329245, -1.1721426561867665],
            true,
        ),
        (
            [-7.771561172376096e-16, -1.554522809648566],
            0.5295228096485661,
            [-0.36630564311329245, -1.1721426561867665],
            [0.366305643113292, -1.172142656186767],
            false,
        ),
        (
            [0.56, -0.9699484522385711],
            0.28,
            [0.366305643113292, -1.172142656186767],
            [0.8319524955604634, -0.9033013205790903],
            true,
        ),
        (
            [1.3462562439180195, -0.7772614048242827],
            0.5295228096485661,
            [0.8319524955604634, -0.9033013205790903],
            [1.198258138673755, -0.26884133560767637],
            false,
        ),
    ];
    source.map(|(center_mm, radius_mm, start_mm, end_mm, convex)| {
        let a = (start_mm[1] - center_mm[1]).atan2(start_mm[0] - center_mm[0]);
        let b = (end_mm[1] - center_mm[1]).atan2(end_mm[0] - center_mm[0]);
        let sweep_rad = if convex {
            (b - a).rem_euclid(TAU)
        } else {
            -(a - b).rem_euclid(TAU)
        };
        Arc {
            center_mm,
            radius_mm,
            start_mm,
            end_mm,
            start_angle_rad: a,
            sweep_rad,
        }
    })
}
fn drive() -> Drive {
    Drive::from_retained(
        fixture(),
        Cone {
            tip_depth_mm: 1.1 + 1.4 / 3_f64.sqrt(),
            depth_per_radius: 1. / 3_f64.sqrt(),
        },
        TransportError::ZERO,
    )
    .unwrap()
}
#[test]
fn actual_twelve_arcs_and_cone_have_distinct_finite_walls_and_bounded_mass() {
    let drive = drive();
    assert_eq!(*drive.arcs(), fixture());
    assert!(drive.source_error_mm() > 0.);
    assert!(drive.maximum_radius_bound_mm() >= 1.4);
    assert!((drive.floor_depth_mm([0., 0.]).unwrap() - drive.cone().tip_depth_mm).abs() < 1e-14);
    for i in 0..720 {
        let t = i as f64 * TAU / 720.;
        let radius = drive.radius_at_angle_mm(t).unwrap();
        assert!((1.025 - 1e-12..=1.4 + 1e-12).contains(&radius));
        let p = [radius * t.cos(), radius * t.sin()];
        let below = drive.cone().tip_depth_mm - radius * drive.cone().depth_per_radius;
        assert!(drive.contains([p[0] * 0.999, p[1] * 0.999, below * 0.999]));
        assert!(!drive.contains([p[0] * 1.001, p[1] * 1.001, 0.1]));
        assert!(!drive.contains([p[0] * 0.5, p[1] * 0.5, drive.cone().tip_depth_mm + 0.01]));
    }
    let bounds = drive.volume_bounds_mm3(1e-7, 100_000).unwrap();
    assert!(bounds[0] > 0. && bounds[1] - bounds[0] <= 1e-7);
    // Independent polar midpoint integration of the actual conical source.
    let n = 200_000;
    let mut volume = 0.;
    for i in 0..n {
        let theta = (i as f64 + 0.5) * TAU / n as f64;
        let r = drive.radius_at_angle_mm(theta).unwrap();
        volume += (drive.cone().tip_depth_mm * r * r / 2.
            - drive.cone().depth_per_radius * r * r * r / 3.)
            * TAU
            / n as f64;
    }
    assert!(bounds[0] <= volume && volume <= bounds[1]);
}

#[test]
fn source_failures_are_explicit_and_leave_the_original_cavity_unchanged() {
    let original = drive();
    for cone in [
        Cone {
            tip_depth_mm: f64::NAN,
            depth_per_radius: 1.,
        },
        Cone {
            tip_depth_mm: 1.,
            depth_per_radius: 0.,
        },
        Cone {
            tip_depth_mm: 1.,
            depth_per_radius: f64::INFINITY,
        },
        Cone {
            tip_depth_mm: 1e8,
            depth_per_radius: 1.,
        },
        Cone {
            tip_depth_mm: 2.,
            depth_per_radius: 1e8,
        },
        Cone {
            tip_depth_mm: 0.1,
            depth_per_radius: 1.,
        },
    ] {
        assert!(Drive::from_retained(fixture(), cone, TransportError::ZERO).is_err());
    }
    for transport in [
        TransportError {
            numeric_conversion_mm: -1.,
            arithmetic_mm: 0.,
        },
        TransportError {
            numeric_conversion_mm: 0.,
            arithmetic_mm: -1.,
        },
        TransportError {
            numeric_conversion_mm: f64::INFINITY,
            arithmetic_mm: 0.,
        },
    ] {
        assert_eq!(
            Drive::from_retained(fixture(), original.cone(), transport),
            Err(Error::Invalid)
        );
    }
    for changed in 0..9 {
        let mut arcs = fixture();
        match changed {
            0 => arcs[0].radius_mm = f64::NAN,
            1 => arcs[0].end_mm[0] += 0.001,
            2 => arcs[0].radius_mm = 0.,
            3 => arcs[0].sweep_rad = 0.,
            4 => arcs[0].sweep_rad = PI,
            5 => arcs[0].sweep_rad = -arcs[0].sweep_rad,
            6 => arcs[0].radius_mm += 0.001,
            7 => arcs[0].start_angle_rad = 100.,
            _ => arcs[0].center_mm[0] = 1e8,
        }
        assert!(Drive::from_retained(arcs, original.cone(), TransportError::ZERO).is_err());
    }
    let translated = fixture().map(|mut a| {
        a.center_mm[0] += 10.;
        a.start_mm[0] += 10.;
        a.end_mm[0] += 10.;
        a
    });
    assert_eq!(
        Drive::from_retained(
            translated,
            Cone {
                tip_depth_mm: 100.,
                depth_per_radius: 1.
            },
            TransportError::ZERO
        ),
        Err(Error::Unresolved)
    );
    for p in [[f64::NAN, 0.], [0., f64::INFINITY]] {
        assert_eq!(original.floor_depth_mm(p), None);
    }
    assert_eq!(original.radius_at_angle_mm(f64::NAN), None);
    assert!(!original.contains([0., 0., -0.1]));
    let disjoint = Arc {
        center_mm: [2., 0.],
        radius_mm: 1.,
        ..fixture()[0]
    };
    assert_eq!(ray_radius(disjoint, [0., 1.]), None);
    assert_eq!(original, drive());
}

#[test]
fn finite_work_cancellation_and_requested_precision_do_not_mutate_retry() {
    use std::cell::Cell;
    let drive = drive();
    for width in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(drive.volume_bounds_mm3(width, 1000), Err(Error::Invalid));
    }
    for work in [0, 1, 1_000_001] {
        assert_eq!(drive.volume_bounds_mm3(1e-7, work), Err(Error::Limit));
    }
    assert_eq!(drive.volume_bounds_mm3(1e-20, 1000), Err(Error::Unresolved));
    assert_eq!(
        drive.volume_bounds_with_cancel(1e-7, 1000, &|| true),
        Err(Error::Cancelled)
    );
    let ticks = Cell::new(0);
    assert_eq!(
        drive.volume_bounds_with_cancel(1e-7, 1000, &|| {
            ticks.set(ticks.get() + 1);
            ticks.get() > 5
        }),
        Err(Error::Cancelled)
    );
    let first = drive.volume_bounds_mm3(1e-7, 100_000).unwrap();
    assert_eq!(first, drive.volume_bounds_mm3(1e-7, 100_000).unwrap());
    let area = drive
        .arcs
        .iter()
        .fold(I::point(0.), |sum, a| sum + a.area());
    let minimum =
        2. * drive.source_volume_error().hi + (area * I::point(drive.cone.tip_depth_mm)).width();
    assert_eq!(
        drive.volume_bounds_mm3(minimum.next_up(), 100_000),
        Err(Error::Unresolved)
    );
    assert!(drive.volume_bounds_mm3(1e100, 1000).unwrap()[0] > 0.);
}

#[test]
fn unresolved_polar_stations_and_mass_enclosures_never_mint_precision() {
    assert_eq!(
        polar_span([1., 0.], [1., 1e-15], 1., 1e-12),
        Err(Error::Unresolved)
    );
    assert_eq!(
        polar_span([1., 0.], [-1., 0.], 1., 1e-12),
        Err(Error::Unresolved)
    );
    for bounds in [
        I::new(f64::NEG_INFINITY, f64::INFINITY),
        I::new(-1., 1.),
        I::new(1., 10.),
    ] {
        assert_eq!(volume_receipt(bounds, 0.1), Err(Error::Unresolved));
    }
}

#[test]
fn a_twice_wound_actual_circular_graph_is_not_a_single_recess() {
    let points: [[f64; 2]; 6] = std::array::from_fn(|i| {
        let t = i as f64 * TAU / 6.;
        [t.cos(), t.sin()]
    });
    let arcs: [Arc; 6] = std::array::from_fn(|i| {
        let start = points[i];
        let end = points[(i + 1) % 6];
        let middle = (i as f64 + 0.5) * TAU / 6.;
        let center = if i.is_multiple_of(2) {
            [0.; 2]
        } else {
            [2. * middle.cos(), 2. * middle.sin()]
        };
        let radius = (start[0] - center[0]).hypot(start[1] - center[1]);
        let a = (start[1] - center[1]).atan2(start[0] - center[0]);
        let b = (end[1] - center[1]).atan2(end[0] - center[0]);
        let sweep = if i.is_multiple_of(2) {
            (b - a).rem_euclid(TAU)
        } else {
            -(a - b).rem_euclid(TAU)
        };
        Arc {
            center_mm: center,
            radius_mm: radius,
            start_mm: start,
            end_mm: end,
            start_angle_rad: a,
            sweep_rad: sweep,
        }
    });
    assert_eq!(
        Drive::from_retained(
            std::array::from_fn(|i| arcs[i % 6]),
            Cone {
                tip_depth_mm: 10.,
                depth_per_radius: 1.
            },
            TransportError::ZERO
        ),
        Err(Error::Unresolved)
    );
}

#[test]
fn retained_rotations_metric_scales_and_transport_keep_actual_circular_sources() {
    let baseline = drive().volume_bounds_mm3(1e-7, 100_000).unwrap();
    for scale in [0.001, 1., 1000.] {
        for rotation in [0_f64, 0.63] {
            let (s, c) = rotation.sin_cos();
            let point = |[x, y]: [f64; 2]| [scale * (c * x - s * y), scale * (s * x + c * y)];
            let arcs = fixture().map(|a| Arc {
                center_mm: point(a.center_mm),
                start_mm: point(a.start_mm),
                end_mm: point(a.end_mm),
                radius_mm: a.radius_mm * scale,
                start_angle_rad: a.start_angle_rad + rotation,
                ..a
            });
            let cone = Cone {
                tip_depth_mm: drive().cone().tip_depth_mm * scale,
                depth_per_radius: drive().cone().depth_per_radius,
            };
            let source = Drive::from_retained(
                arcs,
                cone,
                TransportError {
                    numeric_conversion_mm: 2e-14 * scale,
                    arithmetic_mm: 1e-14 * scale,
                },
            )
            .unwrap();
            assert_eq!(*source.arcs(), arcs);
            let bounds = source
                .volume_bounds_mm3(1e-6 * scale.powi(3), 100_000)
                .unwrap();
            assert!(
                bounds[0] <= baseline[1] * scale.powi(3)
                    && baseline[0] * scale.powi(3) <= bounds[1]
            );
            let p = point([0.5, 0.3]);
            assert!(source.contains([p[0], p[1], scale]));
            assert!(source.source_error_mm() > drive().source_error_mm() * scale.min(1.));
        }
    }
}

#[test]
fn unresolved_near_tangent_ray_and_polar_trims_are_refused_without_radius_clamping() {
    let mut arcs = fixture();
    let start = 2. * PI / 3. - 1e-9;
    let sweep = 1e-10;
    let point = |t: f64| [2. + t.cos(), t.sin()];
    arcs[0] = Arc {
        center_mm: [2., 0.],
        radius_mm: 1.,
        start_mm: point(start),
        end_mm: point(start + sweep),
        start_angle_rad: start,
        sweep_rad: sweep,
    };
    arcs[1].start_mm = arcs[0].end_mm;
    arcs[11].end_mm = arcs[0].start_mm;
    assert_eq!(
        Drive::from_retained(
            arcs,
            Cone {
                tip_depth_mm: 10.,
                depth_per_radius: 1.
            },
            TransportError::ZERO
        ),
        Err(Error::Unresolved)
    );
    let points: [[f64; 2]; 12] = std::array::from_fn(|i| {
        let t = i as f64 * TAU / 12.;
        [t.cos(), t.sin()]
    });
    let mut unresolved = 0;
    for offset in [1e-11, 2e-11, 4e-11] {
        let contour = std::array::from_fn(|i: usize| {
            let middle = (i as f64 + 0.5) * TAU / 12.;
            let distance = (PI / 12.).cos().recip() + offset;
            let center = if i.is_multiple_of(2) {
                [0.; 2]
            } else {
                [distance * middle.cos(), distance * middle.sin()]
            };
            let from = points[i];
            let to = points[(i + 1) % 12];
            let a = (from[1] - center[1]).atan2(from[0] - center[0]);
            let b = (to[1] - center[1]).atan2(to[0] - center[0]);
            let radius = (from[0] - center[0]).hypot(from[1] - center[1]);
            Arc {
                center_mm: center,
                radius_mm: radius,
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
            contour,
            Cone {
                tip_depth_mm: 10.,
                depth_per_radius: 1.,
            },
            TransportError::ZERO,
        )
        .unwrap();
        for p in points {
            let angle = p[1].atan2(p[0]);
            if source.radius_at_angle_mm(angle).is_none() {
                assert_eq!(source.floor_depth_mm(p), None);
                unresolved += 1;
            }
        }
    }
    assert!(unresolved > 0);
}

#[test]
fn authored_tangent_authority_is_explicit_and_never_inferred_from_retained_arcs() {
    let spec = TangentSpec {
        major_diameter_mm: 2.8,
        minor_diameter_mm: 2.05,
        outer_lobe_radius_mm: 0.28,
    };
    let source = Drive::from_authored_tangent(spec, drive().cone()).unwrap();
    assert_eq!(
        source.boundary_source(),
        BoundarySource::AuthoredTangent(spec)
    );
    assert_eq!(drive().boundary_source(), BoundarySource::RetainedFamily);
    let replay = Drive::from_retained(*source.arcs(), source.cone(), TransportError::ZERO).unwrap();
    assert_eq!(replay.boundary_source(), BoundarySource::RetainedFamily);
    for (arc, old) in source.arcs().iter().zip(fixture()) {
        for (a, b) in [arc.center_mm, arc.start_mm, arc.end_mm].into_iter().zip([
            old.center_mm,
            old.start_mm,
            old.end_mm,
        ]) {
            assert!((a[0] - b[0]).hypot(a[1] - b[1]) < source.source_error_mm());
        }
        assert!((arc.radius_mm - old.radius_mm).abs() < source.source_error_mm());
    }
    let actual = source.volume_bounds_mm3(1e-7, 100_000).unwrap();
    let previous = drive().volume_bounds_mm3(1e-7, 100_000).unwrap();
    assert!(actual[0] <= previous[1] && previous[0] <= actual[1]);
    for bad in [
        TangentSpec {
            major_diameter_mm: f64::NAN,
            ..spec
        },
        TangentSpec {
            minor_diameter_mm: 3.,
            ..spec
        },
        TangentSpec {
            outer_lobe_radius_mm: 1.4,
            ..spec
        },
        TangentSpec {
            major_diameter_mm: 1e8,
            minor_diameter_mm: 1e7,
            ..spec
        },
        TangentSpec {
            minor_diameter_mm: 0.1,
            ..spec
        },
        TangentSpec {
            minor_diameter_mm: 2.8 * 2. * (0.4 * (3_f64.sqrt() / 2.) + 0.1),
            ..spec
        },
    ] {
        assert!(Drive::from_authored_tangent(bad, source.cone()).is_err());
    }
    assert_eq!(
        Drive::from_authored_tangent(
            spec,
            Cone {
                tip_depth_mm: -1.,
                depth_per_radius: 1.
            }
        ),
        Err(Error::Invalid)
    );
}

#[test]
fn finite_support_lower_bounds_cover_outward_nodes_at_adjacent_endpoint_cells() {
    for (a, sweep) in fixture()
        .map(|arc| (arc.start_angle_rad, arc.sweep_rad))
        .into_iter()
        .chain([
            (0., 1e-10),
            (-f64::MIN_POSITIVE, 2. * f64::MIN_POSITIVE),
            (TAU, PI - 1e-6),
            (-TAU, -PI + 1e-6),
        ])
    {
        let domain = calculation_domain(a, sweep);
        let b = a + sweep;
        let lo = a.min(b);
        let hi = a.max(b);
        let end = I::point(a) + I::point(sweep);
        assert!(domain.lo <= end.lo.min(b) && end.hi.max(b) <= domain.hi);
        let mid = lo + (hi - lo) * 0.5;
        for (lo, hi) in [
            (lo, hi),
            (lo, lo.next_up()),
            (hi.next_down(), hi),
            (mid, mid.next_up()),
        ] {
            let (_, nodes) = gauss_nodes(lo, hi);
            for node in nodes {
                assert!(domain.lo <= node.lo && node.hi <= domain.hi);
            }
        }
    }
}
