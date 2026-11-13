use super::*;
fn lifted(index: usize, posed: bool) -> (usize, usize) {
    let text = super::super::tests::text(index, posed);
    let b = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let endpoints = b
        .endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 1_000_000,
            },
            None,
        )
        .unwrap();
    let analytic = b
        .analytic_support_bounds(
            &endpoints,
            AnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    let lifts = analytic
        .finite_edge_lifts(endpoints, LiftOptions::default(), None)
        .unwrap_or_else(|e| panic!("metric{index} posed{posed} liftrefused{e:?}"));
    assert!(std::ptr::eq(lifts.boundary(), &b));
    assert!(std::ptr::eq(lifts.endpoint_bounds().boundary(), &b));
    assert!(lifts.work_used() < WORK);
    assert!(
        lifts.retained_capacity_bytes() > std::mem::size_of_val(&lifts)
            && lifts.retained_capacity_bytes() < BYTES
    );
    let mut lines = 0;
    let mut arcs = 0;
    eprintln!("first retained typed analytic lift: {:?}", lifts.edges()[0]);
    for lift in lifts.edges() {
        let edge = &b.edges()[lift.edge_index()];
        match lift.kind() {
            LiftKind::Line { parameters: p } => {
                lines += 1;
                assert!(if edge.same_sense {
                    p[1][0] > p[0][1]
                } else {
                    p[1][1] < p[0][0]
                });
            }
            LiftKind::Arc {
                start_angle_rad: a,
                signed_sweep_rad: s,
            } => {
                arcs += 1;
                assert!(a[0].is_finite() && a[0] <= a[1]);
                assert!(if edge.same_sense {
                    s[0] > 0. && s[1] < tau().lo
                } else {
                    s[1] < 0. && s[0] > -tau().lo
                });
            }
        }
    }
    assert_eq!(lines + arcs, 44);
    (lines, arcs)
}
#[test]
fn actual_all_metric_finite_analytic_edges_keep_vector_magnitude_senses_and_periods() {
    for index in 0..6 {
        for posed in [false, true] {
            let counts = lifted(index, posed);
            eprintln!("lift metric{index} posed{posed}: {counts:?}");
        }
    }
}
#[test]
fn angle_charts_keep_all_quadrants_antipodes_and_ambiguous_small_arcs_explicit() {
    let xy = |x, y| [I::point(x), I::point(y)];
    for (p, expected) in [
        (xy(1., 1.), std::f64::consts::FRAC_PI_4),
        (xy(-1., 1.), 3. * std::f64::consts::FRAC_PI_4),
        (xy(-1., -1.), 5. * std::f64::consts::FRAC_PI_4),
        (xy(1., -1.), -std::f64::consts::FRAC_PI_4),
        (xy(0., 1.), std::f64::consts::FRAC_PI_2),
        (xy(0., -1.), -std::f64::consts::FRAC_PI_2),
    ] {
        let a = angle(p[1], p[0]).unwrap();
        assert!(a.lo <= expected && expected <= a.hi);
    }
    assert!(angle(I::new(-1., 1.), I::new(-1., 1.)).is_err());
    assert!(scaled(xy(0., 0.)).is_err());
    let first = xy(1., 0.);
    for forward in [false, true] {
        let half = directed_sweep(first, xy(-1., 0.), forward).unwrap();
        let expected = if forward {
            std::f64::consts::PI
        } else {
            -std::f64::consts::PI
        };
        assert!(half.lo <= expected && expected <= half.hi);
        let long = directed_sweep(
            first,
            if forward { xy(0., -1.) } else { xy(0., 1.) },
            forward,
        )
        .unwrap();
        let expected = if forward {
            1.5 * std::f64::consts::PI
        } else {
            -1.5 * std::f64::consts::PI
        };
        assert!(long.lo <= expected && expected <= long.hi);
        assert!(directed_sweep(first, first, forward).is_err());
        assert!(
            directed_sweep(
                first,
                [I::new(1. - 1e-15, 1. + 1e-15), I::new(-1e-15, 1e-15)],
                forward
            )
            .is_err()
        );
    }
    assert!(directed_sweep(first, xy(1., -16. * f64::EPSILON * 1.1), true).is_err());
}

#[test]
fn independent_vector_magnitude_and_ellipse_parameter_oracles_use_original_geometry() {
    let mut budget = Budget {
        remaining: WORK,
        interrupt: None,
    };
    let origin = [3., -4., 5.];
    let points = [origin, [3., -4., 9.]];
    for forward in [false, true] {
        let p = if forward {
            points
        } else {
            [points[1], points[0]]
        };
        let LiftKind::Line { parameters } =
            line(origin, [0., 0., 10.], 2., p, forward, &mut budget)
                .unwrap()
                .kind()
        else {
            panic!()
        };
        for (i, expected) in if forward { [0., 2.] } else { [2., 0.] }
            .into_iter()
            .enumerate()
        {
            assert!(parameters[i][0] <= expected && expected <= parameters[i][1]);
        }
    }
    assert!(
        line(
            origin,
            [0., 0., 10.],
            2.,
            [origin, origin],
            true,
            &mut budget
        )
        .is_err()
    );
    let axis = Axis {
        origin: [0.; 3],
        normal: [0., 0., 7.],
        reference: [5., 0., 0.],
    };
    let LiftKind::Arc {
        start_angle_rad: a,
        signed_sweep_rad: s,
    } = arc(
        axis,
        Some([2., 1.]),
        [[2., 0., 0.], [0., 1., 0.]],
        true,
        &mut budget,
    )
    .unwrap()
    .kind()
    else {
        panic!()
    };
    assert!(a[0] <= 0. && a[1] >= 0.);
    assert!(s[0] <= std::f64::consts::FRAC_PI_2 && std::f64::consts::FRAC_PI_2 <= s[1]);
    let tiny = scaled([I::point(1e-300), I::point(2e-300)]).unwrap();
    let huge = scaled([I::point(1e300), I::point(2e300)]).unwrap();
    assert!(angle(tiny[1], tiny[0]).unwrap().lo <= angle(huge[1], huge[0]).unwrap().hi);
    assert!(scaled([I::point(f64::INFINITY), I::point(1.)]).is_err());
}
#[test]
fn borrowed_finite_lifts_keep_foreign_arena_resources_and_all_cancel_stages_closed() {
    let b = read_phillips_screw_boundary(&super::super::tests::text(0, false), 1e-7, 1_000_000)
        .unwrap();
    let foreign = b.clone();
    fn endpoint(b: &Boundary) -> EndpointCouplingBounds<'_> {
        b.endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 1_000_000,
            },
            None,
        )
        .unwrap()
    }
    let query = |o: LiftOptions, interrupt: Option<&dyn Fn() -> bool>| {
        let e = endpoint(&b);
        let a = b
            .analytic_support_bounds(
                &e,
                AnalyticSupportOptions {
                    maximum_distance_mm: 1e-7,
                    work_limit: 10_000_000,
                },
                None,
            )
            .unwrap();
        a.finite_edge_lifts(e, o, interrupt)
    };
    let a = b
        .analytic_support_bounds(
            &endpoint(&b),
            AnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    assert!(matches!(
        a.finite_edge_lifts(endpoint(&foreign), LiftOptions::default(), None),
        Err(PhillipsReadError::Invalid)
    ));
    let o = LiftOptions::default();
    for v in [
        LiftOptions {
            maximum_work: 0,
            ..o
        },
        LiftOptions {
            maximum_work: WORK + 1,
            ..o
        },
        LiftOptions {
            maximum_bytes: 0,
            ..o
        },
        LiftOptions {
            maximum_bytes: BYTES + 1,
            ..o
        },
    ] {
        assert!(matches!(query(v, None), Err(PhillipsReadError::Invalid)));
    }
    for v in [
        LiftOptions {
            maximum_work: 1,
            ..o
        },
        LiftOptions {
            maximum_bytes: 1,
            ..o
        },
    ] {
        assert!(matches!(query(v, None), Err(PhillipsReadError::Limit)));
    }
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let complete = query(o, Some(&count)).unwrap();
    let total = calls.get();
    for stop in 0..total {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() > stop
        };
        assert!(matches!(
            query(o, Some(&cancel)),
            Err(PhillipsReadError::Cancelled)
        ));
    }
    assert_eq!(
        query(o, None).unwrap().edges().len(),
        complete.edges().len()
    );
}

#[test]
fn a_real_line_same_sense_flip_cannot_change_the_original_vector_to_repair_order() {
    let text = super::super::tests::text(0, false);
    let original = read_phillips_screw_boundary(&text, 1e-7, 1_000_000).unwrap();
    let e = original
        .edges()
        .iter()
        .find(|e| matches!(e.curve, Curve::Line { .. }))
        .unwrap();
    let prefix = format!("#{}=", e.entity);
    let (old, new) = if e.same_sense {
        (".T.", ".F.")
    } else {
        (".F.", ".T.")
    };
    let changed = text
        .lines()
        .map(|line| {
            if line.starts_with(&prefix) {
                line.replacen(old, new, 1)
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert_ne!(text, changed);
    let b = read_phillips_screw_boundary(&changed, 1e-7, 1_000_000).unwrap();
    let endpoints = b
        .endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 1_000_000,
            },
            None,
        )
        .unwrap();
    let analytic = b
        .analytic_support_bounds(
            &endpoints,
            AnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .unwrap();
    assert!(matches!(
        analytic.finite_edge_lifts(endpoints, LiftOptions::default(), None),
        Err(PhillipsReadError::Unsupported)
    ));
}

#[test]
fn primitive_frame_evaluation_obeys_actual_cancellation_before_parameter_publication() {
    let axis = Axis {
        origin: [0.; 3],
        normal: [0., 0., 1.],
        reference: [1., 0., 0.],
    };
    let cancel = || true;
    let mut budget = Budget {
        remaining: WORK,
        interrupt: Some(&cancel),
    };
    assert!(matches!(
        arc(axis, None, [[1., 0., 0.], [0., 1., 0.]], true, &mut budget),
        Err(PhillipsReadError::Cancelled)
    ));
}
