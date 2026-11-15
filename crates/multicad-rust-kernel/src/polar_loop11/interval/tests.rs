use super::*;
fn registry(clockwise: bool, width: f64) -> (Vec<[[f64; 2]; 2]>, Vec<Vec<usize>>) {
    let poles = super::super::tests::diamond(1, clockwise, 0.);
    let box_ = |p: P| [[p[0] - width, p[0] + width], [p[1] - width, p[1] + width]];
    let mut points: Vec<_> = poles.iter().map(|p| box_(p[0])).collect();
    let indices = poles
        .iter()
        .enumerate()
        .map(|(curve, poles)| {
            let mut indices = vec![curve];
            for &p in &poles[1..11] {
                indices.push(points.len());
                points.push(box_(p));
            }
            indices.push((curve + 1) % 4);
            indices
        })
        .collect();
    (points, indices)
}
fn curve_views(indices: &[Vec<usize>]) -> Vec<IntervalCurve<'_>> {
    indices
        .iter()
        .enumerate()
        .map(|(i, c)| IntervalCurve {
            curve_owner: 100 + i,
            spans: 1,
            point_indices: c,
        })
        .collect()
}
#[test]
fn shared_interval_variables_certify_every_assignment_without_a_midpoint_or_closure_proxy() {
    for clockwise in [false, true] {
        for width in [0., 1e-12] {
            let (points, indices) = registry(clockwise, width);
            let curves = curve_views(&indices);
            let input = IntervalInput {
                points: &points,
                curves: &curves,
                counter_clockwise: !clockwise,
            };
            let result = certify_intervals(input, Default::default(), None).unwrap();
            assert!(std::ptr::eq(
                result.input().points.as_ptr(),
                points.as_ptr()
            ));
            assert_eq!(result.winding(), if clockwise { -1 } else { 1 });
            assert!(
                result.angular_advance()[0] < std::f64::consts::TAU
                    && result.angular_advance()[1] > std::f64::consts::TAU
            );
            assert!(result.nodes() >= result.leaves() && result.leaves() > 4);
            assert_eq!(
                result.input_budget_bytes(),
                std::mem::size_of_val(points.as_slice())
                    + 48 * std::mem::size_of::<usize>()
                    + std::mem::size_of_val(curves.as_slice())
            );
            assert_eq!(result.scratch_budget_bytes(), SCRATCH);
        }
    }
    let (mut points, mut indices) = registry(false, 1e-12);
    points.push(points[0]);
    *indices[3].last_mut().unwrap() = points.len() - 1;
    let curves = curve_views(&indices);
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &points,
                curves: &curves,
                counter_clockwise: true
            },
            Default::default(),
            None
        ),
        Err(Error::Invalid)
    ));
    let (points, indices) = registry(false, 0.);
    let curves = curve_views(&indices);
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &points,
                curves: &curves,
                counter_clockwise: false
            },
            Default::default(),
            None
        ),
        Err(Error::Unresolved)
    ));
    let wide = vec![[[-20., 20.]; 2]; points.len()];
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &wide,
                curves: &curves,
                counter_clockwise: true
            },
            Options {
                maximum_depth: 0,
                ..Default::default()
            },
            None
        ),
        Err(Error::Unresolved)
    ));
}
#[test]
fn interval_registry_queries_honor_real_work_node_byte_depth_and_cancellation_frontiers() {
    let (points, indices) = registry(false, 1e-12);
    let curves = curve_views(&indices);
    let input = IntervalInput {
        points: &points,
        curves: &curves,
        counter_clockwise: true,
    };
    let options = Options::default();
    let result = certify_intervals(input, options, None).unwrap();
    for bad in [
        Options {
            maximum_work: result.work_used() - 1,
            ..options
        },
        Options {
            maximum_nodes: result.nodes() - 1,
            ..options
        },
        Options {
            maximum_bytes: result.input_budget_bytes() + SCRATCH - 1,
            ..options
        },
    ] {
        assert!(matches!(
            certify_intervals(input, bad, None),
            Err(Error::Limit)
        ));
    }
    let exact = certify_intervals(
        input,
        Options {
            maximum_work: result.work_used(),
            maximum_nodes: result.nodes(),
            maximum_bytes: result.input_budget_bytes() + SCRATCH,
            ..options
        },
        None,
    )
    .unwrap();
    assert_eq!(exact.work_used(), result.work_used());
    let polls = std::cell::Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    certify_intervals(input, options, Some(&count)).unwrap();
    for at in 1..=polls.get() {
        let seen = std::cell::Cell::new(0usize);
        let cancel = || {
            seen.set(seen.get() + 1);
            seen.get() == at
        };
        assert!(matches!(
            certify_intervals(input, options, Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    for bad in [
        Options {
            maximum_work: WORK + 1,
            ..options
        },
        Options {
            maximum_nodes: NODES + 1,
            ..options
        },
        Options {
            maximum_bytes: BYTES + 1,
            ..options
        },
        Options {
            maximum_depth: 33,
            ..options
        },
    ] {
        assert!(matches!(
            certify_intervals(input, bad, None),
            Err(Error::Invalid)
        ));
    }
}
#[test]
fn malformed_or_oversized_interval_inputs_refuse_without_a_geometry_claim() {
    let (points, indices) = registry(false, 0.);
    let valid_curves = curve_views(&indices);
    let options = Options::default();
    for input in [
        IntervalInput {
            points: &points,
            curves: &[],
            counter_clockwise: true,
        },
        IntervalInput {
            points: &[],
            curves: &valid_curves,
            counter_clockwise: true,
        },
    ] {
        assert!(matches!(
            certify_intervals(input, options, None),
            Err(Error::Invalid)
        ));
    }
    let many = vec![valid_curves[0]; CURVES + 1];
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &points,
                curves: &many,
                counter_clockwise: true
            },
            options,
            None
        ),
        Err(Error::Invalid)
    ));
    let many_points = vec![[[1., 1.], [0., 0.]]; POINTS + 1];
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &many_points,
                curves: &valid_curves,
                counter_clockwise: true
            },
            options,
            None
        ),
        Err(Error::Limit)
    ));
    let too_many = vec![0; POINTS + 1];
    let huge = [IntervalCurve {
        curve_owner: 1,
        spans: 1,
        point_indices: &too_many,
    }];
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &points,
                curves: &huge,
                counter_clockwise: true
            },
            options,
            None
        ),
        Err(Error::Limit)
    ));
    let half = vec![0; 11 * (POINTS / 22 + 1) + 1];
    let huge = [
        IntervalCurve {
            curve_owner: 1,
            spans: half.len() / 11,
            point_indices: &half,
        },
        IntervalCurve {
            curve_owner: 2,
            spans: half.len() / 11,
            point_indices: &half,
        },
    ];
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &points,
                curves: &huge,
                counter_clockwise: true
            },
            options,
            None
        ),
        Err(Error::Limit)
    ));
    for bad_point in [
        [[f64::NAN, 0.], [0., 0.]],
        [[2., 1.], [0., 0.]],
        [[0., 1e8], [0., 0.]],
    ] {
        let mut bad = points.clone();
        bad[0] = bad_point;
        assert!(matches!(
            certify_intervals(
                IntervalInput {
                    points: &bad,
                    curves: &valid_curves,
                    counter_clockwise: true
                },
                options,
                None
            ),
            Err(Error::Invalid)
        ));
    }
    let mut bad = valid_curves.clone();
    bad[0].spans = 0;
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &points,
                curves: &bad,
                counter_clockwise: true
            },
            options,
            None
        ),
        Err(Error::Invalid)
    ));
    let mut bad_indices = indices.clone();
    bad_indices[0][1] = points.len();
    let bad = curve_views(&bad_indices);
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &points,
                curves: &bad,
                counter_clockwise: true
            },
            options,
            None
        ),
        Err(Error::Invalid)
    ));
    let mut bad = valid_curves.clone();
    bad[2].curve_owner = bad[0].curve_owner;
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &points,
                curves: &bad,
                counter_clockwise: true
            },
            options,
            None
        ),
        Err(Error::Invalid)
    ));
    let mut bad_indices = indices.clone();
    *bad_indices[1].last_mut().unwrap() = 0;
    bad_indices[2][0] = 0;
    let bad = curve_views(&bad_indices);
    assert!(matches!(
        certify_intervals(
            IntervalInput {
                points: &points,
                curves: &bad,
                counter_clockwise: true
            },
            options,
            None
        ),
        Err(Error::Invalid)
    ));
}
