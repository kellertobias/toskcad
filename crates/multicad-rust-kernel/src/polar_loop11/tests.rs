use super::*;
pub(super) fn diamond(turns: usize, clockwise: bool, z: f64) -> Vec<Vec<P>> {
    let vertices = if clockwise {
        [[11., 0.], [0., -11.], [-11., 0.], [0., 11.]]
    } else {
        [[11., 0.], [0., 11.], [-11., 0.], [0., -11.]]
    };
    (0..4 * turns)
        .map(|edge| {
            let (a, b) = (vertices[edge % 4], vertices[(edge + 1) % 4]);
            (0..12)
                .map(|i| {
                    let t = i as f64;
                    [
                        (a[0] * (11. - t) + b[0] * t) / 11.,
                        (a[1] * (11. - t) + b[1] * t) / 11.,
                        z,
                    ]
                })
                .collect()
        })
        .collect()
}
pub(super) fn curves(poles: &[Vec<P>]) -> Vec<Curve<'_>> {
    poles
        .iter()
        .enumerate()
        .map(|(i, p)| Curve {
            curve_owner: 100 + i,
            start_owner: i,
            end_owner: (i + 1) % poles.len(),
            spans: 1,
            poles: p,
        })
        .collect()
}
#[test]
fn literal_diamond_both_orientations_proves_closed_onefold_radial_jordan() {
    for clockwise in [false, true] {
        let poles = diamond(1, clockwise, 3.125);
        let curves = curves(&poles);
        let input = Input {
            curves: &curves,
            counter_clockwise: !clockwise,
        };
        let c = certify(input, Default::default(), None).unwrap();
        assert!(std::ptr::eq(c.input().curves.as_ptr(), curves.as_ptr()));
        assert_eq!(c.plane_z(), 3.125);
        assert_eq!(c.winding(), if clockwise { -1 } else { 1 });
        assert!(
            c.angular_advance()[0] < 2. * std::f64::consts::PI
                && c.angular_advance()[1] > 2. * std::f64::consts::PI
        );
        assert!(c.nodes() >= c.leaves() && c.leaves() > 4);
        assert_eq!(
            c.input_budget_bytes(),
            4 * 12 * std::mem::size_of::<P>() + std::mem::size_of_val(curves.as_slice())
        );
        assert_eq!(c.scratch_budget_bytes(), SCRATCH);
    }
    let poles = diamond(2, false, 0.);
    let curves = curves(&poles);
    assert!(matches!(
        certify(
            Input {
                curves: &curves,
                counter_clockwise: true
            },
            Default::default(),
            None
        ),
        Err(Error::Unresolved)
    ));
}
#[test]
fn arctangent_has_an_independent_exact_rational_series_oracle() {
    use num_bigint::BigInt as Z;
    use num_rational::BigRational as Q;
    for value in [0.01, 0.25, 0.5] {
        let r = Q::from_float(value).unwrap();
        let square = &r * &r;
        let mut power = r.clone();
        let mut sum = r.clone();
        for k in 1..=40 {
            power *= &square;
            let term = &power / Q::from_integer(Z::from(2 * k + 1));
            if k % 2 == 0 {
                sum += term;
            } else {
                sum -= term;
            }
        }
        let remainder = &power * &square / Q::from_integer(Z::from(83));
        let actual = atan_small(I::point(value));
        assert!(Q::from_float(actual.lo).unwrap() <= &sum - &remainder);
        assert!(Q::from_float(actual.hi).unwrap() >= sum);
    }
}
#[test]
fn work_node_bytes_depth_and_every_cancellation_frontier_are_real() {
    let poles = diamond(1, false, 0.);
    let curves = curves(&poles);
    let input = Input {
        curves: &curves,
        counter_clockwise: true,
    };
    let defaults = Options::default();
    let baseline = certify(input, defaults, None).unwrap();
    for o in [
        Options {
            maximum_work: baseline.work_used() - 1,
            ..defaults
        },
        Options {
            maximum_nodes: baseline.nodes() - 1,
            ..defaults
        },
        Options {
            maximum_bytes: baseline.input_budget_bytes() + SCRATCH - 1,
            ..defaults
        },
    ] {
        assert!(matches!(certify(input, o, None), Err(Error::Limit)));
    }
    let exact = certify(
        input,
        Options {
            maximum_work: baseline.work_used(),
            maximum_nodes: baseline.nodes(),
            maximum_bytes: baseline.input_budget_bytes() + SCRATCH,
            ..defaults
        },
        None,
    )
    .unwrap();
    assert_eq!(exact.work_used(), baseline.work_used());
    assert!(matches!(
        certify(
            input,
            Options {
                maximum_depth: 0,
                ..defaults
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    let polls = std::cell::Cell::new(0usize);
    let poll = || {
        polls.set(polls.get() + 1);
        false
    };
    certify(input, defaults, Some(&poll)).unwrap();
    for at in 1..=polls.get() {
        let calls = std::cell::Cell::new(0usize);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == at
        };
        assert!(
            matches!(
                certify(input, defaults, Some(&cancel)),
                Err(Error::Cancelled)
            ),
            "checkpoint{at}"
        );
    }
}
#[test]
fn invalid_inputs_and_actual_graph_numeric_closure_refuse() {
    let defaults = Options::default();
    let poles = diamond(1, false, 0.);
    let original = curves(&poles);
    for options in [
        Options {
            maximum_work: WORK + 1,
            ..defaults
        },
        Options {
            maximum_nodes: NODES + 1,
            ..defaults
        },
        Options {
            maximum_bytes: BYTES + 1,
            ..defaults
        },
        Options {
            maximum_depth: 33,
            ..defaults
        },
    ] {
        assert!(matches!(
            certify(
                Input {
                    curves: &original,
                    counter_clockwise: true
                },
                options,
                None
            ),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        certify(
            Input {
                curves: &[],
                counter_clockwise: true
            },
            defaults,
            None
        ),
        Err(Error::Invalid)
    ));
    let excessive = vec![original[0]; CURVES + 1];
    assert!(matches!(
        certify(
            Input {
                curves: &excessive,
                counter_clockwise: true
            },
            defaults,
            None
        ),
        Err(Error::Invalid)
    ));
    for case in 0..6 {
        let mut changed = original.clone();
        match case {
            0 => changed[1].curve_owner = changed[0].curve_owner,
            1 => changed[1].start_owner = changed[0].start_owner,
            2 => changed[0].end_owner = 99,
            3 => changed[0].spans = 0,
            4 => changed[0].spans = usize::MAX,
            _ => changed[0].poles = &[],
        }
        assert!(
            matches!(
                certify(
                    Input {
                        curves: &changed,
                        counter_clockwise: true
                    },
                    defaults,
                    None
                ),
                Err(Error::Invalid)
            ),
            "case{case}"
        );
    }
    for bad in [f64::NAN, f64::INFINITY, 1e7_f64.next_up()] {
        let mut p = poles.clone();
        p[0][3][0] = bad;
        let c = curves(&p);
        assert!(matches!(
            certify(
                Input {
                    curves: &c,
                    counter_clockwise: true
                },
                defaults,
                None
            ),
            Err(Error::Invalid)
        ));
    }
    let mut p = poles.clone();
    p[0][3][2] = 0.01;
    let c = curves(&p);
    assert!(matches!(
        certify(
            Input {
                curves: &c,
                counter_clockwise: true
            },
            defaults,
            None
        ),
        Err(Error::Invalid)
    ));
    let mut p = poles.clone();
    p[0][11][0] = 0.01;
    let c = curves(&p);
    assert!(matches!(
        certify(
            Input {
                curves: &c,
                counter_clockwise: true
            },
            defaults,
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        certify(
            Input {
                curves: &original,
                counter_clockwise: false
            },
            Options {
                maximum_depth: 0,
                ..defaults
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    let huge = vec![[0.; 3]; POINTS + 1];
    let c = [Curve {
        poles: &huge,
        ..original[0]
    }];
    assert!(matches!(
        certify(
            Input {
                curves: &c,
                counter_clockwise: true
            },
            defaults,
            None
        ),
        Err(Error::Limit)
    ));
    let part = vec![[0.; 3]; POINTS / 2 + 1];
    let c = [
        Curve {
            poles: &part,
            ..original[0]
        },
        Curve {
            poles: &part,
            ..original[1]
        },
    ];
    assert!(matches!(
        certify(
            Input {
                curves: &c,
                counter_clockwise: true
            },
            defaults,
            None
        ),
        Err(Error::Limit)
    ));
}
#[test]
fn rotated_literal_square_refines_the_closed_half_plane_boundary() {
    let vertices = [[-100., 1.], [1., 100.], [100., -1.], [-1., -100.]];
    let poles = (0..4)
        .map(|e| {
            let (a, b) = (vertices[e], vertices[(e + 1) % 4]);
            (0..12)
                .map(|i| {
                    std::array::from_fn(|k| {
                        if k == 2 {
                            0.
                        } else {
                            (a[k] * (11 - i) as f64 + b[k] * i as f64) / 11.
                        }
                    })
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let curves = curves(&poles);
    let c = certify(
        Input {
            curves: &curves,
            counter_clockwise: false,
        },
        Default::default(),
        None,
    )
    .unwrap();
    assert_eq!(c.winding(), -1);
    assert!(c.angular_advance()[0] > 6. && c.angular_advance()[1] < 7.);
}
#[test]
fn composite_literal_curve_retains_all_original_cells_and_single_owner() {
    let pieces = diamond(1, false, 2.);
    let mut poles = pieces[0].clone();
    for p in &pieces[1..] {
        poles.extend_from_slice(&p[1..]);
    }
    let curves = [Curve {
        curve_owner: 99,
        start_owner: 1,
        end_owner: 1,
        spans: 4,
        poles: &poles,
    }];
    let c = certify(
        Input {
            curves: &curves,
            counter_clockwise: true,
        },
        Default::default(),
        None,
    )
    .unwrap();
    assert_eq!(c.winding(), 1);
    assert_eq!(c.plane_z(), 2.);
}
#[test]
fn nearly_orthogonal_split_endpoints_keep_real_float_enclosures() {
    let mut axial = -98_f64;
    for _ in 0..16 {
        axial = axial.next_up();
        let vertices = [[-201., axial], [1., 100.], [201., -axial], [-1., -100.]];
        let poles = (0..4)
            .map(|e| {
                let (a, b) = (vertices[e], vertices[(e + 1) % 4]);
                (0..12)
                    .map(|i| {
                        std::array::from_fn(|k| {
                            if k == 2 {
                                0.
                            } else {
                                (a[k] * (11 - i) as f64 + b[k] * i as f64) / 11.
                            }
                        })
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let curves = curves(&poles);
        let c = certify(
            Input {
                curves: &curves,
                counter_clockwise: false,
            },
            Default::default(),
            None,
        )
        .unwrap();
        assert_eq!(c.winding(), -1);
        assert!(c.angular_advance()[0] > 6. && c.angular_advance()[1] < 7.);
    }
}
#[test]
fn positive_quotients_enclose_exact_binary_rational_division_at_extreme_scales() {
    use num_rational::BigRational as Q;
    for (numerator, denominator) in [
        (I::new(1e-320, 2e-320), I::new(8e-320, 1e-319)),
        (I::new(1e-308, 2e-308), I::new(8e-308, 1e-307)),
        (I::new(1e300, 2e300), I::new(8e300, 1e301)),
    ] {
        let ratio = positive_ratio(numerator, denominator);
        assert!(
            Q::from_float(ratio.lo).unwrap()
                <= Q::from_float(numerator.lo).unwrap() / Q::from_float(denominator.hi).unwrap()
        );
        assert!(
            Q::from_float(ratio.hi).unwrap()
                >= Q::from_float(numerator.hi).unwrap() / Q::from_float(denominator.lo).unwrap()
        );
        assert!(ratio.valid() && ratio.lo > 0. && ratio.hi < 0.5);
    }
    for clockwise in [false, true] {
        let mut poles = diamond(1, clockwise, 0.);
        for point in poles.iter_mut().flatten() {
            point[0] *= 1e-160;
            point[1] *= 1e-160;
        }
        let curves = curves(&poles);
        let proof = certify(
            Input {
                curves: &curves,
                counter_clockwise: !clockwise,
            },
            Default::default(),
            None,
        )
        .unwrap();
        assert_eq!(proof.winding(), if clockwise { -1 } else { 1 });
    }
}
