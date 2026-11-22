//! Calculation oracles only. No constructed polynomial grants Source admission.
use super::*;
use num_traits::Signed;

fn choose(n: usize, k: usize) -> BigInt {
    let mut value = BigInt::from(1);
    for i in 1..=k {
        value = value * BigInt::from(n - i + 1) / BigInt::from(i);
    }
    value
}
fn power_q(value: &BigRational, n: usize) -> BigRational {
    (0..n).fold(BigRational::from_integer(1.into()), |p, _| p * value)
}
fn contains(value: I, expected: &BigRational) {
    assert!(value.valid());
    assert!(
        BigRational::from_float(value.lo).unwrap() <= *expected
            && *expected <= BigRational::from_float(value.hi).unwrap(),
        "{value:?} misses {expected}"
    );
}
#[test]
fn every_centered_power_coefficient_encloses_independent_exact_binomial_formula() {
    for degree in 1..=11 {
        let mut controls = Vec::new();
        for i in 0..2 * (degree + 1) {
            controls.push(std::array::from_fn(|j| {
                I::point(((i * 17 + j * 11) % 31) as f64 - 15.)
            }));
        }
        let grid = Grid {
            shape: [2, degree + 1],
            values: controls,
        };
        for center in [0., 0.375, 1.] {
            let c = BigRational::from_float(center).unwrap();
            let one = BigRational::from_integer(1.into());
            for row in 0..2 {
                let result =
                    row_power(&grid, row, center, &mut Budget::new(512_000_000, None)).unwrap();
                for (k, coefficients) in result.iter().enumerate().take(degree + 1) {
                    for (axis, &bound) in coefficients.iter().enumerate() {
                        let mut exact = BigRational::from_integer(0.into());
                        for j in 0..=degree - k {
                            let mut difference = BigRational::from_integer(0.into());
                            for i in 0..=k {
                                let value = BigRational::from_float(
                                    grid.values[row * (degree + 1) + j + i][axis].lo,
                                )
                                .unwrap();
                                let sign = if (k - i) % 2 == 0 {
                                    BigInt::from(1)
                                } else {
                                    BigInt::from(-1)
                                };
                                difference +=
                                    BigRational::from_integer(sign * choose(k, i)) * value;
                            }
                            exact += BigRational::from_integer(
                                choose(degree, k) * choose(degree - k, j),
                            ) * power_q(&c, j)
                                * power_q(&(&one - &c), degree - k - j)
                                * difference;
                        }
                        contains(bound, &exact);
                    }
                }
            }
        }
    }
}
fn multiply(a: &[BigRational], b: &[BigRational]) -> Vec<BigRational> {
    let mut result = vec![BigRational::from_integer(0.into()); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            result[i + j] += x * y;
        }
    }
    result
}
#[test]
fn exact_ruled_degree70_composition_encloses_every_independent_power_coefficient() {
    let mut values = Vec::new();
    for i in 0..20 {
        values.push(std::array::from_fn(|j| {
            I::point(((i * 7 + j * 13) % 23) as f64 - 11.)
        }));
    }
    let grid = Grid {
        shape: [2, 10],
        values,
    };
    let partial = [
        grid.derivative(0, I::point(1.), &mut Budget::new(512_000_000, None))
            .unwrap(),
        grid.derivative(1, I::point(1.), &mut Budget::new(512_000_000, None))
            .unwrap(),
    ];
    let second_v = grid
        .second_v_norm(&mut Budget::new(512_000_000, None))
        .unwrap();
    let patch = Prepared {
        domain: [[I::point(0.), I::point(1.)]; 2],
        grid,
        partial,
        maximum: [1000.; 2],
        second_v,
    };
    let uv: Vec<H> = (0..8)
        .map(|i| {
            [
                I::point((i % 3) as f64 / 4.),
                I::point(((i * 3 + 1) % 7) as f64 / 8.),
                I::point(0.),
                I::point(1.),
            ]
        })
        .collect();
    let result = compose(&patch, &uv, &mut Budget::new(512_000_000, None)).unwrap();
    let to_power = |axis: usize| -> Vec<BigRational> {
        (0..8)
            .map(|k| {
                let mut coefficient = BigRational::from_integer(0.into());
                for (i, point) in uv.iter().enumerate().take(k + 1) {
                    let value = BigRational::from_float(point[axis].lo).unwrap();
                    let sign = if (k - i) % 2 == 0 {
                        BigInt::from(1)
                    } else {
                        BigInt::from(-1)
                    };
                    coefficient +=
                        BigRational::from_integer(sign * choose(7, k) * choose(k, i)) * value;
                }
                coefficient
            })
            .collect()
    };
    let u = to_power(0);
    let v = to_power(1);
    let one = BigRational::from_integer(1.into());
    let mut vp = vec![vec![one.clone()]];
    let mut cp = vec![vec![one.clone()]];
    let mut complement = v.iter().map(|x| -x).collect::<Vec<_>>();
    complement[0] += one.clone();
    for i in 1..=9 {
        vp.push(multiply(&vp[i - 1], &v));
        cp.push(multiply(&cp[i - 1], &complement));
    }
    for (axis, polynomial) in result.iter().enumerate() {
        let mut exact = vec![BigRational::from_integer(0.into()); 71];
        for j in 0..10 {
            let b = multiply(&vp[j], &cp[9 - j]);
            let a = BigRational::from_float(patch.grid.values[j][axis].lo).unwrap();
            let delta = BigRational::from_float(patch.grid.values[10 + j][axis].lo).unwrap() - &a;
            let mut coefficient = u.iter().map(|x| x * &delta).collect::<Vec<_>>();
            coefficient[0] += a;
            for (k, value) in multiply(&b, &coefficient).into_iter().enumerate() {
                exact[k] += BigRational::from_integer(choose(9, j)) * value;
            }
        }
        for k in 0..=70 {
            let mut bernstein = BigRational::from_integer(0.into());
            for (i, value) in exact.iter().enumerate().take(k + 1) {
                bernstein += value * BigRational::new(choose(k, i), choose(70, i));
            }
            contains(polynomial.c[k], &bernstein);
        }
    }
}

#[test]
fn auxiliary_residual_covers_entire_restricted_box_even_outside_original_patch() {
    let base = [[0_i64, 2, -3, 5], [2, -1, 4, 7]];
    let mut values = Vec::new();
    for row in &base {
        for k in 0_usize..10 {
            let mut point = [I::point(0.); 3];
            for (axis, slot) in point.iter_mut().enumerate() {
                let mut q = BigRational::from_integer(0.into());
                let first = k.saturating_sub(6);
                for (i, &b) in row.iter().enumerate().take(k.min(3) + 1).skip(first) {
                    q += BigRational::new(choose(3, i) * choose(6, k - i), choose(9, k))
                        * BigRational::from_integer((b + axis as i64).into());
                }
                let f = q.to_f64().unwrap();
                *slot = I::new(f.next_down(), f.next_up());
            }
            values.push(point);
        }
    }
    let original = Grid {
        shape: [2, 10],
        values,
    };
    let mut work = Budget::new(512_000_000, None);
    let restricted = restricted_rows(&original, [-0.125, 1.125], &mut work).unwrap();
    let auxiliary = proposal(&restricted, 3, &mut work).unwrap();
    let error = elevation_error(&restricted, &auxiliary, &mut work).unwrap();
    assert!(error < 1e-10);
    // Independent exact cubic evaluation includes both endpoints OUTSIDE
    // the old patch. No old-domain convexity is used for this check.
    for (row, row_base) in base.iter().enumerate() {
        for numerator in 0..=16 {
            let r = BigRational::new(numerator.into(), 16.into());
            let t =
                BigRational::new((-1).into(), 8.into()) + BigRational::new(5.into(), 4.into()) * &r;
            for axis in 0..3 {
                let mut actual = BigRational::from_integer(0.into());
                let mut proxy = actual.clone();
                for (i, &b) in row_base.iter().enumerate() {
                    actual += BigRational::from_integer(choose(3, i) * (b + axis as i64))
                        * power_q(&t, i)
                        * power_q(&(BigRational::from_integer(1.into()) - &t), 3 - i);
                    proxy += BigRational::from_float(auxiliary.values[row * 4 + i][axis].lo)
                        .unwrap()
                        * BigRational::from_integer(choose(3, i))
                        * power_q(&r, i)
                        * power_q(&(BigRational::from_integer(1.into()) - &r), 3 - i);
                }
                let delta = (&actual - &proxy).abs();
                assert!(delta <= BigRational::from_float(error).unwrap());
            }
        }
    }
    let mut wrong = auxiliary;
    wrong.values[1][0] = wrong.values[1][0] + I::point(0.01);
    assert!(elevation_error(&restricted, &wrong, &mut work).unwrap() > 1e-4);
}

#[test]
fn entire_degree33_target_composition_encloses_independent_exact_substitution() {
    let mut values = Vec::new();
    for i in 0..12 {
        values.push(std::array::from_fn(|j| {
            I::point(((i * 7 + j * 3) % 19) as f64 - 9.)
        }));
    }
    let grid = Grid {
        shape: [1, 12],
        values,
    };
    let mut work = Budget::new(512_000_000, None);
    let derivative = grid.derivative(1, I::point(1.), &mut work).unwrap();
    let second_v = grid.second_v_norm(&mut work).unwrap();
    let target = TargetPatch {
        domain: [I::point(0.), I::point(1.)],
        grid,
        derivative,
        maximum: 1000.,
        second_v,
    };
    let phi = [
        [I::point(0.), I::point(0.), I::point(0.), I::point(1.)],
        [I::point(0.), I::point(0.), I::point(0.), I::point(1.)],
        [I::point(0.), I::point(0.), I::point(0.), I::point(1.)],
        [I::point(1.); 4],
    ];
    let result = compose_target(&target, &phi, &mut work).unwrap();
    // phi(t)=t^3. First convert the ORIGINAL degree11 Bernstein vector
    // directly to exact power coefficients, then substitute and convert
    // the sparse degree33 power vector back to Bernstein independently.
    for (axis, polynomial) in result.iter().enumerate() {
        let mut coefficients = vec![BigRational::from_integer(0.into()); 34];
        for k in 0..12 {
            for i in 0..=k {
                let sign = if (k - i) % 2 == 0 {
                    BigInt::from(1)
                } else {
                    BigInt::from(-1)
                };
                coefficients[3 * k] +=
                    BigRational::from_integer(sign * choose(11, k) * choose(k, i))
                        * BigRational::from_float(target.grid.values[i][axis].lo).unwrap();
            }
        }
        for k in 0..34 {
            let mut exact = BigRational::from_integer(0.into());
            for (i, value) in coefficients.iter().enumerate().take(k + 1) {
                exact += value * BigRational::new(choose(k, i), choose(33, i));
            }
            contains(polynomial.c[k], &exact);
        }
    }
}

#[test]
fn whole_vector_chord_bound_covers_extended_interval_and_metered_failure_retry() {
    let mut values = Vec::new();
    for row in 1..=2 {
        for i in 0..3 {
            values.push(std::array::from_fn(|axis| {
                I::point(if i == 2 {
                    (row * (axis + 1)) as f64
                } else {
                    0.
                })
            }));
        }
    }
    let grid = Grid {
        shape: [2, 3],
        values,
    };
    let range = [-0.125, 1.125];
    let second = grid
        .second_v_norm(&mut Budget::new(512_000_000, None))
        .unwrap();
    let mut work = Budget::new(512_000_000, None);
    let (model, error) = chord_model(&grid, range, second, 40., &mut work)
        .unwrap()
        .unwrap();
    let used = work.used();
    let limit = BigRational::from_float(error).unwrap();
    for row in 0..2 {
        for numerator in 0..=16 {
            let r = BigRational::new(numerator.into(), 16.into());
            let t =
                BigRational::new((-1).into(), 8.into()) + BigRational::new(5.into(), 4.into()) * &r;
            let mut squared = BigRational::from_integer(0.into());
            for axis in 0..3 {
                let actual = BigRational::from_integer(((row + 1) * (axis + 1)).into()) * &t * &t;
                let a = BigRational::from_float(model.values[2 * row][axis].lo).unwrap();
                let b = BigRational::from_float(model.values[2 * row + 1][axis].lo).unwrap();
                let chord = (&BigRational::from_integer(1.into()) - &r) * a + &r * b;
                squared += (&actual - &chord) * (&actual - &chord);
            }
            assert!(squared <= &limit * &limit);
        }
    }
    assert!(
        chord_model(
            &grid,
            range,
            second,
            error * 4.,
            &mut Budget::new(512_000_000, None)
        )
        .unwrap()
        .is_none()
    );
    assert!(matches!(
        chord_model(&grid, range, second, 40., &mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    let retry = chord_model(&grid, range, second, 40., &mut Budget::new(used, None))
        .unwrap()
        .unwrap();
    assert_eq!(retry.1.to_bits(), error.to_bits());
    let calls = std::cell::Cell::new(0_usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    chord_model(
        &grid,
        range,
        second,
        40.,
        &mut Budget::new(used, Some(&count)),
    )
    .unwrap()
    .unwrap();
    let total = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == total
    };
    assert!(matches!(
        chord_model(
            &grid,
            range,
            second,
            40.,
            &mut Budget::new(used, Some(&cancel))
        ),
        Err(ReadError::Limit)
    ));
    let retry = chord_model(&grid, range, second, 40., &mut Budget::new(used, None))
        .unwrap()
        .unwrap();
    for (a, b) in model
        .values
        .iter()
        .flatten()
        .zip(retry.0.values.iter().flatten())
    {
        assert_eq!(
            (a.lo.to_bits(), a.hi.to_bits()),
            (b.lo.to_bits(), b.hi.to_bits())
        );
    }
}

#[test]
fn original_c1_collar_uses_the_global_second_derivative_chart() {
    // Left original polynomial g(v)=100v² on[0,.5]. Its right C1
    // continuation is f(v)=25+100(v-.5)-100(v-.5)². Thus the exact
    // difference is -200(v-.5)². Both actual global second norms are200.
    // Omitting the .5 patch-span conversion would underbound this example.
    let mut values = Vec::new();
    for row in 0..2 {
        for z in [0., 0., 25.] {
            values.push([I::point(row as f64), I::point(0.), I::point(z)]);
        }
    }
    let grid = Grid {
        shape: [2, 3],
        values,
    };
    let mut work = Budget::new(512_000_000, None);
    let second_v = grid.second_v_norm(&mut work).unwrap();
    let partial = [
        grid.derivative(0, I::point(1.), &mut work).unwrap(),
        grid.derivative(1, I::point(0.5), &mut work).unwrap(),
    ];
    let patch = Prepared {
        domain: [[I::point(0.), I::point(1.)], [I::point(0.), I::point(0.5)]],
        grid,
        partial,
        maximum: [1., 100.],
        second_v,
    };
    let bounds = SurfaceBounds {
        first: [1., 100.],
        second_v: 200.,
    };
    let error = extension_error(
        &patch,
        [I::new(0.2, 0.8), I::new(0.55, 0.6)],
        bounds,
        &mut work,
    )
    .unwrap();
    let delta = BigRational::from_float(0.6).unwrap() - BigRational::new(1.into(), 2.into());
    let exact = BigRational::from_integer(200.into()) * &delta * &delta;
    assert!(BigRational::from_float(error).unwrap() >= exact);
    assert!(error < 2.0000000001);
}
