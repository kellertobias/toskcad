use super::*;
use num_rational::BigRational as Q;
fn q(n: i64) -> Q {
    Q::from_integer(n.into())
}
fn contains(i: I, expected: &Q) {
    assert!(
        Q::from_float(i.lo).unwrap() <= *expected && *expected <= Q::from_float(i.hi).unwrap(),
        "{i:?} misses {expected}"
    );
}
fn cell() -> TensorCell11 {
    RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| (0..12).map(move |j| [i as f64, j as f64, (i * j) as f64]))
            .collect(),
    )
    .unwrap()
    .cell([0, 0])
    .unwrap()
}
#[test]
fn actual_d11_normal_coefficients_match_independent_exact_product_convolution() {
    let points: Vec<[f64; 3]> = (0..12)
        .flat_map(|i| {
            (0..12).map(move |j| {
                [
                    (i * i + j) as f64,
                    (j * j - i) as f64,
                    (i * j + 2 * j) as f64,
                ]
            })
        })
        .collect();
    let original = RetainedTensor11::new([1, 1], points.clone()).unwrap();
    let source = original.cell([0, 0]).unwrap();
    let normal = NormalPolynomial11::new(&source).unwrap();
    let exact: Vec<[Q; 3]> = points
        .iter()
        .map(|p| p.map(|v| Q::from_float(v).unwrap()))
        .collect();
    let mut expected = vec![std::array::from_fn::<_, 3, _>(|_| q(0)); 484];
    for i in 0..11 {
        for j in 0..12 {
            for a in 0..12 {
                for b in 0..11 {
                    let weight = Q::new(
                        (choose(10, i) * choose(11, a) * choose(11, j) * choose(10, b)).into(),
                        (choose(21, i + a) * choose(21, j + b)).into(),
                    );
                    let du = std::array::from_fn::<_, 3, _>(|k| {
                        (&exact[12 * (i + 1) + j][k] - &exact[12 * i + j][k]) * q(11)
                    });
                    let dv = std::array::from_fn::<_, 3, _>(|k| {
                        (&exact[12 * a + b + 1][k] - &exact[12 * a + b][k]) * q(11)
                    });
                    for k in 0..3 {
                        expected[22 * (i + a) + j + b][k] += &weight
                            * (&du[(k + 1) % 3] * &dv[(k + 2) % 3]
                                - &du[(k + 2) % 3] * &dv[(k + 1) % 3]);
                    }
                }
            }
        }
    }
    for (actual, expected) in normal.polynomial.controls.iter().flatten().zip(expected) {
        for k in 0..3 {
            contains(actual[k], &expected[k]);
        }
    }
    assert!(normal.retained_budget_bytes() <= NORMAL11_BYTES);
    assert_eq!(NORMAL11_WORK, 400_000);
}
#[test]
fn actual_correlated_normal_hulls_and_child_taylor_scaling_enclose_independent_geometry() {
    let normal = NormalPolynomial11::new(&cell()).unwrap();
    let domain = [[0.25, 0.75], [0.5, 1.]];
    for bounds in [
        normal.bounds(domain).unwrap(),
        normal.lipschitz_bounds(domain).unwrap(),
    ] {
        for u in [0.25, 0.5, 0.75] {
            for v in [0.5, 0.75, 1.] {
                let expected = [-1331. * v, -1331. * u, 121.];
                for k in 0..3 {
                    assert!(bounds[k][0] <= expected[k] && expected[k] <= bounds[k][1]);
                }
                let at = normal.at([[u, u], [v, v]]).unwrap();
                for k in 0..3 {
                    assert!(at[k][0] <= expected[k] && expected[k] <= at[k][1]);
                }
            }
        }
    }
    for along_v in [false, true] {
        let jet = normal.taylor8(domain, along_v).unwrap();
        let expected = if along_v {
            [-1331. / 8., 0., 0.]
        } else {
            [0., -1331. / 8., 0.]
        };
        for k in 0..3 {
            assert!(jet[1][k][0] <= expected[k] && expected[k] <= jet[1][k][1]);
        }
        for order in jet.iter().skip(2) {
            for b in order {
                assert!(b[0] <= 0. && b[1] >= 0.);
            }
        }
        assert!(jet[0][2][0] <= 121. / 4. && 121. / 4. <= jet[0][2][1]);
    }
    for bad in [
        [[f64::NAN, 1.], [0., 1.]],
        [[-0.1, 1.], [0., 1.]],
        [[0., 1.1], [0., 1.]],
        [[1., 0.], [0., 1.]],
    ] {
        assert_eq!(normal.bounds(bad), Err(Error::Invalid));
        assert_eq!(normal.at(bad), Err(Error::Invalid));
        assert_eq!(normal.lipschitz_bounds(bad), Err(Error::Invalid));
        assert_eq!(normal.taylor8(bad, false), Err(Error::Invalid));
    }
    assert!(normal.bounds(domain).is_ok());
}
