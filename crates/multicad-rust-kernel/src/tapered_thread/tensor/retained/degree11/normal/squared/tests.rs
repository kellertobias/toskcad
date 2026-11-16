use super::*;
use crate::tapered_thread::{Scalar, Series};
use num_rational::BigRational as Q;
fn q(n: i64) -> Q {
    Q::from_integer(n.into())
}
fn contains(bounds: I, exact: Q) {
    assert!(
        Q::from_float(bounds.lo).unwrap() <= exact && exact <= Q::from_float(bounds.hi).unwrap(),
        "{bounds:?} misses {exact}"
    );
}
fn choose_q(n: i64, k: i64) -> Q {
    (0..k).fold(q(1), |x, i| x * q(n - i) / q(i + 1))
}
#[test]
fn actual_squared_normal_all_coefficients_enclose_independent_bilinear_polynomial() {
    let tensor = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| (0..12).map(move |j| [i as f64, j as f64, (i * j) as f64]))
            .collect(),
    )
    .unwrap();
    let source = Arc::new(NormalPolynomial11::new(&tensor.cell([0, 0]).unwrap()).unwrap());
    let squared = SquaredNormal11::new(source.clone()).unwrap();
    assert!(std::ptr::eq(squared.source(), &*source));
    drop(source);
    drop(tensor);
    for (i, row) in squared.controls.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            let ratio = |k: usize| q((k * (k.saturating_sub(1))) as i64) / q(42 * 41);
            contains(*c, q(14641) + q(1771561) * (ratio(i) + ratio(j)));
        }
    }
    assert!(squared.retained_budget_bytes() <= NORMAL11_BYTES);
    let domain = [[0.25, 0.75], [0.5, 1.]];
    let bound = squared.bounds(domain).unwrap();
    for u in [0.25, 0.75] {
        for v in [0.5, 1.] {
            contains(
                bound,
                q(14641)
                    + q(1771561) * (Q::from_float(u * u).unwrap() + Q::from_float(v * v).unwrap()),
            );
        }
    }
    for along_v in [false, true] {
        let jet = squared.taylor8(domain, along_v).unwrap();
        for t in if along_v { [0.5, 1.] } else { [0.25, 0.75] } {
            contains(jet[1], q(1771561) * Q::from_float(t).unwrap() / q(16));
        }
        contains(jet[2], q(1771561) / q(64));
        for c in jet.iter().skip(3) {
            contains(*c, q(0));
        }
    }
    for bad in [
        [[f64::NAN, 1.], [0., 1.]],
        [[0., 1.1], [0., 1.]],
        [[1., 0.], [0., 1.]],
    ] {
        assert_eq!(
            squared.bounds(bad).map(|x| [x.lo, x.hi]),
            Err(Error::Invalid)
        );
        assert_eq!(
            squared.taylor8(bad, false).map(|x| x.map(|x| [x.lo, x.hi])),
            Err(Error::Invalid)
        );
    }
    let retry = squared.bounds(domain).unwrap();
    assert_eq!([retry.lo, retry.hi], [bound.lo, bound.hi]);
}
#[test]
fn actual_curved_squared_normal_retains_polynomial_cancellation_and_sqrt_child_chain() {
    let tensor = RetainedTensor11::new(
        [1, 1],
        (0i64..12)
            .flat_map(|i| {
                (0i64..12).map(move |j| {
                    [
                        (45 * i - i * (i - 1) * (i - 2) / 6) as f64 / 1024.,
                        45. * j as f64 / 1024.,
                        (9 * i * (i - 1) / 2) as f64 / 1024.,
                    ]
                })
            })
            .collect(),
    )
    .unwrap();
    let scalar = SquaredNormal11::new(Arc::new(
        NormalPolynomial11::new(&tensor.cell([0, 0]).unwrap()).unwrap(),
    ))
    .unwrap();
    let k = q(245025) / q(1 << 20);
    for (i, row) in scalar.controls.iter().enumerate() {
        let c2 = if i >= 2 {
            choose_q(i as i64, 2) / choose_q(42, 2)
        } else {
            q(0)
        };
        let c4 = if i >= 4 {
            choose_q(i as i64, 4) / choose_q(42, 4)
        } else {
            q(0)
        };
        let expected = &k * &k * (q(1) + q(2) * c2 + c4);
        for c in row {
            contains(*c, expected.clone());
        }
    }
    let domain = [[0.25, 0.75], [0.5, 1.]];
    for along_v in [false, true] {
        let jet = Series {
            c: scalar.taylor8(domain, along_v).unwrap(),
        }
        .sqrt()
        .checked()
        .unwrap();
        for u in [Q::new(1.into(), 4.into()), Q::new(3.into(), 4.into())] {
            contains(jet.c[0], &k * (q(1) + &u * &u) / q(4));
            contains(jet.c[1], if along_v { q(0) } else { &k * &u / q(4) });
        }
        contains(jet.c[2], if along_v { q(0) } else { &k / q(16) });
        for c in jet.c.iter().skip(3) {
            contains(*c, q(0));
        }
    }
}

#[test]
fn actual_coupled_three_coordinate_squared_coefficients_match_independent_power_basis() {
    let source = RetainedTensor11::new(
        [1, 1],
        (0i64..12)
            .flat_map(|i| {
                (0i64..12).map(move |j| {
                    [
                        (i * i + j) as f64,
                        (j * j - i) as f64,
                        (i * j + 2 * j) as f64,
                    ]
                })
            })
            .collect(),
    )
    .unwrap();
    let squared = SquaredNormal11::new(Arc::new(
        NormalPolynomial11::new(&source.cell([0, 0]).unwrap()).unwrap(),
    ))
    .unwrap();
    // Independent power-basis differentiation gives these actual normal
    // polynomials. Square in that basis, then elevate directly to42x42;
    // this does not reproduce the implementation's Bernstein convolution.
    let normal: [[((usize, usize), i64); 4]; 3] = [
        [
            ((0, 0), -242),
            ((1, 0), -1331),
            ((0, 1), -1331),
            ((0, 2), -26620),
        ],
        [
            ((0, 0), -242),
            ((1, 0), -6171),
            ((2, 0), -26620),
            ((0, 1), 1331),
        ],
        [
            ((0, 0), 242),
            ((1, 0), 2420),
            ((0, 1), 2420),
            ((1, 1), 48400),
        ],
    ];
    let mut power = std::collections::BTreeMap::<(usize, usize), i64>::new();
    for p in normal {
        for ((u, v), a) in p {
            for ((x, y), b) in p {
                *power.entry((u + x, v + y)).or_default() += a * b;
            }
        }
    }
    let basis = |n: usize, k: usize| {
        if n < k {
            q(0)
        } else {
            choose_q(n as i64, k as i64) / choose_q(42, k as i64)
        }
    };
    for (i, row) in squared.controls.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            let exact = power.iter().fold(q(0), |sum, (&(u, v), &coefficient)| {
                sum + q(coefficient) * basis(i, u) * basis(j, v)
            });
            contains(*c, exact);
        }
    }
}
