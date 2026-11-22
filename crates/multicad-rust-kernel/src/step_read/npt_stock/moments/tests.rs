use super::*;
use num_bigint::BigInt as Z;
use num_rational::BigRational as Q;
fn choose(n: usize, k: usize) -> Z {
    (0..k.min(n - k)).fold(Z::from(1), |v, i| v * Z::from(n - i) / Z::from(i + 1))
}
fn product(n: usize, i: usize, m: usize, j: usize) -> Q {
    Q::new(
        choose(n, i) * choose(m, j),
        Z::from(n + m + 1) * choose(n + m, i + j),
    )
}
fn derivative(i: usize, j: usize) -> Q {
    let a = if j == 0 {
        Q::from_integer(Z::from(0))
    } else {
        product(11, i, 10, j - 1)
    };
    let b = if j == 11 {
        Q::from_integer(Z::from(0))
    } else {
        product(11, i, 10, j)
    };
    Q::from_integer(Z::from(11)) * (a - b)
}
fn encloses(bound: I, exact: Q) {
    assert!(Q::from_float(bound.lo).unwrap() <= exact && exact <= Q::from_float(bound.hi).unwrap());
}
#[test]
fn every_fixed_weight_matches_independent_exact_bernstein_product_integrals() {
    let mut budget = Budget::new(Default::default(), None).unwrap();
    let matrix = Moments::new(&mut budget).unwrap();
    for i in 0..12 {
        for j in 0..12 {
            encloses(matrix.pair[i][j], product(11, i, 11, j));
            encloses(matrix.skew[i][j], derivative(i, j) - derivative(j, i));
        }
    }
    let curve = std::array::from_fn(|i| {
        [
            I::point(11. + i as f64),
            I::point((i * i) as f64),
            I::point(0.),
        ]
    });
    encloses(
        matrix.curve_angular(&curve, &mut budget).unwrap(),
        Q::new(Z::from(5203), Z::from(3)),
    );
    let tensor = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            [
                I::point(11. + i as f64),
                I::point((i * j) as f64),
                I::point(0.),
            ]
        })
    });
    encloses(
        matrix.tensor_angular(&tensor, &mut budget).unwrap(),
        Q::new(Z::from(1331), Z::from(2)),
    );
}
#[test]
fn moment_primitives_charge_real_creation_and_query_work_and_cancellation() {
    let options = Options::default();
    let mut budget = Budget::new(
        Options {
            maximum_work: 4095,
            ..options
        },
        None,
    )
    .unwrap();
    assert!(matches!(Moments::new(&mut budget), Err(Error::Limit)));
    let mut budget = Budget::new(options, None).unwrap();
    let matrix = Moments::new(&mut budget).unwrap();
    let controls = [[[I::point(1.); 3]; 12]; 12];
    let mut budget = Budget::new(
        Options {
            maximum_work: 8191,
            ..options
        },
        None,
    )
    .unwrap();
    assert!(matches!(
        matrix.tensor_angular(&controls, &mut budget),
        Err(Error::Limit)
    ));
    let mut budget = Budget::new(
        Options {
            maximum_work: 1023,
            ..options
        },
        None,
    )
    .unwrap();
    assert!(matches!(
        matrix.curve_angular(&controls[0], &mut budget),
        Err(Error::Limit)
    ));
    let mut budget = Budget::new(options, Some(&|| true)).unwrap();
    assert!(matches!(
        matrix.tensor_angular(&controls, &mut budget),
        Err(Error::Cancelled)
    ));
}
