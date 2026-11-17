//! Exact ORIGINAL decimal tensor/Jacobian oracles, not endpoint enum mirrors.
use super::*;
use num_bigint::BigInt;
use num_rational::BigRational as Q;
use num_traits::{One, Zero};
fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn decimal(s: &str) -> Q {
    let (m, e) = s.split_once('E').unwrap_or((s, "0"));
    let exponent = e.parse::<i32>().unwrap();
    let fraction = m.split_once('.').map_or(0, |(_, v)| v.len()) as i32;
    let n = m.replace('.', "").parse::<BigInt>().unwrap();
    let power = exponent - fraction;
    if power >= 0 {
        Q::from_integer(n * BigInt::from(10).pow(power as u32))
    } else {
        Q::new(n, BigInt::from(10).pow((-power) as u32))
    }
}
fn contains(v: I, x: &Q) {
    assert!(q(v.lo) <= *x && *x <= q(v.hi), "{v:?} excludes {x}");
}
fn weights(n: usize, t: Q) -> Vec<Q> {
    (0..=n)
        .map(|i| {
            let choose = (0..i.min(n - i)).fold(1u64, |s, j| s * (n - j) as u64 / (j + 1) as u64);
            Q::from_integer(choose.into()) * t.pow(i as i32) * (Q::one() - &t).pow((n - i) as i32)
        })
        .collect()
}
fn full<const N: usize>(p: &[[I; 3]], uv: [f64; 2]) -> [I; 3] {
    let weights = uv.map(basis::<N>);
    let anchor = p[0];
    let mut value = anchor;
    for i in 0..N {
        for j in 0..N {
            let w = weights[0][i] * weights[1][j];
            for k in 0..3 {
                value[k] = value[k] + w * (p[i * N + j][k] - anchor[k]);
            }
        }
    }
    value
}
fn budget(remaining: usize) -> Budget<'static> {
    Budget {
        used: 0,
        remaining,
        bytes: SCRATCH,
        peak: SCRATCH,
        limit: BYTES,
        interrupt: None,
    }
}
#[test]
fn actual_endpoint_restrictions_match_full_original_scalar_tensor_and_exact_jacobian() {
    for quarter in [false, true] {
        let source = super::tests::source(quarter, false);
        for tensor_index in [0, 5] {
            let tensor = &source.tensors()[tensor_index];
            for root in [[0, 0], [tensor.spans()[0] - 1, tensor.spans()[1] - 1]] {
                let height = 11 * tensor.spans()[1] + 1;
                let p = std::array::from_fn::<_, 144, _>(|i| {
                    source
                        .point(
                            tensor.point_indices()
                                [(11 * root[0] + i / 12) * height + 11 * root[1] + i % 12],
                        )
                        .unwrap()
                        .coordinates()
                        .map(|s| decimal(s.spelling()))
                });
                let pi = std::array::from_fn::<_, 144, _>(|i| {
                    source
                        .point(
                            tensor.point_indices()
                                [(11 * root[0] + i / 12) * height + 11 * root[1] + i % 12],
                        )
                        .unwrap()
                        .coordinates()
                        .map(|s| interval(s.bounds()))
                });
                let normal =
                    normal::Normal::new(&source, tensor, root[0], root[1], &mut budget(16_000_000))
                        .unwrap();
                let mut coordinates = Vec::new();
                for t in [0., 0.125, 0.5, 0.875, 1.] {
                    coordinates.extend([[0., t], [1., t], [t, 0.], [t, 1.]]);
                }
                coordinates.extend([
                    [0.375, 0.625],
                    [f64::from_bits(1), 0.5],
                    [(1f64).next_down(), 0.5],
                ]);
                for uv in coordinates {
                    let w = uv.map(|t| weights(11, q(t)));
                    let d = uv.map(|t| weights(10, q(t)));
                    let mut position: [Q; 3] = std::array::from_fn(|_| Q::zero());
                    let mut du = position.clone();
                    let mut dv = position.clone();
                    for i in 0..12 {
                        for j in 0..12 {
                            for k in 0..3 {
                                position[k] += &p[i * 12 + j][k] * &w[0][i] * &w[1][j];
                                if i < 11 {
                                    du[k] += Q::from_integer(11.into())
                                        * (&p[(i + 1) * 12 + j][k] - &p[i * 12 + j][k])
                                        * &d[0][i]
                                        * &w[1][j];
                                }
                                if j < 11 {
                                    dv[k] += Q::from_integer(11.into())
                                        * (&p[i * 12 + j + 1][k] - &p[i * 12 + j][k])
                                        * &w[0][i]
                                        * &d[1][j];
                                }
                            }
                        }
                    }
                    let exact_normal: [Q; 3] = std::array::from_fn(|i| {
                        &dv[(i + 1) % 3] * &du[(i + 2) % 3] - &dv[(i + 2) % 3] * &du[(i + 1) % 3]
                    });
                    let mut bp = budget(16_000_000);
                    let point =
                        angular::registry::root_point(&source, tensor_index, root, uv, &mut bp)
                            .unwrap();
                    let mut bn = budget(16_000_000);
                    let direction =
                        normal::point_values(normal.coefficients(), uv, &mut bn).unwrap();
                    let full_point = full::<12>(&pi, uv);
                    let full_normal = full::<22>(normal.coefficients(), uv);
                    for k in 0..3 {
                        contains(point[k], &position[k]);
                        contains(full_point[k], &position[k]);
                        contains(direction[k], &exact_normal[k]);
                        contains(full_normal[k], &exact_normal[k]);
                    }
                    let ends = uv.into_iter().filter(|&t| t == 0. || t == 1.).count();
                    assert_eq!(
                        bp.used,
                        match ends {
                            2 => 32,
                            1 => 1056,
                            _ => 4096,
                        }
                    );
                    assert_eq!(
                        bn.used,
                        match ends {
                            2 => 32,
                            1 => 2048,
                            _ => 12288,
                        }
                    );
                    let mut narrow = budget(bp.used - 1);
                    assert!(matches!(
                        angular::registry::root_point(&source, tensor_index, root, uv, &mut narrow),
                        Err(Error::Limit)
                    ));
                    let mut narrow = budget(bn.used - 1);
                    assert!(matches!(
                        normal::point_values(normal.coefficients(), uv, &mut narrow),
                        Err(Error::Limit)
                    ));
                }
            }
        }
    }
}
