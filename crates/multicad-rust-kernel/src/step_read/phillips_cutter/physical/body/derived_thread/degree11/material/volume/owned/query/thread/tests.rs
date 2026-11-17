use super::*;
use crate::tapered_thread::RetainedTensor11;
use num_rational::BigRational as Q;
fn choose(n: usize, k: usize) -> u64 {
    (0..k).fold(1, |a, i| a * (n - i) as u64 / (i + 1) as u64)
}
fn actual_value(poles: &[Point], u: f64, v: f64) -> [Q; 3] {
    let basis = |x: f64| {
        let x = Q::from_float(x).unwrap();
        let one = Q::from_integer(1.into());
        std::array::from_fn::<_, 12, _>(|i| {
            Q::from_integer(choose(11, i).into())
                * x.pow(i as i32)
                * (&one - &x).pow((11 - i) as i32)
        })
    };
    let a = basis(u);
    let b = basis(v);
    std::array::from_fn(|k| {
        (0..12)
            .flat_map(|i| (0..12).map(move |j| (i, j)))
            .fold(Q::from_integer(0.into()), |sum, (i, j)| {
                sum + Q::from_float(poles[12 * i + j][k]).unwrap() * &a[i] * &b[j]
            })
    })
}
#[test]
fn actual_nonlinear_height_cell_hulls_exclude_coupled_rays_without_equal_v() {
    // Raw numerical tensor only, NOT an admitted Body/Source or certificate.
    // Degree-elevated coupled z=v+2u+uv, x=2+u+v²/4,y=-1/2+u+v/8.
    let poles = (0..12)
        .flat_map(|i| {
            (0..12).map(move |j| {
                let u = i as f64 / 11.;
                let v = j as f64 / 11.;
                let v2 = if j < 2 {
                    0.
                } else {
                    (j * (j - 1)) as f64 / 110.
                };
                [2. + u + 0.25 * v2, -0.5 + u + 0.125 * v, v + 2. * u + u * v]
            })
        })
        .collect::<Vec<_>>();
    let source = RetainedTensor11::new([1, 1], poles).unwrap();
    let root = source.cell([0, 0]).unwrap();
    let [low, _] = root.subdivide(0).unwrap();
    let [low, _] = low.subdivide(0).unwrap();
    let [_, high] = low.subdivide(1).unwrap();
    let [leaf, _] = high.subdivide(1).unwrap();
    let hull = leaf.bounds().unwrap().map(|b| I::new(b[0], b[1]));
    let mut values = Vec::new();
    for u in [0.125, 0.875] {
        let uv = [[u; 2], [0.625; 2]];
        let p = root.point_bounds(uv).unwrap().map(|b| I::new(b[0], b[1]));
        let exact = actual_value(source.poles(), u, 0.625);
        for k in 0..3 {
            assert!(Q::from_float(p[k].lo).unwrap() <= exact[k]);
            assert!(Q::from_float(p[k].hi).unwrap() >= exact[k]);
        }
        values.push((p, exact));
    }
    assert!(values[1].1[2] > values[0].1[2]);
    assert!(!outside_projection(values[0].0, hull));
    assert!(values[1].0[2].lo > hull[2].hi);
    assert!(outside_projection(values[1].0, hull));
}
