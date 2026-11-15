use super::*;
fn contains(bound: [f64; 2], value: f64) {
    assert!(
        bound[0] <= value && value <= bound[1],
        "{bound:?} missing{value}"
    );
}
fn bilinear(spans: [usize; 2]) -> RetainedTensor11 {
    let [width, height] = spans.map(|n| 11 * n + 1);
    RetainedTensor11::new(
        spans,
        (0..width)
            .flat_map(|i| (0..height).map(move |j| [i as f64, j as f64, (i * j) as f64]))
            .collect(),
    )
    .unwrap()
}
#[test]
fn actual_degree11_and_independent_bilinear_values_derivatives_cover_complete_cells() {
    let source = bilinear([2, 1]);
    assert_eq!(source.spans(), [2, 1]);
    assert_eq!(source.dimensions(), [23, 12]);
    assert_eq!(source.poles().len(), 276);
    let shared = source.clone();
    assert!(std::ptr::eq(
        source.poles().as_ptr(),
        shared.poles().as_ptr()
    ));
    let cell = source.cell([1, 0]).unwrap();
    assert_eq!(cell.domain(), [[0., 1.]; 2]);
    let bounds = cell.bounds().unwrap();
    contains(bounds[0], 11.);
    contains(bounds[0], 22.);
    for (u, v) in [(0., 0.), (0.25, 0.5), (1., 1.)] {
        let (x, y) = (11. + 11. * u, 11. * v);
        let d = cell.derivative_bounds_at([[u, u], [v, v]]).unwrap();
        for (k, value) in [x, y, x * y].into_iter().enumerate() {
            contains(d.position[k], value);
            contains(cell.point_bounds([[u, u], [v, v]]).unwrap()[k], value);
        }
        contains(d.du[0], 11.);
        contains(d.dv[1], 11.);
        contains(d.du[2], 11. * y);
        contains(d.dv[2], 11. * x);
        contains(d.duv[2], 121.);
        contains(d.duu[2], 0.);
        contains(d.dvv[2], 0.);
    }
    let p = (0..12)
        .flat_map(|i| {
            (0..12).map(move |j| {
                [
                    i as f64,
                    j as f64,
                    (if i == 11 { 1. } else { 0. }) + (if j == 11 { 1. } else { 0. }),
                ]
            })
        })
        .collect();
    let high = RetainedTensor11::new([1, 1], p)
        .unwrap()
        .cell([0, 0])
        .unwrap();
    let d = high.derivative_bounds_at([[0.5; 2], [0.25; 2]]).unwrap();
    contains(d.position[2], 0.5_f64.powi(11) + 0.25_f64.powi(11));
    contains(d.du[2], 11. * 0.5_f64.powi(10));
    contains(d.dv[2], 11. * 0.25_f64.powi(10));
    contains(d.duu[2], 110. * 0.5_f64.powi(9));
    contains(d.dvv[2], 110. * 0.25_f64.powi(9));
}
#[test]
fn subdivision_preserves_actual_source_and_leaf_chain_scaling() {
    let source = bilinear([1, 1]);
    let root = source.cell([0, 0]).unwrap();
    let [left, right] = root.subdivide(0).unwrap();
    assert_eq!(left.domain(), [[0., 0.5], [0., 1.]]);
    assert_eq!(right.domain(), [[0.5, 1.], [0., 1.]]);
    let d = right.derivative_bounds_at([[0.5; 2], [0.25; 2]]).unwrap();
    contains(d.position[0], 8.25);
    contains(d.position[1], 2.75);
    contains(d.position[2], 8.25 * 2.75);
    contains(d.du[0], 5.5);
    contains(d.duv[2], 60.5);
    let child = right.subdivide(1).unwrap().into_iter().next().unwrap();
    let d = child.derivative_bounds().unwrap();
    contains(d.du[0], 5.5);
    contains(d.dv[1], 5.5);
    contains(d.duv[2], 30.25);
    assert_eq!(root.subdivide(2).unwrap_err(), Error::Invalid);
    let mut bounded = root.clone();
    for i in 0..64 {
        bounded = bounded
            .subdivide(i % 2)
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
    }
    assert_eq!(bounded.subdivide(0).unwrap_err(), Error::Limit);
    let mut rightmost = root;
    loop {
        match rightmost.subdivide(0) {
            Ok([_, right]) => rightmost = right,
            Err(error) => {
                assert_eq!(error, Error::Unresolved);
                break;
            }
        }
    }
}
#[test]
fn dimensions_numeric_inputs_and_interval_calculation_have_real_refusals() {
    for spans in [[0, 1], [1, 0]] {
        assert_eq!(
            RetainedTensor11::new(spans, vec![]).unwrap_err(),
            Error::Invalid
        );
    }
    for spans in [[usize::MAX, 1], [usize::MAX / 11 - 1, usize::MAX / 11 - 1]] {
        assert_eq!(
            RetainedTensor11::new(spans, vec![]).unwrap_err(),
            Error::Limit
        );
    }
    assert_eq!(
        RetainedTensor11::new([1, 1], vec![[0.; 3]; 143]).unwrap_err(),
        Error::Invalid
    );
    assert_eq!(
        RetainedTensor11::new([1, 4000], vec![[0.; 3]; 12 * 44001]).unwrap_err(),
        Error::Limit
    );
    for value in [f64::NAN, f64::INFINITY, 1e7_f64.next_up()] {
        let mut points = vec![[0.; 3]; 144];
        points[17][1] = value;
        assert_eq!(
            RetainedTensor11::new([1, 1], points).unwrap_err(),
            Error::Invalid
        );
    }
    let source = bilinear([1, 1]);
    assert_eq!(source.cell([1, 0]).unwrap_err(), Error::Invalid);
    assert_eq!(source.cell([0, 1]).unwrap_err(), Error::Invalid);
    let cell = source.cell([0, 0]).unwrap();
    assert_eq!(
        cell.point_bounds([[-0.1, 0.], [0., 1.]]).unwrap_err(),
        Error::Invalid
    );
    let controls = std::array::from_fn(|i| source.poles()[i].map(|v| I::new(v - 1e-4, v + 1e-4)));
    let calc = calculation_cell11(controls).unwrap();
    let d = calc.derivative_bounds_at([[0.5; 2]; 2]).unwrap();
    contains(d.position[2], 30.25);
    for value in [
        I::new(f64::NAN, 0.),
        I::new(1e7, 1e7_f64.next_up()),
        I::new(f64::NEG_INFINITY, 0.),
    ] {
        let mut invalid = controls;
        invalid[0][0] = value;
        assert_eq!(calculation_cell11(invalid).unwrap_err(), Error::Invalid);
    }
}
