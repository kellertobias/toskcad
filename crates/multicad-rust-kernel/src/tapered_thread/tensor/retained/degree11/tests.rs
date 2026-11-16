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

#[test]
fn first_order_only_whole_leaf_bounds_match_independent_bilinear_and_chain_scaling() {
    let source = bilinear([1, 1]);
    let root = source.cell([0, 0]).unwrap();
    for cell in [
        root.clone(),
        root.subdivide(0).unwrap()[1].subdivide(1).unwrap()[0].clone(),
    ] {
        let first = cell.first_derivative_bounds().unwrap();
        let full = cell.derivative_bounds().unwrap();
        assert_eq!(first.position, full.position);
        assert_eq!(first.du, full.du);
        assert_eq!(first.dv, full.dv);
    }
    let first = root.first_derivative_bounds().unwrap();
    contains(first.position[2], 121.);
    contains(first.du[2], 121.);
    contains(first.dv[2], 121.);
}

#[test]
fn actual_corner_coefficients_follow_independent_rational_tensor_polynomial() {
    use num_rational::BigRational as Q;
    let q = |x: f64| Q::from_float(x).unwrap();
    let source = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| {
                (0..12).map(move |j| [i as f64, j as f64, (i * i + j * j * j + i * j) as f64])
            })
            .collect(),
    )
    .unwrap();
    // E[i^2]=110u^2+11u, E[j^3]=990v^3+330v^2+11v,
    // E[ij]=121uv for Bernstein degree11. This rational polynomial is
    // independent of the implementation's interval de Casteljau subdivision.
    let mut cells = vec![source.cell([0, 0]).unwrap()];
    for axis in [0, 1, 0] {
        cells = cells
            .into_iter()
            .flat_map(|cell| cell.subdivide(axis).unwrap())
            .collect();
    }
    for cell in cells {
        let [[u0, u1], [v0, v1]] = cell.domain();
        for (actual, [u, v]) in
            cell.corner_bounds()
                .into_iter()
                .zip([[u0, v0], [u1, v0], [u1, v1], [u0, v1]])
        {
            let (u, v) = (q(u), q(v));
            let expected = [
                q(11.) * &u,
                q(11.) * &v,
                q(110.) * &u * &u
                    + q(11.) * &u
                    + q(990.) * &v * &v * &v
                    + q(330.) * &v * &v
                    + q(11.) * &v
                    + q(121.) * &u * &v,
            ];
            for (actual, expected) in actual.into_iter().zip(expected) {
                assert!(q(actual[0]) <= expected && expected <= q(actual[1]));
            }
        }
    }
}

#[test]
fn conditioned_actual_derivatives_enclose_independent_coefficients_and_normals() {
    use num_rational::BigRational as Q;
    let q = |x: f64| Q::from_float(x).unwrap();
    let source = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| {
                (0..12).map(move |j| [i as f64, j as f64, (i * i + j * j * j + i * j) as f64])
            })
            .collect(),
    )
    .unwrap();
    let root = source.cell([0, 0]).unwrap();
    for k in [-0.5, 0., 1.75] {
        let [u, r] = root.conditioned_derivative_bounds(k).unwrap();
        let [uc, rc] = conditioned_controls(&root.0.controls, k).unwrap();
        // Independent exact Bernstein coefficients after degree elevation:
        // Pu_z=20i+11+11j; Pv_z=27j(j-1)+60j+11+11i.
        for i in 0..12 {
            for j in 0usize..12 {
                let uz = q((20 * i + 11 + 11 * j) as f64);
                let vz = q((27 * j * j.saturating_sub(1) + 60 * j + 11 + 11 * i) as f64);
                let expected = [q(11.), q(0.), uz.clone()];
                let residual = [-q(k) * q(11.), q(11.), vz - q(k) * uz];
                for axis in 0..3 {
                    assert!(q(u[axis][0]) <= expected[axis] && expected[axis] <= q(u[axis][1]));
                    assert!(q(r[axis][0]) <= residual[axis] && residual[axis] <= q(r[axis][1]));
                    let a = uc[12 * i + j][axis];
                    let b = rc[12 * i + j][axis];
                    assert!(q(a.lo) <= expected[axis] && expected[axis] <= q(a.hi));
                    assert!(q(b.lo) <= residual[axis] && residual[axis] <= q(b.hi));
                }
            }
        }
        let child = root.subdivide(0).unwrap()[1].subdivide(1).unwrap()[0].clone();
        let [u, r] = child.conditioned_derivative_bounds(k).unwrap();
        let u = u.map(|[a, b]| I::new(a, b));
        let r = r.map(|[a, b]| I::new(a, b));
        let normal: [I; 3] = std::array::from_fn(|a| {
            u[(a + 1) % 3] * r[(a + 2) % 3] - u[(a + 2) % 3] * r[(a + 1) % 3]
        });
        for su in [0.5, 0.75, 1.] {
            for sv in [0., 0.25, 0.5] {
                let zu = 220. * su + 11. + 121. * sv;
                let zv = 2970. * sv * sv + 660. * sv + 11. + 121. * su;
                let expected = [-11. * zu * 0.25, -11. * zv * 0.25, 121. * 0.25];
                for axis in 0..3 {
                    contains([normal[axis].lo, normal[axis].hi], expected[axis]);
                }
            }
        }
    }
    for k in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            root.conditioned_derivative_bounds(k).unwrap_err(),
            Error::Invalid
        );
    }
    assert_eq!(
        root.conditioned_derivative_bounds(f64::MAX).unwrap_err(),
        Error::Unresolved
    );
    assert!(root.conditioned_derivative_bounds(1.75).is_ok());
}

#[test]
fn actual_bilinear_residual_and_triangle_warp_have_independent_rational_oracles() {
    use num_rational::BigRational as Q;
    let q = |x: f64| Q::from_float(x).unwrap();
    let nonlinear = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| {
                (0..12).map(move |j| [i as f64, j as f64, (i * i + j * j * j + i * j) as f64])
            })
            .collect(),
    )
    .unwrap();
    let cell = nonlinear.cell([0, 0]).unwrap();
    let bound = cell.bilinear_error_bound_mm().unwrap();
    let mut largest = 0i64;
    for i in 0i64..12 {
        for j in 0i64..12 {
            // Exact degree-elevated bilinear z coefficient is11i+121j+ij.
            let residual = (i * i + j * j * j + i * j) - (11 * i + 121 * j + i * j);
            largest = largest.max(residual.abs());
            assert!(q(bound) >= Q::from_integer(residual.abs().into()));
        }
    }
    assert_eq!(largest, 540);
    assert!(bound < 540. + 1e-8);
    let actual = bilinear([1, 1]).cell([0, 0]).unwrap();
    assert!(actual.bilinear_error_bound_mm().unwrap() < 1e-9);
    // For P=(11u,11v,121uv), either literal diagonal triangle has exact
    // signed residual -121v(1-u) or -121u(1-v). The exact maximum30.25 is
    // attained at the common midpoint, independently of the tensor helper.
    let warp = q(121.) / q(4.);
    for i in 0..=16 {
        for j in 0..=16 {
            let (u, v) = (q(i as f64) / q(16.), q(j as f64) / q(16.));
            let error = if u >= v {
                q(121.) * &v * (q(1.) - &u)
            } else {
                q(121.) * &u * (q(1.) - &v)
            };
            assert!(error >= q(0.) && error <= warp);
        }
    }
    assert_eq!(q(121.) * q(0.5) * (q(1.) - q(0.5)), warp);
    for child in actual.subdivide(0).unwrap() {
        for child in child.subdivide(1).unwrap() {
            assert!(child.bilinear_error_bound_mm().unwrap() < 1e-8);
        }
    }
}

#[test]
fn actual_endpoint_normals_include_independent_polynomial_and_leaf_scaling() {
    use num_rational::BigRational as Q;
    let q = |x: f64| Q::from_float(x).unwrap();
    let root = bilinear([1, 1]).cell([0, 0]).unwrap();
    let mut cells = vec![root.clone()];
    cells.extend(
        root.subdivide(0)
            .unwrap()
            .into_iter()
            .flat_map(|c| c.subdivide(1).unwrap()),
    );
    for cell in cells {
        let [[u0, u1], [v0, v1]] = cell.domain();
        let factor = (q(u1) - q(u0)) * (q(v1) - q(v0));
        for (normal, [u, v]) in
            cell.corner_normal_bounds()
                .into_iter()
                .zip([[u0, v0], [u1, v0], [u1, v1], [u0, v1]])
        {
            let expected = [
                -q(1331.) * q(v) * &factor,
                -q(1331.) * q(u) * &factor,
                q(121.) * &factor,
            ];
            for (actual, expected) in normal.into_iter().zip(expected) {
                assert!(q(actual[0]) <= expected && expected <= q(actual[1]));
            }
        }
    }
}
