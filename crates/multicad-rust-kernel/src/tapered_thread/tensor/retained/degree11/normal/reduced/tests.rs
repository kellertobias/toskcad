use super::*;
use num_rational::BigRational as Q;
#[test]
fn exact_degree_elevation_preserves_independent_bernstein_monomials() {
    let w = elevation();
    let q = |x: f64| Q::from_float(x).unwrap();
    for power in 0..=7 {
        let small: [I; 8] = std::array::from_fn(|j| {
            if j < power {
                I::point(0.)
            } else {
                I::point(choose(j, power) as f64) / I::point(choose(7, power) as f64)
            }
        });
        for (i, row) in w.iter().enumerate() {
            let actual = (0..8).fold(I::point(0.), |s, j| s + row[j] * small[j]);
            let exact = if i < power {
                Q::from_integer(0.into())
            } else {
                Q::new(choose(i, power).into(), choose(21, power).into())
            };
            assert!(q(actual.lo) <= exact && exact <= q(actual.hi));
        }
    }
}
#[test]
fn auxiliary_normal_keeps_actual_owner_and_encloses_high_degree_field() {
    let tensor = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| {
                (0..12).map(move |j| {
                    [
                        i as f64,
                        j as f64,
                        if i == 11 { 1. } else { 0. } + if j == 11 { 1. } else { 0. },
                    ]
                })
            })
            .collect(),
    )
    .unwrap();
    let normal = Arc::new(NormalPolynomial11::new(&tensor.cell([0, 0]).unwrap()).unwrap());
    let address = Arc::as_ptr(&normal);
    let reduced = ReducedNormal11::new(normal).unwrap();
    assert_eq!(reduced.source() as *const _, address);
    drop(tensor);
    assert!(
        reduced
            .coefficient_residual_bounds()
            .iter()
            .any(|b| b[0] < -1e-5 || b[1] > 1e-5)
    );
    for u in [0., 0.25, 0.5, 0.75, 1.] {
        for v in [0., 0.25, 0.5, 0.75, 1.] {
            let bounds = reduced.bounds([[u, u], [v, v]]).unwrap();
            for (b, x) in bounds
                .into_iter()
                .zip([-121. * u.powi(10), -121. * v.powi(10), 121.])
            {
                assert!(b[0] <= x && x <= b[1]);
            }
        }
    }
    assert_eq!(reduced.bounds([[1., 0.], [0., 1.]]), Err(Error::Invalid));
    assert_eq!(
        reduced.bounds([[0., f64::NAN], [0., 1.]]),
        Err(Error::Invalid)
    );
    assert!(reduced.bounds([[0., 1.]; 2]).is_ok());
}

#[test]
fn reduced_cell_restrictions_keep_root_owner_residual_and_independent_normals() {
    let tensor = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| (0..12).map(move |j| [i as f64, j as f64, (i * j) as f64]))
            .collect(),
    )
    .unwrap();
    let reduced = Arc::new(
        ReducedNormal11::new(Arc::new(
            NormalPolynomial11::new(&tensor.cell([0, 0]).unwrap()).unwrap(),
        ))
        .unwrap(),
    );
    let root = reduced.root_cell();
    assert!(std::mem::size_of_val(&root) <= ReducedNormal11::CELL_BYTES);
    let residual = reduced.coefficient_residual_bounds();
    let mut cells = vec![root.clone()];
    for axis in [0, 1, 0] {
        cells = cells
            .into_iter()
            .flat_map(|c| c.subdivide(axis).unwrap())
            .collect();
    }
    for cell in cells {
        assert!(std::ptr::eq(cell.source(), reduced.as_ref()));
        assert_eq!(cell.source().coefficient_residual_bounds(), residual);
        let [[u0, u1], [v0, v1]] = cell.domain();
        let b = cell.bounds().unwrap();
        for u in [u0, (u0 + u1) / 2., u1] {
            for v in [v0, (v0 + v1) / 2., v1] {
                for (b, x) in b.into_iter().zip([-1331. * v, -1331. * u, 121.]) {
                    assert!(b[0] <= x && x <= b[1]);
                }
            }
        }
    }
    assert_eq!(root.subdivide(2).err(), Some(Error::Invalid));
    let mut deep = root.clone();
    for i in 0..64 {
        deep = deep.subdivide(i % 2).unwrap()[0].clone();
    }
    assert_eq!(deep.subdivide(0).err(), Some(Error::Limit));
    let mut right = root;
    loop {
        match right.subdivide(0) {
            Ok([_, r]) => right = r,
            Err(e) => {
                assert_eq!(e, Error::Unresolved);
                break;
            }
        }
    }
}

#[test]
fn convex_coefficient_cone_preserves_coupled_vector_direction_and_refuses_real_escape() {
    let source = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| (0..12).map(move |j| [(i * i) as f64, j as f64, (i * i + j) as f64]))
            .collect(),
    )
    .unwrap();
    let owner = Arc::new(
        ReducedNormal11::new(Arc::new(
            NormalPolynomial11::new(&source.cell([0, 0]).unwrap()).unwrap(),
        ))
        .unwrap(),
    );
    let cell = owner.root_cell();
    // This ACTUAL tilted plane has N=(121+2420u)*(-1,-1,1).
    // Its component box also contains vectors far outside the narrow cone.
    let bounds = cell.bounds().unwrap();
    let spurious = [-2541., -121., 121.];
    for (b, x) in bounds.into_iter().zip(spurious) {
        assert!(b[0] <= x && x <= b[1]);
    }
    let axis = [-1., -1., 1.];
    let cosine = spurious
        .into_iter()
        .zip(axis)
        .map(|(a, b)| a * b)
        .sum::<f64>()
        / (spurious.into_iter().map(|x| x * x).sum::<f64>().sqrt() * 3f64.sqrt());
    assert!(cosine < 0.1f64.cos());
    assert!(
        cell.certifies_normal_cone(axis.map(|x| [x, x]), 0.1)
            .unwrap()
    );
    assert!(
        !cell
            .certifies_normal_cone(axis.map(|x| [-x, -x]), 0.1)
            .unwrap()
    );
    let direction = axis.map(|x| [x - 1e-12, x + 1e-12]);
    let prepared = cell.prepare_normal_cones(direction).unwrap();
    assert!(std::mem::size_of_val(&prepared) <= 4096);
    assert!(prepared.hull_cone(0.1).unwrap());
    assert!(prepared.grouped_cones(0.1).unwrap());
    assert!(prepared.coefficient_cones(0.1).unwrap());
    for angle in [0.01, 0.1, 0.3, 0.7] {
        assert_eq!(
            prepared.coefficient_cones(angle).unwrap(),
            cell.certifies_normal_cone(direction, angle).unwrap()
        );
    }
    for angle in [0., f64::NAN, std::f64::consts::FRAC_PI_2] {
        assert_eq!(prepared.hull_cone(angle), Err(Error::Invalid));
        assert_eq!(prepared.grouped_cones(angle), Err(Error::Invalid));
        assert_eq!(prepared.coefficient_cones(angle), Err(Error::Invalid));
    }
    assert!(matches!(
        cell.prepare_normal_cones([[f64::MAX; 2]; 3]),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        cell.prepare_normal_cones([[f64::NAN; 2]; 3]),
        Err(Error::Invalid)
    ));

    assert!(
        cell.certifies_projected_normal_cone(direction, 0.1)
            .unwrap()
    );
    assert!(
        !cell
            .certifies_projected_normal_cone(axis.map(|x| [-x, -x]), 0.1)
            .unwrap()
    );
    let residual = cell.source().coefficient_residual_bounds();
    for axis in 0..2 {
        for child in cell.subdivide(axis).unwrap() {
            assert_eq!(child.source().coefficient_residual_bounds(), residual);
            assert!(
                child
                    .certifies_projected_normal_cone(direction, 0.1)
                    .unwrap()
            );
        }
    }
    assert_eq!(
        cell.certifies_projected_normal_cone([[f64::MAX; 2]; 3], 0.3),
        Err(Error::Unresolved)
    );
    assert_eq!(
        cell.certifies_projected_normal_cone([[f64::NAN; 2]; 3], 0.3),
        Err(Error::Invalid)
    );
    assert_eq!(
        cell.certifies_projected_normal_cone([[1., -1.]; 3], 0.3),
        Err(Error::Invalid)
    );
    for angle in [0., f64::NAN, std::f64::consts::FRAC_PI_2] {
        assert_eq!(
            cell.certifies_projected_normal_cone(direction, angle),
            Err(Error::Invalid)
        );
    }
    assert!(matches!(
        cell.certifies_projected_normal_cone(direction, f64::from_bits(1)),
        Ok(false) | Err(Error::Unresolved)
    ));
    assert!(
        cell.certifies_projected_normal_cone(direction, 0.1)
            .unwrap()
    );
    let saddle = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| (0..12).map(move |j| [i as f64, j as f64, (i * j) as f64]))
            .collect(),
    )
    .unwrap();
    let owner = Arc::new(
        ReducedNormal11::new(Arc::new(
            NormalPolynomial11::new(&saddle.cell([0, 0]).unwrap()).unwrap(),
        ))
        .unwrap(),
    );
    let saddle = owner.root_cell();
    let prepared = saddle
        .prepare_normal_cones([[0.; 2], [0.; 2], [1.; 2]])
        .unwrap();
    assert!(!prepared.hull_cone(0.3).unwrap());
    assert!(!prepared.grouped_cones(0.3).unwrap());
    assert!(!prepared.coefficient_cones(0.3).unwrap());
    assert_eq!(
        prepared.coefficient_cones(0.3).unwrap(),
        saddle
            .certifies_normal_cone([[0.; 2], [0.; 2], [1.; 2]], 0.3)
            .unwrap()
    );

    // Positive dot everywhere, but actual normals leave the Z cone.
    assert!(
        !saddle
            .certifies_projected_normal_cone([[0.; 2], [0.; 2], [1.; 2]], 0.3)
            .unwrap()
    );
    assert!(
        !saddle
            .certifies_normal_cone([[0.; 2], [0.; 2], [1.; 2]], 0.3)
            .unwrap()
    );
    assert_eq!(
        cell.certifies_normal_cone([[f64::NAN; 2]; 3], 0.3),
        Err(Error::Invalid)
    );
    assert_eq!(
        cell.certifies_normal_cone([[1., -1.]; 3], 0.3),
        Err(Error::Invalid)
    );
    assert!(matches!(
        cell.certifies_normal_cone(axis.map(|x| [x, x]), f64::from_bits(1)),
        Ok(false) | Err(Error::Unresolved)
    ));
    for angle in [0., f64::NAN, std::f64::consts::FRAC_PI_2] {
        assert_eq!(
            cell.certifies_normal_cone([[1.; 2]; 3], angle),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        cell.certifies_normal_cone([[f64::MAX; 2]; 3], 0.3),
        Err(Error::Unresolved)
    );
    assert!(
        cell.certifies_normal_cone(axis.map(|x| [x, x]), 0.1)
            .unwrap()
    );
}

#[test]
fn cached_actual_normal_midpoints_preserve_all_original_proposal_row_bits() {
    let tensor = RetainedTensor11::new(
        [1, 1],
        (0..12)
            .flat_map(|i| {
                (0..12).map(move |j| [i as f64, j as f64, (i * i * j + i * j * j) as f64])
            })
            .collect(),
    )
    .unwrap();
    let source = NormalPolynomial11::new(&tensor.cell([0, 0]).unwrap()).unwrap();
    assert!(
        source
            .polynomial
            .controls
            .iter()
            .flatten()
            .flatten()
            .any(|v| v.lo.to_bits() != v.hi.to_bits())
    );
    let weights = std::array::from_fn::<_, SMALL, _>(|i| basis(i as f64 / 7.));
    let cached = proposal_rows(&source, &weights);
    // Frozen original contraction order, independently recomputing every
    // directed-source midpoint at each node. This oracle checks exact bits,
    // not merely closeness or samples of the final geometry.
    let mut original = [[[0.; 3]; LARGE]; SMALL];
    for (node, values) in original.iter_mut().enumerate() {
        for (j, value) in values.iter_mut().enumerate() {
            for (i, row) in source.polynomial.controls.iter().enumerate() {
                for k in 0..3 {
                    value[k] += weights[node][i] * row[j][k].midpoint();
                }
            }
        }
    }
    for (actual, expected) in cached
        .iter()
        .flatten()
        .flatten()
        .zip(original.iter().flatten().flatten())
    {
        assert_eq!(actual.to_bits(), expected.to_bits());
    }
}
