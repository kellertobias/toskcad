use super::*;
use num_rational::BigRational as Q;
fn q(v: f64) -> Q {
    Q::from_float(v).unwrap()
}
fn identity() -> [[f64; 3]; 3] {
    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
}
#[test]
fn exact_native_chain_preserves_literal_leaves_and_encloses_independent_rational_composition() {
    let rotations = [
        identity(),
        [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
        [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
    ];
    let mut p =
        LiteralOccurrencePlacement::new(rotations[0], [11., -5., 3.], Default::default(), None)
            .unwrap();
    let mut exact: [[Q; 3]; 3] = identity().map(|r| r.map(q));
    let mut translation = [11., -5., 3.].map(q);
    for r in rotations.into_iter().skip(1) {
        let t = [17., -8., 13.];
        p = p.transformed(r, t, Default::default(), None).unwrap();
        translation = std::array::from_fn(|i| {
            (0..3)
                .map(|j| q(r[i][j]) * translation[j].clone())
                .sum::<Q>()
                + q(t[i])
        });
        exact = std::array::from_fn(|i| {
            std::array::from_fn(|j| (0..3).map(|k| q(r[i][k]) * exact[k][j].clone()).sum())
        });
    }
    for (i, row) in exact.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            let b = p.linear_rows_bounds()[i][j];
            assert!(q(b[0]) <= *value && *value <= q(b[1]));
        }
    }
    let absolute = |value: Q| if value < q(0.) { -value } else { value };
    let norm = q(p.linear_l1_norm());
    let exact_norm = (0..3)
        .map(|j| (0..3).map(|i| absolute(exact[i][j].clone())).sum::<Q>())
        .max()
        .unwrap();
    assert!(norm >= exact_norm && norm <= exact_norm + q(1e-12));
    for vector in [[1., -2., 4.], [-0.3, 0.7, -0.9], [0., 0., 1.]] {
        let image: [Q; 3] =
            std::array::from_fn(|i| (0..3).map(|j| exact[i][j].clone() * q(vector[j])).sum());
        let length = vector.into_iter().map(|v| absolute(q(v))).sum::<Q>();
        assert!(image.into_iter().map(absolute).sum::<Q>() <= norm.clone() * length);
    }
    let at = [2., 3., 4.];
    let world = p.map_box_mm(at.map(|v| [v, v])).unwrap();
    for i in 0..3 {
        let wanted =
            (0..3).map(|j| exact[i][j].clone() * q(at[j])).sum::<Q>() + translation[i].clone();
        assert!(q(world[i][0]) <= wanted && wanted <= q(world[i][1]));
    }
    let back = p.inverse_box_mm(world).unwrap();
    assert!(
        back.into_iter()
            .zip(at)
            .all(|(b, x)| b[0] <= x && x <= b[1] && b[1] - b[0] < 1e-10)
    );
    assert_eq!(p.transforms().len(), 3);
    assert_eq!(p.transforms()[1].matrix(), rotations[1]);
    assert_eq!(p.transforms()[1].translation_mm(), [17., -8., 13.]);
    assert!(p.determinant_bounds()[0] > 0. && p.metric_error() < 1e-6);
    assert!(p.area_scale_bounds()[0] > 0. && p.normal_condition_upper() >= 1.);
    assert!(p.work_used() > 0 && p.retained_budget_bytes() < p.accounted_peak_bytes());
    let normal = p.normal_bounds([[0., 0.], [0., 0.], [1., 1.]]).unwrap();
    assert!(normal.iter().flatten().all(|x| x.is_finite()));
    assert!(p.vector_bounds(at.map(|v| [v, v])).is_ok());
    assert!(p.translation_bounds_mm().iter().all(|b| b[0] <= b[1]));
}
#[test]
fn exact_native_chain_refuses_real_invalid_resource_and_precision_requests() {
    for matrix in [
        [[2., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[-1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    ] {
        assert_eq!(
            LiteralOccurrencePlacement::new(matrix, [0.; 3], Default::default(), None).err(),
            Some(Error::Unsupported)
        );
    }
    assert_eq!(
        LiteralOccurrencePlacement::new(identity(), [f64::NAN, 0., 0.], Default::default(), None)
            .err(),
        Some(Error::Invalid)
    );
    assert_eq!(
        LiteralOccurrencePlacement::new(identity(), [1e8, 0., 0.], Default::default(), None).err(),
        Some(Error::Limit)
    );
    assert_eq!(
        LiteralOccurrencePlacement::new(
            identity(),
            [0.; 3],
            Options {
                maximum_work: 1,
                ..Default::default()
            },
            None
        )
        .err(),
        Some(Error::Limit)
    );
    assert_eq!(
        LiteralOccurrencePlacement::new(identity(), [0.; 3], Default::default(), Some(&|| true))
            .err(),
        Some(Error::Cancelled)
    );
    let p = LiteralOccurrencePlacement::new(identity(), [0.; 3], Default::default(), None).unwrap();
    assert_eq!(
        p.normal_bounds([[0., 0.]; 3]).err(),
        Some(Error::Unsupported)
    );
    assert_eq!(p.map_box_mm([[1., 0.]; 3]).err(), Some(Error::Invalid));
    assert_eq!(p.map_box_mm([[1e8, 1e8]; 3]).err(), Some(Error::Limit));
}
