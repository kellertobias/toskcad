use super::*;
#[test]
fn direct_normalization_handles_positive_subnormal_denominators_without_a_reciprocal() {
    let a = I::point(0.);
    let b = I::point(1e-310);
    let q = I::point(5e-311);
    let u = unit(q, a, b);
    assert!(u.valid() && u.lo <= 0.5 && u.hi >= 0.5);
    let u = unit(I::new(-1e-310, 2e-310), a, b);
    assert_eq!((u.lo, u.hi), (0., 1.));
}
#[test]
fn real_periodic_profile_encloses_all_bands_translates_and_coarse_numeric_domains() {
    let owner = super::super::tests::owned(false, 0.);
    let s = owner.thread();
    let nodes = s.phases();
    let global = I::new(s.radius_nodes_mm()[1], s.radius_nodes_mm()[0]);
    let ledger = || Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    for band in 0..4 {
        let mid = (nodes[band] + nodes[band + 1]) / I::point(2.);
        let nodes_r = s.radius_nodes_mm();
        let expected = (num_rational::BigRational::from_float(nodes_r[band]).unwrap()
            + num_rational::BigRational::from_float(nodes_r[band + 1]).unwrap())
            / num_bigint::BigInt::from(2);
        for k in [-100., 0., 100.] {
            let phase = mid + I::point(k) * I::point(s.pitch_mm());
            let r = bounds(s, phase, &mut ledger()).unwrap();
            assert!(
                num_rational::BigRational::from_float(r.lo).unwrap() <= expected
                    && num_rational::BigRational::from_float(r.hi).unwrap() >= expected
            );
        }
    }
    for phase in [I::new(-1., 1.), I::point(1e18), I::point(-1e308)] {
        let r = bounds(s, phase, &mut ledger()).unwrap();
        assert_eq!((r.lo, r.hi), (global.lo, global.hi));
    }
    assert!(matches!(
        bounds(s, I::point(f64::NAN), &mut ledger()),
        Err(PhillipsReadError::Invalid)
    ));
    let mut w = ledger();
    w.remaining = 0;
    assert!(matches!(
        bounds(s, I::point(0.), &mut w),
        Err(PhillipsReadError::Limit)
    ));
}
