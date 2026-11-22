use super::*;
#[test]
fn actual_cavity_void_and_original_outer_material_remain_distinct() {
    for index in 0..3 {
        for plain in [0., 0.2, 1.] {
            for left in [false, true] {
                let screw = super::super::tests::fixture(index, plain, left);
                assert!(screw.contains([0., 0., 4.], 0.));
                assert!(!screw.contains([0., 0., screw.head_top_z() + 0.1], 0.));
                assert!(!screw.contains([0., 0., screw.head_top_z()], 0.));
                assert!(!screw.contains([20., 0., 0.], 0.));
                assert!(!screw.contains([0., 0., 9.], 0.));
                let radius = screw.profile.top()[0];
                assert!(screw.contains([radius - 0.001, 0., screw.head_top_z() + 0.0001], 0.));
                for direction in [
                    [0.; 3],
                    [1., 0., 0.],
                    [0., 0., -1.],
                    [0., 0., 1.],
                    [0.3, -0.7, 0.5],
                ] {
                    let support = screw.support(direction);
                    assert!(support.is_finite());
                    assert_eq!(screw.support(direction.map(|v| 2. * v)), 2. * support);
                    for p in [
                        [0., 0., 4.],
                        [radius - 0.001, 0., screw.head_top_z() + 0.0001],
                    ] {
                        let projection: f64 =
                            p.into_iter().zip(direction).map(|(p, d)| p * d).sum();
                        assert!(projection <= support + 1e-12);
                    }
                }
                assert!(screw.support([f64::NAN, 0., 0.]).is_nan());
                assert!(!screw.contains([f64::NAN, 0., 0.], 0.));
                assert!(!screw.contains([0., 0., 4.], -1.));
                assert!(!screw.contains([0., 0., 4.], f64::NAN));
            }
        }
    }
}
#[test]
fn literal_translated_coordinate_enclosures_keep_boundary_ambiguity() {
    let screw = super::super::tests::fixture(1, 0.2, false);
    let cutter = screw.cavity.source();
    assert_eq!(
        classify::enclosure(cutter, [I::point(0.), I::point(0.), I::new(-1e-12, 1e-12)]).unwrap(),
        Classification::UnresolvedBoundary
    );
    assert_eq!(
        classify::enclosure(cutter, [I::point(0.), I::point(0.), I::point(-1.)]).unwrap(),
        Classification::Outside
    );
    assert!(classify::enclosure(cutter, [I::new(f64::NEG_INFINITY, f64::INFINITY); 3]).is_err());
}
