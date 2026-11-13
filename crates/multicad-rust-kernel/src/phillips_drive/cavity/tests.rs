use super::*;

pub(super) fn cavity(index: usize) -> Cavity {
    let (spec, domain) = crate::phillips_drive::test_fixture(index);
    Cavity::new(Cutter::new(spec, domain).unwrap(), 1e-5, 100_000).unwrap()
}
pub(super) fn oracle(s: Spec, depth: f64) -> f64 {
    let a = s.outer_flank_angle_deg.to_radians();
    let c = s.floor_angle_deg.to_radians();
    let r = s.mouth_radius_mm;
    let join = r * (1. - a.sin());
    let flank = (s.mouth_diameter_mm - s.g_diameter_mm) / (2. * a.tan());
    let tip = flank + s.g_diameter_mm / 2. * c.tan();
    if depth <= join {
        s.mouth_diameter_mm / 2. + r * (1. - a.sin()) / a.cos()
            - (r * r - (depth - r).powi(2)).max(0.).sqrt()
    } else if depth <= flank {
        s.mouth_diameter_mm / 2. - depth * a.tan()
    } else {
        ((tip - depth) / c.tan()).max(0.)
    }
}
#[test]
fn all_metric_cavities_retain_true_source_volume_and_depth_dependent_envelope() {
    for index in 0..6 {
        let c = cavity(index);
        let source = c.source();
        let (s, domain) = crate::phillips_drive::test_fixture(index);
        assert_eq!(source.spec(), s);
        assert_eq!(source.domain(), domain);
        let v = c.volume();
        assert!(v.volume_mm3[0] > 0. && v.volume_mm3[1] - v.volume_mm3[0] <= 1e-5);
        let e = source.envelope();
        let tip = e.tip_depth_mm_bounds()[1];
        let mut previous = f64::INFINITY;
        for i in 0..=200 {
            let d = tip * (i as f64 / 200.);
            let bound = c.radial_suffix_bound_mm(d).unwrap();
            let expected = oracle(s, d);
            assert!(
                bound + 2e-13 >= expected,
                "{index} d={d} {bound} {expected}"
            );
            assert!(bound <= previous + 1e-12);
            previous = bound;
            // Diagonal points stay outside all four corner prisms. The
            // independent classifier checks real H0 ellipse/V cutouts too.
            if bound > 1e-5 && d > 0. {
                let t = 0.99 * expected / 2f64.sqrt();
                assert_eq!(
                    source.classify_point_mm([t, t, d]).unwrap(),
                    Classification::Inside
                );
            }
        }
        for join in [
            e.roundover_join_depth_mm_bounds(),
            e.flank_join_depth_mm_bounds(),
        ] {
            for d in [join[0], join[1]] {
                assert!(c.radial_suffix_bound_mm(d).unwrap() >= oracle(s, d) - 2e-13);
            }
        }
        assert_eq!(c.radial_suffix_bound_mm(tip + 1.).unwrap(), 0.);
    }
}
#[test]
fn countersunk_h1_narrows_before_the_original_metric_thread_entry() {
    let c = cavity(2);
    let mouth = c.radial_suffix_bound_mm(0.).unwrap();
    let entry = c.radial_suffix_bound_mm(1.65).unwrap();
    // Actual M3 countersunk H1 cannot pass a constant mouth-radius core test.
    // Its true floor/flank suffix does fit the nominal M3 minor-radius core.
    assert!(mouth > 1.8 && entry < 0.8 && entry > 0.7);
    assert!(mouth > 1.5 && entry + 0.05 < 1.15);
}
#[test]
fn failures_cancel_and_cached_clone_preserve_actual_finite_source() {
    let c = cavity(0);
    let clone = c.clone();
    assert!(Arc::ptr_eq(&c.source, &clone.source));
    for d in [-1., f64::NAN, f64::INFINITY] {
        assert_eq!(c.radial_suffix_bound_mm(d), Err(Error::Invalid));
    }
    assert!(matches!(
        Cavity::with_cancel(c.source().clone(), 1e-5, 100_000, &|| true),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        Cavity::new(c.source().clone(), -1., 100_000),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        Cavity::new(c.source().clone(), 1e-5, 0),
        Err(Error::Limit)
    ));
    let mut domain = c.source().domain();
    domain.corner_reach_mm = 0.001;
    let source = Cutter::new(c.source().spec(), domain).unwrap();
    assert!(matches!(
        Cavity::new(source, 1e-5, 100_000),
        Err(Error::Unresolved)
    ));
    assert_eq!(c.source(), clone.source());
}
