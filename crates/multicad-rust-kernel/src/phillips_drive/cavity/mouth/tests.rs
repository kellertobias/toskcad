use super::*;
#[test]
fn all_metric_mouths_have_actual_positive_front_cross_derivatives() {
    for index in 0..6 {
        let c = super::super::tests::cavity(index);
        let proof = c.mouth_sweep();
        let p = c.source().spec();
        let beta = p.tilt_angle_deg.to_radians();
        let theta = (p.outer_corner_angle_deg / 2.).to_radians();
        let flank = (p.mouth_diameter_mm - p.g_diameter_mm)
            / (2. * p.outer_flank_angle_deg.to_radians().tan());
        let tb =
            flank + (p.g_diameter_mm - p.b_width_mm) / 2. * p.floor_angle_deg.to_radians().tan();
        let rb = p.b_width_mm / 2. + tb * beta.tan();
        let cb = beta.cos();
        let (join, join_s) = match p.inner {
            Inner::Rounded { radius_mm: r } => (r * theta.cos(), r * (1. - theta.sin())),
            Inner::ShallowV {
                angle_deg,
                e_width_mm: e,
            } => (e / 2., e / (2. * (angle_deg / 2.).to_radians().tan())),
        };
        assert!(proof.minimum_front_tangent_intercept_mm() > 0.);
        assert_eq!(proof.maximum_radius_mm(), c.source.envelope.entry.hi);
        assert_eq!(
            proof.source_edge_count(),
            c.source
                .boundary_topology()
                .unwrap()
                .edges()
                .iter()
                .filter(|e| matches!(e.source_identity(), CurveKey::Mouth { .. }))
                .count()
        );
        for i in 0..=100 {
            let y = join * i as f64 / 100.;
            let (f, df) = match p.inner {
                Inner::ShallowV { angle_deg, .. } => {
                    let slope = 1. / (angle_deg / 2.).to_radians().tan();
                    (y * slope, slope)
                }
                Inner::Rounded { radius_mm: r } => {
                    let q = (r * r - y * y).sqrt();
                    (r - q, y / q)
                }
            };
            let numerator = rb + (f - y * df) / cb;
            assert!(numerator >= proof.minimum_front_tangent_intercept_mm());
            // The outer straight front's actual tangent intercept agrees at
            // the H0 join; no circular XY replacement is used for H0.
            let outer =
                rb + (join_s - join / (p.outer_corner_angle_deg / 2.).to_radians().tan()) / cb;
            assert!(outer >= proof.minimum_front_tangent_intercept_mm());
        }
    }
}
