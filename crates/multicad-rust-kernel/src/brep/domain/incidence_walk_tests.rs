//! Numerical first-pass contract tests; no solid or imported-source authority.
use super::*;
use crate::brep::Frame;
fn sphere() -> Surface {
    Surface::Sphere {
        frame: Frame {
            origin: [0.; 3],
            x: [1., 0., 0.],
            y: [0., 1., 0.],
            z: [0., 0., 1.],
        },
        radius: 3.,
    }
}
fn samples(surface: &Surface, parameters: &[[f64; 2]]) -> Vec<Raw> {
    parameters
        .iter()
        .enumerate()
        .map(|(i, &uv)| Raw {
            position: surface.evaluate(uv).0,
            key: Key::Vertex(i),
            segment: Segment {
                edge: i,
                t: [0., 1.],
            },
        })
        .collect()
}
#[test]
fn analytic_sphere_angles_and_pole_metadata_follow_the_actual_sample_order() {
    let surface = sphere();
    let h = std::f64::consts::FRAC_PI_2;
    let expected_uv = [[0.9, -0.4], [0.9, h], [1.3, 0.8], [1.3, -h], [0.9, -0.4]];
    let expected_poles = [false, true, false, true, false];
    let points = samples(&surface, &expected_uv);
    let mut actual = Vec::new();
    walk_incidence(&surface, &points, 1e-7, |uv| {
        let pole = surface.is_pole(uv, 1e-7);
        let seam = seams(&surface, relative(&surface, uv));
        actual.push((uv, pole, pole || seam[0] || seam[1], seam));
    })
    .unwrap();
    assert_eq!(actual.len(), expected_uv.len());
    for ((uv, pole, boundary, seam), (expected, expected_pole)) in actual
        .into_iter()
        .zip(expected_uv.into_iter().zip(expected_poles))
    {
        // In the identity sphere frame, atan2 recovers the input angles.
        // At each pole the previous sample's longitude is inherited: 0.9
        // at the north pole and 1.3 at the south pole. Neither is a seam.
        for (value, expected) in uv.into_iter().zip(expected) {
            let rounding = 8. * f64::EPSILON * (1. + expected.abs());
            assert!((value - expected).abs() <= rounding);
        }
        assert_eq!(pole, expected_pole);
        assert_eq!(boundary, expected_pole);
        assert_eq!(seam, [false, false]);
    }
}
#[test]
fn actual_off_sphere_point_refuses_before_the_visitor_is_called() {
    let surface = sphere();
    let mut points = samples(&surface, &[[0.8, 0.2], [1.2, 0.3]]);
    points[1].position = points[1].position.map(|x| x * 2.);
    let mut visited = 0;
    assert_eq!(
        walk_incidence(&surface, &points, 1e-7, |_| visited += 1),
        Err(Error::Invalid)
    );
    assert_eq!(visited, 1);
}
#[test]
fn full_unroll_retains_the_late_lune_guard_while_incidence_grants_no_domain_claim() {
    let surface = sphere();
    let h = std::f64::consts::FRAC_PI_2;
    let points = samples(
        &surface,
        &[[0.7, 0.], [0.7, h], [0.7 + 1e-10, 0.], [0.7 + 1e-10, -h]],
    );
    // Every point is genuinely on this sphere. The almost coincident
    // meridians still fail the independent full parameter-domain guard.
    walk_incidence(&surface, &points, 1e-7, |_| {}).unwrap();
    assert!(matches!(
        unroll(&surface, points, 1e-7),
        Err(Error::Invalid)
    ));
}
