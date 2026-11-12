use super::*;
#[test]
fn actual_rotated_corners_drive_area_rays_edges_and_centroid() {
    for angle in [0_f64, 0.37, 2.9] {
        for offset in [0., 1e-14] {
            let (s, c) = angle.sin_cos();
            let corners = headed_screw::hex_vertices(2.)
                .map(|[x, y]| [c * x - s * y + offset, s * x + c * y]);
            let drive = HexDrive::from_retained_corners(corners, 2., 1e-12).unwrap();
            assert_eq!(drive.corners_mm(), corners);
            assert_eq!(drive.is_canonical(), angle == 0. && offset == 0.);
            assert!((drive.area_mm2() - 3_f64.sqrt() * 2.).abs() < 1e-12);
            let center = drive.center_mm();
            assert!((center[0] - offset).abs() < 1e-15 && center[1].abs() < 1e-15);
            for side in 0..6 {
                assert!((drive.side_length_mm(side) - 2. / 3_f64.sqrt()).abs() < 1e-12);
            }
            assert_eq!(drive.edge_mm(usize::MAX), drive.edge_mm(usize::MAX % 6));
            assert!(
                drive.maximum_radius_mm()
                    >= corners
                        .into_iter()
                        .map(|p| p[0].hypot(p[1]))
                        .fold(0., f64::max)
            );
            let [lo, hi] = drive.area_bounds_mm2(1e-10).unwrap();
            assert!(lo < drive.area_mm2() && drive.area_mm2() < hi);
            assert!(drive.area_bounds_mm2(100.).unwrap()[0] >= 0.);
            for i in 0..72 {
                let theta = std::f64::consts::TAU * f64::from(i) / 72.;
                let r = drive.radius_at_angle_mm(theta);
                let side = drive.side_at_angle(theta).unwrap();
                let [a, b] = drive.edge_mm(side);
                let point = [r * theta.cos(), r * theta.sin()];
                assert!(cross(difference(b, a), difference(point, a)).abs() < 1e-12);
            }
            assert!(drive.radius_at_angle_mm(f64::NAN).is_nan());
            assert_eq!(drive.side_at_angle(f64::NAN), None);
            for error in [-1., f64::NAN, f64::INFINITY] {
                assert_eq!(drive.area_bounds_mm2(error), None);
            }
        }
    }
    let canonical = HexDrive::canonical(2., 1e-12).unwrap();
    assert_eq!(canonical.area_bounds_mm2(f64::MAX), None);
    assert!(canonical.is_canonical());
    assert_eq!(canonical.area_mm2(), headed_screw::hex_area(2.));
    assert_eq!(
        canonical.radius_at_angle_mm(0.4),
        headed_screw::hex_radius(2., 0.4)
    );
    assert_eq!(
        canonical.area_bounds_mm2(0.).unwrap()[0].min(canonical.area_mm2()),
        canonical.area_bounds_mm2(0.).unwrap()[0]
    );
}

#[test]
fn drive_admission_refuses_unresolved_or_unproved_shapes_without_fitting() {
    let corners = headed_screw::hex_vertices(2.);
    for af in [0., -1., f64::NAN, f64::INFINITY, 1e8] {
        assert!(HexDrive::canonical(af, 1e-12).is_err());
    }
    for error in [0., -1., f64::NAN, f64::INFINITY, 1., 1e100] {
        assert!(HexDrive::canonical(2., error).is_err());
    }
    for bad in [f64::NAN, f64::INFINITY, 1e8, 1e-3] {
        let mut changed = corners;
        changed[0][0] = bad;
        assert!(HexDrive::from_retained_corners(changed, 2., 1e-12).is_err());
    }
    let mut reverse = corners;
    reverse.reverse();
    assert!(HexDrive::from_retained_corners(reverse, 2., 1e-12).is_err());
    let shifted = corners.map(|[x, y]| [x + 8.5e-10, y]);
    assert!(HexDrive::from_retained_corners(shifted, 2., 1e-10).is_err());
    let untouched = HexDrive::canonical(2., 1e-12).unwrap();
    let mut changed = untouched.corners_mm();
    changed[0][0] += 1.;
    assert!(HexDrive::from_retained_corners(changed, 2., 1e-12).is_err());
    assert_eq!(untouched.corners_mm(), corners);
}
