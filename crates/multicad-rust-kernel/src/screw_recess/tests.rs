use super::*;
use crate::hexalobular_drive::{BoundarySource, TangentSpec};
fn torx() -> Drive {
    Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.1 + 1.4 / 3_f64.sqrt(),
            depth_per_radius: 1. / 3_f64.sqrt(),
        },
    )
    .unwrap()
}
#[test]
fn typed_recesses_preserve_real_floors_contours_mass_and_cheap_immutable_clones() {
    let hex = Recess::canonical_flat_hex(2., 0.6).unwrap();
    let curved = Recess::conical_hexalobular(torx(), 1e-7, 100_000).unwrap();
    assert_eq!(hex.flat_hex_geometry().unwrap().1, 0.6);
    assert!(hex.hexalobular_geometry().is_none());
    assert!(curved.flat_hex_geometry().is_none());
    assert!(matches!(
        curved.hexalobular_geometry().unwrap().boundary_source(),
        BoundarySource::AuthoredTangent(_)
    ));
    for (source, count) in [(&hex, 6), (&curved, 12)] {
        let cloned = source.clone();
        assert_eq!(&cloned, source);
        assert_eq!(source.wall_count(), count);
        assert!(source.maximum_radius_bound_mm() > 1.);
        let [lo, hi] = source.volume_bounds_mm3();
        assert!(lo <= source.nominal_volume_mm3() && source.nominal_volume_mm3() <= hi);
        assert_eq!(source.contains_void([0., 0., 0.1], 0.), Some(true));
        assert_eq!(source.contains_void([0., 0., -1.], 0.), Some(true));
        assert_eq!(source.contains_void([3., 0., 0.1], 0.), Some(false));
        assert_eq!(
            source.contains_void([0., 0., source.maximum_depth_mm()], 0.),
            Some(false)
        );
        assert_eq!(source.contains_void([0., 0., 0.1], 0.2), Some(true));
        assert!(source.radius_at_angle_mm(0.).unwrap() >= 1.);
        assert_eq!(source.radius_at_angle_mm(f64::NAN), None);
        for bad in [f64::NAN, -1., f64::INFINITY] {
            assert_eq!(source.contains_void([0., 0., 0.1], bad), None);
        }
        assert_eq!(source.contains_void([f64::NAN, 0., 0.1], 0.), None);
    }
    assert_eq!(hex.source_error_mm(), 0.);
    assert!(curved.source_error_mm() > 0.);
    let cloned = curved.clone();
    assert!(std::ptr::eq(
        curved.hexalobular_geometry().unwrap(),
        cloned.hexalobular_geometry().unwrap()
    ));
    assert_eq!(
        hex.nominal_volume_mm3(),
        crate::headed_screw::hex_area(2.) * 0.6
    );
    assert_eq!(hex.contains_void([0., 0., 1.3], 0.), Some(false));
    assert_eq!(curved.contains_void([0., 0., 1.3], 0.), Some(true));
    assert_eq!(curved.contains_void([1.25, 0., 1.3], 0.), Some(false));
}
#[test]
fn failed_recess_work_and_invalid_dimensions_leave_the_previous_source_usable() {
    for af in [0., f64::NAN, f64::INFINITY] {
        assert!(Recess::canonical_flat_hex(af, 0.6).is_err());
    }
    for depth in [0., -1., f64::NAN, f64::INFINITY, 1e8] {
        assert!(Recess::canonical_flat_hex(2., depth).is_err());
    }
    let source = torx();
    assert_eq!(
        Recess::conical_hexalobular(source.clone(), 1e-7, 0),
        Err(Error::Limit)
    );
    assert_eq!(
        Recess::conical_hexalobular(source.clone(), 1e-30, 1000),
        Err(Error::Unresolved)
    );
    assert_eq!(
        Recess::conical_hexalobular_with_cancel(source.clone(), 1e-7, 1000, &|| true),
        Err(Error::Cancelled)
    );
    assert!(Recess::conical_hexalobular(source, 1e-7, 100_000).is_ok());
}
