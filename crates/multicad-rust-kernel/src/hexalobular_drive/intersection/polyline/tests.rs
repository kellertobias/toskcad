use super::*;
use std::cell::Cell;
fn source(scale: f64) -> Drive {
    Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8 * scale,
            minor_diameter_mm: 2.05 * scale,
            outer_lobe_radius_mm: 0.28 * scale,
        },
        Cone {
            tip_depth_mm: (1.1 + 1.4 / 3_f64.sqrt()) * scale,
            depth_per_radius: 1. / 3_f64.sqrt(),
        },
    )
    .unwrap()
}
fn distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let d: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
    let v: [f64; 3] = std::array::from_fn(|i| p[i] - a[i]);
    let dot = |a: [f64; 3], b: [f64; 3]| (0..3).map(|i| a[i] * b[i]).sum::<f64>();
    let t = (dot(d, v) / dot(d, d)).clamp(0., 1.);
    (0..3)
        .map(|i| (v[i] - t * d[i]).powi(2))
        .sum::<f64>()
        .sqrt()
}
#[test]
fn whole_finite_circle_cone_polylines_keep_contacts_and_actual_deflection() {
    for scale in [0.125, 1., 8.] {
        let drive = source(scale);
        let original = drive.clone();
        for i in 0..12 {
            let error = 0.001 * scale;
            let points = drive.intersection_polyline(i, error, 10_000).unwrap();
            let arc = drive.arcs[i];
            let contact = |[x, y]: [f64; 2]| {
                [
                    x,
                    y,
                    drive.cone.tip_depth_mm - drive.cone.depth_per_radius * x.hypot(y),
                ]
            };
            assert_eq!(points[0], contact(arc.start_mm));
            assert_eq!(*points.last().unwrap(), contact(arc.end_mm));
            for sample in 0..=1000 {
                let t = arc.start_angle_rad + arc.sweep_rad * sample as f64 / 1000.;
                let actual = contact([
                    arc.center_mm[0] + arc.radius_mm * t.cos(),
                    arc.center_mm[1] + arc.radius_mm * t.sin(),
                ]);
                let delta = points
                    .windows(2)
                    .map(|p| distance(actual, p[0], p[1]))
                    .fold(f64::INFINITY, f64::min);
                assert!(delta <= error);
            }
        }
        assert_eq!(drive, original);
    }
}
#[test]
fn malformed_precision_and_point_limits_refuse_without_partial_result_then_recover() {
    let drive = source(1.);
    for error in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            drive.intersection_polyline(0, error, 100).unwrap_err(),
            Error::Invalid
        );
    }
    for limit in [0, 1, 1_000_001] {
        assert_eq!(
            drive.intersection_polyline(0, 0.001, limit).unwrap_err(),
            Error::Limit
        );
    }
    assert_eq!(
        drive.intersection_polyline(12, 0.001, 100).unwrap_err(),
        Error::Invalid
    );
    assert_eq!(
        drive.intersection_polyline(0, 1e-30, 100).unwrap_err(),
        Error::Unresolved
    );
    // Every carrier fits; line subdivision must account for its finer chord demand.
    let carrier_count = drive
        .intersection_carrier(0, 0.001 / 4., 10_000)
        .unwrap()
        .len();
    assert_eq!(
        drive
            .intersection_polyline(0, 0.001, carrier_count + 1)
            .unwrap_err(),
        Error::Limit
    );
    assert!(drive.intersection_polyline(0, 0.001, 10_000).is_ok());
}
#[test]
fn every_source_and_subdivision_cancellation_checkpoint_is_atomic() {
    let drive = source(1.);
    let original = drive.clone();
    let count = Cell::new(0);
    let expected = drive
        .intersection_polyline_with_cancel(0, 0.001, 10_000, &|| {
            count.set(count.get() + 1);
            false
        })
        .unwrap();
    for cutoff in 1..=count.get() {
        let seen = Cell::new(0);
        assert_eq!(
            drive
                .intersection_polyline_with_cancel(0, 0.001, 10_000, &|| {
                    seen.set(seen.get() + 1);
                    seen.get() == cutoff
                })
                .unwrap_err(),
            Error::Cancelled
        );
        assert_eq!(drive, original);
    }
    assert_eq!(
        drive.intersection_polyline(0, 0.001, 10_000).unwrap(),
        expected
    );
}
