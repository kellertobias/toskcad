use super::*;

fn drive() -> Drive {
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
fn contains([lo, hi]: [f64; 2], value: f64) -> bool {
    lo <= value && value <= hi
}
fn midpoint([lo, hi]: [f64; 2]) -> f64 {
    lo.midpoint(hi)
}

#[test]
fn nominal_torx_properties_keep_actual_circular_walls_and_sloping_floor() {
    let source = drive();
    let proof = source.nominal_surface_evaluations(1e-7, 100_000).unwrap();
    assert!(contains(proof.mouth.center_mm[0], 0.));
    assert!(contains(proof.mouth.center_mm[1], 0.));
    assert!(contains(proof.mouth.center_mm[2], 0.));
    let area = I::new(proof.mouth.area_mm2[0], proof.mouth.area_mm2[1]);
    let slope = (I::point(1.) + I::point(source.cone.depth_per_radius).square()).sqrt();
    let floor = area * slope;
    assert!(proof.floor.area_mm2[0] <= floor.hi && floor.lo <= proof.floor.area_mm2[1]);
    let [lo, hi] = source.volume_bounds_mm3(1e-7, 100_000).unwrap();
    let volume = I::new(lo, hi) / area;
    assert!(proof.floor.center_mm[2][0] <= volume.hi && volume.lo <= proof.floor.center_mm[2][1]);
    assert!(midpoint(proof.floor.center_mm[2]) < source.cone.tip_depth_mm);
    assert_eq!(proof.source_error_mm, source.source_error_mm());
    for i in 0..6 {
        let a = proof.walls[i];
        let b = proof.walls[i + 6];
        assert!(a.area_mm2[0] <= b.area_mm2[1] && b.area_mm2[0] <= a.area_mm2[1]);
        for axis in 0..2 {
            assert!(a.center_mm[axis][0] <= -b.center_mm[axis][0]);
            assert!(-b.center_mm[axis][1] <= a.center_mm[axis][1]);
        }
        assert!(a.center_mm[2][0] <= b.center_mm[2][1] && b.center_mm[2][0] <= a.center_mm[2][1]);
    }
}

#[test]
fn shifted_actual_contour_has_nonzero_centers_and_independent_wall_integrals() {
    let original = drive();
    let mut arcs = *original.arcs();
    for arc in &mut arcs {
        for point in [&mut arc.center_mm, &mut arc.start_mm, &mut arc.end_mm] {
            point[0] += 0.04;
            point[1] -= 0.03;
        }
    }
    let source = Drive::from_retained(arcs, original.cone(), TransportError::ZERO).unwrap();
    let proof = source.nominal_surface_evaluations(1e-7, 100_000).unwrap();
    assert!(contains(proof.mouth.center_mm[0], 0.04));
    assert!(contains(proof.mouth.center_mm[1], -0.03));
    assert_eq!(proof.mouth.center_mm[0], proof.floor.center_mm[0]);
    assert_eq!(proof.mouth.center_mm[1], proof.floor.center_mm[1]);
    for (arc, wall) in source.arcs.iter().zip(proof.walls) {
        let mut area = 0.;
        let mut moments = [0.; 3];
        let n = 100_000;
        let dt = arc.sweep_rad.abs() / f64::from(n);
        for i in 0..n {
            let t = arc.start_angle_rad + arc.sweep_rad * (f64::from(i) + 0.5) / f64::from(n);
            let x = arc.center_mm[0] + arc.radius_mm * t.cos();
            let y = arc.center_mm[1] + arc.radius_mm * t.sin();
            let depth = source.cone.tip_depth_mm - source.cone.depth_per_radius * x.hypot(y);
            let da = arc.radius_mm * depth * dt;
            area += da;
            moments[0] += da * x;
            moments[1] += da * y;
            moments[2] += da * depth / 2.;
        }
        assert!(wall.area_mm2[0] - 1e-9 <= area && area <= wall.area_mm2[1] + 1e-9);
        for (center, moment) in wall.center_mm.into_iter().zip(moments) {
            let independent = moment / area;
            assert!(center[0] - 1e-9 <= independent && independent <= center[1] + 1e-9);
        }
    }
}

#[test]
fn property_precision_limits_cancellation_and_failed_retry_are_explicit() {
    use std::cell::Cell;
    let source = drive();
    let snapshot = source.clone();
    for width in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            source.nominal_surface_evaluations(width, 100),
            Err(Error::Invalid)
        );
    }
    for limit in [0, 1_000_001, 1] {
        assert_eq!(
            source.nominal_surface_evaluations(1e-7, limit),
            Err(Error::Limit)
        );
    }
    assert_eq!(
        source.nominal_surface_evaluations_with_cancel(1e-7, 100_000, &|| true),
        Err(Error::Cancelled)
    );
    let calls = Cell::new(0usize);
    source
        .nominal_surface_evaluations_with_cancel(1e-7, 100_000, &|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    let mut threshold = 1;
    while threshold < calls.get() {
        let seen = Cell::new(0);
        assert_eq!(
            source.nominal_surface_evaluations_with_cancel(1e-7, 100_000, &|| {
                seen.set(seen.get() + 1);
                seen.get() > threshold
            }),
            Err(Error::Cancelled)
        );
        threshold *= 2;
    }
    assert_eq!(
        source.nominal_surface_evaluations(1e-300, 100_000),
        Err(Error::Unresolved)
    );
    assert_eq!(
        evaluation(I::new(-1., 1.), [I::point(0.); 3], 1.),
        Err(Error::Unresolved)
    );
    assert_eq!(
        evaluation(I::point(1.), [I::point(f64::INFINITY); 3], 1.),
        Err(Error::Unresolved)
    );
    assert_eq!(source, snapshot);
    assert!(source.nominal_surface_evaluations(1e-7, 100_000).is_ok());
}
