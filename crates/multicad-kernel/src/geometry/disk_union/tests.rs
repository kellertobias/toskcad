use super::*;
use std::f64::consts::{FRAC_PI_4, PI};

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-7 * expected.abs().max(1.),
        "{actual} != {expected}"
    );
}

fn rod(center: [f64; 2], radius: f64, z: f64, height: f64) -> Solid {
    Solid::cylinder(
        PlaneFrame {
            origin: [center[0], center[1], z],
            ..Default::default()
        },
        radius,
        height,
        policy(),
    )
    .unwrap()
}

/// Area common to two disks whose centres are `d` apart.
fn lens(r1: f64, r2: f64, d: f64) -> f64 {
    let a = r1 * r1 * ((d * d + r1 * r1 - r2 * r2) / (2. * d * r1)).portable_acos();
    let b = r2 * r2 * ((d * d + r2 * r2 - r1 * r1) / (2. * d * r2)).portable_acos();
    let k = ((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2)).sqrt();
    a + b - 0.5 * k
}

fn join(a: &Solid, b: &Solid) -> Solid {
    a.boolean(b, BooleanOperation::Union, policy())
        .unwrap()
        .unwrap()
}

/// The drill guide's thumb knob: a radius-10 hub and eight radius-2 lobes
/// on a radius-9 circle, all 9 mm tall, joined one lobe at a time.
#[test]
fn a_hub_takes_its_lobes_one_at_a_time_as_one_exact_profile() {
    let mut knob = rod([12., 30.], 10., 20., 9.);
    for n in 0..8 {
        let (sin, cos) = (n as f64 * FRAC_PI_4).portable_sin_cos();
        let lobe = rod([12. + 9. * cos, 30. + 9. * sin], 2., 20., 9.);
        knob = if n % 2 == 0 {
            join(&knob, &lobe)
        } else {
            join(&lobe, &knob)
        };
        assert!(knob.normal_sweep_body().is_some());
    }
    let outside = PI * 4. - lens(10., 2., 9.);
    near(knob.volume_mm3().unwrap(), 9. * (100. * PI + 8. * outside));
    let bounds = knob.bounds().unwrap();
    near(bounds.max[0], 12. + 11.);
    near(bounds.min[2], 20.);
    near(bounds.max[2], 29.);
    let (_, _, _, rings) = knob.normal_sweep_section().unwrap();
    assert_eq!(rings.len(), 1);
    assert_eq!(rings[0].len(), 16);
}

#[test]
fn a_ring_of_disks_encloses_a_hole() {
    let mut ring = rod([3., 0.], 2., 0., 5.);
    for n in 1..6 {
        let (sin, cos) = (n as f64 * PI / 3.).portable_sin_cos();
        ring = join(&ring, &rod([3. * cos, 3. * sin], 2., 0., 5.));
    }
    // Neighbours overlap in a lens; every other pair is apart.
    near(
        ring.volume_mm3().unwrap(),
        5. * (6. * PI * 4. - 6. * lens(2., 2., 3.)),
    );
    let (_, _, _, rings) = ring.normal_sweep_section().unwrap();
    assert_eq!(rings.len(), 2);
}

#[test]
fn tangent_apart_skew_or_offset_cylinders_are_left_alone() {
    let declined = |a: &Solid, b: &Solid| {
        a.parallel_disk_union(b, BooleanOperation::Union, policy())
            .is_none()
    };
    let hub = rod([0., 0.], 2., 0., 5.);
    // Externally and internally tangent, and two separate regions.
    assert!(declined(&hub, &rod([3., 0.], 1., 0., 5.)));
    assert!(declined(&hub, &rod([1., 0.], 1., 0., 5.)));
    assert!(declined(&hub, &rod([5., 0.], 1., 0., 5.)));
    // A different span, or a different axis.
    assert!(declined(&hub, &rod([2., 0.], 1., 0.5, 5.)));
    let skew = Solid::cylinder(
        PlaneFrame {
            origin: [2., 0., 0.],
            u: [1., 0., 0.],
            v: [0., 0., 1.],
        },
        1.,
        5.,
        policy(),
    )
    .unwrap();
    assert!(declined(&hub, &skew));
    assert!(
        hub.parallel_disk_union(
            &rod([2., 0.], 1., 0., 5.),
            BooleanOperation::Subtract,
            policy()
        )
        .is_none()
    );
    // A nested lobe leaves the hub alone; the union is the hub itself.
    let nested = hub
        .parallel_disk_union(
            &rod([0.5, 0.], 1., 0., 5.),
            BooleanOperation::Union,
            policy(),
        )
        .unwrap()
        .unwrap()
        .unwrap();
    near(nested.volume_mm3().unwrap(), 20. * PI);
    // A sweep with a straight side is not a disk union.
    let slab = Solid::box_solid([4., 4., 5.], policy()).unwrap();
    assert!(declined(&slab, &hub));
}
