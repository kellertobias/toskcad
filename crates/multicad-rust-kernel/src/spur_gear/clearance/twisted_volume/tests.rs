use super::*;
use crate::spur_gear::clearance::{DiscSide, ResultKind, certify_herringbone, classify_spur_disc};
use crate::spur_gear::herringbone::tests::fixture;
use std::f64::consts::PI;

fn pair(helical: bool) -> [Herringbone; 2] {
    [24, 36].map(|teeth| {
        let (spec, twist) = fixture(teeth);
        if helical {
            Herringbone::new_helical(spec, twist, 1e-5, 1e7)
        } else {
            Herringbone::new(spec, twist, 1e-5, 1e7)
        }
        .unwrap()
    })
}

/// The 36-tooth gear on the +x axis turned by `phase`.
fn placements(phase: f64) -> [Placement; 2] {
    let (sin, cos) = phase.sin_cos();
    [
        Placement::default(),
        Placement {
            center_mm: [60. / 15_f64.to_radians().cos(), 0.],
            axes: [[cos, sin], [-sin, cos]],
            uncertainty_mm: 0.,
        },
    ]
}

/// Certified bounds of the common section area at `z_mm` from a grid of
/// cells, each classified by the disc that covers it.
fn grid_area(gears: [&Herringbone; 2], at: [Placement; 2], z_mm: f64, cell: f64) -> [f64; 2] {
    let placed = placed_at(gears, at, z_mm, 0.);
    let tips = gears.map(Herringbone::tip_radius_mm);
    let reach = placed[1].center_mm[0];
    let [x0, x1] = [reach - tips[1], tips[0]];
    let y1 = tips[0];
    let radius = cell * std::f64::consts::FRAC_1_SQRT_2;
    let [mut inside, mut boundary] = [0usize; 2];
    let columns = ((x1 - x0) / cell).ceil() as usize;
    let rows = (2. * y1 / cell).ceil() as usize;
    for i in 0..columns {
        for j in 0..rows {
            let point = [x0 + (i as f64 + 0.5) * cell, -y1 + (j as f64 + 0.5) * cell];
            let sides = [0, 1].map(|s| {
                let p = placed[s];
                let d = sub(point, p.center_mm);
                let [u, v] = p.axes;
                let local = [d[0] * u[0] + d[1] * u[1], d[0] * v[0] + d[1] * v[1]];
                classify_spur_disc(gears[s].profile(), local, radius, 1e-7).unwrap()
            });
            if sides.contains(&DiscSide::Outside) {
                continue;
            }
            if sides == [DiscSide::Inside; 2] {
                inside += 1;
            } else {
                boundary += 1;
            }
        }
    }
    let unit = cell * cell;
    [inside as f64 * unit, (inside + boundary) as f64 * unit]
}

#[test]
fn misphased_twisted_pairs_have_the_volume_of_their_exact_common_sections() {
    for helical in [true, false] {
        let gears = pair(helical);
        let gears = [&gears[0], &gears[1]];
        // Tooth tip against tooth tip on the centre line.
        let at = placements(PI);
        assert_eq!(
            certify_herringbone(gears, at, 1e-7, None).unwrap(),
            ResultKind::CommonInterior
        );
        let found = twisted_interference_volume(gears, at, 1e-7, 1e-7, None)
            .unwrap()
            .unwrap();
        let [low, high] = found.enclosure_mm3;
        assert!(found.error_mm3 <= 1e-7);
        assert!(
            low <= found.volume_mm3
                && found.volume_mm3 <= high
                && high - low <= 0.01 * high
                && low > 1.,
            "{found:?}"
        );
        // A section area lies inside its certified cell bounds.
        let z = if helical { 3.41 } else { 0.37 };
        let exact = area(gears, at, z, 1e-7, None).unwrap().unwrap();
        let [a, b] = grid_area(gears, at, z, 0.05);
        assert!(
            a <= exact && exact <= b && b - a < 0.2 * exact,
            "{a} {exact} {b}"
        );
        // A dense midpoint rule over the same exact sections.
        let n = 200;
        let h = 4. / f64::from(n);
        let midpoint: f64 = (0..n)
            .map(|k| {
                area(gears, at, (f64::from(k) + 0.5) * h, 1e-7, None)
                    .unwrap()
                    .unwrap()
                    * h
            })
            .sum();
        assert!(
            (midpoint - found.volume_mm3).abs() < 1e-4 * found.volume_mm3,
            "{midpoint} {found:?}"
        );
        // Either order gives the same common material.
        let swapped =
            twisted_interference_volume([gears[1], gears[0]], [at[1], at[0]], 1e-7, 1e-7, None)
                .unwrap()
                .unwrap();
        assert!((swapped.volume_mm3 - found.volume_mm3).abs() < 1e-6);
    }
}

#[test]
fn mating_pose_has_no_common_volume_and_invalid_pairs_are_refused() {
    let gears = pair(true);
    let gears = [&gears[0], &gears[1]];
    let at = placements(PI * (1. + 1. / 36.));
    let found = twisted_interference_volume(gears, at, 1e-7, 1e-7, None)
        .unwrap()
        .unwrap();
    assert_eq!(found.volume_mm3, 0.);
    assert_eq!(found.enclosure_mm3[0], 0.);
    let (spec, twist) = fixture(24);
    let gapped = Herringbone::new_gapped(spec, twist, 1., 1e-5, 1e7).unwrap();
    for (gears, volume) in [([&gapped, gears[1]], 1e-7), (gears, 0.), (gears, f64::NAN)] {
        assert_eq!(
            twisted_interference_volume(gears, at, 1e-7, volume, None),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        twisted_interference_volume(gears, placements(PI), 1e-7, 1e-7, Some(&|| true)),
        Err(Error::Cancelled)
    );
}

#[test]
fn sections_through_a_thin_flank_sliver_follow_their_crossings() {
    // The 20/30 normal-system pair of the kernel gear tests, tip against
    // tip. Near 5.7 mm the drive passes a flank contact and the common
    // region gains a sliver thinner than the disc proof's chord tube.
    use crate::spur_gear::Spec;
    let normal = |teeth: u32, beta: f64| {
        let module = 2.;
        let pressure = (20_f64.to_radians().tan() / beta.cos()).atan();
        let pitch = module * f64::from(teeth) / (2. * beta.cos());
        let spec = Spec {
            teeth,
            pitch_radius_mm: pitch,
            base_radius_mm: pitch * pressure.cos(),
            tip_radius_mm: pitch + module,
            root_radius_mm: pitch - 1.25 * module,
            root_fillet_radius_mm: 0.2 * module,
            half_base_angle_rad: PI / (2. * f64::from(teeth)) + pressure.tan()
                - pressure
                - 0.05 / (2. * pitch),
            face_width_mm: 8.,
            bore_radius_mm: 2.,
        };
        (
            Herringbone::new_helical(spec, beta.tan() / pitch, 1e-5, 1e7).unwrap(),
            pitch,
        )
    };
    let beta = 20_f64.to_radians();
    let ((a, ra), (b, rb)) = (normal(20, beta), normal(30, -beta));
    let at = [
        Placement::default(),
        Placement {
            center_mm: [ra + rb, 0.],
            axes: [[-1., 0.], [0., -1.]],
            uncertainty_mm: 0.,
        },
    ];
    let gears = [&a, &b];
    let areas: Vec<f64> = (0..=40)
        .map(|k| {
            let z = 5.68 + 0.001 * f64::from(k);
            area(gears, at, z, 1e-7, None).unwrap().unwrap()
        })
        .collect();
    // The area is continuously differentiable through the contact.
    for w in areas.windows(3) {
        let bend = (w[0] - 2. * w[1] + w[2]).abs();
        assert!(bend < 1e-4, "{areas:?}");
        assert!(w[1] < w[0], "{areas:?}");
    }
    let found = twisted_interference_volume(gears, at, 1e-7, 4e-7, None)
        .unwrap()
        .unwrap();
    let [low, high] = found.enclosure_mm3;
    assert!(low <= found.volume_mm3 && found.volume_mm3 <= high && low > 100.);
}
