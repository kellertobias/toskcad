use multicad_core::document::EntityId;
use multicad_kernel::{
    GeometryError, GeometryPolicy, KnuckleHelicalThread, PlaneFrame, Solid,
    rust_kernel::rounded_thread,
};
/// Owned-kernel knuckle (DIN 405) thread: tangent crest/root circles on
/// 15-degree flanks, whose minor radius is exactly the major less half a pitch.
/// Exact STEP serialization of the owned solid; the Rust analogue of the
/// former OCCT BRep text used to prove that meshing never mutates geometry.
fn step_bytes(solid: &Solid) -> Vec<u8> {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("snapshot.step");
    solid.write_step(&path).unwrap();
    std::fs::read(path).unwrap()
}

fn rust_knuckle(p: KnuckleHelicalThread, feature: EntityId) -> Result<Solid, GeometryError> {
    Solid::rust_rounded_thread(
        rounded_thread::Spec {
            major_radius_mm: p.major_radius_mm,
            pitch_mm: p.pitch_mm,
            length_mm: p.length_mm,
            depth_mm: p.pitch_mm / 2.,
            crest_radius_mm: p.crest_radius_mm,
            root_radius_mm: p.root_radius_mm,
            flank_angle_deg: 15.,
            starts: 1,
            left_hand: p.left_hand,
        },
        feature,
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
}
fn spec(internal: bool, left_hand: bool) -> KnuckleHelicalThread {
    let pitch = 3.175;
    let a = 15_f64.to_radians();
    let r = pitch * (a.cos() - a.sin()) / (4. * (1. - a.sin()));
    let delta = 0.05 * pitch * a.sin() / (1. - a.sin());
    KnuckleHelicalThread {
        major_radius_mm: 10. + if internal { 0.05 * pitch } else { 0. },
        pitch_mm: pitch,
        length_mm: 2.37 * pitch,
        crest_radius_mm: r - if internal { delta } else { 0. },
        root_radius_mm: r + if internal { delta } else { 0. },
        left_hand,
    }
}
/// Fraction of the probe cube's centre and eight corners that the owned
/// kernel's exact point classification places in material. Owned rounded
/// bodies have no exact Boolean or overlap volume (the former OCCT probe
/// intersected the cube). Each probe cube sits farther from the expected
/// boundary than its half-diagonal, so the former volume fraction was 0 or 1;
/// nine exact classifications detect the same side unless a material feature
/// smaller than the cube lies wholly inside it. This is a proxy, not an
/// exact replacement for the cube volume.
fn material(shape: &Solid, point: [f64; 3]) -> f64 {
    let edge = 0.00002;
    let samples = std::iter::once([0.; 3]).chain((0..8).map(|corner| {
        std::array::from_fn(|axis| {
            if corner >> axis & 1 == 1 {
                edge / 2.
            } else {
                -edge / 2.
            }
        })
    }));
    samples
        .filter(|offset: &[f64; 3]| {
            shape
                .rust_thread_contains(std::array::from_fn(|i| point[i] + offset[i]))
                .unwrap()
        })
        .count() as f64
        / 9.
}
fn analytic_volume(p: KnuckleHelicalThread) -> f64 {
    let a = 15_f64.to_radians();
    let rc = p.crest_radius_mm;
    let rr = p.root_radius_mm;
    let major = p.major_radius_mm;
    let minor = major - p.pitch_mm / 2.;
    let circle = |c: f64, r: f64, z: f64, sign: f64| {
        (c * c + r * r) * z - z.powi(3) / 3.
            + sign * c * (z * (r * r - z * z).sqrt() + r * r * (z / r).asin())
    };
    let zc = rc * a.cos();
    let zr = p.pitch_mm / 2. - rr * a.cos();
    let k = (1. - a.sin()) / a.cos();
    let b = 1. / a.tan();
    let aa = major + rc * k * b;
    let flank = |z: f64| aa * aa * z - aa * b * z * z + b * b * z.powi(3) / 3.;
    let moment = 2.
        * (circle(major - rc, rc, zc, 1.) + flank(zr) - flank(zc)
            + circle(minor + rr, rr, rr * a.cos(), -1.));
    std::f64::consts::PI * p.length_mm * moment / p.pitch_mm
}
#[test]
fn native_knuckle_circles_flanks_volume_and_step_keep_distinct_cutter_radii_for_both_hands() {
    let dir = tempfile::tempdir().unwrap();
    for internal in [false, true] {
        for left in [false, true] {
            let p = spec(internal, left);
            let id = EntityId::new();
            let solid = rust_knuckle(p, id).unwrap();
            assert_eq!(solid.solid_count().unwrap(), 1);
            let bounds = solid.bounds().unwrap();
            assert!(bounds.min[2].abs() < 16e-7);
            assert!((bounds.max[2] - p.length_mm).abs() < 16e-7);
            let topology = solid.topology_snapshot(id, 0).unwrap();
            assert!(topology.elements.iter().all(|e| !e.origins.is_empty()));
            // Owned rounded threads name each exact meridian band; the root
            // circle is one periodic band (OCCT split it into root/lower and
            // root/upper faces) that also carries the shared root identity.
            for role in [
                "cap.start",
                "cap.end",
                "thread/rounded/crest.start",
                "thread/rounded/lower",
                "thread/rounded/root",
                "thread/rounded/upper",
                "thread/rounded/crest.end",
            ] {
                assert!(
                    topology
                        .elements
                        .iter()
                        .any(|e| e.origins.iter().any(|o| o.source == id && o.role == role)),
                    "missing {role}"
                );
            }
            assert!(topology.elements.iter().any(|e| {
                ["thread/rounded/root", "thread.root"]
                    .iter()
                    .all(|role| e.origins.iter().any(|o| o.source == id && o.role == *role))
            }));
            let path = dir.path().join(format!("knuckle-{internal}-{left}.step"));
            solid.write_step(&path).unwrap();
            let back = Solid::read_step(&path, Default::default()).unwrap();
            let area: f64 = solid.faces().unwrap().iter().map(|f| f.area_mm2).sum();
            let exact = analytic_volume(p);
            let budget = GeometryPolicy::default().linear_tolerance_mm * area + 1e-7 * exact;
            let before = step_bytes(&solid);
            solid.mesh(Default::default()).unwrap();
            assert_eq!(step_bytes(&solid), before);
            for shape in [&solid, &back] {
                assert!(
                    (shape.volume_mm3().unwrap() - exact).abs() < budget,
                    "volume {internal}/{left}"
                );
                let a = 15_f64.to_radians();
                let hand = if left { -1. } else { 1. };
                for valley in [false, true] {
                    let radius = if valley {
                        p.root_radius_mm
                    } else {
                        p.crest_radius_mm
                    };
                    let center = if valley {
                        p.major_radius_mm - p.pitch_mm / 2. + radius
                    } else {
                        p.major_radius_mm - radius
                    };
                    for fraction in [-0.8, -0.4, 0., 0.4, 0.8] {
                        let dz = fraction * radius * a.cos();
                        let radial = (radius * radius - dz * dz).sqrt();
                        let boundary = center + if valley { -radial } else { radial };
                        let z = 1.37 * p.pitch_mm;
                        let phi = hand
                            * std::f64::consts::TAU
                            * (z - dz - if valley { p.pitch_mm / 2. } else { 0. })
                            / p.pitch_mm;
                        for (dr, inside) in [(-0.001, true), (0.001, false)] {
                            let r = boundary + dr;
                            let v = material(shape, [r * phi.cos(), r * phi.sin(), z]);
                            assert!(
                                if inside { v > 0.999 } else { v < 1e-5 },
                                "circle {internal}/{left}/{valley}/{fraction}/{dr}: {v}"
                            );
                        }
                    }
                }
                let k = (1. - a.sin()) / a.cos();
                for drop in [0.22 * p.pitch_mm, 0.28 * p.pitch_mm] {
                    for direction in [-1., 1.] {
                        let boundary =
                            2. * p.pitch_mm + direction * (drop * a.tan() + p.crest_radius_mm * k);
                        for (offset, inside) in
                            [(-direction * 0.001, true), (direction * 0.001, false)]
                        {
                            let v =
                                material(shape, [p.major_radius_mm - drop, 0., boundary + offset]);
                            assert!(
                                if inside { v > 0.999 } else { v < 1e-5 },
                                "flank {internal}/{left}/{drop}/{direction}: {v}"
                            );
                        }
                    }
                }
            }
        }
    }
}
#[test]
fn malformed_and_swapped_internal_circle_inputs_are_refused() {
    let p = spec(true, false);
    for bad in [
        KnuckleHelicalThread {
            crest_radius_mm: 0.,
            ..p
        },
        KnuckleHelicalThread {
            crest_radius_mm: p.root_radius_mm,
            root_radius_mm: p.crest_radius_mm,
            ..p
        },
        KnuckleHelicalThread {
            length_mm: p.pitch_mm / 2.,
            ..p
        },
        KnuckleHelicalThread {
            major_radius_mm: 0.1,
            ..p
        },
    ] {
        assert!(
            Solid::knuckle_helical_thread_named(bad, EntityId::new(), Default::default()).is_err()
        );
    }
}

fn mesh_material_at(mesh: &multicad_kernel::Mesh, point: [f64; 3]) -> bool {
    let [x, y, z] = point;
    let mut crossings = 0;
    for indices in mesh.indices.as_chunks::<3>().0 {
        let t: [[f64; 3]; 3] = std::array::from_fn(|i| {
            let n = indices[i] as usize * 3;
            [mesh.vertices[n], mesh.vertices[n + 1], mesh.vertices[n + 2]]
        });
        let mut points = [[0.; 2]; 2];
        let mut count = 0;
        for i in 0..3 {
            let a = t[i];
            let b = t[(i + 1) % 3];
            if (a[2] <= z && z < b[2]) || (b[2] <= z && z < a[2]) {
                let f = (z - a[2]) / (b[2] - a[2]);
                if count < 2 {
                    points[count] = [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])];
                }
                count += 1;
            }
        }
        if count == 2 {
            let [a, b] = points;
            if (a[1] > y) != (b[1] > y) && x < (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]) + a[0] {
                crossings += 1;
            }
        }
    }
    crossings % 2 == 1
}
#[test]
fn small_knuckle_fixed_meridian_keeps_its_actual_crest_circle_with_radial_offset() {
    let pitch = 2.54;
    let a = 15_f64.to_radians();
    let radius = pitch * (a.cos() - a.sin()) / (4. * (1. - a.sin()));
    let p = KnuckleHelicalThread {
        major_radius_mm: 4.995,
        pitch_mm: pitch,
        length_mm: 3.37 * pitch,
        crest_radius_mm: radius,
        root_radius_mm: radius,
        left_hand: false,
    };
    let shape = rust_knuckle(p, EntityId::new()).unwrap();
    let dz = -0.65 * radius * a.cos();
    let boundary = p.major_radius_mm - radius + (radius * radius - dz * dz).sqrt();
    let before = step_bytes(&shape);
    let mesh = shape
        .mesh(multicad_kernel::MeshOptions {
            linear_deflection_mm: 0.001,
            angular_deflection_rad: 0.04,
            maximum_triangles: 2_000_000,
        })
        .unwrap();
    assert_eq!(step_bytes(&shape), before);
    let area: f64 = shape.faces().unwrap().iter().map(|f| f.area_mm2).sum();
    let expected = analytic_volume(p);
    assert!(
        (shape.volume_mm3().unwrap() - expected).abs()
            < GeometryPolicy::default().linear_tolerance_mm * area + 1e-7 * expected
    );
    for phi in [0_f64, 0.43] {
        let z = 2. * pitch + dz + phi * pitch / std::f64::consts::TAU;
        for (dr, inside) in [(-0.001, true), (0.001, false)] {
            let r = boundary + dr;
            let point = [r * phi.cos(), r * phi.sin(), z];
            let fraction = material(&shape, point);
            assert!(
                if inside {
                    fraction > 0.999
                } else {
                    fraction < 1e-5
                },
                "native phi={phi}, dr={dr}: {fraction}"
            );
            assert_eq!(
                mesh_material_at(&mesh, point),
                inside,
                "mesh phi={phi}, dr={dr}"
            );
        }
    }
}
