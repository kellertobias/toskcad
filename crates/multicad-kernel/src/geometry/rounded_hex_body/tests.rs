use super::*;
use multicad_rust_kernel::{display::DisplayMesh, manufacturing::certify};

/// DIN 405 Rd 20 x 1/8": external shaft and the internal cutter with its
/// distinct crest/root circles and radial clearance.
fn rd(internal: bool, left_hand: bool) -> Spec {
    let p = 3.175;
    let a = 15_f64.to_radians();
    let r = p * (a.cos() - a.sin()) / (4. * (1. - a.sin()));
    let d = if internal {
        0.05 * p * a.sin() / (1. - a.sin())
    } else {
        0.
    };
    Spec {
        major_radius_mm: 10. + if internal { 0.05 * p } else { 0. },
        pitch_mm: p,
        length_mm: 3.37 * p,
        depth_mm: p / 2.,
        crest_radius_mm: r - d,
        root_radius_mm: r + d,
        flank_angle_deg: 15.,
        starts: 1,
        left_hand,
    }
}
fn frame() -> PlaneFrame {
    PlaneFrame {
        origin: [4., -3., 7.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    }
}

#[test]
fn rounded_hex_nut_and_screw_have_exact_volume_closed_mesh_faces_and_topology() {
    let policy = GeometryPolicy::default();
    for internal in [true, false] {
        for left in [false, true] {
            let feature = EntityId::new();
            let spec = rd(internal, left);
            let solid = if internal {
                Solid::rust_rounded_hex_nut(spec, 30., 1., feature, frame(), policy)
            } else {
                Solid::rust_rounded_hex_head_screw(spec, 30., 3., feature, frame(), policy)
            }
            .unwrap();
            let (stock, _) = solid.rust_rounded_hex_data().unwrap();
            let thread = multicad_rust_kernel::rounded_thread::Thread::new(spec).unwrap();
            let prism = 3_f64.sqrt() / 2. * 900. * if internal { spec.length_mm } else { 3. };
            let expected = if internal {
                prism - thread.volume_mm3()
            } else {
                prism + thread.volume_mm3()
            };
            let volume = solid.volume_mm3().unwrap();
            assert!(
                (volume - expected).abs() <= 1e-9 * expected,
                "{volume} {expected}"
            );
            assert_eq!(stock.volume_mm3(), volume);

            // One closed certified shell whose signed volume meets the bound.
            let options = MeshOptions {
                linear_deflection_mm: 0.02,
                angular_deflection_rad: 0.3,
                maximum_triangles: 2_000_000,
            };
            let mesh = solid.manufacturing_mesh(options).unwrap();
            let proof = certify(&DisplayMesh {
                vertices: mesh.vertices.clone(),
                indices: mesh.indices.clone(),
                face_indices: mesh.face_indices.clone(),
            })
            .unwrap();
            assert!(
                (proof.signed_volume_mm3 - volume).abs() < 0.05 * volume.cbrt(),
                "{} {volume}",
                proof.signed_volume_mm3
            );
            let faces = solid.faces().unwrap();
            assert_eq!(faces.len(), stock.face_count());
            let used: BTreeSet<u32> = mesh.face_indices.iter().copied().collect();
            assert_eq!(used, faces.iter().map(|f| f.index).collect());
            // Every planar face's mesh triangles lie in its stated plane, and
            // the stated areas agree with the mesh to its deflection.
            let v = mesh.vertices.as_chunks::<3>().0;
            let mut areas = vec![0.; faces.len()];
            for (t, &f) in mesh
                .indices
                .as_chunks::<3>()
                .0
                .iter()
                .zip(&mesh.face_indices)
            {
                let [a, b, c] = t.map(|i| v[i as usize]);
                let n = cross(
                    std::array::from_fn(|i| b[i] - a[i]),
                    std::array::from_fn(|i| c[i] - a[i]),
                );
                areas[f as usize] += 0.5 * dot(n, n).sqrt();
                if let Some((origin, normal)) = faces[f as usize].plane {
                    for p in [a, b, c] {
                        let d = dot(std::array::from_fn(|i| p[i] - origin[i]), normal);
                        assert!(d.abs() < 1e-9, "face {f} off plane {d}");
                    }
                    assert!(dot(n, normal) > 0., "face {f} normal");
                }
            }
            for face in &faces {
                let relative = (areas[face.index as usize] - face.area_mm2).abs() / face.area_mm2;
                assert!(
                    relative < 0.02,
                    "face {} {} {}",
                    face.index,
                    areas[face.index as usize],
                    face.area_mm2
                );
            }

            // Stable semantic topology, including the hex walls and corners.
            let topology = solid.topology_snapshot(feature, 1).unwrap();
            let roles: BTreeSet<String> = topology
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Face)
                .flat_map(|e| e.origins.iter().map(|o| o.role.clone()))
                .collect();
            let prefix = if internal { "nut" } else { "head" };
            for k in 0..6 {
                assert!(roles.contains(&format!("{prefix}.wall/{k}")), "{roles:?}");
            }
            assert_eq!(
                roles
                    .iter()
                    .any(|r| r.starts_with("receiving.inward/thread/")),
                internal
            );
            assert_eq!(roles.contains("head.seat"), !internal);
            assert_eq!(
                topology
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == TopologyKind::Face)
                    .count(),
                faces.len()
            );

            // Bounds hold the hex corners (and the screw tip); membership.
            let bounds = solid.bounds().unwrap();
            for p in v {
                for (i, x) in p.iter().enumerate() {
                    assert!(bounds.min[i] <= x + 1e-9 && *x <= bounds.max[i] + 1e-9);
                }
            }
            let local = |p: Point| position(frame(), p);
            assert_eq!(
                solid.rust_thread_contains(local([0., 0., 1.])).unwrap(),
                !internal
            );
            assert!(
                solid
                    .rust_thread_contains(local([14.9, 0., if internal { 1. } else { -1. }]))
                    .unwrap()
            );
            assert!(!solid.rust_thread_contains(local([15.1, 0., -1.])).unwrap());
            solid.display_edges(0.05).unwrap();
            assert!(!solid.supports_step_export());

            let moved = solid
                .transformed(RigidTransform {
                    rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
                    translation_mm: [1., 2., 3.],
                })
                .unwrap();
            assert_eq!(moved.volume_mm3().unwrap(), volume);
            moved.topology_snapshot(feature, 2).unwrap();
        }
    }
}

#[test]
fn rounded_hex_hardware_refuses_thin_walls_and_seats() {
    let policy = GeometryPolicy::default();
    let spec = rd(true, false);
    for (width, wall) in [(21., 1.), (30., 6.)] {
        assert_eq!(
            Solid::rust_rounded_hex_nut(spec, width, wall, EntityId::new(), frame(), policy)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Input
        );
    }
    assert!(
        Solid::rust_rounded_hex_head_screw(
            rd(false, false),
            20.,
            3.,
            EntityId::new(),
            frame(),
            policy
        )
        .is_err()
    );
    assert!(
        Solid::rust_rounded_hex_head_screw(
            rd(false, false),
            30.,
            0.,
            EntityId::new(),
            frame(),
            policy
        )
        .is_err()
    );
}
