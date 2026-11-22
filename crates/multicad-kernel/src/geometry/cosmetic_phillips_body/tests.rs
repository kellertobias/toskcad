use super::*;
use crate::geometry::helical_body::phillips::tests::fixture;
use multicad_rust_kernel::manufacturing::certify;

fn frames() -> [PlaneFrame; 2] {
    [
        PlaneFrame::default(),
        PlaneFrame {
            origin: [10., 20., 30.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        },
    ]
}

#[test]
fn cosmetic_phillips_keeps_physical_head_and_recess_on_one_closed_smooth_shaft() {
    let policy = GeometryPolicy::default();
    for (index, plain) in [(2, 0.), (1, 0.), (1, 1.5)] {
        for frame in frames() {
            let feature = EntityId::new();
            let source = fixture(index, plain, false);
            let physical =
                Solid::rust_phillips_screw_from_source(source.clone(), feature, frame, policy)
                    .unwrap();
            let cosmetic =
                Solid::rust_cosmetic_phillips_screw(source.clone(), feature, frame, policy)
                    .unwrap();
            let body = cosmetic.cosmetic_phillips_body().unwrap();
            let roles = body.geometry.faces().to_vec();
            assert_eq!(
                roles
                    .iter()
                    .filter(|r| matches!(r, Role::ThreadWall(_)))
                    .count(),
                1
            );

            // Smooth envelope material: the head and recess are unchanged, the
            // helical grooves are filled up to the major radius.
            let volume = cosmetic.volume_mm3().unwrap();
            let threaded = physical.volume_mm3().unwrap();
            let [a, b] = source.visible_thread_interval_mm();
            let r = source.thread().spec().major_radius_mm;
            assert!(volume > threaded, "{volume} <= {threaded}");
            assert!(
                volume - threaded < PI * r * r * (b - a),
                "{volume} {threaded}"
            );

            // One closed display shell owning every cosmetic face.
            let mesh = cosmetic.mesh(Default::default()).unwrap();
            let options = multicad_rust_kernel::display::Tessellation {
                deflection: 0.01,
                angle: 0.5,
                triangle_limit: 2_000_000,
            };
            let local = body.geometry.mesh(options, &|| false).unwrap();
            certify(&local).unwrap();
            let used: BTreeSet<u32> = local.face_indices.iter().copied().collect();
            assert_eq!(used, (0..roles.len() as u32).collect());
            let placed: BTreeSet<u32> = mesh.face_indices.iter().copied().collect();
            assert_eq!(placed, used);

            // Faces: physical head/recess records, analytic cosmetic shaft.
            let faces = cosmetic.faces().unwrap();
            let physical_faces = physical.faces().unwrap();
            assert_eq!(faces.len(), roles.len());
            for (i, (face, role)) in faces.iter().zip(&roles).enumerate() {
                assert_eq!(face.index, i as u32);
                if !is_thread(*role) {
                    let original = &physical_faces[body.physical_index(*role)];
                    assert_eq!(face.area_mm2, original.area_mm2);
                    assert_eq!(face.plane, original.plane);
                }
            }
            let wall = &faces[body.index(Role::ThreadWall(0))];
            assert!((wall.area_mm2 - TAU * r * (b - a)).abs() < 1e-12 * wall.area_mm2);
            let bounds = cosmetic.bounds().unwrap();
            let tip = position(frame, [0., 0., b]);
            for (i, value) in tip.into_iter().enumerate() {
                assert!(bounds.min[i] <= value && value <= bounds.max[i]);
            }

            // Topology: physical head/recess references retained; the helical
            // references are replaced by the explicit cosmetic envelope.
            let topology = cosmetic.topology_snapshot(feature, 1).unwrap();
            let original = physical.topology_snapshot(feature, 1).unwrap();
            let names = |t: &TopologySnapshot, kind| -> BTreeSet<String> {
                t.elements
                    .iter()
                    .filter(|e| e.key.kind == kind)
                    .flat_map(|e| e.origins.iter().map(|o| o.role.clone()))
                    .collect()
            };
            let face_names = names(&topology, TopologyKind::Face);
            assert!(face_names.contains("thread.cosmetic.envelope"));
            assert!(face_names.iter().all(|n| !n.starts_with("thread/")));
            let edge_names = names(&topology, TopologyKind::Edge);
            assert!(edge_names.contains("cosmetic.shaft/tip"));
            assert!(edge_names.contains("cosmetic.shaft/start"));
            for name in names(&original, TopologyKind::Face) {
                if name.starts_with("screw.head/") || name.starts_with("screw.drive/") {
                    assert!(face_names.contains(&name), "{name}");
                }
            }
            for name in names(&original, TopologyKind::Edge) {
                if name.starts_with("screw.drive/") {
                    assert!(edge_names.contains(&name), "{name}");
                }
            }
            let faces_in_topology: Vec<_> = topology
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Face)
                .collect();
            assert_eq!(faces_in_topology.len(), roles.len());
            for (i, face) in faces_in_topology.iter().enumerate() {
                assert_eq!(face.key.ordinal, i as u32);
                assert!(!face.origins.is_empty());
                assert!(!face.boundary.is_empty());
            }
            // Every face boundary names a published edge.
            for face in &faces_in_topology {
                for edge in &face.boundary {
                    assert!(edge_names.contains(&edge.role), "{}", edge.role);
                }
            }

            // Cosmetic geometry never earns manufacturing or STEP exchange.
            assert_eq!(
                cosmetic
                    .manufacturing_mesh(Default::default())
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Unsupported
            );
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("cosmetic.step");
            std::fs::write(&path, b"keep destination").unwrap();
            assert!(cosmetic.write_step(&path).is_err());
            assert_eq!(std::fs::read(&path).unwrap(), b"keep destination");

            // A rigid move keeps the same local envelope and volume.
            let moved = cosmetic
                .transformed(RigidTransform {
                    rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
                    translation_mm: [3., -4., 5.],
                })
                .unwrap();
            assert_eq!(moved.volume_mm3().unwrap(), volume);
            assert_eq!(moved.faces().unwrap().len(), roles.len());
        }
    }
}
