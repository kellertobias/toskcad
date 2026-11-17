use multicad_core::{document::EntityId, references::TopologyKind};
use multicad_kernel::{
    GeometryPolicy, MeshOptions, PlaneFrame, RoundedBlindMeshBudget, Solid,
    rust_kernel::{Rational, rounded_thread::Spec},
};
fn spec(left_hand: bool) -> Spec {
    let p = 1.27;
    let a = 27.5_f64.to_radians();
    let r = p * a.cos() / (12. * (1. - a.sin()));
    Spec {
        major_radius_mm: 3.175,
        pitch_mm: p,
        length_mm: 3.81,
        depth_mm: p / (3. * a.tan()),
        crest_radius_mm: r,
        root_radius_mm: r,
        flank_angle_deg: 27.5,
        starts: 1,
        left_hand,
    }
}
fn make(left: bool, frame: PlaneFrame) -> Solid {
    Solid::rust_rounded_blind_bore(
        spec(left),
        Rational::new(127, 20).unwrap(),
        Rational::new(127, 25).unwrap(),
        EntityId::new(),
        frame,
        GeometryPolicy::default(),
    )
    .unwrap()
}
fn quality() -> MeshOptions {
    MeshOptions {
        linear_deflection_mm: 0.05,
        angular_deflection_rad: 0.3,
        maximum_triangles: 500_000,
    }
}
#[test]
fn authored_blind_retains_full_cutter_true_floor_and_stock_top_with_fresh_world_proof() {
    for left in [false, true] {
        for posed in [false, true] {
            let frame = if posed {
                PlaneFrame {
                    origin: [17., -8., 13.],
                    u: [0.8, 0.6, 0.],
                    v: [-0.36, 0.48, 0.8],
                }
            } else {
                PlaneFrame::default()
            };
            let solid = make(left, frame);
            let data = solid.rust_rounded_blind_bore_data().unwrap();
            assert_eq!(data.bore.cutter().spec(), spec(left));
            let expected =
                std::f64::consts::PI * 6.35_f64.powi(2) * 5.08 - data.bore.cutter().volume_mm3();
            assert!((solid.volume_mm3().unwrap() - expected).abs() < 1e-9);
            let faces = solid.faces().unwrap();
            assert_eq!(faces.len(), 9);
            assert!(faces.iter().all(|f| f.area_mm2 > 0.));
            let floor = faces.iter().find(|f| f.index == 1).unwrap();
            let top = faces.iter().find(|f| f.index == 8).unwrap();
            for k in 0..3 {
                assert!(
                    (floor.plane.unwrap().0[k] - (frame.origin[k] + frame.normal()[k] * 3.81))
                        .abs()
                        < 1e-12
                );
                assert!(
                    (top.plane.unwrap().0[k] - (frame.origin[k] + frame.normal()[k] * 5.08)).abs()
                        < 1e-12
                );
                assert!((floor.plane.unwrap().1[k] + frame.normal()[k]).abs() < 1e-12);
            }
            let topology = solid.topology_snapshot(EntityId::new(), 1).unwrap();
            assert_eq!(
                topology
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == TopologyKind::Face)
                    .count(),
                9
            );
            let (mesh, receipt) = solid
                .rust_rounded_blind_bore_mesh_with_receipt(
                    quality(),
                    RoundedBlindMeshBudget::default(),
                    None,
                )
                .unwrap();
            assert_eq!(receipt.source(), &data.bore);
            assert!(receipt.embedding_pairs() > 0);
            assert!(mesh.face_indices.iter().all(|&f| f < 9));
            assert!((0..9).all(|i| mesh.face_indices.contains(&i)));
            for (&f, ids) in mesh
                .face_indices
                .iter()
                .zip(mesh.indices.as_chunks::<3>().0)
            {
                if f == 1 || f == 8 {
                    for &id in ids {
                        let p = &mesh.vertices[3 * id as usize..3 * id as usize + 3];
                        let z = (0..3)
                            .map(|k| (p[k] - frame.origin[k]) * frame.normal()[k])
                            .sum::<f64>();
                        assert!((z - if f == 1 { 3.81 } else { 5.08 }).abs() < 1e-10);
                    }
                }
            }
            let manufactured = solid.manufacturing_mesh_with_receipt(quality()).unwrap();
            assert_eq!(manufactured.mesh.indices, mesh.indices);
            assert_eq!(manufactured.mesh.vertices, mesh.vertices);
            assert!(manufactured.edge_pairs > 0);
        }
    }
}
#[test]
fn actual_query_resource_frontiers_cancel_and_unchanged_retry() {
    use multicad_kernel::GeometryErrorKind;
    use std::cell::Cell;
    let solid = make(false, PlaneFrame::default());
    let original = solid.rust_rounded_blind_bore_data().unwrap();
    let calls = Cell::new(0usize);
    let poll = || {
        calls.set(calls.get() + 1);
        false
    };
    let (mesh, proof) = solid
        .rust_rounded_blind_bore_mesh_with_receipt(quality(), Default::default(), Some(&poll))
        .unwrap();
    let used = RoundedBlindMeshBudget::default().maximum_work - proof.remaining_work();
    let bytes = proof.native_buffer_peak_bytes();
    let polls = calls.get();
    assert!(polls > 10 && used > proof.producer_reserved_work());
    for budget in [
        RoundedBlindMeshBudget {
            maximum_work: used - 1,
            ..Default::default()
        },
        RoundedBlindMeshBudget {
            maximum_native_buffer_bytes: bytes - 1,
            ..Default::default()
        },
    ] {
        let error = solid
            .rust_rounded_blind_bore_mesh_with_receipt(quality(), budget, None)
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Limit);
    }
    let exact = RoundedBlindMeshBudget {
        maximum_work: used,
        maximum_native_buffer_bytes: bytes,
    };
    let (retry, receipt) = solid
        .rust_rounded_blind_bore_mesh_with_receipt(quality(), exact, None)
        .unwrap();
    assert_eq!(receipt.remaining_work(), 0);
    assert_eq!(retry.vertices, mesh.vertices);
    assert_eq!(retry.indices, mesh.indices);
    for stop in [1, polls / 2, polls] {
        let observed = Cell::new(0usize);
        let cancel = || {
            observed.set(observed.get() + 1);
            observed.get() >= stop
        };
        let error = solid
            .rust_rounded_blind_bore_mesh_with_receipt(quality(), Default::default(), Some(&cancel))
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Cancelled);
        assert!(observed.get() >= stop);
        assert_eq!(solid.rust_rounded_blind_bore_data().unwrap(), original);
    }
    assert_eq!(
        solid.rust_rounded_blind_bore_contains([0., 0., 4.5]),
        Some(true)
    );
    assert_eq!(
        solid.rust_rounded_blind_bore_contains([0., 0., 2.]),
        Some(false)
    );
    assert_eq!(
        solid.rust_rounded_blind_bore_contains([7., 0., 4.5]),
        Some(false)
    );
    let (retry, _) = solid
        .rust_rounded_blind_bore_mesh_with_receipt(quality(), exact, None)
        .unwrap();
    assert_eq!(retry.vertices, mesh.vertices);
}
#[test]
fn source_known_step_and_owned_lifecycle_preserve_the_authored_floor() {
    let source = make(true, PlaneFrame::default());
    let shared = source.share().unwrap();
    let data = source.rust_rounded_blind_bore_data().unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("blind.step");
    source.write_step(&path).unwrap();
    let input = std::fs::read_to_string(path).unwrap();
    assert!(
        (source
            .verify_source_known_rounded_blind_step(&input)
            .unwrap()
            - source.volume_mm3().unwrap())
        .abs()
            < 1e-6
    );
    drop(source);
    assert_eq!(shared.rust_rounded_blind_bore_data().unwrap(), data);
    assert!(
        (shared
            .verify_source_known_rounded_blind_step(&input)
            .unwrap()
            - shared.volume_mm3().unwrap())
        .abs()
            < 1e-6
    );
}
