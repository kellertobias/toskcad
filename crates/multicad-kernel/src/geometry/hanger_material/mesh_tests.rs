use super::source_tests::{source, volume_options};
use super::*;

#[test]
fn actual_uncapped_metric_wall_is_not_a_closed_manufacturing_solid() {
    let source = source(false, false);
    let before = source;
    let shell = source
        .metric()
        .tessellate_shell(
            multicad_rust_kernel::display::Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            16_000_000,
            None,
        )
        .unwrap();
    assert!(shell.source_chord_bound_mm() <= 0.1);
    assert!(!shell.mesh().indices.is_empty());
    let failure = mesh::certify(shell.mesh()).unwrap_err();
    assert_eq!(failure.kind, GeometryErrorKind::InvalidResult);
    assert_eq!(failure.operation, "Rust authored hanger manufacturing mesh");
    assert_eq!(source, before);
}

#[test]
fn connected_hanger_solid_retains_source_through_share_and_actual_pose() {
    let source = source(true, true);
    let feature = EntityId::new();
    let policy = GeometryPolicy::default();
    let solid = Solid::rust_hanger(source, feature, PlaneFrame::default(), policy, None).unwrap();
    assert_eq!(solid.solid_count().unwrap(), 1);
    assert_eq!(
        solid.rust_hanger_data(),
        Some(HangerData {
            source,
            frame: PlaneFrame::default()
        })
    );
    let shared = solid.share().unwrap();
    assert_eq!(shared.rust_hanger_data(), solid.rust_hanger_data());
    let original_volume = solid.volume_mm3().unwrap();
    let transform = RigidTransform {
        rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
        translation_mm: [11., -7., 5.],
    };
    let placed = solid.transformed(transform).unwrap();
    assert_eq!(
        placed.rust_hanger_data(),
        Some(HangerData {
            source,
            frame: transform.frame(PlaneFrame::default())
        })
    );
    assert!((placed.volume_mm3().unwrap() - original_volume).abs() < 1e-7);
    let proof = placed
        .manufacturing_mesh_with_receipt(MeshOptions::default())
        .unwrap();
    assert!(proof.edge_pairs > 0 && proof.vertex_links > 0);
    let bounds = placed.bounds().unwrap();
    for point in proof.mesh.vertices.as_chunks::<3>().0 {
        for (axis, coordinate) in point.iter().enumerate() {
            assert!(bounds.min[axis] <= *coordinate && *coordinate <= bounds.max[axis]);
        }
    }
    assert_eq!(
        solid.rust_hanger_data().unwrap().frame,
        PlaneFrame::default()
    );
    assert_eq!(solid.rust_hanger_data().unwrap().source, source);
}

#[test]
fn incomplete_hanger_interfaces_refuse_without_occt_or_partial_files() {
    let source = source(false, false);
    let feature = EntityId::new();
    let policy = GeometryPolicy::default();
    assert_eq!(
        Solid::rust_hanger(source, EntityId::nil(), PlaneFrame::default(), policy, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        Solid::rust_hanger(
            source,
            feature,
            PlaneFrame::default(),
            policy,
            Some(&|| true)
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    let solid = Solid::rust_hanger(source, feature, PlaneFrame::default(), policy, None).unwrap();
    let failure = solid.faces().unwrap_err();
    assert_eq!(failure.kind, GeometryErrorKind::Unsupported);
    assert_eq!(failure.operation, "Rust hanger references");
    assert_eq!(
        solid.topology_snapshot(feature, 0).unwrap_err().operation,
        "Rust hanger references"
    );
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("existing.step");
    std::fs::write(&path, b"retained previous file").unwrap();
    assert!(solid.write_step(&path).is_err());
    assert_eq!(std::fs::read(&path).unwrap(), b"retained previous file");
}
#[test]
fn authored_hanger_world_manufacturing_retains_source_quality_and_closed_material() {
    for (left, nonbinary, frame) in [
        (false, false, PlaneFrame::default()),
        (
            true,
            true,
            PlaneFrame {
                origin: [11., -7., 5.],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
        ),
    ] {
        let source = source(left, nonbinary);
        let before = source;
        let material = HangerMaterial::new(source, frame, GeometryPolicy::default(), None).unwrap();
        let proof = material
            .manufacturing_mesh(MeshOptions::default(), None)
            .unwrap();
        assert!(proof.edge_pairs > 0 && proof.vertex_links > 0);
        assert_eq!(proof.removed_exact_zero_triangles, 0);
        assert_eq!(proof.maximum_displacement_mm, 0.);
        let bounds = material.bounds();
        for point in proof.mesh.vertices.as_chunks::<3>().0 {
            for (axis, coordinate) in point.iter().enumerate() {
                assert!(bounds.min[axis] <= *coordinate && *coordinate <= bounds.max[axis]);
            }
        }
        let volume = material.volume_bounds(volume_options(), None).unwrap();
        assert!(volume[0] > 0. && volume[1] >= volume[0] && volume[1] - volume[0] <= 1e-7);
        assert_eq!(*material.source(), before);
        assert_eq!(material.frame(), frame);
        assert_eq!(source, before);
    }
}

#[test]
fn hanger_exact_quarter_placement_retains_actual_source_volume_and_failure_atomicity() {
    let source = source(false, true);
    let before = source;
    let frame = PlaneFrame {
        origin: [10., -7., 3.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    };
    let material = HangerMaterial::new(source, frame, GeometryPolicy::default(), None).unwrap();
    let local = source.volume(volume_options(), None).unwrap().material_mm3;
    let world = material.volume_bounds(volume_options(), None).unwrap();
    assert!(world[0] <= local[0] && world[1] >= local[1]);
    assert_eq!(
        material
            .mesh(MeshOptions::default(), Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert_eq!(
        HangerMaterial::new(source, frame, GeometryPolicy::default(), Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(
        HangerMaterial::new(
            source,
            PlaneFrame {
                u: [1.01, 0., 0.],
                ..frame
            },
            GeometryPolicy::default(),
            None
        )
        .is_err()
    );
    assert!(
        material
            .mesh(
                MeshOptions {
                    maximum_triangles: 32,
                    ..Default::default()
                },
                None
            )
            .is_err()
    );
    assert_eq!(
        material
            .volume_bounds(
                QueryOptions {
                    volume_error_mm3: 0.,
                    ..volume_options()
                },
                None
            )
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(*material.source(), before);
}

#[test]
fn cancelled_hanger_final_manufacturing_checkpoint_cannot_publish_a_partial_mesh() {
    let source = source(true, false);
    let material = HangerMaterial::new(
        source,
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    material
        .manufacturing_mesh(MeshOptions::default(), Some(&count))
        .unwrap();
    let last = calls.get();
    assert!(last > 10);
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == last
    };
    assert_eq!(
        material
            .manufacturing_mesh(MeshOptions::default(), Some(&cancel))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert_eq!(calls.get(), last);
    assert_eq!(*material.source(), source);
}
