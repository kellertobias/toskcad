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
