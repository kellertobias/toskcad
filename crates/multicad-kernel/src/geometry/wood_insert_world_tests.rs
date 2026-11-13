//! Private integration probe: these candidates are not admitted as Solid geometry.
use super::*;
use multicad_rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    helical_thread::Spec,
    wood_insert::Insert,
};
fn source(left: bool) -> Insert {
    let quality = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let outer = Spec {
        major_radius_mm: 5.,
        pitch_mm: 2.,
        length_mm: 4.,
        depth_mm: 1.,
        crest_width_mm: 0.4,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    };
    let pitch_mid = (2.675 + 2.775) / 2.;
    let minor_mid = (2.459 + 2.599) / 2.;
    let inner = Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 4.,
        depth_mm: (3. - minor_mid) / 2.,
        crest_width_mm: 0.25 - (3. - pitch_mid) / 3_f64.sqrt(),
        ..outer
    };
    Insert::new(
        AuthoredThread::new(outer, quality).unwrap(),
        AuthoredThread::new(inner, quality).unwrap(),
        0.5,
        None,
    )
    .unwrap()
}
#[test]
fn actual_paired_insert_world_mesh_retains_quality_and_literal_embedding() {
    for (left, frame) in [
        (false, PlaneFrame::default()),
        (
            true,
            PlaneFrame {
                origin: [11., -7., 5.],
                u: [0.6, 0.8, 0.],
                v: [-0.8, 0.6, 0.],
            },
        ),
        (
            false,
            PlaneFrame {
                origin: [1e4, -2e4, 3e4],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
        ),
    ] {
        let source = source(left);
        let before = source;
        let timer = std::time::Instant::now();
        let material =
            WoodInsertMaterial::new(source, frame, GeometryPolicy::default(), None).unwrap();
        let proof = material
            .manufacturing_mesh(MeshOptions::default(), None)
            .unwrap();
        assert!(proof.edge_pairs > 0 && proof.vertex_links > 0);
        assert_eq!(proof.removed_exact_zero_triangles, 0);
        assert_eq!(material.source(), &before);
        assert_eq!(material.frame(), frame);
        let bounds = material.bounds().unwrap();
        let volume = material
            .volume_bounds(
                multicad_rust_kernel::wood_insert::QueryOptions {
                    volume_error_mm3: 1e-7,
                    maximum_cells: 10,
                },
                None,
            )
            .unwrap();
        assert!(volume[0] > 0. && volume[1] >= volume[0]);
        let mesh = proof.mesh;
        for p in mesh.vertices.as_chunks::<3>().0 {
            for (i, coordinate) in p.iter().enumerate() {
                assert!(bounds.min[i] <= *coordinate && *coordinate <= bounds.max[i]);
            }
        }
        assert!(!mesh.indices.is_empty());
        assert_eq!(source, before);
        eprintln!(
            "actual insert world LH{left} {} facets {:?}",
            mesh.indices.len() / 3,
            timer.elapsed()
        );
    }
}

#[test]
fn paired_insert_native_failures_and_cancellation_leave_source_unchanged() {
    let source = source(false);
    let before = source;
    let policy = GeometryPolicy::default();
    assert_eq!(
        WoodInsertMaterial::new(source, PlaneFrame::default(), policy, Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(
        WoodInsertMaterial::new(
            source,
            PlaneFrame {
                u: [1.01, 0., 0.],
                ..PlaneFrame::default()
            },
            policy,
            None
        )
        .is_err()
    );
    let material = WoodInsertMaterial::new(source, PlaneFrame::default(), policy, None).unwrap();
    use multicad_rust_kernel::wood_insert::QueryOptions;
    assert_eq!(
        material
            .volume_bounds(
                QueryOptions {
                    volume_error_mm3: 0.,
                    maximum_cells: 10
                },
                None
            )
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        material
            .volume_bounds(
                QueryOptions {
                    volume_error_mm3: 1e-7,
                    maximum_cells: 9
                },
                None
            )
            .unwrap_err()
            .kind,
        GeometryErrorKind::Limit
    );
    assert_eq!(
        material
            .mesh(MeshOptions::default(), Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(
        material
            .mesh(
                MeshOptions {
                    maximum_triangles: 32,
                    ..MeshOptions::default()
                },
                None
            )
            .is_err()
    );
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
    assert_eq!(*material.source(), before);
    assert_eq!(source, before);
}

#[test]
fn true_source_numeric_floor_and_open_shell_refuse_native_publication() {
    let source = source(false);
    let options = multicad_rust_kernel::wood_insert::QueryOptions {
        volume_error_mm3: 1e-7,
        maximum_cells: 10,
    };
    let material = WoodInsertMaterial::new(
        source,
        PlaneFrame {
            origin: [11., -7., 5.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        },
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let local = source.volume(options, None).unwrap().material_mm3;
    let tight = multicad_rust_kernel::wood_insert::QueryOptions {
        volume_error_mm3: (local[1] - local[0]).next_up(),
        ..options
    };
    assert!(source.volume(tight, None).is_ok());
    assert_eq!(
        material.volume_bounds(tight, None).unwrap_err().kind,
        GeometryErrorKind::Unsupported
    );
    let shell = source
        .outer()
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
    assert_eq!(
        wood_insert_material::certify(shell.mesh())
            .unwrap_err()
            .kind,
        GeometryErrorKind::InvalidResult
    );
    let huge_policy = GeometryPolicy {
        maximum_coordinate_mm: f64::MAX,
        ..GeometryPolicy::default()
    };
    assert!(
        WoodInsertMaterial::new(
            source,
            PlaneFrame {
                origin: [f64::MAX, 0., 0.],
                ..PlaneFrame::default()
            },
            huge_policy,
            None
        )
        .is_err()
    );
}

#[test]
fn large_sources_construct_without_an_unrequested_absolute_volume_precision() {
    let base = source(false);
    assert!(
        wood_insert_material::placed_bounds(
            &base,
            PlaneFrame {
                origin: [f64::INFINITY, 0., 0.],
                ..PlaneFrame::default()
            }
        )
        .is_err()
    );
    let scaled = |spec: Spec| Spec {
        major_radius_mm: spec.major_radius_mm * 100.,
        pitch_mm: spec.pitch_mm * 100.,
        length_mm: spec.length_mm * 100.,
        depth_mm: spec.depth_mm * 100.,
        crest_width_mm: spec.crest_width_mm * 100.,
        ..spec
    };
    let options = SourceOptions {
        position_limit_mm: 1e-2,
        normal_angle_limit_rad: 1e-2,
    };
    let source = Insert::new(
        AuthoredThread::new(scaled(base.outer().spec()), options).unwrap(),
        AuthoredThread::new(scaled(base.inner().spec()), options).unwrap(),
        50.,
        None,
    )
    .unwrap();
    let material = WoodInsertMaterial::new(
        source,
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let options = multicad_rust_kernel::wood_insert::QueryOptions {
        volume_error_mm3: 1e-7,
        maximum_cells: 10,
    };
    assert_eq!(
        material.volume_bounds(options, None).unwrap_err().kind,
        GeometryErrorKind::Unsupported
    );
    assert!(
        material
            .volume_bounds(
                multicad_rust_kernel::wood_insert::QueryOptions {
                    volume_error_mm3: 0.01,
                    ..options
                },
                None
            )
            .unwrap()[0]
            > 0.
    );
}
