use super::*;
const BYTES: usize = 64 * 1024 * 1024;
const WORK: usize = 512_000_000;
fn quality() -> Tessellation {
    Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    }
}
fn candidate(quarter: bool) -> StockClipMeshCandidate {
    let source = super::super::super::super::tests::canonical(quarter, 1.);
    let atlas = source.stock_clip_plan(Default::default(), None).unwrap();
    atlas
        .cell_plan(128, 128, WORK, BYTES, None)
        .unwrap()
        .mesh_candidate(quality(), WORK, BYTES, None)
        .unwrap()
}
fn identity() -> [[f64; 3]; 3] {
    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
}
#[test]
fn actual_source_placed_mesh_retains_connectivity_source_quality_and_independent_affine_geometry() {
    for quarter in [false, true] {
        for matrix in [
            identity(),
            [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
            [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
        ] {
            let c = candidate(quarter);
            let old = c.mesh().clone();
            let original_volume = c.source_volume_mm3();
            let pose = StockClipPlacement::new(matrix, [17., -8., 13.]).unwrap();
            assert_eq!(pose.matrix(), matrix);
            assert_eq!(pose.translation_mm(), [17., -8., 13.]);
            let determinant = pose.determinant_bounds();
            assert!(determinant[0] > 0.);
            assert!(pose.metric_error() >= 0. && pose.metric_error() <= 1e-6);
            assert!(pose.area_scale_bounds()[0] > 0.);
            let p = c.placed(pose, quality(), WORK, BYTES, None).unwrap();
            assert_eq!(p.mesh().indices, old.indices);
            assert_eq!(p.mesh().face_indices, old.face_indices);
            assert_eq!(p.mesh().vertices.len(), old.vertices.len());
            assert!(p.source_distance_mm() <= quality().deflection);
            assert!(p.normal_chord() <= 2. * (quality().angle / 2.).sin());
            for (local, world) in old
                .vertices
                .as_chunks::<3>()
                .0
                .iter()
                .zip(p.mesh().vertices.as_chunks::<3>().0)
            {
                for i in 0..3 {
                    let oracle =
                        (0..3).map(|j| matrix[i][j] * local[j]).sum::<f64>() + [17., -8., 13.][i];
                    assert!((world[i] - oracle).abs() < 1e-12);
                }
            }
            let v = p.source_volume_mm3();
            assert!(v[0] <= original_volume[0] && v[1] >= original_volume[1]);
            assert!(p.literal_volume_mm3()[0] > 0.);
            assert!(p.work_used() <= WORK && p.retained_budget_bytes() <= BYTES);
            assert_eq!(
                manufacturing::certify(p.mesh()).unwrap().edge_pairs,
                p.connectivity().edge_pairs
            );
            let mesh = p.into_mesh();
            assert_eq!(mesh.indices, old.indices);
        }
    }
}
#[test]
fn actual_placement_refusals_honor_quality_total_work_memory_and_cancellation() {
    for value in [f64::NAN, f64::INFINITY] {
        let mut matrix = identity();
        matrix[1][0] = value;
        assert_eq!(
            StockClipPlacement::new(matrix, [0.; 3]).unwrap_err(),
            Error::Invalid
        );
        assert_eq!(
            StockClipPlacement::new(identity(), [value, 0., 0.]).unwrap_err(),
            Error::Invalid
        );
    }
    assert_eq!(
        StockClipPlacement::new(identity(), [1e7_f64.next_up(), 0., 0.]).unwrap_err(),
        Error::Limit
    );
    let mut matrix = identity();
    matrix[0][0] = 3.;
    assert_eq!(
        StockClipPlacement::new(matrix, [0.; 3]).unwrap_err(),
        Error::Limit
    );
    matrix[0][0] = 1.01;
    assert_eq!(
        StockClipPlacement::new(matrix, [0.; 3]).unwrap_err(),
        Error::Unsupported
    );
    matrix[0][0] = -1.;
    assert_eq!(
        StockClipPlacement::new(matrix, [0.; 3]).unwrap_err(),
        Error::Unsupported
    );
    let pose = StockClipPlacement::new(identity(), [17., -8., 13.]).unwrap();
    for (work, bytes) in [(0, BYTES), (WORK + 1, BYTES), (WORK, 0), (WORK, BYTES + 1)] {
        assert_eq!(
            candidate(false)
                .placed(pose, quality(), work, bytes, None)
                .unwrap_err(),
            Error::Invalid
        );
    }
    let c = candidate(false);
    let used = c.work_used();
    assert_eq!(
        c.placed(pose, quality(), used + 1, BYTES, None)
            .unwrap_err(),
        Error::Limit
    );
    let c = candidate(false);
    let bytes = c.retained_budget_bytes();
    assert_eq!(
        c.placed(pose, quality(), WORK, bytes + 1, None)
            .unwrap_err(),
        Error::Limit
    );
    assert_eq!(
        candidate(false)
            .placed(
                pose,
                Tessellation {
                    deflection: 1e-12,
                    ..quality()
                },
                WORK,
                BYTES,
                None
            )
            .unwrap_err(),
        Error::Unresolved
    );
    assert_eq!(
        candidate(false)
            .placed(
                pose,
                Tessellation {
                    angle: 0.01,
                    ..quality()
                },
                WORK,
                BYTES,
                None
            )
            .unwrap_err(),
        Error::Unresolved
    );
    let far = StockClipPlacement::new(identity(), [1e7, 0., 0.]).unwrap();
    assert_eq!(
        candidate(false)
            .placed(far, quality(), WORK, BYTES, None)
            .unwrap_err(),
        Error::Limit
    );
    let calls = std::cell::Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    candidate(false)
        .placed(pose, quality(), WORK, BYTES, Some(&count))
        .unwrap();
    for stop in [0, 1, calls.get() - 1] {
        let c = candidate(false);
        let calls = std::cell::Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() > stop
        };
        assert_eq!(
            c.placed(pose, quality(), WORK, BYTES, Some(&cancel))
                .unwrap_err(),
            Error::Cancelled
        );
    }
}

#[test]
fn interval_normal_and_literal_matrix_roundoff_are_independent_numeric_obligations() {
    let matrix = [[1. + 4e-12, 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let pose = StockClipPlacement::new(matrix, [0.; 3]).unwrap();
    let determinant = pose.determinant_bounds();
    assert!(
        determinant[0] > 1. && determinant[0] <= matrix[0][0] && matrix[0][0] <= determinant[1]
    );
    assert_eq!(
        pose.normal([I::point(0.); 3]).unwrap_err(),
        Error::Unresolved
    );
    assert_eq!(
        pose.normal([I::new(-1., 1.); 3]).unwrap_err(),
        Error::Unresolved
    );
    let c = candidate(false);
    let original = c.source_volume_mm3();
    let placed = c.placed(pose, quality(), WORK, BYTES, None).unwrap();
    let volume = placed.source_volume_mm3();
    assert!(volume[0] <= original[0] * matrix[0][0] && volume[1] >= original[1] * matrix[0][0]);
    assert_eq!(
        candidate(false)
            .placed(
                pose,
                Tessellation {
                    angle: 0.,
                    ..quality()
                },
                WORK,
                BYTES,
                None
            )
            .unwrap_err(),
        Error::Invalid
    );
}

#[test]
fn actual_affine_inverse_preserves_boxes_without_transpose_approximation() {
    for matrix in [
        identity(),
        [[1. + 4e-12, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
    ] {
        let pose = StockClipPlacement::new(matrix, [17., -8., 13.]).unwrap();
        let local = [[1., 1.], [2., 2.], [3., 3.]];
        let world = pose.point_bounds(local).unwrap();
        let back = pose.inverse_point_bounds(world).unwrap();
        for i in 0..3 {
            assert!(back[i][0] <= local[i][0] && back[i][1] >= local[i][1]);
            assert!(back[i][1] - back[i][0] < 1e-11);
        }
    }
    let pose = StockClipPlacement::new(identity(), [17., -8., 13.]).unwrap();
    for input in [[[f64::NAN, 0.]; 3], [[0., f64::INFINITY]; 3], [[1., 0.]; 3]] {
        assert_eq!(pose.point_bounds(input).unwrap_err(), Error::Invalid);
        assert_eq!(
            pose.inverse_point_bounds(input).unwrap_err(),
            Error::Invalid
        );
    }
    assert_eq!(
        pose.point_bounds([[1e7_f64.next_up(); 2]; 3]).unwrap_err(),
        Error::Limit
    );
    assert_eq!(
        pose.inverse_point_bounds([[1e7_f64.next_up(); 2]; 3])
            .unwrap_err(),
        Error::Limit
    );
    assert_eq!(pose.point_bounds([[1e7; 2]; 3]).unwrap_err(), Error::Limit);
    assert_eq!(
        pose.inverse_point_bounds([[-1e7; 2]; 3]).unwrap_err(),
        Error::Limit
    );
}
