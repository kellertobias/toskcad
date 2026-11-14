//! Public accepted-source export regressions and occurrence ownership failures.
use super::*;
use multicad_rust_kernel::authored_flat_thread::{AuthoredThread, SourceOptions};
fn actual_insert(length: f64, left: bool, frame: PlaneFrame) -> Solid {
    let base = crate::geometry::wood_insert_world_tests::source(left);
    let quality = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let law = |original: &AuthoredThread| {
        let mut spec = original.spec();
        spec.length_mm = length;
        AuthoredThread::new(spec, quality).unwrap()
    };
    let paired = source::Insert::new(
        law(base.outer()),
        law(base.inner()),
        base.minimum_wall_mm(),
        None,
    )
    .unwrap();
    Solid::rust_wood_insert(
        paired,
        EntityId::new(),
        frame,
        GeometryPolicy::default(),
        None,
    )
    .unwrap()
}
#[test]
fn actual_eight_mm_inserts_export_with_unchanged_selected_member_limits() {
    for (left, frame) in [
        (false, PlaneFrame::default()),
        (
            true,
            PlaneFrame {
                origin: [31., -17., 9.],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
        ),
    ] {
        let solid = actual_insert(8., left, frame);
        let before = solid.rust_wood_insert_data().unwrap();
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("actual-eight-mm.step");
        solid
            .write_step(&path)
            .expect("standard8mm insert must export through normal Native path");
        let text = std::fs::read_to_string(path).unwrap();
        let proof = solid
            .verify_wood_insert_step_source(&text, 1e-7, None)
            .unwrap()
            .unwrap();
        assert!(proof.decoded_volume_bounds_mm3[0] > 0.);
        assert!(proof.decoded_volume_bounds_mm3[0] <= 462.74635728646706);
        assert!(proof.decoded_volume_bounds_mm3[1] >= 462.74635728646706);
        assert!(proof.maximum_position_error_mm <= 1e-7);
        assert!(proof.maximum_row_derivative_error_mm <= 1e-6);
        let compound_proof = solid
            .verify_partitioned_step_source(&text, 1e-7, None)
            .unwrap()
            .unwrap();
        assert_eq!(compound_proof.members.solid_count, 1);
        assert!(compound_proof.members.decoded_scalar_members.is_empty());
        assert_eq!(solid.rust_wood_insert_data().unwrap(), before);
    }
}
#[test]
fn single_source_verifier_rejects_multiple_real_occurrences_then_retries_unchanged() {
    let insert = actual_insert(4., false, PlaneFrame::default());
    let before = insert.rust_wood_insert_data().unwrap();
    let compound = Solid::compound_named(
        &[&insert, &insert],
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("occurrences.step");
    compound.write_step(&path).unwrap();
    let multiple = std::fs::read_to_string(&path).unwrap();
    let error = insert
        .verify_wood_insert_step_source(&multiple, 1e-7, None)
        .unwrap()
        .unwrap_err();
    assert_eq!(error.kind, GeometryErrorKind::InvalidResult);
    assert_eq!(
        error.message,
        "Insert candidate must have one actual occurrence"
    );
    insert.write_step(&path).unwrap();
    let single = std::fs::read_to_string(path).unwrap();
    insert
        .verify_wood_insert_step_source(&single, 1e-7, None)
        .unwrap()
        .unwrap();
    assert_eq!(insert.rust_wood_insert_data().unwrap(), before);
    let sphere = Solid::rust_sphere_named(
        1.,
        PlaneFrame::default(),
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(
        sphere
            .verify_wood_insert_step_source(&single, 1e-7, None)
            .is_none()
    );
}
