//! Backend-neutral exact text witnesses and their independent readback.
use multicad_kernel::{BooleanOperation, GeometryPolicy, RigidTransform, Solid};

#[test]
fn rust_exact_text_is_deterministic_and_reads_back_as_an_independent_solid() {
    let policy = GeometryPolicy::default();
    let solid = Solid::box_solid([10., 20., 30.], policy).unwrap();
    let text = solid.exact_text().unwrap();
    assert!(text.starts_with("ISO-10303-21;"));
    assert_eq!(solid.exact_text().unwrap(), text);
    let independent = Solid::from_exact_text(&text, policy).unwrap();
    assert!((independent.volume_mm3().unwrap() - 6000.).abs() < 1e-9);
    assert_eq!(independent.bounds().unwrap().max, [10., 20., 30.]);

    // A moved solid has different exact text.
    let moved = solid
        .transformed(RigidTransform {
            translation_mm: [1., 0., 0.],
            ..Default::default()
        })
        .unwrap();
    assert_ne!(moved.exact_text().unwrap(), text);
    // Malformed text is refused, never approximated.
    assert!(Solid::from_exact_text("ISO-10303-21;\nEND-ISO-10303-21;\n", policy).is_err());
}

#[test]
fn rust_parallel_intersection_matches_the_serial_exact_result() {
    let policy = GeometryPolicy::default();
    let a = Solid::box_solid([10.; 3], policy).unwrap();
    let b = a
        .transformed(RigidTransform {
            translation_mm: [5., 0., 0.],
            ..Default::default()
        })
        .unwrap();
    let serial = a
        .boolean(&b, BooleanOperation::Intersect, policy)
        .unwrap()
        .unwrap();
    let parallel = a.intersect_parallel(&b, policy).unwrap().unwrap();
    assert!((parallel.volume_mm3().unwrap() - 500.).abs() < 1e-9);
    assert_eq!(parallel.volume_mm3().unwrap(), serial.volume_mm3().unwrap());
}
