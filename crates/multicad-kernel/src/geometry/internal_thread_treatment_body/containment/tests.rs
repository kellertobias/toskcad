use super::*;
use crate::rust_kernel::{helical_thread::BoxClassification, thread_end_relief};
fn male(left: bool, tail: bool) -> Solid {
    let request = super::super::tests::cut(PlaneFrame::default(), left, 1);
    Solid::rust_helical_thread_entry_chamfer(
        request.profile,
        None,
        0.4,
        tail.then_some(thread_end_relief::Intent {
            width_mm: 0.5,
            radial_depth_mm: 0.05,
        }),
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}
fn receiver(left: bool, end: Rational, entry: bool, tail: bool) -> Solid {
    let frame = PlaneFrame::default();
    let stock = super::super::tests::fixture_stock(end, frame);
    let mut request = super::super::tests::cut(frame, left, 1);
    if !entry {
        request.intent.entry_depth_mm = None;
    }
    if !tail {
        request.intent.terminal = None;
    }
    if end == Rational::new(4, 1).unwrap() {
        stock
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                request,
                EntityId::new(),
                GeometryPolicy::default(),
                None,
            )
            .unwrap()
    } else {
        stock
            .rust_cut_authored_internal_thread_treatments_with_exact_extent_and_interrupt(
                request,
                end,
                EntityId::new(),
                GeometryPolicy::default(),
                None,
            )
            .unwrap()
    }
}
fn box_at(center: [f64; 3], dims: [f64; 3]) -> Solid {
    Solid::box_solid(dims, GeometryPolicy::default())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: std::array::from_fn(|i| center[i] - dims[i] / 2.),
            ..Default::default()
        })
        .unwrap()
}
#[test]
fn whole_actual_entry_and_terminal_boxes_retain_unchanged_probe_volume_in_both_poses() {
    let pose = RigidTransform {
        rotation: [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
        translation_mm: [11., -7., 5.],
    };
    for left in [false, true] {
        let a = male(left, true);
        let b = receiver(left, Rational::new(4, 1).unwrap(), true, true);
        let original = (
            a.rust_thread_entry_chamfer_data().unwrap(),
            b.rust_internal_thread_treatment_data().unwrap(),
        );
        for (source, z, radius, inside) in [
            (&a, 0.01, 1.109, true),
            (&a, 0.01, 1.111, false),
            (&a, 3.8, 1.199, true),
            (&a, 3.8, 1.201, false),
            (&b, 0.01, 1.741, true),
            (&b, 0.01, 1.739, false),
            (&b, 3.8, 1.626, true),
            (&b, 3.8, 1.624, false),
            (&b, 3.8, 2.6, false),
        ] {
            let probe = box_at([radius, 0., z], [0.00002; 3]);
            let moved_source = source.transformed(pose).unwrap();
            let moved_probe = probe.transformed(pose).unwrap();
            for (source, probe) in [(source, &probe), (&moved_source, &moved_probe)] {
                for (first, second) in [(source, probe), (probe, source)] {
                    let result = first
                        .boolean(
                            second,
                            BooleanOperation::Intersect,
                            GeometryPolicy::default(),
                        )
                        .unwrap();
                    assert_eq!(result.is_some(), inside, "left={left} z={z} r={radius}");
                    if let Some(result) = result {
                        assert_eq!(result.volume_mm3().unwrap(), probe.volume_mm3().unwrap());
                        assert_eq!(result.bounds().unwrap(), probe.bounds().unwrap());
                    }
                }
            }
        }
        assert_eq!(
            (
                a.rust_thread_entry_chamfer_data().unwrap(),
                b.rust_internal_thread_treatment_data().unwrap()
            ),
            original
        );
    }
}
#[test]
fn actual_intersecting_boxes_cannot_be_decided_from_corners_or_nominal_endpoints() {
    use BoxClassification::*;
    let a = male(false, true);
    let b = receiver(false, Rational::new(4, 1).unwrap(), true, true);
    for (source, probe) in [
        // All corners are outside the cone, but the center is inside material.
        (&a, box_at([0., 0., 0.01], [4., 4., 0.00002])),
        // Corners AND center are outside receiver material, but its annulus
        // crosses this complete box. Neither set of samples proves Outside.
        (&b, box_at([0., 0., 0.01], [8., 8., 0.00002])),
        (&a, box_at([1.11, 0., 0.01], [0.02, 0.02, 0.00002])),
        (&b, box_at([1.625, 0., 3.8], [0.02, 0.02, 0.00002])),
        (&b, box_at([2., 0., 4.], [0.00002; 3])),
        // Entry no longer dominates the actual original helical field.
        (&a, box_at([1.2, 0., 0.2], [0.00002; 3])),
        (&b, box_at([1.4, 0., 0.4], [0.00002; 3])),
        (&a, box_at([1.2, 0., 1.], [0.00002; 3])),
    ] {
        assert!(
            source
                .boolean(
                    &probe,
                    BooleanOperation::Intersect,
                    GeometryPolicy::default()
                )
                .is_err()
        );
    }
    assert_eq!(
        a.rust_thread_entry_chamfer_data()
            .unwrap()
            .source
            .classify_treatment_box([[-2., 2.], [-2., 2.], [0.00999, 0.01001]])
            .unwrap(),
        UnresolvedBoundary
    );
    assert_eq!(
        b.rust_internal_thread_treatment_data()
            .unwrap()
            .material
            .classify_treatment_box([[-4., 4.], [-4., 4.], [0.00999, 0.01001]])
            .unwrap(),
        UnresolvedBoundary
    );
    let stock =
        super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), PlaneFrame::default());
    assert!(
        a.boolean(
            &stock,
            BooleanOperation::Intersect,
            GeometryPolicy::default()
        )
        .is_err()
    );
    let probe = box_at([1.109, 0., 0.01], [0.00002; 3]);
    assert!(
        a.boolean(&probe, BooleanOperation::Union, GeometryPolicy::default())
            .is_err()
    );
    assert!(
        a.boolean(
            &probe,
            BooleanOperation::Subtract,
            GeometryPolicy::default()
        )
        .is_err()
    );
}
#[test]
fn exact_extended_end_and_absent_treatments_preserve_actual_finite_domains() {
    use BoxClassification::*;
    let b = receiver(false, Rational::new(41, 10).unwrap(), true, true);
    // This occupied terminal cylinder is ABOVE the original4mm thread domain.
    let inside = box_at([1.626, 0., 4.05], [0.00002; 3]);
    assert_eq!(
        b.boolean(
            &inside,
            BooleanOperation::Intersect,
            GeometryPolicy::default()
        )
        .unwrap()
        .unwrap()
        .volume_mm3()
        .unwrap(),
        inside.volume_mm3().unwrap()
    );
    let beyond = box_at([1.626, 0., 4.11], [0.00002; 3]);
    assert!(
        b.boolean(
            &beyond,
            BooleanOperation::Intersect,
            GeometryPolicy::default()
        )
        .unwrap()
        .is_none()
    );
    let touching = box_at([1.626, 0., 4.1], [0.00002; 3]);
    assert!(
        b.boolean(
            &touching,
            BooleanOperation::Intersect,
            GeometryPolicy::default()
        )
        .is_err()
    );
    let a = male(false, false);
    assert_eq!(
        a.rust_thread_entry_chamfer_data()
            .unwrap()
            .source
            .classify_treatment_box([[1.19, 1.20], [0., 0.00001], [3.79, 3.80]])
            .unwrap(),
        UnresolvedBoundary
    );
    let entry_only = receiver(false, Rational::new(4, 1).unwrap(), true, false);
    assert_eq!(
        entry_only
            .rust_internal_thread_treatment_data()
            .unwrap()
            .material
            .classify_treatment_box([[1.62, 1.63], [0., 0.00001], [3.79, 3.80]])
            .unwrap(),
        UnresolvedBoundary
    );
    let tail_only = receiver(false, Rational::new(4, 1).unwrap(), false, true);
    assert_eq!(
        tail_only
            .rust_internal_thread_treatment_data()
            .unwrap()
            .material
            .classify_treatment_box([[1.7, 1.71], [0., 0.00001], [0.01, 0.02]])
            .unwrap(),
        UnresolvedBoundary
    );
    for bounds in [
        [[f64::NAN, 1.], [0., 1.], [0., 1.]],
        [[1., 0.], [0., 1.], [0., 1.]],
    ] {
        assert!(
            a.rust_thread_entry_chamfer_data()
                .unwrap()
                .source
                .classify_treatment_box(bounds)
                .is_err()
        );
        assert!(
            b.rust_internal_thread_treatment_data()
                .unwrap()
                .material
                .classify_treatment_box(bounds)
                .is_err()
        );
    }
}
