use super::*;

fn source() -> RadialBoundary {
    BoundaryOwnership::new(crate::radial_tensor::tests::cam::wedge(&[1, 2]), 1e-7, None)
        .unwrap()
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap()
}

#[test]
fn proper_source_pose_requalifies_actual_values_and_preserves_atomic_cancellation() {
    use std::cell::Cell;
    let source = source();
    let original = source.input().clone();
    let frame = [[100., -30., 50.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]];
    let calls = Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    let placed = source
        .transformed_with_interrupt(frame, 1e-7, 1e-6, Some(&observe))
        .unwrap();
    assert_eq!(source.input(), &original);
    assert_ne!(placed.input(), &original);
    let n = [0.48, -0.64, 0.6];
    for (p, q) in original.vertices.iter().zip(&placed.input().vertices) {
        let expected = std::array::from_fn::<_, 3, _>(|k| {
            frame[0][k] + frame[1][k] * p[0] + frame[2][k] * p[1] + n[k] * p[2]
        });
        for k in 0..3 {
            assert!((expected[k] - q[k]).abs() < 1e-12);
        }
    }
    let total = calls.get();
    calls.set(0);
    let mut stages = vec![1];
    let owned = BoundaryOwnership::new(placed.input().clone(), 1e-7, Some(&observe)).unwrap();
    stages.push(calls.get() + 1);
    let trimmed = owned.bind_surfaces(1e-7, Some(&observe)).unwrap();
    stages.push(calls.get() + 1);
    let axial = trimmed.certify_axial(1e-7, 1e-6, Some(&observe)).unwrap();
    stages.push(calls.get() + 1);
    let caps = axial.certify_caps(Some(&observe)).unwrap();
    stages.push(calls.get() + 1);
    let ordered = caps.certify_ordering(Some(&observe)).unwrap();
    stages.push(calls.get() + 1);
    ordered.certify_material(1e-7, Some(&observe)).unwrap();
    let preparation = total - calls.get();
    // Every new source-copy checkpoint, then the first checkpoint in each
    // existing cold proof stage. Deep stage cancellation has its own tests.
    let checkpoints = (1..=preparation)
        .chain(stages.into_iter().map(|at| preparation + at))
        .chain([total]);
    for at in checkpoints {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= at
        };
        assert_eq!(
            source
                .transformed_with_interrupt(frame, 1e-7, 1e-6, Some(&cancel))
                .err(),
            Some(Error::Cancelled)
        );
        assert_eq!(source.input(), &original);
    }
    let text = crate::step::write_radial_boundary(&placed, 1e-7).unwrap();
    let restored = crate::step_read::read_radial_boundary(&text, 1e-7)
        .unwrap()
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    let mesh = restored
        .tessellate_with_interrupt(
            crate::display::Tessellation {
                deflection: 0.05,
                angle: 0.4,
                triangle_limit: 50000,
            },
            None,
        )
        .unwrap();
    assert!(crate::manufacturing::certify(&mesh).is_ok());
}

#[test]
fn invalid_or_unresolved_source_poses_never_replace_the_original_boundary() {
    let source = source();
    let original = source.input().clone();
    let identity = [[0.; 3], [1., 0., 0.], [0., 1., 0.]];
    for (tolerance, minimum) in [
        (0., 1e-6),
        (f64::NAN, 1e-6),
        (1e-7, 0.),
        (1e-7, f64::INFINITY),
    ] {
        assert_eq!(
            source
                .transformed_with_interrupt(identity, tolerance, minimum, None)
                .err(),
            Some(Error::Invalid)
        );
    }
    for frame in [
        [[f64::NAN, 0., 0.], identity[1], identity[2]],
        [identity[0], [f64::NAN, 0., 0.], identity[2]],
        [identity[0], [1.01, 0., 0.], identity[2]],
        [identity[0], identity[1], [0., 1.01, 0.]],
        [identity[0], identity[1], identity[1]],
        [[1e7, 0., 0.], identity[1], identity[2]],
    ] {
        assert_eq!(
            source
                .transformed_with_interrupt(frame, 1e-7, 1e-6, None)
                .err(),
            Some(Error::Invalid)
        );
    }
    assert!(
        source
            .transformed_with_interrupt(identity, 1e-20, 1e-6, None)
            .is_err()
    );
    assert!(
        source
            .transformed_with_interrupt(identity, 1e-7, 1e7, None)
            .is_err()
    );
    assert_eq!(source.input(), &original);
    assert!(
        source
            .transformed_with_interrupt(identity, 1e-7, 1e-6, None)
            .is_ok()
    );
}
