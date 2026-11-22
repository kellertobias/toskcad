use super::*;
#[test]
fn numeric_world_transfer_guard_refuses_overflow_before_publishing_a_bound() {
    // Pure numeric helper frontier, not an admitted Body or a file fixture.
    // Preserve the defensive refusal for future callers with different scales.
    let g = Factors {
        a: vec![[1e-8; 2]],
        au: vec![1.],
        b: [1e-8; 2],
        bv: 1.,
        z: [0.; 2],
        major: 1.,
    };
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let frame = [[0.; 3], [f64::MAX, 0., 0.], [0., 1., 0.]];
    assert!(world_transfer(0, &Target::raw(frame), [0.; 3], &g, &mut work).is_err());
    work.remaining = 0;
    assert!(world_transfer(0, &Target::raw([[0.; 3]; 3]), [0.; 3], &g, &mut work).is_err());
}

#[test]
fn untrusted_reference_grid_shape_contract_refuses_before_tensor_indexing() {
    for (w, h, coherent, ew, eh) in [
        (0, 12, true, 0, 12),
        (12, 0, true, 12, 0),
        (12, 12, false, 12, 12),
        (12, 12, true, 23, 12),
        (12, 12, true, 12, 23),
    ] {
        assert!(matches!(
            require_matching_grid(w, h, coherent, ew, eh),
            Err(ReadError::Unsupported)
        ));
    }
    assert!(require_matching_grid(12, 23, true, 12, 23).is_ok());
}

#[test]
fn genuine_independent_reference_and_out_of_range_index_refuse_before_pole_access() {
    let text = super::super::super::degree11::tests::actual_m3(true, false);
    super::super::super::super::common_mouth::tests::derived_thread_probe::with_head(
        &text,
        &mut |head| {
            let body = head.body();
            let mut shell = Shell::new(&body.entities, body.tolerance, true);
            let mut work = Work {
                remaining: WORK,
                interrupt: None,
                cancelled: Cell::new(false),
            };
            let raw = raw(head, &mut shell, &mut work, BYTES).unwrap();
            let reference = crate::step::helical::prepare_retained_reference(
                &body.thread,
                body.visible[0],
                body.tolerance,
            )
            .unwrap();
            let target = Target::raw(body.inferred_thread_frame());
            assert!(
                literal_transfer(
                    &raw[0],
                    &reference,
                    reference.bands.len(),
                    &target,
                    &shell,
                    &mut work
                )
                .is_err()
            );
            let mut spec = body.thread.spec();
            spec.length_mm *= 2.;
            let different = crate::helical_thread::Thread::new(spec).unwrap();
            let other = crate::step::helical::prepare_retained_reference(
                &different,
                body.visible[0],
                body.tolerance,
            )
            .unwrap();
            assert!(matches!(
                literal_transfer(&raw[0], &other, 0, &target, &shell, &mut work),
                Err(ReadError::Unsupported)
            ));
            assert!(
                literal_transfer(&raw[0], &reference, 0, &target, &shell, &mut work)
                    .unwrap()
                    .into_iter()
                    .all(|d| d <= 1e-7)
            );
        },
    );
}
