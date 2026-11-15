use super::*;
pub(super) fn owned(left: bool, plain: f64) -> OwnedPeriodicMaterial {
    let source = crate::step_read::phillips_cutter::physical::tests::fixture(0, left, true, 1);
    let mut spec = source.spec();
    spec.shank_length_mm = plain;
    let profile =
        crate::head_profile::Profile::new_for_thread_union(spec.head.unwrap(), 1., plain).unwrap();
    OwnedPeriodicMaterial::from_authored_parts(
        source.thread(),
        &profile,
        source.cavity().source(),
        spec,
        Default::default(),
        None,
    )
    .unwrap()
}
fn around(p: [f64; 3], half: f64) -> [[f64; 2]; 3] {
    p.map(|v| {
        [
            (I::point(v) - I::point(half)).lo,
            (I::point(v) + I::point(half)).hi,
        ]
    })
}
#[test]
fn actual_periodic_owner_queries_material_void_surfaces_and_real_head_shift() {
    for left in [false, true] {
        for plain in [0., 0.03, 0.3] {
            let source = owned(left, plain);
            let query = |p| {
                source
                    .classify_local_box(around(p, 1e-6), Default::default(), None)
                    .unwrap()
                    .location()
            };
            assert_eq!(query([0., 0., 1.]), C::Inside);
            assert_eq!(query([1.3, 0., -0.5 - plain]), C::Inside);
            let mouth = expression(source.cavity_mouth_expression_mm()).midpoint();
            assert_eq!(query([0., 0., mouth + 0.2]), C::Outside);
            assert_eq!(query([3., 0., 1.]), C::Outside);
            assert_eq!(query([0.9, 0.9, 1.]), C::Outside);
            assert_eq!(query([0., 0., 3.]), C::Outside);
            assert_eq!(query([0., 0., 2.]), C::UnresolvedBoundary);
            let wide = source
                .classify_local_box([[-2., 2.], [-2., 2.], [0.5, 0.6]], Default::default(), None)
                .unwrap();
            assert_eq!(wide.location(), C::UnresolvedBoundary);
            let huge = source
                .classify_local_box([[-f64::MAX, f64::MAX]; 3], Default::default(), None)
                .unwrap();
            assert_eq!(huge.location(), C::UnresolvedBoundary);
            assert!(std::ptr::eq(wide.source(), &source));
            assert_eq!(wide.bounds_mm(), [[-2., 2.], [-2., 2.], [0.5, 0.6]]);
            assert_eq!(
                wide.work_used() + wide.remaining_work(),
                Options::default().maximum_work
            );
        }
    }
}
#[test]
fn all_actual_periodic_radial_bands_classify_whole_nearby_boxes_in_both_hands() {
    for left in [false, true] {
        let source = owned(left, 0.);
        for band in 0..4 {
            for u in [0.25, 0.75] {
                let b = source
                    .thread()
                    .bounds(band, [u, u], [1., 1.], 16_000_000, None)
                    .unwrap();
                let point = b.position_mm.map(|v| interval(v).midpoint());
                for (factor, expected) in [
                    (0.98, C::Inside),
                    (1.02, C::Outside),
                    (1., C::UnresolvedBoundary),
                ] {
                    let p = [point[0] * factor, point[1] * factor, 1.];
                    let query = source
                        .classify_local_box(around(p, 1e-7), Default::default(), None)
                        .unwrap();
                    assert_eq!(
                        query.location(),
                        expected,
                        "left{left} band{band} u{u} p{p:?}"
                    );
                }
            }
        }
    }
}
#[test]
fn actual_periodic_material_query_frontiers_keep_one_balance_and_cancel_each_handoff() {
    let source = owned(false, 0.);
    let b = source
        .thread()
        .bounds(0, [0.5, 0.5], [1., 1.], 16_000_000, None)
        .unwrap();
    let p = b.position_mm.map(|v| interval(v).midpoint());
    let box_ = around([p[0] * 0.98, p[1] * 0.98, 1.], 1e-7);
    let o = Options::default();
    let good = source.classify_local_box(box_, o, None).unwrap();
    for bad in [
        Options {
            maximum_work: MAX + 1,
            ..o
        },
        Options {
            maximum_head_cells: 0,
            ..o
        },
        Options {
            maximum_head_cells: 1_000_001,
            ..o
        },
    ] {
        assert!(matches!(
            source.classify_local_box(box_, bad, None),
            Err(E::Invalid)
        ));
    }
    for bad in [0, 64, 4096, good.work_used() - 1] {
        assert!(matches!(
            source.classify_local_box(
                box_,
                Options {
                    maximum_work: bad,
                    ..o
                },
                None
            ),
            Err(E::Limit)
                | Err(E::Thread(crate::tapered_thread::Error::Limit))
                | Err(E::Head(crate::head_profile::RepairedHeadError::Limit))
        ));
    }
    assert!(
        source
            .classify_local_box(
                box_,
                Options {
                    maximum_work: good.work_used(),
                    ..o
                },
                None
            )
            .is_ok()
    );
    assert!(matches!(
        source.classify_local_box(
            around([1.3, 0., -0.5], 1e-7),
            Options {
                maximum_head_cells: 1,
                ..o
            },
            None
        ),
        Err(E::Head(crate::head_profile::RepairedHeadError::Limit))
    ));
    let mut bad = box_;
    bad[0] = [f64::NAN, 0.];
    assert!(matches!(
        source.classify_local_box(bad, o, None),
        Err(E::Invalid)
    ));
    for bounds in [
        box_,
        around([1.3, 0., -0.5], 1e-7),
        around([0., 0., -1.4], 1e-7),
        around([3., 0., 1.], 1e-7),
        around([0., 0., 2.], 1e-7),
    ] {
        let calls = Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        source.classify_local_box(bounds, o, Some(&count)).unwrap();
        let complete = calls.get();
        assert!(complete >= 2);
        for stop in 1..=complete {
            calls.set(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() == stop
            };
            let result = source.classify_local_box(bounds, o, Some(&cancel));
            assert!(
                matches!(
                    result,
                    Err(E::Cancelled)
                        | Err(E::Head(crate::head_profile::RepairedHeadError::Cancelled))
                ),
                "callback{stop}/{complete}"
            );
        }
    }
    assert!(source.classify_local_box(box_, o, None).is_ok());
}
#[test]
fn standalone_periodic_numeric_query_checks_caps_and_conservatively_handles_wide_phases() {
    let owner = owned(false, 0.);
    let source = owner.thread();
    let caps = owner.cap_expressions_mm();
    let p = [I::new(0.89, 0.91), I::new(-0.01, 0.01), I::new(0.5, 1.5)];
    assert_eq!(
        thread_location(source, caps, p, THREAD_WORK, None)
            .unwrap()
            .0,
        C::UnresolvedBoundary
    );
    for bad in [
        [[0., 0.], [0., 0.]],
        [[-1., 0.], [2., 0.]],
        [[0., 0.], [3., 0.]],
        [[f64::NAN, 0.], [2., 0.]],
    ] {
        assert!(matches!(
            thread_location(source, bad, p, THREAD_WORK, None),
            Err(E::Thread(crate::tapered_thread::Error::Invalid))
        ));
    }
    assert!(matches!(
        thread_location(source, caps, p, THREAD_WORK + 1, None),
        Err(E::Thread(crate::tapered_thread::Error::Invalid))
    ));
    let mut bad = p;
    bad[0] = I::point(f64::NAN);
    assert!(matches!(
        thread_location(source, caps, bad, THREAD_WORK, None),
        Err(E::Thread(crate::tapered_thread::Error::Invalid))
    ));
    for e in [
        PhillipsReadError::Invalid,
        PhillipsReadError::Limit,
        PhillipsReadError::Cancelled,
        PhillipsReadError::Unsupported,
    ] {
        let mapped = thread_error(e);
        assert!(matches!(mapped, E::Thread(_) | E::Cancelled));
    }
}
