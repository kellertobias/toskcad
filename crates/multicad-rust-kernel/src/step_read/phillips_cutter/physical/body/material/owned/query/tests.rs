use super::*;
use std::cell::Cell;
fn box_at(p: [f64; 3], half: f64) -> [[f64; 2]; 3] {
    p.map(|x| [x - half, x + half])
}
fn owned(scale: f64, posed: bool, left: bool) -> OwnedDerivedMaterial {
    super::super::mesh_probe::owned_neck(scale, posed, left, 2. * scale, 0.1 * scale)
}
fn cancelled(e: OwnedMaterialQueryError) -> bool {
    matches!(
        e,
        OwnedMaterialQueryError::Cancelled
            | OwnedMaterialQueryError::Head(RepairedHeadError::Cancelled)
            | OwnedMaterialQueryError::Thread(crate::tapered_thread::Error::Cancelled)
    )
}
#[test]
fn actual_owned_material_boxes_keep_cavity_shift_and_both_source_components() {
    for (scale, posed, left) in [(1., false, false), (1.137, true, true)] {
        let source = owned(scale, posed, left);
        let original = source.original_input().to_owned();
        let head = source.head().view();
        let top = head.station(0).unwrap();
        assert_eq!(top[1], source.chart().outer_top_z);
        let caps = source.thread().cap_axes();
        let thread_z = (caps[0].origin[2] + caps[1].origin[2]) * 0.5;
        let r = source.local_enclosure_bounds()[0][1];
        for (point, expected) in [
            ([0., 0., top[1] + 0.1 * scale], C::Outside),
            ([0.99 * top[0], 0., top[1] + 0.001 * scale], C::Inside),
            ([0., 0., thread_z], C::Inside),
            ([0.99 * r, 0.99 * r, thread_z], C::Outside),
            ([r + scale, 0., thread_z], C::Outside),
        ] {
            let bounds = box_at(point, 1e-7 * scale);
            let q = source
                .classify_local_box(bounds, Default::default(), None)
                .unwrap();
            assert_eq!(q.location(), expected, "scale={scale}, point={point:?}");
            assert!(std::ptr::eq(q.source(), &source));
            assert_eq!(q.bounds_mm(), bounds);
            assert_eq!(q.work_used() + q.remaining_work(), MAXIMUM_WORK);
            let tight = OwnedMaterialQueryOptions {
                maximum_work: q.work_used(),
                ..Default::default()
            };
            let exact = source.classify_local_box(bounds, tight, None).unwrap();
            assert_eq!(exact.location(), expected);
            assert_eq!(exact.remaining_work(), 0);
            let short = OwnedMaterialQueryOptions {
                maximum_work: q.work_used() - 1,
                ..tight
            };
            assert!(source.classify_local_box(bounds, short, None).is_err());
        }
        let whole = source.local_enclosure_bounds();
        assert_eq!(
            source
                .classify_local_box(whole, Default::default(), None)
                .unwrap()
                .location(),
            C::UnresolvedBoundary
        );
        // One endpoint is in the real cavity, the other in retained head
        // material. The complete head slab cannot be called Inside or Outside.
        let cut_crossing = [[0., 0.99 * top[0]], [0., 0.], [top[1] + 0.1 * scale; 2]];
        assert_eq!(
            source
                .classify_local_box(cut_crossing, Default::default(), None)
                .unwrap()
                .location(),
            C::UnresolvedBoundary
        );
        assert_eq!(source.original_input(), original);
    }
}
#[test]
fn actual_owned_box_resources_and_every_checkpoint_cancel_preserve_source() {
    let source = owned(1., false, false);
    let caps = source.thread().cap_axes();
    let bounds = box_at(
        [0., 0., (caps[0].origin[2] + caps[1].origin[2]) * 0.5],
        1e-7,
    );
    let original = source.original_input().to_owned();
    let polls = Cell::new(0);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    let baseline = source
        .classify_local_box(bounds, Default::default(), Some(&count))
        .unwrap();
    let checkpoints = polls.get();
    assert!(checkpoints > 4);
    for stop in 1..=checkpoints {
        polls.set(0);
        let cancel = || {
            polls.set(polls.get() + 1);
            polls.get() >= stop
        };
        let error = source
            .classify_local_box(bounds, Default::default(), Some(&cancel))
            .err()
            .unwrap();
        assert!(cancelled(error), "checkpoint={stop}: {error:?}");
    }
    for options in [
        OwnedMaterialQueryOptions {
            maximum_work: MAXIMUM_WORK + 1,
            ..Default::default()
        },
        OwnedMaterialQueryOptions {
            maximum_head_cells: 0,
            ..Default::default()
        },
        OwnedMaterialQueryOptions {
            maximum_head_cells: 1_000_001,
            ..Default::default()
        },
    ] {
        assert_eq!(
            source.classify_local_box(bounds, options, None).err(),
            Some(OwnedMaterialQueryError::Invalid)
        );
    }
    assert_eq!(
        source
            .classify_local_box(
                bounds,
                OwnedMaterialQueryOptions {
                    maximum_work: 0,
                    ..Default::default()
                },
                None
            )
            .err(),
        Some(OwnedMaterialQueryError::Limit)
    );
    // The actual thread-midpoint lies beyond the head's finite axial slab.
    // With no provider balance the head refuses; two further units finish
    // that head exclusion but leave no balance for the retained thread.
    for (maximum_work, expected) in [
        (
            4160,
            OwnedMaterialQueryError::Head(RepairedHeadError::Limit),
        ),
        (
            4162,
            OwnedMaterialQueryError::Thread(crate::tapered_thread::Error::Limit),
        ),
    ] {
        assert_eq!(
            source
                .classify_local_box(
                    bounds,
                    OwnedMaterialQueryOptions {
                        maximum_work,
                        ..Default::default()
                    },
                    None,
                )
                .err(),
            Some(expected)
        );
    }
    for bad in [[1., -1.], [f64::NAN, 0.], [0., f64::INFINITY]] {
        let mut invalid = bounds;
        invalid[1] = bad;
        assert_eq!(
            source
                .classify_local_box(invalid, Default::default(), None)
                .err(),
            Some(OwnedMaterialQueryError::Invalid)
        );
    }
    assert_eq!(
        source
            .classify_local_box(bounds, Default::default(), None)
            .unwrap()
            .location(),
        baseline.location()
    );
    assert_eq!(source.original_input(), original);
}
