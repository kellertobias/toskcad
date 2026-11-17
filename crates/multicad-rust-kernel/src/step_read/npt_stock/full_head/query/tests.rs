use super::*;
use std::{cell::Cell, sync::Arc};
fn source(quarter: bool) -> OwnedActualHead {
    OwnedActualHead::read(
        Arc::from(super::super::source_tests::candidate(quarter, 18.)),
        Default::default(),
        None,
    )
    .unwrap()
}
fn point(p: [f64; 3]) -> [[f64; 2]; 3] {
    p.map(|v| [v, v])
}
#[test]
fn actual_sewn_union_queries_head_thread_air_and_contact_with_retained_source() {
    for quarter in [false, true] {
        let source = source(quarter);
        let [lo, hi] = source.head_stations_mm();
        let options = HeadQueryOptions::default();
        for (box_, expected) in [
            (point([0., 0., (lo + hi) / 2.]), Location::Inside),
            (point([0., 0., lo / 2.]), Location::Inside),
            (point([10., 8., (lo + hi) / 2.]), Location::Outside),
            (point([30., 0., lo / 2.]), Location::Outside),
            (point([0., 0., hi + 1.]), Location::Outside),
            (point([0., 0., -1.]), Location::Outside),
            (
                [[0., 0.], [0., 0.], [lo - 0.01, lo + 0.01]],
                Location::Unresolved,
            ),
        ] {
            let query = source.classify_local_box(box_, options, None).unwrap();
            assert!(std::ptr::eq(query.source(), &source));
            assert_eq!(query.local_bounds_mm(), box_);
            assert_eq!(query.location(), expected, "{box_:?}");
            assert!(query.work_used() > 0 && query.query_budget_bytes() <= MAX_BYTES);
        }
        // The air corner is INSIDE the head AABB, but outside actual hex halfspaces.
        let bounds = source.local_sewn_enclosure_mm();
        assert!(bounds[0][0] < 10. && bounds[0][1] > 10. && bounds[1][0] < 8. && bounds[1][1] > 8.);
    }
}
#[test]
fn actual_file_then_extra_native_inverse_queries_match_independent_physical_points() {
    for quarter in [false, true] {
        let source = Arc::new(source(quarter));
        let lo = source.head_stations_mm()[0];
        let intrinsic = LiteralOccurrencePlacement::new(
            [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
            [0., 0., lo],
            Default::default(),
            None,
        )
        .unwrap();
        let output =
            transport::encoding::propose_mapped(&source, &intrinsic, Default::default(), None)
                .unwrap();
        let mapped = transport::MappedActualHead::read_reusing_local(
            output.actual_bytes(),
            source,
            Default::default(),
            None,
        )
        .unwrap();
        let a = LiteralOccurrencePlacement::new(
            [[0.6, -0.8, 0.], [0.8, 0.6, 0.], [0., 0., 1.]],
            [20., -30., 11.],
            Default::default(),
            None,
        )
        .unwrap();
        for (local, expected) in [
            ([0., 0., lo + 1.], Location::Inside),
            ([0., 0., lo / 2.], Location::Inside),
            ([10., 8., lo + 1.], Location::Outside),
        ] {
            // Independent explicit original F, then A arithmetic, never reader's map.
            let world = [local[0], -local[1], lo - local[2]];
            let placed = [
                0.6 * world[0] - 0.8 * world[1] + 20.,
                0.8 * world[0] + 0.6 * world[1] - 30.,
                world[2] + 11.,
            ];
            let query = mapped
                .classify_box(point(world), Default::default(), None)
                .unwrap();
            assert_eq!(query.location(), expected);
            assert!(std::ptr::eq(query.source(), &mapped));
            assert!(query.occurrence().is_none());
            let query = mapped
                .classify_placed_box(point(placed), &a, Default::default(), None)
                .unwrap();
            assert_eq!(query.location(), expected);
            assert!(std::ptr::eq(query.occurrence().unwrap(), &a));
            assert!(query.work_used() > 256);
            let work = query.work_used();
            let bytes = query.query_budget_bytes();
            for (maximum_work, maximum_bytes, ok) in [
                (work - 1, MAX_BYTES, false),
                (work, MAX_BYTES, true),
                (QUERY_WORK, bytes - 1, false),
                (QUERY_WORK, bytes, true),
            ] {
                let result = mapped.classify_placed_box(
                    point(placed),
                    &a,
                    HeadQueryOptions {
                        query: QueryOptions {
                            maximum_work,
                            maximum_bytes,
                            ..Default::default()
                        },
                        ..Default::default()
                    },
                    None,
                );
                assert_eq!(result.is_ok(), ok);
                if !ok {
                    assert!(matches!(result, Err(Error::Limit)));
                }
            }
            let extra = mapped
                .classify_placed_box(
                    point(placed),
                    &a,
                    HeadQueryOptions {
                        additional_resident_bytes: 4096,
                        ..Default::default()
                    },
                    None,
                )
                .unwrap();
            assert_eq!(extra.query_budget_bytes(), bytes + 4096);
            assert_eq!(extra.work_used(), work);
            let calls = Cell::new(0usize);
            let observe = || {
                calls.set(calls.get() + 1);
                false
            };
            mapped
                .classify_placed_box(point(placed), &a, Default::default(), Some(&observe))
                .unwrap();
            let at = calls.get();
            let seen = Cell::new(0usize);
            let cancel = || {
                seen.set(seen.get() + 1);
                seen.get() == at
            };
            assert!(matches!(
                mapped.classify_placed_box(point(placed), &a, Default::default(), Some(&cancel)),
                Err(Error::Cancelled)
            ));
            assert_eq!(
                mapped
                    .classify_placed_box(point(placed), &a, Default::default(), None)
                    .unwrap()
                    .location(),
                expected
            );
        }
    }
}
#[test]
fn actual_head_query_work_bytes_cancel_and_unchanged_retry() {
    let source = source(false);
    let z = source.head_stations_mm()[0];
    for input in [point([0., 0., z + 1.]), point([0., 0., z / 2.])] {
        let good = source
            .classify_local_box(input, Default::default(), None)
            .unwrap();
        let work = good.work_used();
        let bytes = good.query_budget_bytes();
        let location = good.location();
        for (maximum_work, maximum_bytes, ok) in [
            (work - 1, MAX_BYTES, false),
            (work, MAX_BYTES, true),
            (QUERY_WORK, bytes - 1, false),
            (QUERY_WORK, bytes, true),
        ] {
            let result = source.classify_local_box(
                input,
                HeadQueryOptions {
                    query: QueryOptions {
                        maximum_work,
                        maximum_bytes,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                None,
            );
            assert_eq!(result.is_ok(), ok);
            if !ok {
                assert!(matches!(result, Err(Error::Limit)));
            }
        }
        let calls = Cell::new(0);
        let observe = || {
            calls.set(calls.get() + 1);
            false
        };
        source
            .classify_local_box(input, Default::default(), Some(&observe))
            .unwrap();
        for at in [1, calls.get() / 2, calls.get()] {
            let seen = Cell::new(0);
            let cancel = || {
                seen.set(seen.get() + 1);
                seen.get() >= at
            };
            assert!(matches!(
                source.classify_local_box(input, Default::default(), Some(&cancel)),
                Err(Error::Cancelled)
            ));
        }
        let retry = source
            .classify_local_box(input, Default::default(), None)
            .unwrap();
        assert_eq!(retry.location(), location);
        assert_eq!(retry.work_used(), work);
        let extra = source
            .classify_local_box(
                input,
                HeadQueryOptions {
                    additional_resident_bytes: 4096,
                    ..Default::default()
                },
                None,
            )
            .unwrap();
        assert_eq!(extra.query_budget_bytes(), bytes + 4096);
    }
    for input in [
        [[f64::NAN, 0.], [0., 0.], [0., 0.]],
        [[1., 0.], [0., 0.], [0., 0.]],
    ] {
        assert!(matches!(
            source.classify_local_box(input, Default::default(), None),
            Err(Error::Invalid)
        ));
    }
}
