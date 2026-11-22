use super::*;
use crate::step_read::npt_stock::full_head::source_tests::with_source;
fn owner(quarter: bool) -> Arc<OwnedActualHead> {
    Arc::new(
        OwnedActualHead::read(
            Arc::from(super::super::tests::candidate(quarter)),
            Default::default(),
            None,
        )
        .unwrap(),
    )
}
fn pose(length: f64, rotated: bool) -> LiteralOccurrencePlacement {
    let p = LiteralOccurrencePlacement::new(
        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
        [0., 0., length],
        Default::default(),
        None,
    )
    .unwrap();
    if rotated {
        p.transformed(
            [[0.6, -0.8, 0.], [0.8, 0.6, 0.], [0., 0., 1.]],
            [2., -3., 0.],
            Default::default(),
            None,
        )
        .unwrap()
    } else {
        p
    }
}
fn intent() -> HexHeadIntent {
    HexHeadIntent {
        across_flats_mm: 18.,
        height_mm: 3.,
    }
}
#[test]
fn complete_verified_mapped_file_has_actual_material_original_source_frame_and_unchanged_cold_import()
 {
    for quarter in [false, true] {
        with_source(quarter, |source| {
            let owner = owner(quarter);
            for rotated in [false, true] {
                let p = pose(source.curves().source().source().spec().length_mm, rotated);
                let verified = encode_verified_mapped(
                    Arc::clone(&owner),
                    source,
                    intent(),
                    &p,
                    Default::default(),
                    None,
                )
                .unwrap();
                assert!(std::ptr::eq(verified.source_plan(), source));
                assert!(std::ptr::eq(verified.original_placement(), &p));
                assert_eq!(verified.head_definition().across_flats_mm, 18.);
                assert_eq!(verified.head_definition().height_mm, 3.);
                let r = verified.receipt();
                assert!(
                    r.aggregate_work <= 272_000_000 && r.accounted_peak_bytes <= 64 * 1024 * 1024
                );
                assert_eq!(
                    r.aggregate_work,
                    r.encoding.fresh_encoding_work + r.source.aggregate_work + r.orchestration_work
                );
                assert_eq!(
                    r.source.aggregate_work,
                    owner.receipt().aggregate_work
                        + r.source.actual_file.original_cold_work_reserved
                        + r.source.actual_file.envelope_work
                        + source.work_used()
                        + p.work_used()
                        + r.source.head_and_frame_work
                        + r.source.fresh_thread_correspondence_work
                );
                let cold = MappedActualHead::read(
                    Arc::from(verified.actual().actual_file()),
                    Default::default(),
                    None,
                )
                .unwrap();
                assert_eq!(
                    verified.actual().world_sewn_enclosure_mm(),
                    cold.world_sewn_enclosure_mm()
                );
                assert_eq!(
                    verified.actual().world_sewn_volume_mm3(),
                    cold.world_sewn_volume_mm3()
                );
                let actual = verified.into_actual();
                assert_eq!(
                    actual.local_geometry().head_polygon_mm(),
                    cold.local_geometry().head_polygon_mm()
                );
                eprintln!("VERIFIED_MAPPED quarter={quarter} rotated={rotated} {r:?}");
            }
        });
    }
}
#[test]
fn composed_mapped_work_bytes_cancel_failure_and_retry_are_cumulative_and_atomic() {
    use std::cell::Cell;
    with_source(false, |source| {
        let owner = owner(false);
        let p = pose(source.curves().source().source().spec().length_mm, true);
        let good = encode_verified_mapped(
            Arc::clone(&owner),
            source,
            intent(),
            &p,
            Default::default(),
            None,
        )
        .unwrap();
        let r = good.receipt();
        let text = good.actual().actual_file().to_owned();
        drop(good);
        let refs = Arc::strong_count(&owner);
        for (work, bytes, ok) in [
            (r.aggregate_work - 1, 64 * 1024 * 1024, false),
            (r.aggregate_work, 64 * 1024 * 1024, true),
            (272_000_000, r.accounted_peak_bytes - 1, false),
            (272_000_000, r.accounted_peak_bytes, true),
        ] {
            let value = encode_verified_mapped(
                Arc::clone(&owner),
                source,
                intent(),
                &p,
                SourceOptions {
                    resources: Options {
                        maximum_work: work,
                        maximum_bytes: bytes,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                None,
            );
            assert_eq!(value.is_ok(), ok);
            if !ok {
                assert!(matches!(value, Err(Error::Limit)));
            }
            drop(value);
            assert_eq!(Arc::strong_count(&owner), refs);
        }
        let n = Cell::new(0);
        let observe = || {
            n.set(n.get() + 1);
            false
        };
        encode_verified_mapped(
            Arc::clone(&owner),
            source,
            intent(),
            &p,
            Default::default(),
            Some(&observe),
        )
        .unwrap();
        for at in [1, n.get() / 2, n.get()] {
            let seen = Cell::new(0);
            let cancel = || {
                seen.set(seen.get() + 1);
                seen.get() >= at
            };
            assert!(matches!(
                encode_verified_mapped(
                    Arc::clone(&owner),
                    source,
                    intent(),
                    &p,
                    Default::default(),
                    Some(&cancel)
                ),
                Err(Error::Cancelled)
            ));
            assert_eq!(Arc::strong_count(&owner), refs);
        }
        assert!(matches!(
            encode_verified_mapped(
                Arc::clone(&owner),
                source,
                HexHeadIntent {
                    across_flats_mm: 20.,
                    height_mm: 3.
                },
                &p,
                Default::default(),
                None
            ),
            Err(Error::Unsupported)
        ));
        assert_eq!(Arc::strong_count(&owner), refs);
        let extra = encode_verified_mapped(
            Arc::clone(&owner),
            source,
            intent(),
            &p,
            SourceOptions {
                resources: Options {
                    additional_resident_bytes: 4096,
                    ..Default::default()
                },
                ..Default::default()
            },
            None,
        )
        .unwrap();
        assert_eq!(extra.receipt().aggregate_work, r.aggregate_work);
        assert_eq!(
            extra.receipt().accounted_peak_bytes,
            r.accounted_peak_bytes + 4096
        );
        drop(extra);
        let retry =
            encode_verified_mapped(owner, source, intent(), &p, Default::default(), None).unwrap();
        assert_eq!(retry.actual().actual_file(), text);
    });
}
