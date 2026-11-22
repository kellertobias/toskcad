use super::*;
use crate::step_read::owned_insert::DocumentOptions;
#[test]
fn same_selected_occupied_capture_matches_public_retain_and_rejects_cross_proof() {
    for quarter in [false, true] {
        let text = super::super::super::tests::fixture(quarter);
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let item = document.represented_members()[0].item_entity();
        let selected = Selected::parse(&document, item, Default::default(), None).unwrap();
        let occupied = selected
            .parsed()
            .projection(Default::default(), None)
            .unwrap()
            .occupied(Default::default(), None)
            .unwrap();
        let reused = selected
            .retain_occupied(&occupied, 64_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        let original = selected.retain(Default::default(), None).unwrap();
        assert_eq!(reused.work_used(), original.work_used() + 128);
        assert_eq!(reused.local_volume_mm3(), original.local_volume_mm3());
        assert_eq!(reused.world_volume_mm3(), original.world_volume_mm3());
        assert_eq!(reused.local_cap_areas_mm2(), original.local_cap_areas_mm2());
        assert_eq!(reused.local_enclosure_mm(), original.local_enclosure_mm());
        assert_eq!(
            reused.linear_columns_bounds(),
            original.linear_columns_bounds()
        );
        assert_eq!(
            reused.original_mapping_identity(),
            original.original_mapping_identity()
        );
        for i in 0..original.point_count() {
            let a = original.point(i).unwrap();
            let z = reused.point(i).unwrap();
            assert_eq!(a.entity(), z.entity());
            for (a, z) in a.coordinates().into_iter().zip(z.coordinates()) {
                assert_eq!(a.spelling(), z.spelling());
                assert_eq!(a.bounds(), z.bounds());
            }
        }
        for id in selected.entity_ids() {
            for component in 0..document.records()[id].len() {
                assert_eq!(
                    reused.original_component(*id, component),
                    original.original_component(*id, component)
                );
            }
        }
        let other = Selected::parse(&document, item, Default::default(), None).unwrap();
        assert!(matches!(
            other.retain_occupied(&occupied, 64_000_000, 64 * 1024 * 1024, None),
            Err(Error::Invalid)
        ));
    }
}
#[test]
fn same_selected_capture_work_storage_cancel_frontiers_are_atomic() {
    use std::cell::Cell;
    let text = super::super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let item = document.represented_members()[0].item_entity();
    let selected = Selected::parse(&document, item, Default::default(), None).unwrap();
    let occupied = selected
        .parsed()
        .projection(Default::default(), None)
        .unwrap()
        .occupied(Default::default(), None)
        .unwrap();
    let calls = Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let good = selected
        .retain_occupied(&occupied, 64_000_000, 64 * 1024 * 1024, Some(&count))
        .unwrap();
    let total = calls.get();
    for work in [good.work_used() - 1, good.work_used()] {
        let r = selected.retain_occupied(&occupied, work, 64 * 1024 * 1024, None);
        if work < good.work_used() {
            assert!(matches!(r, Err(Error::Limit)));
        } else {
            assert_eq!(r.unwrap().local_volume_mm3(), good.local_volume_mm3());
        }
    }
    for bytes in [
        good.construction_query_peak_bytes() - 1,
        good.construction_query_peak_bytes(),
    ] {
        let r = selected.retain_occupied(&occupied, 64_000_000, bytes, None);
        if bytes < good.construction_query_peak_bytes() {
            assert!(matches!(r, Err(Error::Limit)));
        } else {
            assert_eq!(
                r.unwrap().retained_budget_bytes(),
                good.retained_budget_bytes()
            );
        }
    }
    for at in [1, total / 2, total] {
        let calls = Cell::new(0usize);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= at
        };
        assert!(matches!(
            selected.retain_occupied(&occupied, 64_000_000, 64 * 1024 * 1024, Some(&cancel)),
            Err(Error::Cancelled)
        ));
        assert_eq!(
            selected
                .retain_occupied(&occupied, 64_000_000, 64 * 1024 * 1024, None)
                .unwrap()
                .local_volume_mm3(),
            good.local_volume_mm3()
        );
    }
    assert!(matches!(
        selected.retain_occupied(&occupied, 128, 64 * 1024 * 1024, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        selected.retain_occupied(
            &occupied,
            64_000_000,
            occupied.query_budget_bytes() - 1,
            None
        ),
        Err(Error::Limit)
    ));
}
