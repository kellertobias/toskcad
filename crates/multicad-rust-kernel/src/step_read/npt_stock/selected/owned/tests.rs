use super::*;
use crate::step_read::owned_insert::DocumentOptions;
#[test]
fn actual_owned_literal_geometry_survives_original_document_and_input_destruction() {
    for (quarter, mapped) in [(false, false), (true, false), (false, true)] {
        let owned = {
            let text = super::super::super::tests::fixture(quarter);
            let text = if mapped {
                super::super::tests::mapped(&text, "(-3.,0.,0.)", "1.")
            } else {
                text
            };
            let text = if mapped {
                text.replace("CARTESIAN_POINT('',", "CARTESIAN_POINT('測点',")
            } else {
                text
            };
            let document =
                CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
            let selected = Selected::parse(
                &document,
                document.represented_members()[0].item_entity(),
                Default::default(),
                None,
            )
            .unwrap();
            let owned = selected.retain(Default::default(), None).unwrap();
            assert_eq!(owned.source_entity_count(), selected.entity_ids().len());
            assert_eq!(owned.point_count(), selected.parsed.points.len());
            for (index, original) in selected.parsed.points.iter().enumerate() {
                let copy = owned.point(index).unwrap();
                assert_eq!(copy.id, original.id);
                for (copy, original) in copy.coordinates.into_iter().zip(original.coordinates) {
                    assert_eq!(copy.spelling(), original.spelling());
                    assert_eq!(copy.value(), original.value());
                    assert_eq!(copy.bounds(), original.bounds());
                }
            }
            for id in selected.entity_ids() {
                for (component, original) in document.records()[id].iter().enumerate() {
                    assert_eq!(owned.original_component(*id, component), original.source);
                }
            }
            eprintln!(
                "NPT OWNED quarter={quarter} mapped={mapped} work={} retained={} querypeak={} docbytes={}",
                owned.work_used(),
                owned.retained_budget_bytes(),
                owned.construction_query_peak_bytes(),
                owned.original_document_budget_bytes()
            );
            owned
        };
        let clone = owned.clone();
        drop(owned);
        assert_eq!(clone.vertices().len(), 12);
        assert_eq!(clone.edges().len(), 20);
        assert_eq!(clone.faces().len(), 10);
        assert_eq!(clone.curves().len(), 20);
        assert_eq!(clone.tensors().len(), 8);
        assert!(clone.source_entity_count() > clone.point_count());
        assert!(clone.point(usize::MAX).is_none());
        assert!(clone.original_component(u32::MAX, 0).is_none());
        assert!(
            clone
                .original_component(clone.point(0).unwrap().id, usize::MAX)
                .is_none()
        );
        assert_eq!(clone.member().mapping_link().is_some(), mapped);
        assert!(clone.local_volume_mm3()[0] > 0.);
        assert!(
            clone.world_volume_mm3()[0] <= clone.local_volume_mm3()[1]
                && clone.world_volume_mm3()[1] >= clone.local_volume_mm3()[0]
        );
        assert!(clone.local_cap_areas_mm2().iter().all(|a| a[0] > 0.));
        assert!(clone.determinant_bounds()[0] <= 1. && clone.determinant_bounds()[1] >= 1.);
        assert!(clone.retained_budget_bytes() < clone.construction_query_peak_bytes());
        assert!(clone.original_document_budget_bytes() > clone.retained_budget_bytes());
        assert!(
            clone
                .linear_columns_bounds()
                .iter()
                .flatten()
                .all(|i| i[0] <= i[1])
        );
        assert!(
            clone
                .inverse_columns_bounds()
                .iter()
                .flatten()
                .all(|i| i[0] <= i[1])
        );
        assert!(
            clone
                .origin_mm()
                .iter()
                .chain(&clone.inverse_origin_mm())
                .all(|i| i[0] <= i[1])
        );
        assert!(
            clone
                .local_enclosure_mm()
                .iter()
                .chain(&clone.world_enclosure_mm())
                .all(|i| i[0] <= i[1])
        );
        let p = [1., 2., 3.];
        let world = if mapped { [9., -4., 6.] } else { p };
        super::super::tests::contains_point(clone.map_box_mm(p.map(|x| [x, x])).unwrap(), world);
        super::super::tests::contains_point(
            clone.inverse_box_mm(world.map(|x| [x, x])).unwrap(),
            p,
        );
        let cell = clone.cell(0, [0, 0]).unwrap();
        assert!(cell.bounds().unwrap().iter().all(|b| b[0] <= b[1]));
        assert!(matches!(
            clone.cell(usize::MAX, [0, 0]),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            clone.cell(0, [usize::MAX, 0]),
            Err(Error::Invalid)
        ));
    }
}
#[test]
fn owned_copy_has_one_aggregate_work_ledger_and_atomic_storage_cancel_frontiers() {
    let text = super::super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let selected = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let standard = MaterialOptions::default();
    for options in [
        MaterialOptions {
            maximum_work: 0,
            ..standard
        },
        MaterialOptions {
            maximum_work: 64_000_001,
            ..standard
        },
        MaterialOptions {
            maximum_bytes: 0,
            ..standard
        },
        MaterialOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..standard
        },
        MaterialOptions {
            maximum_nodes: 65_537,
            ..standard
        },
        MaterialOptions {
            maximum_depth: 33,
            ..standard
        },
    ] {
        assert!(matches!(
            selected.retain(options, None),
            Err(Error::Invalid)
        ));
    }
    let baseline = selected.retain(standard, None).unwrap();
    let exact = selected
        .retain(
            MaterialOptions {
                maximum_work: baseline.work_used(),
                maximum_bytes: baseline.construction_query_peak_bytes(),
                ..standard
            },
            None,
        )
        .unwrap();
    assert_eq!(exact.work_used(), baseline.work_used());
    assert_eq!(
        exact.retained_budget_bytes(),
        baseline.retained_budget_bytes()
    );
    assert!(matches!(
        selected.retain(
            MaterialOptions {
                maximum_work: baseline.work_used() - 1,
                ..standard
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        selected.retain(
            MaterialOptions {
                maximum_bytes: baseline.construction_query_peak_bytes() - 1,
                ..standard
            },
            None
        ),
        Err(Error::Limit)
    ));
    let occupied = selected
        .parsed
        .projection(standard, None)
        .unwrap()
        .occupied(standard, None)
        .unwrap();
    let fixed = size_of::<OwnedMaterial>()
        + size_of::<Source>()
        + 2 * size_of::<usize>()
        + 2 * std::mem::align_of::<usize>();
    let at_fixed = occupied.query_budget_bytes() + size_of::<Placement<'_, '_, '_>>() + fixed - 1;
    assert!(matches!(
        selected.retain(
            MaterialOptions {
                maximum_bytes: at_fixed,
                ..standard
            },
            None
        ),
        Err(Error::Limit)
    ));
    let calls = std::cell::Cell::new(0usize);
    selected
        .retain(
            standard,
            Some(&|| {
                calls.set(calls.get() + 1);
                false
            }),
        )
        .unwrap();
    for checkpoint in [1, calls.get() / 2, calls.get() - 1] {
        let count = std::cell::Cell::new(0usize);
        let stop = || {
            count.set(count.get() + 1);
            count.get() >= checkpoint
        };
        assert!(matches!(
            selected.retain(standard, Some(&stop)),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(selected.parsed.document.literal_input(), text);
    assert!(selected.retain(standard, None).is_ok());
}
#[test]
fn actual_mapped_world_volume_must_remain_finite_after_literal_retention() {
    let local = super::super::super::tests::fixture(false);
    let text = super::super::tests::mapped(&local, "(-3.,0.,0.)", "1.E102");
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let selected = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    assert!(selected.placement(Default::default(), None).is_ok());
    assert!(
        selected
            .parsed
            .projection(Default::default(), None)
            .unwrap()
            .occupied(Default::default(), None)
            .is_ok()
    );
    assert!(matches!(
        selected.retain(Default::default(), None),
        Err(Error::Unsupported)
    ));
}
