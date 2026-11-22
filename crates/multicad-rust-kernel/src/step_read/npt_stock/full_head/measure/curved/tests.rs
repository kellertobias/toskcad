use super::*;
use std::cell::Cell;
fn source(quarter: bool) -> Arc<OwnedActualHead> {
    Arc::new(
        OwnedActualHead::read(
            Arc::from(super::super::super::source_tests::candidate(quarter, 18.)),
            Default::default(),
            None,
        )
        .unwrap(),
    )
}
#[test]
fn explicit_actual_plan_preserves_original_owner_and_curved_integrals_without_forcing_unselected_descriptors()
 {
    for quarter in [false, true] {
        let source = source(quarter);
        let before = Arc::strong_count(&source);
        let head = source
            .planar_face_measure(9, Default::default(), None)
            .unwrap();
        assert_eq!(head.inherited_plan_work_used(), 0);
        assert_eq!(head.aggregate_work_used(), head.work_used());
        assert_eq!(Arc::strong_count(&source), before);
        let full = source
            .auxiliary_thread()
            .geometry_plan(Default::default(), None)
            .unwrap();
        for ordinal in 0..8 {
            let plan =
                HeadGeometryPlan::new(Arc::clone(&source), ordinal, Default::default(), None)
                    .unwrap();
            assert_eq!(plan.selected_face_ordinal(), ordinal);
            assert!(std::ptr::eq(plan.source(), source.as_ref()));
            assert!(plan.mapped_source().is_none());
            assert_eq!(Arc::strong_count(&source), before + 1);

            let value = plan
                .face_measure(ordinal, Default::default(), None)
                .unwrap_or_else(|e| panic!("quarter={quarter} ordinal={ordinal}: {e:?}"));
            eprintln!(
                "ACTUAL_SELECTED_CURVED quarter={quarter} ordinal={ordinal} prep={} fresh={} aggregate={} peak={}",
                plan.preparation_work_used(),
                value.work_used(),
                value.aggregate_work_used(),
                value.accounted_peak_bytes()
            );
            assert!(matches!(
                plan.face_measure((ordinal + 1) % 8, Default::default(), None),
                Err(Error::Unsupported)
            ));

            assert!(std::ptr::eq(value.source(), source.as_ref()));
            assert_eq!(
                value.inherited_plan_work_used(),
                plan.preparation_work_used()
            );
            assert_eq!(
                value.aggregate_work_used(),
                value.work_used() + plan.preparation_work_used()
            );
            assert!(
                value.aggregate_work_used() <= MEASURE_WORK
                    && value.accounted_peak_bytes() <= MAX_BYTES
            );

            let entity = source.physical_roles().thread_face_entities[ordinal];
            let index = source
                .auxiliary_thread()
                .faces()
                .iter()
                .position(|f| f.entity() == entity)
                .unwrap();
            // Independently rerun the old actual raw full-surface integral,
            // not an authored radius/nominal area or a cached descriptor.
            let raw = full
                .face_measure(
                    index,
                    FaceMeasureOptions {
                        area_width_mm2: 1e-6,
                        centroid_width_mm: 1e-7,
                        ..Default::default()
                    },
                    None,
                )
                .unwrap();
            let factor = source.axial_sewing_transfer().positive_axial_scale_bounds();
            let raw_area = raw.area_mm2();
            let area = value.area_mm2();
            assert!(
                area[0] <= raw_area[1] * factor[1].max(1.)
                    && area[1] >= raw_area[0] * factor[0].min(1.)
            );
            let center = source
                .axial_sewing_transfer()
                .map_raw_box_mm(raw.centroid_mm())
                .unwrap();
            for (got, old) in value.centroid_mm().iter().zip(center) {
                assert!(got[0] <= old[1] && got[1] >= old[0]);
            }
        }
        assert_eq!(Arc::strong_count(&source), before);
        assert!(matches!(
            HeadGeometryPlan::new(Arc::clone(&source), 8, Default::default(), None),
            Err(Error::Unsupported)
        ));
        assert!(matches!(
            HeadGeometryPlan::new(Arc::clone(&source), 17, Default::default(), None),
            Err(Error::Invalid)
        ));
    }
}
#[test]
fn actual_plan_and_curved_measure_resource_cancel_retry_and_mapped_owner_are_bounded() {
    let source = source(false);
    let plan = HeadGeometryPlan::new(Arc::clone(&source), 0, Default::default(), None).unwrap();
    let work = plan.preparation_work_used();
    let bytes = plan.construction_peak_bytes();
    for (maximum_work, maximum_bytes, ok) in [
        (work - 1, MAX_BYTES, false),
        (work, MAX_BYTES, true),
        (MEASURE_WORK, bytes - 1, false),
        (MEASURE_WORK, bytes, true),
    ] {
        let result = HeadGeometryPlan::new(
            Arc::clone(&source),
            0,
            HeadGeometryOptions {
                plan: LiteralGeometryOptions {
                    maximum_work,
                    maximum_bytes,
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
    HeadGeometryPlan::new(Arc::clone(&source), 0, Default::default(), Some(&observe)).unwrap();
    for at in [1, calls.get()] {
        let seen = Cell::new(0);
        let cancel = || {
            seen.set(seen.get() + 1);
            seen.get() >= at
        };
        assert!(matches!(
            HeadGeometryPlan::new(Arc::clone(&source), 0, Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    let other = HeadGeometryPlan::new(
        Arc::clone(&source),
        1,
        HeadGeometryOptions {
            additional_resident_bytes: plan.retained_budget_bytes()
                - source.receipt().retained_bytes,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    let other_extra = other.retained_budget_bytes() - source.receipt().retained_bytes;
    let good = plan.face_measure(0, Default::default(), None).unwrap();
    let work = good.aggregate_work_used();
    let bytes = good.accounted_peak_bytes();
    let area = good.area_mm2();
    let shared = plan
        .face_measure(
            0,
            HeadMeasureOptions {
                additional_resident_bytes: other_extra,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    assert_eq!(shared.accounted_peak_bytes(), bytes + other_extra);

    for (maximum_work, maximum_bytes, ok) in [
        (work - 1, MAX_BYTES, false),
        (work, MAX_BYTES, true),
        (MEASURE_WORK, bytes - 1, false),
        (MEASURE_WORK, bytes, true),
    ] {
        let result = plan.face_measure(
            0,
            HeadMeasureOptions {
                measure: FaceMeasureOptions {
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
    plan.face_measure(0, Default::default(), Some(&observe))
        .unwrap();
    for at in [1, calls.get()] {
        let seen = Cell::new(0);
        let cancel = || {
            seen.set(seen.get() + 1);
            seen.get() >= at
        };
        assert!(matches!(
            plan.face_measure(0, Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(
        plan.face_measure(0, Default::default(), None)
            .unwrap()
            .area_mm2(),
        area
    );
    assert_eq!(
        plan.face_measure(
            0,
            HeadMeasureOptions {
                additional_resident_bytes: 4096,
                ..Default::default()
            },
            None
        )
        .unwrap()
        .accounted_peak_bytes(),
        bytes + 4096
    );
    let intrinsic = LiteralOccurrencePlacement::new(
        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
        [0., 0., source.head_stations_mm()[0]],
        Default::default(),
        None,
    )
    .unwrap();
    let output =
        transport::encoding::propose_mapped(&source, &intrinsic, Default::default(), None).unwrap();
    let mapped = Arc::new(
        transport::MappedActualHead::read_reusing_local(
            output.actual_bytes(),
            Arc::clone(&source),
            Default::default(),
            None,
        )
        .unwrap(),
    );
    let before = Arc::strong_count(&mapped);
    let mapped_plan =
        HeadGeometryPlan::new_mapped(Arc::clone(&mapped), 0, Default::default(), None).unwrap();
    assert!(std::ptr::eq(
        mapped_plan.mapped_source().unwrap(),
        mapped.as_ref()
    ));
    assert!(std::ptr::eq(mapped_plan.source(), source.as_ref()));
    assert_eq!(Arc::strong_count(&mapped), before + 1);
    let m = mapped_plan
        .face_measure(0, Default::default(), None)
        .unwrap();
    assert_eq!(m.area_mm2(), area);
    assert!(mapped_plan.retained_budget_bytes() > plan.retained_budget_bytes());
    drop(mapped_plan);
    assert_eq!(Arc::strong_count(&mapped), before);
}
