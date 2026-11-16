use super::*;
use crate::step_read::npt_stock::full_head::source_tests::with_source;
#[test]
fn actual_mapped_candidate_preserves_physical_records_then_independently_earns_source_and_frame() {
    for quarter in [false, true] {
        with_source(quarter, |plan| {
            let local = OwnedActualHead::read(
                Arc::from(super::super::tests::candidate(quarter)),
                Default::default(),
                None,
            )
            .unwrap();
            let length = plan.curves().source().source().spec().length_mm;
            for pose in 0..3 {
                let (m, t) = if pose == 0 {
                    ([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3])
                } else {
                    (
                        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
                        [0., 0., length],
                    )
                };
                let mut original =
                    LiteralOccurrencePlacement::new(m, t, Default::default(), None).unwrap();
                if pose == 2 {
                    original = original
                        .transformed(
                            [[0.6, -0.8, 0.], [0.8, 0.6, 0.], [0., 0., 1.]],
                            [2., -3., 0.],
                            Default::default(),
                            None,
                        )
                        .unwrap();
                }
                let candidate =
                    propose_mapped(&local, &original, Default::default(), None).unwrap();
                // Separate ACTUAL-byte read earns material and actual mapping;
                // candidate production itself supplies no material authority.
                let decoded =
                    MappedActualHead::read(Arc::from(candidate.as_str()), Default::default(), None)
                        .unwrap();
                let qualified = decoded
                    .qualify_original_source(
                        plan,
                        HexHeadIntent {
                            across_flats_mm: 18.,
                            height_mm: 3.,
                        },
                        &original,
                        Default::default(),
                        None,
                    )
                    .unwrap();
                assert!(qualified.receipt().world_sewn_reference_position_error_mm <= 1e-7);
                let actual =
                    CompoundDocument::parse(candidate.as_str(), Default::default(), None).unwrap();
                for i in 0..local.original_record_count() {
                    let id = local.original_entity_at(i).unwrap();
                    let mut p = 0;
                    while let Some(component) = local.original_component_at(i, p) {
                        if !component.starts_with("ADVANCED_BREP_SHAPE_REPRESENTATION") {
                            assert_eq!(Some(component), actual.records()[&id][p].source);
                        }
                        p += 1;
                    }
                }
                let encoding = candidate.receipt();
                // Diagnostic accounting only: duplicate local qualification makes
                // this simple composed path exceed272M. Do not reset/relax it or
                // call this a single admitted272M request.
                let composed = encoding.aggregate_work + qualified.receipt().aggregate_work
                    - original.work_used();
                assert!(composed > 272_000_000);
                eprintln!(
                    "MAPPED_ENCODING quarter={quarter} pose={pose} {encoding:?} independent_output={:?} duplicate_composite={composed}",
                    qualified.receipt()
                );
            }
        });
    }
}
#[test]
fn mapped_candidate_encoding_frontiers_cancel_and_retry_are_atomic() {
    use std::cell::Cell;
    let local = OwnedActualHead::read(
        Arc::from(super::super::tests::candidate(false)),
        Default::default(),
        None,
    )
    .unwrap();
    let p = LiteralOccurrencePlacement::new(
        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
        [0., 0., 6.70306],
        Default::default(),
        None,
    )
    .unwrap();
    let good = propose_mapped(&local, &p, Default::default(), None).unwrap();
    let r = good.receipt();
    for (work, bytes, ok) in [
        (r.aggregate_work - 1, 64 * 1024 * 1024, false),
        (r.aggregate_work, 64 * 1024 * 1024, true),
        (272_000_000, r.accounted_peak_bytes - 1, false),
        (272_000_000, r.accounted_peak_bytes, true),
    ] {
        let value = propose_mapped(
            &local,
            &p,
            EncodingOptions {
                maximum_work: work,
                maximum_bytes: bytes,
                ..Default::default()
            },
            None,
        );
        assert_eq!(value.is_ok(), ok);
        if !ok {
            assert!(matches!(value, Err(Error::Limit)));
        }
    }
    let n = Cell::new(0);
    let observe = || {
        n.set(n.get() + 1);
        false
    };
    propose_mapped(&local, &p, Default::default(), Some(&observe)).unwrap();
    for at in [1, n.get() / 2, n.get()] {
        let seen = Cell::new(0);
        let cancel = || {
            seen.set(seen.get() + 1);
            seen.get() >= at
        };
        assert!(matches!(
            propose_mapped(&local, &p, Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    let retry = propose_mapped(&local, &p, Default::default(), None).unwrap();
    assert_eq!(retry.as_str(), good.as_str());
    assert_eq!(retry.receipt().aggregate_work, r.aggregate_work);
}
