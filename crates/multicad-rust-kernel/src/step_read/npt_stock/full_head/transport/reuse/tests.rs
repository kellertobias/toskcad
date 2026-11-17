use super::*;
use crate::step_read::npt_stock::{
    LiteralOccurrencePlacement, full_head::source_tests::with_source,
};
fn local_owner(quarter: bool) -> Arc<OwnedActualHead> {
    Arc::new(
        OwnedActualHead::read(
            Arc::from(super::super::tests::candidate(quarter)),
            Default::default(),
            None,
        )
        .unwrap(),
    )
}
fn pose(length: f64) -> LiteralOccurrencePlacement {
    LiteralOccurrencePlacement::new(
        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
        [0., 0., length],
        Default::default(),
        None,
    )
    .unwrap()
    .transformed(
        [[0.6, -0.8, 0.], [0.8, 0.6, 0.], [0., 0., 1.]],
        [2., -3., 0.],
        Default::default(),
        None,
    )
    .unwrap()
}
#[test]
fn actual_output_full_component_binding_transfers_material_once_and_complete_pipeline_fits_unchanged_cap()
 {
    for quarter in [false, true] {
        with_source(quarter, |source| {
            let local = local_owner(quarter);
            let p = pose(source.curves().source().source().spec().length_mm);
            let candidate = encoding::propose_mapped(
                &local,
                &p,
                encoding::EncodingOptions {
                    maximum_work: 272_000_000 - source.work_used(),
                    additional_resident_bytes: source.retained_budget_bytes(),
                    ..Default::default()
                },
                None,
            )
            .unwrap();
            let encoding = candidate.receipt();
            let input = candidate.actual_bytes();
            drop(candidate);
            let decoded = MappedActualHead::read_reusing_local(
                Arc::clone(&input),
                Arc::clone(&local),
                OwnedHeadOptions {
                    resources: Options {
                        maximum_work: 272_000_000
                            - source.work_used()
                            - p.work_used()
                            - encoding.fresh_encoding_work,
                        additional_resident_bytes: source.retained_budget_bytes()
                            + p.retained_budget_bytes(),
                        ..Default::default()
                    },
                    ..Default::default()
                },
                None,
            )
            .unwrap();
            assert!(decoded.receipt().exact_local_material_reused);
            assert!(std::ptr::eq(decoded.local_geometry(), &*local));
            let qualified = decoded
                .qualify_original_source(
                    source,
                    HexHeadIntent {
                        across_flats_mm: 18.,
                        height_mm: 3.,
                    },
                    &p,
                    SourceOptions {
                        resources: Options {
                            maximum_work: 272_000_000 - encoding.fresh_encoding_work,
                            ..Default::default()
                        },
                        ..Default::default()
                    },
                    None,
                )
                .unwrap();
            let aggregate = qualified.receipt().aggregate_work + encoding.fresh_encoding_work;
            assert!(aggregate <= 272_000_000);
            assert!(
                qualified
                    .receipt()
                    .accounted_peak_bytes
                    .max(encoding.accounted_peak_bytes)
                    <= 64 * 1024 * 1024
            );
            // Cold public import independently re-earns material, using no owner,
            // trace, source recipe or pointer association from the export path.
            let cold = MappedActualHead::read(input, Default::default(), None).unwrap();
            assert!(!cold.receipt().exact_local_material_reused);
            assert_eq!(
                decoded.world_sewn_enclosure_mm(),
                cold.world_sewn_enclosure_mm()
            );
            assert_eq!(
                decoded.world_sewn_volume_mm3(),
                cold.world_sewn_volume_mm3()
            );
            assert_eq!(
                decoded.local_geometry().head_polygon_mm(),
                cold.local_geometry().head_polygon_mm()
            );
            let cold_source = cold
                .qualify_original_source(
                    source,
                    HexHeadIntent {
                        across_flats_mm: 18.,
                        height_mm: 3.,
                    },
                    &p,
                    Default::default(),
                    None,
                )
                .unwrap();
            assert_eq!(
                qualified.receipt().original_raw_thread.position_error_mm(),
                cold_source
                    .receipt()
                    .original_raw_thread
                    .position_error_mm()
            );
            assert_eq!(
                qualified.receipt().world_sewn_reference_position_error_mm,
                cold_source.receipt().world_sewn_reference_position_error_mm
            );
            eprintln!(
                "MAPPED_REUSE quarter={quarter} encoding={encoding:?} source={:?} aggregate={aggregate}",
                qualified.receipt()
            );
        });
    }
}
fn change(input: &str, id: u32, body: &str) -> String {
    let start = input.find(&format!("#{id}=")).unwrap();
    let end = start + input[start..].find(';').unwrap() + 1;
    format!("{}#{id}={body};{}", &input[..start], &input[end..])
}
#[test]
fn actual_output_binding_rejects_changed_scalar_topology_metadata_and_orphans_and_is_not_pointer_trust()
 {
    let local = local_owner(false);
    let p = pose(6.70306);
    let candidate = encoding::propose_mapped(&local, &p, Default::default(), None).unwrap();
    let input = candidate.actual_bytes();
    drop(candidate);
    let doc = CompoundDocument::parse(&input, Default::default(), None).unwrap();
    let point = local.physical_roles().head_vertex_entities[0];
    let vertex = &doc.records()[&point][0];
    let Value::Reference(coord) = vertex.args[1] else {
        panic!()
    };
    // Identical numeric value with DIFFERENT actual scalar spelling still
    // refuses literal identity. It can cold-import independently as geometry.
    let original = doc.records()[&coord][0].source.unwrap();
    let scalar = crate::step_read::part21::numeric_spellings(&doc.records()[&coord][0]).unwrap()[0];
    let spelling = if let Some(unsigned) = scalar.strip_prefix('-') {
        format!("-0{unsigned}")
    } else {
        format!("+{scalar}")
    };
    let changed = original.replacen(scalar, &spelling, 1);
    assert_ne!(changed, original);
    let wrong = change(&input, coord, &changed);
    assert!(MappedActualHead::read(Arc::from(wrong.as_str()), Default::default(), None).is_ok());
    assert!(matches!(
        MappedActualHead::read_reusing_local(
            Arc::from(wrong),
            Arc::clone(&local),
            Default::default(),
            None
        ),
        Err(Error::Unsupported)
    ));
    // Identical geometry but DIFFERENT point-reference ownership. Reuse an
    // existing equal-valued LINE origin and remove the now-unreferenced old
    // vertex point record, so no orphan is smuggled into this positive file.
    let shared = doc
        .records()
        .iter()
        .find_map(|(&id, p)| {
            (id != coord && matches!(p.as_slice(),[e] if e.source==Some(original))).then_some(id)
        })
        .unwrap();
    let mut sharing = change(&input, point, &format!("VERTEX_POINT('',#{shared})"));
    let start = sharing.find(&format!("#{coord}=")).unwrap();
    let end = start + sharing[start..].find(';').unwrap() + 1;
    sharing.replace_range(start..end, "");
    MappedActualHead::read(Arc::from(sharing.as_str()), Default::default(), None).unwrap();
    assert!(matches!(
        MappedActualHead::read_reusing_local(
            Arc::from(sharing),
            Arc::clone(&local),
            Default::default(),
            None
        ),
        Err(Error::Unsupported)
    ));
    let (&product, parts) = doc
        .records()
        .iter()
        .find(|(_, p)| matches!(p.as_slice(),[e] if e.name=="PRODUCT"))
        .unwrap();
    let source = parts[0].source.unwrap();
    let metadata = change(
        &input,
        product,
        &source.replacen("PRODUCT('", "PRODUCT('changed", 1),
    );
    assert!(matches!(
        MappedActualHead::read_reusing_local(
            Arc::from(metadata),
            Arc::clone(&local),
            Default::default(),
            None
        ),
        Err(Error::Unsupported)
    ));
    let shell = doc
        .records()
        .values()
        .find(|p| matches!(p.as_slice(),[e] if e.name=="CLOSED_SHELL"))
        .unwrap()[0]
        .source
        .unwrap();
    let fresh = *doc.records().keys().next_back().unwrap() + 1;
    let mut orphan = input.to_string();
    let end = orphan.rfind("ENDSEC;").unwrap();
    orphan.insert_str(end, &format!("#{fresh}={shell};\n"));
    assert!(matches!(
        MappedActualHead::read_reusing_local(
            Arc::from(orphan),
            Arc::clone(&local),
            Default::default(),
            None
        ),
        Err(Error::Unsupported)
    ));
    // A fresh independently earned owner over the SAME original components
    // succeeds despite distinct arena/data addresses. Complete equality is
    // the authority transfer, not pointer/length or writer provenance.
    let independent = local_owner(false);
    assert!(!Arc::ptr_eq(&local, &independent));
    let good =
        MappedActualHead::read_reusing_local(input, independent, Default::default(), None).unwrap();
    assert_eq!(
        good.local_geometry().head_polygon_mm(),
        local.head_polygon_mm()
    );
}
#[test]
fn source_known_output_work_bytes_cancel_and_failed_retry_preserve_the_prior_owner() {
    use std::cell::Cell;
    let local = local_owner(false);
    let p = pose(6.70306);
    let candidate = encoding::propose_mapped(&local, &p, Default::default(), None).unwrap();
    let input = candidate.actual_bytes();
    drop(candidate);
    let good = MappedActualHead::read_reusing_local(
        Arc::clone(&input),
        Arc::clone(&local),
        Default::default(),
        None,
    )
    .unwrap();
    let r = good.receipt();
    drop(good);
    let references = Arc::strong_count(&local);
    for (work, bytes, ok) in [
        (r.aggregate_work - 1, 64 * 1024 * 1024, false),
        (r.aggregate_work, 64 * 1024 * 1024, true),
        (272_000_000, r.accounted_peak_bytes - 1, false),
        (272_000_000, r.accounted_peak_bytes, true),
    ] {
        let value = MappedActualHead::read_reusing_local(
            Arc::clone(&input),
            Arc::clone(&local),
            OwnedHeadOptions {
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
        assert_eq!(Arc::strong_count(&local), references);
    }
    let n = Cell::new(0);
    let observe = || {
        n.set(n.get() + 1);
        false
    };
    MappedActualHead::read_reusing_local(
        Arc::clone(&input),
        Arc::clone(&local),
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
            MappedActualHead::read_reusing_local(
                Arc::clone(&input),
                Arc::clone(&local),
                Default::default(),
                Some(&cancel)
            ),
            Err(Error::Cancelled)
        ));
        assert_eq!(Arc::strong_count(&local), references);
    }
    let retry =
        MappedActualHead::read_reusing_local(input, Arc::clone(&local), Default::default(), None)
            .unwrap();
    assert_eq!(retry.receipt().aggregate_work, r.aggregate_work);
}
