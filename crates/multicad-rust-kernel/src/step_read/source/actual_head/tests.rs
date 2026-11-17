use super::*;
use crate::step_read::npt_stock::{
    Error, LiteralOccurrencePlacement,
    full_head::{OwnedActualHead, source_tests::candidate, transport::encoding},
};
use std::cell::Cell;
fn mapped() -> Arc<str> {
    let local = Arc::new(
        OwnedActualHead::read(Arc::from(candidate(false, 18.)), Default::default(), None).unwrap(),
    );
    let frame = LiteralOccurrencePlacement::new(
        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
        [0., 0., local.head_stations_mm()[0]],
        Default::default(),
        None,
    )
    .unwrap()
    .transformed(
        [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
        [20., -30., 11.],
        Default::default(),
        None,
    )
    .unwrap();
    encoding::propose_mapped(&local, &frame, Default::default(), None)
        .unwrap()
        .actual_bytes()
}
fn read(
    input: Arc<str>,
    options: OwnedHeadOptions,
    stop: Option<&dyn Fn() -> bool>,
) -> std::result::Result<Arc<MappedActualHead>, Error> {
    match read_actual_mapped_head_source(input, options, stop)? {
        SourceImport::ActualMappedHead(actual) => Ok(actual),
        _ => panic!("actual reader returned a recipe"),
    }
}
#[test]
fn actual_public_descriptor_whole_file_resource_frontiers_cancel_and_retry() {
    let input = mapped();
    let calls = Cell::new(0usize);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    let actual = read(input.clone(), Default::default(), Some(&observe)).unwrap();
    assert_eq!(actual.actual_file(), input.as_ref());
    assert!(!actual.receipt().exact_local_material_reused);
    let work = actual.receipt().aggregate_work + DESCRIPTOR_WORK;
    let bytes = actual.receipt().accounted_peak_bytes;
    let external = actual.receipt().additional_resident_bytes;
    assert_eq!(external, DESCRIPTOR_BYTES);
    let options = OwnedHeadOptions {
        resources: crate::step_read::npt_stock::full_head::Options {
            maximum_work: work,
            maximum_bytes: bytes,
            ..Default::default()
        },
        ..Default::default()
    };
    assert_eq!(
        read(input.clone(), options, None)
            .unwrap()
            .world_sewn_volume_mm3(),
        actual.world_sewn_volume_mm3()
    );
    for (w, b) in [
        (work - 1, bytes),
        (work, bytes - 1),
        (DESCRIPTOR_WORK, bytes),
    ] {
        let mut low = options;
        low.resources.maximum_work = w;
        low.resources.maximum_bytes = b;
        assert!(
            matches!(read(input.clone(), low, None), Err(Error::Limit)),
            "w={w} b={b}"
        );
    }
    for (w, b) in [
        (0, bytes),
        (work, 0),
        (272_000_001, bytes),
        (work, 64 * 1024 * 1024 + 1),
    ] {
        let mut bad = options;
        bad.resources.maximum_work = w;
        bad.resources.maximum_bytes = b;
        assert!(matches!(
            read(input.clone(), bad, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        read(input.clone(), options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let total = calls.get();
    assert!(total > 3);
    let late = Cell::new(0usize);
    assert!(matches!(
        read(
            input.clone(),
            Default::default(),
            Some(&|| {
                let n = late.get() + 1;
                late.set(n);
                n >= total
            })
        ),
        Err(Error::Cancelled)
    ));
    let retry = read(input.clone(), options, None).unwrap();
    assert_eq!(retry.actual_file(), actual.actual_file());
    assert_eq!(
        retry.world_sewn_enclosure_mm(),
        actual.world_sewn_enclosure_mm()
    );
    assert_eq!(retry.receipt().aggregate_work + DESCRIPTOR_WORK, work);
    eprintln!("ACTUAL_PUBLIC_DESCRIPTOR work={work} peak={bytes} checkpoints={total}");
}
fn change(input: &str, id: u32, body: &str) -> String {
    let a = input.find(&format!("#{id}=")).unwrap();
    let b = a + input[a..].find(';').unwrap() + 1;
    format!("{}#{id}={body};{}", &input[..a], &input[b..])
}
#[test]
fn normal_public_route_refuses_malformed_operator_and_mixed_whole_files_without_recipe_fallback() {
    let input = mapped();
    let doc = CompoundDocument::parse(&input, Default::default(), None).unwrap();
    let member = doc.represented_members()[0];
    let op = member.mapping_link().unwrap().operator_entity();
    let old = doc.records()[&op]
        .iter()
        .find(|p| p.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
        .unwrap()
        .source
        .unwrap();
    // The actual operator remains recognized by topology, but a negative scale
    // cannot become a recipe or a guessed proper intrinsic placement.
    let mut parts = old.to_owned();
    let comma = parts.rfind(',').unwrap();
    parts.replace_range(comma + 1..parts.len() - 1, "-1.");
    let bad = input.replacen(old, &parts, 1);
    assert!(read_source_compound(&bad, 1e-7).is_err());
    assert!(read_actual_mapped_head_source(Arc::from(bad), Default::default(), None).is_err());
    let item = member.item_entity();
    let fresh = doc.records().keys().next_back().unwrap() + 1;
    let item_parts = &doc.records()[&item];
    assert_eq!(item_parts.len(), 1);
    let duplicate = format!("#{fresh}={};\n", item_parts[0].source.unwrap());
    let (&top, top_parts) = doc
        .records()
        .iter()
        .find(|(_, p)| matches!(p.as_slice(),[e] if e.name=="ADVANCED_BREP_SHAPE_REPRESENTATION"))
        .unwrap();
    let [
        Value::String(_),
        Value::List(items),
        Value::Reference(context),
    ] = top_parts[0].args.as_slice()
    else {
        panic!()
    };
    let axis = items
        .iter()
        .find_map(|v| match v {
            Value::Reference(n) if *n != item => Some(*n),
            _ => None,
        })
        .unwrap();
    let mut mixed = change(
        &input,
        top,
        &format!("ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{item},#{fresh},#{axis}),#{context})"),
    );
    let end = mixed.rfind("ENDSEC;").unwrap();
    mixed.insert_str(end, &duplicate);
    let neutral = CompoundDocument::parse(&mixed, Default::default(), None).unwrap();
    assert_eq!(
        neutral.represented_members().len(),
        2,
        "genuine mixed represented file"
    );
    assert!(recognized(&neutral, 16_000_000).unwrap().0);
    drop(neutral);
    assert!(matches!(
        read_source_compound(&mixed, 1e-7),
        Err(ReadError::Unsupported)
    ));
    assert!(read_actual_mapped_head_source(Arc::from(mixed), Default::default(), None).is_err());
    drop(doc);
    assert!(matches!(
        read_source_compound(&input, 1e-7).unwrap().as_slice(),
        [SourceImport::ActualMappedHead(_)]
    ));
}

mod discovery;
