use super::*;
use crate::{Rational, step::Writer};
use std::cell::Cell;

fn boxes(count: usize) -> String {
    let mut writer = Writer::new();
    for _ in 0..count {
        writer
            .append_box(
                [Rational::new(4, 1).unwrap(); 3],
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
    }
    writer.finish().unwrap()
}
fn document(text: &str) -> CompoundDocument<'_> {
    CompoundDocument::parse(text, Default::default(), None).unwrap()
}
fn compare(
    a: &CompoundDocument<'_>,
    b: &CompoundDocument<'_>,
    index: usize,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Binding> {
    bind_local_members(
        a,
        a.represented_members()[0].item_entity(),
        b,
        b.represented_members()[index].item_entity(),
        options,
        interrupt,
    )
}
#[test]
fn actual_written_members_bind_with_rebased_ids_and_separate_document_arenas() {
    let single = boxes(1);
    let mixed = boxes(2);
    let a = document(&single);
    let b = document(&mixed);
    let first = compare(&a, &b, 0, Default::default(), None).unwrap();
    let second = compare(&a, &b, 1, Default::default(), None).unwrap();
    assert_eq!(first.matched_entities(), second.matched_entities());
    assert!(second.matched_entities() > 50);
    assert!(second.work_used() > second.matched_entities());
    assert!(second.query_peak_bytes() >= second.matched_entities() * PAIR_BYTES);
    assert_ne!(
        a.represented_members()[0].item_entity(),
        b.represented_members()[1].item_entity()
    );
    assert_eq!(a.literal_input(), single);
    assert_eq!(b.literal_input(), mixed);
}
#[test]
fn original_decimal_leaves_and_owner_sharing_cannot_be_replaced_by_float_equality() {
    let text = boxes(1);
    let a = document(&text);
    let (point, source) = a
        .records()
        .iter()
        .flat_map(|(id, records)| records.iter().map(move |r| (*id, r)))
        .find(|(_, record)| {
            record.name == "CARTESIAN_POINT"
                && part21::numeric_spellings(record)
                    .unwrap()
                    .iter()
                    .any(|n| n.parse::<f64>().unwrap() == 4.)
        })
        .unwrap();
    let source = source.source.unwrap();
    let number = *part21::numeric_spellings(&a.records()[&point][0])
        .unwrap()
        .iter()
        .find(|n| n.parse::<f64>().unwrap() == 4.)
        .unwrap();
    let alias = "4.0000000000000000000000000000000001";
    assert_eq!(
        number.parse::<f64>().unwrap(),
        alias.parse::<f64>().unwrap()
    );
    let changed = text.replacen(source, &source.replacen(number, alias, 1), 1);
    assert!(matches!(
        compare(&a, &document(&changed), 0, Default::default(), None),
        Err(Error::Different)
    ));

    // Preserve every coordinate but split one actual shared point owner.
    let duplicated = format!("#1100000={source};\nENDSEC;\nEND-ISO-10303-21;");
    let with_duplicate = text.replacen("ENDSEC;\nEND-ISO-10303-21;", &duplicated, 1);
    let reference = format!("#{point}");
    let owner_line = text
        .lines()
        .find(|line| line.contains("=VERTEX_POINT(") && line.ends_with(&format!("{reference});")))
        .unwrap();
    let changed =
        with_duplicate.replacen(owner_line, &owner_line.replace(&reference, "#1100000"), 1);
    assert!(matches!(
        compare(&a, &document(&changed), 0, Default::default(), None),
        Err(Error::Different)
    ));
}
#[test]
fn actual_query_resource_and_cancellation_frontiers_preserve_original_documents() {
    let text = boxes(1);
    let a = document(&text);
    let good = compare(&a, &a, 0, Default::default(), None).unwrap();
    for options in [
        Options {
            maximum_work: good.work_used() - 1,
            ..Default::default()
        },
        Options {
            maximum_bytes: good.query_peak_bytes() - 1,
            ..Default::default()
        },
        Options {
            maximum_entities: good.matched_entities() - 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            compare(&a, &a, 0, options, None),
            Err(Error::Limit)
        ));
    }
    assert!(
        compare(
            &a,
            &a,
            0,
            Options {
                maximum_work: good.work_used(),
                maximum_bytes: good.query_peak_bytes(),
                maximum_entities: good.matched_entities()
            },
            None
        )
        .is_ok()
    );
    let calls = Cell::new(0usize);
    compare(
        &a,
        &a,
        0,
        Default::default(),
        Some(&|| {
            calls.set(calls.get() + 1);
            false
        }),
    )
    .unwrap();
    for stop in [1, calls.get() / 2, calls.get()] {
        let count = Cell::new(0usize);
        assert!(matches!(
            compare(
                &a,
                &a,
                0,
                Default::default(),
                Some(&|| {
                    count.set(count.get() + 1);
                    count.get() == stop
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    for options in [
        Options {
            maximum_work: 0,
            ..Default::default()
        },
        Options {
            maximum_work: 64_000_001,
            ..Default::default()
        },
        Options {
            maximum_bytes: 0,
            ..Default::default()
        },
        Options {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
        Options {
            maximum_entities: 0,
            ..Default::default()
        },
        Options {
            maximum_entities: 1_048_577,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            compare(&a, &a, 0, options, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        bind_local_members(
            &a,
            u32::MAX,
            &a,
            a.represented_members()[0].item_entity(),
            Default::default(),
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        bind_local_members(
            &a,
            a.represented_members()[0].item_entity(),
            &a,
            u32::MAX,
            Default::default(),
            None
        ),
        Err(Error::Invalid)
    ));
    assert_eq!(a.literal_input(), text);
}
