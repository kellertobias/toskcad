use super::super::tests::{argument_change, record_change};
use super::*;

#[track_caller]
fn refused(text: &str, expected: Error) {
    let actual = File::parse(text, Options::default(), None).err();
    assert_eq!(actual, Some(expected));
}
#[test]
fn actual_represented_headers_require_one_linked_nonempty_member_list() {
    let (text, _, origin, rep) = tests::mapped_file();
    let file = File::parse(&text, Options::default(), None).unwrap();
    let context = reference(
        &kind(&file.entities, rep, "ADVANCED_BREP_SHAPE_REPRESENTATION")
            .unwrap()
            .args[2],
    )
    .unwrap();
    let link = *file
        .entities
        .iter()
        .find(|(_, r)| r[0].name == "SHAPE_DEFINITION_REPRESENTATION")
        .unwrap()
        .0;
    refused(
        &record_change(&text, rep, |l| {
            l.replace("ADVANCED_BREP_SHAPE_REPRESENTATION", "SHAPE_REPRESENTATION")
        }),
        Error::Unsupported,
    );
    let duplicate = format!(
        "#1100000=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#1000007),#{context});\nENDSEC;\nEND-ISO-10303-21;"
    );
    refused(
        &text.replacen("ENDSEC;\nEND-ISO-10303-21;", &duplicate, 1),
        Error::Unsupported,
    );
    refused(&argument_change(&text, rep, 0, "$"), Error::Invalid);
    refused(&argument_change(&text, link, 1, "#1000005"), Error::Invalid);
    refused(
        &argument_change(&text, rep, 1, "(#1000007,#1000007)"),
        Error::Invalid,
    );
    refused(
        &argument_change(&text, rep, 1, "(#1000000)"),
        Error::Unsupported,
    );
    refused(
        &argument_change(&text, rep, 1, &format!("(#{origin})")),
        Error::Unsupported,
    );
    let mut records = String::new();
    let mut ids = Vec::new();
    for i in 0..=crate::step::MAX_COMPOUND_SOLIDS {
        let id = 1100000 + i;
        records.push_str(&format!("#{id}=MAPPED_ITEM('',#1000006,#1000004);\n"));
        ids.push(format!("#{id}"));
    }
    records.push_str("ENDSEC;\nEND-ISO-10303-21;");
    let many = text.replacen("ENDSEC;\nEND-ISO-10303-21;", &records, 1);
    refused(
        &argument_change(&many, rep, 1, &format!("({})", ids.join(","))),
        Error::Limit,
    );
    assert_eq!(
        File::parse(&text, Options::default(), None)
            .unwrap()
            .represented_members()
            .len(),
        1
    );
}
#[test]
fn actual_map_records_require_complete_local_identity_ownership() {
    let (text, solid, origin, _) = tests::mapped_file();
    for items in [
        format!("(#{solid})"),
        format!("(#{solid},#{solid})"),
        "(#1000000,#1000001)".to_owned(),
    ] {
        refused(
            &argument_change(&text, 1000005, 1, &items),
            Error::Unsupported,
        );
    }
    for (id, arg, value, expected) in [
        (1000007, 0, "$", Error::Invalid),
        (1000005, 0, "$", Error::Invalid),
        (1000004, 0, "$", Error::Unsupported),
    ] {
        let changed = if id == 1000004 {
            record_change(&text, id, |l| {
                l.replace("CARTESIAN_TRANSFORMATION_OPERATOR_3D", "OTHER_OPERATOR_3D")
            })
        } else {
            argument_change(&text, id, arg, value)
        };
        refused(&changed, expected);
    }
    refused(
        &record_change(&text, 1000006, |_| {
            format!("#1000006=REPRESENTATION_MAP(#{origin});")
        }),
        Error::Invalid,
    );
    refused(
        &argument_change(&text, 1000007, 2, "#2000000"),
        Error::Invalid,
    );
    refused(&argument_change(&text, origin, 0, "$"), Error::Invalid);
    let file = File::parse(&text, Options::default(), None).unwrap();
    let axis = kind(&file.entities, origin, "AXIS2_PLACEMENT_3D").unwrap();
    let direction = reference(&axis.args[2]).unwrap();
    refused(&argument_change(&text, direction, 0, "$"), Error::Invalid);
    refused(
        &argument_change(&text, direction, 1, "(0.,1.)"),
        Error::Unsupported,
    );
    refused(
        &argument_change(&text, direction, 1, "(0.,0.,-1.)"),
        Error::Unsupported,
    );
    refused(
        &argument_change(&text, origin, 2, ".T."),
        Error::Unsupported,
    );
    // Default source axes have the same exact datum semantics as positive axis ratios.
    let omitted = argument_change(&argument_change(&text, origin, 2, "$"), origin, 3, "$");
    assert_eq!(
        File::parse(&omitted, Options::default(), None)
            .unwrap()
            .represented_members()
            .len(),
        1
    );
}
// Real authored geometry, with only product/representation metadata changed.
// Wrapper datums do not place the BREP; no transformed/extra body is discarded.
fn related_box(wrappers: usize, reverse: bool) -> (String, String, u32, u32, u32) {
    let q = crate::Rational::new(4, 1).unwrap();
    let original = crate::step::write_box([q; 3], [0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    let file = File::parse(&original, Options::default(), None).unwrap();
    let (&link, parts) = file
        .entities
        .iter()
        .find(|(_, p)| p[0].name == "SHAPE_DEFINITION_REPRESENTATION")
        .unwrap();
    let [Value::Reference(shape), Value::Reference(top)] = parts[0].args.as_slice() else {
        panic!("actual SDR")
    };
    let (shape, top) = (*shape, *top);
    let context = reference(
        &kind(&file.entities, top, "ADVANCED_BREP_SHAPE_REPRESENTATION")
            .unwrap()
            .args[2],
    )
    .unwrap();
    let datum = *file
        .entities
        .iter()
        .find(|(_, p)| p[0].name == "AXIS2_PLACEMENT_3D")
        .unwrap()
        .0;
    let mut records = String::new();
    for i in 0..wrappers {
        let id = 1_100_000 + i as u32;
        let next = if i + 1 == wrappers { top } else { id + 1 };
        let (a, b) = if reverse { (next, id) } else { (id, next) };
        records.push_str(&format!("#{id}=SHAPE_REPRESENTATION('',(#{datum}),#{context});\n#{}=SHAPE_REPRESENTATION_RELATIONSHIP('',$,#{a},#{b});\n",1_200_000+i));
    }
    let linked = argument_change(&original, link, 1, "#1100000");
    let end = linked.rfind("ENDSEC;").unwrap();
    let linked = format!("{}{}{}", &linked[..end], records, &linked[end..]);
    (original, linked, shape, top, context)
}
#[test]
fn legal_related_representation_paths_keep_actual_members_and_selected_geometry() {
    for reverse in [false, true] {
        for count in [1, 3, 63] {
            let (original, text, _, _, _) = related_box(count, reverse);
            let baseline = File::parse(&original, Options::default(), None).unwrap();
            let file = File::parse(&text, Options::default(), None).unwrap();
            assert_eq!(file.represented_members(), baseline.represented_members());
            let member = file.represented_members()[0];
            let view = file
                .local_member_document(member.item_entity(), Options::default(), None)
                .unwrap();
            assert!(std::ptr::eq(view.source(), &file));
            assert_eq!(view.member(), member);
            let actual = crate::step_read::read_analytic_solid(view.as_str(), 1e-7).unwrap();
            assert_eq!(actual.solid_count, 1);
            assert_eq!(actual.face_count, 6);
            assert!((actual.volume_mm3 - 64.).abs() < 1e-10);
            assert_eq!(file.literal_input(), text);
            assert_eq!(
                crate::step_read::read_analytic_solid(&text, 1e-7)
                    .unwrap()
                    .face_count,
                6
            );
        }
    }
}
#[test]
fn related_graphs_refuse_ambiguous_cycles_extra_geometry_and_bad_schema_atomically() {
    let (original, text, shape, top, context) = related_box(2, false);
    let add = |extra: &str| {
        let end = text.rfind("ENDSEC;").unwrap();
        format!("{}{}{}", &text[..end], extra, &text[end..])
    };
    // A second path gives no unique metadata tree; no branch is silently picked.
    refused(
        &add(&format!(
            "#1300000=SHAPE_REPRESENTATION_RELATIONSHIP('','',#1100000,#{top});\n"
        )),
        Error::Unsupported,
    );
    // Two product shapes claim the same related represented graph.
    let base = File::parse(&original, Options::default(), None).unwrap();
    let product = kind(&base.entities, shape, "PRODUCT_DEFINITION_SHAPE")
        .unwrap()
        .source
        .unwrap();
    refused(
        &add(&format!(
            "#1300000={product};\n#1300001=SHAPE_DEFINITION_REPRESENTATION(#1300000,#1100000);\n"
        )),
        Error::Unsupported,
    );
    refused(
        &argument_change(&text, 1_200_001, 2, &format!("#{top}")),
        Error::Unsupported,
    );
    refused(&argument_change(&text, 1_200_000, 0, "$"), Error::Invalid);
    refused(&argument_change(&text, 1_200_000, 1, ".T."), Error::Invalid);
    refused(
        &argument_change(&text, 1_200_000, 2, "#2000000000"),
        Error::Invalid,
    );
    refused(&argument_change(&text, 1_100_000, 0, "$"), Error::Invalid);
    refused(
        &argument_change(&text, 1_100_000, 1, "(#1100001)"),
        Error::Unsupported,
    );
    let datum = base.represented_members()[0].local_solid_entity();
    // Duplicate datum SET items are malformed, even when individually valid.
    let actual = File::parse(&text, Options::default(), None).unwrap();
    let axis = refs(
        &kind(&actual.entities, 1_100_000, "SHAPE_REPRESENTATION")
            .unwrap()
            .args[1],
    )
    .unwrap()[0];
    refused(
        &argument_change(&text, 1_100_000, 1, &format!("(#{axis},#{axis})")),
        Error::Invalid,
    );
    refused(
        &argument_change(&text, 1_100_000, 1, &format!("(#{datum})")),
        Error::Unsupported,
    );
    let bad_dimension = argument_change(&text, 1_100_000, 2, "#1300000");
    let end = bad_dimension.rfind("ENDSEC;").unwrap();
    let bad_dimension = format!(
        "{}#1300000=GEOMETRIC_REPRESENTATION_CONTEXT(2.);{}",
        &bad_dimension[..end],
        &bad_dimension[end..]
    );
    assert!(File::parse(&bad_dimension, Options::default(), None).is_err());
    let transformed = record_change(&text, 1_200_000, |_| {
        "#1200000=(REPRESENTATION_RELATIONSHIP('','',#1100000,#1100001)SHAPE_REPRESENTATION_RELATIONSHIP());".to_owned()
    });
    refused(&transformed, Error::Unsupported);
    let (_, too_many, _, _, _) = related_box(64, false);
    refused(&too_many, Error::Limit);
    assert_eq!(
        File::parse(&text, Options::default(), None)
            .unwrap()
            .represented_members()
            .len(),
        1
    );
    let _ = context;
}
#[test]
fn related_metadata_does_not_allocate_unbounded_datums_or_claim_unrelated_transform_graphs() {
    let (_, text, _, _, context) = related_box(2, false);
    let actual = File::parse(&text, Options::default(), None).unwrap();
    let axis = refs(
        &kind(&actual.entities, 1_100_000, "SHAPE_REPRESENTATION")
            .unwrap()
            .args[1],
    )
    .unwrap()[0];
    let datum = kind(&actual.entities, axis, "AXIS2_PLACEMENT_3D")
        .unwrap()
        .source
        .unwrap();
    let with_datums = |count: usize| {
        let records: String = (0..count)
            .map(|i| format!("#{}={datum};\n", 1_300_000 + i))
            .collect();
        let end = text.rfind("ENDSEC;").unwrap();
        let extended = format!("{}{}{}", &text[..end], records, &text[end..]);
        let ids: Vec<String> = (0..count).map(|i| format!("#{}", 1_300_000 + i)).collect();
        argument_change(&extended, 1_100_000, 1, &format!("({})", ids.join(",")))
    };
    assert!(File::parse(&with_datums(64), Options::default(), None).is_ok());
    refused(&with_datums(65), Error::Limit);
    let end = text.rfind("ENDSEC;").unwrap();
    let unrelated = format!(
        "{}#1400000=SHAPE_REPRESENTATION('',(#{axis}),#{context});\n#1400001=SHAPE_REPRESENTATION('',(#{axis}),#{context});\n#1400002=(REPRESENTATION_RELATIONSHIP('','',#1400000,#1400001)SHAPE_REPRESENTATION_RELATIONSHIP());\n{}",
        &text[..end],
        &text[end..]
    );
    let document = File::parse(&unrelated, Options::default(), None).unwrap();
    assert_eq!(document.represented_members(), actual.represented_members());
    assert_eq!(
        crate::step_read::read_analytic_solid(&unrelated, 1e-7)
            .unwrap()
            .face_count,
        6
    );
    let connected = unrelated.replace("#1400000,#1400001)", "#1100000,#1400001)");
    assert_ne!(connected, unrelated);
    refused(&connected, Error::Unsupported);
    assert!(File::parse(&unrelated, Options::default(), None).is_ok());
}
#[test]
fn related_representation_queries_keep_exact_work_frontier_and_observed_cancellation_retry() {
    use std::cell::Cell;
    let (_, text, _, _, _) = related_box(3, true);
    let calls = Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let file = File::parse(&text, Options::default(), Some(&count)).unwrap();
    let work = file.work_used();
    let polls = calls.get();
    assert!(polls > 10);
    refused_with_work(&text, work - 1, Error::Limit);
    let exact = Options {
        maximum_work: work,
        ..Options::default()
    };
    assert!(File::parse(&text, exact, None).is_ok());
    for stop in [1, polls / 2, polls] {
        let seen = Cell::new(0usize);
        let cancel = || {
            seen.set(seen.get() + 1);
            seen.get() >= stop
        };
        assert!(matches!(
            File::parse(&text, Options::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
        assert!(seen.get() >= stop);
        assert_eq!(
            File::parse(&text, exact, None)
                .unwrap()
                .represented_members(),
            file.represented_members()
        );
    }
}
fn refused_with_work(text: &str, work: usize, expected: Error) {
    assert!(
        matches!(File::parse(text,Options{maximum_work:work,..Options::default()},None),Err(e) if e==expected)
    );
}
