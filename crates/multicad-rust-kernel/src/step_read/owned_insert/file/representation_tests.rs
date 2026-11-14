use super::super::tests::{argument_change, record_change};
use super::*;

fn refused(text: &str, expected: Error) {
    assert!(matches!(File::parse(text, Options::default(), None), Err(e) if e == expected));
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
