//! Mutate the actual emitted records; retain the ordinary neutral/read pipeline.
use super::*;

fn replace_field(text: &str, id: u32, index: usize, replacement: &str) -> String {
    change_record(text, id, |record| {
        let mut values = fields(record);
        values[index] = replacement;
        with_fields(record, &values)
    })
}

#[test]
fn actual_tensor_rows_and_curve_subdivision_counts_must_match_the_neutral_schema() {
    let text = fixture(false);
    let tensor = record_id(&text, "=B_SPLINE_SURFACE_WITH_KNOTS(", 0);
    let ragged = change_record(&text, tensor, |record| {
        let mut args = fields(record);
        let mut rows = fields(args[3]);
        let mut first = fields(rows[0]);
        first.pop().unwrap();
        let shortened = format!("({})", first.join(","));
        rows[0] = &shortened;
        let grid = format!("({})", rows.join(","));
        args[3] = &grid;
        with_fields(record, &args)
    });
    assert_eq!(rejected(&ragged), Error::Invalid);
    let curve = record_id(&text, "=B_SPLINE_CURVE_WITH_KNOTS('',11,", 0);
    let non_dyadic = change_record(&text, curve, |record| {
        let mut args = fields(record);
        args[6] = "(12,11,11,12)";
        args[7] = "(0.,0.3333333333333333,0.6666666666666667,1.)";
        with_fields(record, &args)
    });
    assert_eq!(rejected(&non_dyadic), Error::Unsupported);
    let mismatched = replace_field(&text, curve, 6, "(12,12)");
    assert_eq!(rejected(&mismatched), Error::Invalid);
    let document = CompoundDocument::parse(&text, Default::default(), None).unwrap();
    assert!(
        Parsed::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .is_ok()
    );
    assert_eq!(document.literal_input(), text);
}

#[test]
fn actual_packet_semantic_mutations_refuse_false_bindings() {
    let text = fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let solid = document.represented_members()[0].item_entity();
    assert!(matches!(
        Parsed::parse(&document, solid + 1, Default::default(), None),
        Err(Error::Unsupported)
    ));
    let mut cases = vec![];
    for (fragment, index, replacements) in [
        ("=B_SPLINE_CURVE_WITH_KNOTS('',11,", 1, vec!["11.0", "-11"]),
        (
            "=B_SPLINE_SURFACE_WITH_KNOTS(",
            8,
            vec!["(12)", "(12,11,11,12)", "(12,10,12)"],
        ),
        ("=B_SPLINE_SURFACE_WITH_KNOTS(", 10, vec!["(0.,0.1,1.)"]),
        ("=DIRECTION(", 1, vec!["(0.,1.)"]),
        ("=ADVANCED_FACE(", 3, vec![".T."]),
        ("=EDGE_CURVE(", 4, vec![".F."]),
        ("=SURFACE_CURVE(", 3, vec![".CURVE_3D."]),
        ("=B_SPLINE_CURVE_WITH_KNOTS('',1,", 1, vec!["2", "1.0"]),
        (
            "=B_SPLINE_CURVE_WITH_KNOTS('',1,",
            6,
            vec!["(2,1,2)", "(1,2)"],
        ),
        ("=B_SPLINE_CURVE_WITH_KNOTS('',1,", 7, vec!["(0.,0.5)"]),
    ] {
        let id = record_id(&text, fragment, 0);
        for replacement in replacements {
            cases.push(replace_field(&text, id, index, replacement));
        }
    }
    let records = part21::parse(&text).unwrap();
    let shell = record_id(&text, "=CLOSED_SHELL(", 0);
    let shell_record = text
        .lines()
        .find(|line| line.starts_with(&format!("#{shell}=")))
        .unwrap();
    let face_list = fields(shell_record)[1];
    let face_ids = face_list
        .trim_matches(['(', ')'])
        .split(',')
        .collect::<Vec<_>>();
    let duplicate_faces = format!(
        "({},{})",
        face_ids[0],
        face_ids[1..]
            .iter()
            .map(|_| face_ids[0])
            .collect::<Vec<_>>()
            .join(",")
    );
    cases.push(replace_field(&text, shell, 1, &duplicate_faces));
    let first_face = first_reference(face_ids[0]);
    let second_face = first_reference(face_ids[1]);
    let first_support = match &records[&first_face][0].args[2] {
        Value::Reference(id) => *id,
        _ => panic!(),
    };
    cases.push(replace_field(
        &text,
        second_face,
        2,
        &format!("#{first_support}"),
    ));
    let point = record_id(&text, "=CARTESIAN_POINT(", 0);
    cases.push(replace_field(&text, first_face, 2, &format!("#{point}")));
    let edge = record_id(&text, "=EDGE_CURVE(", 0);
    let edge_record = text
        .lines()
        .find(|line| line.starts_with(&format!("#{edge}=")))
        .unwrap();
    let edge_fields = fields(edge_record);
    cases.push(replace_field(&text, edge, 2, edge_fields[1]));
    let surface_curve = record_id(&text, "=SURFACE_CURVE(", 0);
    let surface_record = text
        .lines()
        .find(|line| line.starts_with(&format!("#{surface_curve}=")))
        .unwrap();
    let pcurves = fields(surface_record)[2]
        .trim_matches(['(', ')'])
        .split(',')
        .collect::<Vec<_>>();
    cases.push(replace_field(
        &text,
        surface_curve,
        2,
        &format!("({0},{0})", pcurves[0]),
    ));
    cases.push(replace_field(
        &text,
        surface_curve,
        2,
        &format!("({0},{1},{0})", pcurves[0], pcurves[1]),
    ));
    let line = record_id(&text, "=B_SPLINE_CURVE_WITH_KNOTS('',1,", 0);
    let line_record = text
        .lines()
        .find(|record| record.starts_with(&format!("#{line}=")))
        .unwrap();
    let poles = fields(line_record)[2]
        .trim_matches(['(', ')'])
        .split(',')
        .collect::<Vec<_>>();
    cases.push(replace_field(
        &text,
        line,
        2,
        &format!("({0},{0})", poles[0]),
    ));
    for (index, changed) in cases.iter().enumerate() {
        assert!(
            refuses_pipeline(changed, solid),
            "semantic mutation {index} accepted"
        );
    }
    assert!(Parsed::parse(&document, solid, Default::default(), None).is_ok());
    assert_eq!(document.literal_input(), text);
}

#[test]
fn original_numeric_records_must_be_consumed_and_have_the_requested_entity_kind() {
    let text = fixture(false);
    let records = part21::parse(&text).unwrap();
    let point = record_id(&text, "=CARTESIAN_POINT(", 0);
    let mut budget = Budget::new(Default::default(), None).unwrap();
    assert!(matches!(
        records::kind(&records, point, "DIRECTION", &mut budget),
        Err(Error::Unsupported)
    ));
    let entity = &records[&point][0];
    let mut numbers = numeric::Numbers::new(entity, &mut budget).unwrap();
    let Value::List(coordinates) = &entity.args[1] else {
        panic!()
    };
    numbers.scalar(&coordinates[0], &mut budget).unwrap();
    assert_eq!(numbers.finish(), Err(Error::Invalid));
}

#[test]
fn unrelated_point_records_do_not_change_selected_geometry_or_query_resources() {
    use std::fmt::Write;
    let mut text = fixture(false);
    let original = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let solid = original.represented_members()[0].item_entity();
    let baseline = Parsed::parse(&original, solid, Default::default(), None).unwrap();
    let expected = (
        baseline.points().len(),
        baseline.work_used(),
        baseline.query_budget_bytes(),
    );
    let original_bytes = original.retained_budget_bytes();
    drop(baseline);
    drop(original);
    let mut extra = String::new();
    for id in 1_000_000..1_000_000 + 262_144 + 4097 {
        writeln!(&mut extra, "#{id}=CARTESIAN_POINT('',(0.,0.,0.));").unwrap();
    }
    let end = text.rfind("ENDSEC;").unwrap();
    text.insert_str(end, &extra);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let selected = Parsed::parse(&document, solid, Default::default(), None).unwrap();
    assert_eq!(
        (
            selected.points().len(),
            selected.work_used(),
            selected.query_budget_bytes()
        ),
        expected
    );
    assert!(document.retained_budget_bytes() > original_bytes);
    assert_eq!(document.literal_input(), text);
}

#[test]
fn actual_graph_query_cancellation_prefixes_preserve_the_original_arena() {
    let text = fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let solid = document.represented_members()[0].item_entity();
    let polls = std::cell::Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    let baseline = Parsed::parse(&document, solid, Default::default(), Some(&count)).unwrap();
    let total = polls.get();
    for stop in (1..=total).step_by((total / 64).max(1)) {
        polls.set(0);
        let cancel = || {
            polls.set(polls.get() + 1);
            polls.get() == stop
        };
        assert!(matches!(
            Parsed::parse(&document, solid, Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
        assert_eq!(polls.get(), stop);
    }
    let recovered = Parsed::parse(&document, solid, Default::default(), None).unwrap();
    assert_eq!(recovered.work_used(), baseline.work_used());
    assert_eq!(
        recovered.query_budget_bytes(),
        baseline.query_budget_bytes()
    );
    assert_eq!(document.literal_input(), text);
}

#[test]
fn distinct_boundary_point_owners_require_exact_decimal_equality() {
    let text = fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let solid = document.represented_members()[0].item_entity();
    let baseline = Parsed::parse(&document, solid, Default::default(), None).unwrap();
    let face = baseline
        .faces
        .iter()
        .find(|f| matches!(f.support, Support::Tensor(_)))
        .unwrap();
    let curve = &baseline.curves[baseline.edges[face.coedges[0].edge].curve];
    let original = baseline.points[curve.poles[1]].id;
    let new_id = document.records().last_key_value().unwrap().0 + 1;
    let Value::List(xyz) = &document.records()[&original][0].args[1] else {
        panic!()
    };
    let axis = xyz
        .iter()
        .position(|v| matches!(v, Value::Number(n) if *n!=0.))
        .unwrap();
    let point_line = text
        .lines()
        .find(|line| line.starts_with(&format!("#{original}=")))
        .unwrap();
    let coordinates = fields(point_line)[1]
        .trim_matches(['(', ')'])
        .split(',')
        .collect::<Vec<_>>();
    let original_token = coordinates[axis];
    let (mantissa, exponent) = original_token
        .split_once('E')
        .unwrap_or((original_token, "0"));
    assert!(mantissa.contains('.'));
    for exact in [true, false] {
        let suffix = if exact {
            "0"
        } else {
            "0000000000000000000000000000000000000001"
        };
        let token = format!("{mantissa}{suffix}E{exponent}");
        assert_eq!(
            token.parse::<f64>().unwrap(),
            original_token.parse::<f64>().unwrap()
        );
        let mut xyz = coordinates.clone();
        xyz[axis] = &token;
        let alias = format!("#{new_id}=CARTESIAN_POINT('',({}));\n", xyz.join(","));
        let mut changed = change_record(&text, curve.id, |line| {
            let mut args = fields(line);
            let mut poles = args[2]
                .trim_matches(['(', ')'])
                .split(',')
                .collect::<Vec<_>>();
            let replacement = format!("#{new_id}");
            poles[1] = &replacement;
            let list = format!("({})", poles.join(","));
            args[2] = &list;
            with_fields(line, &args)
        });
        let end = changed.rfind("ENDSEC;").unwrap();
        changed.insert_str(end, &alias);
        let actual = CompoundDocument::parse(&changed, DocumentOptions::default(), None).unwrap();
        let result = Parsed::parse(&actual, solid, Default::default(), None);
        if exact {
            let parsed = result.unwrap();
            assert!(
                parsed
                    .curves
                    .iter()
                    .any(|c| c.poles.iter().any(|p| parsed.points[*p].id == new_id))
            );
            assert!(
                parsed
                    .tensors
                    .iter()
                    .any(|t| t.poles.iter().any(|p| parsed.points[*p].id == original))
            );
        } else {
            assert!(matches!(result, Err(Error::Unsupported)));
        }
        assert_eq!(actual.literal_input(), changed);
    }
    assert_eq!(document.literal_input(), text);
}

#[test]
fn actual_cap_and_edge_ownership_mutations_cannot_publish_a_closed_shell() {
    let text = fixture(false);
    let document = CompoundDocument::parse(&text, Default::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let records = document.records();
    let caps = parsed
        .faces
        .iter()
        .filter(|f| matches!(f.support, Support::Plane { .. }))
        .collect::<Vec<_>>();
    assert_eq!(caps.len(), 2);
    let cap_line = text
        .lines()
        .find(|line| line.starts_with(&format!("#{}=", caps[0].id)))
        .unwrap();
    let cap_fields = fields(cap_line);
    let duplicate = change_record(&text, caps[1].id, |line| {
        let mut args = fields(line);
        args[1] = cap_fields[1];
        args[2] = cap_fields[2];
        with_fields(line, &args)
    });
    assert_eq!(rejected(&duplicate), Error::Invalid);
    let Value::List(bounds) = &records[&caps[0].id][0].args[1] else {
        panic!()
    };
    let Value::Reference(bound) = bounds[0] else {
        panic!()
    };
    let Value::Reference(wire) = records[&bound][0].args[1] else {
        panic!()
    };
    let Value::List(coedges) = &records[&wire][0].args[1] else {
        panic!()
    };
    let ids = coedges
        .iter()
        .map(|c| match c {
            Value::Reference(id) => *id,
            _ => panic!(),
        })
        .collect::<Vec<_>>();
    let repeated = replace_field(
        &text,
        wire,
        1,
        &format!("(#{0},#{0},#{1},#{2})", ids[0], ids[2], ids[3]),
    );
    assert_eq!(rejected(&repeated), Error::Invalid);
    let mut reversed = text.clone();
    for id in &ids {
        reversed = change_record(&reversed, *id, |line| {
            let mut args = fields(line);
            args[4] = if args[4] == ".T." { ".F." } else { ".T." };
            with_fields(line, &args)
        });
    }
    let order = format!(
        "({})",
        ids.iter()
            .rev()
            .map(|id| format!("#{id}"))
            .collect::<Vec<_>>()
            .join(",")
    );
    reversed = replace_field(&reversed, wire, 1, &order);
    assert_eq!(rejected(&reversed), Error::Invalid);
    assert!(matches!(parsed.faces[0].support, Support::Tensor(_)));
    let [a, b] = [0, 1].map(|i| parsed.edges[parsed.faces[0].coedges[i].edge].id);
    let Value::Reference(geometry) = records[&a][0].args[3] else {
        panic!()
    };
    let reused = replace_field(&text, b, 3, &format!("#{geometry}"));
    assert_eq!(rejected(&reused), Error::Unsupported);
    let swapped = change_record(&text, a, |line| {
        let mut args = fields(line);
        args.swap(1, 2);
        with_fields(line, &args)
    });
    assert_eq!(rejected(&swapped), Error::Unsupported);
    assert_eq!(document.literal_input(), text);
}
