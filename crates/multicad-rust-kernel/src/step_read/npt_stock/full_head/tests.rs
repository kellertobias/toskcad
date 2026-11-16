use super::*;
fn candidate(quarter: bool) -> String {
    let source = super::super::tests::fixture(quarter);
    let doc = CompoundDocument::parse(&source, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let cap = parsed
        .faces
        .iter()
        .filter_map(|f| match f.support {
            Support::Plane { height } => Some((parsed.points[height].coordinates[2].value(), f.id)),
            _ => None,
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap()
        .1;
    let head = HeadInput::regular_hex(&parsed, cap, 18., 3.).unwrap();
    let proof = head::certify(&parsed, head, 1e-7, 1e-6, Default::default(), None).unwrap();
    crate::step::npt_stock::full_head::propose(&proof, Default::default(), None)
        .unwrap()
        .as_str()
        .to_owned()
}
#[test]
fn actual_full_clipped_d11_file_earns_independent_occupied_head_with_derived_cap() {
    for quarter in [false, true] {
        let actual = candidate(quarter);
        let packet = DecodedFile::parse(&actual, Default::default(), None).unwrap();
        assert!(packet.derived_auxiliary_text().contains("DERIVED"));
        assert_eq!(packet.actual_file(), actual);
        for b in packet.component_bindings() {
            assert_eq!(
                &actual[b.actual_range()],
                &packet.derived_auxiliary_text()[b.derived_range()]
            );
        }
        let aux = CompoundDocument::parse(
            packet.derived_auxiliary_text(),
            DocumentOptions::default(),
            None,
        )
        .unwrap();
        let parsed = Parsed::parse(
            &aux,
            aux.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let proof = packet.certify(&parsed, 1e-7, 1e-6, None).unwrap();
        assert!(proof.material().capped_thread().volume_mm3()[0] > 0.);
        assert!(proof.material().minimum_wall_mm() > 1e-6);
        assert_eq!(proof.material().head().faces.len(), 8);
        assert!(proof.material().join_sewing_mm() <= 1e-7 / 8.);
        assert_ne!(
            packet.cap_chart().physical_support_entity,
            packet.cap_chart().derived_support_entity
        );
        let actual_doc =
            CompoundDocument::parse(&actual, DocumentOptions::default(), None).unwrap();
        assert!(
            !actual_doc
                .records()
                .contains_key(&proof.material().head().auxiliary_cap_entity)
        );
        assert!(proof.receipt().aggregate_work <= 272_000_000);
        assert!(proof.receipt().accounted_peak_bytes <= 64 * 1024 * 1024);
        eprintln!(
            "ACTUAL_CLIPPED_FULL_HEAD quarter={quarter} packet={:?} occupied={:?} chart={:?}",
            packet.receipt(),
            proof.receipt(),
            packet.cap_chart()
        );
    }
}

fn occupied(
    input: &str,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Receipt> {
    let packet = DecodedFile::parse(input, options, interrupt)?;
    let aux = CompoundDocument::parse(
        packet.derived_auxiliary_text(),
        DocumentOptions::default(),
        interrupt,
    )
    .map_err(neutral_error)?;
    let parsed = Parsed::parse(
        &aux,
        aux.represented_members()[0].item_entity(),
        Default::default(),
        interrupt,
    )?;
    Ok(packet.certify(&parsed, 1e-7, 1e-6, interrupt)?.receipt())
}
fn change(input: &str, id: u32, body: &str) -> String {
    let prefix = format!("#{id}=");
    let start = input.find(&prefix).unwrap();
    let end = start + input[start..].find(';').unwrap() + 1;
    format!("{}{}{};{}", &input[..start], prefix, body, &input[end..])
}
#[test]
fn actual_full_clipped_file_resource_frontiers_and_cancel_retry_preserve_input() {
    use std::cell::Cell;
    let input = candidate(false);
    let before = input.clone();
    let good = occupied(&input, Default::default(), None).unwrap();
    for maximum_work in [good.aggregate_work - 1, good.aggregate_work] {
        let result = occupied(
            &input,
            Options {
                maximum_work,
                ..Default::default()
            },
            None,
        );
        if maximum_work < good.aggregate_work {
            assert!(matches!(result, Err(Error::Limit)));
        } else {
            assert_eq!(result.unwrap().aggregate_work, good.aggregate_work);
        }
    }
    for maximum_bytes in [good.accounted_peak_bytes - 1, good.accounted_peak_bytes] {
        let result = occupied(
            &input,
            Options {
                maximum_bytes,
                ..Default::default()
            },
            None,
        );
        if maximum_bytes < good.accounted_peak_bytes {
            assert!(matches!(result, Err(Error::Limit)));
        } else {
            assert_eq!(
                result.unwrap().accounted_peak_bytes,
                good.accounted_peak_bytes
            );
        }
    }
    let packet = DecodedFile::parse(&input, Default::default(), None).unwrap();
    let aux = CompoundDocument::parse(
        packet.derived_auxiliary_text(),
        DocumentOptions::default(),
        None,
    )
    .unwrap();
    let parsed = Parsed::parse(
        &aux,
        aux.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let count = Cell::new(0usize);
    let counter = || {
        count.set(count.get() + 1);
        false
    };
    packet.certify(&parsed, 1e-7, 1e-6, Some(&counter)).unwrap();
    let total = count.get();
    for at in [1, total / 2, total] {
        let calls = Cell::new(0usize);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= at
        };
        assert!(matches!(
            packet.certify(&parsed, 1e-7, 1e-6, Some(&cancel)),
            Err(Error::Cancelled)
        ));
        assert_eq!(
            packet
                .certify(&parsed, 1e-7, 1e-6, None)
                .unwrap()
                .receipt()
                .aggregate_work,
            good.aggregate_work
        );
    }
    assert!(matches!(
        DecodedFile::parse(
            &input,
            Options {
                additional_resident_bytes: usize::MAX,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        DecodedFile::parse(
            &input,
            Options {
                additional_resident_bytes: 64 * 1024 * 1024,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(input, before);
}
#[test]
fn actual_full_clipped_file_refuses_subulp_head_alias_and_nonlinear_edge() {
    use super::super::records::{list, reference};
    use crate::step::pipe_hardware::numbers as exact;
    let input = candidate(false);
    let doc = CompoundDocument::parse(&input, DocumentOptions::default(), None).unwrap();
    let records = doc.records();
    let (&line_id, line) = records
        .iter()
        .find_map(|(id, e)| matches!(e.as_slice(),[e] if e.name=="LINE").then_some((id, &e[0])))
        .unwrap();
    let [_, origin, vector] = line.args.as_slice() else {
        panic!()
    };
    let origin = reference(origin).unwrap();
    let vector = reference(vector).unwrap();
    let [_, direction, _] = records[&vector][0].args.as_slice() else {
        panic!()
    };
    let direction = reference(direction).unwrap();
    let e = &records[&direction][0];
    let tokens = crate::step_read::part21::numeric_spellings(e).unwrap();
    let token = tokens.iter().find(|s| s.len() > 20).unwrap();
    let (mantissa, exponent) = token.split_once('E').unwrap();
    let exponent = exponent.parse::<i32>().unwrap();
    let mantissa = mantissa.strip_suffix('.').unwrap();
    let alias = format!("{mantissa}1.E{}", exponent - 1);
    let body = e.source.unwrap().replacen(*token, &alias, 1);
    let altered = change(&input, direction, &body);
    let alias_doc = CompoundDocument::parse(&altered, DocumentOptions::default(), None).unwrap();
    assert_eq!(
        e.args,
        alias_doc.records()[&direction][0].args,
        "actual new decimal must alias the same binary cache"
    );
    assert!(matches!(
        occupied(&altered, Default::default(), None),
        Err(Error::Unsupported)
    ));
    // A genuine bowed cubic has the SAME actual finite edge endpoint IDs.
    // Its interior controls differ; checking endpoints alone must not admit it.
    let edge = records
        .values()
        .find_map(|e| match e.as_slice() {
            [e] if e.name == "EDGE_CURVE"
                && matches!(&e.args[3],Value::Reference(c) if *c==line_id) =>
            {
                Some(e)
            }
            _ => None,
        })
        .unwrap();
    let end_vertex = reference(&edge.args[2]).unwrap();
    let end = reference(&records[&end_vertex][0].args[1]).unwrap();
    let coordinates = |id: u32| {
        let p = list(&records[&id][0].args[1]).unwrap();
        p.iter()
            .map(|v| {
                if let Value::Number(n) = v {
                    *n
                } else {
                    panic!()
                }
            })
            .collect::<Vec<_>>()
    };
    let a = coordinates(origin);
    let z = coordinates(end);
    let max = *records.keys().next_back().unwrap();
    let c1 = max + 1;
    let c2 = max + 2;
    let mut bowed = change(
        &input,
        line_id,
        &format!(
            "B_SPLINE_CURVE_WITH_KNOTS('',3,(#{origin},#{c1},#{c2},#{end}),.UNSPECIFIED.,.F.,.F.,(4,4),(0.,1.),.UNSPECIFIED.)"
        ),
    );
    let mut extra = String::new();
    for (id, t) in [(c1, 1. / 3.), (c2, 2. / 3.)] {
        let p = std::array::from_fn::<_, 3, _>(|i| {
            a[i] * (1. - t) + z[i] * t + if i == 2 { 0.1 } else { 0. }
        });
        extra.push_str(&format!(
            "#{id}=CARTESIAN_POINT('',({},{},{}));\n",
            exact::literal(p[0]).unwrap(),
            exact::literal(p[1]).unwrap(),
            exact::literal(p[2]).unwrap()
        ));
    }
    let at = bowed.rfind("ENDSEC;").unwrap();
    bowed.insert_str(at, &extra);
    let bowed_doc = CompoundDocument::parse(&bowed, DocumentOptions::default(), None).unwrap();
    assert_eq!(
        records[&end_vertex][0].source,
        bowed_doc.records()[&end_vertex][0].source
    );
    assert!(matches!(
        occupied(&bowed, Default::default(), None),
        Err(Error::Unsupported)
    ));
    assert!(occupied(&input, Default::default(), None).is_ok());
}
#[test]
fn actual_full_clipped_file_rejects_detached_topology_foreign_packet_and_unit_change() {
    use super::super::records::reference;
    let input = candidate(false);
    let doc = CompoundDocument::parse(&input, DocumentOptions::default(), None).unwrap();
    let records = doc.records();
    let annulus = records
        .values()
        .find_map(|e| match e.as_slice() {
            [e] if e.name == "ADVANCED_FACE"
                && matches!(&e.args[1],Value::List(b) if b.len()==2) =>
            {
                Some(e)
            }
            _ => None,
        })
        .unwrap();
    let Value::List(bounds) = &annulus.args[1] else {
        panic!()
    };
    let bound = reference(&bounds[1]).unwrap();
    let wire = reference(&records[&bound][0].args[1]).unwrap();
    let Value::List(uses) = &records[&wire][0].args[1] else {
        panic!()
    };
    let use_id = reference(&uses[0]).unwrap();
    let e = &records[&use_id][0];
    let edge = reference(&e.args[3]).unwrap();
    let sense = if e.args[4] == Value::Enum("T") {
        ".F."
    } else {
        ".T."
    };
    let reversed = change(
        &input,
        use_id,
        &format!("ORIENTED_EDGE('',*,*,#{edge},{sense})"),
    );
    assert!(occupied(&reversed, Default::default(), None).is_err());
    let unit = input.replacen("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT($,.METRE.)", 1);
    assert!(occupied(&unit, Default::default(), None).is_err());
    let packet = DecodedFile::parse(&input, Default::default(), None).unwrap();
    let foreign_text = packet.derived_auxiliary_text().to_owned();
    let foreign = CompoundDocument::parse(&foreign_text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &foreign,
        foreign.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    assert!(matches!(
        packet.certify(&parsed, 1e-7, 1e-6, None),
        Err(Error::Invalid)
    ));
}
