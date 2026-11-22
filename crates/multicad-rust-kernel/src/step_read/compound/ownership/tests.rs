use super::*;
#[test]
fn actual_neutral_flat_topology_has_one_owner_and_bounded_work() {
    let mut writer = crate::step::Writer::new();
    for x in [0., 4.] {
        writer
            .append_sphere(
                crate::sphere::Sphere::new(1.).unwrap(),
                [x, 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
    }
    let input = writer.finish().unwrap();
    let doc = CompoundDocument::parse(&input, Default::default(), None).unwrap();
    let proof = validate_document_topology_ownership(&doc, 1e-7, Default::default(), None).unwrap();
    assert!(std::ptr::eq(proof.document(), &doc));
    assert_eq!(proof.represented_member_count(), 2);
    assert!(proof.topology_entity_count() > 0);
    let exact = DocumentOwnershipOptions {
        maximum_work: proof.work_used(),
        maximum_query_bytes: proof.query_peak_bytes(),
        ..Default::default()
    };
    validate_document_topology_ownership(&doc, 1e-7, exact, None).unwrap();
    for o in [
        DocumentOwnershipOptions {
            maximum_work: exact.maximum_work - 1,
            ..Default::default()
        },
        DocumentOwnershipOptions {
            maximum_query_bytes: exact.maximum_query_bytes - 1,
            ..Default::default()
        },
        DocumentOwnershipOptions {
            maximum_members: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            validate_document_topology_ownership(&doc, 1e-7, o, None),
            Err(DocumentOwnershipError::Limit)
        ));
    }
    assert!(matches!(
        validate_document_topology_ownership(&doc, 1e-7, Default::default(), Some(&|| true)),
        Err(DocumentOwnershipError::Cancelled)
    ));
}

#[test]
fn actual_document_ownership_rejects_orphans_shared_topology_and_invalid_query_contracts() {
    let mut writer = crate::step::Writer::new();
    for x in [0., 4.] {
        writer
            .append_sphere(
                crate::sphere::Sphere::new(1.).unwrap(),
                [x, 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
    }
    let input = writer.finish().unwrap();
    let doc = CompoundDocument::parse(&input, Default::default(), None).unwrap();
    for tolerance in [0., f64::NAN] {
        assert!(matches!(
            validate_document_topology_ownership(&doc, tolerance, Default::default(), None),
            Err(DocumentOwnershipError::Invalid)
        ));
    }
    for o in [
        DocumentOwnershipOptions {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        DocumentOwnershipOptions {
            maximum_query_bytes: 0,
            ..Default::default()
        },
        DocumentOwnershipOptions {
            maximum_query_bytes: BYTES + 1,
            ..Default::default()
        },
        DocumentOwnershipOptions {
            maximum_members: 0,
            ..Default::default()
        },
        DocumentOwnershipOptions {
            maximum_members: 1025,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            validate_document_topology_ownership(&doc, 1e-7, o, None),
            Err(DocumentOwnershipError::Invalid)
        ));
    }
    let mut orphan = input.clone();
    let at = orphan.rfind("ENDSEC;").unwrap();
    orphan.insert_str(at, "#999999=VERTEX_POINT('',#1);\n");
    let extra = CompoundDocument::parse(&orphan, Default::default(), None).unwrap();
    assert!(matches!(
        validate_document_topology_ownership(&extra, 1e-7, Default::default(), None),
        Err(DocumentOwnershipError::Invalid)
    ));
    let roots = doc.represented_members();
    let a = &doc.records()[&roots[0].local_solid_entity()][0];
    let b = &doc.records()[&roots[1].local_solid_entity()][0];
    let Value::Reference(a_shell) = a.args[1] else {
        panic!("actual first shell");
    };
    let Value::Reference(b_shell) = b.args[1] else {
        panic!("actual second shell");
    };
    let original = b.source.unwrap();
    let changed = original.replacen(&format!("#{b_shell}"), &format!("#{a_shell}"), 1);
    let shared = input.replacen(original, &changed, 1);
    let shared = CompoundDocument::parse(&shared, Default::default(), None).unwrap();
    assert!(matches!(
        validate_document_topology_ownership(&shared, 1e-7, Default::default(), None),
        Err(DocumentOwnershipError::Invalid)
    ));
    let vertex_loop = doc
        .records()
        .values()
        .flatten()
        .find(|e| e.name == "VERTEX_LOOP")
        .unwrap()
        .source
        .unwrap();
    let complex = input.replacen(
        vertex_loop,
        &format!("({vertex_loop}REPRESENTATION_ITEM(''))"),
        1,
    );
    let syntactic = CompoundDocument::parse(&complex, Default::default(), None).unwrap();
    assert!(matches!(
        validate_document_topology_ownership(&syntactic, 1e-7, Default::default(), None),
        Err(DocumentOwnershipError::Unsupported)
    ));
}
