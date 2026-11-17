use super::*;
fn check_tokens(document: &CompoundDocument<'_>, view: &SelectedSourceDocument<'_, '_>) {
    let entities = part21::parse(view.as_str()).unwrap();
    for (id, after) in &entities {
        if let Some(before) = document.source.entities.get(id) {
            assert_eq!(before.len(), after.len());
            for (a, b) in before.iter().zip(after) {
                assert_eq!(a.source, b.source);
            }
        }
    }
    let parsed = File::parse(view.as_str(), Options::default(), None).unwrap();
    assert_eq!(parsed.represented_members().len(), 1);
    let member = parsed.represented_members()[0];
    assert_eq!(member.item_entity(), view.member().local_solid_entity());
    assert_eq!(
        member.local_solid_entity(),
        view.member().local_solid_entity()
    );
    assert!(member.mapping_link().is_none());
}
#[test]
fn selected_source_document_retains_actual_mapping_and_only_local_tokens() {
    let (text, solid, _, _) = super::super::tests::mapped_file();
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let original = document.represented_members()[0];
    let source = document
        .local_source_document(original.item_entity(), Options::default(), None)
        .unwrap();
    assert!(std::ptr::eq(source.document(), &document));
    assert_eq!(source.member(), original);
    assert_eq!(source.mapping_link(), original.mapping_link());
    assert_eq!(source.member().local_solid_entity(), solid);
    assert!(source.mapping_link().is_some());
    check_tokens(&document, &source);
    assert!(matches!(
        document.local_member_document(original.item_entity(), Options::default(), None),
        Err(Error::Unsupported)
    ));
    let child = parse(source.as_str(), Options::default(), None).unwrap();
    let bindings = child.support_bindings(Options::default(), None).unwrap();
    let local = bindings.material_bounds(Options::default(), None).unwrap();
    assert_eq!(local.bounds_mm()[2], [0., 4.]);
    let packet = document
        .member_packet(original.item_entity(), Options::default(), None)
        .unwrap();
    let mapping = packet
        .parsed()
        .physical_mapping(
            source.mapping_link().unwrap().operator_entity(),
            Options::default(),
            None,
        )
        .unwrap();
    for ([lo, hi], exact) in mapping.origin_mm().into_iter().zip([10., 20., 30.]) {
        assert!(lo <= exact && exact <= hi && hi - lo < 1e-12);
    }
    let world = mapping.map_box_mm(local.bounds_mm()).unwrap();
    assert!(world[2][0] > 29. && world[2][1] > 33.);
}
#[test]
fn selected_source_document_ordinary_mixed_members_keep_exact_old_payload() {
    let text = super::super::tests::compound();
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    assert_eq!(document.represented_members().len(), 3);
    for &member in document.represented_members() {
        let old = document
            .local_member_document(member.item_entity(), Options::default(), None)
            .unwrap();
        let source = document
            .local_source_document(member.item_entity(), Options::default(), None)
            .unwrap();
        assert_eq!(source.as_str(), old.as_str());
        assert_eq!(source.member(), old.member());
        assert!(source.mapping_link().is_none());
        check_tokens(&document, &source);
    }
}
#[test]
fn selected_source_document_limits_cancel_and_retry_preserve_real_document() {
    let (text, _, _, _) = super::super::tests::mapped_file();
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let member = document.represented_members()[0];
    let options = Options::default();
    let polls = std::cell::Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    let source = document
        .local_source_document(member.item_entity(), options, Some(&count))
        .unwrap();
    eprintln!(
        "SELECTED_LOCAL_SOURCE work={} query_bytes={} doc_work={} doc_bytes={}",
        source.work_used(),
        source.query_budget_bytes(),
        document.work_used(),
        document.retained_budget_bytes()
    );
    for bad in [
        Options {
            maximum_work: source.work_used() - 1,
            ..options
        },
        Options {
            retained_budget_bytes: source.query_budget_bytes() - 1,
            ..options
        },
    ] {
        assert!(matches!(
            document.local_source_document(member.item_entity(), bad, None),
            Err(Error::Limit)
        ));
    }
    let exact = document
        .local_source_document(
            member.item_entity(),
            Options {
                maximum_work: source.work_used(),
                retained_budget_bytes: source.query_budget_bytes(),
            },
            None,
        )
        .unwrap();
    assert_eq!(exact.as_str(), source.as_str());
    assert!(matches!(
        document.local_source_document(u32::MAX, options, None),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        document.local_source_document(
            member.item_entity(),
            Options {
                maximum_work: 0,
                ..options
            },
            None
        ),
        Err(Error::Invalid)
    ));
    for stop in [1, polls.get()] {
        let at = std::cell::Cell::new(0usize);
        let cancel = || {
            at.set(at.get() + 1);
            at.get() >= stop
        };
        assert!(matches!(
            document.local_source_document(member.item_entity(), options, Some(&cancel)),
            Err(Error::Cancelled)
        ));
        let retry = document
            .local_source_document(member.item_entity(), options, None)
            .unwrap();
        assert_eq!(retry.as_str(), source.as_str());
        assert_eq!(retry.member(), source.member());
    }
    assert_eq!(document.literal_input(), text);
}

#[test]
fn selected_source_document_genuine_mixed_mapped_member_has_no_unrelated_geometry() {
    let text = super::tests::mixed(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    assert_eq!(document.represented_members().len(), 2);
    let member = document
        .represented_members()
        .iter()
        .copied()
        .find(|m| m.mapping_link().is_some())
        .unwrap();
    let source = document
        .local_source_document(member.item_entity(), Options::default(), None)
        .unwrap();
    assert_eq!(source.member(), member);
    assert_eq!(source.mapping_link(), member.mapping_link());
    check_tokens(&document, &source);
    let child = parse(source.as_str(), Options::default(), None).unwrap();
    let bindings = child.support_bindings(Options::default(), None).unwrap();
    assert_eq!(
        bindings
            .material_bounds(Options::default(), None)
            .unwrap()
            .bounds_mm()[2],
        [0., 4.]
    );
    assert!(matches!(
        document.local_member_document(member.item_entity(), Options::default(), None),
        Err(Error::Unsupported)
    ));
}
