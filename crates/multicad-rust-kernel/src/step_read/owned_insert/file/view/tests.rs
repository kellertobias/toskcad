use super::*;
#[test]
fn actual_ordinary_member_view_preserves_definitions_and_uses_generic_geometry_reader() {
    let original = super::super::tests::compound();
    // This real header includes a quoted fake marker, escaped quote and layout
    // comments. The view must preserve it through the ACTUAL section boundary.
    let text = original
        .replace("MultiCAD analytic solid", "before \nDATA; '' after")
        .replace("\nDATA;", "\n/* marker DATA; */ DATA /* layout */ ;");
    let file = File::parse(&text, Options::default(), None).unwrap();
    let shell = Shell::new(&file.entities, 1e-7, false);
    let sphere = file
        .represented_members()
        .iter()
        .find(|m| shell.member_faces(m.local_solid_entity()).unwrap().len() < 3)
        .unwrap()
        .item_entity();
    let view = file
        .local_member_document(sphere, Options::default(), None)
        .unwrap();
    assert!(std::ptr::eq(view.source(), &file));
    assert_eq!(view.member().item_entity(), sphere);
    assert!(view.work_used() <= 16_000_000 && view.query_budget_bytes() <= 64 * 1024 * 1024);
    assert!(view.as_str().len() < text.len() / 10);
    let entities = part21::parse(view.as_str()).unwrap();
    for (&id, after) in &entities {
        if let Some(before) = file.entities.get(&id) {
            assert_eq!(before.len(), after.len());
            for (a, b) in before.iter().zip(after) {
                assert_eq!(a.source, b.source);
            }
        }
    }
    let (decoded, _) = crate::step_read::read_import_solid(view.as_str(), 1e-7).unwrap();
    assert_eq!(decoded.solid_count, 1);
    assert_eq!(decoded.bounds_mm, [[19., -1., -1.], [21., 1., 1.]]);
    assert!((decoded.volume_mm3 - 4. * std::f64::consts::PI / 3.).abs() < 1e-12);
    let retry = file
        .local_member_document(sphere, Options::default(), None)
        .unwrap();
    assert_eq!(view.as_str(), retry.as_str());
}
#[test]
fn mapped_view_resource_cancellation_and_high_original_ids_never_drop_a_transform() {
    let (text, _, _, _) = super::super::tests::mapped_file();
    let file = File::parse(&text, Options::default(), None).unwrap();
    assert!(matches!(
        file.local_member_document(1000007, Options::default(), None),
        Err(Error::Unsupported)
    ));
    let text = super::super::tests::compound();
    let file = File::parse(&text, Options::default(), None).unwrap();
    let id = file.represented_members()[0].item_entity();
    assert!(matches!(
        file.local_member_document(u32::MAX, Options::default(), None),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        file.local_member_document(
            id,
            Options {
                maximum_work: 0,
                ..Options::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        file.local_member_document(
            id,
            Options {
                maximum_work: file.work_used(),
                ..Options::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        file.local_member_document(
            id,
            Options {
                retained_budget_bytes: file.retained_budget_bytes(),
                ..Options::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let polls = std::cell::Cell::new(0);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    let expected = file
        .local_member_document(id, Options::default(), Some(&count))
        .unwrap();
    let total = polls.get();
    for at in [1, total] {
        let polls = std::cell::Cell::new(0);
        let stop = || {
            polls.set(polls.get() + 1);
            polls.get() == at
        };
        assert!(matches!(
            file.local_member_document(id, Options::default(), Some(&stop)),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(
        file.local_member_document(id, Options::default(), None)
            .unwrap()
            .as_str(),
        expected.as_str()
    );
    let high = text.replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        "#4294967295=APPLICATION_CONTEXT('unused high ID');\nENDSEC;\nEND-ISO-10303-21;",
        1,
    );
    let file = File::parse(&high, Options::default(), None).unwrap();
    let view = file
        .local_member_document(id, Options::default(), None)
        .unwrap();
    assert!(crate::step_read::read_import_solid(view.as_str(), 1e-7).is_ok());
}
