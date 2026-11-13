use super::*;
fn fixture() -> String {
    let source = super::super::tests::fixture(1, false, false, 1);
    crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap()
}
#[test]
fn exact_scalar_identity_is_distinct_from_binary64_and_resolves_axis_references() {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let b = body.cavity();
    assert!(b.exact_support_key(29).is_none());
    assert!(b.exact_support_keys_capacity_bytes() > 29 * std::mem::size_of::<ExactSupportKey>());
    for i in 0..29 {
        let key = b.exact_support_key(i).unwrap();
        assert!(key.retained_capacity_bytes() > std::mem::size_of::<ExactSupportKey>());
        assert_eq!(key, b.exact_support_key(i).unwrap());
    }
    // Actual file mutation: exact decimal radius differs while its binary64
    // value remains identical. Neither float equality nor entity IDs may merge it.
    let id = b
        .faces()
        .iter()
        .find(|f| matches!(f.support, Surface::Torus { .. }))
        .unwrap()
        .surface_entity;
    let parsed = crate::step_read::part21::parse(&text).unwrap();
    let entity = &parsed[&id][0];
    let token = crate::step_read::part21::numeric_spellings(entity).unwrap()[0];
    let value: f64 = token.parse().unwrap();
    let (mantissa, exponent) = token.split_once('E').unwrap();
    let different = format!("{mantissa}000000000000000000000000000000000000001E{exponent}");
    assert_eq!(value, different.parse::<f64>().unwrap());
    let raw = entity.source.unwrap();
    let changed = raw.replacen(token, &different, 1);
    let altered = text.replace(raw, &changed);
    let after = read_phillips_body(&altered, 1e-7, 1_000_000).unwrap();
    for i in 0..29 {
        let affected = b.faces()[i].surface_entity == id;
        assert_eq!(
            b.exact_support_key(i) != after.cavity().exact_support_key(i),
            affected
        );
    }
    // Labels and reference IDs themselves are not numeric support parameters.
    let mut entities = crate::step_read::part21::parse(&text).unwrap();
    let mut work = Work {
        remaining: 1_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let original = read(&Shell::new(&entities, 1e-7, true), id, &mut work).unwrap();
    entities.get_mut(&id).unwrap()[0].source = None;
    assert!(
        read(&Shell::new(&entities, 1e-7, true), id, &mut work)
            .unwrap()
            .is_none()
    );
    assert!(original.is_some());
}
#[test]
fn raw_numeric_tokens_exclude_labels_references_and_comments_and_preserve_exact_spelling() {
    let input = "THING('label 123 ''456',#789, /* 321 */ (+1.0E-2,-0.,2.))";
    let entity = Entity {
        source: Some(input),
        name: "THING",
        args: vec![],
    };
    assert_eq!(
        crate::step_read::part21::numeric_spellings(&entity).unwrap(),
        ["+1.0E-2", "-0.", "2."]
    );
    for input in [
        "THING('unclosed",
        "THING(#0)",
        "THING(/*unclosed",
        "THING(#)",
    ] {
        let entity = Entity {
            source: Some(input),
            name: "THING",
            args: vec![],
        };
        assert!(crate::step_read::part21::numeric_spellings(&entity).is_none());
    }
    let entity = Entity {
        source: None,
        name: "THING",
        args: vec![],
    };
    assert!(crate::step_read::part21::numeric_spellings(&entity).is_none());
}
#[test]
fn unavailable_or_malformed_authority_and_work_are_never_replaced_by_float_keys() {
    let text = fixture();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let id = body.cavity().faces()[0].surface_entity;
    let mut entities = crate::step_read::part21::parse(&text).unwrap();
    let axis = reference(&entities[&id][0].args[1]).unwrap();
    let coordinates: Vec<_> = entities[&axis][0].args[1..]
        .iter()
        .map(|v| reference(v).unwrap())
        .collect();
    let mut work = Work {
        remaining: 1_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let original = read(&Shell::new(&entities, 1e-7, true), id, &mut work)
        .unwrap()
        .unwrap();
    for entity in [axis, coordinates[0], coordinates[1], coordinates[2]] {
        let old = entities.get_mut(&entity).unwrap()[0].source.take();
        assert!(
            read(&Shell::new(&entities, 1e-7, true), id, &mut work)
                .unwrap()
                .is_none()
        );
        entities.get_mut(&entity).unwrap()[0].source = old;
    }
    for entity in coordinates {
        let old = entities[&entity][0].source.unwrap();
        let mut altered = old.to_owned();
        let token = crate::step_read::part21::numeric_spellings(&entities[&entity][0]).unwrap()[0];
        altered = altered.replacen(token, &format!("+{token}"), 1);
        let changed_text = text.replace(old, &altered);
        let changed_entities = crate::step_read::part21::parse(&changed_text).unwrap();
        assert_ne!(
            read(&Shell::new(&changed_entities, 1e-7, true), id, &mut work)
                .unwrap()
                .unwrap(),
            original
        );
    }
    let mut zero = Work {
        remaining: 0,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    assert!(matches!(
        read(&Shell::new(&entities, 1e-7, true), id, &mut zero),
        Err(ReadError::Limit)
    ));
    let stop = || true;
    zero.remaining = 1_000_000;
    zero.interrupt = Some(&stop);
    assert!(read(&Shell::new(&entities, 1e-7, true), id, &mut zero).is_err());
    entities.get_mut(&axis).unwrap()[0].args.pop();
    assert!(matches!(
        read(&Shell::new(&entities, 1e-7, true), id, &mut work),
        Err(ReadError::Invalid)
    ));
}
