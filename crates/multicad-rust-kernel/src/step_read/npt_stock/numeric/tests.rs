use super::*;
fn with_values<T>(values: [&str; 3], query: impl FnOnce([Scalar<'_>; 3]) -> T) -> T {
    let text = format!(
        "ISO-10303-21;HEADER;FILE_DESCRIPTION(('Test'),'2;1');FILE_NAME('x','1970-01-01T00:00:00',(''),(''),'','','');FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;#1=CARTESIAN_POINT('',({},{},{}));ENDSEC;END-ISO-10303-21;",
        values[0], values[1], values[2]
    );
    let records = part21::parse(&text).unwrap();
    let entity = &records[&1][0];
    let Value::List(coords) = &entity.args[1] else {
        panic!()
    };
    let mut budget = Budget::new(Default::default(), None).unwrap();
    let mut numbers = Numbers::new(entity, &mut budget).unwrap();
    let values = std::array::from_fn(|i| numbers.scalar(&coords[i], &mut budget).unwrap());
    numbers.finish().unwrap();
    query(values)
}
#[test]
fn actual_decimal_chart_translation_is_exact_even_when_caches_coincide() {
    for (values, expected) in [
        (["-5.1964789195629741", "1.", "-4.1964789195629741"], true),
        (
            ["-1.21110570427484565", "2.", "7.8889429572515435E-1"],
            true,
        ),
        (["9999999999999999999999999", "1", "1E25"], true),
        (["1E-30", "-1E-30", "0."], true),
        (["-1E25", "-1", "-10000000000000000000000001"], true),
        (
            ["1.0000000000000000001", "1", "2.0000000000000000002"],
            false,
        ),
        (["0.1", "0.", "0.10"], true),
        (["1E-320", "0.", "0."], false),
        (["0.", "1E-320", "1E-320"], true),
    ] {
        with_values(values, |[a, b, c]| {
            let mut budget = Budget::new(Default::default(), None).unwrap();
            assert_eq!(
                translated_equal(a, b, c, &mut budget).unwrap(),
                expected,
                "{values:?}"
            );
            assert_eq!(budget.bytes, SCRATCH);
        });
    }
}
#[test]
fn decimal_chart_translation_honors_resource_cancel_and_retry() {
    with_values(
        ["-5.1964789195629741", "1.", "-4.1964789195629741"],
        |[a, b, c]| {
            let mut good = Budget::new(Default::default(), None).unwrap();
            assert!(translated_equal(a, b, c, &mut good).unwrap());
            let work = good.used;
            let bytes = good.peak;
            for options in [
                Options {
                    maximum_work: work - 1,
                    ..Default::default()
                },
                Options {
                    maximum_bytes: bytes - 1,
                    ..Default::default()
                },
            ] {
                let mut limited = Budget::new(options, None).unwrap();
                assert_eq!(translated_equal(a, b, c, &mut limited), Err(Error::Limit));
            }
            let polls = std::cell::Cell::new(0);
            let cancel = || {
                polls.set(polls.get() + 1);
                polls.get() == 8
            };
            let mut budget = Budget::new(Default::default(), Some(&cancel)).unwrap();
            assert_eq!(
                translated_equal(a, b, c, &mut budget),
                Err(Error::Cancelled)
            );
            let mut retry = Budget::new(Default::default(), None).unwrap();
            assert!(translated_equal(a, b, c, &mut retry).unwrap());
            assert_eq!(retry.used, work);
        },
    );
    with_values(["1E-100000", "1", "1"], |[a, b, c]| {
        assert_eq!(
            translated_equal(a, b, c, &mut Budget::new(Default::default(), None).unwrap()),
            Err(Error::Limit)
        );
    });
}
