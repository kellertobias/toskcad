use super::*;
#[test]
fn actual_metric_candidates_append_and_bind_actual_rebased_complex_capable_payloads() {
    let source = super::super::tests::source(false, 0.);
    let frame = super::super::tests::occurrence(false);
    let candidate = export(&source, frame, Default::default(), None).unwrap();
    let mut writer = step::Writer::new();
    let first = writer
        .append_qualified_periodic_phillips(&candidate, Default::default(), None)
        .unwrap();
    eprintln!("Phi actual rebased append {first:?}");
    assert!(first.matched_entities > 0 && first.work_used > 0);
    let second = writer
        .append_qualified_periodic_phillips(&candidate, Default::default(), None)
        .unwrap();
    assert_eq!(first.matched_entities, second.matched_entities);
    let text = writer.finish().unwrap();
    let doc = read::CompoundDocument::parse(&text, Default::default(), None).unwrap();
    assert_eq!(doc.represented_members().len(), 2);
    for member in doc.represented_members() {
        let selected = doc
            .local_member_document(member.item_entity(), Default::default(), None)
            .unwrap();
        qualify_actual_file(&source, selected.as_str(), frame, Default::default(), None).unwrap();
    }
}
fn prior() -> step::Writer {
    let mut writer = step::Writer::new();
    writer
        .append_sphere(
            crate::sphere::Sphere::new(1.).unwrap(),
            [20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    writer
}
#[test]
fn actual_replay_work_byte_child_and_cancellation_frontiers_restore_prior_members() {
    let source = super::super::tests::source(false, 0.);
    let candidate = export(
        &source,
        super::super::tests::occurrence(false),
        Default::default(),
        None,
    )
    .unwrap();
    let original = prior().finish().unwrap();
    let polls = std::cell::Cell::new(0);
    let observe = || {
        polls.set(polls.get() + 1);
        false
    };
    let proof = prior()
        .append_qualified_periodic_phillips(&candidate, Default::default(), Some(&observe))
        .unwrap();
    let exact = AppendOptions {
        maximum_work: proof.work_used,
        maximum_bytes: proof.accounted_peak_bytes,
        maximum_encoding_work: proof.encoding_work,
    };
    prior()
        .append_qualified_periodic_phillips(&candidate, exact, None)
        .unwrap();
    for options in [
        AppendOptions {
            maximum_work: 1,
            ..Default::default()
        },
        AppendOptions {
            maximum_work: proof.work_used - 1,
            ..Default::default()
        },
        AppendOptions {
            maximum_bytes: 1,
            ..Default::default()
        },
        AppendOptions {
            maximum_bytes: proof.accounted_peak_bytes - 1,
            ..Default::default()
        },
        AppendOptions {
            maximum_encoding_work: 1,
            ..Default::default()
        },
        AppendOptions {
            maximum_encoding_work: 16_000_000,
            ..Default::default()
        },
        AppendOptions {
            maximum_encoding_work: proof.encoding_work - 1,
            ..Default::default()
        },
    ] {
        let mut writer = prior();
        assert_eq!(
            writer
                .append_qualified_periodic_phillips(&candidate, options, None)
                .unwrap_err(),
            Error::Limit
        );
        assert_eq!(writer.finish().unwrap(), original);
    }
    for options in [
        AppendOptions {
            maximum_work: 0,
            ..Default::default()
        },
        AppendOptions {
            maximum_work: 128_000_001,
            ..Default::default()
        },
        AppendOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        AppendOptions {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
        AppendOptions {
            maximum_encoding_work: 0,
            ..Default::default()
        },
        AppendOptions {
            maximum_encoding_work: 32_000_001,
            ..Default::default()
        },
    ] {
        let mut writer = prior();
        assert_eq!(
            writer
                .append_qualified_periodic_phillips(&candidate, options, None)
                .unwrap_err(),
            Error::Invalid
        );
        assert_eq!(writer.finish().unwrap(), original);
    }
    for boundary in [1, 100, polls.get()] {
        let mut writer = prior();
        let calls = std::cell::Cell::new(0);
        let stop = || {
            calls.set(calls.get() + 1);
            calls.get() >= boundary
        };
        assert_eq!(
            writer
                .append_qualified_periodic_phillips(&candidate, Default::default(), Some(&stop))
                .unwrap_err(),
            Error::Cancelled
        );
        assert_eq!(writer.finish().unwrap(), original);
    }
}
