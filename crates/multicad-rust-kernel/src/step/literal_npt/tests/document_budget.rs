use super::*;

#[test]
fn actual_four_member_file_keeps_document_child_limit_under_larger_aggregate() {
    let source = source(false, false);
    let verified = verify_exact_replay(&source, Default::default(), None).unwrap();
    let mut writer = Writer::new();
    for _ in 0..4 {
        writer
            .append_verified_literal_npt(&verified, Default::default(), None)
            .unwrap();
    }
    let packet = writer.finish().unwrap();
    let document = CompoundDocument::parse(
        &packet,
        DocumentOptions {
            maximum_work: WORK,
            retained_budget_bytes: 128 * 1024 * 1024,
        },
        None,
    )
    .unwrap();
    assert_eq!(document.represented_members().len(), 4);
    let bytes = document.retained_budget_bytes();
    // The actual packet and neutral AST fit the new aggregate together with
    // the real retained subject; the independent document child exceeds64MiB.
    assert!(bytes > BYTES);
    assert!(source.retained_budget_bytes() + bytes + FIXED < 128 * 1024 * 1024);
    drop(document);
    assert!(matches!(
        CompoundDocument::parse(
            &packet,
            DocumentOptions {
                maximum_work: WORK,
                retained_budget_bytes: BYTES,
            },
            None
        ),
        Err(crate::step_read::owned_insert::Error::Limit)
    ));
    assert!(matches!(
        verify_readback(&source, None, &packet, 0, Default::default(), None),
        Err(Error::Limit)
    ));
    // The same independent reader still accepts the actual single-member
    // file under the unchanged document ceiling.
    let receipt = verify_readback(
        &source,
        None,
        verified.standalone_step(),
        0,
        Default::default(),
        None,
    )
    .unwrap();
    assert_eq!(receipt.represented_member_count(), 1);
    assert!(receipt.actual_volume_mm3()[0] > 0.);
}
