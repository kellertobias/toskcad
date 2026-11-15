use super::*;
#[test]
fn actual_verified_candidates_append_and_reopen_without_an_authored_import_recipe() {
    for quarter in [false, true] {
        let source = super::super::tests::source(quarter);
        let curves = CurvePlan::new(&source, Default::default(), None).unwrap();
        let plan = SurfacePlan::new(curves, Default::default(), None).unwrap();
        let verified = plan.verify_candidate(Default::default(), None).unwrap();
        assert!(std::ptr::eq(verified.source_plan(), &plan));
        assert!(verified.literal_material().world_volume_mm3()[0] > 0.);
        let receipt = verified.receipt();
        assert_eq!(
            receipt.aggregate_work_used(),
            receipt.encoding_work_used()
                + receipt.document_work_used()
                + receipt.material_work_used()
                + receipt.correspondence_work_used()
                + receipt.orchestration_work_used()
        );
        assert!(
            receipt.encoding_work_used() <= 16_000_000
                && receipt.material_work_used() <= 64_000_000
                && receipt.aggregate_work_used() <= 96_000_000
        );
        assert!(receipt.accounted_peak_bytes() <= 64 * 1024 * 1024);
        assert!(
            receipt.position_error_mm() <= 1e-7
                && receipt
                    .original_derivative_error_mm()
                    .iter()
                    .all(|e| *e <= 1e-6)
        );
        assert!(receipt.edge_derivative_error_mm() <= 1e-6);
        assert_eq!(receipt.append_work_used(), 0);
        assert!(receipt.literal_retained_bytes() > 0);
        assert!(
            receipt.encoding_peak_bytes() > 0
                && receipt.document_peak_bytes() > 0
                && receipt.material_peak_bytes() > 0
                && receipt.correspondence_peak_bytes() > 0
        );
        assert!(receipt.local_correspondence().query_work_used() > 0);
        let mut writer = super::super::super::Writer::new();
        let appended = writer
            .append_verified_source_npt_candidate(&verified, Default::default(), None)
            .unwrap();
        assert!(appended.append_work_used() > 0);
        assert_eq!(
            appended.aggregate_work_used(),
            receipt.aggregate_work_used() + appended.append_work_used() + 256
        );
        assert_eq!(
            appended.aggregate_work_used(),
            appended.encoding_work_used()
                + appended.document_work_used()
                + appended.material_work_used()
                + appended.correspondence_work_used()
                + appended.orchestration_work_used()
        );
        let text = writer.finish().unwrap();
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let selected = Selected::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let owned = selected.retain(Default::default(), None).unwrap();
        let correspondence = plan
            .correspond_literal(&owned, Default::default(), None)
            .unwrap();
        assert!(correspondence.position_error_mm() <= 1e-7);
        eprintln!(
            "NPT VERIFIED quarter={quarter} verifywork={} appendwork={} total={} peak={} encoding={} document={} material={} correspondence={}",
            receipt.aggregate_work_used(),
            appended.append_work_used(),
            appended.aggregate_work_used(),
            appended.accounted_peak_bytes(),
            receipt.encoding_work_used(),
            receipt.document_work_used(),
            receipt.material_work_used(),
            receipt.correspondence_work_used()
        );
    }
}
fn prefix() -> super::super::super::Writer {
    let mut writer = super::super::super::Writer::new();
    writer
        .append_box(
            [crate::Rational::new(5, 1).unwrap(); 3],
            [40., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    writer
}
#[test]
fn actual_verified_replay_obeys_explicit_64m_and_cancel_limits_atomically() {
    let source = super::super::tests::source(false);
    let plan = SurfacePlan::new(
        CurvePlan::new(&source, Default::default(), None).unwrap(),
        Default::default(),
        None,
    )
    .unwrap();
    let candidate = plan.verify_candidate(Default::default(), None).unwrap();
    let expected = prefix().finish().unwrap();
    for options in [
        VerifiedOptions {
            maximum_work: 64_000_000,
            ..Default::default()
        },
        VerifiedOptions {
            maximum_work: candidate.receipt().aggregate_work_used() - 1,
            ..Default::default()
        },
        VerifiedOptions {
            maximum_bytes: candidate.receipt().accounted_peak_bytes() - 1,
            ..Default::default()
        },
    ] {
        let mut writer = prefix();
        assert!(matches!(
            writer.append_verified_source_npt_candidate(&candidate, options, None),
            Err(Error::Limit)
        ));
        assert_eq!(writer.finish().unwrap(), expected);
    }
    let mut writer = prefix();
    assert!(matches!(
        writer.append_verified_source_npt_candidate(
            &candidate,
            VerifiedOptions {
                position_error_mm: 1e-12,
                ..Default::default()
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    assert_eq!(writer.finish().unwrap(), expected);
    let calls = std::cell::Cell::new(0usize);
    let mut writer = prefix();
    writer
        .append_verified_source_npt_candidate(
            &candidate,
            Default::default(),
            Some(&|| {
                calls.set(calls.get() + 1);
                false
            }),
        )
        .unwrap();
    for checkpoint in [1, 2, calls.get() / 2, calls.get()] {
        let count = std::cell::Cell::new(0usize);
        let stop = || {
            count.set(count.get() + 1);
            count.get() >= checkpoint
        };
        let mut writer = prefix();
        assert!(matches!(
            writer.append_verified_source_npt_candidate(
                &candidate,
                Default::default(),
                Some(&stop)
            ),
            Err(Error::Cancelled)
        ));
        assert_eq!(writer.finish().unwrap(), expected);
    }
    // Reuse the same immutable proof after failed/cancelled transactions.
    let mut writer = prefix();
    writer
        .append_verified_source_npt_candidate(&candidate, Default::default(), None)
        .unwrap();
    let text = writer.finish().unwrap();
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    assert_eq!(document.represented_members().len(), 2);
    let member = document.represented_members()[1];
    let selected =
        Selected::parse(&document, member.item_entity(), Default::default(), None).unwrap();
    let owned = selected.retain(Default::default(), None).unwrap();
    assert!(
        plan.correspond_literal(&owned, Default::default(), None)
            .is_ok()
    );
}
#[test]
fn verified_pipeline_keeps_actual_phase_and_tolerance_frontiers_explicit() {
    let source = super::super::tests::source(false);
    let plan = SurfacePlan::new(
        CurvePlan::new(&source, Default::default(), None).unwrap(),
        Default::default(),
        None,
    )
    .unwrap();
    for options in [
        VerifiedOptions {
            maximum_work: 0,
            ..Default::default()
        },
        VerifiedOptions {
            maximum_work: 96_000_001,
            ..Default::default()
        },
        VerifiedOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        VerifiedOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
        VerifiedOptions {
            position_error_mm: f64::NAN,
            ..Default::default()
        },
        VerifiedOptions {
            original_derivative_error_mm: 0.,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            plan.verify_candidate(options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        VerifiedOptions {
            maximum_work: 2048,
            ..Default::default()
        },
        VerifiedOptions {
            maximum_work: 1024,
            ..Default::default()
        },
        VerifiedOptions {
            maximum_work: 2048 + plan.work_used() - 1,
            ..Default::default()
        },
        VerifiedOptions {
            maximum_bytes: PIPELINE_SCRATCH,
            ..Default::default()
        },
        VerifiedOptions {
            maximum_bytes: 8 * 1024 * 1024,
            ..Default::default()
        },
        VerifiedOptions {
            maximum_bytes: 30 * 1024 * 1024,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            plan.verify_candidate(options, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        plan.verify_candidate(
            VerifiedOptions {
                position_error_mm: 1e-12,
                ..Default::default()
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        plan.verify_candidate(Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
}
#[test]
fn two_verified_npt_members_with_an_ordinary_box_share_actual_encoder_budget() {
    let source = super::super::tests::source(false);
    let plan = SurfacePlan::new(
        CurvePlan::new(&source, Default::default(), None).unwrap(),
        Default::default(),
        None,
    )
    .unwrap();
    let candidate = plan.verify_candidate(Default::default(), None).unwrap();
    let mut writer = super::super::super::Writer::new();
    let first = writer
        .append_verified_source_npt_candidate(&candidate, Default::default(), None)
        .unwrap();
    writer
        .append_box(
            [crate::Rational::new(5, 1).unwrap(); 3],
            [40., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    let second = writer
        .append_verified_source_npt_candidate(&candidate, Default::default(), None)
        .unwrap();
    assert!(first.append_work_used() <= 16_000_000 && second.append_work_used() <= 16_000_000);
    assert!(second.aggregate_work_used() <= 96_000_000);
    let text = writer.finish().unwrap();
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    assert_eq!(document.represented_members().len(), 3);
    for index in [0, 2] {
        let selected = Selected::parse(
            &document,
            document.represented_members()[index].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let owned = selected.retain(Default::default(), None).unwrap();
        let receipt = plan
            .correspond_literal(&owned, Default::default(), None)
            .unwrap();
        assert!(receipt.position_error_mm() <= 1e-7);
    }
    eprintln!(
        "NPT COMPOUND VERIFIED first_append={} second_append={} second_total={} bytes={}",
        first.append_work_used(),
        second.append_work_used(),
        second.aggregate_work_used(),
        second.accounted_peak_bytes()
    );
}
#[test]
fn verified_cancellation_reaches_actual_document_and_decoder_phase_boundaries() {
    let source = super::super::tests::source(false);
    let plan = SurfacePlan::new(
        CurvePlan::new(&source, Default::default(), None).unwrap(),
        Default::default(),
        None,
    )
    .unwrap();
    let encoded_calls = std::cell::Cell::new(0usize);
    let file = write_source_candidate(
        &plan,
        16_000_000,
        64 * 1024 * 1024,
        Some(&|| {
            encoded_calls.set(encoded_calls.get() + 1);
            false
        }),
    )
    .unwrap();
    let document_calls = std::cell::Cell::new(0usize);
    CompoundDocument::parse(
        file.as_str(),
        DocumentOptions::default(),
        Some(&|| {
            document_calls.set(document_calls.get() + 1);
            false
        }),
    )
    .unwrap();
    // Pipeline orchestration, encoding, then its publication charge precede
    // the actual document. Its publication charge precedes the selected reader.
    for checkpoint in [
        encoded_calls.get() + 3,
        encoded_calls.get() + document_calls.get() + 4,
    ] {
        let count = std::cell::Cell::new(0usize);
        let stop = || {
            count.set(count.get() + 1);
            count.get() >= checkpoint
        };
        assert!(matches!(
            plan.verify_candidate(Default::default(), Some(&stop)),
            Err(Error::Cancelled)
        ));
        assert_eq!(count.get(), checkpoint);
    }
}
#[test]
fn genuine_unrepresentable_flat_land_cannot_promote_a_generated_packet() {
    use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredThread};
    let standard = super::super::tests::source(false);
    let mut spec = standard.source().spec();
    let AuthoredSection::NptFlat { depth_mm, .. } = spec.section else {
        unreachable!()
    };
    spec.section = AuthoredSection::NptFlat {
        depth_mm,
        flat_width_mm: 1e-16,
    };
    let source = AuthoredThread::new(spec, 1e-7, 1e-6)
        .unwrap()
        .stock_clip_plan(Default::default(), None)
        .unwrap();
    let plan = SurfacePlan::new(
        CurvePlan::new(&source, Default::default(), None).unwrap(),
        Default::default(),
        None,
    )
    .unwrap();
    // This is a real accepted analytic source, not a forged SurfacePlan. Its
    // positive crest land is below representable physical pole separation.
    // Independent literal geometry must refuse the degenerate candidate.
    assert!(matches!(
        plan.verify_candidate(Default::default(), None),
        Err(Error::Unresolved)
    ));
}
