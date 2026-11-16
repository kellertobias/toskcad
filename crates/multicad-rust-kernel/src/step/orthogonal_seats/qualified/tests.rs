use super::*;
use std::cell::Cell;
#[test]
fn closed_qualified_bytes_preserve_original_and_healed_mass_for_actual_poses() {
    let source = super::super::tests::source(true);
    for placement in [
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap(),
        StockClipPlacement::new(
            [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
            [11., -7., 5.],
        )
        .unwrap(),
        StockClipPlacement::new(
            [[1. + 2f64.powi(-42), 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            [11., -7., 5.],
        )
        .unwrap(),
    ] {
        let start = std::time::Instant::now();
        let q = qualify(&source, placement, QualificationOptions::default(), &|| {
            false
        })
        .unwrap();
        eprintln!(
            "qualified Angle elapsed={:?} total_work={} encoding={} peak={} bytes={}",
            start.elapsed(),
            q.work_used(),
            q.encoding_work_used(),
            q.peak_accounted_bytes(),
            q.as_str().len()
        );
        assert!(q.matches_source(&source, placement));
        assert_eq!(q.placement().matrix(), placement.matrix());
        assert_eq!(q.source_to_raw_correspondence().original_face_count(), 13);
        assert_eq!(q.raw_to_healed_correspondence().original_face_count(), 13);
        assert!(
            q.source_to_healed_position_error_mm() < 1e-7
                && q.source_to_healed_normalized_derivative_error_mm() < 1e-7
        );
        let local = source.volume_bounds_mm3();
        let determinant = placement.determinant_bounds();
        let volume = crate::tapered_thread::interval::I::new(local[0], local[1])
            * crate::tapered_thread::interval::I::new(determinant[0], determinant[1]);
        assert_eq!(
            q.original_source_volume_bounds_mm3(),
            [volume.lo, volume.hi]
        );
        assert!(q.healed_volume_bounds_mm3()[0] > 0.);
        assert_eq!(q.qualification_material_options().minimum_feature_mm, 1e-5);
        let other = StockClipPlacement::new(placement.matrix(), [99., 3., 4.]).unwrap();
        assert!(!q.matches_source(&source, other));
        assert!(!q.matches_source(&super::super::tests::source(false), placement));
        let owned = q.into_string();
        assert_eq!(owned.matches("CONICAL_SURFACE(").count(), 3);
    }
}
#[test]
fn global_qualification_frontiers_and_late_cancellation_are_atomic() {
    let source = super::super::tests::source(false);
    let placement =
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap();
    let q = qualify(&source, placement, QualificationOptions::default(), &|| {
        false
    })
    .unwrap();
    let exact = QualificationOptions {
        maximum_work: q.work_used(),
        maximum_bytes: q.peak_accounted_bytes(),
        ..QualificationOptions::default()
    };
    assert!(qualify(&source, placement, exact, &|| false).is_ok());
    assert!(matches!(
        qualify(
            &source,
            placement,
            QualificationOptions {
                maximum_work: exact.maximum_work - 1,
                ..exact
            },
            &|| false
        ),
        Err(QualificationError::Limit)
    ));
    assert!(matches!(
        qualify(
            &source,
            placement,
            QualificationOptions {
                maximum_bytes: exact.maximum_bytes - 1,
                ..exact
            },
            &|| false
        ),
        Err(QualificationError::Limit)
    ));
    assert!(matches!(
        qualify(
            &source,
            placement,
            QualificationOptions {
                maximum_work: 0,
                ..exact
            },
            &|| false
        ),
        Err(QualificationError::Invalid)
    ));
    let calls = Cell::new(0usize);
    qualify(&source, placement, QualificationOptions::default(), &|| {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    let stop = calls.get() - 5;
    calls.set(0);
    assert!(matches!(
        qualify(&source, placement, QualificationOptions::default(), &|| {
            calls.set(calls.get() + 1);
            calls.get() > stop
        }),
        Err(QualificationError::Cancelled)
    ));
    assert_eq!(
        qualify(&source, placement, exact, &|| false)
            .unwrap()
            .as_str(),
        q.as_str()
    );
}
#[test]
fn source_strain_and_reflection_cannot_be_published_as_qualified_material() {
    let source = super::super::tests::source(false);
    let strained = StockClipPlacement::new(
        [[1. + 2f64.powi(-24), 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [0.; 3],
    )
    .unwrap();
    assert!(matches!(
        qualify(&source, strained, QualificationOptions::default(), &|| {
            false
        }),
        Err(QualificationError::Unsupported)
    ));
    assert!(StockClipPlacement::new([[-1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).is_err());
    let identity =
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap();
    let invalid = QualificationOptions {
        material: read::MaterialOptions {
            transfer: read::TransferOptions {
                maximum_work: 16_000_001,
                ..read::TransferOptions::default()
            },
            ..read::MaterialOptions::default()
        },
        ..QualificationOptions::default()
    };
    assert!(matches!(
        qualify(&source, identity, invalid, &|| false),
        Err(QualificationError::Invalid)
    ));
}

fn prefix() -> Writer {
    let mut writer = Writer::new();
    writer
        .append_box(
            [crate::Rational::new(4, 1).unwrap(); 3],
            [-100.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    writer
}
#[test]
fn actual_compound_members_bind_exact_literals_and_independent_healed_material() {
    let source = super::super::tests::source(false);
    let placement =
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap();
    let q = qualify(&source, placement, QualificationOptions::default(), &|| {
        false
    })
    .unwrap();
    let other = StockClipPlacement::new(
        [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
        [100., 20., 4.],
    )
    .unwrap();
    let second = qualify(&source, other, QualificationOptions::default(), &|| false).unwrap();
    let mut writer = prefix();
    let a = writer
        .append_qualified_orthogonal_seats(&q, AppendOptions::default(), &|| false)
        .unwrap();
    let b = writer
        .append_qualified_orthogonal_seats(&second, AppendOptions::default(), &|| false)
        .unwrap();
    assert_ne!(a.solid_entity(), b.solid_entity());
    assert!(a.append_work_used() > 0);
    assert_eq!(a.work_used(), q.work_used() + a.append_work_used());
    assert!(a.peak_accounted_bytes() > 0);
    assert_eq!(
        a.original_source_volume_bounds_mm3(),
        q.original_source_volume_bounds_mm3()
    );
    assert_eq!(a.healed_volume_bounds_mm3(), q.healed_volume_bounds_mm3());
    assert!(
        a.source_to_healed_position_error_mm() < 1e-7
            && a.source_to_healed_normalized_derivative_error_mm() < 1e-7
            && a.raw_support_sewing_error_mm() < 1e-7
    );
    let text = writer.finish().unwrap();
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    assert_eq!(document.represented_members().len(), 3);
    for (qualified, receipt) in [(&q, &a), (&second, &b)] {
        let proof = qualified
            .verify_actual_member(
                &document,
                receipt.solid_entity(),
                VerificationOptions::default(),
                &|| false,
            )
            .unwrap();
        assert_eq!(proof.item_entity(), receipt.solid_entity());
        assert_eq!(proof.solid_entity(), receipt.solid_entity());
        assert!(proof.matched_entities() > 100);
        assert_eq!(
            proof.original_source_volume_bounds_mm3(),
            qualified.original_source_volume_bounds_mm3()
        );
        assert!(proof.healed_volume_bounds_mm3()[0] > 0.);
        assert!(
            proof
                .source_to_raw_correspondence()
                .maximum_position_error_mm()
                < 1e-7
                && proof
                    .raw_to_healed_correspondence()
                    .maximum_position_error_mm()
                    < 1e-7
        );
        assert!(
            proof.source_to_healed_position_error_mm() < 1e-7
                && proof.source_to_healed_normalized_derivative_error_mm() < 1e-7
        );
        assert_eq!(
            proof.aggregate_known_work(),
            qualified.work_used() + document.work_used() + proof.work_used()
        );
        assert_eq!(
            proof.actual_document_retained_bytes(),
            document.retained_budget_bytes()
        );
        assert!(proof.query_peak_accounted_bytes() > 0);
    }
    assert!(matches!(
        q.verify_actual_member(
            &document,
            b.solid_entity(),
            VerificationOptions::default(),
            &|| false
        ),
        Err(QualificationError::Unsupported)
    ));
    assert!(matches!(
        q.verify_actual_member(
            &document,
            document.represented_members()[0].item_entity(),
            VerificationOptions::default(),
            &|| false
        ),
        Err(QualificationError::Unsupported)
    ));
    let original = "4.9000000000000004E0";
    let alias = "4.9000000000000004000000000000001E0";
    assert_eq!(
        original.parse::<f64>().unwrap(),
        alias.parse::<f64>().unwrap()
    );
    let changed = text.replacen(original, alias, 1);
    assert_ne!(changed, text);
    let changed_document =
        CompoundDocument::parse(&changed, DocumentOptions::default(), None).unwrap();
    assert!(matches!(
        q.verify_actual_member(
            &changed_document,
            a.solid_entity(),
            VerificationOptions::default(),
            &|| false
        ),
        Err(QualificationError::Unsupported)
    ));
}
#[test]
fn actual_append_and_verification_frontiers_cancel_without_partial_publication() {
    let source = super::super::tests::source(false);
    let placement =
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap();
    let q = qualify(&source, placement, QualificationOptions::default(), &|| {
        false
    })
    .unwrap();
    let mut writer = prefix();
    let receipt = writer
        .append_qualified_orthogonal_seats(&q, AppendOptions::default(), &|| false)
        .unwrap();
    let text = writer.finish().unwrap();
    let append_exact = AppendOptions {
        maximum_work: receipt.work_used(),
        maximum_bytes: receipt.peak_accounted_bytes(),
    };
    let mut writer = prefix();
    assert!(
        writer
            .append_qualified_orthogonal_seats(&q, append_exact, &|| false)
            .is_ok()
    );
    assert_eq!(writer.finish().unwrap(), text);
    let baseline = prefix().finish().unwrap();
    let mut writer = prefix();
    assert!(matches!(
        writer.append_qualified_orthogonal_seats(
            &q,
            AppendOptions {
                maximum_work: append_exact.maximum_work - 1,
                ..append_exact
            },
            &|| false
        ),
        Err(QualificationError::Limit)
    ));
    assert_eq!(writer.finish().unwrap(), baseline);
    let mut writer = prefix();
    assert!(matches!(
        writer.append_qualified_orthogonal_seats(
            &q,
            AppendOptions {
                maximum_bytes: append_exact.maximum_bytes - 1,
                ..append_exact
            },
            &|| false
        ),
        Err(QualificationError::Limit)
    ));
    assert_eq!(writer.finish().unwrap(), baseline);
    let calls = Cell::new(0usize);
    let mut writer = prefix();
    writer
        .append_qualified_orthogonal_seats(&q, AppendOptions::default(), &|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    let stop = calls.get() - 1;
    calls.set(0);
    let mut writer = prefix();
    assert!(matches!(
        writer.append_qualified_orthogonal_seats(&q, AppendOptions::default(), &|| {
            calls.set(calls.get() + 1);
            calls.get() > stop
        }),
        Err(QualificationError::Cancelled)
    ));
    assert_eq!(writer.finish().unwrap(), baseline);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let proof = q
        .verify_actual_member(
            &document,
            receipt.solid_entity(),
            VerificationOptions::default(),
            &|| false,
        )
        .unwrap();
    let exact = VerificationOptions {
        maximum_work: proof.work_used(),
        maximum_bytes: proof.query_peak_accounted_bytes(),
        ..VerificationOptions::default()
    };
    assert!(
        q.verify_actual_member(&document, receipt.solid_entity(), exact, &|| false)
            .is_ok()
    );
    assert!(matches!(
        q.verify_actual_member(
            &document,
            receipt.solid_entity(),
            VerificationOptions {
                maximum_work: exact.maximum_work - 1,
                ..exact
            },
            &|| false
        ),
        Err(QualificationError::Limit)
    ));
    assert!(matches!(
        q.verify_actual_member(
            &document,
            receipt.solid_entity(),
            VerificationOptions {
                maximum_bytes: exact.maximum_bytes - 1,
                ..exact
            },
            &|| false
        ),
        Err(QualificationError::Limit)
    ));
    calls.set(0);
    q.verify_actual_member(
        &document,
        receipt.solid_entity(),
        VerificationOptions::default(),
        &|| {
            calls.set(calls.get() + 1);
            false
        },
    )
    .unwrap();
    let stop = calls.get() - 3;
    calls.set(0);
    assert!(matches!(
        q.verify_actual_member(
            &document,
            receipt.solid_entity(),
            VerificationOptions::default(),
            &|| {
                calls.set(calls.get() + 1);
                calls.get() > stop
            }
        ),
        Err(QualificationError::Cancelled)
    ));
    assert!(
        q.verify_actual_member(&document, receipt.solid_entity(), exact, &|| false)
            .is_ok()
    );
}
