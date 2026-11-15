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
