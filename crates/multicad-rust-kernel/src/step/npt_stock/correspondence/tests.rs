use super::*;
use crate::step_read::{
    npt_stock::Selected,
    owned_insert::{CompoundDocument, DocumentOptions},
};
#[test]
fn actual_serialized_literal_nets_edges_and_owner_graph_correspond_to_original_source() {
    for quarter in [false, true] {
        let source = super::super::tests::source(quarter);
        let curves = CurvePlan::new(&source, Default::default(), None).unwrap();
        let plan = SurfacePlan::new(curves, Default::default(), None).unwrap();
        let file = write_source_candidate(&plan, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let document =
            CompoundDocument::parse(file.as_str(), DocumentOptions::default(), None).unwrap();
        let selected = Selected::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let owned = selected.retain(Default::default(), None).unwrap();
        let receipt = plan
            .correspond_literal(&owned, Default::default(), None)
            .unwrap();
        eprintln!(
            "NPT LOCAL CORRESPONDENCE quarter={quarter} position={} derivative={:?} edgepos={} edged={} query={} aggregate={} bytes={}",
            receipt.position_error_mm(),
            receipt.original_derivative_error_mm(),
            receipt.edge_position_error_mm(),
            receipt.edge_derivative_error_mm(),
            receipt.query_work_used(),
            receipt.aggregate_work_used(),
            receipt.accounted_peak_bytes()
        );
        assert!(receipt.position_error_mm() <= 1e-7);
        assert!(
            receipt
                .original_derivative_error_mm()
                .into_iter()
                .all(|e| e <= 1e-6)
        );
        assert!(receipt.edge_derivative_error_mm() <= 1e-6);
        for mapping in [
            receipt.tensor_mapping().as_slice(),
            receipt.edge_mapping().as_slice(),
            receipt.vertex_mapping().as_slice(),
        ] {
            let unique: std::collections::BTreeSet<_> = mapping.iter().copied().collect();
            assert_eq!(unique.len(), mapping.len());
        }
    }
}
