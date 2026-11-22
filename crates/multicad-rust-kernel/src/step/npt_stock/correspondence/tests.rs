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
fn retained(plan: &SurfacePlan<'_>) -> OwnedMaterial {
    let file = write_source_candidate(plan, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let document =
        CompoundDocument::parse(file.as_str(), DocumentOptions::default(), None).unwrap();
    Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap()
}
#[test]
fn actual_correspondence_keeps_tolerance_resource_and_cancellation_controls() {
    let source = super::super::tests::source(false);
    let plan = SurfacePlan::new(
        CurvePlan::new(&source, Default::default(), None).unwrap(),
        Default::default(),
        None,
    )
    .unwrap();
    let owned = retained(&plan);
    let standard = CorrespondenceOptions::default();
    for options in [
        CorrespondenceOptions {
            position_error_mm: 0.,
            ..standard
        },
        CorrespondenceOptions {
            position_error_mm: f64::NAN,
            ..standard
        },
        CorrespondenceOptions {
            original_derivative_error_mm: 0.,
            ..standard
        },
        CorrespondenceOptions {
            original_derivative_error_mm: f64::INFINITY,
            ..standard
        },
        CorrespondenceOptions {
            maximum_work: 0,
            ..standard
        },
        CorrespondenceOptions {
            maximum_work: 96_000_001,
            ..standard
        },
        CorrespondenceOptions {
            maximum_bytes: 0,
            ..standard
        },
        CorrespondenceOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..standard
        },
    ] {
        assert!(matches!(
            plan.correspond_literal(&owned, options, None),
            Err(Error::Invalid)
        ));
    }
    let baseline = plan.correspond_literal(&owned, standard, None).unwrap();
    assert!(
        plan.correspond_literal(
            &owned,
            CorrespondenceOptions {
                maximum_work: baseline.aggregate_work_used(),
                maximum_bytes: baseline.accounted_peak_bytes(),
                ..standard
            },
            None
        )
        .is_ok()
    );
    for options in [
        CorrespondenceOptions {
            maximum_work: plan.work_used() + owned.work_used(),
            ..standard
        },
        CorrespondenceOptions {
            maximum_work: baseline.aggregate_work_used() - 1,
            ..standard
        },
        CorrespondenceOptions {
            maximum_bytes: baseline.accounted_peak_bytes() - 1,
            ..standard
        },
    ] {
        assert!(matches!(
            plan.correspond_literal(&owned, options, None),
            Err(Error::Limit)
        ));
    }
    for options in [
        CorrespondenceOptions {
            position_error_mm: 1e-12,
            ..standard
        },
        CorrespondenceOptions {
            original_derivative_error_mm: 1e-12,
            ..standard
        },
    ] {
        assert!(matches!(
            plan.correspond_literal(&owned, options, None),
            Err(Error::Unresolved)
        ));
    }
    // A genuinely smaller balance also stops during owner assignment, rather
    // than only its final publication. The continuing inherited ledger remains.
    let inherited = plan.work_used() + owned.work_used();
    assert!(matches!(
        plan.correspond_literal(
            &owned,
            CorrespondenceOptions {
                maximum_work: inherited + baseline.query_work_used() * 93 / 100,
                ..standard
            },
            None
        ),
        Err(Error::Limit)
    ));
    let calls = std::cell::Cell::new(0usize);
    plan.correspond_literal(
        &owned,
        standard,
        Some(&|| {
            calls.set(calls.get() + 1);
            false
        }),
    )
    .unwrap();
    for checkpoint in [1, calls.get() / 2, calls.get()] {
        let count = std::cell::Cell::new(0usize);
        assert!(matches!(
            plan.correspond_literal(
                &owned,
                standard,
                Some(&|| {
                    count.set(count.get() + 1);
                    count.get() >= checkpoint
                })
            ),
            Err(Error::Cancelled)
        ));
    }
}
#[test]
fn actual_shell_face_order_does_not_supply_source_patch_identity() {
    let source = super::super::tests::source(false);
    let plan = SurfacePlan::new(
        CurvePlan::new(&source, Default::default(), None).unwrap(),
        Default::default(),
        None,
    )
    .unwrap();
    let file = write_source_candidate(&plan, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let document =
        CompoundDocument::parse(file.as_str(), DocumentOptions::default(), None).unwrap();
    let parsed = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let line = file
        .as_str()
        .lines()
        .find(|line| line.contains("=CLOSED_SHELL("))
        .unwrap();
    let prefix = line.split_once('=').unwrap().0;
    let mut faces = parsed
        .parsed()
        .faces()
        .iter()
        .map(|face| face.entity())
        .collect::<Vec<_>>();
    let (first, second) = (0..8)
        .find_map(|first| {
            ((first + 1)..8)
                .find(|second| {
                    plan.patches()[first].spans() == plan.patches()[*second].spans()
                        && plan.patches()[first].source_patch().slab
                            == plan.patches()[*second].source_patch().slab
                })
                .map(|second| (first, second))
        })
        .expect("real repeated source-grid dimensions");
    // Change one slab's proposal order, rather than applying the same atlas
    // symmetry to both slabs. Broad numeric proposals must genuinely backtrack.
    faces.swap(first, second);
    let replacement = format!(
        "{prefix}=CLOSED_SHELL('',({}));",
        faces
            .iter()
            .map(|id| format!("#{id}"))
            .collect::<Vec<_>>()
            .join(",")
    );
    let reordered = file.as_str().replacen(line, &replacement, 1);
    let document = CompoundDocument::parse(&reordered, DocumentOptions::default(), None).unwrap();
    let owned = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap();
    let exact = plan
        .correspond_literal(&owned, Default::default(), None)
        .unwrap();
    assert_ne!(exact.tensor_mapping(), [0, 1, 2, 3, 4, 5, 6, 7]);
    // Deliberately broad requested error admits several numeric proposals;
    // the real owner graph must still select a consistent bijection.
    let broad = plan
        .correspond_literal(
            &owned,
            CorrespondenceOptions {
                position_error_mm: 1e6,
                original_derivative_error_mm: 1e6,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    assert_eq!(
        broad
            .edge_mapping()
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        20
    );
}
#[test]
fn actual_rotated_wire_does_not_forge_the_direct_source_parameter_witness() {
    let source = super::super::tests::source(false);
    let plan = SurfacePlan::new(
        CurvePlan::new(&source, Default::default(), None).unwrap(),
        Default::default(),
        None,
    )
    .unwrap();
    let file = write_source_candidate(&plan, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let line = file
        .as_str()
        .lines()
        .find(|line| line.contains("=EDGE_LOOP("))
        .unwrap();
    let open = line.find(",(").unwrap() + 2;
    let close = line.rfind("))").unwrap();
    let mut edges = line[open..close].split(',').collect::<Vec<_>>();
    edges.rotate_left(1);
    let changed = format!("{}{}{}", &line[..open], edges.join(","), &line[close..]);
    let text = file.as_str().replacen(line, &changed, 1);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let owned = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap();
    // The independent reader still earns actual material. The direct replay
    // witness additionally preserves the original UV-side/coedge origin;
    // equivalence through a different wire origin is not inferred here.
    assert!(owned.world_volume_mm3()[0] > 0.);
    assert!(matches!(
        plan.correspond_literal(&owned, Default::default(), None),
        Err(Error::Unresolved)
    ));
}
