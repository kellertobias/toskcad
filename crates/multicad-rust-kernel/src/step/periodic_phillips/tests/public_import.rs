//! Exercise actual public compound dispatch, not an already selected reader.
use super::*;
use crate::step_read::{ReadError, SourceImport, owned_insert::CompoundDocument};

fn ordinary(writer: &mut crate::step::Writer) {
    writer
        .append_box(
            [2, 3, 4].map(|v| crate::Rational::new(v, 1).unwrap()),
            [-20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
}

fn candidate(left: bool) -> (Material, [[f64; 3]; 4], QualifiedStep) {
    let material = source(left, 0.);
    let frame = occurrence(left);
    let file = export(&material, frame, Default::default(), None).unwrap();
    (material, frame, file)
}

#[test]
fn public_mixed_periodic_import_retains_actual_material_after_all_input_owners_drop() {
    let (members, expected) = {
        let mut writer = crate::step::Writer::new();
        let mut expected = Vec::new();
        for left in [false, true] {
            let (material, frame, file) = candidate(left);
            writer
                .append_qualified_periodic_phillips(&file, Default::default(), None)
                .unwrap();
            expected.push((material, frame, file.qualification().verification()));
            if !left {
                ordinary(&mut writer);
            }
        }
        let input = writer.finish().unwrap();
        let members = crate::step_read::read_source_compound(&input, 1e-7).unwrap();
        assert_eq!(members.len(), 3);
        (members, expected)
    };
    // Only the returned actual source owners survive; no parser, writer,
    // candidate packet, or mixed input String remains alive here.
    let mut expected = expected.into_iter();
    let mut owners = Vec::new();
    for (index, member) in members.into_iter().enumerate() {
        if index == 1 {
            let SourceImport::Analytic(value) = member else {
                panic!("The box must keep its independent analytic receipt")
            };
            assert_eq!(value.0.volume_mm3, 24.);
            assert_eq!(value.0.solid_count, 1);
            continue;
        }
        let SourceImport::PeriodicPhillips(material) = member else {
            panic!("The actual periodic member was not retained")
        };
        assert!(matches!(
            material.provenance(),
            step_read::PhillipsPeriodicMaterialProvenance::Decoded { .. }
        ));
        let input = material.original_input().unwrap();
        let document = CompoundDocument::parse(input, Default::default(), None).unwrap();
        assert_eq!(document.represented_members().len(), 1);
        assert!(material.original_curve_owners().is_some());
        assert!(material.thread_transfers().is_some());
        assert!(material.retained_budget_bytes() > input.len());
        let (source, frame, prior) = expected.next().unwrap();
        let proof = source
            .match_material(&material, frame, Default::default(), None)
            .unwrap();
        assert!(proof.maximum_world_distance_mm() <= 1e-7);
        assert!(proof.maximum_world_chart_derivative_mm() <= 1e-7);
        let volume = material.world_volume_mm3();
        assert!(volume[0] > 0. && volume[0] <= volume[1]);
        assert!(volume[0] <= prior.decoded_volume_bounds_mm3[1]);
        assert!(prior.decoded_volume_bounds_mm3[0] <= volume[1]);
        owners.push(material);
    }
    assert!(expected.next().is_none());
    assert_eq!(owners.len(), 2);
    // Equal nominal size/mass does not justify losing the physical handedness.
    let wrong = source(false, 0.);
    assert!(
        wrong
            .match_material(&owners[1], occurrence(true), Default::default(), None)
            .is_err()
    );
}

#[test]
fn public_mixed_periodic_import_rejects_orphan_owners_without_rewriting_input() {
    let (_, _, file) = candidate(false);
    let mut writer = crate::step::Writer::new();
    ordinary(&mut writer);
    writer
        .append_qualified_periodic_phillips(&file, Default::default(), None)
        .unwrap();
    let original = writer.finish().unwrap();
    let mut malformed = original.clone();
    malformed.insert_str(
        malformed.rfind("ENDSEC;").unwrap(),
        "#999999=VERTEX_POINT('',#1);\n",
    );
    let before = malformed.clone();
    assert!(matches!(
        crate::step_read::read_source_compound(&malformed, 1e-7),
        Err(ReadError::Invalid)
    ));
    assert_eq!(malformed, before);
    assert_eq!(
        crate::step_read::read_source_compound(&original, 1e-7)
            .unwrap()
            .len(),
        2
    );
    for tolerance in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(crate::step_read::read_source_compound(&original, tolerance).is_err());
    }
}

#[test]
fn public_periodic_import_reserves_all_actual_member_stages_before_material_publication() {
    let material = source(false, -8.);
    let file = export(&material, occurrence(false), Default::default(), None).unwrap();
    let mut writer = crate::step::Writer::new();
    for _ in 0..14 {
        writer
            .append_qualified_periodic_phillips(&file, Default::default(), None)
            .expect("the actual short-source append must fit existing encoding limits");
    }
    let input = writer.finish().unwrap();
    let before = input.clone();
    let document = CompoundDocument::parse(&input, Default::default(), None).unwrap();
    assert_eq!(document.represented_members().len(), 14);
    drop(document);
    // Each member is an actual independently qualified encoded solid with its
    // own rebased topology. Their staged sum exceeds this route's existing 4B
    // budget; a cheap discovery success must not publish partial material.
    assert!(matches!(
        crate::step_read::read_source_compound(&input, 1e-7),
        Err(ReadError::Limit)
    ));
    assert_eq!(input, before);
}

#[test]
fn public_periodic_import_keeps_the_existing_void_shell_member_route() {
    use crate::orthogonal::{Box3, Operation, boolean};
    let stock = boolean(
        &[Box3 {
            min: [0.; 3],
            max: [15., 10., 8.],
        }],
        &[Box3 {
            min: [2., 3., 1.],
            max: [7., 7., 4.],
        }],
        Operation::Subtract,
        1e-7,
    )
    .unwrap()
    .unwrap();
    let boundary = crate::orthogonal_boundary::build(stock.quads()).unwrap();
    let (_, _, file) = candidate(false);
    let mut writer = crate::step::Writer::new();
    writer.append_orthogonal(&boundary, |p| p).unwrap();
    writer
        .append_qualified_periodic_phillips(&file, Default::default(), None)
        .unwrap();
    let input = writer.finish().unwrap();
    assert!(input.contains("BREP_WITH_VOIDS"));
    let mut imported = crate::step_read::read_source_compound(&input, 1e-7).unwrap();
    assert_eq!(imported.len(), 2);
    assert!(matches!(
        imported.pop(),
        Some(SourceImport::PeriodicPhillips(_))
    ));
    let SourceImport::Analytic(voided) = imported.pop().unwrap() else {
        panic!("The void shell must preserve its existing analytic importer")
    };
    assert_eq!(voided.0.volume_mm3, 1140.);
    assert_eq!(voided.0.face_count, 12);
}

#[test]
fn public_periodic_import_does_not_admit_a_changed_interior_tensor_from_discovery() {
    let (_, _, file) = candidate(false);
    let mut writer = crate::step::Writer::new();
    ordinary(&mut writer);
    writer
        .append_qualified_periodic_phillips(&file, Default::default(), None)
        .unwrap();
    let input = writer.finish().unwrap();
    let document = CompoundDocument::parse(&input, Default::default(), None).unwrap();
    let local = document
        .local_member_document(
            document.represented_members()[1].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
    let body = step_read::read_phillips_body(local.as_str(), 1e-7, 1_000_000).unwrap();
    let point = body.literal_tensor_control_points(0).unwrap()[1][1];
    let old = input
        .lines()
        .find(|line| {
            let Some((_, args)) = line.split_once("=CARTESIAN_POINT('',(") else {
                return false;
            };
            let Some(args) = args.strip_suffix("));") else {
                return false;
            };
            let values = args
                .split(',')
                .map(str::parse::<f64>)
                .collect::<Result<Vec<_>, _>>();
            values.is_ok_and(|values| values.as_slice() == point)
        })
        .unwrap();
    let entity = old.split_once('=').unwrap().0;
    let shifted = [point[0] + 0.001, point[1], point[2]];
    let changed = input.replacen(
        old,
        &format!(
            "{entity}=CARTESIAN_POINT('',({:.16E},{:.16E},{:.16E}));",
            shifted[0], shifted[1], shifted[2]
        ),
        1,
    );
    let altered = CompoundDocument::parse(&changed, Default::default(), None).unwrap();
    let local = altered
        .local_member_document(
            altered.represented_members()[1].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
    let discovered = step_read::read_phillips_body(local.as_str(), 1e-7, 1_000_000).unwrap();
    assert_eq!(
        discovered.literal_tensor_control_points(0).unwrap()[1][1],
        shifted
    );
    assert_eq!(
        discovered.inferred_thread_spec(),
        body.inferred_thread_spec()
    );
    // Discovery retained the actual changed coefficient, but it cannot grant
    // the independent complete material/closure/mass proof required to import.
    assert!(matches!(
        crate::step_read::read_source_compound(&changed, 1e-7),
        Err(ReadError::Invalid | ReadError::Unsupported)
    ));
    assert_eq!(
        crate::step_read::read_source_compound(&input, 1e-7)
            .unwrap()
            .len(),
        2
    );
}
