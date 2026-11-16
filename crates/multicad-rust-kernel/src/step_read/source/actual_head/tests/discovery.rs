//! Real neutral arenas and physical file mutations exercise discovery only;
//! an identifying cue never substitutes for the independent material reader.
use super::*;
fn shell(document: &CompoundDocument<'_>) -> (u32, Vec<u32>) {
    let solid = document.represented_members()[0].local_solid_entity();
    let [Value::String(_), Value::Reference(shell)] = document.records()[&solid][0].args.as_slice()
    else {
        panic!("actual solid");
    };
    let [Value::String(_), Value::List(faces)] = document.records()[shell][0].args.as_slice()
    else {
        panic!("actual shell");
    };
    (
        *shell,
        faces
            .iter()
            .map(|v| {
                if let Value::Reference(id) = v {
                    *id
                } else {
                    panic!("actual face reference")
                }
            })
            .collect(),
    )
}
#[test]
fn actual_descriptor_external_resident_overflow_and_full_resident_refuse_before_retry() {
    let input = mapped();
    for external in [usize::MAX, 64 * 1024 * 1024] {
        let options = OwnedHeadOptions {
            resources: crate::step_read::npt_stock::full_head::Options {
                additional_resident_bytes: external,
                ..Default::default()
            },
            ..Default::default()
        };
        assert!(matches!(
            read(input.clone(), options, None),
            Err(Error::Limit)
        ));
    }
    let actual = read(input.clone(), Default::default(), None).unwrap();
    assert_eq!(actual.actual_file(), input.as_ref());
    assert_eq!(actual.receipt().additional_resident_bytes, DESCRIPTOR_BYTES);
    assert!(actual.world_sewn_volume_mm3()[0] > 0.);
}
#[test]
fn actual_discovery_work_frontier_and_unmapped_whole_file_refuse_without_losing_original() {
    let local = candidate(false, 18.);
    let document = CompoundDocument::parse(&local, Default::default(), None).unwrap();
    assert!(document.represented_members()[0].mapping_link().is_none());
    let (cue, work) = recognized(&document, 16_000_000).unwrap();
    assert!(cue && work > 128);
    assert!(matches!(
        recognized(&document, work - 1),
        Err(ReadError::Limit)
    ));
    assert_eq!(recognized(&document, work).unwrap(), (true, work));
    assert!(matches!(
        read_source_compound(&local, 1e-7),
        Err(ReadError::Unsupported)
    ));
    assert!(
        read_actual_mapped_head_source(Arc::from(local.as_str()), Default::default(), None)
            .is_err()
    );
    // The mapping requirement is explicit; the same physical local file still
    // earns its independent local material through its actual supported API.
    let owned = OwnedActualHead::read(Arc::from(local.as_str()), Default::default(), None).unwrap();
    assert_eq!(owned.actual_file(), local);
    assert!(owned.local_sewn_volume_mm3()[0] > 0.);
    let mapped = mapped();
    assert!(matches!(
        read_source_compound(&mapped, 1e-7).unwrap().as_slice(),
        [SourceImport::ActualMappedHead(_)]
    ));
}
#[derive(Clone, Copy, Debug)]
enum Expected {
    Invalid,
    Unsupported,
    Cue(bool),
}
#[test]
fn actual_parsed_schema_discovery_refusals_and_false_cues_preserve_whole_reader_retry() {
    let input = mapped();
    let document = CompoundDocument::parse(&input, Default::default(), None).unwrap();
    let solid = document.represented_members()[0].local_solid_entity();
    let (shell, faces) = shell(&document);
    assert_eq!(faces.len(), 17);
    let (face, surface) = faces
        .iter()
        .find_map(|&id| {
            let f = &document.records()[&id][0];
            let [Value::String(_), _, Value::Reference(surface), _] = f.args.as_slice() else {
                return None;
            };
            (document.records()[surface][0].name == "B_SPLINE_SURFACE_WITH_KNOTS")
                .then_some((id, *surface))
        })
        .unwrap();
    let old_shell = document.records()[&shell][0].source.unwrap();
    let old_face = document.records()[&face][0].source.unwrap();
    let old_surface = document.records()[&surface][0].source.unwrap();
    let mut numbers = faces.iter().map(|id| format!("#{id}")).collect::<Vec<_>>();
    numbers[0] = "1.".into();
    let non_reference = format!("CLOSED_SHELL('',({}))", numbers.join(","));
    let wrong_degrees = input.replace(
        "B_SPLINE_SURFACE_WITH_KNOTS('',11,11",
        "B_SPLINE_SURFACE_WITH_KNOTS('',10,11",
    );
    assert_ne!(wrong_degrees, input.as_ref());
    let cases = [
        (
            "solid argument shape",
            change(&input, solid, "MANIFOLD_SOLID_BREP('broken')"),
            Expected::Invalid,
        ),
        (
            "shell kind",
            change(
                &input,
                shell,
                &old_shell.replacen("CLOSED_SHELL", "OPEN_SHELL", 1),
            ),
            Expected::Cue(false),
        ),
        (
            "shell argument shape",
            change(&input, shell, "CLOSED_SHELL('broken')"),
            Expected::Invalid,
        ),
        (
            "shell non-reference",
            change(&input, shell, &non_reference),
            Expected::Invalid,
        ),
        (
            "face argument shape",
            change(&input, face, "ADVANCED_FACE('broken')"),
            Expected::Invalid,
        ),
        (
            "non-advanced face",
            change(
                &input,
                face,
                &old_face.replacen("ADVANCED_FACE", "FACE_SURFACE", 1),
            ),
            Expected::Cue(true),
        ),
        (
            "complex surface",
            change(
                &input,
                surface,
                &format!("({old_surface}GEOMETRIC_REPRESENTATION_ITEM())"),
            ),
            Expected::Unsupported,
        ),
        ("other tensor degree", wrong_degrees, Expected::Cue(false)),
    ];
    let mut reached = 0;
    for (name, bytes, expected) in cases {
        assert_ne!(bytes, input.as_ref(), "{name}");
        // No parser record, source span, entity owner or certificate is forged.
        // Earlier independent parser refusals remain honest outcomes.
        match CompoundDocument::parse(&bytes, Default::default(), None) {
            Ok(arena) => {
                reached += 1;
                let actual = recognized(&arena, 16_000_000);
                assert!(
                    match (actual, expected) {
                        (Err(ReadError::Invalid), Expected::Invalid)
                        | (Err(ReadError::Unsupported), Expected::Unsupported) => true,
                        (Ok((cue, _)), Expected::Cue(want)) => cue == want,
                        _ => false,
                    },
                    "{name}: expected {expected:?}"
                );
            }
            Err(error) => {
                eprintln!("ACTUAL_DISCOVERY_EARLIER_PARSER_REFUSAL name={name} error={error:?}")
            }
        }
        assert!(
            read_actual_mapped_head_source(Arc::from(bytes), Default::default(), None).is_err(),
            "{name}"
        );
    }
    assert!(reached > 0);
    assert!(recognized(&document, 16_000_000).unwrap().0);
    drop(document);
    let retry = read(input.clone(), Default::default(), None).unwrap();
    assert_eq!(retry.actual_file(), input.as_ref());
    assert!(retry.world_sewn_volume_mm3()[0] > 0.);
}

#[test]
fn ordinary_analytic_box_and_actual_void_shell_are_not_head_discovery_cues() {
    use crate::orthogonal::{Box3, Operation, boolean};
    let mut writer = crate::step::Writer::new();
    writer
        .append_box(
            [2, 3, 4].map(|n| crate::Rational::new(n, 1).unwrap()),
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    let box_file = writer.finish().unwrap();
    let solid = boolean(
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
    let boundary = crate::orthogonal_boundary::build(solid.quads()).unwrap();
    let void_file = crate::step::write_orthogonal(&boundary, |p| p).unwrap();
    for (bytes, kind, expected_volume) in [
        (box_file, "MANIFOLD_SOLID_BREP", 24.),
        (void_file, "BREP_WITH_VOIDS", 1140.),
    ] {
        let document = CompoundDocument::parse(&bytes, Default::default(), None).unwrap();
        let member = document.represented_members()[0];
        assert_eq!(
            document.records()[&member.local_solid_entity()][0].name,
            kind
        );
        let (cue, work) = recognized(&document, 16_000_000).unwrap();
        assert!(!cue && work >= 128);
        drop(document);
        let mut members = read_source_compound(&bytes, 1e-7).unwrap();
        assert_eq!(members.len(), 1);
        let SourceImport::Analytic(member) = members.pop().unwrap() else {
            panic!("ordinary family must remain analytic");
        };
        assert_eq!(member.0.volume_mm3, expected_volume);
        assert_eq!(member.0.solid_count, 1);
    }
}

#[test]
fn actual_clipped_stock_without_head_stays_literal_after_discovery_and_owner_drop() {
    use crate::step_read::npt_stock::{Selected, full_head::source_tests::with_source};
    for quarter in [false, true] {
        let bytes = with_source(quarter, |source| {
            crate::step::npt_stock::write_source_candidate(
                source,
                16_000_000,
                64 * 1024 * 1024,
                None,
            )
            .unwrap()
            .into_string()
        });
        let document = CompoundDocument::parse(&bytes, Default::default(), None).unwrap();
        let member = document.represented_members()[0];
        assert!(!recognized(&document, 16_000_000).unwrap().0);
        let solid = member.local_solid_entity();
        let original = document.records()[&solid][0].source.unwrap().to_owned();
        let selected =
            Selected::parse(&document, member.item_entity(), Default::default(), None).unwrap();
        let independent = selected.retain(Default::default(), None).unwrap();
        let mut imported = read_source_compound(&bytes, 1e-7).unwrap();
        assert_eq!(imported.len(), 1);
        let SourceImport::LiteralNpt(actual) = imported.pop().unwrap() else {
            panic!("actual stock must remain literal, not become a head or recipe");
        };
        assert_eq!(actual.world_volume_mm3(), independent.world_volume_mm3());
        assert_eq!(
            actual.world_enclosure_mm(),
            independent.world_enclosure_mm()
        );
        drop(selected);
        drop(document);
        drop(bytes);
        drop(independent);
        assert_eq!(actual.original_component(solid, 0), Some(original.as_str()));
        assert!(actual.world_volume_mm3()[0] > 0.);
        assert!(
            actual
                .world_enclosure_mm()
                .into_iter()
                .flatten()
                .all(f64::is_finite)
        );
    }
    let input = mapped();
    assert!(matches!(
        read_source_compound(&input, 1e-7).unwrap().as_slice(),
        [SourceImport::ActualMappedHead(_)]
    ));
}
