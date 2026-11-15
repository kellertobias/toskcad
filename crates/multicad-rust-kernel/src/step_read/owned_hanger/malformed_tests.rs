use super::*;
use crate::authored_hanger_bolt::Hanger;

// Mutate actual serialized entities, so every rejection traverses the public
// Part21 parser and decoder rather than an invented admitted intermediate.
fn replace_component(input: &str, name: &str, replacement: &str) -> String {
    let arena = part21::parse(input).unwrap();
    let original = arena
        .values()
        .flatten()
        .find(|entity| entity.name == name)
        .unwrap()
        .source
        .unwrap();
    assert_ne!(original, replacement);
    input.replacen(original, replacement, 1)
}

// Select literal argument ranges without reserializing untouched coordinates.
// Input here is the already-parsed writer fixture; nested lists and quoted
// labels retain their original spellings.
fn argument(text: &str, index: usize) -> std::ops::Range<usize> {
    let open = text.find('(').unwrap();
    let bytes = text.as_bytes();
    let mut depth = 0;
    let mut quoted = false;
    let mut start = open + 1;
    (start..bytes.len())
        .filter_map(|i| {
            let b = bytes[i];
            if b == b'\'' {
                quoted = !quoted;
            } else if !quoted {
                if b == b'(' {
                    depth += 1;
                } else if b == b')' && depth > 0 {
                    depth -= 1;
                } else if (b == b',' || b == b')') && depth == 0 {
                    let range = start..i;
                    start = i + 1;
                    return Some(range);
                }
            }
            None
        })
        .nth(index)
        .unwrap()
}
fn replace_argument(text: &str, path: &[usize], replacement: &str) -> String {
    let selected = argument(text, path[0]);
    let changed = if path.len() == 1 {
        replacement.to_owned()
    } else {
        replace_argument(&text[selected.clone()], &path[1..], replacement)
    };
    let mut result = text.to_owned();
    result.replace_range(selected, &changed);
    result
}

#[test]
fn actual_hanger_schema_fields_refuse_invalid_geometry_and_topology() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let original = file.as_str();
    let arena = part21::parse(original).unwrap();
    let cases: &[(&str, &[usize], &str)] = &[
        ("B_SPLINE_CURVE", &[0], "1.5"),
        ("B_SPLINE_CURVE", &[0], "33"),
        ("B_SPLINE_CURVE", &[1], "0."),
        ("B_SPLINE_CURVE", &[1, 0], "0."),
        ("B_SPLINE_CURVE", &[1], "()"),
        ("B_SPLINE_CURVE", &[2], ".POLYLINE_FORM."),
        ("B_SPLINE_CURVE_WITH_KNOTS", &[0], "()"),
        ("B_SPLINE_CURVE_WITH_KNOTS", &[0, 0], "0"),
        ("B_SPLINE_CURVE_WITH_KNOTS", &[1, 1], "-1."),
        ("RATIONAL_B_SPLINE_CURVE", &[0, 0], "0."),
        ("B_SPLINE_SURFACE", &[0], "1.5"),
        ("B_SPLINE_SURFACE", &[0], "3"),
        ("B_SPLINE_SURFACE", &[2], "()"),
        ("B_SPLINE_SURFACE", &[2, 1], "()"),
        ("B_SPLINE_SURFACE", &[2, 0, 0], "0."),
        ("B_SPLINE_SURFACE_WITH_KNOTS", &[0], "()"),
        ("B_SPLINE_SURFACE_WITH_KNOTS", &[0, 0], "0"),
        ("B_SPLINE_SURFACE_WITH_KNOTS", &[1, 0], "0"),
        ("B_SPLINE_SURFACE_WITH_KNOTS", &[2, 0], "0.0000000000001"),
        ("B_SPLINE_SURFACE_WITH_KNOTS", &[3, 1], "-1."),
        ("RATIONAL_B_SPLINE_SURFACE", &[0, 0, 0], "0."),
        ("RATIONAL_B_SPLINE_SURFACE", &[0, 0], "()"),
        ("PLANE", &[1], "$"),
        ("PLANE", &[1], "#4294967294"),
        ("GEOMETRIC_REPRESENTATION_CONTEXT", &[0], "2.5"),
        ("GLOBAL_UNIT_ASSIGNED_CONTEXT", &[0], "()"),
        ("ADVANCED_BREP_SHAPE_REPRESENTATION", &[1], "()"),
        ("ADVANCED_BREP_SHAPE_REPRESENTATION", &[1], "0."),
        ("PRODUCT_DEFINITION_SHAPE", &[2], "0."),
        ("CLOSED_SHELL", &[1], "()"),
        ("ADVANCED_FACE", &[1], "()"),
        ("ADVANCED_FACE", &[3], ".U."),
        ("FACE_OUTER_BOUND", &[2], ".F."),
        ("EDGE_LOOP", &[1], "()"),
        ("ORIENTED_EDGE", &[1], "$"),
        ("ORIENTED_EDGE", &[4], ".U."),
        ("EDGE_CURVE", &[4], ".U."),
        ("SURFACE_CURVE", &[2], "()"),
        ("SURFACE_CURVE", &[3], ".PCURVE_S1."),
        ("DEFINITIONAL_REPRESENTATION", &[1], "()"),
        ("PCURVE", &[1], "#4294967294"),
        ("VERTEX_POINT", &[1], "$"),
    ];
    for &(name, path, replacement) in cases {
        let component = arena
            .values()
            .flatten()
            .find(|entity| entity.name == name)
            .unwrap()
            .source
            .unwrap();
        let changed = replace_argument(component, path, replacement);
        let input = original.replacen(component, &changed, 1);
        assert!(
            part21::parse(&input).is_ok(),
            "field mutation must remain valid Part21: {name}{path:?}"
        );
        let error = decode(&input, ReadOptions::default(), None).err();
        assert!(
            matches!(error, Some(Error::Invalid | Error::Unsupported)),
            "{name}{path:?}={replacement}: {error:?}"
        );
    }
    let shell = arena
        .values()
        .flatten()
        .find(|entity| entity.name == "CLOSED_SHELL")
        .unwrap();
    let Value::List(faces) = &shell.args[1] else {
        unreachable!()
    };
    let face_ids = faces
        .iter()
        .map(|value| match value {
            Value::Reference(id) => *id,
            _ => unreachable!(),
        })
        .collect::<Vec<_>>();
    let mut duplicate = face_ids.clone();
    duplicate.push(face_ids[0]);
    for ids in [&duplicate[..], &face_ids[..face_ids.len() - 1]] {
        let references = format!(
            "({})",
            ids.iter()
                .map(|id| format!("#{id}"))
                .collect::<Vec<_>>()
                .join(",")
        );
        let component = shell.source.unwrap();
        let changed = replace_argument(component, &[1], &references);
        let input = original.replacen(component, &changed, 1);
        assert!(part21::parse(&input).is_ok());
        assert!(matches!(
            decode(&input, Default::default(), None),
            Err(Error::Invalid)
        ));
    }
    let edge = arena
        .values()
        .flatten()
        .find(|entity| entity.name == "EDGE_CURVE")
        .unwrap();
    let Value::Reference(first_vertex) = edge.args[1] else {
        unreachable!()
    };
    let component = edge.source.unwrap();
    let changed = replace_argument(component, &[2], &format!("#{first_vertex}"));
    let input = original.replacen(component, &changed, 1);
    assert!(part21::parse(&input).is_ok());
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Unsupported)
    ));
    let (point_id, _) = arena
        .iter()
        .find(|(_, entities)| entities[0].name == "CARTESIAN_POINT")
        .unwrap();
    let plane = arena
        .values()
        .flatten()
        .find(|entity| entity.name == "PLANE")
        .unwrap()
        .source
        .unwrap();
    let changed = replace_argument(plane, &[1], &format!("#{point_id}"));
    let input = original.replacen(plane, &changed, 1);
    assert!(part21::parse(&input).is_ok());
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Unsupported)
    ));

    let complex = arena
        .values()
        .find(|entities| {
            entities.len() > 1 && entities.iter().any(|e| e.name == "RATIONAL_B_SPLINE_CURVE")
        })
        .unwrap();
    let rational = complex
        .iter()
        .find(|e| e.name == "RATIONAL_B_SPLINE_CURVE")
        .unwrap()
        .source
        .unwrap();
    let missing = original.replacen(rational, "", 1);
    assert!(part21::parse(&missing).is_ok());
    assert!(matches!(
        decode(&missing, Default::default(), None),
        Err(Error::Unsupported)
    ));
    let duplicate = rational.replacen("RATIONAL_B_SPLINE_CURVE", "B_SPLINE_CURVE", 1);
    let duplicate = original.replacen(rational, &duplicate, 1);
    // Duplicate complex component names are rejected by the shared parser,
    // before the Hanger-specific semantic layer can consume the record.
    assert!(part21::parse(&duplicate).is_err());
    assert!(matches!(
        decode(&duplicate, Default::default(), None),
        Err(Error::Invalid)
    ));

    let decoded = decode(original, Default::default(), None).unwrap();
    let point_2d = decoded
        .points()
        .values()
        .find(|p| p.coordinates().len() == 2)
        .unwrap()
        .entity();
    let curve = decoded
        .curves()
        .values()
        .find(|curve| {
            decoded.points()[&curve.point_entities()[0]]
                .coordinates()
                .len()
                == 3
        })
        .unwrap();
    let entities = &arena[&curve.entity()];
    let (name, controls) = if entities.len() == 1 {
        ("B_SPLINE_CURVE_WITH_KNOTS", 2)
    } else {
        ("B_SPLINE_CURVE", 1)
    };
    let component = entities
        .iter()
        .find(|e| e.name == name)
        .unwrap()
        .source
        .unwrap();
    let changed = replace_argument(component, &[controls, 0], &format!("#{point_2d}"));
    let input = original.replacen(component, &changed, 1);
    assert!(part21::parse(&input).is_ok());
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Invalid)
    ));
    let tensor = arena
        .values()
        .flatten()
        .find(|e| e.name == "B_SPLINE_SURFACE")
        .unwrap()
        .source
        .unwrap();
    let changed = replace_argument(tensor, &[2, 0, 0], &format!("#{point_2d}"));
    let input = original.replacen(tensor, &changed, 1);
    assert!(part21::parse(&input).is_ok());
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Invalid)
    ));
    assert!(
        !decode(original, ReadOptions::default(), None)
            .unwrap()
            .faces()
            .is_empty()
    );
}

#[test]
fn actual_hanger_material_binds_real_decoder_owners_and_rejects_invalid_options() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let decoded = decode(file.as_str(), Default::default(), None).unwrap();
    let other = decode(file.as_str(), Default::default(), None).unwrap();
    for options in [
        GridOptions {
            maximum_preparation_work: 0,
            ..Default::default()
        },
        GridOptions {
            maximum_preparation_work: 256_000_001,
            ..Default::default()
        },
        GridOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        GridOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            ExactPoleGrid::prepare(&decoded, options, None),
            Err(Error::Invalid)
        ));
    }
    let grid = ExactPoleGrid::prepare(&decoded, Default::default(), None).unwrap();
    for options in [
        PolarOptions {
            maximum_work: 0,
            ..Default::default()
        },
        PolarOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            PolarSupports::prepare(&grid, options, None),
            Err(Error::Invalid)
        ));
    }
    let polar = PolarSupports::prepare(&grid, Default::default(), None).unwrap();
    for options in [
        LiftOptions {
            position_limit_mm: f64::NAN,
            ..Default::default()
        },
        LiftOptions {
            derivative_limit_mm: 0.,
            ..Default::default()
        },
        LiftOptions {
            maximum_work: 0,
            ..Default::default()
        },
        LiftOptions {
            maximum_edge_work: 16_000_001,
            ..Default::default()
        },
        LiftOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        LiftOptions {
            maximum_subdivisions: 8193,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            decoded.surface_lifts(options, None),
            Err(Error::Invalid)
        ));
    }
    let lifts = decoded.surface_lifts(Default::default(), None).unwrap();
    let other_lifts = other.surface_lifts(Default::default(), None).unwrap();
    assert!(!std::ptr::eq(lifts.source(), other_lifts.source()));
    assert!(matches!(
        HealedMaterial::prepare(&polar, &other_lifts, Default::default(), None),
        Err(Error::Invalid)
    ));
    for options in [
        HealingOptions {
            maximum_work: 0,
            ..Default::default()
        },
        HealingOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        HealingOptions {
            maximum_sewing_error_mm: f64::INFINITY,
            ..Default::default()
        },
        HealingOptions {
            maximum_sewing_error_mm: 0.,
            ..Default::default()
        },
        HealingOptions {
            maximum_sewing_error_mm: 1e-6,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            HealedMaterial::prepare(&polar, &lifts, options, None),
            Err(Error::Invalid)
        ));
    }
    let valid = HealedMaterial::prepare(&polar, &lifts, Default::default(), None).unwrap();
    assert!(std::ptr::eq(valid.finite_lifts().source(), &decoded));
    assert!(valid.healed_volume_bounds_mm3()[0] > 0.);
    assert_eq!(decoded.literal_input(), file.as_str());
}

#[test]
fn malformed_actual_hanger_entities_never_publish_decoded_material() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let original = file.as_str();
    let admitted = decode(original, ReadOptions::default(), None).unwrap();
    assert!(!admitted.faces().is_empty());
    for options in [
        ReadOptions {
            maximum_work: 0,
            ..Default::default()
        },
        ReadOptions {
            retained_budget_bytes: 0,
            ..Default::default()
        },
        ReadOptions {
            maximum_work: 256_000_001,
            ..Default::default()
        },
        ReadOptions {
            retained_budget_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            decode(original, options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        ReadOptions {
            maximum_work: 1,
            ..Default::default()
        },
        ReadOptions {
            retained_budget_bytes: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(decode(original, options, None), Err(Error::Limit)));
    }
    assert!(matches!(
        decode(original, Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        decode(&original[..original.len() - 5], Default::default(), None),
        Err(Error::Invalid)
    ));
    let binary = replace_component(
        original,
        "CARTESIAN_POINT",
        "CARTESIAN_POINT('',(\"00\",0.,0.))",
    );
    assert!(matches!(
        decode(&binary, Default::default(), None),
        Err(Error::Unsupported)
    ));
    let nested = format!("CARTESIAN_POINT('',{}0{})", "(".repeat(80), ")".repeat(80));
    let nested = replace_component(original, "CARTESIAN_POINT", &nested);
    assert!(matches!(
        decode(&nested, Default::default(), None),
        Err(Error::Limit)
    ));
    // Structural admission retains arbitrary pole coordinates. A nonzero
    // decimal that caches as binary64 zero must retain its exact spelling and
    // must not acquire the exact-zero identity used by geometric validators.
    let underflow = replace_component(
        original,
        "CARTESIAN_POINT",
        "CARTESIAN_POINT('',(1.E-400,0.,0.))",
    );
    let changed = decode(&underflow, ReadOptions::default(), None).unwrap();
    let tiny = changed
        .points()
        .values()
        .flat_map(|point| point.coordinates())
        .find(|value| value.spelling() == "1.E-400")
        .unwrap();
    assert_eq!(tiny.value, 0.);
    assert!(!tiny.is_zero());

    // Keep the actual direction and change only its unit coordinate to a
    // distinct decimal with the same binary64 cache.
    let arena = part21::parse(original).unwrap();
    let direction = arena
        .values()
        .flatten()
        .find(|entity| entity.name == "DIRECTION")
        .unwrap();
    let unit = part21::numeric_spellings(direction)
        .unwrap()
        .into_iter()
        .find(|number| number.parse::<f64>().unwrap() == 1.)
        .unwrap();
    let direction_text = direction.source.unwrap();
    let alias = direction_text.replacen(unit, "1.0000000000000000000000001", 1);
    let alias = original.replacen(direction_text, &alias, 1);
    assert!(matches!(
        decode(&alias, ReadOptions::default(), None),
        Err(Error::Unsupported)
    ));
    for (name, replacement) in [
        ("CARTESIAN_POINT", "CARTESIAN_POINT('',(0.))"),
        ("CARTESIAN_POINT", "CARTESIAN_POINT()"),
        ("CARTESIAN_POINT", "CARTESIAN_POINT('',(.T.,0.,0.))"),
        ("B_SPLINE_CURVE", "B_SPLINE_CURVE()"),
        (
            "B_SPLINE_CURVE",
            "B_SPLINE_CURVE(0,(),.UNSPECIFIED.,.F.,.F.)",
        ),
        ("B_SPLINE_CURVE_WITH_KNOTS", "B_SPLINE_CURVE_WITH_KNOTS()"),
        ("RATIONAL_B_SPLINE_CURVE", "RATIONAL_B_SPLINE_CURVE((0.))"),
        (
            "B_SPLINE_SURFACE",
            "B_SPLINE_SURFACE(3,10,(),.UNSPECIFIED.,.F.,.F.,.F.)",
        ),
        (
            "RATIONAL_B_SPLINE_SURFACE",
            "RATIONAL_B_SPLINE_SURFACE((()))",
        ),
        ("PCURVE", "PCURVE('',$,$)"),
        ("PCURVE", "PCURVE()"),
        (
            "DEFINITIONAL_REPRESENTATION",
            "DEFINITIONAL_REPRESENTATION()",
        ),
        ("EDGE_CURVE", "EDGE_CURVE('',$,$,$,.F.)"),
        ("EDGE_CURVE", "EDGE_CURVE()"),
        ("SURFACE_CURVE", "SURFACE_CURVE()"),
        ("PLANE", "PLANE()"),
        ("PRODUCT_DEFINITION_SHAPE", "PRODUCT_DEFINITION_SHAPE()"),
        (
            "SHAPE_DEFINITION_REPRESENTATION",
            "SHAPE_DEFINITION_REPRESENTATION()",
        ),
        (
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            "ADVANCED_BREP_SHAPE_REPRESENTATION()",
        ),
        ("ADVANCED_FACE", "ADVANCED_FACE()"),
        (
            "GEOMETRIC_REPRESENTATION_CONTEXT",
            "GEOMETRIC_REPRESENTATION_CONTEXT()",
        ),
        ("AXIS2_PLACEMENT_3D", "AXIS2_PLACEMENT_3D('')"),
        (
            "GEOMETRIC_REPRESENTATION_CONTEXT",
            "GEOMETRIC_REPRESENTATION_CONTEXT(4)",
        ),
        ("SI_UNIT", "SI_UNIT(.CENTI.,.METRE.)"),
        ("DIRECTION", "DIRECTION('',(0.,0.,0.))"),
    ] {
        let changed = replace_component(original, name, replacement);
        let error = decode(&changed, ReadOptions::default(), None).err();
        assert!(
            matches!(error, Some(Error::Invalid | Error::Unsupported)),
            "{name}: {error:?}"
        );
    }
    // Failure never invalidates an earlier immutable result or its literal data.
    assert_eq!(admitted.literal_input(), original);
    assert_eq!(
        decode(original, ReadOptions::default(), None)
            .unwrap()
            .faces()
            .len(),
        admitted.faces().len()
    );
}
