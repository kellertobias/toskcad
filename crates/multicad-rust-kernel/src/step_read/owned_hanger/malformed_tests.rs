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
        ("EDGE_CURVE", "EDGE_CURVE('',$,$,$,.F.)"),
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
