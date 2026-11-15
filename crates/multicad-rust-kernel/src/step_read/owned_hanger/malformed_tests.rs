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
fn actual_tensor_weight_shapes_refuse_even_when_all_numeric_weights_survive() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let original = file.as_str();
    let arena = part21::parse(original).unwrap();
    let rational = arena
        .values()
        .flatten()
        .find(|e| e.name == "RATIONAL_B_SPLINE_SURFACE")
        .unwrap();
    let component = rational.source.unwrap();
    let matrix = &component[argument(component, 0)];
    let Value::List(weight_rows) = &rational.args[0] else {
        unreachable!()
    };
    let rows = weight_rows
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let Value::List(weights) = row else {
                unreachable!()
            };
            let row = &matrix[argument(matrix, i)];
            (0..weights.len())
                .map(|j| row[argument(row, j)].to_owned())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert!(rows.len() > 1 && rows[0].len() > 1);
    let count = rows.iter().map(Vec::len).sum::<usize>();
    for kind in 0..3 {
        let mut changed = rows.clone();
        let matrix = if kind == 0 {
            let moved = changed[0].pop().unwrap();
            changed[1].push(moved);
            format!(
                "({})",
                changed
                    .iter()
                    .map(|r| format!("({})", r.join(",")))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        } else if kind == 1 {
            let remainder = changed[0].split_off(1);
            changed[1].extend(remainder);
            let scalar = changed[0][0].clone();
            format!(
                "({scalar},{})",
                changed[1..]
                    .iter()
                    .map(|r| format!("({})", r.join(",")))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        } else {
            changed[0][0] = format!("({})", changed[0][0]);
            format!(
                "({})",
                changed
                    .iter()
                    .map(|r| format!("({})", r.join(",")))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        let changed = replace_argument(component, &[0], &matrix);
        let input = original.replacen(component, &changed, 1);
        let changed_arena = part21::parse(&input).unwrap();
        let entity = changed_arena
            .values()
            .flatten()
            .find(|e| e.name == "RATIONAL_B_SPLINE_SURFACE")
            .unwrap();
        fn numbers(v: &Value<'_>) -> usize {
            match v {
                Value::Number(_) => 1,
                Value::List(v) => v.iter().map(numbers).sum(),
                _ => 0,
            }
        }
        assert_eq!(
            numbers(&entity.args[0]),
            count,
            "numeric count must not be the refusal reason"
        );
        assert!(
            matches!(
                decode(&input, Default::default(), None),
                Err(Error::Invalid)
            ),
            "shape variant{kind}"
        );
    }
    assert!(decode(original, Default::default(), None).is_ok());
}

#[test]
fn actual_rectangular_tensor_grid_must_match_its_knot_dimensions() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let original = file.as_str();
    let arena = part21::parse(original).unwrap();
    let surface = arena
        .values()
        .flatten()
        .find(|e| e.name == "B_SPLINE_SURFACE")
        .unwrap();
    let Value::List(rows) = &surface.args[2] else {
        unreachable!()
    };
    let matrix = format!(
        "({})",
        rows.iter()
            .map(|row| {
                let Value::List(columns) = row else {
                    unreachable!()
                };
                assert!(columns.len() > 2);
                format!(
                    "({})",
                    columns[..columns.len() - 1]
                        .iter()
                        .map(|v| match v {
                            Value::Reference(id) => format!("#{id}"),
                            _ => unreachable!(),
                        })
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    );
    let component = surface.source.unwrap();
    let changed = replace_argument(component, &[2], &matrix);
    let input = original.replacen(component, &changed, 1);
    assert!(part21::parse(&input).is_ok());
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Invalid)
    ));
    assert!(decode(original, Default::default(), None).is_ok());
}

#[test]
fn actual_hanger_retained_memory_frontiers_refuse_and_recover() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let input = file.as_str();
    let full = decode(input, Default::default(), None).unwrap();
    let check = |bytes| {
        decode(
            input,
            ReadOptions {
                retained_budget_bytes: bytes,
                ..Default::default()
            },
            None,
        )
    };
    // Public requests traverse the actual parser, curves, tensors and final
    // retained graph. This budget is the documented retained-data accounting,
    // not a claim about the allocator's complete peak memory.
    let mut refused = 0;
    for step in 1..=32 {
        let bytes = ReadOptions::default().retained_budget_bytes * step / 32;
        match check(bytes) {
            Ok(decoded) => {
                assert_eq!(decoded.work_used(), full.work_used());
                assert_eq!(
                    decoded.retained_budget_bytes(),
                    full.retained_budget_bytes()
                );
            }
            Err(Error::Limit) => refused += 1,
            Err(error) => panic!("Unexpected byte-budget failure: {error:?}"),
        }
    }
    assert!(refused > 0);
    let mut lower = 1;
    let mut upper = ReadOptions::default().retained_budget_bytes;
    while upper - lower > 1 {
        let middle = lower + (upper - lower) / 2;
        match check(middle) {
            Ok(_) => upper = middle,
            Err(Error::Limit) => lower = middle,
            Err(error) => panic!("Unexpected frontier failure: {error:?}"),
        }
    }
    assert!(matches!(check(lower), Err(Error::Limit)));
    assert!(check(upper).is_ok());
    assert_eq!(
        decode(input, Default::default(), None)
            .unwrap()
            .retained_budget_bytes(),
        full.retained_budget_bytes()
    );
}

#[test]
fn actual_coedge_orientation_and_vertex_point_dimensions_are_checked() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let original = file.as_str();
    let arena = part21::parse(original).unwrap();
    let oriented = arena
        .values()
        .flatten()
        .find(|e| e.name == "ORIENTED_EDGE")
        .unwrap();
    let flipped = match oriented.args[4] {
        Value::Enum("T") => ".F.",
        Value::Enum("F") => ".T.",
        _ => unreachable!(),
    };
    let component = oriented.source.unwrap();
    let changed = replace_argument(component, &[4], flipped);
    let input = original.replacen(component, &changed, 1);
    assert!(part21::parse(&input).is_ok());
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Invalid)
    ));
    let Value::Reference(vertex_id) = oriented.args[3] else {
        unreachable!()
    };
    let edge = &arena[&vertex_id][0];
    assert_eq!(edge.name, "EDGE_CURVE");
    let Value::Reference(vertex_id) = edge.args[1] else {
        unreachable!()
    };
    // Preserve the complete coedge cycle while dangling its shared vertex:
    // every occurrence of that same vertex reference changes together, and
    // its original entity remains as unused input rather than being invented.
    let input = original
        .replace(&format!("#{vertex_id},"), "#4294967294,")
        .replace(&format!("#{vertex_id})"), "#4294967294)");
    assert_ne!(input, original);
    assert!(part21::parse(&input).is_ok());
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Invalid)
    ));

    // Keep every curve's original points and dimensions. Only the vertex's
    // binding is changed to a genuine existing pcurve point in this file.
    let planar_point = arena
        .iter()
        .find_map(|(&id, entities)| {
            let e = &entities[0];
            (e.name == "CARTESIAN_POINT" && matches!(&e.args[1], Value::List(v) if v.len() == 2))
                .then_some(id)
        })
        .unwrap();
    let vertex = arena
        .values()
        .flatten()
        .find(|e| e.name == "VERTEX_POINT")
        .unwrap();
    let component = vertex.source.unwrap();
    let changed = replace_argument(component, &[1], &format!("#{planar_point}"));
    let input = original.replacen(component, &changed, 1);
    assert!(part21::parse(&input).is_ok());
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Invalid)
    ));
    assert!(decode(original, Default::default(), None).is_ok());
}

#[test]
fn actual_bare_curve_record_cannot_replace_a_qualified_complex_curve() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let original = file.as_str();
    let arena = part21::parse(original).unwrap();
    let (&id, entities) = arena
        .iter()
        .find(|(_, entities)| entities.iter().any(|e| e.name == "B_SPLINE_CURVE"))
        .unwrap();
    assert!(entities.len() > 1);
    let base = entities
        .iter()
        .find(|e| e.name == "B_SPLINE_CURVE")
        .unwrap();
    let prefix = format!("#{id}=");
    let record = original
        .lines()
        .find(|line| line.starts_with(&prefix))
        .unwrap();
    let input = original.replacen(record, &format!("{prefix}{};", base.source.unwrap()), 1);
    let parsed = part21::parse(&input).unwrap();
    assert_eq!(parsed[&id].len(), 1);
    assert_eq!(parsed[&id][0].name, "B_SPLINE_CURVE");
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Unsupported)
    ));
    assert!(decode(original, Default::default(), None).is_ok());
}

#[test]
fn actual_two_closed_sources_cannot_be_claimed_as_one_spherical_shell() {
    let mut spec = crate::authored_hanger_bolt::tests::fixture(false, false);
    spec.metric = crate::authored_metric_thread::AuthoredMetricThread::new(
        crate::authored_metric_thread::Spec {
            length_mm: 1.,
            ..spec.metric.spec()
        },
    )
    .unwrap();
    spec.wood = crate::authored_flat_thread::AuthoredThread::new(
        crate::helical_thread::Spec {
            length_mm: 4.,
            ..spec.wood.spec()
        },
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-6,
        },
    )
    .unwrap();
    // Two genuine short authored sources fit the existing reader budget;
    // the default long pair is correctly refused earlier by its 64MiB cap.
    let source = Hanger::new(spec, None).unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let mut writer = crate::step::Writer::new();
    for _ in 0..2 {
        writer
            .append_qualified_hanger(&file, Default::default(), None)
            .unwrap();
    }
    let original = writer.finish().unwrap();
    let arena = part21::parse(&original).unwrap();
    let (&solid, entities) = arena
        .iter()
        .find(|(_, es)| es[0].name == "MANIFOLD_SOLID_BREP")
        .unwrap();
    let Value::Reference(shell) = entities[0].args[1] else {
        unreachable!()
    };
    let mut faces = Vec::new();
    let mut shell_count = 0;
    for entity in arena
        .values()
        .flatten()
        .filter(|e| e.name == "CLOSED_SHELL")
    {
        shell_count += 1;
        let Value::List(ids) = &entity.args[1] else {
            unreachable!()
        };
        faces.extend(ids.iter().map(|v| match v {
            Value::Reference(id) => format!("#{id}"),
            _ => unreachable!(),
        }));
    }
    assert_eq!(shell_count, 2);
    let component = arena[&shell][0].source.unwrap();
    let changed = replace_argument(component, &[1], &format!("({})", faces.join(",")));
    let input = original.replacen(component, &changed, 1);
    let representation = arena
        .values()
        .flatten()
        .find(|e| e.name == "ADVANCED_BREP_SHAPE_REPRESENTATION")
        .unwrap();
    let component = representation.source.unwrap();
    let changed = replace_argument(component, &[1], &format!("(#{solid})"));
    let input = input.replacen(component, &changed, 1);
    assert!(part21::parse(&input).is_ok());
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Unsupported)
    ));
    assert!(decode(file.as_str(), Default::default(), None).is_ok());
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
        ("DIRECTION", &[1], "(0.,1.)"),
        ("GEOMETRIC_REPRESENTATION_CONTEXT", &[0], "2.5"),
        ("GLOBAL_UNIT_ASSIGNED_CONTEXT", &[0], "()"),
        ("ADVANCED_BREP_SHAPE_REPRESENTATION", &[1], "()"),
        ("ADVANCED_BREP_SHAPE_REPRESENTATION", &[1], "0."),
        ("PRODUCT_DEFINITION_SHAPE", &[2], "0."),
        ("SHAPE_DEFINITION_REPRESENTATION", &[1], "#4294967294"),
        ("CLOSED_SHELL", &[1], "()"),
        ("ADVANCED_FACE", &[1], "()"),
        ("ADVANCED_FACE", &[1], "(#4294967294)"),
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

#[test]
fn actual_modified_supports_and_extreme_decimals_cannot_acquire_material_proofs() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let original = file.as_str();
    let decoded = decode(original, Default::default(), None).unwrap();
    let arena = part21::parse(original).unwrap();
    let first_plane = &decoded.planes()[0];
    let second_plane = &decoded.planes()[1];
    let origin = arena[&first_plane.origin_entity()][0].source.unwrap();
    let second_height = decoded.points()[&second_plane.origin_entity()].coordinates()[2].spelling();
    for (path, value, structural_refusal) in [
        (&[1, 0][..], "0.001", true),
        (&[1, 2][..], second_height, false),
    ] {
        let changed = replace_argument(origin, path, value);
        let input = original.replacen(origin, &changed, 1);
        let parsed = decode(&input, Default::default(), None);
        if structural_refusal {
            assert!(matches!(parsed, Err(Error::Unsupported)));
            continue;
        }
        let parsed = parsed.unwrap();
        let grid = ExactPoleGrid::prepare(&parsed, Default::default(), None).unwrap();
        assert!(matches!(
            PolarSupports::prepare(&grid, Default::default(), None),
            Err(Error::Unsupported)
        ));
    }
    let point_id = decoded.tensors()[0].point_entities()[0];
    let pole = arena[&point_id][0].source.unwrap();
    for (value, expected) in [
        ("1.E-16385", Error::Limit),
        ("1.E-2147483649", Error::Unsupported),
    ] {
        let changed = replace_argument(pole, &[1, 0], value);
        let input = original.replacen(pole, &changed, 1);
        let parsed = decode(&input, Default::default(), None).unwrap();
        assert!(
            matches!(ExactPoleGrid::prepare(&parsed, Default::default(), None), Err(error) if error == expected)
        );
    }
    let grid = ExactPoleGrid::prepare(&decoded, Default::default(), None).unwrap();
    let polar = PolarSupports::prepare(&grid, Default::default(), None).unwrap();
    assert!(std::ptr::eq(polar.grid().source(), &decoded));
}

#[test]
fn actual_scalar_tokens_and_representation_metadata_preserve_decoder_contract() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let file =
        crate::step::hanger::export_qualified_hanger(&source, Default::default(), None).unwrap();
    let original = file.as_str();
    let admitted = decode(original, Default::default(), None).unwrap();
    let arena = part21::parse(original).unwrap();
    let curve = arena
        .values()
        .flatten()
        .find(|e| e.name == "B_SPLINE_CURVE")
        .unwrap();
    let Value::List(points) = &curve.args[1] else {
        unreachable!()
    };
    let Value::Reference(point_id) = &points[0] else {
        unreachable!()
    };
    let point = arena[point_id][0].source.unwrap();
    for (path, replacement) in [(vec![0], "7."), (vec![1, 0], "1.E400")] {
        let changed = replace_argument(point, &path, replacement);
        let input = original.replacen(point, &changed, 1);
        assert!(matches!(
            decode(&input, Default::default(), None),
            Err(Error::Invalid)
        ));
    }
    for (name, path, replacement) in [
        ("B_SPLINE_CURVE_WITH_KNOTS", vec![1, 0], "'zero'"),
        ("B_SPLINE_SURFACE_WITH_KNOTS", vec![3, 0], "'zero'"),
        ("B_SPLINE_SURFACE", vec![6], ".T."),
    ] {
        let e = arena.values().flatten().find(|e| e.name == name).unwrap();
        let text = e.source.unwrap();
        let changed = replace_argument(text, &path, replacement);
        let input = original.replacen(text, &changed, 1);
        assert!(
            matches!(
                decode(&input, Default::default(), None),
                Err(Error::Invalid | Error::Unsupported)
            ),
            "{name}"
        );
    }
    let with = arena
        .values()
        .flatten()
        .find(|e| e.name == "B_SPLINE_SURFACE_WITH_KNOTS")
        .unwrap()
        .source
        .unwrap();
    let truncated = format!("{})", &with[..with.rfind(',').unwrap()]);
    let input = original.replacen(with, &truncated, 1);
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Invalid)
    ));

    let rep = arena
        .values()
        .flatten()
        .find(|e| e.name == "ADVANCED_BREP_SHAPE_REPRESENTATION")
        .unwrap();
    let Value::Reference(context) = &rep.args[2] else {
        unreachable!()
    };
    let dimension = arena[context]
        .iter()
        .find(|e| e.name == "GEOMETRIC_REPRESENTATION_CONTEXT")
        .unwrap()
        .source
        .unwrap();
    assert_eq!(original.matches(dimension).count(), 1);
    let changed = replace_argument(dimension, &[0], "4");
    let input = original.replacen(dimension, &changed, 1);
    assert!(matches!(
        decode(&input, Default::default(), None),
        Err(Error::Unsupported)
    ));

    let first = arena.keys().copied().max().unwrap() + 1;
    let axis = first + 3;
    let extra = format!(
        "#{first}=CARTESIAN_POINT('',(0.,0.,0.));\n#{}=DIRECTION('',(0.,0.,1.));\n#{}=DIRECTION('',(1.,0.,0.));\n#{axis}=AXIS2_PLACEMENT_3D('',#{first},#{},#{});\n",
        first + 1,
        first + 2,
        first + 1,
        first + 2
    );
    let text = rep.source.unwrap();
    let range = argument(text, 1);
    let items = format!("{},#{axis})", &text[range.start..range.end - 1]);
    let changed = replace_argument(text, &[1], &items);
    let extended = original.replacen(text, &changed, 1).replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("{extra}ENDSEC;\nEND-ISO-10303-21;"),
        1,
    );
    let retained = decode(&extended, Default::default(), None).unwrap();
    assert_eq!(retained.faces().len(), admitted.faces().len());
    assert_eq!(retained.edges().len(), admitted.edges().len());
    assert_eq!(retained.curves().len(), admitted.curves().len());
    assert_eq!(retained.solid_entity(), admitted.solid_entity());
    let axis_text = format!(
        "AXIS2_PLACEMENT_3D('',#{first},#{},#{})",
        first + 1,
        first + 2
    );
    assert!(extended.contains(&axis_text));
    let invalid = extended.replacen(&axis_text, "AXIS2_PLACEMENT_3D('')", 1);
    assert!(matches!(
        decode(&invalid, Default::default(), None),
        Err(Error::Unsupported)
    ));
    let target = format!("#{first}=CARTESIAN_POINT('',(0.,0.,0.));");
    let invalid = extended.replacen(
        &target,
        &format!("#{first}=CARTESIAN_POINT('',(1.,0.,0.));"),
        1,
    );
    assert!(matches!(
        decode(&invalid, Default::default(), None),
        Err(Error::Unsupported)
    ));
    assert_eq!(admitted.literal_input(), original);
}
